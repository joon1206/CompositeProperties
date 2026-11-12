"""Reduced-order thermal warpage screen for a square core, die and edge ring.

This is a free, small-deflection plate model, not a package finite-element
model. It minimizes plane-stress thermoelastic energy with smooth symmetric
Ritz displacement fields. Lengths are mm, moduli GPa, CTEs ppm/K.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import Inputs, evaluate


def isotropic_q(e: float, nu: float) -> np.ndarray:
    if e <= 0 or not -1 < nu < 0.5:
        raise ValueError("invalid isotropic elastic constants")
    return e / (1 - nu**2) * np.array(
        [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]], dtype=float
    )


def core_q_alpha(p: Inputs = Inputs(), matrix_method: str = "legacy"):
    """Map repository X/Z in-plane directions to the package x/y plane."""
    result = evaluate(p, method="exact", matrix_method=matrix_method)
    indices = [0, 2, 4]  # xx, zz, xz in the repository's 3D Voigt order
    compliance = np.linalg.solve(result["C"], np.eye(6))
    q = np.linalg.solve(compliance[np.ix_(indices, indices)], np.eye(3))
    alpha = result["alpha"][indices]
    return q, alpha


@dataclass(frozen=True)
class ScreenInputs:
    package_side: float = 100.0
    interposer_side: float = 70.0
    die_area_fraction: float = 0.95
    core_thickness: float = 1.6
    build_up_thickness: float = 0.0  # each side; zero until effective properties are known
    build_up_e: float | None = None
    build_up_nu: float | None = None
    build_up_cte: float | None = None
    die_thickness: float = 0.775
    die_e: float = 130.0  # illustrative effective value; not in project deck
    die_nu: float = 0.28
    die_cte: float = 4.4  # project deck's effective die CTE
    ring_width: float = 10.0
    ring_thickness: float = 2.0
    ring_e: float = 200.0  # illustrative steel value
    ring_nu: float = 0.30
    ring_cte: float = 12.0  # illustrative steel value
    stress_free_temp: float = 150.0  # project deck assumption
    observation_temp: float = 25.0

    def checked(self):
        if min(self.package_side, self.interposer_side, self.core_thickness,
               self.die_thickness, self.ring_width, self.ring_thickness) <= 0:
            raise ValueError("lengths must be positive")
        if self.build_up_thickness < 0:
            raise ValueError("build-up thickness cannot be negative")
        if self.build_up_thickness and (self.build_up_e is None or self.build_up_nu is None
                                        or self.build_up_cte is None):
            raise ValueError("effective build-up E, nu, and CTE are required")
        if self.build_up_thickness:
            isotropic_q(self.build_up_e, self.build_up_nu)
        if not (0 < self.die_area_fraction <= 1):
            raise ValueError("die area fraction must be in (0, 1]")
        if self.interposer_side > self.package_side or 2 * self.ring_width >= self.package_side:
            raise ValueError("die and ring must fit inside package")
        isotropic_q(self.die_e, self.die_nu)
        isotropic_q(self.ring_e, self.ring_nu)
        return self


# Each basis is a pair of in-plane polynomials u/L, v/L, or a polynomial w/L.
# A monomial dictionary uses (power of X, power of Y): coefficient.
_MEMBRANE = (
    ({(1, 0): 1}, {(0, 1): 1}),
    ({(3, 0): 1}, {(0, 3): 1}),
    ({(1, 2): 1}, {(2, 1): 1}),
    ({(5, 0): 1}, {(0, 5): 1}),
    ({(3, 2): 1}, {(2, 3): 1}),
    ({(1, 4): 1}, {(4, 1): 1}),
    ({(7, 0): 1}, {(0, 7): 1}),
    ({(5, 2): 1}, {(2, 5): 1}),
    ({(3, 4): 1}, {(4, 3): 1}),
    ({(1, 6): 1}, {(6, 1): 1}),
)
_BENDING = (
    {(2, 0): 1, (0, 2): 1},
    {(4, 0): 1, (0, 4): 1},
    {(2, 2): 1},
    {(6, 0): 1, (0, 6): 1},
    {(4, 2): 1, (2, 4): 1},
    {(8, 0): 1, (0, 8): 1},
    {(6, 2): 1, (2, 6): 1},
    {(4, 4): 1},
)


def _poly_derivative(poly, dx, dy, x, y):
    total = 0.0
    for (px, py), coefficient in poly.items():
        if px < dx or py < dy:
            continue
        factor = coefficient
        for i in range(dx):
            factor *= px - i
        for i in range(dy):
            factor *= py - i
        total += factor * x ** (px - dx) * y ** (py - dy)
    return total


def _strain_matrices(x: float, y: float, side: float, membrane, bending):
    X, Y = x / side, y / side
    n = len(membrane) + len(bending)
    h, k = np.zeros((3, n)), np.zeros((3, n))
    for j, (u, v) in enumerate(membrane):
        h[:, j] = (_poly_derivative(u, 1, 0, X, Y),
                   _poly_derivative(v, 0, 1, X, Y),
                   _poly_derivative(u, 0, 1, X, Y) + _poly_derivative(v, 1, 0, X, Y))
    for j, w in enumerate(bending, start=len(membrane)):
        k[:, j] = (_poly_derivative(w, 2, 0, X, Y) / side,
                   _poly_derivative(w, 0, 2, X, Y) / side,
                   2 * _poly_derivative(w, 1, 1, X, Y) / side)
    return h, k


def _layer_resultants(layers, delta_temp):
    a, b, d = (np.zeros((3, 3)) for _ in range(3))
    nt, mt = np.zeros(3), np.zeros(3)
    for z0, z1, q, alpha in layers:
        dz = z1 - z0
        dz2 = (z1*z1 - z0*z0)/2
        dz3 = (z1**3 - z0**3)/3
        a += q * dz
        b += q * dz2
        d += q * dz3
        thermal = q @ alpha * (delta_temp * 1e-6)
        nt += thermal * dz
        mt += thermal * dz2
    return a, b, d, nt, mt


def _quadrature_breaks(p: ScreenInputs):
    die_half = p.interposer_side * np.sqrt(p.die_area_fraction) / 2
    ring_inner = p.package_side / 2 - p.ring_width
    half = p.package_side / 2
    return sorted(set([-half, -ring_inner, -die_half, die_half, ring_inner, half]))


def solve_screen(core_cte: float, core_q: np.ndarray | None = None,
                 ring: bool = False, inputs: ScreenInputs = ScreenInputs(),
                 order: int = 8):
    """Return free thermal shape and bow for the stated illustrative stack.

    The die patch is bonded to the core or to an optional symmetric effective
    build-up. Interposer, explicit copper, underfill, EMC, cure shrinkage and
    local attachments are omitted. A steel ring occupies the outer frame.
    """
    p = inputs.checked()
    if order not in (2, 4, 6, 8):
        raise ValueError("order must be 2, 4, 6, or 8")
    n_membrane, n_bending = {2: (1, 1), 4: (3, 3), 6: (6, 5), 8: (10, 8)}[order]
    membrane, bending = _MEMBRANE[:n_membrane], _BENDING[:n_bending]
    if ring and p.interposer_side*np.sqrt(p.die_area_fraction) > p.package_side-2*p.ring_width:
        raise ValueError("die footprint overlaps ring")
    q_core = core_q_alpha()[0] if core_q is None else np.asarray(core_q, dtype=float)
    if q_core.shape != (3, 3) or not np.allclose(q_core, q_core.T) or np.linalg.eigvalsh(q_core).min() <= 0:
        raise ValueError("core_q must be a positive-definite symmetric 3x3 matrix")
    q_die, q_ring = isotropic_q(p.die_e, p.die_nu), isotropic_q(p.ring_e, p.ring_nu)
    delta_temp = p.observation_temp - p.stress_free_temp
    zc = p.core_thickness / 2
    core = (-zc, zc, q_core, np.array([core_cte]*2+[0.0]))
    common = [core]
    if p.build_up_thickness:
        qb = isotropic_q(p.build_up_e, p.build_up_nu)
        ab = np.array([p.build_up_cte]*2+[0.0])
        common += [(-zc-p.build_up_thickness, -zc, qb, ab),
                   (zc, zc+p.build_up_thickness, qb, ab)]
    top = zc+p.build_up_thickness
    die = (top, top+p.die_thickness, q_die, np.array([p.die_cte]*2+[0.0]))
    rim = (top, top+p.ring_thickness, q_ring, np.array([p.ring_cte]*2+[0.0]))
    local = {
        "bare": _layer_resultants(common, delta_temp),
        "die": _layer_resultants([*common, die], delta_temp),
        "ring": _layer_resultants([*common, rim], delta_temp),
    }
    n = len(membrane) + len(bending)
    stiffness, force = np.zeros((n, n)), np.zeros(n)
    breaks = _quadrature_breaks(p)
    points, weights = np.polynomial.legendre.leggauss(order//2+3)
    die_half = p.interposer_side*np.sqrt(p.die_area_fraction)/2
    ring_inner = p.package_side/2-p.ring_width
    for xa, xb in zip(breaks[:-1], breaks[1:]):
        for ya, yb in zip(breaks[:-1], breaks[1:]):
            cx, cy = (xa+xb)/2, (ya+yb)/2
            if abs(cx) < die_half and abs(cy) < die_half:
                region = "die"
            elif ring and (abs(cx) > ring_inner or abs(cy) > ring_inner):
                region = "ring"
            else:
                region = "bare"
            a, b, d, nt, mt = local[region]
            for xi, wi in zip(points, weights):
                x = cx+(xb-xa)*xi/2
                for yi, wj in zip(points, weights):
                    y = cy+(yb-ya)*yi/2
                    h, k = _strain_matrices(x, y, p.package_side, membrane, bending)
                    weight = wi*wj*(xb-xa)*(yb-ya)/4
                    stiffness += weight*(h.T@a@h+h.T@b@k+k.T@b@h+k.T@d@k)
                    force += weight*(h.T@nt+k.T@mt)
    coefficients = np.linalg.solve(stiffness, force)
    grid = np.linspace(-p.package_side/2, p.package_side/2, 101)
    xx, yy = np.meshgrid(grid, grid)
    shape = np.zeros_like(xx)
    for value, poly in zip(coefficients[len(membrane):], bending):
        shape += p.package_side*value*sum(c*(xx/p.package_side)**i*(yy/p.package_side)**j
                                           for (i, j), c in poly.items())
    shape -= shape[50, 50]  # center as display reference; bow is translation invariant
    bow = float(np.ptp(shape))
    corner_minus_center = float(shape[0, 0])
    return {"bow_mm": bow, "corner_minus_center_mm": corner_minus_center,
            "center_curvature_per_mm": float(2*coefficients[len(membrane)]/p.package_side),
            "coefficients": coefficients, "x_mm": grid, "shape_mm": shape}
