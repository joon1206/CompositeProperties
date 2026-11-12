"""Free thermal response of a bare, uniformly layered PCB laminate.

Classical laminate plate theory with a through-thickness stack. Units: mm,
GPa, ppm/K. A copper/dielectric build-up is represented by an in-plane
isostrain effective layer; its copper coverage is a screening input.
No die, interposer, component, stiffener or ring is present.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import Inputs, evaluate


def isotropic_q(e: float, nu: float) -> np.ndarray:
    if e <= 0 or not -1 < nu < 0.5:
        raise ValueError("invalid isotropic elastic constants")
    return e/(1-nu**2)*np.array([[1, nu, 0], [nu, 1, 0],
                                 [0, 0, (1-nu)/2]], dtype=float)


def core_q_alpha(p: Inputs = Inputs(), matrix_method: str = "legacy"):
    """Map repository X/Z in-plane directions to the PCB x/y plane."""
    result = evaluate(p, method="exact", matrix_method=matrix_method)
    indices = [0, 2, 4]  # xx, zz, xz in the repository's 3D Voigt order
    compliance = np.linalg.solve(result["C"], np.eye(6))
    q = np.linalg.solve(compliance[np.ix_(indices, indices)], np.eye(3))
    return q, result["alpha"][indices]


@dataclass(frozen=True)
class Layer:
    name: str
    thickness_mm: float
    q_gpa: np.ndarray
    alpha_ppm_k: np.ndarray

    def checked(self):
        q = np.asarray(self.q_gpa, dtype=float)
        alpha = np.asarray(self.alpha_ppm_k, dtype=float)
        if self.thickness_mm <= 0 or q.shape != (3, 3) or alpha.shape != (3,):
            raise ValueError("layer needs positive thickness, 3x3 Q and 3-vector alpha")
        if not np.allclose(q, q.T) or np.linalg.eigvalsh(q).min() <= 0:
            raise ValueError("layer Q must be symmetric positive definite")
        return self


def solve_laminate(layers: list[Layer], delta_temp_k: float, side_mm: float):
    """Solve free midplane strain and curvature of a uniform flat stack.

    Layers are ordered bottom to top. Curvature order is kxx, kyy, kxy where
    kxy=2*d²w/dxdy. Bow is peak-to-valley of its quadratic free shape.
    Initial bow, residual stress, cure shrinkage and supports are excluded.
    """
    if not layers or side_mm <= 0:
        raise ValueError("need layers and positive board side")
    for layer in layers:
        layer.checked()
    total = sum(layer.thickness_mm for layer in layers)
    z = -total/2
    a, b, d = (np.zeros((3, 3)) for _ in range(3))
    nt, mt = np.zeros(3), np.zeros(3)
    for layer in layers:
        z0, z1 = z, z+layer.thickness_mm
        q, alpha = layer.q_gpa, layer.alpha_ppm_k
        dz, dz2, dz3 = z1-z0, (z1*z1-z0*z0)/2, (z1**3-z0**3)/3
        a += q*dz
        b += q*dz2
        d += q*dz3
        thermal = q@alpha*(delta_temp_k*1e-6)
        nt += thermal*dz
        mt += thermal*dz2
        z = z1
    system = np.block([[a, b], [b, d]])
    response = np.linalg.solve(system, np.r_[nt, mt])
    eps0, kappa = response[:3], response[3:]
    # An unconstrained square plate has quadratic shape under uniform
    # curvature. Sample its corners, side centers and center; extrema of a
    # general quadratic may also lie on an edge, so use a dense grid.
    coords = np.linspace(-side_mm/2, side_mm/2, 301)
    xx, yy = np.meshgrid(coords, coords)
    w = .5*kappa[0]*xx**2+.5*kappa[1]*yy**2+.5*kappa[2]*xx*yy
    return {"midplane_strain": eps0, "curvature_per_mm": kappa,
            "bow_mm": float(np.ptp(w)), "shape_mm": w,
            "thermal_force": nt, "thermal_moment": mt,
            "A": a, "B": b, "D": d}


@dataclass(frozen=True)
class PCBScreenInputs:
    side_mm: float = 100.0  # illustrative panel size; replace with actual outline
    core_thickness_mm: float = 1.6
    build_up_pairs_each_side: int = 10  # project deck's 10 or 14 pairs
    dielectric_ply_mm: float = .030  # project deck
    copper_ply_mm: float = .012  # project deck
    top_copper_coverage: float = .75  # project deck's nominal metal ratio
    bottom_copper_coverage: float = .75
    dielectric_e_gpa: float = 20.0  # hypothetical effective value, not measured
    dielectric_nu: float = .30
    dielectric_cte_ppm_k: float = 10.0  # hypothetical, not measured
    copper_e_gpa: float = 110.0  # illustrative, not measured for this board
    copper_nu: float = .34
    copper_cte_ppm_k: float = 17.0
    delta_temp_k: float = -100.0  # normalized cooling increment, not process history

    def checked(self):
        if self.side_mm <= 0 or self.core_thickness_mm <= 0:
            raise ValueError("positive panel size and core thickness required")
        if self.build_up_pairs_each_side < 0 or int(self.build_up_pairs_each_side) != self.build_up_pairs_each_side:
            raise ValueError("build-up pair count must be nonnegative integer")
        if min(self.dielectric_ply_mm, self.copper_ply_mm) <= 0:
            raise ValueError("build-up ply thicknesses must be positive")
        if not (0 <= self.top_copper_coverage <= 1 and 0 <= self.bottom_copper_coverage <= 1):
            raise ValueError("copper coverage must lie in [0,1]")
        isotropic_q(self.dielectric_e_gpa, self.dielectric_nu)
        isotropic_q(self.copper_e_gpa, self.copper_nu)
        return self


def effective_build_up(inputs: PCBScreenInputs, coverage: float, name: str) -> Layer:
    """In-plane Voigt mixture across copper patterns and dielectric.

    The assumed copper volume share is copper thickness * area coverage /
    total build-up thickness. This is a screening approximation, not a
    validated trace-pattern homogenization or a local stress model.
    """
    p = inputs.checked()
    if p.build_up_pairs_each_side == 0:
        raise ValueError("no build-up pairs specified")
    total = p.build_up_pairs_each_side*(p.dielectric_ply_mm+p.copper_ply_mm)
    copper_share = p.copper_ply_mm*coverage/(p.dielectric_ply_mm+p.copper_ply_mm)
    qd = isotropic_q(p.dielectric_e_gpa, p.dielectric_nu)
    qc = isotropic_q(p.copper_e_gpa, p.copper_nu)
    q = (1-copper_share)*qd+copper_share*qc
    beta = ((1-copper_share)*qd@np.array([p.dielectric_cte_ppm_k]*2+[0.])
            + copper_share*qc@np.array([p.copper_cte_ppm_k]*2+[0.]))
    return Layer(name, total, q, np.linalg.solve(q, beta))


def screen_pcb(core_cte_ppm_k: float, inputs: PCBScreenInputs = PCBScreenInputs(),
               core_q: np.ndarray | None = None):
    p = inputs.checked()
    q = core_q_alpha()[0] if core_q is None else np.asarray(core_q, dtype=float)
    layers = [Layer("core", p.core_thickness_mm, q,
                    np.array([core_cte_ppm_k]*2+[0.]))]
    if p.build_up_pairs_each_side:
        layers = [effective_build_up(p, p.bottom_copper_coverage, "bottom build-up"),
                  *layers,
                  effective_build_up(p, p.top_copper_coverage, "top build-up")]
    result = solve_laminate(layers, p.delta_temp_k, p.side_mm)
    result["layers"] = layers
    return result
