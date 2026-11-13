"""Synthetic periodic voxel weave and small-strain thermoelastic FE homogenization.

Coordinates follow the project convention: X and Z in plane, Y normal.
Stiffness is GPa, CTE is ppm/K, and strain uses engineering shear.
These cells are geometry experiments, not reconstructions of a real fabric.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, cg, splu

from .model import Inputs, filled_matrix, ply
from .woven import rotate_thermoelastic, _rz, _rx, _ry


def isotropic_phase(e: float, nu: float, alpha: float):
    lam = e * nu / ((1 + nu) * (1 - 2 * nu))
    mu = e / (2 * (1 + nu))
    c = np.diag([lam + 2 * mu] * 3 + [mu] * 3)
    c[:3, :3] += lam * (np.ones((3, 3)) - np.eye(3))
    return c, np.array([alpha, alpha, alpha, 0., 0., 0.])


@lru_cache(maxsize=None)
def _element_quadrature(n: int):
    # Unit cell, n elements per side. Local node order is binary (x,y,z).
    signs = np.array([(x, y, z) for x in (-1, 1) for y in (-1, 1)
                      for z in (-1, 1)], dtype=float)
    quadrature = []
    for qx in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
        for qy in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
            for qz in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
                q = np.array([qx, qy, qz])
                grad = np.empty((8, 3))
                for a, s in enumerate(signs):
                    for d in range(3):
                        grad[a, d] = n * s[d] / 4 * np.prod(
                            1 + s[np.arange(3) != d] * q[np.arange(3) != d])
                b = np.zeros((6, 24))
                for a, (dx, dy, dz) in enumerate(grad):
                    k = 3 * a
                    b[0, k], b[1, k + 1], b[2, k + 2] = dx, dy, dz
                    b[3, k + 1], b[3, k + 2] = dz, dy
                    b[4, k], b[4, k + 2] = dz, dx
                    b[5, k], b[5, k + 1] = dy, dx
                quadrature.append(b)
    return np.array(quadrature)


def _dofs(n: int, i: int, j: int, k: int):
    nodes = [(((i + x) % n * n + (j + y) % n) * n + (k + z) % n)
             for x in (0, 1) for y in (0, 1) for z in (0, 1)]
    return (3 * np.repeat(nodes, 3) + np.tile(np.arange(3), 8)).astype(int)


def homogenize(cells_c: np.ndarray, cells_a: np.ndarray):
    """Periodic Q1 hexahedron FE; six mechanical solves and one thermal solve.

    Cellwise C and alpha must be (n,n,n,6,6)/(n,n,n,6), respectively.
    Returns homogenized C, alpha and diagnostics. The thermal RHS is built
    from C*alpha; average stress is evaluated by volume integration.
    """
    c = np.asarray(cells_c, float)
    a = np.asarray(cells_a, float)
    n = c.shape[0]
    if c.shape != (n, n, n, 6, 6) or a.shape != (n, n, n, 6):
        raise ValueError("expected cubic voxel arrays of C and alpha")
    if not np.allclose(c, c.swapaxes(-1, -2), atol=1e-9):
        raise ValueError("cell stiffness must be symmetric")
    bs = _element_quadrature(n)
    vol = 1 / n**3
    bmean = bs.mean(axis=0)
    # All B matrices have equal Gauss weights; element volume already includes detJ.
    entries = n**3 * 24 * 24
    rows = np.empty(entries, dtype=np.int32)
    cols = np.empty(entries, dtype=np.int32)
    vals = np.empty(entries, dtype=float)
    cursor = 0
    rhs = np.zeros((3 * n**3, 7))
    elements = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                ce, ae = c[i, j, k], a[i, j, k] * 1e-6
                dof = _dofs(n, i, j, k)
                ke = vol * np.einsum('qai,ab,qbj->ij', bs, ce, bs) / 8
                fe = -vol * (bmean.T @ ce)
                loc = np.empty((24, 7))
                loc[:, :6] = fe
                loc[:, 6] = -fe @ ae
                np.add.at(rhs, dof, loc)
                rows[cursor:cursor+576] = np.repeat(dof, 24)
                cols[cursor:cursor+576] = np.tile(dof, 24)
                vals[cursor:cursor+576] = ke.ravel()
                cursor += 576
                elements.append((dof, ce, ae))
    ndof = 3 * n**3
    kk = coo_matrix((vals, (rows, cols)), shape=(ndof, ndof)).tocsc()
    free = np.arange(3, ndof)
    u = np.zeros((ndof, 7))
    reduced = kk[free, :][:, free]
    solver_info = []
    if n <= 8:
        u[free] = splu(reduced).solve(rhs[free])
        solver_info = [0] * 7
    else:
        # Sparse LU can exhaust Windows' limited virtual-memory headroom even
        # for 12^3 cells. CG uses the SPD periodic stiffness without fill-in.
        reduced = reduced.tocsr()
        diagonal = reduced.diagonal()
        preconditioner = LinearOperator(reduced.shape, matvec=lambda x: x/diagonal)
        for col in range(7):
            sol, info = cg(reduced, rhs[free, col], M=preconditioner,
                           rtol=1e-9, atol=1e-12, maxiter=10000)
            if info != 0:
                raise RuntimeError(f"periodic CG did not converge: RHS {col}, info={info}")
            u[free, col] = sol
            solver_info.append(info)
    avg = np.zeros((6, 7))
    for dof, ce, ae in elements:
        fluctuation = bmean @ u[dof]
        strain = fluctuation + np.column_stack((np.eye(6), np.zeros(6)))
        strain[:, 6] -= ae
        avg += vol * ce @ strain
    ceff = (avg[:, :6] + avg[:, :6].T) / 2
    alpha = -np.linalg.solve(ceff, avg[:, 6]) * 1e6
    voigt = c.mean(axis=(0, 1, 2))
    reuss = np.linalg.inv(np.linalg.inv(c).mean(axis=(0, 1, 2)))
    return ceff, alpha, {
        "n": n, "symmetry_error": float(np.max(np.abs(avg[:, :6] - avg[:, :6].T))),
        "min_eigenvalue": float(np.linalg.eigvalsh(ceff).min()),
        "voigt_gap_min_eigenvalue": float(np.linalg.eigvalsh(voigt - ceff).min()),
        "reuss_gap_min_eigenvalue": float(np.linalg.eigvalsh(ceff - reuss).min()),
        "mean_phase_cte": a.mean(axis=(0, 1, 2)).tolist(),
        "solver": "sparse_LU" if n <= 8 else "CG_Jacobi",
        "max_relative_equilibrium_residual": float(max(
            np.linalg.norm(reduced @ u[free, j] - rhs[free, j]) /
            max(np.linalg.norm(rhs[free, j]), 1e-30) for j in range(7))),
    }


def properties(c: np.ndarray, a: np.ndarray):
    s = np.linalg.inv(c)
    return {"EX": 1/s[0, 0], "EY": 1/s[1, 1], "EZ": 1/s[2, 2],
            "GXY": 1/s[5, 5], "GYZ": 1/s[3, 3], "GXZ": 1/s[4, 4],
            "NUXY": -s[0, 1]/s[0, 0], "NUYZ": -s[1, 2]/s[1, 1],
            "NUXZ": -s[0, 2]/s[0, 0],
            "AX": a[0], "AY": a[1], "AZ": a[2]}


@dataclass(frozen=True)
class CellSpec:
    n: int = 8
    amplitude: float = 0.10
    tow_half_height: float = 0.12
    tow_half_width: float = 0.23
    local_fiber_fraction: float = 0.65
    global_matrix_silica_fraction: float = 0.40
    silica_in_tows: float | None = None  # None means uniform matrix silica
    matrix_method: str = "mori_tanaka"


def synthetic_cell(p: Inputs, spec: CellSpec):
    """Alternating plain weave; voxel centers classify elliptical tow volumes.

    The two tow families weave in Y. They are continuous across X/Z periodic
    boundaries. Silica conservation is enforced using the *actual* voxel tow
    fraction; included/excluded silica is balanced in the outside matrix.
    """
    n = spec.n
    xyz = (np.arange(n) + .5) / n
    phase = np.zeros((n, n, n), dtype=np.int8)  # 0 outside, 1 warp, 2 weft
    slopes = np.zeros((n, n, n))
    for i, x in enumerate(xyz):
        for j, y in enumerate(xyz):
            for k, z in enumerate(xyz):
                warp_options = []
                for center, sign in ((.25, 1), (.75, -1)):
                    dz = abs(z - center)
                    yc = .5 + spec.amplitude * np.sin(2*np.pi*x) * sign
                    metric = ((y-yc)/spec.tow_half_height)**2 + (dz/spec.tow_half_width)**2
                    warp_options.append((metric, spec.amplitude*2*np.pi*np.cos(2*np.pi*x)*sign))
                weft_options = []
                for center, sign in ((.25, 1), (.75, -1)):
                    dx = abs(x - center)
                    yc = .5 - spec.amplitude * np.sin(2*np.pi*z) * sign
                    metric = ((y-yc)/spec.tow_half_height)**2 + (dx/spec.tow_half_width)**2
                    weft_options.append((metric, -spec.amplitude*2*np.pi*np.cos(2*np.pi*z)*sign))
                wm, ws = min(warp_options)
                fm, fs = min(weft_options)
                if min(wm, fm) < 1:
                    phase[i, j, k] = 1 if wm <= fm else 2
                    slopes[i, j, k] = ws if wm <= fm else fs
    vt = float(np.count_nonzero(phase) / phase.size)
    vf = vt * spec.local_fiber_fraction
    if not 0 < vt < 1:
        raise ValueError("voxel resolution missed all tows")
    target_silica = spec.global_matrix_silica_fraction * (1-vf)
    tow_matrix = vt * (1-spec.local_fiber_fraction)
    phi_in = (spec.global_matrix_silica_fraction if spec.silica_in_tows is None
              else spec.silica_in_tows)
    phi_out = (target_silica - tow_matrix*phi_in)/(1-vt)
    if not (0 <= phi_in <= 1 and 0 <= phi_out <= 1):
        raise ValueError("silica partition is not volume feasible")
    cp, ap, _ = ply(replace(p, vf=spec.local_fiber_fraction, phi=phi_in),
                    matrix_method=spec.matrix_method)
    em, nu, am = filled_matrix(replace(p, phi=phi_out), spec.matrix_method)
    cm, at = isotropic_phase(em, nu, am)
    cells_c = np.empty((n, n, n, 6, 6))
    cells_a = np.empty((n, n, n, 6))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                typ = phase[i, j, k]
                if typ == 0:
                    cells_c[i,j,k], cells_a[i,j,k] = cm, at
                else:
                    theta = np.arctan(slopes[i,j,k])
                    rotation = (_rz(theta) if typ == 1 else
                                _rx(-theta) @ _ry(-np.pi/2))
                    cells_c[i,j,k], cells_a[i,j,k] = rotate_thermoelastic(cp, ap, rotation)
    return cells_c, cells_a, {"tow_fraction": vt, "fiber_fraction": vf,
                             "silica_fraction": target_silica,
                             "resin_fraction": 1-vf-target_silica,
                             "silica_in_tows": phi_in, "silica_outside_tows": phi_out,
                             "warp_voxels": int(np.count_nonzero(phase == 1)),
                             "weft_voxels": int(np.count_nonzero(phase == 2)),
                             "rms_tow_angle_deg": float(np.rad2deg(np.sqrt(
                                 np.mean(np.arctan(slopes[phase > 0])**2))))}
