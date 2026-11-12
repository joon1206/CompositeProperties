"""Orientation-based surrogate for woven composite elastic and CTE properties.

The two tow families are approximated by transversely isotropic UD material
segments tilted by symmetric crimp angles out of the laminate plane. Their
3D stiffness and thermal stress are Voigt-averaged. This is a transparent
upper-bound-style orientation model, not a weave-unit-cell homogenization.
Repository axes: X/Z in plane, Y through thickness. GPa and ppm/K.
"""
from __future__ import annotations

import numpy as np

from .model import Inputs, ply


_PAIRS = ((0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1))


def rotate_thermoelastic(c: np.ndarray, alpha: np.ndarray, rotation: np.ndarray):
    """Rotate engineering-Voigt stiffness and eigenstrain consistently."""
    c = np.asarray(c, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    r = np.asarray(rotation, dtype=float)
    if c.shape != (6, 6) or alpha.shape != (6,) or r.shape != (3, 3):
        raise ValueError("expected 6x6 stiffness, 6-vector CTE and 3x3 rotation")
    if not np.allclose(r.T@r, np.eye(3), atol=1e-12) or np.linalg.det(r) < .999999999:
        raise ValueError("rotation must be proper orthogonal")
    c4 = np.empty((3, 3, 3, 3))
    for i, (a, b) in enumerate(_PAIRS):
        for j, (d, e) in enumerate(_PAIRS):
            c4[a, b, d, e] = c4[b, a, d, e] = c4[a, b, e, d] = c4[b, a, e, d] = c[i, j]
    rotated4 = np.einsum("ia,jb,kc,ld,abcd->ijkl", r, r, r, r, c4,
                         optimize=True)
    rotated_c = np.array([[rotated4[a, b, d, e] for d, e in _PAIRS]
                          for a, b in _PAIRS])
    rotated_c = (rotated_c+rotated_c.T)/2
    a2 = np.zeros((3, 3))
    for i, (a, b) in enumerate(_PAIRS):
        a2[a, b] = a2[b, a] = alpha[i] if a == b else alpha[i]/2
    rotated_a2 = r@a2@r.T
    rotated_alpha = np.array([rotated_a2[a, b] if a == b
                              else 2*rotated_a2[a, b] for a, b in _PAIRS])
    return rotated_c, rotated_alpha


def _rz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def _rx(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1., 0, 0], [0, c, -s], [0, s, c]])


def _ry(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1., 0], [-s, 0, c]])


def woven_orientation(p: Inputs = Inputs(), *, matrix_method="legacy",
                      warp_crimp_deg=0.0, weft_crimp_deg=0.0,
                      warp_share=0.5, filled_matrix_properties=None):
    """Compute 9 engineering elastic constants and 3 CTEs.

    Crimp angles tilt tow axes into Y (thickness). Equal positive/negative
    segments cancel net coupling to the thickness direction. `warp_share`
    weights the X family; the Z family receives the remainder. At zero
    crimp and 50/50 share, this equals the existing consistent Voigt
    cross-ply. Fiber and matrix inputs are otherwise unchanged.
    """
    if not (0 <= warp_share <= 1):
        raise ValueError("warp_share must be in [0,1]")
    if not (0 <= warp_crimp_deg < 45 and 0 <= weft_crimp_deg < 45):
        raise ValueError("crimp angles must be in [0,45) degrees")
    if matrix_method == "measured":
        if filled_matrix_properties is None:
            raise ValueError("measured matrix method requires (E, nu, CTE)")
        c, alpha, ply_details = ply(p, matrix_override=filled_matrix_properties)
    else:
        if filled_matrix_properties is not None:
            raise ValueError("matrix override requires matrix_method='measured'")
        c, alpha, ply_details = ply(p, matrix_method)
    tw, tf = np.deg2rad([warp_crimp_deg, weft_crimp_deg])
    rotations = ((warp_share/2, _rz(tw)),
                 (warp_share/2, _rz(-tw)),
                 ((1-warp_share)/2, _rx(tf)@_ry(-np.pi/2)),
                 ((1-warp_share)/2, _rx(-tf)@_ry(-np.pi/2)))
    ce, beta = np.zeros((6, 6)), np.zeros(6)
    for weight, rotation in rotations:
        cr, ar = rotate_thermoelastic(c, alpha, rotation)
        ce += weight*cr
        beta += weight*(cr@ar)
    ae = np.linalg.solve(ce, beta)
    s = np.linalg.solve(ce, np.eye(6))
    return {"EX": 1/s[0, 0], "EY": 1/s[1, 1], "EZ": 1/s[2, 2],
            "GXY": 1/s[5, 5], "GYZ": 1/s[3, 3], "GXZ": 1/s[4, 4],
            "NUXY": -s[0, 1]/s[0, 0],
            "NUYZ": -s[1, 2]/s[1, 1],
            "NUXZ": -s[0, 2]/s[0, 0],
            "AX": ae[0], "AY": ae[1], "AZ": ae[2],
            "C": ce, "alpha": ae, "ply": ply_details}
