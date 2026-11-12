import numpy as np
import pytest

from composite_research.model import Inputs, evaluate, filled_matrix
from composite_research.woven import rotate_thermoelastic, woven_orientation


def test_zero_crimp_reproduces_consistent_cross_ply():
    expected = evaluate(Inputs(), method="voigt")
    actual = woven_orientation(Inputs())
    for key in ("EX", "EY", "EZ", "GXY", "GYZ", "GXZ",
                "NUXY", "NUYZ", "NUXZ", "AX", "AY", "AZ"):
        assert actual[key] == pytest.approx(expected[key], rel=1e-12)


def test_rotations_preserve_isotropic_material_and_thermal_strain():
    e, nu, a = 25., .3, 12.
    g = e/(2*(1+nu))
    lam = e*nu/((1+nu)*(1-2*nu))
    c = np.zeros((6, 6))
    c[:3, :3] = lam
    c[np.diag_indices(3)] = lam+2*g
    c[3, 3] = c[4, 4] = c[5, 5] = g
    alpha = np.array([a, a, a, 0., 0., 0.])
    theta = np.deg2rad(17)
    rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0], [0, 0, 1.]])
    cr, ar = rotate_thermoelastic(c, alpha, rotation)
    np.testing.assert_allclose(cr, c, atol=1e-13)
    np.testing.assert_allclose(ar, alpha, atol=1e-13)


def test_symmetric_crimp_has_no_net_shear_cte_and_is_positive_definite():
    for method in ("legacy", "mori_tanaka"):
        for warp, weft in ((0, 0), (5, 15), (15, 5), (20, 20)):
            r = woven_orientation(Inputs(), matrix_method=method,
                                  warp_crimp_deg=warp, weft_crimp_deg=weft)
            assert np.linalg.eigvalsh(r["C"]).min() > 0
            np.testing.assert_allclose(r["C"], r["C"].T, atol=1e-12)
            np.testing.assert_allclose(r["alpha"][3:], 0, atol=1e-12)
            assert r["EX"] > 0 and r["EZ"] > 0


def test_warp_weft_crimp_exchange_swaps_in_plane_properties():
    a = woven_orientation(Inputs(), matrix_method="mori_tanaka",
                          warp_crimp_deg=5, weft_crimp_deg=15)
    b = woven_orientation(Inputs(), matrix_method="mori_tanaka",
                          warp_crimp_deg=15, weft_crimp_deg=5)
    assert a["EX"] == pytest.approx(b["EZ"])
    assert a["EZ"] == pytest.approx(b["EX"])
    assert a["AX"] == pytest.approx(b["AZ"])
    assert a["AZ"] == pytest.approx(b["AX"])


def test_measured_filled_matrix_route_matches_identical_estimate():
    p = Inputs()
    estimated = woven_orientation(p, matrix_method="legacy", warp_crimp_deg=7,
                                  weft_crimp_deg=12)
    direct = woven_orientation(p, matrix_method="measured",
                               filled_matrix_properties=filled_matrix(p, "legacy"),
                               warp_crimp_deg=7, weft_crimp_deg=12)
    for key in ("EX", "EY", "EZ", "AX", "AY", "AZ"):
        assert direct[key] == pytest.approx(estimated[key], rel=1e-12)
