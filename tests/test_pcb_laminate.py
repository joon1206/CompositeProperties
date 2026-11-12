from dataclasses import replace

import numpy as np
import pytest

from composite_research.pcb_laminate import (
    Layer, PCBScreenInputs, core_q_alpha, isotropic_q, screen_pcb,
    solve_laminate,
)


def test_bare_symmetric_board_is_flat_for_both_catalog_scenarios():
    for core_cte in (1.75, 5.5, 6.53):
        result = screen_pcb(core_cte)
        assert result["bow_mm"] < 1e-11
        np.testing.assert_allclose(result["curvature_per_mm"], 0, atol=1e-13)
        assert result["midplane_strain"][0] != 0


def test_asymmetric_copper_reverses_curvature_when_faces_swap():
    p = replace(PCBScreenInputs(), top_copper_coverage=.9,
                bottom_copper_coverage=.6)
    first = screen_pcb(1.75, p)
    second = screen_pcb(1.75, replace(p, top_copper_coverage=.6,
                                     bottom_copper_coverage=.9))
    np.testing.assert_allclose(first["curvature_per_mm"],
                               -second["curvature_per_mm"], rtol=1e-10, atol=1e-14)
    assert first["bow_mm"] == pytest.approx(second["bow_mm"])


def test_uniform_cte_and_zero_temperature_change_have_no_bow():
    p = replace(PCBScreenInputs(), top_copper_coverage=1.,
                bottom_copper_coverage=.2, dielectric_cte_ppm_k=5.,
                copper_cte_ppm_k=5.)
    assert screen_pcb(5., p)["bow_mm"] < 1e-11
    assert screen_pcb(1.75, replace(p, delta_temp_k=0.))["bow_mm"] < 1e-11


def test_equal_stiffness_bimetal_analytic_curvature():
    q = isotropic_q(30, .25)
    bottom = Layer("bottom", .5, q, np.array([2., 2., 0.]))
    top = Layer("top", .5, q, np.array([8., 8., 0.]))
    result = solve_laminate([bottom, top], -100., 100.)
    expected = 3*(8-2)*(-100)*1e-6/(2*1.0)
    np.testing.assert_allclose(result["curvature_per_mm"], [expected, expected, 0],
                               rtol=1e-12, atol=1e-14)


def test_bow_scales_with_board_size_and_core_stiffness_is_valid():
    q, alpha = core_q_alpha()
    assert np.linalg.eigvalsh(q).min() > 0
    assert alpha[0] == pytest.approx(alpha[1])
    p = replace(PCBScreenInputs(), top_copper_coverage=.9,
                bottom_copper_coverage=.6)
    b100 = screen_pcb(1.75, p)["bow_mm"]
    b250 = screen_pcb(1.75, replace(p, side_mm=250.))["bow_mm"]
    assert b250/b100 == pytest.approx((250/100)**2)
