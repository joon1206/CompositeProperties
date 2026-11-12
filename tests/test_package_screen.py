from dataclasses import replace

import numpy as np
import pytest

from composite_research.package_screen import (
    ScreenInputs, _layer_resultants, core_q_alpha, isotropic_q, solve_screen,
)


def test_repository_axis_mapping_and_positive_stiffness():
    q, alpha = core_q_alpha()
    assert np.linalg.eigvalsh(q).min() > 0
    assert alpha[0] == pytest.approx(alpha[1])
    assert alpha[2] == pytest.approx(0)


def test_zero_thermal_mismatch_is_flat_even_with_ring():
    q = isotropic_q(32, .25)
    p = replace(ScreenInputs(), die_cte=4.4, ring_cte=4.4,
                build_up_thickness=.42, build_up_e=20, build_up_nu=.3,
                build_up_cte=4.4)
    for ring in (False, True):
        result = solve_screen(4.4, q, ring, p)
        assert result["bow_mm"] < 1e-11
    assert solve_screen(1.75, q, False, replace(p, observation_temp=p.stress_free_temp))["bow_mm"] < 1e-11
    with pytest.raises(ValueError, match="build-up"):
        solve_screen(4.4, q, False, replace(p, build_up_e=None))


def test_full_coverage_bilayer_matches_exact_abd_solution():
    p = replace(ScreenInputs(), interposer_side=100.0, die_area_fraction=1.0)
    q_core = isotropic_q(32.0, .25)
    numerical = solve_screen(1.75, q_core, False, p)
    zc = p.core_thickness/2
    layers = [(-zc, zc, q_core, np.array([1.75, 1.75, 0.])),
              (zc, zc+p.die_thickness, isotropic_q(p.die_e, p.die_nu),
               np.array([p.die_cte, p.die_cte, 0.]))]
    a, b, d, nt, mt = _layer_resultants(layers, p.observation_temp-p.stress_free_temp)
    exact = np.linalg.solve(np.block([[a, b], [b, d]]), np.r_[nt, mt])
    kappa = exact[3]
    assert exact[4] == pytest.approx(kappa)
    assert numerical["center_curvature_per_mm"] == pytest.approx(kappa, rel=1e-9)
    assert numerical["corner_minus_center_mm"] == pytest.approx(kappa*p.package_side**2/4, rel=1e-9)


def test_ring_can_change_preferred_cte_in_this_screen():
    q = isotropic_q(32, .25)
    p = ScreenInputs()
    no_ring = [solve_screen(a, q, False, p)["bow_mm"] for a in (4.4, 10.0)]
    with_ring = [solve_screen(a, q, True, p)["bow_mm"] for a in (4.4, 10.0)]
    assert no_ring[0] < no_ring[1]
    assert with_ring[1] < with_ring[0]


def test_symmetric_build_up_can_reverse_catalog_scenario_ranking():
    q, _ = core_q_alpha()
    bare = ScreenInputs()
    build_up = replace(bare, build_up_thickness=.42, build_up_e=20,
                       build_up_nu=.3, build_up_cte=10)
    bare_bow = [solve_screen(a, q, False, bare)["bow_mm"] for a in (1.75, 5.5)]
    built_bow = [solve_screen(a, q, False, build_up)["bow_mm"] for a in (1.75, 5.5)]
    assert bare_bow[1] < bare_bow[0]
    assert built_bow[0] < built_bow[1]


def test_polynomial_order_convergence_for_screening_cases():
    q, _ = core_q_alpha()
    p = ScreenInputs()
    for ring in (False, True):
        for cte in (1.75, 5.5):
            b6 = solve_screen(cte, q, ring, p, order=6)["bow_mm"]
            b8 = solve_screen(cte, q, ring, p, order=8)["bow_mm"]
            assert abs(b8-b6)/b8 < .08
