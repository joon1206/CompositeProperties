import numpy as np
import pytest

from composite_research.model import Inputs
from composite_research.weave_cell import CellSpec, homogenize, isotropic_phase, synthetic_cell


def test_periodic_homogenization_recovers_homogeneous_thermoelasticity():
    c, a = isotropic_phase(10., .3, 20.)
    ce, ae, diag = homogenize(np.broadcast_to(c, (4,4,4,6,6)),
                              np.broadcast_to(a, (4,4,4,6)))
    np.testing.assert_allclose(ce, c, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ae, a, rtol=1e-12, atol=1e-12)
    assert diag["symmetry_error"] < 1e-10


@pytest.mark.parametrize("phi_in", [None, 0., .7])
def test_synthetic_weave_conserves_composition(phi_in):
    spec = CellSpec(n=8, silica_in_tows=phi_in)
    c, a, info = synthetic_cell(Inputs(), spec)
    assert c.shape == (8,8,8,6,6)
    assert a.shape == (8,8,8,6)
    assert info["fiber_fraction"] + info["silica_fraction"] + info["resin_fraction"] == pytest.approx(1.)
    expected_silica = (info["tow_fraction"]*(1-spec.local_fiber_fraction)*info["silica_in_tows"]
                       + (1-info["tow_fraction"])*info["silica_outside_tows"])
    assert expected_silica == pytest.approx(info["silica_fraction"])


def test_nonuniform_periodic_cell_obeys_energy_bounds():
    c, a, _ = synthetic_cell(Inputs(), CellSpec(n=8, silica_in_tows=0.))
    ce, ae, info = homogenize(c, a)
    assert np.linalg.eigvalsh(ce).min() > 0
    assert info["symmetry_error"] < 1e-9
    assert info["voigt_gap_min_eigenvalue"] > -1e-8
    assert info["reuss_gap_min_eigenvalue"] > -1e-8
    assert np.all(np.isfinite(ae))
