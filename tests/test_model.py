from dataclasses import replace
import numpy as np
import pytest
from composite_research.model import Inputs, fractions, ply, rotated_cross_ply, exact_layers, voigt_layers, evaluate, filled_matrix, matrix_hs_bounds
from composite_research.reference import numerical_effective, constrained_response


def layers(p=Inputs()):
    c,a,_=ply(p);cr,ar=rotated_cross_ply(c,a)
    return [(0.5,c,a),(0.5,cr,ar)]


def test_phase_fractions_and_weight():
    f=fractions(Inputs())
    assert (f['fiber'],f['silica'],f['resin'])==pytest.approx((.42,.406,.174))
    assert sum(f[k] for k in ('fiber_weight','silica_weight','resin_weight'))==pytest.approx(1)


@pytest.mark.parametrize('p',[
    replace(Inputs(),vf=0),replace(Inputs(),vf=1),
    replace(Inputs(),phi=0),replace(Inputs(),phi=1),Inputs()])
def test_symmetry_spd_and_reference(p):
    lay=layers(p);c,a=exact_layers(lay)
    for n in (1,2,4):
        cn,an=numerical_effective(lay,n)
        np.testing.assert_allclose(cn,c,rtol=2e-14,atol=2e-14)
        np.testing.assert_allclose(an,a,rtol=2e-14,atol=2e-14)
    np.testing.assert_allclose(c,c.T,rtol=1e-14,atol=1e-14)
    assert np.linalg.eigvalsh(c).min()>0


def test_identical_layers_and_isotropic_limits():
    for p in (replace(Inputs(),vf=0,phi=0),replace(Inputs(),vf=1)):
        c,a,_=ply(p);cr,ar=rotated_cross_ply(c,a)
        ce,ae=exact_layers([(0.5,c,a),(0.5,cr,ar)])
        np.testing.assert_allclose(ce,c,rtol=1e-14,atol=1e-14)
        np.testing.assert_allclose(ae,a,rtol=1e-14,atol=1e-14)


def test_voigt_energy_bound_and_in_plane_equality():
    lay=layers();ce,ae=exact_layers(lay);cv,av=voigt_layers(lay)
    assert np.linalg.eigvalsh(cv-ce).min()>-1e-12
    # All admissible layers share tangential strains: their reduced plane-stress
    # response is equal even though full 3D uniaxial moduli need not be.
    np.testing.assert_allclose(ce[np.ix_([0,2,4],[0,2,4])]-ce[np.ix_([0,2,4],[1,3,5])]@np.linalg.solve(ce[np.ix_([1,3,5],[1,3,5])],ce[np.ix_([1,3,5],[0,2,4])]),
                               sum(w*(c[np.ix_([0,2,4],[0,2,4])]-c[np.ix_([0,2,4],[1,3,5])]@np.linalg.solve(c[np.ix_([1,3,5],[1,3,5])],c[np.ix_([1,3,5],[0,2,4])])) for w,c,a in lay))


def test_cross_ply_symmetry():
    r=evaluate(Inputs())
    assert r['EX']==pytest.approx(r['EZ'])
    assert r['GXY']==pytest.approx(r['GYZ'])
    assert r['AX']==pytest.approx(r['AZ'])
    legacy=evaluate(Inputs(),'legacy')
    assert legacy['NUYZ']!=pytest.approx(r['NUYZ'])
    assert legacy['NUYZ']/legacy['EZ']==pytest.approx(r['NUYZ']/r['EY'],rel=1e-4)


def test_thermal_free_expansion_and_reciprocity():
    c,a=exact_layers(layers())
    np.testing.assert_allclose(constrained_response(layers(),a,temperature=1),np.zeros(6),atol=1e-12)
    s=np.linalg.solve(c,np.eye(6))
    assert s[0,1]==pytest.approx(s[1,0])


def test_matrix_closure_limits_and_bounds():
    for phi in (0,.1,.4,.7,1):
        p=replace(Inputs(),phi=phi)
        em,nu,alpha=filled_matrix(p,'mori_tanaka')
        b=matrix_hs_bounds(p)
        assert em==pytest.approx(b['E_lower'])
        assert b['E_lower']<=b['E_upper']+1e-12
        assert filled_matrix(p,'legacy')[2]==pytest.approx(alpha)
    assert filled_matrix(replace(Inputs(),phi=0),'mori_tanaka')==pytest.approx((2.5,.37,50))
    assert filled_matrix(replace(Inputs(),phi=1),'mori_tanaka')==pytest.approx((70,.17,.5))
