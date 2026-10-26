import numpy as np
import pytest
from dataclasses import replace
from composite_research.model import Inputs,evaluate

jax=pytest.importorskip('jax')
from composite_research.jax_model import BASE,compiled_response,compiled_jacobian,compiled_study


def test_jax_matches_numpy_and_finite_difference():
    p=Inputs();r=evaluate(p)
    target=np.array([r[k] for k in ('EX','EY','GXY','AX','AY')])
    np.testing.assert_allclose(np.asarray(compiled_response(BASE)),target,rtol=1e-13)
    jac=np.asarray(compiled_jacobian(BASE))
    keys=['vf','phi','ef','af','er','ar','ep','ap']
    for i,k in enumerate(keys):
        x=float(BASE[i]);h=max(abs(x)*1e-5,1e-6)
        hi=evaluate(replace(p,**{k:x+h}))
        lo=evaluate(replace(p,**{k:x-h}))
        fd=np.array([(hi[z]-lo[z])/(2*h) for z in ('EX','EY','GXY','AX','AY')])
        np.testing.assert_allclose(jac[:,i],fd,rtol=2e-7,atol=2e-8)


def test_batched_study_matches_four_numpy_models():
    cases=[Inputs(),replace(Inputs(),vf=.25,phi=.1),
           replace(Inputs(),vf=.55,phi=.55,er=3.4,ar=44.)]
    keys=('vf','phi','ef','af','er','ar','ep','ap')
    batch=np.array([[getattr(p,k) for k in keys] for p in cases])
    actual=np.asarray(compiled_study(batch))
    for i,p in enumerate(cases):
        for j,(method,matrix) in enumerate((('exact','legacy'),('voigt','legacy'),
                                             ('legacy','legacy'),('exact','mori_tanaka'))):
            r=evaluate(p,method,matrix)
            target=np.array([r[k] for k in ('EX','EY','GXY','AX','AY')])
            np.testing.assert_allclose(actual[i,j],target,rtol=3e-13,atol=3e-13)
