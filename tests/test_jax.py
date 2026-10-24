import numpy as np
import pytest
from dataclasses import replace
from composite_research.model import Inputs,evaluate

jax=pytest.importorskip('jax')
from composite_research.jax_model import BASE,compiled_response,compiled_jacobian


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
