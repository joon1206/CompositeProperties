"""Independent constrained-energy numerical reference for flat laminates.

Piecewise-constant normal/shear strains are the exact solution in flat layers.
Subdivision supplies an assembly/refinement check; it is not a 3D fiber RVE.
"""
import numpy as np

AA = [0,2,4]
BB = [1,3,5]


def constrained_response(layers, macro_strain, subdivide=1, temperature=0.0):
    if subdivide < 1:
        raise ValueError("subdivide must be positive")
    pieces = [(w/subdivide,c,a) for w,c,a in layers for _ in range(subdivide)]
    n = len(pieces)
    h = np.zeros((3*n,3*n)); rhs=np.zeros(3*n); weights=np.zeros((3,3*n))
    ea=np.asarray(macro_strain)[AA]
    eb=np.asarray(macro_strain)[BB]
    for i,(w,c,a) in enumerate(pieces):
        b=slice(3*i,3*i+3)
        h[b,b]=w*c[np.ix_(BB,BB)]
        rhs[b]=w*(temperature*(c@a)[BB]-c[np.ix_(BB,AA)]@ea)
        weights[:,b]=w*np.eye(3)
    kkt=np.block([[h,weights.T],[weights,np.zeros((3,3))]])
    sol=np.linalg.solve(kkt,np.r_[rhs,eb])[:3*n]
    mean_stress=np.zeros(6)
    for i,(w,c,a) in enumerate(pieces):
        strain=np.zeros(6);strain[AA]=ea;strain[BB]=sol[3*i:3*i+3]
        mean_stress+=w*(c@strain-temperature*c@a)
    return mean_stress


def numerical_effective(layers, subdivide=1):
    zero=np.zeros(6)
    columns=[constrained_response(layers,np.eye(6)[i],subdivide) for i in range(6)]
    c=np.column_stack(columns)
    beta=-constrained_response(layers,zero,subdivide,temperature=1)
    return c,np.linalg.solve(c,beta)
