"""Small-strain thermoelastic micromechanics and exact flat-layer homogenization.

Engineering Voigt order: xx, yy, zz, yz, xz, xy. Engineering shear strain
is gamma_ij=2 epsilon_ij. Y is the laminate normal; X and Z are in-plane.
All stiffnesses are GPa, CTEs ppm/K, densities g/cm^3, and fractions
dimensionless.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np


@dataclass(frozen=True)
class Inputs:
    ef: float = 82.0
    nuf: float = 0.22
    af: float = 2.49
    rhof: float = 2.5
    vf: float = 0.42
    ep: float = 70.0
    nup: float = 0.17
    ap: float = 0.5
    rhop: float = 2.2
    phi: float = 0.70  # silica fraction OF NON-FIBER MATRIX
    er: float = 2.5
    nur: float = 0.37
    ar: float = 50.0
    rhor: float = 1.2

    def checked(self):
        if not (0 <= self.vf <= 1 and 0 <= self.phi <= 1):
            raise ValueError("volume fractions must be in [0,1]")
        if any(x <= 0 for x in (self.ef, self.ep, self.er, self.rhof, self.rhop, self.rhor)):
            raise ValueError("elastic moduli and densities must be positive")
        if any(not (-1 < x < .5) for x in (self.nuf, self.nup, self.nur)):
            raise ValueError("isotropic Poisson ratios must be in (-1,0.5)")
        return self


def fractions(p: Inputs):
    p.checked()
    f, s, r = p.vf, (1-p.vf)*p.phi, (1-p.vf)*(1-p.phi)
    mass = f*p.rhof + s*p.rhop + r*p.rhor
    return {"fiber": f, "silica": s, "resin": r,
            "fiber_weight": f*p.rhof/mass, "silica_weight": s*p.rhop/mass,
            "resin_weight": r*p.rhor/mass, "density": mass}


def halpin_tsai(em, ei, fraction, xi):
    eta = (ei/em-1)/(ei/em+xi)
    return em*(1+xi*eta*fraction)/(1-eta*fraction)


def filled_matrix(p: Inputs, method="legacy"):
    """Two alternative isotropic filled-resin closures.

    'legacy' preserves the script's empirical formulas. 'mori_tanaka'
    assumes bonded spherical silica in an isotropic resin and uses the exact
    two-phase isotropic bulk-modulus/CTE identity. Neither is experimental
    validation, especially near high filler packing.
    """
    p.checked()
    q, r = p.phi, 1-p.phi
    if method == "mori_tanaka":
        kr,gr=p.er/(3*(1-2*p.nur)),p.er/(2*(1+p.nur))
        kp,gp=p.ep/(3*(1-2*p.nup)),p.ep/(2*(1+p.nup))
        z=gr*(9*kr+8*gr)/(6*(kr+2*gr))
        km=kr+q*(kp-kr)/(1+r*(kp-kr)/(kr+4*gr/3))
        gm=gr+q*(gp-gr)/(1+r*(gp-gr)/(gr+z))
        em=9*km*gm/(3*km+gm)
        nu=(3*km-2*gm)/(2*(3*km+gm))
        alpha=p.ap+(p.ar-p.ap)*(1/km-1/kp)/(1/kr-1/kp)
        return em,nu,alpha
    if method != "legacy":
        raise ValueError(f"unknown filled-matrix method: {method}")
    em = halpin_tsai(p.er, p.ep, q, 2)
    nu = p.nup*q+p.nur*r
    numerator = 2*(p.ep*(1-2*p.nur)-p.er*(1-2*p.nup))
    denominator = 2*p.er*r*(1-2*p.nup)+p.ep*(2*q*(1-2*p.nur)+1+p.nur)
    alpha = p.ar*r+p.ap*q+(p.ap-p.ar)*q*r*numerator/denominator
    return em, nu, alpha


def matrix_hs_bounds(p: Inputs):
    """Hashin-Shtrikman isotropic K/G and implied E bounds (GPa).

    Assumes stiff particles (silica) in a softer resin; the bounds do not
    certify the empirical Halpin-Tsai E/linear-nu pair as a microstructure.
    """
    p.checked()
    kr,gr=p.er/(3*(1-2*p.nur)),p.er/(2*(1+p.nur))
    kp,gp=p.ep/(3*(1-2*p.nup)),p.ep/(2*(1+p.nup))
    if kp <= kr or gp <= gr:
        raise ValueError("HS implementation requires silica stiffer in K and G")
    q,r=p.phi,1-p.phi
    z=lambda k,g:g*(9*k+8*g)/(6*(k+2*g))
    kl=kr+q/(1/(kp-kr)+r/(kr+4*gr/3))
    ku=kp+r/(1/(kr-kp)+q/(kp+4*gp/3))
    gl=gr+q/(1/(gp-gr)+r/(gr+z(kr,gr)))
    gu=gp+r/(1/(gr-gp)+q/(gp+z(kp,gp)))
    young=lambda k,g:9*k*g/(3*k+g)
    return {"K_lower":kl,"K_upper":ku,"G_lower":gl,"G_upper":gu,
            "E_lower":young(kl,gl),"E_upper":young(ku,gu)}


def ply(p: Inputs, matrix_method="legacy"):
    """Reference script's transversely isotropic ply approximation."""
    em, num, am = filled_matrix(p,matrix_method)
    vm = 1-p.vf
    el = p.vf*p.ef+vm*em
    nult = p.vf*p.nuf+vm*num
    et = halpin_tsai(em, p.ef, p.vf, 2)
    gf, gm = p.ef/(2*(1+p.nuf)), em/(2*(1+num))
    glt = halpin_tsai(gm, gf, p.vf, 1)
    ntt = p.nuf*p.vf+num*vm*(1+num-nult*em/el)/(1-num*num+num*nult*em/el)
    gtt = et/(2*(1+ntt))
    a1 = (p.vf*p.ef*p.af+vm*em*am)/el
    a2 = p.af*p.vf*(1+p.nuf)+am*vm*(1+num)-a1*nult
    s = np.zeros((6,6))
    s[0,0],s[1,1],s[2,2] = 1/el,1/et,1/et
    s[0,1]=s[1,0]=s[0,2]=s[2,0]=-nult/el
    s[1,2]=s[2,1]=-ntt/et
    s[3,3],s[4,4],s[5,5] = 1/gtt,1/glt,1/glt
    if np.linalg.eigvalsh(s).min() <= 0:
        raise ValueError("ply compliance is not positive definite")
    c = np.linalg.solve(s,np.eye(6))
    return c, np.array([a1,a2,a2,0.,0.,0.]), {"Em":em,"num":num,"am":am,"EL":el,"ET":et,"GLT":glt,"GTT":gtt,"nuLT":nult,"nuTT":ntt,"a1":a1,"a2":a2}


def rotated_cross_ply(c, alpha):
    # 90 degrees about Y; sign of shear coupling is immaterial for this TI ply.
    idx = [2,1,0,5,4,3]
    return c[np.ix_(idx,idx)],alpha[idx]


def voigt_layers(layers):
    """Uniform-strain stiffness upper bound and CONSISTENT thermal stress mix."""
    c = sum(w*ck for w,ck,ak in layers)
    beta = sum(w*(ck@ak) for w,ck,ak in layers)
    return c,np.linalg.solve(c,beta)


def exact_layers(layers):
    """Exact infinite, flat, perfectly bonded layer homogenization along Y.

    Tangential strains (xx,zz,xz) and normal tractions (yy,yz,xy) are
    continuous. A Schur complement eliminates the layer normal strains.
    """
    aa = [0,2,4]
    bb = [1,3,5]
    kbar=np.zeros((3,3)); hbar=np.zeros((3,3))
    mbar=np.zeros((3,3)); nbar=np.zeros((3,3))
    qbar=np.zeros(3); rbar=np.zeros(3)
    if not np.isclose(sum(w for w,_,_ in layers),1):
        raise ValueError("layer fractions must sum to one")
    for w,c,a in layers:
        caa,cab,cba,cbb=c[np.ix_(aa,aa)],c[np.ix_(aa,bb)],c[np.ix_(bb,aa)],c[np.ix_(bb,bb)]
        beta=c@a
        invb=np.linalg.solve(cbb,np.eye(3))
        h=cab@invb
        kbar+=w*(caa-h@cba)
        hbar+=w*h
        mbar+=w*invb
        nbar+=w*(invb@cba)
        qbar+=w*(beta[aa]-h@beta[bb])
        rbar+=w*(invb@beta[bb])
    mb=np.linalg.solve(mbar,np.eye(3))
    c=np.zeros((6,6)); beta=np.zeros(6)
    c[np.ix_(aa,aa)]=kbar+hbar@mb@nbar
    c[np.ix_(aa,bb)]=hbar@mb
    c[np.ix_(bb,aa)]=mb@nbar
    c[np.ix_(bb,bb)]=mb
    beta[aa]=qbar+hbar@mb@rbar
    beta[bb]=mb@rbar
    return c,np.linalg.solve(c,beta)


def legacy(p: Inputs):
    """Reproduce only the meaningful algebra of the pasted MATLAB script.

    Its explicit CTE equations are retained, even though they do not equal
    the thermal eigenstrain implied by its own averaged stiffness.
    """
    c,a,d=ply(p)
    cr,ar=rotated_cross_ply(c,a)
    cv,_=voigt_layers([(0.5,c,a),(0.5,cr,ar)])
    ax=(d["a1"]*d["EL"]+d["a2"]*d["ET"])/(d["EL"]+d["ET"])
    ay=d["a2"]+(d["nuTT"]+d["nuLT"])*(ax-d["a1"])
    return cv,np.array([ax,ay,ax,0.,0.,0.])


def evaluate(p: Inputs, method="exact", matrix_method="legacy"):
    if method=="legacy" and matrix_method!="legacy":
        raise ValueError("legacy laminate requires legacy matrix model")
    c,a,d=ply(p,matrix_method)
    cr,ar=rotated_cross_ply(c,a)
    layers=[(0.5,c,a),(0.5,cr,ar)]
    ce,ae={"exact":exact_layers,"voigt":voigt_layers,"legacy":lambda _:legacy(p)}[method](layers)
    s=np.linalg.solve(ce,np.eye(6))
    return {"EX":1/s[0,0],"EY":1/s[1,1],"EZ":1/s[2,2],
            "GXY":1/s[5,5],"GYZ":1/s[3,3],"GXZ":1/s[4,4],
            "NUXY":-s[0,1]/s[0,0],
            # MATLAB labels -S_zy/S_zz as NUYZ, but this is nu_ZY.
            "NUYZ":(-s[2,1]/s[2,2] if method=="legacy" else -s[1,2]/s[1,1]),
            "NUXZ":-s[0,2]/s[0,0],
            "AX":ae[0],"AY":ae[1],"AZ":ae[2],"C":ce,"alpha":ae,
            "ply":d,"fractions":fractions(p)}
