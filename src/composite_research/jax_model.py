"""JAX x64 mirror of the reference-script closure and exact layer solution.

Vector x = [fiber Vf, nested silica phi, Ef, alpha_f, Er, alpha_r, Ep,
alpha_p]; other constituent properties are fixed at the illustrative baseline.
"""
import jax
jax.config.update("jax_enable_x64",True)
import jax.numpy as j

BASE=j.array([.42,.70,82.,2.49,2.5,50.,70.,.5])


def _ht(em,ei,f,xi):
    eta=(ei/em-1)/(ei/em+xi)
    return em*(1+xi*eta*f)/(1-eta*f)


def _ply_detail(x, matrix_method="legacy"):
    vf,phi,ef,af,er,ar,ep,ap=x
    nf,np_,nr=.22,.17,.37
    vm=1-vf;fr=1-phi
    if matrix_method == "legacy":
        em=_ht(er,ep,phi,2)
        nm=np_*phi+nr*fr
        num=2*(ep*(1-2*nr)-er*(1-2*np_))
        den=2*er*fr*(1-2*np_)+ep*(2*phi*(1-2*nr)+1+nr)
        am=ar*fr+ap*phi+(ap-ar)*phi*fr*num/den
    elif matrix_method == "mori_tanaka":
        kr,gr=er/(3*(1-2*nr)),er/(2*(1+nr))
        kp,gp=ep/(3*(1-2*np_)),ep/(2*(1+np_))
        z=gr*(9*kr+8*gr)/(6*(kr+2*gr))
        km=kr+phi*(kp-kr)/(1+fr*(kp-kr)/(kr+4*gr/3))
        gm=gr+phi*(gp-gr)/(1+fr*(gp-gr)/(gr+z))
        em=9*km*gm/(3*km+gm)
        nm=(3*km-2*gm)/(2*(3*km+gm))
        am=ap+(ar-ap)*(1/km-1/kp)/(1/kr-1/kp)
    else:
        raise ValueError(matrix_method)
    el=vf*ef+vm*em
    nlt=vf*nf+vm*nm
    et=_ht(em,ef,vf,2)
    gf=ef/(2*(1+nf));gm=em/(2*(1+nm))
    glt=_ht(gm,gf,vf,1)
    ntt=nf*vf+nm*vm*(1+nm-nlt*em/el)/(1-nm*nm+nm*nlt*em/el)
    gtt=et/(2*(1+ntt))
    a1=(vf*ef*af+vm*em*am)/el
    a2=af*vf*(1+nf)+am*vm*(1+nm)-a1*nlt
    s=j.array([[1/el,-nlt/el,-nlt/el,0,0,0],
               [-nlt/el,1/et,-ntt/et,0,0,0],
               [-nlt/el,-ntt/et,1/et,0,0,0],
               [0,0,0,1/gtt,0,0],
               [0,0,0,0,1/glt,0],
               [0,0,0,0,0,1/glt]])
    return j.linalg.solve(s,j.eye(6)),j.array([a1,a2,a2,0,0,0]),(el,et,nlt,ntt,a1,a2)


def _ply(x):
    c,a,_=_ply_detail(x)
    return c,a


def _reorder(c,a):
    idx=j.array([2,1,0,5,4,3])
    return c[idx[:,None],idx[None,:]],a[idx]


def _effective(c0,a0,method="exact"):
    c1,a1=_reorder(c0,a0)
    if method == "voigt":
        ce=(c0+c1)/2
        alpha=j.linalg.solve(ce,(c0@a0+c1@a1)/2)
        return ce,alpha
    if method != "exact":
        raise ValueError(method)
    aa=j.array([0,2,4]);bb=j.array([1,3,5])
    def blocks(c,a):
        caa=c[aa[:,None],aa[None,:]]
        cab=c[aa[:,None],bb[None,:]]
        cba=c[bb[:,None],aa[None,:]]
        cbb=c[bb[:,None],bb[None,:]]
        ib=j.linalg.solve(cbb,j.eye(3));h=cab@ib;beta=c@a
        return (caa-h@cba,h,ib,ib@cba,beta[aa]-h@beta[bb],ib@beta[bb])
    b0=blocks(c0,a0);b1=blocks(c1,a1)
    k,h,m,n,q,r=tuple((u+v)/2 for u,v in zip(b0,b1))
    mb=j.linalg.solve(m,j.eye(3))
    caa=k+h@mb@n;cab=h@mb;cba=mb@n;cbb=mb
    beta_a=q+h@mb@r;beta_b=mb@r
    # Permutation into [aa,bb] order avoids functional indexed assignments.
    cb=j.concatenate([j.concatenate([caa,cab],axis=1),j.concatenate([cba,cbb],axis=1)],axis=0)
    betab=j.concatenate([beta_a,beta_b])
    order=j.array([0,3,1,4,2,5])
    ce=cb[order[:,None],order[None,:]]
    beta=betab[order]
    alpha=j.linalg.solve(ce,beta)
    return ce,alpha


def _summary(ce,alpha):
    s=j.linalg.solve(ce,j.eye(6))
    return j.array([1/s[0,0],1/s[1,1],1/s[5,5],alpha[0],alpha[1]])


def response(x):
    c,a=_ply(x)
    return _summary(*_effective(c,a))


def study_response(x):
    """[exact, consistent Voigt, supplied-script, Mori-Tanaka] x 5 outputs."""
    c,a,(el,et,nlt,ntt,a1,a2)=_ply_detail(x)
    ce,ae=_effective(c,a)
    cv,av=_effective(c,a,"voigt")
    ax=(a1*el+a2*et)/(el+et)
    ay=a2+(ntt+nlt)*(ax-a1)
    legacy=_summary(cv,j.array([ax,ay,ax,0.,0.,0.]))
    cm,am,_=_ply_detail(x,"mori_tanaka")
    return j.stack((_summary(ce,ae),_summary(cv,av),legacy,
                    _summary(*_effective(cm,am))))


compiled_response=jax.jit(response)
compiled_jacobian=jax.jit(jax.jacfwd(response))
compiled_study=jax.jit(jax.vmap(study_response))
