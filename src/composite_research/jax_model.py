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


def _ply(x):
    vf,phi,ef,af,er,ar,ep,ap=x
    nf,np_,nr=.22,.17,.37
    vm=1-vf;fr=1-phi
    em=_ht(er,ep,phi,2)
    nm=np_*phi+nr*fr
    num=2*(ep*(1-2*nr)-er*(1-2*np_))
    den=2*er*fr*(1-2*np_)+ep*(2*phi*(1-2*nr)+1+nr)
    am=ar*fr+ap*phi+(ap-ar)*phi*fr*num/den
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
    return j.linalg.solve(s,j.eye(6)),j.array([a1,a2,a2,0,0,0])


def _reorder(c,a):
    idx=j.array([2,1,0,5,4,3])
    return c[idx[:,None],idx[None,:]],a[idx]


def response(x):
    c0,a0=_ply(x);c1,a1=_reorder(c0,a0)
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
    s=j.linalg.solve(ce,j.eye(6))
    return j.array([1/s[0,0],1/s[1,1],1/s[5,5],alpha[0],alpha[1]])


compiled_response=jax.jit(response)
compiled_jacobian=jax.jit(jax.jacfwd(response))
