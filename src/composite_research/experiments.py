"""Executable model audit, parameter sweep, uncertainty and convergence study."""
from __future__ import annotations

from dataclasses import asdict, replace
from pathlib import Path
import csv, hashlib, json, platform, time
import numpy as np
from scipy.stats import qmc, spearmanr
from scipy.optimize import brentq
from .model import Inputs, evaluate, ply, rotated_cross_ply, exact_layers, voigt_layers
from .reference import numerical_effective

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/'results'/'raw'; FIG=ROOT/'results'/'figures'
OUTPUTS=['EX','EY','EZ','GXY','GYZ','GXZ','NUXY','NUYZ','NUXZ','AX','AY','AZ']


def save_csv(path,rows):
    rows=list(rows)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)


def model_comparison():
    p=Inputs();rows=[]
    for method in ('legacy','voigt','exact'):
        r=evaluate(p,method)
        rows.append({'method':method,**{k:float(r[k]) for k in OUTPUTS}})
    r=evaluate(p,'exact','mori_tanaka')
    rows.append({'method':'exact_mori_tanaka_matrix',**{k:float(r[k]) for k in OUTPUTS}})
    save_csv(RAW/'baseline.csv',rows)
    return rows


def grid_study():
    rows=[]
    for vf in np.linspace(.2,.6,17):
        for phi in np.linspace(0,.7,29):
            p=replace(Inputs(),vf=float(vf),phi=float(phi))
            e=evaluate(p,'exact');v=evaluate(p,'voigt');l=evaluate(p,'legacy')
            mt=evaluate(p,'exact','mori_tanaka')
            rows.append({'vf':vf,'phi':phi,'silica_total':(1-vf)*phi,
                         'EX_exact':e['EX'],'EY_exact':e['EY'],'GXY_exact':e['GXY'],
                         'AX_exact':e['AX'],'AY_exact':e['AY'],
                         'EX_mt':mt['EX'],'EY_mt':mt['EY'],'AX_mt':mt['AX'],'AY_mt':mt['AY'],
                         'EY_voigt':v['EY'],'GXY_voigt':v['GXY'],
                         'AX_legacy':l['AX'],'AY_legacy':l['AY'],
                         'AY_voigt':v['AY'],
                         'EX_error_pct':100*(v['EX']/e['EX']-1),
                         'matrix_EX_error_pct':100*(e['EX']/mt['EX']-1),
                         'matrix_EY_error_pct':100*(e['EY']/mt['EY']-1),
                         'GXY_error_pct':100*(v['GXY']/e['GXY']-1),
                         'AX_error':l['AX']-e['AX'],'AY_error':l['AY']-e['AY']})
    save_csv(RAW/'grid.csv',rows)
    return rows


# Illustrative intervals, not empirically calibrated confidence intervals.
RANGES={'vf':(.32,.52),'phi':(.4,.7),'ef':(76.,86.),'af':(1.8,5.6),
        'er':(2.,4.),'ar':(40.,60.),'ep':(65.,75.),'ap':(.3,.7)}


def uncertainty(seed=20261024,n=512):
    rng=qmc.Sobol(d=len(RANGES),scramble=True,seed=seed)
    sample=rng.random_base2(int(np.log2(n)))
    keys=list(RANGES);rows=[]
    for i,x in enumerate(sample):
        kwargs={k:lo+(hi-lo)*x[j] for j,(k,(lo,hi)) in enumerate(RANGES.items())}
        p=replace(Inputs(),**kwargs)
        e=evaluate(p,'exact');v=evaluate(p,'voigt');l=evaluate(p,'legacy')
        mt=evaluate(p,'exact','mori_tanaka')
        rows.append({'sample':i,**kwargs,
                     'EX_exact':e['EX'],'EY_exact':e['EY'],'GXY_exact':e['GXY'],
                     'AX_exact':e['AX'],'AY_exact':e['AY'],
                     'EX_mt':mt['EX'],'EY_mt':mt['EY'],'AX_mt':mt['AX'],'AY_mt':mt['AY'],
                     'matrix_EX_bias_pct':100*(e['EX']/mt['EX']-1),
                     'matrix_EY_bias_pct':100*(e['EY']/mt['EY']-1),
                     'matrix_AY_bias':e['AY']-mt['AY'],
                     'EX_bias_pct':100*(v['EX']/e['EX']-1),
                     'GXY_bias_pct':100*(v['GXY']/e['GXY']-1),
                     'AX_bias':l['AX']-e['AX'],'AY_bias':l['AY']-e['AY']})
    save_csv(RAW/'uncertainty.csv',rows)
    summary={}
    for k in ('EX_exact','EY_exact','GXY_exact','AX_exact','AY_exact','matrix_EX_bias_pct','matrix_EY_bias_pct','matrix_AY_bias','EX_bias_pct','GXY_bias_pct','AX_bias','AY_bias'):
        z=np.array([row[k] for row in rows]);summary[k]={"median":float(np.median(z)),"p05":float(np.quantile(z,.05)),"p95":float(np.quantile(z,.95)),"min":float(z.min()),"max":float(z.max())}
    for outcome in ('EY_exact','AX_exact','AY_exact','matrix_EX_bias_pct','matrix_AY_bias','AY_bias'):
        summary[outcome]['spearman']={k:float(spearmanr([r[k] for r in rows],[r[outcome] for r in rows]).statistic) for k in keys}
    (RAW/'uncertainty_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return rows,summary


def derivatives():
    p=Inputs();rows=[]
    for name in RANGES:
        x=getattr(p,name);h=max(abs(x)*1e-4,1e-5)
        hi=evaluate(replace(p,**{name:x+h}),'exact')
        lo=evaluate(replace(p,**{name:x-h}),'exact')
        rows.append({'parameter':name,'step':h,**{k:(hi[k]-lo[k])/(2*h) for k in ('EX','EY','AX','AY')}})
    save_csv(RAW/'derivatives.csv',rows)


def convergence():
    p=Inputs();c,a,_=ply(p);cr,ar=rotated_cross_ply(c,a)
    lay=[(.5,c,a),(.5,cr,ar)]
    ce,ae=exact_layers(lay);rows=[]
    for n in (1,2,4,8,16):
        cn,an=numerical_effective(lay,n)
        rows.append({'subdivision_per_layer':n,'relative_C_error':np.linalg.norm(cn-ce)/np.linalg.norm(ce),
                     'absolute_alpha_error_ppmK':np.linalg.norm(an-ae),
                     'min_eigenvalue_GPa':np.linalg.eigvalsh(cn).min()})
    save_csv(RAW/'convergence.csv',rows)
    return rows


def sign_reversal():
    """Post-hoc falsification of the idea that legacy CTE has fixed bias."""
    rows=[]
    for vf in np.linspace(.2,.6,17):
        def bias(phi):
            p=replace(Inputs(),vf=float(vf),phi=float(phi))
            return evaluate(p,'legacy')['AY']-evaluate(p,'exact')['AY']
        roots=[];grid=np.linspace(0,.7,71)
        for lo,hi in zip(grid[:-1],grid[1:]):
            if bias(lo)*bias(hi)<0:roots.append(brentq(bias,lo,hi,xtol=1e-13))
        rows.append({'vf':vf,'AY_bias_phi0':bias(0),'AY_bias_phi07':bias(.7),
                     'zero_crossing_phi':roots[0] if len(roots)==1 else '',
                     'root_count':len(roots)})
    save_csv(RAW/'cte_bias_sign_reversal.csv',rows)
    return rows


def benchmark():
    p=Inputs();n=1000
    t=time.perf_counter()
    for _ in range(n):evaluate(p,'exact')
    dt=time.perf_counter()-t
    out={'iterations':n,'seconds':dt,'microseconds_per_evaluation':dt/n*1e6,
         'cpu':platform.processor(),'python':platform.python_version(),'numpy':np.__version__}
    (RAW/'benchmark.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out


def jax_benchmark(samples):
    """Measure x64 CPU JIT compile and steady-state batch time with synchronization."""
    try:
        import jax
        import jax.numpy as jnp
        from .jax_model import response,compiled_jacobian
    except ImportError:
        return {'status':'JAX unavailable'}
    keys=['vf','phi','ef','af','er','ar','ep','ap']
    batch=jnp.asarray([[r[k] for k in keys] for r in samples])
    batched=jax.jit(jax.vmap(response))
    t=time.perf_counter();result=batched(batch);result.block_until_ready();compile_seconds=time.perf_counter()-t
    times=[]
    for _ in range(5):
        t=time.perf_counter();batched(batch).block_until_ready();times.append(time.perf_counter()-t)
    t=time.perf_counter();jac=compiled_jacobian(batch[0]);jac.block_until_ready();jac_seconds=time.perf_counter()-t
    np_results=np.array([[r[k] for k in ('EX_exact','EY_exact','GXY_exact','AX_exact','AY_exact')] for r in samples])
    out={'status':'run','jax_version':jax.__version__,'backend':jax.default_backend(),
         'devices':[str(d) for d in jax.devices()],'batch_size':len(samples),
         'compile_and_first_seconds':compile_seconds,'steady_median_seconds':float(np.median(times)),
         'jacobian_first_seconds':jac_seconds,
         'max_absolute_vs_numpy':float(np.max(np.abs(np.asarray(result)-np_results)))}
    (RAW/'jax_benchmark.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
    return out


def figures(grid, samples):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'figure.dpi':160,'savefig.dpi':200})
    vs=sorted(set(r['vf'] for r in grid));ps=sorted(set(r['phi'] for r in grid))
    fig,axs=plt.subplots(1,3,figsize=(12,3.6),layout='constrained')
    for ax,k,title,unit in zip(axs,['matrix_EX_error_pct','GXY_error_pct','AY_error'],
                                ['In-plane matrix closure','Interlaminar shear','Thickness CTE'],
                                ['Legacy vs MT difference (%)','Voigt error (%)','Legacy minus exact (ppm/K)']):
        z=np.array([[next(r[k] for r in grid if r['vf']==v and r['phi']==p) for p in ps] for v in vs])
        im=ax.pcolormesh(ps,vs,z,shading='nearest',cmap='coolwarm')
        ax.set(xlabel='Silica fraction within matrix',ylabel='Fiber volume fraction',title=title)
        fig.colorbar(im,ax=ax,label=unit)
    fig.savefig(FIG/'model_error_map.png');fig.savefig(FIG/'model_error_map.pdf');plt.close(fig)
    fig,axs=plt.subplots(1,2,figsize=(8.5,3.5),layout='constrained')
    for ax,k,label in zip(axs,['matrix_EX_bias_pct','AY_bias'],['Legacy vs MT EX difference (%)','Legacy AY bias (ppm/K)']):
        ax.hist([r[k] for r in samples],bins=35,color='#276a85',edgecolor='white')
        ax.set(xlabel=label,ylabel='Illustrative parameter draws')
    fig.savefig(FIG/'uncertainty_bias.png');fig.savefig(FIG/'uncertainty_bias.pdf');plt.close(fig)


def main():
    RAW.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    b=model_comparison();g=grid_study();u,s=uncertainty();derivatives();c=convergence();sr=sign_reversal();perf=benchmark();jb=jax_benchmark(u);figures(g,u)
    source=Path(__file__).with_name('model.py')
    manifest={'model_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
              'inputs':asdict(Inputs()),'uncertainty_ranges':RANGES,'sobol_seed':20261024,
              'grid_points':len(g),'uncertainty_samples':len(u),'python':platform.python_version(),
              'numpy':np.__version__,'baseline':b,'convergence':c,'sign_reversal':sr,'benchmark':perf,'jax_benchmark':jb}
    (RAW/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps({'baseline':b,'uncertainty_summary':{k:{x:y for x,y in v.items() if x!='spearman'} for k,v in s.items()},'benchmark':perf,'jax_benchmark':jb},indent=2))


if __name__=='__main__':main()
