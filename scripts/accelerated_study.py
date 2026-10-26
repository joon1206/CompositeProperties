"""Reproducible four-model, five-output CPU/GPU study and timing protocol.

Run inside WSL with LD_LIBRARY_PATH unset to avoid OpenFOAM's old libstdc++.
Set JAX_PLATFORMS=cpu or cuda explicitly; a failed CUDA setup must not fall back.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
from pathlib import Path
import platform
import time

import jax
import numpy as np
from scipy.stats import qmc

from composite_research.experiments import RANGES
from composite_research.jax_model import compiled_study
from composite_research.model import Inputs, evaluate

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'results'/'raw'
KEYS=tuple(RANGES)
OUT=('EX','EY','GXY','AX','AY')
MODELS=(('exact','legacy'),('voigt','legacy'),('legacy','legacy'),
        ('exact','mori_tanaka'))


def samples(n=8192):
    assert n > 0 and n & (n-1) == 0
    u=qmc.Sobol(d=len(KEYS),scramble=True,seed=20261024).random_base2(n.bit_length()-1)
    lo=np.array([RANGES[k][0] for k in KEYS]);hi=np.array([RANGES[k][1] for k in KEYS])
    return lo+(hi-lo)*u


def numpy_study(x):
    out=[]
    for row in x:
        p=replace(Inputs(),**dict(zip(KEYS,map(float,row))))
        values=[]
        for m,mm in MODELS:
            r=evaluate(p,m,mm)
            values.append([r[k] for k in OUT])
        out.append(values)
    return np.asarray(out)


def stats(z):
    return {k:float(v) for k,v in zip(('min','p05','median','p95','max'),
            (np.min(z),np.quantile(z,.05),np.median(z),np.quantile(z,.95),np.max(z)))}


def cpu_name():
    info=Path('/proc/cpuinfo')
    if info.exists():
        for line in info.read_text().splitlines():
            if line.startswith('model name'):
                return line.split(':',1)[1].strip()
    return platform.processor()


def bench(fun, reps=9):
    times=[]
    for _ in range(reps):
        t=time.perf_counter();fun();times.append(time.perf_counter()-t)
    return {'median_seconds':float(np.median(times)),
            'min_seconds':float(min(times)),'max_seconds':float(max(times)),
            'repetitions':reps}


def validate_saved_sweep(filename,device):
    with (RAW/filename).open(newline='') as f:
        rows=list(csv.DictReader(f))
    default=Inputs()
    x=np.array([[float(r[k]) if k in r else getattr(default,k) for k in KEYS]
                for r in rows])
    actual=np.asarray(compiled_study(jax.device_put(x,device)))
    names={0:'exact',1:'voigt',2:'legacy',3:'mt'}
    errors=[]
    for i,r in enumerate(rows):
        for model,label in names.items():
            for output,kind in enumerate(OUT):
                key=f'{kind}_{label}'
                if key in r:
                    errors.append(abs(actual[i,model,output]-float(r[key])))
    return {'n':len(rows),'compared_values':len(errors),
            'max_absolute_difference':float(max(errors))}


def run(backend, extended=False, analysis_only=False):
    if jax.default_backend()!=backend:
        raise RuntimeError(f'requested {backend}, got {jax.default_backend()}')
    RAW.mkdir(parents=True,exist_ok=True)
    x=samples()
    device=jax.devices(backend)[0]
    records=[]
    for n in (() if analysis_only else (64,512,2048,8192)):
        host=x[:n]
        batch=jax.device_put(host,device)
        t=time.perf_counter();first=compiled_study(batch);first.block_until_ready()
        first_seconds=time.perf_counter()-t
        # Warm device/driver before measuring steady state.
        for _ in range(3):compiled_study(batch).block_until_ready()
        resident=bench(lambda:compiled_study(batch).block_until_ready())
        # Includes host input transfer and host output copy, but excludes JIT compilation.
        roundtrip=bench(lambda:np.asarray(compiled_study(jax.device_put(host,device))),5)
        row={'n':n,'compile_and_first_seconds':first_seconds,
             'device_resident':resident,'host_roundtrip':roundtrip}
        if n <= 2048:
            ref=numpy_study(host)
            row['max_abs_vs_numpy']=float(np.max(np.abs(np.asarray(first)-ref)))
            row['numpy_scalar']=bench(lambda:numpy_study(host),3)
        records.append(row)
        print(f'{backend} {n}: resident {resident["median_seconds"]:.6f}s, '
              f'roundtrip {roundtrip["median_seconds"]:.6f}s',flush=True)
    result={'backend':backend,'jax_version':jax.__version__,
            'device':str(device),'cpu':cpu_name(),
            'python':platform.python_version(),'outputs':list(OUT),
            'models':[f'{m}/{mm}' for m,mm in MODELS],
            'timing_protocol':'median of synchronized repetitions; first call separately',
            'batch_seed':20261024,'sizes':records}
    if not analysis_only:
        (RAW/f'performance_{backend}.json').write_text(json.dumps(result,indent=2))
    if extended:
        # This is a larger illustrative scenario sweep, not a measured confidence interval.
        y=np.asarray(compiled_study(jax.device_put(x,device)))
        exact,voigt,legacy,mt=np.moveaxis(y,1,0)
        ex_gap=100*(exact[:,0]/mt[:,0]-1)
        ey_gap=100*(exact[:,1]/mt[:,1]-1)
        shear_gap=100*(voigt[:,2]/exact[:,2]-1)
        ay_formula=legacy[:,4]-exact[:,4]
        ay_matrix=exact[:,4]-mt[:,4]
        ay_net=legacy[:,4]-mt[:,4]
        hidden=(np.abs(ay_net)<=.2)&(np.abs(ay_formula)>=.5)&(np.abs(ay_matrix)>=.5)
        summary={'n':len(x),'seed':20261024,'ranges':RANGES,
                 'backend':backend,'EX_closure_gap_pct':stats(ex_gap),
                 'EY_closure_gap_pct':stats(ey_gap),
                 'GXY_layer_gap_pct':stats(shear_gap),
                 'AY_legacy_formula_gap_ppmK':stats(ay_formula),
                 'AY_closure_gap_ppmK':stats(ay_matrix),
                 'AY_legacy_vs_MT_net_gap_ppmK':stats(ay_net),
                 'fraction_net_AY_gap_within_0p2_ppmK':float(np.mean(np.abs(ay_net)<=.2)),
                 'fraction_hidden_opposing_gaps':float(np.mean(hidden)),
                 'fraction_EX_closure_gap_gt_5pct':float(np.mean(ex_gap>5)),
                 'fraction_EY_closure_gap_gt_10pct':float(np.mean(ey_gap>10)),
                 'fraction_AY_formula_gap_positive':float(np.mean(ay_formula>0)),
                 'fraction_AY_closure_gap_negative':float(np.mean(ay_matrix<0)),
                 'max_abs_vs_numpy_32':float(np.max(np.abs(y[:32]-numpy_study(x[:32])))),
                 'saved_grid_check':validate_saved_sweep('grid.csv',device),
                 'saved_uncertainty_check':validate_saved_sweep('uncertainty.csv',device)}
        (RAW/f'extended_study_{backend}.json').write_text(json.dumps(summary,indent=2))
        print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--backend',choices=('gpu','cpu'),required=True)
    ap.add_argument('--extended',action='store_true')
    ap.add_argument('--analysis-only',action='store_true')
    args=ap.parse_args()
    run(args.backend,args.extended,args.analysis_only)
