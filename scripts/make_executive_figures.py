"""Build presentation-ready charts from measured benchmark and study JSON."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'results'/'raw'
FIG=ROOT/'results'/'figures'
FIG.mkdir(parents=True,exist_ok=True)

NAVY='#132943';TEAL='#087E8B';ORANGE='#D36B38';GRAY='#586779';LIGHT='#EDF2F5'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,
                     'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':220})


def label_box(ax,x,y,w,h,title,value,detail,color=TEAL):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02',
                                fc=LIGHT,ec='none',transform=ax.transAxes))
    ax.text(x+.035,y+h-.18,title,color=GRAY,fontsize=11,
            transform=ax.transAxes,va='top')
    ax.text(x+.035,y+.39,value,color=color,fontsize=27,fontweight='bold',
            transform=ax.transAxes,va='center')
    ax.text(x+.035,y+.11,detail,color=NAVY,fontsize=10.5,
            transform=ax.transAxes,va='bottom')


def performance():
    cpu=json.loads((RAW/'performance_cpu.json').read_text())
    gpu=json.loads((RAW/'performance_gpu.json').read_text())
    ns=[r['n'] for r in cpu['sizes']]
    tc=np.array([r['host_roundtrip']['median_seconds']*1000 for r in cpu['sizes']])
    tg=np.array([r['host_roundtrip']['median_seconds']*1000 for r in gpu['sizes']])
    speed=tc/tg
    fig=plt.figure(figsize=(13.33,7.5),facecolor='white')
    fig.text(.055,.95,'GPU value appears at larger study sizes',fontsize=22,
             fontweight='bold',color=NAVY,va='top')
    fig.text(.055,.89,'Same four models and five outputs per scenario • JAX on CPU vs RTX 2080 Ti',
             fontsize=12.5,color=GRAY)
    cards=fig.add_axes([.055,.64,.89,.20]);cards.axis('off')
    label_box(cards,.00,.02,.30,.91,'8,192 scenarios',f'{speed[-1]:.1f}× faster',
              'GPU vs optimized CPU, warm run')
    label_box(cards,.34,.02,.30,.91,'GPU batch time',f'{tg[-1]:.1f} ms',
              f'CPU: {tc[-1]:.1f} ms; data round trip included')
    label_box(cards,.68,.02,.30,.91,'Measured large batches','≥2,048',
              'GPU wins at measured larger batches',ORANGE)
    ax=fig.add_axes([.12,.23,.79,.34])
    colors=[GRAY if s<1 else TEAL for s in speed]
    bars=ax.barh(np.arange(len(ns)),speed,color=colors,height=.54)
    ax.axvline(1,color=ORANGE,linestyle='--',linewidth=1.5)
    ax.set_yticks(np.arange(len(ns)),[f'{n:,} scenarios' for n in ns])
    ax.invert_yaxis();ax.set_xlim(0,max(speed)*1.22)
    ax.set_xlabel('Speedup = CPU time ÷ GPU time  (above 1× favors GPU)',color=NAVY)
    ax.set_title('Measured speedup by batch size',loc='left',fontweight='bold',color=NAVY,pad=12)
    ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    for bar,s,c,g in zip(bars,speed,tc,tg):
        ax.text(bar.get_width()+.08,bar.get_y()+bar.get_height()/2,
                f'{s:.2f}×   ({c:.1f} vs {g:.1f} ms)',va='center',fontsize=10.5,color=NAVY)
    fig.text(.055,.145,'TIMING SCOPE',fontsize=9.5,fontweight='bold',color=TEAL)
    fig.text(.055,.115,'Median synchronized runs; input transfer and output copy included. Compilation excluded from bars.',
             fontsize=10,color=GRAY)
    fig.text(.055,.087,'First call took 0.6–0.7 s on CPU and 0.7–1.0 s on GPU across measured sizes.',
             fontsize=10,color=GRAY)
    fig.text(.055,.058,'Machine: Ryzen 5 5600X, RTX 2080 Ti, WSL2; exact timings vary with load and hardware.',
             fontsize=10,color=GRAY)
    for ext in ('png','pdf'):fig.savefig(FIG/f'executive_gpu_speedup.{ext}',bbox_inches='tight')
    plt.close(fig)


def science():
    summary=json.loads((RAW/'extended_study_gpu.json').read_text())
    with (RAW/'baseline.csv').open(newline='') as f:
        base={r['method']:r for r in csv.DictReader(f)}
    ay_legacy=float(base['legacy']['AY'])
    ay_exact=float(base['exact']['AY'])
    ay_mt=float(base['exact_mori_tanaka_matrix']['AY'])
    vals=np.array([ay_legacy-ay_exact,ay_exact-ay_mt,ay_legacy-ay_mt])
    fig=plt.figure(figsize=(13.33,7.5),facecolor='white')
    fig.text(.055,.95,'Model errors can cancel and hide disagreement',fontsize=22,
             fontweight='bold',color=NAVY,va='top')
    fig.text(.055,.89,'Useful review priority: validate the filled-matrix model before refining layer averaging',
             fontsize=12.5,color=GRAY)
    ax=fig.add_axes([.24,.49,.58,.29])
    names=['Supplied CTE formula vs exact layers','Matrix closure: exact vs Mori–Tanaka',
           'Net: supplied model vs Mori–Tanaka']
    y=np.arange(3)
    ax.barh(y,vals,color=[ORANGE,TEAL,NAVY],height=.5)
    ax.axvline(0,color=NAVY,lw=1)
    ax.set_yticks(y,names);ax.invert_yaxis();ax.set_xlim(-.83,.88)
    ax.set_xlabel('Difference in thickness CTE (ppm/K)',color=NAVY)
    ax.set_title('Illustrative baseline: two sizable errors leave a small net gap',
                 loc='left',fontweight='bold',color=NAVY,pad=13)
    ax.grid(axis='x',alpha=.15);ax.set_axisbelow(True)
    for yy,v in zip(y,vals):
        ax.text(v+(.025 if v>=0 else -.025),yy,f'{v:+.3f}',
                ha='left' if v>=0 else 'right',va='center',fontweight='bold',color=NAVY)
    frac=summary['fraction_hidden_opposing_gaps']*100
    fig.text(.055,.35,f'{frac:.1f}%',fontsize=37,color=ORANGE,fontweight='bold')
    fig.text(.28,.363,'of 8,192 illustrative cases had |net gap| ≤0.2 ppm/K',fontsize=13,color=NAVY)
    fig.text(.28,.328,'while each opposing component was ≥0.5 ppm/K in magnitude.',fontsize=12,color=GRAY)
    fig.text(.055,.235,'MODEL CHOICE IMPACT',fontsize=10,fontweight='bold',color=TEAL)
    fig.text(.055,.202,f'Thickness modulus: {summary["EY_closure_gap_pct"]["median"]:.1f}% median matrix-closure gap; '
             f'in-plane modulus: {summary["EX_closure_gap_pct"]["median"]:.1f}%.',fontsize=12,color=NAVY)
    fig.text(.055,.166,f'Layer-averaging gap in interlaminar shear: {summary["GXY_layer_gap_pct"]["median"]:.2f}% median.',
             fontsize=12,color=NAVY)
    fig.text(.055,.09,'Scenario ranges are illustrative, not calibrated probabilities or measured material behavior.',
             fontsize=10,color=GRAY)
    fig.text(.055,.06,'The 0.2 and 0.5 ppm/K thresholds are exploratory screens; model agreement is not experimental validation.',
             fontsize=10,color=GRAY)
    for ext in ('png','pdf'):fig.savefig(FIG/f'executive_model_risk.{ext}',bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':
    performance();science()
