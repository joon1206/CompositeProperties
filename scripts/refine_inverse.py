"""Constrained local refinement of the inverse-screen low-CTE boundary."""
from __future__ import annotations
import csv
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from composite_research.model import Inputs
from composite_research.woven import woven_orientation
from inverse_weave_study import BOUNDS, FIBERS, KEYS, RAW, write_csv


def main():
    with (RAW / "inverse_all_samples.csv").open(newline="", encoding="utf-8") as f:
        samples = list(csv.DictReader(f))
    names = list(BOUNDS)
    bounds = list(BOUNDS.values())
    out = []
    for fiber, (ef, nuf, af) in FIBERS.items():
        for method in ("legacy", "mori_tanaka"):
            group = [r for r in samples if r["fiber"]==fiber and r["matrix_method"]==method]
            # Include feasible starts and near-boundary starts; local optimization
            # does not certify global infeasibility.
            group.sort(key=lambda r: max(float(r["AX"]),float(r["AZ"]))
                + 5*sum(max(31-float(r[k]), 0, float(r[k])-33)**2 for k in ("EX","EZ")))
            starts = group[:5]
            for start in starts:
                x0 = np.array([float(start[k]) for k in names])
                def calc(x):
                    d = dict(zip(names, x))
                    p = replace(Inputs(), ef=ef, nuf=nuf, af=af,
                                **{k:d[k] for k in ("vf","phi","er","ar","ep","ap")})
                    return woven_orientation(p, matrix_method=method,
                                             warp_crimp_deg=d["warp_crimp_deg"],
                                             weft_crimp_deg=d["weft_crimp_deg"],
                                             warp_share=d["warp_share"])
                def objective(x):
                    r=calc(x)
                    return (r["AX"]+r["AZ"])/2 + .02*(r["AX"]-r["AZ"])**2
                constraints = [
                    {"type":"ineq", "fun": lambda x,k=k: calc(x)[k]-31} for k in ("EX","EZ")
                ] + [
                    {"type":"ineq", "fun": lambda x,k=k: 33-calc(x)[k]} for k in ("EX","EZ")
                ]
                fit = minimize(objective, x0, method="SLSQP", bounds=bounds,
                               constraints=constraints, options={"maxiter":120,"ftol":1e-8})
                r=calc(fit.x)
                feasible = all(31-1e-5 <= r[k] <= 33+1e-5 for k in ("EX","EZ"))
                out.append({"fiber":fiber,"matrix_method":method,"start_sample_id":start["sample_id"],
                            "optimizer_success":fit.success,"feasible_modulus_XZ":feasible,
                            "objective":float(fit.fun),"iterations":fit.nit,
                            **dict(zip(names, map(float,fit.x))),
                            **{k:float(r[k]) for k in KEYS}})
            print("refined",fiber,method,flush=True)
    write_csv(RAW / "inverse_local_refinements.csv",out)
    feasible = [r for r in out if r["feasible_modulus_XZ"]]
    best = min(feasible,key=lambda r:max(r["AX"],r["AZ"]))
    (RAW / "inverse_refined_best.json").write_text(json.dumps(best,indent=2),encoding="utf-8")
    print(json.dumps(best,indent=2))


if __name__=="__main__":
    main()
