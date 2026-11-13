"""Compare repository filled-resin estimates with isotropic two-phase bounds."""
from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from composite_research.model import (Inputs, filled_matrix, matrix_hs_bounds,
                                      matrix_cte_hs_bounds)

ROOT = Path(__file__).resolve().parents[1]


def main():
    raw, figdir = ROOT / "results" / "raw", ROOT / "results" / "figures"
    raw.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for phi in np.linspace(0, .75, 151):
        p = replace(Inputs(), phi=float(phi))
        legacy_e, legacy_nu, legacy_a = filled_matrix(p, "legacy")
        mt_e, mt_nu, mt_a = filled_matrix(p, "mori_tanaka")
        elastic = matrix_hs_bounds(p)
        thermal = matrix_cte_hs_bounds(p)
        rows.append({"phi":float(phi), "legacy_E_GPa":legacy_e,
                     "mori_tanaka_E_GPa":mt_e, "HS_E_lower_GPa":elastic["E_lower"],
                     "HS_E_upper_GPa":elastic["E_upper"],
                     "legacy_nu":legacy_nu, "mori_tanaka_nu":mt_nu,
                     "legacy_CTE_ppmK":legacy_a, "mori_tanaka_CTE_ppmK":mt_a,
                     "HS_CTE_lower_ppmK":thermal["CTE_lower"],
                     "HS_CTE_upper_ppmK":thermal["CTE_upper"]})
    with (raw / "matrix_model_bounds.csv").open("w",newline="",encoding="utf-8") as f:
        writer = csv.DictWriter(f,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    x = [r["phi"] for r in rows]
    fig, axes = plt.subplots(1,2,figsize=(10.3,4.9))
    ax=axes[0]
    ax.fill_between(x,[r["HS_E_lower_GPa"] for r in rows],
                    [r["HS_E_upper_GPa"] for r in rows],color="#d9e7ed",
                    label="Hashin–Shtrikman admissible interval")
    ax.plot(x,[r["mori_tanaka_E_GPa"] for r in rows],color="#126782",lw=2,
            label="Mori–Tanaka = lower bound")
    ax.plot(x,[r["legacy_E_GPa"] for r in rows],color="#b65f2d",lw=2,
            label="Supplied Halpin–Tsai E")
    ax.set(ylabel="Filled-resin Young's modulus (GPa)",xlabel="Silica fraction of nonfiber matrix")
    ax.legend(frameon=False,fontsize=8)
    ax=axes[1]
    ax.fill_between(x,[r["HS_CTE_lower_ppmK"] for r in rows],
                    [r["HS_CTE_upper_ppmK"] for r in rows],color="#e9e1cf",
                    label="CTE interval from bulk-modulus bounds")
    ax.plot(x,[r["mori_tanaka_CTE_ppmK"] for r in rows],color="#126782",lw=2,
            label="Mori–Tanaka = supplied formula")
    ax.set(ylabel="Filled-resin CTE (ppm/K)",xlabel="Silica fraction of nonfiber matrix")
    ax.legend(frameon=False,fontsize=8)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle("Filled-resin model choice affects stiffness and thermal expansion")
    fig.text(.5,.015,"Two-phase isotropic bounds are not a calibrated particle morphology; illustrative constituent inputs.",
             ha="center",fontsize=8,color="#555")
    fig.tight_layout(rect=(0,.04,1,.94))
    fig.savefig(figdir / "matrix_model_bounds.png",dpi=180)
    fig.savefig(figdir / "matrix_model_bounds.pdf")
    plt.close(fig)
    for phi in (.4,.7,.75):
        r=next(row for row in rows if np.isclose(row["phi"],phi))
        print(phi,"E MT/legacy/HS upper",r["mori_tanaka_E_GPa"],
              r["legacy_E_GPa"],r["HS_E_upper_GPa"],
              "CTE HS lower/MT",r["HS_CTE_lower_ppmK"],r["mori_tanaka_CTE_ppmK"])


if __name__=="__main__":
    main()
