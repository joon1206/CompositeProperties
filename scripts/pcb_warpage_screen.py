"""Regenerate no-ring, bare-PCB thermal laminate screening figures and CSVs."""
from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from composite_research.pcb_laminate import PCBScreenInputs, core_q_alpha, screen_pcb


ROOT = Path(__file__).resolve().parents[1]
RAW, FIG = ROOT/"results"/"raw", ROOT/"results"/"figures"


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    q, baseline_alpha = core_q_alpha()
    baseline_cte = float(baseline_alpha[0])
    scenarios = (("Illustrative lower-CTE core: 2.0", 2.0, "#6a51a3"),
                 (f"Repository baseline core: {baseline_cte:.2f}", baseline_cte, "#31a354"))
    # Equal top and bottom coverage gives zero thermal bow in this model.
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.9), sharey=True)
    rows = []
    for ax, thickness in zip(axes, (1.6, 3.0)):
        imbalance = np.linspace(0, .5, 26)
        for label, core_cte, color in scenarios:
            values = []
            for difference in imbalance:
                p = replace(PCBScreenInputs(), core_thickness_mm=thickness,
                            top_copper_coverage=.75+difference/2,
                            bottom_copper_coverage=.75-difference/2)
                result = screen_pcb(core_cte, p, q)
                bow = result["bow_mm"]*1e3
                values.append(bow)
                rows.append({"side_mm": p.side_mm, "core_thickness_mm": thickness,
                             "core_CTE_ppm_K": core_cte,
                             "top_copper_coverage": p.top_copper_coverage,
                             "bottom_copper_coverage": p.bottom_copper_coverage,
                             "coverage_difference": difference, "delta_temp_K": p.delta_temp_k,
                             "dielectric_E_GPa": p.dielectric_e_gpa,
                             "dielectric_CTE_ppm_K": p.dielectric_cte_ppm_k,
                             "bow_um": bow,
                             "board_midplane_CTE_ppm_K":
                             result["midplane_strain"][0]/p.delta_temp_k*1e6})
            ax.plot(imbalance, values, color=color, lw=2.5, label=label)
            marker = screen_pcb(core_cte, replace(PCBScreenInputs(),
                                core_thickness_mm=thickness,
                                top_copper_coverage=.90,
                                bottom_copper_coverage=.60), q)["bow_mm"]*1e3
            ax.scatter([.30], [marker], color=color, s=38, zorder=3)
        ax.set_title(f"Core thickness {thickness:.1f} mm")
        ax.set_xlabel("Top minus bottom copper coverage")
        ax.set_xlim(0, .5)
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.2)
    axes[0].set_ylabel("Peak-to-valley thermal bow (µm per 100 K)")
    axes[0].legend(loc="upper left", frameon=False, fontsize=8)
    fig.suptitle("Bare PCB: copper imbalance creates thermal bow", fontsize=14, weight="bold")
    fig.text(.5, .015, "100 mm free panel, 10 build-up pairs per side. E and CTE of dielectric/copper are illustrative.\n"
             "Zero imbalance gives zero bow in the ideal symmetric stack; dots mark a 90% top / 60% bottom example.",
             ha="center", va="bottom", fontsize=8.3, color="#444")
    fig.tight_layout(rect=(0, .09, 1, .94))
    fig.savefig(FIG/"pcb_copper_imbalance.png", dpi=180)
    fig.savefig(FIG/"pcb_copper_imbalance.pdf")
    with (RAW/"pcb_copper_imbalance.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    # Show where the core-CTE ranking is sensitive to unknown dielectric
    # properties; this is an assumed-value grid, not a calibrated interval.
    e_values = np.linspace(10, 40, 13)
    a_values = np.linspace(5, 30, 21)
    advantage = np.zeros((len(e_values), len(a_values)))
    sensitivity_rows = []
    for i, dielectric_e in enumerate(e_values):
        for j, dielectric_cte in enumerate(a_values):
            p = replace(PCBScreenInputs(), core_thickness_mm=1.6,
                        top_copper_coverage=.90, bottom_copper_coverage=.60,
                        dielectric_e_gpa=float(dielectric_e),
                        dielectric_cte_ppm_k=float(dielectric_cte))
            low = screen_pcb(2.0, p, q)["bow_mm"]*1e3
            baseline = screen_pcb(baseline_cte, p, q)["bow_mm"]*1e3
            advantage[i, j] = low-baseline
            sensitivity_rows.append({"dielectric_E_GPa": dielectric_e,
                                     "dielectric_CTE_ppm_K": dielectric_cte,
                                     "bow_core_2p0_um": low,
                                     "bow_core_baseline_um": baseline,
                                     "low_minus_baseline_um": low-baseline})
    fig2, ax = plt.subplots(figsize=(8.8, 5.9))
    bound = float(np.max(np.abs(advantage)))
    mesh = ax.pcolormesh(a_values, e_values, advantage, shading="auto",
                         cmap="PiYG", vmin=-bound, vmax=bound)
    ax.contour(a_values, e_values, advantage, levels=[0], colors="#333", linewidths=1.2)
    ax.scatter([10], [20], color="#111", s=45, marker="x", zorder=4)
    ax.text(10.6, 20.4, "illustrative default", fontsize=8, color="#111")
    ax.set_xlabel("Assumed dielectric CTE (ppm/K)")
    ax.set_ylabel("Assumed dielectric modulus (GPa)")
    ax.set_title("Core choice also depends on unknown dielectric properties", weight="bold")
    cb = fig2.colorbar(mesh, ax=ax)
    cb.set_label("Bow with 2.0 core − bow with baseline core (µm)")
    fig2.text(.5, .018, "100 mm board, 1.6 mm core, 90% top / 60% bottom copper, 100 K cooling.\n"
              "Green: baseline core gives less bow. Purple: 2.0 core gives less bow. Assumed-value screen only.",
              ha="center", va="bottom", fontsize=8.5, color="#444")
    fig2.tight_layout(rect=(0, .09, 1, 1))
    fig2.savefig(FIG/"pcb_dielectric_sensitivity.png", dpi=180)
    fig2.savefig(FIG/"pcb_dielectric_sensitivity.pdf")
    with (RAW/"pcb_dielectric_sensitivity.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=sensitivity_rows[0].keys())
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    for thickness in (1.6, 3.0):
        for _, core_cte, _ in scenarios:
            p = replace(PCBScreenInputs(), core_thickness_mm=thickness,
                        top_copper_coverage=.90, bottom_copper_coverage=.60)
            bow = screen_pcb(core_cte, p, q)["bow_mm"]*1e3
            print(f"core {thickness:.1f} mm, CTE {core_cte:.2f}: "
                  f"{bow:.1f} µm per 100 K under illustrative 90/60 copper coverage")


if __name__ == "__main__":
    main()
