"""Generate the conditional core/die/ring warpage screening study.

Run: python scripts/package_warpage_screen.py
"""
from __future__ import annotations

import csv
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from composite_research.model import Inputs, evaluate
from composite_research.package_screen import ScreenInputs, core_q_alpha, solve_screen


ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw"
FIG = ROOT / "results" / "figures"
WINDOWS = {"Detailed table (30–120 °C)": 1.75,
           "Overview table (below Tg)": 5.5}


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    q, _ = core_q_alpha()
    original = evaluate(Inputs())
    rows = []
    fig, axs = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True, sharey=True)
    colors = {False: "#2171b5", True: "#d95f0e"}
    for row_index, thickness in enumerate((1.6, 3.0)):
        for col_index, ring in enumerate((False, True)):
            settings = replace(ScreenInputs(), core_thickness=thickness)
            ax = axs[row_index, col_index]
            cte_values = np.linspace(0, 10, 41)
            responses = [solve_screen(float(cte), q, ring, settings)["bow_mm"]*1e3
                         for cte in cte_values]
            ax.plot(cte_values, responses, color=colors[ring], lw=2.3)
            ax.axhline(0, color="#555", lw=.8)
            ax.axvspan(.5, 3.0, color="#6a51a3", alpha=.11)
            ax.axvspan(4.5, 6.5, color="#31a354", alpha=.10)
            ax.axvline(original["AX"], color="#444", ls="--", lw=1)
            ax.set_title(f"Core {thickness:.1f} mm · {'steel ring' if ring else 'no ring'}")
            ax.grid(alpha=.18)
            for source, cte in WINDOWS.items():
                r = solve_screen(cte, q, ring, settings)
                rows.append({"core_thickness_mm": thickness, "ring": ring,
                             "CTE_source": source, "core_CTE_ppm_K": cte,
                             "corner_minus_center_um": r["corner_minus_center_mm"]*1e3,
                             "bow_um": r["bow_mm"]*1e3,
                             "die_E_GPa": settings.die_e, "die_CTE_ppm_K": settings.die_cte,
                             "ring_E_GPa": settings.ring_e if ring else "",
                             "ring_CTE_ppm_K": settings.ring_cte if ring else "",
                             "stress_free_C": settings.stress_free_temp,
                             "observation_C": settings.observation_temp})
                ax.scatter([cte], [r["bow_mm"]*1e3],
                           color="#6a51a3" if cte < 4 else "#31a354", zorder=3, s=35)
    for ax in axs[-1]:
        ax.set_xlabel("Assumed in-plane core CTE (ppm/K)")
    fig.supylabel("Peak-to-valley bow (µm; lower is better)", x=.012, fontsize=11)
    fig.suptitle("Illustrative core CTE screen: material choice depends on package design", fontsize=13, weight="bold")
    fig.text(.5, .025, "Free-plate screen: core + effective die patch + optional edge ring only.\n"
             "Purple/green shading = two Resonac catalog CTE windows; dashed line = repository baseline. Not a package forecast.",
             ha="center", va="bottom", fontsize=8.5, color="#444")
    fig.tight_layout(rect=(.04, .075, 1, .95))
    fig.savefig(FIG/"package_warpage_screen.png", dpi=180)
    fig.savefig(FIG/"package_warpage_screen.pdf")
    with (RAW/"package_warpage_screen.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    sensitivity = []
    for thickness in (1.6, 3.0):
        for ring in (False, True):
            for die_cte in (3.5, 4.4, 5.3):
                for ring_cte in ((10.0, 12.0, 14.0) if ring else (12.0,)):
                    settings = replace(ScreenInputs(), core_thickness=thickness,
                                       die_cte=die_cte, ring_cte=ring_cte)
                    for source, core_cte in WINDOWS.items():
                        r = solve_screen(core_cte, q, ring, settings)
                        sensitivity.append({"core_thickness_mm": thickness, "ring": ring,
                                            "die_CTE_ppm_K": die_cte,
                                            "ring_CTE_ppm_K": ring_cte if ring else "",
                                            "CTE_source": source, "core_CTE_ppm_K": core_cte,
                                            "bow_um": r["bow_mm"]*1e3})
    with (RAW/"package_warpage_sensitivity.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=sensitivity[0].keys())
        writer.writeheader()
        writer.writerows(sensitivity)
    # The deck specifies build-up thickness but not an effective E or CTE.
    # Sweep those missing inputs explicitly; the plotted E=20 GPa is a
    # hypothetical sensitivity value, not a measured package property.
    fig2, axs2 = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex=True, sharey=True)
    build_rows = []
    build_ctes = np.linspace(0, 20, 41)
    for row_index, thickness in enumerate((1.6, 3.0)):
        for col_index, ring in enumerate((False, True)):
            ax = axs2[row_index, col_index]
            curves = {}
            for label, core_cte, color in (("Detailed range midpoint: 1.75", 1.75, "#6a51a3"),
                                           ("Overview range midpoint: 5.5", 5.5, "#31a354")):
                bow = []
                for build_cte in build_ctes:
                    settings = replace(ScreenInputs(), core_thickness=thickness,
                                       build_up_thickness=.42, build_up_e=20,
                                       build_up_nu=.3, build_up_cte=float(build_cte))
                    value = solve_screen(core_cte, q, ring, settings)["bow_mm"]*1e3
                    bow.append(value)
                    build_rows.append({"core_thickness_mm": thickness, "ring": ring,
                                       "core_CTE_ppm_K": core_cte,
                                       "build_up_thickness_each_side_mm": .42,
                                       "build_up_E_GPa": 20, "build_up_nu": .3,
                                       "build_up_CTE_ppm_K": float(build_cte),
                                       "bow_um": value})
                curves[core_cte] = np.array(bow)
                ax.plot(build_ctes, bow, color=color, lw=2.3, label=label)
            difference = curves[1.75]-curves[5.5]
            switches = np.where(np.sign(difference[:-1]) != np.sign(difference[1:]))[0]
            if len(switches):
                i = switches[0]
                switch = float(build_ctes[i]-difference[i]*(build_ctes[i+1]-build_ctes[i])/
                               (difference[i+1]-difference[i]))
                ax.axvline(switch, color="#444", ls="--", lw=1)
                ax.text(.97, .96, f"Ranking switches near {switch:.1f} ppm/K",
                        transform=ax.transAxes, ha="right", va="top", fontsize=8,
                        bbox={"facecolor": "white", "edgecolor": "none", "alpha": .8})
                print(f"Core {thickness:.1f} mm, ring={ring}: ranking switch at build-up CTE ≈{switch:.1f} ppm/K")
            ax.set_title(f"Core {thickness:.1f} mm · {'steel ring' if ring else 'no ring'}")
            ax.grid(alpha=.18)
    axs2[0, 0].legend(loc="upper left", frameon=False, fontsize=8)
    for ax in axs2[-1]:
        ax.set_xlabel("Assumed effective build-up CTE (ppm/K)")
    fig2.supylabel("Peak-to-valley bow (µm; lower is better)", x=.012, fontsize=11)
    fig2.suptitle("Missing build-up properties can reverse the core choice", fontsize=14, weight="bold")
    fig2.text(.5, .025, "Hypothetical symmetric build-up: 0.42 mm per side, E=20 GPa, ν=0.30.\n"
              "Die and ring assumptions match the first chart. Screening result only; build-up properties are unmeasured.",
              ha="center", va="bottom", fontsize=8.5, color="#444")
    fig2.tight_layout(rect=(.04, .075, 1, .95))
    fig2.savefig(FIG/"package_buildup_sensitivity.png", dpi=180)
    fig2.savefig(FIG/"package_buildup_sensitivity.pdf")
    with (RAW/"package_buildup_sensitivity.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=build_rows[0].keys())
        writer.writeheader()
        writer.writerows(build_rows)
    thresholds = []
    for thickness in (1.6, 3.0):
        for build_e in (10.0, 20.0, 40.0):
            def bow_difference(build_cte):
                settings = replace(ScreenInputs(), core_thickness=thickness,
                                   build_up_thickness=.42, build_up_e=build_e,
                                   build_up_nu=.3, build_up_cte=build_cte)
                return (solve_screen(1.75, q, False, settings)["bow_mm"]
                        - solve_screen(5.5, q, False, settings)["bow_mm"])
            threshold = brentq(bow_difference, 0, 20, xtol=.002)
            thresholds.append({"core_thickness_mm": thickness,
                               "build_up_E_GPa": build_e,
                               "equal_bow_build_up_CTE_ppm_K": threshold})
    with (RAW/"package_buildup_thresholds.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=thresholds[0].keys())
        writer.writeheader()
        writer.writerows(thresholds)
    for row in rows:
        print(f"{row['core_thickness_mm']:.1f} mm, ring={row['ring']}, "
              f"CTE={row['core_CTE_ppm_K']:.2f}: signed={row['corner_minus_center_um']:.1f} µm, "
              f"bow={row['bow_um']:.1f} µm")


if __name__ == "__main__":
    main()
