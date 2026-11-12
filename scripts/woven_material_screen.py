"""Regenerate material-property sensitivity to matrix closure and tow crimp."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from composite_research.model import Inputs
from composite_research.woven import woven_orientation


ROOT = Path(__file__).resolve().parents[1]
RAW, FIG = ROOT/"results"/"raw", ROOT/"results"/"figures"
OUTPUTS = ("EX", "EY", "EZ", "GXY", "GYZ", "GXZ",
           "NUXY", "NUYZ", "NUXZ", "AX", "AY", "AZ")


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    p = Inputs()
    angles = np.arange(0, 21, 1)
    rows = []
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8))
    for method, label, color in (("legacy", "Supplied filled-resin formula", "#2171b5"),
                                 ("mori_tanaka", "Bonded-sphere alternative", "#d95f0e")):
        values = []
        for angle in angles:
            result = woven_orientation(p, matrix_method=method,
                                       warp_crimp_deg=float(angle),
                                       weft_crimp_deg=float(angle))
            values.append(result)
            rows.append({"matrix_method": method, "warp_crimp_deg": angle,
                         "weft_crimp_deg": angle, "warp_share": .5,
                         **{key: result[key] for key in OUTPUTS}})
        axes[0].plot(angles, [r["EX"] for r in values], lw=2.5, color=color, label=label)
        axes[1].plot(angles, [r["AX"] for r in values], lw=2.5, color=color, label=label)
    axes[0].set_ylabel("In-plane Young's modulus, EX (GPa)")
    axes[1].set_ylabel("In-plane thermal expansion, AX (ppm/K)")
    for ax in axes:
        ax.set_xlabel("Assumed equal warp/weft crimp angle (degrees)")
        ax.grid(alpha=.2)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("Material model: filled resin matters more than modest crimp", weight="bold", fontsize=13)
    fig.text(.5, .015, "Illustrative repository constituent inputs; crimp is an orientation surrogate, not measured fabric geometry.\n"
             "Solid lines are model predictions, not a fit to Resonac material.",
             ha="center", va="bottom", fontsize=8.2, color="#444")
    fig.tight_layout(rect=(0, .10, 1, .94))
    fig.savefig(FIG/"woven_model_comparison.png", dpi=180)
    fig.savefig(FIG/"woven_model_comparison.pdf")
    write_csv(RAW/"woven_material_screen.csv", rows)

    # Equal tow volume in warp and weft does not force equal directional
    # properties when their out-of-plane crimp differs.
    anisotropy = []
    weft_angles = np.arange(0, 21, 1)
    for angle in weft_angles:
        result = woven_orientation(p, matrix_method="mori_tanaka",
                                   warp_crimp_deg=5., weft_crimp_deg=float(angle))
        anisotropy.append({"warp_crimp_deg": 5., "weft_crimp_deg": angle,
                           "matrix_method": "mori_tanaka", "warp_share": .5,
                           **{key: result[key] for key in OUTPUTS}})
    fig2, axes2 = plt.subplots(1, 2, figsize=(10.4, 4.8))
    for key, color, label in (("EX", "#2171b5", "Warp X"),
                              ("EZ", "#6a51a3", "Weft Z")):
        axes2[0].plot(weft_angles, [r[key] for r in anisotropy], lw=2.5,
                      color=color, label=label)
    for key, color, label in (("AX", "#2171b5", "Warp X"),
                              ("AZ", "#6a51a3", "Weft Z")):
        axes2[1].plot(weft_angles, [r[key] for r in anisotropy], lw=2.5,
                      color=color, label=label)
    axes2[0].set_ylabel("Young's modulus (GPa)")
    axes2[1].set_ylabel("Thermal expansion (ppm/K)")
    for ax in axes2:
        ax.set_xlabel("Assumed weft crimp angle (degrees)")
        ax.grid(alpha=.2)
    axes2[0].legend(frameon=False)
    fig2.suptitle("Unequal warp/weft crimp creates directional properties", weight="bold", fontsize=13)
    fig2.text(.5, .015, "Warp crimp fixed at 5°; Mori–Tanaka filled-resin alternative; equal family share.\n"
              "Illustrative orientation model, not a measured fabric or a complete woven unit cell.",
              ha="center", va="bottom", fontsize=8.2, color="#444")
    fig2.tight_layout(rect=(0, .10, 1, .94))
    fig2.savefig(FIG/"woven_directional_properties.png", dpi=180)
    fig2.savefig(FIG/"woven_directional_properties.pdf")
    write_csv(RAW/"woven_anisotropy.csv", anisotropy)

    for method in ("legacy", "mori_tanaka"):
        for angle in (0, 10, 20):
            r = woven_orientation(p, matrix_method=method,
                                  warp_crimp_deg=angle, weft_crimp_deg=angle)
            print(f"{method}, {angle:2}° crimp: EX={r['EX']:.3f} GPa, "
                  f"AX={r['AX']:.3f} ppm/K, EY={r['EY']:.3f} GPa, AY={r['AY']:.3f} ppm/K")


if __name__ == "__main__":
    main()
