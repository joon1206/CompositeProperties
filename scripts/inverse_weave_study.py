"""Reproducible inverse screen, synthetic periodic weave study, and material cards.

Run from repository root: python scripts/inverse_weave_study.py
All input ranges and temperature scenarios are assumptions, not measurements.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import qmc

from composite_research.model import Inputs
from composite_research.woven import woven_orientation
from composite_research.weave_cell import CellSpec, homogenize, properties, synthetic_cell

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "results" / "raw"
FIG = ROOT / "results" / "figures"
KEYS = ("EX", "EY", "EZ", "GXY", "GYZ", "GXZ", "NUXY", "NUYZ", "NUXZ", "AX", "AY", "AZ")
FIBERS = {
    "E_glass_illustrative": (76., .21, 5.6),
    "S_glass_illustrative": (87., .22, 2.8),
    "original_fiber_illustrative": (82., .22, 2.49),
}
BOUNDS = {"vf": (.25, .65), "phi": (.25, .75), "er": (1.5, 5.0),
          "ar": (25., 65.), "ep": (65., 75.), "ap": (.3, .7),
          "warp_crimp_deg": (0., 20.), "weft_crimp_deg": (0., 20.),
          "warp_share": (.40, .60)}
WINDOWS = {"overview": ((31., 33.), (4.5, 6.5)),
           "detail": ((31., 33.), (.5, 3.0))}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def inverse_screen():
    # Same 2048 Sobol points across six fiber/closure groups allow paired comparisons.
    unit = qmc.Sobol(d=len(BOUNDS), scramble=True, seed=20261113).random_base2(11)
    names = list(BOUNDS)
    rows = []
    for fiber_name, (ef, nuf, af) in FIBERS.items():
        for method in ("legacy", "mori_tanaka"):
            for idx, uu in enumerate(unit):
                settings = {name: lo + uu[k] * (hi-lo)
                            for k, (name, (lo, hi)) in enumerate(BOUNDS.items())}
                p = replace(Inputs(), ef=ef, nuf=nuf, af=af,
                            **{k: settings[k] for k in ("vf", "phi", "er", "ar", "ep", "ap")})
                r = woven_orientation(p, matrix_method=method,
                                      warp_crimp_deg=settings["warp_crimp_deg"],
                                      weft_crimp_deg=settings["weft_crimp_deg"],
                                      warp_share=settings["warp_share"])
                row = {"sample_id": idx, "fiber": fiber_name, "matrix_method": method,
                       "ef": ef, "nuf": nuf, "af": af, **settings,
                       **{k: float(r[k]) for k in KEYS}}
                for label, (ew, aw) in WINDOWS.items():
                    row[f"meets_{label}_X"] = ew[0] <= row["EX"] <= ew[1] and aw[0] <= row["AX"] <= aw[1]
                    row[f"meets_{label}_XZ"] = (ew[0] <= row["EX"] <= ew[1] and
                                                ew[0] <= row["EZ"] <= ew[1] and
                                                aw[0] <= row["AX"] <= aw[1] and
                                                aw[0] <= row["AZ"] <= aw[1])
                rows.append(row)
            print("inverse", fiber_name, method, flush=True)
    write_csv(RAW / "inverse_all_samples.csv", rows)
    for label in WINDOWS:
        feasible = [r for r in rows if r[f"meets_{label}_XZ"]]
        if feasible:
            write_csv(RAW / f"inverse_feasible_{label}.csv", feasible)
    return rows


def weave_screen():
    rows = []
    # Three silica partitions at each of two crimp amplitudes and two meshes.
    # Cell parameters are synthetic and chosen to avoid tow intersection.
    for amplitude in (0.06, 0.10):
        for label, phi_in in (("uniform", None), ("between_tows", 0.),
                              ("inside_tows", .7)):
            for n in ((8, 12) if amplitude == .06 else (8, 12, 16, 20)):
                spec = CellSpec(n=n, amplitude=amplitude, silica_in_tows=phi_in)
                c, a, meta = synthetic_cell(Inputs(), spec)
                ce, ae, diag = homogenize(c, a)
                rve = properties(ce, ae)
                # Match actual discretized global composition and RMS tow slope.
                # This removes a trivial composition difference, but the surrogate
                # still cannot represent resin pockets or silica location.
                p = replace(Inputs(), vf=meta["fiber_fraction"],
                            phi=spec.global_matrix_silica_fraction)
                fast = woven_orientation(p, matrix_method=spec.matrix_method,
                                         warp_crimp_deg=meta["rms_tow_angle_deg"],
                                         weft_crimp_deg=meta["rms_tow_angle_deg"])
                row = {"amplitude": amplitude, "silica_location": label, "n": n,
                       **meta, **diag,
                       **{f"rve_{k}": float(rve[k]) for k in KEYS},
                       **{f"fast_{k}": float(fast[k]) for k in KEYS}}
                rows.append(row)
                print("weave", amplitude, label, n, "EX", rve["EX"],
                      "AX", rve["AX"], flush=True)
    write_csv(RAW / "synthetic_weave_cells.csv", rows)
    return rows


def material_cards(rows):
    # Three actual 25 C samples, chosen by low/median/high X CTE among
    # designs that meet the overview window in both in-plane directions.
    feasible = sorted((r for r in rows if r["meets_overview_XZ"]), key=lambda r: r["AX"])
    if len(feasible) < 3:
        # This fallback still gives coherent scenarios and is explicitly marked.
        feasible = sorted((r for r in rows if 31 <= r["EX"] <= 33), key=lambda r: r["AX"])
        selection = "modulus_only_fallback"
    else:
        selection = "overview_XZ_feasible"
    picked = {"low": feasible[0], "central": feasible[len(feasible)//2],
              "high": feasible[-1]}
    profiles = {25: (1., 0.), 100: (.85, 5.), 150: (.65, 12.)}
    cards = []
    for label, base in picked.items():
        for temp, (er_scale, ar_add) in profiles.items():
            p = replace(Inputs(), ef=base["ef"], nuf=base["nuf"], af=base["af"],
                        vf=base["vf"], phi=base["phi"], er=base["er"]*er_scale,
                        ar=base["ar"]+ar_add, ep=base["ep"], ap=base["ap"])
            result = woven_orientation(p, matrix_method=base["matrix_method"],
                                       warp_crimp_deg=base["warp_crimp_deg"],
                                       weft_crimp_deg=base["weft_crimp_deg"],
                                       warp_share=base["warp_share"])
            c = result["C"]
            s = np.linalg.inv(c)
            card = {"card": label, "temperature_C": temp, "selection": selection,
                    "sample_id": base["sample_id"], "fiber": base["fiber"],
                    "matrix_method": base["matrix_method"],
                    "er_scale_vs_25C": er_scale, "ar_add_vs_25C_ppmK": ar_add,
                    **{k: base[k] for k in ("ef", "nuf", "af", "vf", "phi", "ep", "ap",
                                           "warp_crimp_deg", "weft_crimp_deg", "warp_share")},
                    "er": p.er, "ar": p.ar,
                    **{k: float(result[k]) for k in KEYS},
                    "min_C_eigenvalue_GPa": float(np.linalg.eigvalsh(c).min()),
                    "max_C_symmetry_error_GPa": float(np.abs(c-c.T).max()),
                    "reciprocity_XY_error": float(abs(s[0,1]-s[1,0])),
                    "reciprocity_YZ_error": float(abs(s[1,2]-s[2,1])),
                    "reciprocity_XZ_error": float(abs(s[0,2]-s[2,0]))}
            cards.append(card)
    write_csv(RAW / "illustrative_12_property_cards.csv", cards)
    # Componentwise ranges across all 25 C overview-feasible designs. Each
    # extremum/median retains the actual design that produced it; these are
    # sensitivity envelopes, not simultaneous physical material cards.
    range_rows = []
    for temp, (er_scale, ar_add) in profiles.items():
        predictions = []
        for base in feasible:
            p = replace(Inputs(), ef=base["ef"], nuf=base["nuf"], af=base["af"],
                        vf=base["vf"], phi=base["phi"], er=base["er"]*er_scale,
                        ar=base["ar"]+ar_add, ep=base["ep"], ap=base["ap"])
            result = woven_orientation(p, matrix_method=base["matrix_method"],
                                       warp_crimp_deg=base["warp_crimp_deg"],
                                       weft_crimp_deg=base["weft_crimp_deg"],
                                       warp_share=base["warp_share"])
            predictions.append((base,result))
        for key in KEYS:
            ranked = sorted(predictions,key=lambda pair:pair[1][key])
            row = {"temperature_C":temp,"property":key,
                   "unit":"GPa" if key.startswith(("E","G")) else
                          "ppm/K" if key.startswith("A") else "dimensionless",
                   "scenario_count":len(ranked)}
            for label, (base,result) in zip(("low","central","high"),
                                            (ranked[0],ranked[len(ranked)//2],ranked[-1])):
                row.update({f"{label}_value":float(result[key]),
                            f"{label}_fiber":base["fiber"],
                            f"{label}_matrix_method":base["matrix_method"],
                            f"{label}_sample_id":base["sample_id"],
                            f"{label}_vf":base["vf"], f"{label}_phi":base["phi"],
                            f"{label}_warp_crimp_deg":base["warp_crimp_deg"],
                            f"{label}_weft_crimp_deg":base["weft_crimp_deg"],
                            f"{label}_warp_share":base["warp_share"]})
            range_rows.append(row)
    write_csv(RAW / "illustrative_property_envelopes.csv", range_rows)
    return cards


def figures(rows, cells, cards):
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    color = {"E_glass_illustrative": "#2a6f9e", "S_glass_illustrative": "#d27627",
             "original_fiber_illustrative": "#7461a8"}
    for name in FIBERS:
        part = [r for r in rows if r["fiber"] == name and r["matrix_method"] == "mori_tanaka"]
        ax.scatter([r["EX"] for r in part], [r["AX"] for r in part], s=5, alpha=.20,
                   color=color[name], label=name.replace("_illustrative", "").replace("_", " "))
    from matplotlib.patches import Rectangle
    for label, (_, aw) in WINDOWS.items():
        ax.add_patch(Rectangle((31, aw[0]), 2, aw[1]-aw[0], fill=False,
                               lw=2, ls="-" if label=="overview" else "--",
                               color="#157f54" if label=="overview" else "#b03d3d"))
        ax.text(33.2, (aw[0]+aw[1])/2, label, va="center", fontsize=9)
    ax.set(xlabel="Predicted warp Young's modulus EX (GPa)",
           ylabel="Predicted warp CTE AX (ppm/K)",
           title="Inverse screen: assumed designs versus comparison windows")
    ax.legend(frameon=False, ncol=3, fontsize=8, loc="upper right")
    ax.grid(alpha=.15)
    fig.text(.5, .01, "Sobol sample; Mori–Tanaka matrix closure shown. Windows are product-level comparisons, not validation.",
             ha="center", fontsize=8, color="#555")
    fig.tight_layout(rect=(0,.04,1,1))
    fig.savefig(FIG / "inverse_design_tradeoffs.png", dpi=180)
    fig.savefig(FIG / "inverse_design_tradeoffs.pdf")
    plt.close(fig)

    feasible = [r for r in rows if r["meets_overview_XZ"]]
    fig, axes = plt.subplots(1, 2, figsize=(9.7, 4.5))
    for name in FIBERS:
        for method, marker in (("legacy", "o"), ("mori_tanaka", "s")):
            subset = [r for r in feasible if r["fiber"]==name and r["matrix_method"]==method]
            label = name.replace("_illustrative", "").replace("_", " ") + ", " + method.replace("_", " ")
            axes[0].scatter([r["vf"] for r in subset], [r["phi"] for r in subset],
                            c=color[name], marker=marker, s=45, alpha=.85, label=label)
            axes[1].scatter([r["AX"] for r in subset], [r["AY"] for r in subset],
                            c=color[name], marker=marker, s=45, alpha=.85)
    axes[0].set(xlabel="Fiber fraction of composite", ylabel="Silica fraction of nonfiber matrix",
                title="Many compositions meet the broad window")
    axes[1].set(xlabel="In-plane CTE AX (ppm/K)", ylabel="Thickness CTE AY (ppm/K)",
                title="Thickness expansion also changes")
    axes[0].legend(frameon=False, fontsize=7, loc="lower left")
    for ax in axes: ax.grid(alpha=.15)
    fig.tight_layout()
    fig.savefig(FIG / "inverse_feasible_tradeoffs.png", dpi=180)
    fig.savefig(FIG / "inverse_feasible_tradeoffs.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.7, 4.7))
    for amp, marker in ((.06, "o"), (.10, "s")):
        for loc, color2 in (("uniform", "#687681"), ("between_tows", "#187a5a"),
                             ("inside_tows", "#b85b31")):
            rr = next(r for r in cells if r["amplitude"]==amp and r["n"]==(12 if amp == .06 else 20)
                      and r["silica_location"]==loc)
            for ax, key, unit in ((axes[0], "EX", "GPa"), (axes[1], "AX", "ppm/K")):
                ax.scatter(rr[f"fast_{key}"], rr[f"rve_{key}"], marker=marker,
                           s=80, color=color2, edgecolor="white", linewidth=.5,
                           label=f"{loc.replace('_',' ')}, a={amp}" if ax is axes[0] else None)
    for ax, key in zip(axes, ("EX", "AX")):
        vals = [r[f"{m}_{key}"] for r in cells for m in ("fast", "rve")]
        lo, hi = min(vals)*.95, max(vals)*1.05
        ax.plot([lo,hi],[lo,hi], color="#333", lw=1, ls="--")
        ax.set(xlim=(lo,hi), ylim=(lo,hi), xlabel=f"Fast orientation {key}",
               ylabel=f"Periodic cell {key}")
        ax.grid(alpha=.15)
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Geometry and silica location change the predicted properties")
    fig.tight_layout()
    fig.savefig(FIG / "synthetic_weave_comparison.png", dpi=180)
    fig.savefig(FIG / "synthetic_weave_comparison.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.4))
    for card, color2 in (("low", "#187a5a"), ("central", "#2a6f9e"), ("high", "#b85b31")):
        cc = [r for r in cards if r["card"]==card]
        for ax, key in zip(axes, ("EX", "AX")):
            ax.plot([r["temperature_C"] for r in cc], [r[key] for r in cc],
                    marker="o", color=color2, label=card)
    axes[0].set_ylabel("Warp Young's modulus (GPa)")
    axes[1].set_ylabel("Warp CTE (ppm/K)")
    for ax in axes:
        ax.set_xlabel("Temperature (°C)")
        ax.grid(alpha=.2)
    axes[0].legend(frameon=False)
    fig.suptitle("Illustrative temperature scenarios; no measured resin curve")
    fig.tight_layout()
    fig.savefig(FIG / "illustrative_card_temperature.png", dpi=180)
    fig.savefig(FIG / "illustrative_card_temperature.pdf")
    plt.close(fig)


def summary(rows, cells, cards):
    paired = {}
    for r in rows:
        paired.setdefault((r["fiber"],r["sample_id"]),{})[r["matrix_method"]] = r
    out = {"sample_count": len(rows), "seed": 20261113, "sobol_points_per_group": 2048,
           "bounds": BOUNDS, "fiber_presets": FIBERS, "windows": WINDOWS,
           "counts": {label: {"X": sum(r[f"meets_{label}_X"] for r in rows),
                              "XZ": sum(r[f"meets_{label}_XZ"] for r in rows)}
                      for label in WINDOWS},
           "paired_designs_meeting_overview_under_both_matrix_models": sum(
               all(v[f"meets_overview_XZ"] for v in pair.values()) for pair in paired.values()),
           "counts_by_group": {f"{f}/{m}": {label: sum(r[f"meets_{label}_XZ"] for r in rows
                    if r["fiber"]==f and r["matrix_method"]==m) for label in WINDOWS}
                    for f in FIBERS for m in ("legacy", "mori_tanaka")},
           "weave_mesh_sensitivity": [],
           "card_selection": {c["card"]: {"fiber": c["fiber"], "sample_id": c["sample_id"],
                                           "matrix_method": c["matrix_method"]}
                              for c in cards if c["temperature_C"]==25}}
    for amp in (.06, .10):
        for loc in ("uniform", "between_tows", "inside_tows"):
            sequence = sorted((r for r in cells if r["amplitude"]==amp and
                               r["silica_location"]==loc), key=lambda r:r["n"])
            for a,b in zip(sequence[:-1],sequence[1:]):
                out["weave_mesh_sensitivity"].append({"amplitude": amp, "silica_location": loc,
                    "mesh_pair": f"{a['n']}→{b['n']}",
                    "EX_percent_change": 100*(b["rve_EX"]-a["rve_EX"])/b["rve_EX"],
                    "AX_percent_change": 100*(b["rve_AX"]-a["rve_AX"])/b["rve_AX"],
                    "fiber_fraction_delta": b["fiber_fraction"]-a["fiber_fraction"]})
    (RAW / "inverse_weave_summary.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reuse-inverse", action="store_true",
                        help="reuse existing inverse_all_samples.csv after an interrupted weave run")
    parser.add_argument("--reuse-cells", action="store_true",
                        help="reuse existing synthetic_weave_cells.csv to refresh figures and summary")
    args = parser.parse_args()
    RAW.mkdir(exist_ok=True, parents=True)
    FIG.mkdir(exist_ok=True, parents=True)
    if args.reuse_inverse:
        with (RAW / "inverse_all_samples.csv").open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        for r in rows:
            for key in ("sample_id",):
                r[key] = int(r[key])
            for key in (*BOUNDS, "ef", "nuf", "af", *KEYS):
                r[key] = float(r[key])
            for label in WINDOWS:
                for suffix in ("X", "XZ"):
                    r[f"meets_{label}_{suffix}"] = r[f"meets_{label}_{suffix}"] == "True"
    else:
        rows = inverse_screen()
    if args.reuse_cells:
        with (RAW / "synthetic_weave_cells.csv").open(newline="", encoding="utf-8") as f:
            cells = list(csv.DictReader(f))
        for r in cells:
            for key in r:
                if key not in ("silica_location", "solver", "mean_phase_cte"):
                    r[key] = float(r[key])
            r["n"] = int(r["n"])
    else:
        cells = weave_screen()
    cards = material_cards(rows)
    figures(rows, cells, cards)
    print(json.dumps(summary(rows, cells, cards)["counts"], indent=2))


if __name__ == "__main__":
    main()
