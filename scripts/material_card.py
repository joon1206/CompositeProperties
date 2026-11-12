"""Turn temperature-specific constituent inputs into a 12-property card.

Example:
  python scripts/material_card.py data/material_card_input_example.csv results/raw/material_card_example.csv --matrix-method mori_tanaka

Each input row is one measured or assumed temperature state. The script does
not interpolate, infer a glass transition, or fit missing measurements.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

from composite_research.model import Inputs
from composite_research.woven import woven_orientation


COMMON_FIELDS = ("temperature_C", "ef", "nuf", "af", "vf",
                 "warp_crimp_deg", "weft_crimp_deg", "warp_share")
ESTIMATE_FIELDS = ("ep", "nup", "ap", "phi", "er", "nur", "ar")
MEASURED_FIELDS = ("em_measured", "num_measured", "am_measured")
PROPERTIES = ("EX", "EY", "EZ", "GXY", "GYZ", "GXZ",
              "NUXY", "NUYZ", "NUXZ", "AX", "AY", "AZ")


def build_card(input_path: Path, output_path: Path, matrix_method: str):
    if matrix_method not in ("legacy", "mori_tanaka", "measured"):
        raise ValueError("unknown matrix method")
    required = COMMON_FIELDS+(MEASURED_FIELDS if matrix_method == "measured"
                              else ESTIMATE_FIELDS)
    with input_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = set(required)-set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"input CSV missing columns: {', '.join(sorted(missing))}")
        input_rows = list(reader)
    if not input_rows:
        raise ValueError("input CSV has no data rows")
    output_rows = []
    for index, row in enumerate(input_rows, start=2):
        try:
            values = {key: float(row[key]) for key in required}
        except (TypeError, ValueError) as exc:
            raise ValueError(f"non-numeric or empty value in CSV row {index}") from exc
        p = Inputs(ef=values["ef"], nuf=values["nuf"], af=values["af"],
                   vf=values["vf"],
                   **{key: values[key] for key in ESTIMATE_FIELDS if key in values})
        result = woven_orientation(
            p, matrix_method=matrix_method,
            warp_crimp_deg=values["warp_crimp_deg"],
            weft_crimp_deg=values["weft_crimp_deg"],
            warp_share=values["warp_share"],
            filled_matrix_properties=(tuple(values[key] for key in MEASURED_FIELDS)
                                      if matrix_method == "measured" else None))
        output_rows.append({"temperature_C": values["temperature_C"],
                            "matrix_method": matrix_method,
                            "warp_crimp_deg": values["warp_crimp_deg"],
                            "weft_crimp_deg": values["weft_crimp_deg"],
                            "warp_share": values["warp_share"],
                            "Em": result["ply"]["Em"],
                            "num": result["ply"]["num"],
                            "am": result["ply"]["am"],
                            **{key: result[key] for key in PROPERTIES}})
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=output_rows[0].keys())
        writer.writeheader()
        writer.writerows(output_rows)
    return output_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--matrix-method", choices=("legacy", "mori_tanaka", "measured"),
                        default="mori_tanaka")
    args = parser.parse_args()
    rows = build_card(args.input_csv, args.output_csv, args.matrix_method)
    print(f"Wrote {len(rows)} temperature state(s) to {args.output_csv}")


if __name__ == "__main__":
    main()
