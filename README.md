# Composite thermoelastic research audit

This repository executes a numerical study of the user's illustrative
glass-fiber / silica-filled resin 0/90 cross-ply. The [research report](manuscript/report.md)
contains the derivation, numerical findings, limitations, and closest prior
work. This is an audited model comparison, not validated material data.

## Run

Python 3.11+ is required. From the repository root:

```powershell
python -m pip install -e ".[test,jax]"
python -m pytest -q
python -m composite_research.experiments
```

The JAX extra is optional; the core study uses NumPy/SciPy. JAX runs on CPU
on the tested Windows installation. A separate WSL2 CUDA run used the
RTX 2080 Ti for a batched four-model extension. See the
[executive brief](EXECUTIVE_BRIEF.md) for measured CPU/GPU speedup and the
new model finding. Output CSV/JSON is under `results/raw`,
figures under `results/figures`, and literature metadata and bibliography
under `research/literature`. See [reproducibility notes](REPRODUCIBILITY.md).

Key source modules: `model.py` (constituent, ply, and exact layer models),
`reference.py` (independent constrained-energy calculation),
`jax_model.py` (differentiable x64 mirror), and `experiments.py` (all runs).
The accelerated study and chart scripts are in `scripts/`.

The original user paste is preserved in `data/source_matlab_paste.txt` for
traceability. The paste contains some missing multiplication operators in
the supplied text; the present implementation reconstructs the algebra
shown rather than claiming to have run the MATLAB file verbatim.

## Woven composite material properties

The [woven-material guide](research/woven_material_model.md) explains how to
calculate nine elastic engineering constants and three CTEs for the
resin/silica/glass composite itself. The orientation surrogate adds separate
warp and weft crimp angles to the existing 0/90 model and compares two
filled-resin descriptions. At zero crimp it reproduces the consistent
cross-ply baseline. This is an unvalidated material screen, not a resolved
weave unit cell.

```powershell
python scripts/woven_material_screen.py
python scripts/material_card.py data/material_card_input_example.csv results/raw/material_card_example.csv --matrix-method mori_tanaka
```

Figures and raw scenario tables are saved in `results/figures` and
`results/raw`. The material-card tool also accepts directly measured
filled-matrix properties; see the input template in `data/`. The model
and checks are in `src/composite_research/woven.py`,
`tests/test_woven.py` and `tests/test_material_card.py`.

## Inverse design and periodic weave cells

The [inverse and weave study](research/inverse_weave_study.md) searches
illustrative composition/property assumptions backward from two in-plane
modulus/CTE windows. It also compares the fast orientation model with
periodic finite-element homogenization of synthetic plain-weave cells,
including silica between or inside tows. Coherent 12-property cards at
25, 100 and 150 °C, all sampled inputs, numerical checks and figures are
saved under `results/`. None of these cards is a measured material.

```powershell
python scripts/inverse_weave_study.py
python scripts/refine_inverse.py
```

The FE cells and inverse screen use CPU NumPy/SciPy; they do not require
JAX or a GPU. A resumed run can use `--reuse-inverse` and `--reuse-cells`
to regenerate cards, figures and summary from existing CSV output.

The [filled-resin model deep dive](research/matrix_model_deep_dive.md)
compares the two existing estimates with Hashin–Shtrikman elastic and
Levin thermal bounds. Its key correction is that the two existing
filled-resin CTE formulas are algebraically identical. Regenerate its
figure and table with `python scripts/matrix_model_deep_dive.py`.
