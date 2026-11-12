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

## Package-level screening study

The [package warpage screen](research/package_warpage_screen.md) connects the
computed core stiffness to a reduced-order thermal plate model using the
project deck's core, die and ring dimensions. It compares the two conflicting
Resonac Type LH in-plane CTE ranges and sweeps the unknown effective build-up
CTE. Its main result is that adding build-up layers can reverse the preferred
core CTE in the no-ring case. It is a conditional sensitivity model, not a
validated package warpage forecast.

```powershell
python scripts/package_warpage_screen.py
```

Figures and raw scenario tables are saved in `results/figures` and
`results/raw`. The plate implementation and checks are in
`src/composite_research/package_screen.py` and `tests/test_package_screen.py`.
