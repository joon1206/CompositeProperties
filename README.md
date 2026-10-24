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
on the tested Windows installation. Output CSV/JSON is under `results/raw`,
figures under `results/figures`, and literature metadata and bibliography
under `research/literature`. See [reproducibility notes](REPRODUCIBILITY.md).

Key source modules: `model.py` (constituent, ply, and exact layer models),
`reference.py` (independent constrained-energy calculation),
`jax_model.py` (differentiable x64 mirror), and `experiments.py` (all runs).

The original user paste is preserved in `data/source_matlab_paste.txt` for
traceability. The paste contains some missing multiplication operators in
the supplied text; the present implementation reconstructs the algebra
shown rather than claiming to have run the MATLAB file verbatim.
