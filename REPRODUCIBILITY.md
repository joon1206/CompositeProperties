# Reproducibility record

1. Install: `python -m pip install -e ".[test,jax]"` from `C:\comp_work`.
2. Verify: `python -m pytest -q` (13 tests on the recorded run).
3. Reproduce all calculations and figures:
   `python -m composite_research.experiments`.
4. Compare `results/raw/manifest.json` with the source hash and inputs;
   all raw CSVs are deterministic except measured run times.
5. Refresh bibliography only if desired:
   `python scripts/build_literature.py` (requires Crossref network access).

Tested environment: Windows, Python 3.14.0, NumPy 2.3.3, SciPy 1.16.2,
JAX 0.12.0.1 CPU backend; AMD Family 25 Model 33 CPU, RTX 2080 Ti detected
but not used. The randomized scenario design is scrambled Sobol, seed
20261024, 512 draws over the eight ranges encoded in `experiments.py`.
The 493-point grid and 17 sign-reversal roots are deterministic.

The data are predictions under linear elasticity, ideal flat bonded layers,
transversely isotropic ply approximations, and two alternative isotropic
filled-matrix closures. No matched experimental data, full fiber/particle
RVE, or GPU/cluster execution is claimed. Subdivision checks apply only to
the analytically exact flat-layer constrained-energy problem.
The LaTeX manuscript source was written and its bibliography generated, but
no TeX compiler was installed here, so a compiled manuscript PDF was not
verified.

Raw result files: `baseline.csv`, `grid.csv`, `uncertainty.csv`,
`uncertainty_summary.json`, `derivatives.csv`, `convergence.csv`,
`cte_bias_sign_reversal.csv`, `benchmark.json`, `jax_benchmark.json`,
and `manifest.json`.
