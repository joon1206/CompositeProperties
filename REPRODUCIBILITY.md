# Reproducibility record

1. Install: `python -m pip install -e ".[test,jax]"` from `C:\comp_work`.
2. Verify: `python -m pytest -q` (13 tests on the recorded run).
3. Reproduce all calculations and figures:
   `python -m composite_research.experiments`.
4. Compare `results/raw/manifest.json` with the source hash and inputs;
   all raw CSVs are deterministic except measured run times.
5. Refresh bibliography only if desired:
   `python scripts/build_literature.py` (requires Crossref network access).

Original tested environment: Windows, Python 3.14.0, NumPy 2.3.3, SciPy 1.16.2,
JAX 0.12.0.1 CPU backend; AMD Family 25 Model 33 CPU. The later accelerated
run used WSL2 Ubuntu, Python 3.14.4, JAX 0.12.0.1, the same Ryzen 5 5600X
CPU and an RTX 2080 Ti (`cuda:0`) with x64 enabled. The randomized scenario
design is scrambled Sobol, seed
20261024, 512 draws over the eight ranges encoded in `experiments.py`.
The 493-point grid and 17 sign-reversal roots are deterministic.

The data are predictions under linear elasticity, ideal flat bonded layers,
transversely isotropic ply approximations, and two alternative isotropic
filled-matrix closures. No matched experimental data, full fiber/particle
RVE, or cluster execution is claimed. Subdivision checks apply only to
the analytically exact flat-layer constrained-energy problem.
The LaTeX manuscript source was written and its bibliography generated, but
no TeX compiler was installed here, so a compiled manuscript PDF was not
verified.

Raw result files: `baseline.csv`, `grid.csv`, `uncertainty.csv`,
`uncertainty_summary.json`, `derivatives.csv`, `convergence.csv`,
`cte_bias_sign_reversal.csv`, `benchmark.json`, `jax_benchmark.json`,
and `manifest.json`.

## Accelerated study

The WSL2 GPU virtual environment is
`/home/joon1206/venvs/composite-gpu-capped`. The interactive shell may source
OpenFOAM's old `libstdc++` through `LD_LIBRARY_PATH`; unset that variable for
JAX. The explicit `JAX_PLATFORMS` selection prevents a silent CPU fallback.
From `C:\comp_work` in PowerShell:

```powershell
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/comp_work && env -u LD_LIBRARY_PATH XLA_PYTHON_CLIENT_PREALLOCATE=false JAX_PLATFORMS=cuda /home/joon1206/venvs/composite-gpu-capped/bin/python -m pytest -q'
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/comp_work && env -u LD_LIBRARY_PATH XLA_PYTHON_CLIENT_PREALLOCATE=false JAX_PLATFORMS=cuda /home/joon1206/venvs/composite-gpu-capped/bin/python scripts/accelerated_study.py --backend gpu --extended'
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/comp_work && env -u LD_LIBRARY_PATH XLA_PYTHON_CLIENT_PREALLOCATE=false JAX_PLATFORMS=cpu /home/joon1206/venvs/composite-gpu-capped/bin/python scripts/accelerated_study.py --backend cpu'
wsl.exe -d Ubuntu -- bash -lc 'cd /mnt/c/comp_work && env -u LD_LIBRARY_PATH /home/joon1206/venvs/composite-gpu-capped/bin/python scripts/make_executive_figures.py'
```

The benchmark computes four model variants and five outputs per scenario.
For each size, it records compilation plus first call, nine warmed device
resident calls, and five warmed host input/output round trips. All calls are
synchronized before timing ends. The comparison chart uses median round-trip
times. The CPU and GPU runs use identical Sobol inputs and JAX code. The
NumPy scalar baseline is also measured for up to 2,048 scenarios. The
8,192-case GPU extension uses the same illustrative ranges as the original
512-case sweep. Its JSON includes maximum differences against the saved
grid and uncertainty CSVs. Source: `scripts/accelerated_study.py`;
chart source: `scripts/make_executive_figures.py`.

New artifacts: `performance_cpu.json`, `performance_gpu.json`, and
`extended_study_gpu.json` under `results/raw`; `executive_gpu_speedup` and
`executive_model_risk` as PNG/PDF under `results/figures`.
