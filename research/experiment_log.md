# Experiment log

| Experiment | Question | Actual execution | Main result | Raw output |
|---|---|---|---|---|
| Baseline | Which algebraic stage changes outputs? | Four model chains | Matrix closure largest stiffness change; legacy thermal equation changes AY | `results/raw/baseline.csv` |
| Grid | Does error persist over composition? | 17×29=493 cases | Matrix EX discrepancy up to 13.53%; Voigt GXY bias up to 1.57% | `results/raw/grid.csv` |
| Scenario sweep | Are rankings stable under illustrative input changes? | 512 scrambled Sobol draws | Median matrix EX discrepancy 7.04%; layer EX 0.0012% | `results/raw/uncertainty.csv` |
| Sign reversal | Is legacy AY bias one-sided? | 17 bracketing/root solves | One crossing for each Vf 0.2–0.6, phi 0.206–0.283 | `results/raw/cte_bias_sign_reversal.csv` |
| Numerical reference | Are block formulas assembled correctly? | 1–16 sublayers | Relative C discrepancy about 10⁻¹⁵ or less | `results/raw/convergence.csv` |
| Differentiation | Do automatic derivatives agree? | JAX jacfwd vs centered differences | All five outputs/eight inputs pass tests | `tests/test_jax.py` |
| Performance | Is JAX useful for batches here? | CPU 512 JIT batch | 0.55 s initial; 1.6 ms steady median; no GPU backend | `results/raw/jax_benchmark.json` |
