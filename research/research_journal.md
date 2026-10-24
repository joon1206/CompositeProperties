# Research journal

## 2026-10-24 — audit and first cycle

Inspected the two supplied text files and an empty workspace. The MATLAB
paste contains apparent lost multiplication symbols, a stiffness average
presented as physical laminate homogenization, independent laminate CTE
equations, and a mislabeled `NUYZ` expression. Chose a layered-model audit
after literature revealed close prior exact and hierarchical methods.
Implemented pure NumPy exact layer blocks and an independent energy
reference. Initial result refuted a large balanced cross-ply Voigt error.

## 2026-10-24 — revised physical comparison

Added bonded-sphere Mori–Tanaka and Hashin–Shtrikman bounds. Ran grid and
scenario sweep. The particulate stiffness closure dominates layer error in
the tested parameter windows. An error sign change in the original
thickness CTE led to a post-hoc root-search experiment. Added JAX x64
autodifferentiation and checked derivatives against finite differences.

## 2026-10-24 — skeptical review

Mathematical criticism: source Poisson label reversed and original CTE was
not a consistent thermal-stress mixture; fixed in reference outputs.
Novelty criticism: exact thermoelastic layers and silica/fiber RVEs already
published; withdrew novelty claims. Validity criticism: constituent data,
packing, interfaces, and actual laminate measurements unavailable; reported
only model comparisons. Publication potential remains unproven.
