# Research status

## 2026-11-12 — Bare-PCB thermal laminate screen

Added a no-component, no-ring PCB laminate screen using the calculated core
stiffness and configurable top/bottom dielectric-copper build-up layers.
The ideal symmetric board has zero thermal bow under uniform temperature
change, while copper imbalance produces bow. A sensitivity grid shows that
dielectric properties can change the relative bow of core-CTE scenarios.
This is an uncalibrated screen; actual board stack-up, process history,
temperature-dependent properties and measured bare-board shapes are needed
for a physical prediction. See `research/pcb_warpage_screen.md`.

## 2026-10-24 — Composite model audit

Completed: source and hardware audit; 32 Crossref-verified scholarly records;
four explicit novelty reviews; exact flat-layer and independent numerical
reference; alternative particulate closure and Hashin–Shtrikman bounds;
13 passing tests; 493-case grid, 512-case scenario sweep, 17 root searches;
JAX x64 derivative checks and synchronized CPU batch benchmark; report,
figures, source manuscript and reproducibility record.

Strongest supported result: in this *idealized balanced cross-ply*, the
choice of filled-matrix closure changes stiffness much more than replacing
consistent Voigt layer averaging by exact flat-layer homogenization. The
source's separate thermal formulas and `NUYZ` label need correction.

Not established: a novel general materials result, experimental predictive
accuracy, processability at 70 vol% nested silica, GPU speedup, or a
submission-ready original research paper. The numerical audit itself is
complete and reproducible. The physical validation condition cannot be
completed without matched measured materials and microstructure data.
