# Four explicit novelty reviews

1. **Broad mapping.** Classical Schapery/Kerner/Halpin–Tsai, Hashin–Shtrikman
bounds, exact flat-layer theory, and particle/fiber RVEs already cover the
major algorithmic pieces. See `annotated_review.md` and 32 DOI records.
2. **Adversarial before implementation.**
[Chen et al. (2002)](https://doi.org/10.1016/S0266-3538(02)00017-9)
already give exact thermoelastic layer formulas;
[Whitcomb and Noh (2000)](https://doi.org/10.1177/002199830003400605)
give 3D sublaminate formulas; the
[2017 glass/polyimide/silica work](https://doi.org/10.1016/j.compositesa.2017.02.015)
already treats the hierarchical material combination.
3. **Post-implementation.** The main observed result is that balanced
cross-ply symmetries make consistent Voigt estimates close to exact, while
the script's CTE formulas and matrix closure matter more. This is a
case-specific consequence of established equations, not a new physical law.
The sign reversal in CTE error is likewise an error pattern of one script.
4. **Pre-manuscript.** Experimental silica/epoxy papers report dependence on
resin, filler morphology, temperature, and high loading; there is no matched
validation set for this formulation. The current work is suitable as a
reproducible technical audit and research note. **A submission-ready original
materials study has not been established.**

Closest competitors and distinctions are in `related_work_matrix.csv`. No
novel method or experimentally confirmed physical effect is claimed.
