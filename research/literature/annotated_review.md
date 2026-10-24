# Annotated review

The [2002 Chen–Yu–Feng–Lu paper](https://doi.org/10.1016/S0266-3538(02)00017-9)
derives exact thermoelastic formulas for perfectly bonded flat layers by
enforcing continuous tangential strain and normal traction. Our laminate
algorithm uses the same physical conditions. **It is an implementation and
audit, not a new homogenization method.**

[Whitcomb and Noh (2000)](https://doi.org/10.1177/002199830003400605)
derive 3D sublaminate homogenization and report finite-element verification.
This is the closest elastic precedent. Our constrained-energy calculation
verifies our code independently but does not improve their method.

[The 2017 glass fiber/polyimide/silica study](https://doi.org/10.1016/j.compositesa.2017.02.015)
uses a 3D unit-cell micromechanics model and an interphase for silica
nanoparticles. It directly defeats a claim that hierarchical glass-fiber
plus silica CTE modeling is new. Its polyimide and unidirectional architecture
are different from our illustrative epoxy cross-ply.

[Chaturvedi and Shen (1998)](https://doi.org/10.1016/S1359-6454(98)00089-5)
studied silica-filled epoxy encapsulant with particle unit cells and a
viscoelastic matrix. They report that CTE was comparatively insensitive to
matrix viscosity for their tested conditions, while local thermal stress was
more sensitive to arrangement and rate. This motivates separating effective
CTE from local stress claims.

[Ito, Suganuma and Wakashima (2000)](https://doi.org/10.1016/S0266-3538(00)00073-7)
measured glass/polypropylene angle-ply CTE and found a discrepancy from a
simple elastic laminate prediction, analyzed with matrix creep. This limits
transfer of elastic-only CTE predictions across temperature/time regimes.

[Jo and Lee (2017)](https://doi.org/10.1016/j.matpr.2017.06.113) report
measured silica/epoxy changes at 70 **weight** percent silica: about 51–55%
increase in Young's modulus and 25% CTE reduction relative to their own
baseline. This is not our 70 **volume** percent nested silica case, nor a
calibration for our assumed resin. The distinction matters greatly.

[Dittanet, Pearson and Kongkachuichay (2017)](https://doi.org/10.1016/j.ijadhadh.2017.06.006)
measured two epoxy systems with silica nanoparticles and compared Kerner CTE
and Halpin–Tsai stiffness; the abstract reports deviations at higher filler
weight fractions. [Ahmad et al. (2008)](https://doi.org/10.1016/j.compscitech.2007.07.015)
found particle shape affected mechanical properties. Neither supports treating
a spherical closure as validated at the illustrative 70 vol% loading.

[Baranau and Tallarek (2014)](https://doi.org/10.1039/C3SM52959B) analyze
monodisperse and polydisperse hard-sphere packing; about 0.64 is a common
random close packing value for monodisperse spheres. A nested silica fraction
of 0.70 therefore needs a polydisperse, ordered, or otherwise specified
geometry. It is retained as a stress-test input, not a proven manufacturable
monodisperse mixture.

[Kerner (1956)](https://doi.org/10.1088/0370-1301/69/8/305),
[Schapery (1968)](https://doi.org/10.1177/002199836800200308),
and [Hashin–Shtrikman (1963)](https://doi.org/10.1016/0022-5096(63)90060-7)
establish classical thermoelastic approximations and bounds. The source
script's formulas should be viewed in this established tradition, with their
particular assumptions stated.

[Recent FR4 measurement/modeling](https://doi.org/10.1016/j.jmrt.2026.01.011)
includes temperature-dependent thermal characterization of woven laminates.
It reinforces that our constant-temperature unidirectional cross-ply case is
not an experimental surrogate for woven FR4.

Only publisher abstracts/pages were available for many entries. The empty
`validated_findings` fields in `papers.csv` deliberately avoid inventing
paper-specific results from metadata.
