# Hypotheses and analysis chronology

The initial working hypothesis, formed after the first literature pass and
before the principal parameter sweep, was that uniform-strain averaging would
materially overpredict through-thickness stiffness and shift CTE relative to
exact layer homogenization. This was **refuted** for the balanced cross-ply
under the supplied transversely isotropic ply model: EY and all principal
CTEs agree to roundoff when the Voigt thermal stress is mixed consistently.
The script's independently prescribed CTE formulas cause the larger error.

The second working hypothesis was that the particulate matrix closure would
have a larger impact than the layer closure, especially at high silica
fraction. Supported *within the two selected idealized models*: the 512-point
illustrative Sobol sweep gave median legacy-vs-Mori–Tanaka EX difference of
7.04%, versus median Voigt-vs-exact EX difference 0.0012%. This is not a
probabilistic statement about real manufactured material; ranges are
illustrative and no probability distribution was calibrated.

The CTE-error sign reversal was observed in the pilot grid and then tested
with one-dimensional root finding. It was **post-hoc**: 17/17 tested fiber
fractions had one root in the silica interval 0–0.7. No formal statistical
significance is claimed for deterministic model output.

Baselines: source-script algebra, thermally consistent Voigt layers, exact
flat layers with the same ply, and exact flat layers with bonded-sphere
Mori–Tanaka matrix. Variables: nested silica fraction 0–0.7, fiber fraction
0.2–0.6, and eight illustrative constituent inputs in `experiments.py`.
Metrics: elastic moduli (GPa), CTE (ppm/K), closure biases, positive
definiteness, numerical reference error, and derivative agreement. A
material claim would require measured constituent and laminate properties;
that criterion was not met.
