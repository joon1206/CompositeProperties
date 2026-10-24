# Literature search protocol (2026-10-24)

Scope: continuous glass fiber / silica-filled polymer thermoelasticity, exact
flat-layer homogenization, particulate micromechanics, packing, experimental
CTE and stiffness, and computational methods. Search primary publisher pages,
DOI records, and Crossref metadata. Inspect accessible abstracts or full text
before making scientific claims. Do not infer results from a citation title.

Search strategies: (1) architecture/material terms, (2) exact layer
homogenization and historical backward citation tracing, (3) particle-scale
models and measured silica/epoxy behavior, (4) adversarial queries for the
same methodology/result, including woven FR4 and temperature dependence.
Results and access level appear in `papers.csv`; 32 DOI records were resolved
via Crossref. `references.bib` is generated from those records. Publisher
abstract pages were inspected for the key works summarized below. The other
records are bibliographic context only until their content is checked.

Inclusion: relevant peer-reviewed primary methods or material measurements.
Exclusion: unverifiable titles, marketing property sheets without specimen
methods, generated material claims, and secondary sources as evidence for a
specific experimental result. No citation-count threshold was used. This is
targeted scoping, not a claim of an exhaustive systematic review.

Metadata refresh: `python scripts/build_literature.py`. The script obtains
metadata from Crossref and preserves a failure log; its hand-written notes
are limited to sources whose abstracts were inspected. Search dates refer to
the client date 2026-10-24.
