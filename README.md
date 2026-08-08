# Adaptive-memory audit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21850352.svg)](https://doi.org/10.5281/zenodo.21850352)
[![GitHub release](https://img.shields.io/github/v/release/dddd1007/adaptive-memory-audit)](https://github.com/dddd1007/adaptive-memory-audit/releases/tag/v1.0.0)

Public R and Julia reproducibility materials for:

> **Auditing psychological inferences about memory from adaptive learning logs: Cognitive-model recovery, policy logging, and randomized validation**

Authors: Muxue Zhou and Xiaokai Xia, School of Educational Science,
Hengyang Normal University.

This repository supports reader inspection and computational reproduction of
the article's simulations, empirical models, sensitivity analyses, figures,
tables, and validation checks. The manuscript, supplementary document,
submission correspondence, reviewer response, and document templates are not
distributed here.

## Scope

The project evaluates when delay effects can be interpreted in adaptive,
closed-loop learning systems. Julia implements the known-truth simulations and
model-recovery experiments. R analyzes the de-identified empirical transition
table, summarizes simulations, builds figures and tables, and runs the
machine-readable validation suite.

Included materials:

- R and Julia analysis source code;
- one minimized de-identified empirical analysis table;
- formal simulation and empirical result files;
- generated figures and tables;
- validation records and runtime information.

Excluded materials:

- the article and all submission documents;
- raw Anki exports and event-level logs;
- names, emails, telephone numbers, device identifiers, free text, source
  filenames, original event identifiers, and linkage keys;
- absolute dates, absolute timestamps, and exact within-day clock times;
- n-back, Stroop, FSRS-state, and participant-summary data that are outside the
  article's analyses.

## Data privacy boundary

The public empirical workflow begins at
[`data/analysis_transitions_deidentified.csv`](data/analysis_transitions_deidentified.csv).
It contains 2,045 cross-day transitions, 10 released fields, 12 pseudonymous
participant codes, and 100 pseudonymous card codes. The codes were created
during private de-identification; their crosswalk is absent from the repository
and release package. No card content is included.

The private event-to-transition preparation was audited before release, but raw
and event-level records are withheld to reduce longitudinal re-identification
risk. Consequently, readers can reproduce the statistical analyses that begin
with the released transition table, while the private de-identification and
event reconstruction step is outside the public computational boundary. See
[`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) and
[`data/PUBLIC_DATA_AUDIT.md`](data/PUBLIC_DATA_AUDIT.md).

## Requirements

The formal results were produced with:

- Julia 1.11.7; `Project.toml` supports Julia 1.10 or later and uses no external
  Julia packages;
- R 4.5.3;
- R package `lme4` 2.0.6 for two crossed participant/card random-intercept
  logistic models; base and recommended R packages are used elsewhere.

The tested dependency versions are also recorded in
[`R_DEPENDENCIES.csv`](R_DEPENDENCIES.csv).

Run commands from the repository root. `Rscript` is required for every mode;
`julia` is additionally required for full and smoke modes.

## Quick start

### Reference reproduction

```sh
Rscript run_all.R --mode=reference
```

This mode starts from the included formal replication-level outputs. It
regenerates summaries, five main and three supplementary figures, seven main and
eight supplementary tables, runtime manifests, and the 38-check validation
report. It does not rerun the longest Monte Carlo studies, empirical bootstrap,
or mixed-effects models.

### Full native rerun

```sh
Rscript run_all.R --mode=full
```

This mode reruns all Julia simulations and both R empirical analyses before
rebuilding and validating the generated assets. The empirical sensitivity grid
contains:

```text
4 planned-interval thresholds x 2 outcome definitions x 3 penalties
= 24 cells x 1,200 participant-cluster bootstraps
= 28,800 estimates
```

The full run also fits 24 leave-one-participant-out specifications and two
crossed random-intercept benchmarks. Runtime depends on available CPU resources.

### Smoke test

```sh
Rscript run_all.R --mode=smoke
```

Smoke mode runs one replication per Julia design cell and writes only to
`validation/smoke_outputs/`. It checks execution and output schemas without
overwriting the formal results.

## Repository map

| Path | Contents |
|---|---|
| `src/julia/` | Closed-loop, cognitive-state, calibration, WCLS, and context-proxy simulations |
| `src/R/` | Empirical analyses, summaries, builders, public-data checks, and global validation |
| `data/` | Minimized de-identified transition data, dictionary, and release audit |
| `outputs/` | Formal simulation and empirical outputs plus run manifests |
| `figures/` | Generated PNG and SVG figures |
| `tables/` | Full-precision CSV and display-rounded Markdown tables |
| `validation/` | Public-data, model-extension, empirical, runtime, and 38-check records |
| `run_all.R` | Main reproduction entry point |
| `Project.toml` | Julia project metadata and compatibility declaration |

## Validation

The project-level gate is implemented in
[`src/R/validate_outputs.R`](src/R/validate_outputs.R). The released reference
results pass all 38 checks. The checks cover formal row counts, known-truth
conclusions, WCLS/LPM numerical equivalence, context-proxy behavior,
observation-model boundaries, calibration robustness, the minimized empirical
schema, sensitivity results, generated asset inventories, and the R/Julia-only
source boundary.

Scientific reproduction is assessed against the declared designs, row counts,
estimands, direction, bias, coverage, recovery criteria, and validation checks.
Floating-point and platform differences can prevent byte-identical Monte Carlo
files across systems.

## Citation

| Item | Link |
|---|---|
| Repository | <https://github.com/dddd1007/adaptive-memory-audit> |
| GitHub release | <https://github.com/dddd1007/adaptive-memory-audit/releases/tag/v1.0.0> |
| Version DOI (v1.0.0) | <https://doi.org/10.5281/zenodo.21850352> |
| Concept DOI (all versions) | <https://doi.org/10.5281/zenodo.21850351> |
| Zenodo record | <https://zenodo.org/records/21850352> |

GitHub renders citation metadata from [`CITATION.cff`](CITATION.cff). Please
cite both this software package and the accompanying article when using the
materials.

Suggested software citation:

```text
Zhou, M., & Xia, X. (2026). Adaptive-memory audit: R and Julia reproducibility
materials (Version 1.0.0) [Computer software].
https://doi.org/10.5281/zenodo.21850352
```

Prefer the **version DOI** (`10.5281/zenodo.21850352`) when citing this exact
release. Use the **concept DOI** (`10.5281/zenodo.21850351`) only when you want
to cite the software as a whole across versions. After the article is published,
add the article DOI to `CITATION.cff`.

## Licenses

| Material | License |
|---|---|
| `run_all.R`, `Project.toml`, and `src/` | MIT License |
| `data/analysis_transitions_deidentified.csv` | Creative Commons Attribution 4.0 International |
| `outputs/`, `figures/`, `tables/`, validation records, and Markdown documentation | Creative Commons Attribution 4.0 International |

See [`LICENSE`](LICENSE), [`DATA_LICENSE.md`](DATA_LICENSE.md), and
[`LICENSES/CC-BY-4.0.txt`](LICENSES/CC-BY-4.0.txt). The software and materials
are provided without warranty. Users remain responsible for ethical and lawful
handling of human-participant-derived data.
