# Reproducibility status for public release v1.0.0

## Verified execution environment

- Julia 1.11.7 executed the formal simulation extensions. The project supports
  Julia 1.10 or later and uses no external Julia packages.
- R 4.5.3 executed the empirical analyses, sensitivity grid, summaries,
  figures, tables, and validation steps.
- `lme4` 2.0.6 fitted two crossed participant/card random-intercept logistic
  benchmarks. Both fits converged and were non-singular.
- Base and recommended R packages are used for the remaining R workflow.

## Included computational evidence

- 19,200 parameter-recovery estimator rows;
- 5,760 unlogged-context estimator rows;
- 12,000 legacy micro-randomized benchmark rows;
- 12,000 WCLS estimator rows across 48 design cells;
- 2,400 context-proxy estimator rows;
- 8,000 cognitive-state delay-recovery estimator rows;
- 800 state-model diagnostic datasets;
- 4,000 observation-model misspecification rows;
- 27 calibration-search candidates and 2,000 independent robustness rows;
- 24 empirical sensitivity cells with 1,200 participant-cluster bootstrap
  estimates per cell;
- five main and three supplementary figures in PNG and SVG;
- seven main and eight supplementary tables in CSV and Markdown.

The public empirical input contains 2,045 cross-day transitions and 199 Again
outcomes. Its eight release checks pass. The global validation record reports
38 of 38 checks passed.

## Public reproduction boundary

The public workflow starts from the minimized de-identified transition table.
The article's empirical estimators, sensitivity analyses, mixed-effects
benchmarks, figures, tables, and validation checks are reproducible from that
point.

Raw exports, event-level timing records, absolute dates, original identifiers,
linkage material, and the private event-to-transition preparation are excluded.
The private preparation was audited before release, but readers cannot rerun
that step from this repository. This boundary reduces longitudinal
re-identification risk and is documented in the data dictionary and public data
audit.

## Reproduction modes

`--mode=reference` rebuilds summaries and generated assets from the included
formal outputs. It does not rerun Monte Carlo simulations, the empirical
bootstrap, or mixed-effects models.

`--mode=full` reruns the Julia simulations and R empirical analyses before
rebuilding and validating all generated assets.

`--mode=smoke` runs one Julia replication per design cell and writes only to
`validation/smoke_outputs/`.

The current public release combines completed formal native runs with a
validated reference build. It does not claim that every component was captured
in a newly archived, single-command full-run console log.

## Interpretation of successful reproduction

Monte Carlo reruns use documented deterministic random streams. Differences in
operating system, numerical libraries, threading, and floating-point behavior
can prevent byte-identical CSV files. Scientific reproduction is defined by the
declared design, row counts, estimands, direction, bias, coverage, recovery
criteria, and successful validation checks.
