# Empirical sensitivity audit

Generated: 2026-08-08

## Public analysis input

The public workflow begins with 2,045 de-identified cross-day transitions and
10 released analysis fields. The table contains 12 pseudonymous participant
codes and 100 pseudonymous card codes. All 8 public-data checks passed. Raw exports, absolute
dates, exact within-day timestamps, original identifiers, linkage keys, and the
private de-identification step are outside this repository.

## Sensitivity design

The audit evaluated four planned-interval thresholds, two outcome definitions,
and three L2 penalties, yielding 24 cells. Every cell retained 1,200 successful
participant-cluster bootstrap estimates. The 0.25-, 0.50-, and 1-day thresholds
selected the same 1,724 transitions because the only observed sub-day plans
were 60, 330, and 600 seconds.

All cells used common participant-resample multiplicities from the independently
reserved sensitivity seed 20260807.

## Focal results

At the 0.5-day threshold and C = 10, the Again outcome gave OR = 0.759,
bootstrap 95% interval [0.348, 2.547]. The Again-or-Hard outcome gave OR =
0.781 [0.655, 2.205].

The leave-one-participant-out range was 0.589 to 1.342 for Again, with 11 of 12
estimates below 1, and 0.721 to 0.983 for Again-or-Hard, with 12 of 12 below 1.

The crossed random-intercept benchmarks gave OR = 0.811 [0.605, 1.089],
p = 0.164 for Again and OR = 0.802 [0.648, 0.992], p = 0.042 for
Again-or-Hard. Both models converged and were non-singular.

## Interpretation boundary

These observational estimates are association diagnostics. Adaptive
scheduling, unlogged learner state, the 12-cluster sample, and dependence on the
outcome definition prevent a causal spacing-effect interpretation.

## Validation

Targeted validation status: 13 of 13 checks passed.
