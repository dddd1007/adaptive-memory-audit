# Public data release audit

Audit date: 2026-08-08

## Released dataset

- File: `analysis_transitions_deidentified.csv`
- Grain: eligible cross-day participant-card transition
- Rows: 2,045
- Released fields: 10
- Pseudonymous participant codes: 12
- Pseudonymous card codes: 100
- Again outcomes: 199

## De-identification and minimization checks

The public export contains no names, email addresses, telephone numbers,
postal addresses, device identifiers, free text, source filenames, original
event IDs, original participant IDs, original card IDs, linkage keys, absolute
dates, absolute timestamps, or within-day clock times.

Only fields used by the article's empirical models were retained. The public
export removed event-level response times, raw scheduler encodings, exact local
clock time, redundant outcome fields, rows without a linked future retrieval,
and all n-back, Stroop, FSRS-state, and participant-summary datasets.

Participant and card codes were created during private de-identification. The
crosswalk to original identifiers is absent from the repository and release
package.

## Automated quality checks

The release validation checks:

1. expected row and field counts;
2. participant and card code formats and cardinalities;
3. absence of direct-identifier and absolute-time field names;
4. uniqueness of the participant-card-study-day key;
5. completeness of all released fields;
6. allowed rating values and relative study-day range;
7. finite positive planned and enacted delays;
8. reproduction of the 199 Again outcomes.

All eight checks passed for this release. The machine-readable result is stored
in `validation/public_data_validation.csv`.

## Residual risk and analytical boundary

The data contain 12 longitudinal learner histories and exact relative delays,
which retain some residual re-identification risk. Those variables are needed
for participant-cluster bootstrap, crossed random effects, scheduler
conditioning, and the focal delay models. The release omits the richer event
stream and unrelated cognitive-task records that would increase linkability.

The raw-to-transition preprocessing audit remains a private provenance step.
The public code begins with the released transition table and reproduces the
article's statistical analyses from that point forward.
