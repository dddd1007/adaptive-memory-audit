# Public de-identified analysis data

## File and analytical grain

`analysis_transitions_deidentified.csv` contains 2,045 eligible cross-day
transitions. One row represents a link from the final scheduler event for a
participant-card pair on one observed study day to the first rated retrieval
for the same pair on its next observed study day.

The table contains 12 pseudonymous participant codes and 100 pseudonymous card
codes created during private de-identification. Original identifiers, linkage
tables, card content, raw exports, absolute dates, and within-day timestamps are
absent.

## Fields

| Field | Type | Definition |
|---|---|---|
| `participant_code` | string | Pseudonymous learner label, P01 to P12 |
| `card_code` | string | Pseudonymous card label, C001 to C100; card content is withheld |
| `study_day` | integer | Relative calendar study day, ranging from 1 to 13 for eligible transitions |
| `daily_review_index` | integer | Order of the observed study day within the participant-card history, ranging from 1 to 8 |
| `first_ease` | integer | First rating on the current observed day: 1 Again, 2 Hard, 3 Good, 4 Easy |
| `last_ease` | integer | Final rating on the current observed day using the same 1-to-4 scale |
| `last_scheduled_interval_days` | numeric | Interval recorded by the final scheduler event on the current day, expressed in days |
| `next_gap_days` | numeric | Enacted elapsed time from that final event to the next observed day's first rated retrieval, expressed in days |
| `next_first_ease` | integer | Rating at the linked next observed day's first retrieval, using the same 1-to-4 scale |
| `scheduler_encoding` | string | Coarse scheduler classification: `fsrs_like`, `legacy_like`, or `mixed` |

The primary Again outcome is derived as
`as.integer(next_first_ease == 1)`. The broader sensitivity outcome is derived
as `as.integer(next_first_ease <= 2)`.

## Completeness and ranges

- 2,045 rows and 10 fields;
- no missing values in released fields;
- 199 next-retrieval Again outcomes;
- relative study days 1 to 13;
- planned intervals 0.000694 to 38 days;
- enacted gaps 0.004211 to 11.619250 days;
- scheduler groups: 407 `fsrs_like`, 708 `legacy_like`, and 930 `mixed`
  transitions.

## Derivation and privacy boundary

The public table is a minimized export from a private preprocessing workflow.
The private workflow ordered review events within participant and card, grouped
them by local study day, carried the final scheduler decision forward, and
linked it to the next observed day's first rated retrieval. That transformation
was independently audited before public release.

The raw event stream and de-identification program are intentionally outside
the repository because they require private source files and contain richer
longitudinal timing information. The released workflow therefore reproduces
all empirical models from the analysis-ready transition table, while the raw
event reconstruction itself remains outside the public reproduction boundary.

Relative gap precision is retained because delay is the focal analytical
variable. The table still represents a small longitudinal sample, so users
should avoid re-identification attempts or linkage to outside records.

## License and citation

The table is licensed under CC BY 4.0 as described in
[`../DATA_LICENSE.md`](../DATA_LICENSE.md). Cite the repository and the
accompanying article when reusing it.
