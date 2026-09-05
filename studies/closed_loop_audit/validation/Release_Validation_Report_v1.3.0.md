# Release validation report, v1.3.0 candidate

Date: 2026-09-03

## Result

All release checks passed.

## Computational checks

- The complete Monte Carlo workflow ran from the archived Julia source under
  Julia 1.12.7.
- Primary observational and randomized cells contain 1,000 replications;
  secondary cells contain the declared 400 or 500 replications.
- The dynamic likelihood check covers spaced retrieval and adaptive tutoring
  under random and observed-history adaptive selection.
- The observation-mechanism stress test contains 400 replications in each of
  five declared conditions.
- R 4.6.1 with `clubSandwich` 0.7.0 independently refitted 51 frozen datasets.
  Maximum absolute differences were `6.11e-11` for coefficients, `1.04e-4`
  for CR2 standard errors, and `2.99e-4` for Satterthwaite degrees of freedom.
- Focal Julia 1.12.7 and Julia 1.10.12 LTS results differed by at most
  `1.42e-14`.

## Manuscript checks

- Abstract: 111 words, within the JEBS 100 to 120 word requirement.
- Keywords: five.
- Anonymous main manuscript: 32 US Letter pages, five tables, four figures,
  and 22 native Word display-equation objects.
- Anonymous supplement: 12 US Letter pages and nine tables.
- Both documents use 1-inch margins, 12-point body text, and double spacing.
- DOCX ZIP integrity, heading structure, image placement, and accessibility
  checks passed.
- Final rendered pages are unchanged by the metadata privacy scrub.
- Scans found no author names, institutional identifiers, email addresses,
  manuscript number, or project DOI in either anonymous file.

## Analysis boundary

Julia is the authoritative implementation for data generation, estimators, and
Monte Carlo results. R provides independent refits, summaries, and figures.
Python is not part of the research-analysis dependency graph; document QA tools
do not generate or alter reported numerical results.
