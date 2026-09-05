# Revision repair audit

## Scope

This repair replaces the provisional Python reimplementation with an
authoritative Julia simulation stack and an independent R validation layer.
The prior manuscript and archive candidate remain frozen for comparison.

## Material corrections

1. Reimplemented the revised observed-belief scheduler, tutor, context stress
   test, micro-randomized benchmark, and short-sequence model-recovery check in
   Julia.
2. Added a correctly specified dynamic likelihood study under random and
   observed-history adaptive selection to test the manuscript's ignorability
   claim directly.
3. Added the recorded tutor-belief estimator and completed the belief-versus-
   plan comparison under the latent-state policy stress test.
4. Standardized all conditional-coefficient and WCLS intervals to CR2
   covariance with coefficient-specific Satterthwaite degrees of freedom.
5. Increased focal observational and randomized cells to 1,000 replications,
   reported Monte Carlo standard-error diagnostics, and restored the 70%
   compliance condition.
6. Reduced the update-class diagnostic to four opportunities and retained it
   in the supplement because adaptive selection did not amplify the tested
   misspecification.
7. Rewrote the manuscript's theory, methods, results, discussion, tables, and
   computational provenance around the frozen Julia results.
8. Added the corresponding dynamic likelihood check for spaced retrieval and
   an observation-mechanism stress test that separates outcome omission from
   informative observation.
9. Re-executed the complete simulation stack in Julia 1.12.7, repeated focal
   fits in Julia 1.10.12 LTS, and regenerated the independent R validation with
   R 4.6.1 and `clubSandwich` 0.7.0.

## Validation finding about the earlier archive

The historical v1.1.0 replication CSV files do not reproduce from the Julia
source distributed beside them. Re-execution of the archived Julia focal cell
also produces lowercase Boolean fields, while the distributed CSV uses
capitalized Boolean fields. This is consistent with the earlier outputs having
passed through a different implementation or post-processing path. Those
historical numerical files are excluded from the repaired manuscript evidence
and must not be presented as a clean Julia rerun.

## Current validation

R independently refitted 51 frozen Julia-generated estimator instances. The
largest absolute differences were `6.11e-11` for coefficients, `1.04e-4` for
standard errors, and `2.99e-4` for Satterthwaite degrees of freedom. Focal
Julia 1.12.7 and Julia 1.10.12 LTS outputs differed by at most `1.42e-14`.

## Diagnostic of the earlier .910 coverage result

An exploratory Python run with 200 replications reported .910 coverage for the
100-learner scheduling cell. Its mean estimate (.01429), bias (.00012), and
mean model SE (.00442) were close to the final Julia values (.01428, .00012,
and .00444). The exploratory empirical SD was .00495, compared with .00452 in
the 1,000-replication Julia run, whose coverage was .946. Under nominal .95
coverage, 182 or fewer covered intervals out of 200 has an exact lower-tail
probability of about .012. The Python development code also used a simplified
degrees-of-freedom expression that was not algebraically identical to the
released CR2/Satterthwaite implementation. The discrepancy is therefore
consistent with Monte Carlo fluctuation plus a non-equivalent development
interval calculation; it does not indicate a change in the point estimand.

## Release consequence

The repaired archive should be published as a new Zenodo version after the
author reviews the code, outputs, and manuscript. It should supersede the
earlier v1.1.0 candidate while preserving the earlier record's version history.
