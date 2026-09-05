# Closed-loop audit computational study

This directory contains the computational study accompanying the September 2026 revision on target-specific ignorability in adaptive learning logs. It is separate from the earlier repository-root study and its published release metadata. No new DOI or GitHub Release is implied.

Julia implements simulation and estimation. R independently checks likelihood fits and summaries and refits the regression benchmarks. The primary numerical findings include conditional-coefficient recovery across spaced retrieval and tutoring, randomized proximal assignment effects, and a matched-rate informative-observation experiment with uncertain initial state. The 12-learner empirical case is descriptive.

## Evidence and version boundaries

- `outputs_julia/` contains retained earlier Monte Carlo outputs. These were subjected to current focal checks and independent regression refits; they were not all freshly simulated on the current machine.
- `studies/informative_observation/outputs/confirmatory/` contains the final 3 × 1,000-replication response-only likelihood study. All 3,000 original seeds were uniformly refitted with clipping-aware piecewise Gauss-Legendre integration. Valid fits number 994, 990 and 987; all other fits retain boundary status and count as noncoverage in unconditional summaries.
- Earlier GH21 results and numerical diagnostics remain explicitly named for provenance. They do not supply the final reported coefficients or intervals.
- Known-observation-mechanism joint fits are fixed-sample numerical references, not another complete 3,000-replication performance study or a proof that unknown observation mechanisms are identified.
- `validation/submission_checks/` adds portable read-only audits of stored fits, coverage denominators and a deterministic logistic non-collapsibility identity. The numerical algorithms and frozen design are unchanged.
- `validation/CR2_COORDINATE_AUDIT.md` describes the Pearson linearization used for independent CR2 agreement. Default GLM adjustment coordinates need not give identical standard errors.

The `planned_not_implemented` string in the original `design.toml` is retained as part of the hashed design specification. Actual completion and numerical refinement are documented by the frozen run manifest, `execution.toml`, replication files, and `validation/DESIGN_AND_RUN_NOTES.md`. Historical validation reports describe their original runs; they are not current blanket certification.

## Requirements

Use Julia 1.12.7 for exact frozen Xoshiro stream reproduction. The Julia study uses standard libraries only. Independent validation was run with native R 4.6.1 and clubSandwich 0.7.0; package versions are listed in `R_DEPENDENCIES_CURRENT.csv`. Base R suffices for the new likelihood and summary checks. Install platform-native tools and keep libraries/caches outside this repository. No GPU is needed.

## Validate without changing distributed results

From this directory, run:

```bash
bash tools/run_public_check.sh /path/to/new-check-directory
```

The destination must not exist. This command copies the study before running checks. It verifies frozen source and configuration hashes, all 3,000 stored fit identities and diagnostic flags, likelihood and integration tests, the coverage denominator identity, 51 independent Pearson-coordinate CR2 refits, and the minimized public-data checks. It does not claim to rerun the full Monte Carlo study. `JULIA_BIN` and `RSCRIPT_BIN` may specify executable paths; otherwise the tools on PATH are used. Julia version must be 1.12.7.

To repeat all informative-observation Monte Carlo fits in a separate directory using the original design and seeds:

```bash
bash tools/reproduce_new.sh fresh /path/to/new-simulation-directory
```

That separate entry point also reruns independent R fixed-data likelihood checks and the quadrature audit. Its `check` mode verifies and retains matching formal checkpoints rather than refitting them. Do not recalibrate the frozen intercepts to seek a preferred result. Failures are retained.

The earlier full-study entry point is `run_all.sh`. Use it only in a separate copy because it regenerates old outputs and retains historical numerical conventions. It is not necessary for the current read-only checks. The empirical ridge implementation is in `empirical_workflow/src/R/`; earlier mixed-effects artifacts remain historical and do not support causal inference for the 12-learner case.

## Public data and licenses

The only participant-level analysis input is the already minimized 2,045-row table in `empirical_workflow/data/analysis_transitions_deidentified.csv`. No raw exports, original identifiers, linkage keys, exact timestamps, learning-item content or manuscript/submission files are included. Synthetic latent-state validation datasets are simulated, not participant disclosures. Existing license and attribution files are preserved; see the repository-root licenses and the empirical workflow licenses.

`SOURCE_SYNC_MANIFEST.csv` records byte-level correspondence with the local research working tree. The repository-root SHA256 manifest covers distributed files. Local manuscript editing, typesetting, correspondence and document-generation scripts remain outside the public research package.
