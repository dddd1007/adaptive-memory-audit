#!/usr/bin/env bash
# Validate an isolated copy; keep distributed results and frozen source unchanged.
set -euo pipefail
source_root="$(cd "$(dirname "$0")/.." && pwd)"
destination="${1:?Usage: run_public_check.sh NEW_OUTPUT_DIRECTORY}"
[[ ! -e "$destination" ]] || { echo 'Output directory already exists.' >&2; exit 2; }
mkdir -p "$destination"
destination="$(cd "$destination" && pwd)"
cp -R "$source_root" "$destination/research"
research="$destination/research"
mkdir -p "$destination/checks"
study="$research/studies/informative_observation"
julia_bin="${JULIA_BIN:-julia}"
rscript_bin="${RSCRIPT_BIN:-Rscript}"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
[[ "$("$julia_bin" --startup-file=no -e 'print(VERSION)')" == '1.12.7' ]] || { echo 'Julia 1.12.7 is required for this frozen study.' >&2; exit 2; }
"$julia_bin" --startup-file=no -e 'using InteractiveUtils; versioninfo()' > "$destination/checks/Julia_version.txt"
"$julia_bin" --startup-file=no "$research/validation/submission_checks/audit_stored_fits.jl" "$study" "$destination/checks"
"$julia_bin" --startup-file=no "$study/src/test_likelihood.jl"
"$julia_bin" --startup-file=no "$study/src/test_piecewise.jl"
"$rscript_bin" --vanilla "$research/validation/submission_checks/check_coverage_identity.R" "$destination/checks" "$study"
"$rscript_bin" --vanilla "$research/validation/validate_pearson_cr2.R" "$research" "$destination/checks/pearson_refits.csv"
"$rscript_bin" --vanilla "$research/validation/submission_checks/check_public_evidence.R" "$research" "$destination/checks"
echo 'Public reproduction checks passed. Stored Monte Carlo results were audited, not fully rerun.'
