#!/usr/bin/env bash
# Reproduce in a separate directory; keep the distributed archive unchanged.
set -euo pipefail
script_dir="$(cd "$(dirname "$0")" && pwd)"
if [[ -d "$script_dir/../research/studies/informative_observation" ]]; then
  study_source="$script_dir/../research/studies/informative_observation"
else
  study_source="$(cd "$script_dir/../studies/informative_observation" && pwd)"
fi
mode="${1:-check}"
destination="${2:?Usage: reproduce_new.sh check|fresh NEW_OUTPUT_DIRECTORY}"
[[ "$mode" == check || "$mode" == fresh ]] || { echo 'Mode must be check or fresh'; exit 2; }
[[ ! -e "$destination" ]] || { echo 'Output directory already exists; choose a new directory.'; exit 2; }
mkdir -p "$destination"
cp -R "$study_source" "$destination/study"
study="$(cd "$destination/study" && pwd)"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
if command -v juliaup >/dev/null 2>&1; then
  julia_cmd=(julia +1.12.7 --startup-file=no)
else
  [[ "$(julia --startup-file=no -e 'print(VERSION)')" == 1.12.7 ]] || { echo 'Julia 1.12.7 required'; exit 2; }
  julia_cmd=(julia --startup-file=no)
fi
"${julia_cmd[@]}" "$study/src/test_likelihood.jl"
"${julia_cmd[@]}" "$study/src/test_piecewise.jl"
if [[ "$mode" == fresh ]]; then
  mv "$study/outputs/confirmatory" "$study/outputs/confirmatory_distributed"
fi
"${julia_cmd[@]}" "$study/src/run_study.jl" confirmatory 4
"${julia_cmd[@]}" "$study/src/summarize.jl"
Rscript --vanilla "$study/validation/summarize_and_plot.R" "$study"
if [[ "$mode" == fresh ]]; then
  "${julia_cmd[@]}" "$study/src/export_validation.jl" confirmatory
  Rscript --vanilla "$study/validation/independent_likelihood.R" "$study" "$study/validation/confirmatory_julia_refits.csv" "$study/validation/confirmatory_r_comparison.csv"
  "${julia_cmd[@]}" "$study/src/check_quadrature.jl"
fi
