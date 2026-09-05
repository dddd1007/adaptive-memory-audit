#!/usr/bin/env bash
set -euo pipefail

project_dir="$(cd "$(dirname "$0")" && pwd)"
julia --startup-file=no "$project_dir/src/julia/run_revision_simulations.jl"
Rscript "$project_dir/src/r/summarize_and_validate.R" "$project_dir"

if [[ "${RUN_EMPIRICAL:-0}" == "1" ]]; then
  Rscript "$project_dir/empirical_workflow/run_empirical.R"
fi
