arguments_full <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", arguments_full, value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[1L]) else "run_all.R"
project_dir <- normalizePath(dirname(script_path), mustWork = TRUE)
Sys.setenv(BRM_PROJECT_DIR = project_dir)

arguments <- commandArgs(trailingOnly = TRUE)
mode_argument <- grep("^--mode=", arguments, value = TRUE)
mode <- if (length(mode_argument)) sub("^--mode=", "", mode_argument[1L]) else "reference"
valid_modes <- c("reference", "full", "smoke")
if (!mode %in% valid_modes) stop("Unknown mode: ", mode, ". Choose: ", paste(valid_modes, collapse = ", "))

run_julia <- function(script) {
  julia <- Sys.which("julia")
  if (!nzchar(julia)) stop("Julia is unavailable. Install Julia 1.10 or later for full/smoke modes.")
  status <- system2(unname(julia), c(paste0("--project=", project_dir), script))
  if (!identical(status, 0L)) stop("Julia script failed: ", script)
}

run_r <- function(script) {
  source(file.path(project_dir, "src", "R", script), echo = FALSE, chdir = FALSE)
}

if (mode == "smoke") {
  smoke_dir <- file.path(project_dir, "validation", "smoke_outputs")
  if (!dir.exists(smoke_dir)) dir.create(smoke_dir, recursive = TRUE)
  Sys.setenv(
    BRM_OUTPUT_DIR = smoke_dir,
    BRM_RECOVERY_REPS = "1",
    BRM_CONTEXT_REPS = "1",
    BRM_MRT_REPS = "1",
    BRM_WCLS_REPS = "1",
    BRM_PROXY_REPS = "1",
    BRM_STATE_REPS = "1",
    BRM_STATE_OBS_REPS = "1",
    BRM_CALIBRATION_SEARCH_REPS = "1",
    BRM_CALIBRATION_ROBUST_REPS = "1"
  )
  run_julia(file.path(project_dir, "src", "julia", "run_closed_loop_simulation.jl"))
  run_julia(file.path(project_dir, "src", "julia", "run_state_model_simulation.jl"))
  run_julia(file.path(project_dir, "src", "julia", "run_calibration_sensitivity.jl"))
  message("Smoke simulation completed in ", smoke_dir)
  quit(save = "no", status = 0L)
}

if (mode == "full") {
  Sys.unsetenv(c("BRM_OUTPUT_DIR", "BRM_RECOVERY_REPS", "BRM_CONTEXT_REPS",
                 "BRM_MRT_REPS", "BRM_WCLS_REPS", "BRM_PROXY_REPS",
                 "BRM_EXTENSIONS_ONLY", "BRM_STATE_REPS", "BRM_STATE_OBS_REPS",
                 "BRM_STATE_OBS_ONLY", "BRM_STATE_OBS_SEED_OFFSET",
                 "BRM_CALIBRATION_SEARCH_REPS", "BRM_CALIBRATION_ROBUST_REPS",
                 "BRM_CALIBRATION_SEARCH_SEED_OFFSET",
                 "BRM_CALIBRATION_ROBUST_SEED_OFFSET",
                 "BRM_EMPIRICAL_SENSITIVITY_REPS"))
  run_julia(file.path(project_dir, "src", "julia", "run_closed_loop_simulation.jl"))
  run_julia(file.path(project_dir, "src", "julia", "run_state_model_simulation.jl"))
  run_julia(file.path(project_dir, "src", "julia", "run_calibration_sensitivity.jl"))
  run_r("analyse_empirical_case.R")
  run_r("audit_empirical_sensitivity.R")
}

if (mode %in% c("reference", "full")) {
  run_r("summarise_simulations.R")
  run_r("validate_model_extensions.R")
  run_r("build_figures.R")
  run_r("build_tables.R")
  run_r("check_r_sources.R")
  run_r("validate_outputs.R")
  run_r("build_manifests.R")
}
message("Public R/Julia reproducibility workflow completed in mode: ", mode)
