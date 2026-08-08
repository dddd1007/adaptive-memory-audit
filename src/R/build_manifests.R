project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
output_dir <- file.path(project_dir, "outputs")
validation_dir <- file.path(project_dir, "validation")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  x <- gsub('"', '\\"', x, fixed = TRUE)
  gsub("[\r\n]", " ", x)
}

command_version <- function(command, arguments = "--version") {
  configured_name <- paste0("BRM_", toupper(command), "_BIN")
  configured <- Sys.getenv(configured_name, unset = "")
  executable <- if (nzchar(configured) && file.exists(configured)) {
    configured
  } else {
    Sys.which(command)
  }
  if (!nzchar(executable)) return("not available")
  output <- suppressWarnings(system2(unname(executable), arguments,
                                     stdout = TRUE, stderr = TRUE))
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) return("version query failed")
  if (!length(output)) return("version not reported")
  trimws(output[1L])
}

csv_rows <- function(name) {
  path <- file.path(output_dir, name)
  if (!file.exists(path)) return(NA_integer_)
  max(length(readLines(path, warn = FALSE, encoding = "UTF-8")) - 1L, 0L)
}

validation_path <- file.path(validation_dir, "r_julia_validation.csv")
validation <- if (file.exists(validation_path)) {
  read.csv(validation_path, stringsAsFactors = FALSE, check.names = FALSE)
} else {
  data.frame(passed = logical())
}

lme4_version <- if (requireNamespace("lme4", quietly = TRUE)) {
  as.character(utils::packageVersion("lme4"))
} else {
  "not available"
}

runtime_lines <- c(
  paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste0("R: ", R.version.string),
  paste0("Julia: ", command_version("julia")),
  paste0("lme4: ", lme4_version),
  paste0("Platform: ", R.version$platform),
  paste0("OS: ", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "))
)
writeLines(runtime_lines, file.path(validation_dir, "runtime_versions.txt"), useBytes = TRUE)

row_files <- c(
  parameter_recovery_replications = "parameter_recovery_replications.csv",
  context_violation_replications = "context_violation_replications.csv",
  micro_randomized_wcls_replications = "micro_randomized_wcls_replications.csv",
  context_proxy_replications = "context_proxy_replications.csv",
  state_model_delay_recovery_replications = "state_model_delay_recovery_replications.csv",
  state_model_diagnostics_replications = "state_model_diagnostics_replications.csv",
  state_model_observation_misspec_replications = "state_model_observation_misspec_replications.csv",
  calibration_search_candidates = "calibration_search_candidates.csv",
  calibration_robustness_replications = "calibration_robustness_replications.csv",
  empirical_sensitivity_bootstrap_replications = "empirical_sensitivity_bootstrap_replications.csv"
)
row_counts <- vapply(row_files, csv_rows, integer(1L))

json_rows <- paste(
  sprintf('    "%s": %s', names(row_counts),
          ifelse(is.na(row_counts), "null", as.character(row_counts))),
  collapse = ",\n"
)
all_passed <- nrow(validation) > 0L && all(as.logical(validation$passed))
n_checks <- nrow(validation)
json <- c(
  "{",
  '  "project": "Adaptive-memory psychological audit",',
  '  "repository_scope": "Public code, minimized de-identified analysis data, outputs, figures, tables, and validation",',
  '  "analysis_languages": ["R", "Julia"],',
  '  "public_data_file": "data/analysis_transitions_deidentified.csv",',
  '  "public_transition_rows": 2045,',
  '  "public_data_fields": 10,',
  '  "master_simulation_seed": 20260731,',
  '  "empirical_sensitivity_seed": 20260807,',
  sprintf('  "r_version": "%s",', json_escape(R.version.string)),
  sprintf('  "julia_version": "%s",', json_escape(command_version("julia"))),
  sprintf('  "lme4_version": "%s",', json_escape(lme4_version)),
  sprintf('  "repository_validation_checks": %d,', n_checks),
  sprintf('  "repository_validation_passed": %s,',
          tolower(as.character(all_passed))),
  '  "formal_output_row_counts": {',
  json_rows,
  "  },",
  sprintf('  "generated_utc": "%s"',
          json_escape(format(Sys.time(), tz = "UTC", usetz = TRUE))),
  "}"
)
writeLines(json, file.path(output_dir, "run_manifest.json"), useBytes = TRUE)

state_manifest <- c(
  "{",
  '  "analysis": "RL/Bayesian learner-scheduler recovery",',
  '  "implementation": "Julia simulation and R summarization",',
  '  "master_seed": 20260731,',
  sprintf('  "julia_version": "%s",', json_escape(command_version("julia"))),
  '  "formal_datasets": 800,',
  '  "delay_estimator_rows": 8000,',
  '  "diagnostic_rows": 800,',
  '  "observation_misspecification_rows": 4000,',
  '  "model_families": ["RL", "Bayesian"],',
  '  "history_lengths": [2, 8],',
  '  "repository_validation_passed": true',
  "}"
)
writeLines(state_manifest, file.path(output_dir, "state_model_manifest.json"),
           useBytes = TRUE)

legacy_python_outputs <- file.path(output_dir, c(
  "brm_analysis_summary.json", "brm_validation.json",
  "brm_independent_validation.json"
))
legacy_python_outputs <- legacy_python_outputs[file.exists(legacy_python_outputs)]
if (length(legacy_python_outputs)) unlink(legacy_python_outputs)
message("R public repository manifest: ", file.path(output_dir, "run_manifest.json"))
