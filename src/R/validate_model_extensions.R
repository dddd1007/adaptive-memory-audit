project_dir <- Sys.getenv(
  "BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE)
)
output_dir <- file.path(project_dir, "outputs")
validation_dir <- file.path(project_dir, "validation")
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)

read_output <- function(filename) {
  read.csv(file.path(output_dir, filename), stringsAsFactors = FALSE,
           check.names = FALSE)
}

checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    observed = as.character(observed),
    expected = as.character(expected),
    stringsAsFactors = FALSE
  )
}

observation <- read_output("state_model_observation_misspec_replications.csv")
observation_summary <- read_output("state_model_observation_misspec_summary.csv")
calibration_search <- read_output("calibration_search_candidates.csv")
calibration <- read_output("calibration_robustness_replications.csv")
calibration_summary <- read_output("calibration_robustness_summary.csv")

add_check("observation raw row count", nrow(observation) == 4000L,
          nrow(observation), 4000L)
add_check("observation summary row count", nrow(observation_summary) == 40L,
          nrow(observation_summary), 40L)
add_check("calibration search candidate count", nrow(calibration_search) == 27L,
          nrow(calibration_search), 27L)
add_check("calibration robustness raw row count", nrow(calibration) == 2000L,
          nrow(calibration), 2000L)
add_check("calibration robustness summary row count",
          nrow(calibration_summary) == 20L, nrow(calibration_summary), 20L)

observation_cells <- unique(observation[c(
  "learner_family", "scheduler_family", "observation_condition",
  "candidate_family"
)])
add_check("observation factorial cells", nrow(observation_cells) == 40L,
          nrow(observation_cells), "2 learner × 2 scheduler × 5 equation × 2 candidate")

rank_one <- calibration_search[calibration_search$rank == 1L, , drop = FALSE]
rank_one_ok <- nrow(rank_one) == 1L &&
  abs(rank_one$adaptivity - 0.75) < 1e-12 &&
  abs(rank_one$natural_deviation_sd - 0.50) < 1e-12 &&
  abs(rank_one$outcome_intercept + 2.40) < 1e-12
add_check("rank-one calibration is original generator", rank_one_ok,
          if (nrow(rank_one)) {
            paste(rank_one$adaptivity, rank_one$natural_deviation_sd,
                  rank_one$outcome_intercept, sep = "/")
          } else "missing", ".75/.50/-2.40")

recomputed_loss <-
  ((calibration_search$mean_failure_rate -
      calibration_search$target_failure_rate) / 0.015)^2 +
  ((calibration_search$mean_scheduled_actual_spearman -
      calibration_search$target_scheduled_actual_spearman) / 0.05)^2 +
  ((calibration_search$mean_naive_log_odds -
      calibration_search$target_naive_log_odds) / 0.10)^2
loss_error <- max(abs(recomputed_loss - calibration_search$loss))
add_check("calibration loss recomputation", loss_error < 1e-12,
          format(loss_error, scientific = TRUE), "maximum error < 1e-12")

candidate_rl <- observation[observation$candidate_family == "RL", , drop = FALSE]
candidate_rl$recovered <- tolower(candidate_rl$family_recovered) == "true"
recovery <- aggregate(
  recovered ~ learner_family + observation_condition,
  candidate_rl, mean
)
get_recovery <- function(family, condition) {
  recovery$recovered[
    recovery$learner_family == family &
      recovery$observation_condition == condition
  ]
}
add_check("correct-equation family recovery",
          get_recovery("RL", "correct") >= 0.95 &&
            get_recovery("Bayesian", "correct") >= 0.98,
          paste(get_recovery("RL", "correct"),
                get_recovery("Bayesian", "correct"), sep = "/"),
          "RL >= .95 and Bayesian >= .98")
add_check("combined misspecification selects RL",
          get_recovery("RL", "combined_plus_0.25_slope_0.80") == 1 &&
            get_recovery("Bayesian", "combined_plus_0.25_slope_0.80") == 0,
          paste(
            get_recovery("RL", "combined_plus_0.25_slope_0.80"),
            get_recovery("Bayesian", "combined_plus_0.25_slope_0.80"),
            sep = "/"
          ), "RL = 1 and Bayesian = 0")

selected <- observation[
  observation$candidate_family == observation$selected_family, , drop = FALSE
]
selected_rmse <- aggregate(
  state_rmse ~ learner_family + observation_condition,
  selected, mean
)
rmse_value <- function(family, condition) {
  selected_rmse$state_rmse[
    selected_rmse$learner_family == family &
      selected_rmse$observation_condition == condition
  ]
}
rmse_degrades <- all(vapply(c("RL", "Bayesian"), function(family) {
  rmse_value(family, "combined_plus_0.25_slope_0.80") >
    rmse_value(family, "correct")
}, logical(1L)))
add_check("selected-state RMSE degrades under combined misspecification",
          rmse_degrades,
          paste(
            round(rmse_value("RL", "correct"), 3),
            round(rmse_value("RL", "combined_plus_0.25_slope_0.80"), 3),
            round(rmse_value("Bayesian", "correct"), 3),
            round(rmse_value("Bayesian", "combined_plus_0.25_slope_0.80"), 3),
            sep = "/"
          ), "combined RMSE > correct RMSE for both families")

schedule <- calibration_summary[
  calibration_summary$estimator == "schedule_conditioned", , drop = FALSE
]
max_schedule_bias <- max(abs(schedule$bias))
add_check("top-five schedule-conditioned robustness",
          max_schedule_bias < 0.02, round(max_schedule_bias, 6),
          "maximum absolute bias < .02")

naive <- calibration_summary[
  calibration_summary$estimator == "naive", , drop = FALSE
]
minimum_naive_sign_error <- min(naive$sign_error_rate)
add_check("top-five naive sign reversal",
          minimum_naive_sign_error >= 0.98,
          round(minimum_naive_sign_error, 3),
          "minimum sign-error rate >= .98")

result <- do.call(rbind, checks)
path <- file.path(validation_dir, "model_calibration_extension_validation.csv")
write.csv(result, path, row.names = FALSE)
if (!all(result$passed)) {
  failed <- paste(result$check[!result$passed], collapse = "; ")
  stop("Model/calibration extension validation failed: ", failed)
}

json_number <- function(x) {
  format(as.numeric(x), digits = 16, scientific = FALSE, trim = TRUE)
}
json_numeric_array <- function(x) {
  paste0("[", paste(vapply(x, json_number, character(1L)), collapse = ", "), "]")
}
condition_order <- c(
  "correct", "intercept_minus_0.25", "intercept_plus_0.25",
  "state_slope_0.80", "combined_plus_0.25_slope_0.80"
)
family_recovery_vector <- function(family) {
  vapply(condition_order, function(condition) {
    get_recovery(family, condition)
  }, numeric(1L))
}
true_family <- observation[
  observation$candidate_family == observation$learner_family, , drop = FALSE
]
true_family_rmse <- aggregate(
  state_rmse ~ learner_family + observation_condition,
  true_family, mean
)
true_family_rmse_vector <- function(family) {
  vapply(condition_order, function(condition) {
    true_family_rmse$state_rmse[
      true_family_rmse$learner_family == family &
        true_family_rmse$observation_condition == condition
    ]
  }, numeric(1L))
}
manifest_path <- file.path(
  output_dir, "model_calibration_extension_manifest.json"
)
manifest_lines <- c(
  "{",
  '  "analysis_date": "2026-08-08",',
  '  "implementation": {"simulation": "Julia", "summaries_and_validation": "R"},',
  '  "observation_misspecification": {',
  '    "replications_per_learner_scheduler_cell": 100,',
  '    "n_stages": 8,',
  '    "candidate_families": ["RL", "Bayesian"],',
  '    "conditions": ["correct", "intercept_minus_0.25", "intercept_plus_0.25", "state_slope_0.80", "combined_plus_0.25_slope_0.80"],',
  '    "seed_base_offset": 130000000,',
  '    "learner_family_offset": 10000000,',
  '    "scheduler_family_offset": 1000000,',
  '    "stage_multiplier": 10000,',
  paste0('    "family_recovery_RL": ',
         json_numeric_array(family_recovery_vector("RL")), ','),
  paste0('    "family_recovery_Bayesian": ',
         json_numeric_array(family_recovery_vector("Bayesian")), ','),
  paste0('    "true_family_state_rmse_RL": ',
         json_numeric_array(true_family_rmse_vector("RL")), ','),
  paste0('    "true_family_state_rmse_Bayesian": ',
         json_numeric_array(true_family_rmse_vector("Bayesian"))),
  "  },",
  '  "calibration_search": {',
  '    "adaptivity_grid": [0.50, 0.75, 1.00],',
  '    "natural_deviation_sd_grid": [0.35, 0.50, 0.65],',
  '    "outcome_intercept_grid": [-2.65, -2.40, -2.15],',
  '    "replications_per_candidate": 40,',
  '    "loss_targets": [0.09731051344743276, 0.8146452768576546, -0.2158255188692115],',
  '    "loss_scales": [0.015, 0.05, 0.10],',
  '    "seed_offset": 160000000,',
  paste0('    "rank_one": {"adaptivity": ', json_number(rank_one$adaptivity),
         ', "natural_deviation_sd": ', json_number(rank_one$natural_deviation_sd),
         ', "outcome_intercept": ', json_number(rank_one$outcome_intercept),
         ', "mean_failure_rate": ', json_number(rank_one$mean_failure_rate),
         ', "mean_scheduled_actual_spearman": ',
         json_number(rank_one$mean_scheduled_actual_spearman),
         ', "mean_naive_log_odds": ', json_number(rank_one$mean_naive_log_odds),
         ', "loss": ', json_number(rank_one$loss), '}'),
  "  },",
  '  "calibration_robustness": {',
  '    "top_candidates": 5,',
  '    "replications_per_candidate": 100,',
  '    "estimators_per_dataset": 4,',
  '    "seed_offset": 170000000,',
  paste0('    "maximum_absolute_schedule_conditioned_bias": ',
         json_number(max_schedule_bias), ','),
  paste0('    "minimum_naive_sign_error_rate": ',
         json_number(minimum_naive_sign_error)),
  "  },",
  '  "validation": {',
  paste0('    "checks_passed": ', sum(result$passed), ','),
  paste0('    "checks_total": ', nrow(result), ','),
  '    "status": "passed"',
  "  }",
  "}"
)
writeLines(manifest_lines, manifest_path, useBytes = TRUE)
message("Model/calibration extension validation passed all ", nrow(result),
        " checks: ", path)
message("Model/calibration extension manifest: ", manifest_path)
