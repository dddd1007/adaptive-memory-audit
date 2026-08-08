project_dir <- Sys.getenv(
  "BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE)
)
source(file.path(project_dir, "src", "R", "theme_brm.R"))
source(file.path(project_dir, "src", "R", "validate_public_data.R"))

output_dir <- file.path(project_dir, "outputs")
figure_dir <- file.path(project_dir, "figures")
table_dir <- file.path(project_dir, "tables")
validation_dir <- file.path(project_dir, "validation")
if (!dir.exists(validation_dir)) dir.create(validation_dir, recursive = TRUE)

checks <- list()
add_check <- function(name, passed, observed, expected) {
  checks[[length(checks) + 1L]] <<- data.frame(
    check = as.character(name),
    passed = isTRUE(passed),
    observed = paste(as.character(observed), collapse = ""),
    expected = paste(as.character(expected), collapse = ""),
    stringsAsFactors = FALSE
  )
}

read_output <- function(filename) {
  brm_read_csv(file.path(output_dir, filename))
}

mean_selected_rmse <- function(data, learner, condition) {
  selected_rows <-
    data$learner_family == learner &
      data$observation_condition == condition &
      data$candidate_family == data$selected_family
  values <- as.numeric(data$state_rmse[selected_rows])
  mean(values[is.finite(values)])
}

# ---------------------------------------------------------------------------
# Original known-truth studies and their central conclusions (checks 1--10).
# ---------------------------------------------------------------------------
parameter_rep <- read_output("parameter_recovery_replications.csv")
context_rep <- read_output("context_violation_replications.csv")
mrt_rep <- read_output("micro_randomized_replications.csv")
state_delay_rep <- read_output("state_model_delay_recovery_replications.csv")
state_diag_rep <- read_output("state_model_diagnostics_replications.csv")

add_check(
  "parameter-recovery replication rows",
  nrow(parameter_rep) == 19200L, nrow(parameter_rep), 19200L
)
add_check(
  "context-violation replication rows",
  nrow(context_rep) == 5760L, nrow(context_rep), 5760L
)
add_check(
  "legacy micro-randomized benchmark rows",
  nrow(mrt_rep) == 12000L, nrow(mrt_rep), 12000L
)
add_check(
  "state-model delay-recovery rows",
  nrow(state_delay_rep) == 8000L, nrow(state_delay_rep), 8000L
)
add_check(
  "state-model diagnostic datasets",
  nrow(state_diag_rep) == 800L, nrow(state_diag_rep), 800L
)

parameter <- read_output("parameter_recovery_summary.csv")
baseline <- parameter[
  parameter$n_participants == 12 &
    abs(parameter$adaptivity - 0.75) < 1e-12 &
    abs(parameter$state_reliability - 0.80) < 1e-12,
  , drop = FALSE
]
naive <- baseline[baseline$estimator == "naive", , drop = FALSE]
schedule <- baseline[
  baseline$estimator == "schedule_conditioned", , drop = FALSE
]
core_closed_loop_ok <-
  nrow(naive) == 1L && nrow(schedule) == 1L &&
  naive$sign_error_rate >= 0.99 &&
  schedule$coverage_95 >= 0.90 && schedule$coverage_95 <= 0.98
add_check(
  "closed-loop sign reversal and logged-plan recovery",
  core_closed_loop_ok,
  sprintf(
    "naive sign error %.3f; plan coverage %.3f",
    naive$sign_error_rate, schedule$coverage_95
  ),
  "sign error >= .990; coverage .900--.980"
)

family <- read_output("state_model_family_recovery_summary.csv")
rl2 <- family$family_recovery_rate[
  family$learner_family == "RL" & family$n_stages == 2
]
rl8 <- family$family_recovery_rate[
  family$learner_family == "RL" & family$n_stages == 8
]
bayes2 <- family$family_recovery_rate[
  family$learner_family == "Bayesian" & family$n_stages == 2
]
bayes8 <- family$family_recovery_rate[
  family$learner_family == "Bayesian" & family$n_stages == 8
]
history_ok <-
  length(c(rl2, rl8, bayes2, bayes8)) == 4L &&
  rl8 > rl2 && bayes8 > bayes2 && rl8 >= 0.90 && bayes8 >= 0.95
add_check(
  "longer histories improve both learner-family recoveries",
  history_ok,
  sprintf("RL %.3f->%.3f; Bayesian %.3f->%.3f", rl2, rl8, bayes2, bayes8),
  "both increase; eight-stage RL >= .90 and Bayesian >= .95"
)

state <- read_output("state_model_state_recovery_summary.csv")
state8 <- state[state$n_stages == 8, , drop = FALSE]
learner_match <- brm_match_bool(state8$learner_family_match, TRUE)
policy_match <- brm_match_bool(state8$scheduler_family_match, TRUE)
state_matched <- mean(state8$state_rmse[learner_match])
state_mismatched <- mean(state8$state_rmse[!learner_match])
policy_matched <- mean(state8$policy_rmse[policy_match])
policy_mismatched <- mean(state8$policy_rmse[!policy_match])
add_check(
  "model-family matching improves state and policy reconstruction",
  state_matched < state_mismatched && policy_matched < policy_mismatched,
  sprintf(
    "state %.3f<%.3f; policy %.3f<%.3f",
    state_matched, state_mismatched, policy_matched, policy_mismatched
  ),
  "matched RMSE < mismatched RMSE for both layers"
)

layer <- read_output("state_model_layer_summary.csv")
logged8 <- layer[
  layer$n_stages == 8 & layer$design == "Logged plan", , drop = FALSE
]
logged_cell_absolute_bias <- mean(abs(logged8$bias))
logged_coverage <- mean(logged8$coverage_95)
add_check(
  "logged plans retain low delay bias across model layers",
  nrow(logged8) > 0L && logged_cell_absolute_bias < 0.015 &&
    logged_coverage >= 0.90,
  sprintf("absolute bias %.4f; coverage %.3f", logged_cell_absolute_bias,
          logged_coverage),
  "absolute cell bias < .015 and coverage >= .900"
)

context <- read_output("context_violation_summary.csv")
context_bad <- context$bias[
  context$n_participants == 30 & context$context_strength == 1 &
    context$estimator == "schedule_conditioned"
]
context_oracle <- context$bias[
  context$n_participants == 30 & context$context_strength == 1 &
    context$estimator == "oracle"
]
add_check(
  "unlogged context breaks plan conditioning while oracle remains stable",
  length(context_bad) == 1L && length(context_oracle) == 1L &&
    context_bad > 0.40 && abs(context_oracle) < 0.02,
  sprintf("plan bias %.3f; oracle bias %.3f", context_bad, context_oracle),
  "plan bias > .400; absolute oracle bias < .020"
)

# ---------------------------------------------------------------------------
# Reviewer-requested WCLS benchmark (checks 11--15).
# ---------------------------------------------------------------------------
wcls_rep <- read_output("micro_randomized_wcls_replications.csv")
wcls <- read_output("micro_randomized_wcls_summary.csv")

add_check(
  "WCLS replication rows",
  nrow(wcls_rep) == 12000L, nrow(wcls_rep), 12000L
)
wcls_cells <- interaction(
  wcls_rep$n_participants, wcls_rep$perturbation_log_days,
  wcls_rep$perturbation_compliance, wcls_rep$estimator,
  drop = TRUE, lex.order = TRUE
)
wcls_cell_counts <- table(wcls_cells)
add_check(
  "WCLS factorial cells contain 250 replications",
  length(wcls_cell_counts) == 48L && all(wcls_cell_counts == 250L),
  sprintf("%d cells; range %d--%d", length(wcls_cell_counts),
          min(wcls_cell_counts), max(wcls_cell_counts)),
  "48 cells; 250 per cell"
)
wcls_target <- wcls[
  wcls$n_participants == 100 &
    abs(wcls$perturbation_log_days - 0.30) < 1e-12 &
    abs(wcls$perturbation_compliance - 1.00) < 1e-12 &
    wcls$estimator == "covariate_adjusted_itt",
  , drop = FALSE
]
add_check(
  "WCLS target participant-t power",
  nrow(wcls_target) == 1L &&
    abs(wcls_target$participant_t_power_positive - 0.872) < 1e-12,
  sprintf("%.3f", wcls_target$participant_t_power_positive),
  ".872"
)
target_reference_ok <-
  nrow(wcls_target) == 1L &&
  abs(wcls_target$normal_coverage_95 - 0.960) < 1e-12 &&
  abs(wcls_target$participant_t_coverage_95 - 0.960) < 1e-12 &&
  abs(wcls_target$normal_power_positive - 0.876) < 1e-12
add_check(
  "WCLS target normal and participant-t operating characteristics",
  target_reference_ok,
  sprintf(
    "coverage %.3f/%.3f; normal power %.3f",
    wcls_target$normal_coverage_95,
    wcls_target$participant_t_coverage_95,
    wcls_target$normal_power_positive
  ),
  "coverage .960/.960; normal power .876"
)
max_wcls_estimate_difference <- max(
  abs(wcls_rep$estimate - wcls_rep$lpm_estimate), na.rm = TRUE
)
max_wcls_se_difference <- max(
  abs(wcls_rep$standard_error - wcls_rep$lpm_standard_error), na.rm = TRUE
)
add_check(
  "fixed-probability WCLS equals its additive LPM benchmark",
  max_wcls_estimate_difference < 1e-12 && max_wcls_se_difference < 1e-12,
  sprintf("estimate %.4e; SE %.4e", max_wcls_estimate_difference,
          max_wcls_se_difference),
  "both maximum absolute differences < 1e-12"
)

# ---------------------------------------------------------------------------
# Partially observed context stress test (checks 16--20).
# ---------------------------------------------------------------------------
proxy_rep <- read_output("context_proxy_replications.csv")
proxy <- read_output("context_proxy_summary.csv")
add_check(
  "context-proxy replication rows",
  nrow(proxy_rep) == 2400L, nrow(proxy_rep), 2400L
)
proxy_cells <- interaction(
  proxy_rep$proxy_r2_target, proxy_rep$estimator,
  drop = TRUE, lex.order = TRUE
)
proxy_cell_counts <- table(proxy_cells)
add_check(
  "context-proxy factorial cells contain 160 replications",
  length(proxy_cell_counts) == 15L && all(proxy_cell_counts == 160L) &&
    length(unique(proxy_rep$estimator)) == 3L,
  sprintf("%d cells; range %d--%d; %d estimators", length(proxy_cell_counts),
          min(proxy_cell_counts), max(proxy_cell_counts),
          length(unique(proxy_rep$estimator))),
  "15 cells; 160 per cell; 3 estimators"
)
proxy_r2_error <- max(
  abs(proxy_rep$proxy_r2_realized - proxy_rep$proxy_r2_target),
  na.rm = TRUE
)
add_check(
  "context proxies attain finite-sample target R-squared",
  proxy_r2_error < 1e-10,
  sprintf("%.4e", proxy_r2_error), "maximum absolute error < 1e-10"
)
proxy_curve <- proxy[
  proxy$estimator == "plan_plus_proxy", , drop = FALSE
]
proxy_curve <- proxy_curve[order(proxy_curve$proxy_r2_target), , drop = FALSE]
add_check(
  "context-proxy bias decreases monotonically with measurement quality",
  nrow(proxy_curve) == 5L && all(diff(proxy_curve$bias) < 0),
  paste(sprintf("%.4f", proxy_curve$bias), collapse = " -> "),
  "strictly decreasing over R2 = 0, .25, .50, .75, 1"
)
proxy_t_coverage <- proxy_curve$participant_t_coverage_95
add_check(
  "context-proxy participant-t coverage improves with measurement quality",
  length(proxy_t_coverage) == 5L && all(diff(proxy_t_coverage) >= 0) &&
    proxy_t_coverage[1L] == 0 && proxy_t_coverage[5L] >= 0.95,
  paste(sprintf("%.4f", proxy_t_coverage), collapse = " -> "),
  "nondecreasing; starts at 0 and ends >= .95"
)

# ---------------------------------------------------------------------------
# Observation-equation and calibration extensions (checks 21--28).
# ---------------------------------------------------------------------------
misspec_rep <- read_output("state_model_observation_misspec_replications.csv")
misspec <- read_output("state_model_observation_misspec_summary.csv")
add_check(
  "observation-misspecification replication rows",
  nrow(misspec_rep) == 4000L, nrow(misspec_rep), 4000L
)
add_check(
  "observation-misspecification summary rows",
  nrow(misspec) == 40L, nrow(misspec), 40L
)
true_family_rows <- misspec[
  misspec$candidate_family == misspec$learner_family, , drop = FALSE
]
correct_rl <- mean(true_family_rows$family_recovery_rate[
  true_family_rows$learner_family == "RL" &
    true_family_rows$observation_condition == "correct"
])
correct_bayes <- mean(true_family_rows$family_recovery_rate[
  true_family_rows$learner_family == "Bayesian" &
    true_family_rows$observation_condition == "correct"
])
correct_rl_rmse <- mean_selected_rmse(misspec_rep, "RL", "correct")
correct_bayes_rmse <- mean_selected_rmse(
  misspec_rep, "Bayesian", "correct"
)
add_check(
  "correct observation equation permits strong family and state recovery",
  correct_rl >= 0.95 && correct_bayes >= 0.98 &&
    correct_rl_rmse < 0.35 && correct_bayes_rmse < 0.30,
  sprintf(
    "recovery RL %.3f/Bayesian %.3f; RMSE %.3f/%.3f",
    correct_rl, correct_bayes, correct_rl_rmse, correct_bayes_rmse
  ),
  "recovery >= .95/.98 and selected-state RMSE < .35/.30"
)
combined_condition <- "combined_plus_0.25_slope_0.80"
combined_rl <- mean(true_family_rows$family_recovery_rate[
  true_family_rows$learner_family == "RL" &
    true_family_rows$observation_condition == combined_condition
])
combined_bayes <- mean(true_family_rows$family_recovery_rate[
  true_family_rows$learner_family == "Bayesian" &
    true_family_rows$observation_condition == combined_condition
])
combined_rl_rmse <- mean_selected_rmse(
  misspec_rep, "RL", combined_condition
)
combined_bayes_rmse <- mean_selected_rmse(
  misspec_rep, "Bayesian", combined_condition
)
add_check(
  "combined observation error exposes the declared recovery boundary",
  abs(combined_rl - 1) < 1e-12 && abs(combined_bayes) < 1e-12 &&
    combined_rl_rmse > correct_rl_rmse &&
    combined_bayes_rmse > correct_bayes_rmse,
  sprintf(
    "recovery RL %.3f/Bayesian %.3f; RMSE %.3f/%.3f",
    combined_rl, combined_bayes, combined_rl_rmse, combined_bayes_rmse
  ),
  "recovery 1/0 and selected-state RMSE worse for both families"
)

calibration_search <- read_output("calibration_search_candidates.csv")
calibration_rep <- read_output("calibration_robustness_replications.csv")
calibration <- read_output("calibration_robustness_summary.csv")
add_check(
  "calibration search contains the full 27-candidate grid",
  nrow(calibration_search) == 27L &&
    identical(sort(as.integer(calibration_search$rank)), 1:27),
  sprintf("%d candidates; ranks %d--%d", nrow(calibration_search),
          min(calibration_search$rank), max(calibration_search$rank)),
  "27 candidates with unique ranks 1--27"
)
add_check(
  "top-five calibration robustness contains 2,000 estimator rows",
  nrow(calibration_rep) == 2000L && nrow(calibration) == 20L,
  sprintf("raw %d; summary %d", nrow(calibration_rep), nrow(calibration)),
  "2,000 raw rows and 20 summary cells"
)
rank_one <- calibration_search[
  calibration_search$rank == 1L, , drop = FALSE
]
recomputed_loss <-
  ((calibration_search$mean_failure_rate -
      calibration_search$target_failure_rate) / 0.015)^2 +
  ((calibration_search$mean_scheduled_actual_spearman -
      calibration_search$target_scheduled_actual_spearman) / 0.05)^2 +
  ((calibration_search$mean_naive_log_odds -
      calibration_search$target_naive_log_odds) / 0.10)^2
loss_error <- max(abs(recomputed_loss - calibration_search$loss))
rank_one_ok <-
  nrow(rank_one) == 1L &&
  abs(rank_one$adaptivity - 0.75) < 1e-12 &&
  abs(rank_one$natural_deviation_sd - 0.50) < 1e-12 &&
  abs(rank_one$outcome_intercept + 2.40) < 1e-12 &&
  rank_one$loss < 0.01 && loss_error < 1e-12
add_check(
  "rank-one calibration and prespecified loss reproduce",
  rank_one_ok,
  sprintf(
    "(%.2f, %.2f, %.2f), loss %.6f; max loss error %.2e",
    rank_one$adaptivity, rank_one$natural_deviation_sd,
    rank_one$outcome_intercept, rank_one$loss, loss_error
  ),
  "(.75, .50, -2.40), loss < .01, recomputation error < 1e-12"
)
top_five <- calibration[calibration$calibration_rank <= 5, , drop = FALSE]
top_five_schedule <- top_five[
  top_five$estimator == "schedule_conditioned", , drop = FALSE
]
top_five_naive <- top_five[
  top_five$estimator == "naive", , drop = FALSE
]
max_schedule_bias <- max(abs(top_five_schedule$bias))
min_naive_sign_error <- min(top_five_naive$sign_error_rate)
add_check(
  "top-five calibrations preserve plan recovery and naive sign reversal",
  nrow(top_five_schedule) == 5L && nrow(top_five_naive) == 5L &&
    max_schedule_bias < 0.02 && min_naive_sign_error >= 0.98,
  sprintf("max plan |bias| %.5f; min naive sign error %.3f",
          max_schedule_bias, min_naive_sign_error),
  "max plan |bias| < .02 and min naive sign error >= .98"
)

# ---------------------------------------------------------------------------
# Minimized public empirical data and sensitivity audit (checks 29--34).
# ---------------------------------------------------------------------------
public_data_audit <- brm_validate_public_data(project_dir)
public_data <- public_data_audit$data
required_public_fields <- public_data_audit$required_fields
add_check(
  "public empirical input has the minimized schema and pseudonymous codes",
  nrow(public_data) == 2045L &&
    identical(names(public_data)[seq_along(required_public_fields)],
              required_public_fields) &&
    length(unique(public_data$participant_code)) == 12L &&
    length(unique(public_data$card_code)) == 100L,
  sprintf(
    "%d rows; %d fields; %d participant codes; %d card codes",
    nrow(public_data), length(required_public_fields),
    length(unique(public_data$participant_code)),
    length(unique(public_data$card_code))
  ),
  "2,045 rows; 10 fields; 12 participant codes; 100 card codes"
)
n_transitions <- nrow(public_data)
n_failures <- sum(public_data$next_first_ease == 1)
add_check(
  "empirical transition and failure counts reproduce",
  n_transitions == 2045L && n_failures == 199L,
  sprintf("transitions %d; failures %d", n_transitions, n_failures),
  "2,045 transitions and 199 Again failures"
)
public_data_validation <- public_data_audit$checks
public_data_passed <- brm_match_bool(public_data_validation$passed, TRUE)
add_check(
  "all public empirical-data checks pass",
  nrow(public_data_validation) == 8L && all(public_data_passed),
  sprintf("%d of %d", sum(public_data_passed), nrow(public_data_validation)),
  "8 of 8"
)
sensitivity <- read_output("empirical_sensitivity_summary.csv")
sensitivity_rep <- read_output(
  "empirical_sensitivity_bootstrap_replications.csv"
)
sensitivity_cells <- interaction(
  sensitivity$planned_interval_threshold_days,
  sensitivity$outcome_definition,
  sensitivity$regularization_C,
  drop = TRUE, lex.order = TRUE
)
add_check(
  "empirical sensitivity grid retains all 24 by 1,200 bootstraps",
  nrow(sensitivity) == 24L && length(unique(sensitivity_cells)) == 24L &&
    all(sensitivity$n_bootstrap_requested == 1200L) &&
    all(sensitivity$n_bootstrap_successful == 1200L) &&
    nrow(sensitivity_rep) == 28800L,
  sprintf("%d cells; successful range %d--%d; raw %d",
          nrow(sensitivity), min(sensitivity$n_bootstrap_successful),
          max(sensitivity$n_bootstrap_successful), nrow(sensitivity_rep)),
  "24 cells, 1,200 successful per cell, 28,800 raw rows"
)
loo <- read_output("empirical_leave_one_participant_out.csv")
loo_again <- loo[loo$outcome_definition == "Again", , drop = FALSE]
loo_broad <- loo[loo$outcome_definition == "Again or Hard", , drop = FALSE]
loo_ok <-
  nrow(loo) == 24L && nrow(loo_again) == 12L && nrow(loo_broad) == 12L &&
  sum(loo_again$odds_ratio < 1) == 11L &&
  sum(loo_broad$odds_ratio < 1) == 12L &&
  min(loo_again$odds_ratio) < 0.60 && max(loo_again$odds_ratio) > 1.30 &&
  max(loo_broad$odds_ratio) < 1
add_check(
  "leave-one-participant-out influence pattern reproduces",
  loo_ok,
  sprintf(
    "Again %.3f--%.3f (%d/12 <1); broad %.3f--%.3f (%d/12 <1)",
    min(loo_again$odds_ratio), max(loo_again$odds_ratio),
    sum(loo_again$odds_ratio < 1), min(loo_broad$odds_ratio),
    max(loo_broad$odds_ratio), sum(loo_broad$odds_ratio < 1)
  ),
  "Again 11/12 below 1 and range crosses 1; broad 12/12 below 1"
)
mixed <- read_output("empirical_mixed_effects_benchmark.csv")
mixed_converged <- brm_match_bool(mixed$converged, TRUE)
mixed_singular <- brm_match_bool(mixed$singular, TRUE)
mixed_broad <- mixed[
  mixed$outcome_definition == "Again or Hard", , drop = FALSE
]
add_check(
  "crossed mixed-effects empirical benchmark is stable and informative",
  nrow(mixed) == 2L && all(mixed_converged) && !any(mixed_singular) &&
    nrow(mixed_broad) == 1L && mixed_broad$odds_ratio < 1 &&
    mixed_broad$ci_high < 1 && mixed_broad$p_value < 0.05,
  sprintf(
    "converged %d/2; singular %d/2; broad OR %.3f [%.3f, %.3f], p %.3f",
    sum(mixed_converged), sum(mixed_singular), mixed_broad$odds_ratio,
    mixed_broad$ci_low, mixed_broad$ci_high, mixed_broad$p_value
  ),
  "2/2 converged, 0/2 singular; broad outcome CI below 1 and p < .05"
)

# ---------------------------------------------------------------------------
# Publication assets and source-language boundary (checks 35--38).
# ---------------------------------------------------------------------------
figure_stems <- c(
  "figure_1_sequential_learning_program",
  "figure_2_parameter_recovery",
  "figure_3_state_model_misspecification",
  "figure_4_robustness_and_experiment",
  "figure_5_empirical_case",
  "supplementary_figure_s1_event_reconstruction",
  "supplementary_figure_s2_calibration_robustness",
  "supplementary_figure_s3_reviewer_robustness"
)
png_paths <- file.path(figure_dir, paste0(figure_stems, ".png"))
svg_paths <- file.path(figure_dir, paste0(figure_stems, ".svg"))
png_sizes <- file.info(png_paths)$size
svg_sizes <- file.info(svg_paths)$size
add_check(
  "five main and three supplementary PNG figures are publication-ready",
  all(file.exists(png_paths)) && all(is.finite(png_sizes)) &&
    all(png_sizes > 1000),
  sprintf("%d of 8 present; minimum size %s",
          sum(file.exists(png_paths)),
          if (any(is.finite(png_sizes))) min(png_sizes, na.rm = TRUE) else "NA"),
  "8 of 8 present and every file > 1,000 bytes"
)
add_check(
  "five main and three supplementary SVG figures are publication-ready",
  all(file.exists(svg_paths)) && all(is.finite(svg_sizes)) &&
    all(svg_sizes > 1000),
  sprintf("%d of 8 present; minimum size %s",
          sum(file.exists(svg_paths)),
          if (any(is.finite(svg_sizes))) min(svg_sizes, na.rm = TRUE) else "NA"),
  "8 of 8 present and every file > 1,000 bytes"
)

main_table_stems <- c(
  "table_1_audit_components", "table_2_calibration",
  "table_3_parameter_recovery", "table_4_state_model_recovery",
  "table_5_context_sensitivity", "table_6_micro_randomized",
  "table_7_empirical_audit"
)
supplementary_table_stems <- c(
  "supplementary_table_s1_calibration_search",
  "supplementary_table_s2_calibration_robustness",
  "supplementary_table_s3_observation_misspecification",
  "supplementary_table_s4_context_proxy",
  "supplementary_table_s5_wcls_grid",
  "supplementary_table_s6_wcls_equivalence",
  "supplementary_table_s7_empirical_sensitivity",
  "supplementary_table_s8_empirical_influence_mixed"
)
table_paths <- c(
  file.path(table_dir, as.vector(outer(main_table_stems, c(".csv", ".md"), paste0))),
  file.path(table_dir, as.vector(outer(
    supplementary_table_stems, c(".csv", ".md"), paste0
  )))
)
table_sizes <- file.info(table_paths)$size
add_check(
  "seven main and eight supplementary tables have CSV and Markdown forms",
  length(table_paths) == 30L && all(file.exists(table_paths)) &&
    all(is.finite(table_sizes)) && all(table_sizes > 50),
  sprintf("%d of 30 present; minimum size %s",
          sum(file.exists(table_paths)),
          if (any(is.finite(table_sizes))) min(table_sizes, na.rm = TRUE) else "NA"),
  "30 of 30 files: 7 main and 8 supplementary CSV/Markdown pairs"
)

source_files <- list.files(
  file.path(project_dir, "src"), recursive = TRUE, full.names = TRUE
)
source_files <- source_files[!file.info(source_files)$isdir]
source_extensions <- tolower(tools::file_ext(source_files))
source_only_ok <-
  length(source_files) > 0L &&
  all(source_extensions %in% c("r", "jl")) &&
  any(source_extensions == "r") && any(source_extensions == "jl")
add_check(
  "computational source tree is R and Julia only",
  source_only_ok,
  paste(sprintf("%s=%d", sort(unique(source_extensions)),
                as.integer(table(source_extensions)[sort(unique(source_extensions))])),
        collapse = "; "),
  "only .R and .jl source files, with both languages represented"
)

# Exactly 38 publication checks are part of the machine-readable contract.
if (length(checks) != 38L) {
  stop("Internal validation specification error: expected 38 checks, found ",
       length(checks))
}

results <- do.call(rbind, checks)
write.csv(
  results, file.path(validation_dir, "r_julia_validation.csv"),
  row.names = FALSE
)

json_escape <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  x <- gsub('"', '\\"', x, fixed = TRUE)
  gsub("[\r\n]", " ", x)
}
json_rows <- vapply(seq_len(nrow(results)), function(index) {
  row <- results[index, ]
  sprintf(
    '    {"check":"%s","passed":%s,"observed":"%s","expected":"%s"}',
    json_escape(row$check), tolower(as.character(row$passed)),
    json_escape(row$observed), json_escape(row$expected)
  )
}, character(1L))
json <- c(
  "{",
  sprintf('  "all_passed": %s,',
          tolower(as.character(all(results$passed)))),
  sprintf('  "n_checks": %d,', nrow(results)),
  '  "engine": "R base public-repository validation over Julia/R outputs",',
  '  "checks": [',
  paste(json_rows, collapse = ",\n"),
  "  ]",
  "}"
)
writeLines(
  json, file.path(validation_dir, "r_julia_validation.json"),
  useBytes = TRUE
)

if (!all(results$passed)) {
  failed <- results$check[!results$passed]
  stop("Validation failed: ", paste(failed, collapse = "; "))
}
message("R validation passed all ", nrow(results), " checks.")
