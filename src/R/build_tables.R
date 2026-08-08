project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source(file.path(project_dir, "src", "R", "theme_brm.R"))
output_dir <- file.path(project_dir, "outputs")
table_dir <- file.path(project_dir, "tables")
if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

fmt <- function(x, digits = 3L, leading_zero = TRUE) {
  value <- formatC(as.numeric(x), format = "f", digits = digits)
  if (!leading_zero) value <- sub("^(-?)0\\.", "\\1.", value)
  value
}

fmt_sci <- function(x, digits = 2L) {
  formatC(as.numeric(x), format = "e", digits = digits)
}

weighted_mean_safe <- function(x, weights) {
  x <- as.numeric(x)
  weights <- as.numeric(weights)
  keep <- is.finite(x) & is.finite(weights) & weights > 0
  if (!any(keep)) return(NA_real_)
  stats::weighted.mean(x[keep], weights[keep])
}

escape_markdown <- function(x) gsub("\\|", "\\\\|", as.character(x))

write_markdown_table <- function(data, title, filename, align = NULL) {
  if (is.null(align)) align <- c("left", rep("right", ncol(data) - 1L))
  markers <- vapply(align, function(item) {
    switch(item, left = ":--", center = ":--:", right = "--:", ":--")
  }, character(1L))
  header <- paste0("| ", paste(names(data), collapse = " | "), " |")
  separator <- paste0("|", paste(markers, collapse = "|"), "|")
  rows <- apply(data, 1L, function(row) {
    paste0("| ", paste(escape_markdown(row), collapse = " | "), " |")
  })
  lines <- c(paste0("**", title, "**"), "", header, separator, rows, "")
  writeLines(lines, file.path(table_dir, filename), useBytes = TRUE)
}

write_table <- function(data, title, stem, align = NULL) {
  write.csv(data, file.path(table_dir, paste0(stem, ".csv")), row.names = FALSE)
  write_markdown_table(data, title, paste0(stem, ".md"), align)
  message("R table: ", stem, " [", nrow(data), " rows]")
}

markdown_table_lines <- function(data, align = NULL) {
  if (is.null(align)) align <- c("left", rep("right", ncol(data) - 1L))
  markers <- vapply(align, function(item) {
    switch(item, left = ":--", center = ":--:", right = "--:", ":--")
  }, character(1L))
  c(
    paste0("| ", paste(names(data), collapse = " | "), " |"),
    paste0("|", paste(markers, collapse = "|"), "|"),
    apply(data, 1L, function(row) paste0("| ", paste(escape_markdown(row), collapse = " | "), " |"))
  )
}

write_panel_markdown <- function(title, panels, labels, filename, aligns = NULL, note = NULL) {
  lines <- c(paste0("**", title, "**"), "")
  for (index in seq_along(panels)) {
    lines <- c(lines, paste0("*", labels[index], "*"), "",
               markdown_table_lines(panels[[index]], if (is.null(aligns)) NULL else aligns[[index]]), "")
  }
  if (!is.null(note)) lines <- c(lines, paste0("*Note.* ", note), "")
  writeLines(lines, file.path(table_dir, filename), useBytes = TRUE)
}

table_1 <- data.frame(
  Component = c(
    "Empirical estimator audit", "Parameter recovery", "Cognitive-model recovery",
    "Context-sensitivity study", "Micro-randomized study", "Computational validation"
  ),
  `Psychological question` = c(
    "Does observed spacing reflect memory or scheduler selection?",
    "Can the within-state delay effect be recovered?",
    "Which memory-update rule, state, and prediction error are recoverable?",
    "Do attention and availability confound enacted timing and recall?",
    "What is the proximal causal effect of an interval perturbation?",
    "Are the psychological conclusions computationally reproducible?"
  ),
  `Evidence provided` = c(
    "Direction and interval stability",
    "Bias, RMSE, 95% coverage, and sign error",
    "Family, parameter, correlation, and RMSE",
    "Bias and coverage across context strength",
    "Risk-difference bias, coverage, and power",
    "Cross-run agreement and structural checks"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_1, "Table 1. Psychological questions addressed by the audit",
            "table_1_audit_components", rep("left", 3L))

calibration <- brm_read_csv(file.path(output_dir, "empirical_simulation_calibration.csv"))
calibration_labels <- c(
  "Transitions per dataset", "Failure rate", "Planned-enacted Spearman correlation",
  "Naive gap odds ratio per SD"
)
table_2 <- data.frame(
  `Calibration target` = calibration_labels,
  Empirical = c(
    formatC(calibration$observed_value[1L], format = "f", digits = 0, big.mark = ","),
    fmt(calibration$observed_value[2L], 4L, FALSE),
    fmt(calibration$observed_value[3L], 4L, FALSE),
    fmt(calibration$observed_value[4L], 3L, TRUE)
  ),
  `Simulated mean` = c(
    formatC(calibration$simulated_mean[1L], format = "f", digits = 0, big.mark = ","),
    fmt(calibration$simulated_mean[2L], 4L, FALSE),
    fmt(calibration$simulated_mean[3L], 4L, FALSE),
    fmt(calibration$simulated_mean[4L], 3L, TRUE)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_2, "Table 2. Empirical-to-simulation calibration",
            "table_2_calibration")

recovery <- brm_read_csv(file.path(output_dir, "parameter_recovery_summary.csv"))
selected <- recovery[
  recovery$n_participants == 12 & abs(recovery$adaptivity - 0.75) < 1e-10 &
    abs(recovery$state_reliability - 0.80) < 1e-10,
  , drop = FALSE
]
estimator_order <- c("naive", "history", "schedule_conditioned", "oracle")
selected <- selected[match(estimator_order, selected$estimator), ]
table_3 <- data.frame(
  Estimator = c("Naive", "History adjusted", "Schedule conditioned", "Oracle"),
  Truth = fmt(selected$mean_truth, 3L, FALSE),
  Estimate = fmt(selected$mean_estimate, 3L, FALSE),
  Bias = fmt(selected$bias, 3L, FALSE),
  RMSE = fmt(selected$rmse, 3L, FALSE),
  Coverage = fmt(selected$coverage_95, 3L, FALSE),
  `Sign error` = fmt(selected$sign_error_rate, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_3, "Table 3. Parameter recovery in the empirical-like condition",
            "table_3_parameter_recovery")
write_panel_markdown(
  "Table 3. Parameter recovery in the empirical-like condition",
  list(table_3[c("Estimator", "Truth", "Estimate", "Bias")],
       table_3[c("Estimator", "RMSE", "Coverage", "Sign error")]),
  c("Panel A. Point recovery", "Panel B. Error and interval performance"),
  "table_3_parameter_recovery.md"
)

family <- brm_read_csv(file.path(output_dir, "state_model_family_recovery_summary.csv"))
family$order <- match(family$learner_family, c("RL", "Bayesian")) * 10 + family$n_stages
family <- family[order(family$order), ]
table_4 <- data.frame(
  Generator = family$learner_family,
  Transitions = family$n_stages,
  `Family recovery` = fmt(family$family_recovery_rate, 3L, FALSE),
  `True par.` = fmt(family$true_parameter, 3L, FALSE),
  `Mean par.` = fmt(family$mean_parameter_within_true_family, 3L, FALSE),
  `Exact recovery` = fmt(family$parameter_exact_recovery_rate, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_4, "Table 4. Cognitive update-family and parameter recovery",
            "table_4_state_model_recovery")
table_4_panel_b <- table_4[c("Generator", "Transitions", "True par.", "Mean par.", "Exact recovery")]
names(table_4_panel_b) <- c("Model", "Transitions", "True", "Mean", "Exact")
table_4_panel_b$Model <- ifelse(table_4_panel_b$Model == "Bayesian", "Bayes.", "RL")
write_panel_markdown(
  "Table 4. Cognitive update-family and parameter recovery",
  list(table_4[c("Generator", "Transitions", "Family recovery")],
       table_4_panel_b),
  c("Panel A. Model-family recovery", "Panel B. Parameter recovery within the true family"),
  "table_4_state_model_recovery.md"
)

context <- brm_read_csv(file.path(output_dir, "context_violation_summary.csv"))
context <- context[context$n_participants == 30 &
                     context$estimator %in% c("schedule_conditioned", "oracle"), , drop = FALSE]
context$order <- match(context$estimator, c("schedule_conditioned", "oracle")) * 10 + context$context_strength
context <- context[order(context$order), ]
table_5 <- data.frame(
  c = fmt(context$context_strength, 2L),
  `Est.` = ifelse(context$estimator == "schedule_conditioned", "Plan", "Oracle"),
  Truth = fmt(context$mean_truth, 3L, FALSE),
  Bias = fmt(context$bias, 3L, FALSE),
  `Cov.` = fmt(context$coverage_95, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_5, "Table 5. Context sensitivity at 30 participants",
            "table_5_context_sensitivity", c("right", "left", "right", "right", "right"))
write_panel_markdown(
  "Table 5. Context sensitivity at 30 participants",
  list(table_5), "Context stress test", "table_5_context_sensitivity.md",
  list(c("right", "left", "right", "right", "right")),
  "c denotes context strength; Plan denotes the schedule-conditioned estimator."
)

wcls <- brm_read_csv(file.path(output_dir, "micro_randomized_wcls_summary.csv"))
mrt <- wcls[wcls$estimator == "covariate_adjusted_itt" &
              abs(wcls$perturbation_log_days - 0.30) < 1e-10 &
              abs(wcls$perturbation_compliance - 1.00) < 1e-10, , drop = FALSE]
mrt <- mrt[order(mrt$n_participants), , drop = FALSE]
table_6 <- data.frame(
  N = mrt$n_participants,
  Truth = fmt(mrt$mean_truth, 4L, FALSE),
  Estimate = fmt(mrt$mean_estimate, 4L, FALSE),
  Bias = fmt(mrt$bias, 4L, FALSE),
  RMSE = fmt(mrt$rmse, 4L, FALSE),
  `Norm. cov.` = fmt(mrt$normal_coverage_95, 3L, FALSE),
  `t cov.` = fmt(mrt$participant_t_coverage_95, 3L, FALSE),
  `t power` = fmt(mrt$participant_t_power_positive, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_6,
            "Table 6. Covariate-adjusted WCLS performance for a ±0.30-log-day perturbation under full compliance",
            "table_6_micro_randomized")
write_panel_markdown(
  "Table 6. Covariate-adjusted WCLS performance for a ±0.30-log-day perturbation under full compliance",
  list(table_6[c("N", "Truth", "Estimate", "Bias")],
       table_6[c("N", "RMSE", "Norm. cov.", "t cov.", "t power")]),
  c("Panel A. Effect recovery", "Panel B. Interval performance"),
  "table_6_micro_randomized.md",
  note = paste0(
    "WCLS is the covariate-adjusted, centered-treatment estimator. ",
    "Participant-t intervals use learner clusters as the degrees-of-freedom unit."
  )
)

empirical <- brm_read_csv(file.path(output_dir, "empirical_model_comparison.csv"))
empirical_sensitivity <- brm_read_csv(
  file.path(output_dir, "empirical_sensitivity_summary.csv")
)
focal_empirical <- empirical_sensitivity[
  abs(empirical_sensitivity$planned_interval_threshold_days - 0.5) < 1e-10 &
    empirical_sensitivity$outcome_definition == "Again" &
    empirical_sensitivity$regularization_C == 10,
  , drop = FALSE
]
if (nrow(focal_empirical) != 1L) {
  stop("Expected exactly one empirical sensitivity row for .5 day, C=10, Again.")
}
# The schedule-conditioned point model is unchanged, but its interval is taken
# from the reviewer-sensitivity run so that the same formal bootstrap appears
# everywhere in the generated tables and validation assets.
empirical$n_transitions[4L] <- focal_empirical$n_transitions
empirical$n_failures[4L] <- focal_empirical$n_outcomes
empirical$log_odds_coefficient[4L] <- focal_empirical$log_odds_coefficient
empirical$odds_ratio[4L] <- focal_empirical$odds_ratio
empirical$ci_low[4L] <- focal_empirical$ci_low
empirical$ci_high[4L] <- focal_empirical$ci_high
table_7 <- data.frame(
  Model = c("Naive", "History", "Learner-item", "Plan-conditioned"),
  `N/fail.` = paste0(formatC(empirical$n_transitions, format = "f", digits = 0, big.mark = ","),
                     "/", formatC(empirical$n_failures, format = "f", digits = 0, big.mark = ",")),
  OR = fmt(empirical$odds_ratio, 3L),
  `95% CI` = paste0(fmt(empirical$ci_low, 2L, FALSE), "–", fmt(empirical$ci_high, 2L, FALSE)),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(table_7, "Table 7. Empirical estimator audit",
            "table_7_empirical_audit", c("left", "right", "right", "left"))

# Supplementary Table S1: calibration search ranking.
calibration_candidates <- brm_read_csv(
  file.path(output_dir, "calibration_search_candidates.csv")
)
calibration_candidates <- calibration_candidates[
  calibration_candidates$rank <= 5,
  , drop = FALSE
]
calibration_candidates <- calibration_candidates[order(calibration_candidates$rank), ]
supplement_s1 <- data.frame(
  Rank = calibration_candidates$rank,
  Adaptivity = fmt(calibration_candidates$adaptivity, 2L),
  `Deviation SD` = fmt(calibration_candidates$natural_deviation_sd, 2L),
  `Outcome intercept` = fmt(calibration_candidates$outcome_intercept, 2L),
  Replications = calibration_candidates$n_replications,
  `Failure rate` = fmt(calibration_candidates$mean_failure_rate, 4L, FALSE),
  Spearman = fmt(calibration_candidates$mean_scheduled_actual_spearman, 4L, FALSE),
  `Naive log OR` = fmt(calibration_candidates$mean_naive_log_odds, 4L),
  Loss = fmt(calibration_candidates$loss, 3L),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s1,
  "Supplementary Table S1. Top five empirical-signature calibration candidates",
  "supplementary_table_s1_calibration_search"
)
write_panel_markdown(
  "Supplementary Table S1. Top five empirical-signature calibration candidates",
  list(
    supplement_s1[c("Rank", "Adaptivity", "Deviation SD", "Outcome intercept", "Loss")],
    supplement_s1[c("Rank", "Replications", "Failure rate", "Spearman", "Naive log OR")]
  ),
  c("Panel A. Candidate parameters and scaled loss",
    "Panel B. Recovered empirical signatures"),
  "supplementary_table_s1_calibration_search.md",
  note = paste0(
    "Candidates were ranked by the prespecified scaled loss; lower values indicate closer joint calibration. ",
    "Targets were failure rate .0973, planned-enacted Spearman correlation .8146, and naive log odds -.2158."
  )
)

# Supplementary Table S2: out-of-search robustness for the top candidates.
calibration_robustness <- brm_read_csv(
  file.path(output_dir, "calibration_robustness_summary.csv")
)
calibration_robustness <- calibration_robustness[
  calibration_robustness$calibration_rank <= 5 &
    calibration_robustness$estimator %in% c("naive", "schedule_conditioned", "oracle"),
  , drop = FALSE
]
calibration_robustness$order <-
  calibration_robustness$calibration_rank * 10 +
  match(calibration_robustness$estimator,
        c("naive", "schedule_conditioned", "oracle"))
calibration_robustness <- calibration_robustness[
  order(calibration_robustness$order), , drop = FALSE
]
supplement_s2 <- data.frame(
  Rank = calibration_robustness$calibration_rank,
  Estimator = c(naive = "Naive", schedule_conditioned = "Plan",
                oracle = "Oracle")[calibration_robustness$estimator],
  Truth = fmt(calibration_robustness$mean_truth, 3L, FALSE),
  Estimate = fmt(calibration_robustness$mean_estimate, 3L, FALSE),
  Bias = fmt(calibration_robustness$bias, 3L, FALSE),
  RMSE = fmt(calibration_robustness$rmse, 3L, FALSE),
  Coverage = fmt(calibration_robustness$coverage_95, 3L, FALSE),
  `Sign error` = fmt(calibration_robustness$sign_error_rate, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s2,
  "Supplementary Table S2. Out-of-search robustness across the top five calibration candidates",
  "supplementary_table_s2_calibration_robustness"
)
write_panel_markdown(
  "Supplementary Table S2. Out-of-search robustness across the top five calibration candidates",
  list(supplement_s2[c("Rank", "Estimator", "Truth", "Estimate", "Bias")],
       supplement_s2[c("Rank", "Estimator", "RMSE", "Coverage", "Sign error")]),
  c("Panel A. Effect recovery", "Panel B. Error and interval performance"),
  "supplementary_table_s2_calibration_robustness.md",
  note = "Plan denotes the schedule-conditioned estimator. Robustness replications used an independent post-search stream."
)

# Supplementary Table S3: observation-equation misspecification, restricted to
# the candidate model family matching the true learner family and aggregated
# over the two scheduler families.
observation <- brm_read_csv(
  file.path(output_dir, "state_model_observation_misspec_summary.csv")
)
observation <- observation[
  as.character(observation$candidate_family) == as.character(observation$learner_family),
  , drop = FALSE
]
observation_key <- paste(observation$learner_family,
                         observation$observation_condition, sep = "\r")
observation_groups <- split(seq_len(nrow(observation)), observation_key)
observation_aggregated <- do.call(rbind, lapply(observation_groups, function(indices) {
  group <- observation[indices, , drop = FALSE]
  weights <- group$n_replications
  data.frame(
    learner_family = group$learner_family[1L],
    observation_condition = group$observation_condition[1L],
    n_replications = sum(weights),
    family_recovery_rate = weighted_mean_safe(group$family_recovery_rate, weights),
    parameter_exact_recovery_rate = weighted_mean_safe(
      group$parameter_exact_recovery_rate, weights
    ),
    mean_state_rmse = weighted_mean_safe(group$mean_state_rmse, weights),
    mean_state_correlation = weighted_mean_safe(group$mean_state_correlation, weights),
    mean_unsigned_pe_correlation = weighted_mean_safe(
      group$mean_unsigned_pe_correlation, weights
    ),
    stringsAsFactors = FALSE
  )
}))
observation_condition_labels <- c(
  correct = "Correct",
  intercept_minus_0.25 = "Intercept -0.25",
  intercept_plus_0.25 = "Intercept +0.25",
  state_slope_0.80 = "State slope ×0.80",
  combined_plus_0.25_slope_0.80 = "Intercept +0.25; slope ×0.80"
)
observation_aggregated$order <-
  match(observation_aggregated$learner_family, c("RL", "Bayesian")) * 10 +
  match(observation_aggregated$observation_condition,
        names(observation_condition_labels))
observation_aggregated <- observation_aggregated[
  order(observation_aggregated$order), , drop = FALSE
]
condition_label <- observation_condition_labels[
  as.character(observation_aggregated$observation_condition)
]
missing_condition_label <- is.na(condition_label)
condition_label[missing_condition_label] <-
  as.character(observation_aggregated$observation_condition[missing_condition_label])
supplement_s3 <- data.frame(
  `True learner` = observation_aggregated$learner_family,
  Condition = unname(condition_label),
  Replications = observation_aggregated$n_replications,
  `Family recovery` = fmt(observation_aggregated$family_recovery_rate, 3L, FALSE),
  `Exact parameter` = fmt(observation_aggregated$parameter_exact_recovery_rate, 3L, FALSE),
  `State RMSE` = fmt(observation_aggregated$mean_state_rmse, 3L, FALSE),
  `State r` = fmt(observation_aggregated$mean_state_correlation, 3L, FALSE),
  `PE r` = fmt(observation_aggregated$mean_unsigned_pe_correlation, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s3,
  "Supplementary Table S3. Recovery under observation-equation misspecification",
  "supplementary_table_s3_observation_misspecification"
)
write_panel_markdown(
  "Supplementary Table S3. Recovery under observation-equation misspecification",
  list(supplement_s3[c("True learner", "Condition", "Replications", "Family recovery", "Exact parameter")],
       supplement_s3[c("True learner", "Condition", "State RMSE", "State r", "PE r")]),
  c("Panel A. Family and parameter recovery", "Panel B. Latent-signal recovery"),
  "supplementary_table_s3_observation_misspecification.md",
  note = "Rows retain the candidate family matching the true learner and average over the two scheduler families. PE denotes unsigned prediction error."
)

# Supplementary Table S4: measured-context proxy gradient.
context_proxy <- brm_read_csv(file.path(output_dir, "context_proxy_summary.csv"))
context_proxy$order <-
  context_proxy$proxy_r2_target * 10 +
  match(context_proxy$estimator,
        c("plan_conditioned", "plan_plus_proxy", "oracle"))
context_proxy <- context_proxy[order(context_proxy$order), , drop = FALSE]
proxy_estimator_labels <- c(
  plan_conditioned = "Plan only",
  plan_plus_proxy = "Plan + proxy",
  oracle = "Oracle"
)
supplement_s4 <- data.frame(
  `Target R2` = fmt(context_proxy$proxy_r2_target, 2L, FALSE),
  `Realized R2` = fmt(context_proxy$mean_realized_proxy_r2, 3L, FALSE),
  Estimator = unname(proxy_estimator_labels[context_proxy$estimator]),
  Truth = fmt(context_proxy$mean_truth, 3L, FALSE),
  Estimate = fmt(context_proxy$mean_estimate, 3L, FALSE),
  Bias = fmt(context_proxy$bias, 3L, FALSE),
  RMSE = fmt(context_proxy$rmse, 3L, FALSE),
  `Normal cov.` = fmt(context_proxy$normal_coverage_95, 3L, FALSE),
  `Participant-t cov.` = fmt(context_proxy$participant_t_coverage_95, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s4,
  "Supplementary Table S4. Context-proxy quality and delay-effect recovery",
  "supplementary_table_s4_context_proxy"
)
write_panel_markdown(
  "Supplementary Table S4. Context-proxy quality and delay-effect recovery",
  list(supplement_s4[c("Target R2", "Realized R2", "Estimator", "Truth", "Estimate", "Bias")],
       supplement_s4[c("Target R2", "Estimator", "RMSE", "Normal cov.", "Participant-t cov.")]),
  c("Panel A. Effect recovery", "Panel B. Error and interval performance"),
  "supplementary_table_s4_context_proxy.md",
  note = "R2 is the squared correlation between the proxy and the otherwise unlogged standardized context."
)

# Supplementary Table S5: complete covariate-adjusted WCLS design grid.
wcls_adjusted <- wcls[wcls$estimator == "covariate_adjusted_itt", , drop = FALSE]
wcls_adjusted$order <-
  wcls_adjusted$n_participants * 100 +
  wcls_adjusted$perturbation_log_days * 10 +
  wcls_adjusted$perturbation_compliance
wcls_adjusted <- wcls_adjusted[order(wcls_adjusted$order), , drop = FALSE]
wcls_labels <- c(unadjusted_itt = "Unadjusted",
                 covariate_adjusted_itt = "Covariate adjusted")
supplement_s5 <- data.frame(
  N = wcls_adjusted$n_participants,
  Delta = fmt(wcls_adjusted$perturbation_log_days, 2L, FALSE),
  Compliance = fmt(wcls_adjusted$perturbation_compliance, 2L, FALSE),
  Truth = fmt(wcls_adjusted$mean_truth, 4L, FALSE),
  Estimate = fmt(wcls_adjusted$mean_estimate, 4L, FALSE),
  Bias = fmt(wcls_adjusted$bias, 4L, FALSE),
  RMSE = fmt(wcls_adjusted$rmse, 4L, FALSE),
  `Normal cov.` = fmt(wcls_adjusted$normal_coverage_95, 3L, FALSE),
  `Participant-t cov.` = fmt(wcls_adjusted$participant_t_coverage_95, 3L, FALSE),
  `Normal power` = fmt(wcls_adjusted$normal_power_positive, 3L, FALSE),
  `Participant-t power` = fmt(wcls_adjusted$participant_t_power_positive, 3L, FALSE),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s5,
  "Supplementary Table S5. Complete covariate-adjusted weighted-and-centered least-squares design grid",
  "supplementary_table_s5_wcls_grid"
)
write_panel_markdown(
  "Supplementary Table S5. Complete covariate-adjusted weighted-and-centered least-squares design grid",
  list(supplement_s5[c("N", "Delta", "Compliance", "Truth", "Estimate", "Bias")],
       supplement_s5[c("N", "Delta", "Compliance", "RMSE", "Normal cov.", "Participant-t cov.")],
       supplement_s5[c("N", "Delta", "Compliance", "Normal power", "Participant-t power")]),
  c("Panel A. Effect recovery", "Panel B. Error and coverage", "Panel C. Detection probability"),
  "supplementary_table_s5_wcls_grid.md",
  note = "Delta is the randomized perturbation in log days. Participant-t intervals use N - 1 degrees of freedom."
)

# Supplementary Table S6: algebraic/numerical equivalence of fixed-probability
# WCLS and the corresponding centered-treatment linear probability model.
wcls_equivalence <- brm_read_csv(
  file.path(output_dir, "micro_randomized_wcls_equivalence_summary.csv")
)
wcls_equivalence <- wcls_equivalence[
  match(c("unadjusted_itt", "covariate_adjusted_itt"), wcls_equivalence$estimator),
  , drop = FALSE
]
supplement_s6 <- data.frame(
  Estimator = unname(wcls_labels[wcls_equivalence$estimator]),
  Fits = wcls_equivalence$n_fits,
  `Max abs estimate diff.` = fmt_sci(
    wcls_equivalence$maximum_absolute_estimate_difference, 2L
  ),
  `Mean abs estimate diff.` = fmt_sci(
    wcls_equivalence$mean_absolute_estimate_difference, 2L
  ),
  `Max abs SE diff.` = fmt_sci(
    wcls_equivalence$maximum_absolute_standard_error_difference, 2L
  ),
  `Mean abs SE diff.` = fmt_sci(
    wcls_equivalence$mean_absolute_standard_error_difference, 2L
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s6,
  "Supplementary Table S6. Numerical equivalence of WCLS and the centered-treatment linear probability model",
  "supplementary_table_s6_wcls_equivalence"
)

# Supplementary Table S7: full empirical threshold/outcome/penalty grid.
empirical_sensitivity$order <-
  empirical_sensitivity$planned_interval_threshold_days * 100 +
  match(empirical_sensitivity$outcome_definition, c("Again", "Again or Hard")) * 10 +
  match(empirical_sensitivity$regularization_C, c(1, 10, 100))
empirical_sensitivity <- empirical_sensitivity[
  order(empirical_sensitivity$order), , drop = FALSE
]
supplement_s7 <- data.frame(
  `Plan threshold` = fmt(
    empirical_sensitivity$planned_interval_threshold_days, 2L, FALSE
  ),
  Outcome = sub("^Again or Hard$", "Again-or-Hard",
                empirical_sensitivity$outcome_definition),
  `C_ridge` = empirical_sensitivity$regularization_C,
  `N/events` = paste0(
    formatC(empirical_sensitivity$n_transitions, format = "f", digits = 0,
            big.mark = ","), "/",
    formatC(empirical_sensitivity$n_outcomes, format = "f", digits = 0,
            big.mark = ",")
  ),
  OR = fmt(empirical_sensitivity$odds_ratio, 3L),
  `95% bootstrap CI` = paste0(
    fmt(empirical_sensitivity$ci_low, 3L), "-",
    fmt(empirical_sensitivity$ci_high, 3L)
  ),
  `Successful bootstraps` = empirical_sensitivity$n_bootstrap_successful,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_table(
  supplement_s7,
  "Supplementary Table S7. Empirical outcome, threshold, and regularization sensitivity grid",
  "supplementary_table_s7_empirical_sensitivity"
)

# Supplementary Table S8: participant influence and mixed-model benchmark.
leave_one_out <- brm_read_csv(
  file.path(output_dir, "empirical_leave_one_participant_out.csv")
)
loo_groups <- split(
  seq_len(nrow(leave_one_out)),
  factor(leave_one_out$outcome_definition,
         levels = c("Again", "Again or Hard"))
)
loo_panel <- do.call(rbind, lapply(loo_groups, function(indices) {
  group <- leave_one_out[indices, , drop = FALSE]
  data.frame(
    Outcome = sub("^Again or Hard$", "Again-or-Hard",
                  group$outcome_definition[1L]),
    `Minimum OR` = fmt(min(group$odds_ratio), 3L),
    `Maximum OR` = fmt(max(group$odds_ratio), 3L),
    `Below 1` = paste0(sum(group$odds_ratio < 1), "/", nrow(group)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}))
mixed <- brm_read_csv(file.path(output_dir, "empirical_mixed_effects_benchmark.csv"))
mixed_panel <- data.frame(
  Outcome = sub("^Again or Hard$", "Again-or-Hard",
                mixed$outcome_definition),
  `N/events` = paste0(
    formatC(mixed$n_transitions, format = "f", digits = 0, big.mark = ","),
    "/", formatC(mixed$n_outcomes, format = "f", digits = 0, big.mark = ",")
  ),
  OR = fmt(mixed$odds_ratio, 3L),
  `95% Wald CI` = paste0(fmt(mixed$ci_low, 3L), "-", fmt(mixed$ci_high, 3L)),
  p = fmt(mixed$p_value, 3L, FALSE),
  Converged = ifelse(brm_match_bool(mixed$converged, TRUE), "Yes", "No"),
  Singular = ifelse(brm_match_bool(mixed$singular, TRUE), "Yes", "No"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
supplement_s8_csv <- rbind(
  data.frame(
    Panel = "Leave-one-participant-out", Outcome = loo_panel$Outcome,
    N = "", OR = "", `95% CI` = "", p = "",
    `LOO minimum OR` = loo_panel$`Minimum OR`,
    `LOO maximum OR` = loo_panel$`Maximum OR`,
    `LOO below 1` = loo_panel$`Below 1`,
    Converged = "", Singular = "", check.names = FALSE,
    stringsAsFactors = FALSE
  ),
  data.frame(
    Panel = "Crossed random-intercept model", Outcome = mixed_panel$Outcome,
    N = mixed_panel$`N/events`, OR = mixed_panel$OR,
    `95% CI` = mixed_panel$`95% Wald CI`, p = mixed_panel$p,
    `LOO minimum OR` = "", `LOO maximum OR` = "", `LOO below 1` = "",
    Converged = mixed_panel$Converged, Singular = mixed_panel$Singular,
    check.names = FALSE, stringsAsFactors = FALSE
  )
)
write.csv(
  supplement_s8_csv,
  file.path(table_dir, "supplementary_table_s8_empirical_influence_mixed.csv"),
  row.names = FALSE
)
write_panel_markdown(
  "Supplementary Table S8. Participant influence and crossed random-intercept benchmarks",
  list(loo_panel, mixed_panel),
  c("Panel A. Leave-one-participant-out ridge-logistic estimates",
    "Panel B. Crossed participant-card random-intercept logistic models"),
  "supplementary_table_s8_empirical_influence_mixed.md",
  note = paste0(
    "Both panels use the .5-day threshold. The inverse ridge-penalty constant is 10; ",
    "the mixed models report Wald intervals and were fit with crossed learner and card intercepts."
  )
)
message("R table: supplementary_table_s8_empirical_influence_mixed [4 rows]")

all_markdown <- unlist(lapply(sprintf("table_%d", 1:7), function(prefix) {
  file <- list.files(table_dir, pattern = paste0("^", prefix, "_.*\\.md$"), full.names = TRUE)
  c(readLines(file, warn = FALSE), "")
}))
writeLines(all_markdown, file.path(table_dir, "all_tables.md"), useBytes = TRUE)
message("R generated seven main and eight supplementary publication tables in ", table_dir)
