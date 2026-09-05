project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source(file.path(project_dir, "src", "R", "theme_brm.R"))
source(file.path(project_dir, "src", "R", "validate_public_data.R"))

data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")
validation_dir <- file.path(project_dir, "validation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)

# Independent seed reserved for the reviewer-requested empirical sensitivity run.
seed <- 20260807L
bootstrap_replications <- as.integer(Sys.getenv(
  "BRM_EMPIRICAL_SENSITIVITY_REPS", unset = "1200"
))
bootstrap_cores <- as.integer(Sys.getenv(
  "BRM_EMPIRICAL_SENSITIVITY_CORES",
  unset = if (.Platform$OS.type == "windows") "1" else "4"
))
if (!is.finite(bootstrap_replications) || bootstrap_replications < 1L) {
  stop("BRM_EMPIRICAL_SENSITIVITY_REPS must be a positive integer.")
}
if (!is.finite(bootstrap_cores) || bootstrap_cores < 1L) bootstrap_cores <- 1L

standardize <- function(x) {
  x <- as.numeric(x)
  deviation <- stats::sd(x)
  if (!is.finite(deviation) || deviation <= 0) return(rep(0, length(x)))
  (x - mean(x)) / deviation
}

design_matrix <- function(frame) {
  columns <- list(
    `(Intercept)` = rep(1, nrow(frame)),
    log_actual_gap = standardize(frame$log_actual_gap),
    log_planned_gap = standardize(frame$log_planned_gap),
    log_review_index = standardize(frame$log_review_index),
    study_day = standardize(frame$study_day),
    current_failure = standardize(frame$current_failure)
  )
  for (name in c("last_ease", "participant_code", "card_code", "scheduler_encoding")) {
    values <- as.character(frame[[name]])
    levels <- sort(unique(values))
    if (length(levels) > 1L) {
      for (level in levels[-1L]) {
        columns[[paste0(name, "_", level)]] <- as.numeric(values == level)
      }
    }
  }
  matrix <- do.call(cbind, columns)
  storage.mode(matrix) <- "double"
  matrix
}

ridge_logistic <- function(x, y, case_weights = NULL, C = 10,
                           max_iter = 80L, tolerance = 1e-8,
                           initial_beta = NULL) {
  if (is.null(case_weights)) case_weights <- rep(1, length(y))
  keep <- is.finite(case_weights) & case_weights > 0
  x <- x[keep, , drop = FALSE]
  y <- as.numeric(y[keep])
  case_weights <- as.numeric(case_weights[keep])
  if (length(unique(y)) < 2L) return(rep(NA_real_, ncol(x)))
  beta <- if (is.null(initial_beta)) rep(0, ncol(x)) else as.numeric(initial_beta)
  if (length(beta) != ncol(x) || any(!is.finite(beta))) beta <- rep(0, ncol(x))
  penalty <- diag(1 / C, ncol(x))
  for (iteration in seq_len(max_iter)) {
    eta <- drop(x %*% beta)
    probability <- plogis(pmax(-30, pmin(30, eta)))
    variance <- pmax(1e-8, probability * (1 - probability))
    working_weight <- variance * case_weights
    working_response <- eta + (y - probability) / variance
    weighted_x <- x * sqrt(working_weight)
    information <- crossprod(weighted_x) + penalty
    score_target <- crossprod(weighted_x, working_response * sqrt(working_weight))
    candidate <- tryCatch({
      cholesky <- chol(information)
      drop(backsolve(cholesky, forwardsolve(t(cholesky), score_target)))
    },
      error = function(condition) {
        tryCatch(
          drop(solve(information, score_target)),
          error = function(condition) drop(qr.solve(information, score_target, tol = 1e-10))
        )
      })
    if (max(abs(candidate - beta)) < tolerance) {
      beta <- candidate
      break
    }
    beta <- candidate
  }
  beta
}

make_bootstrap_weights <- function(participants, n_replications, cell_seed) {
  participant_levels <- sort(unique(as.character(participants)))
  set.seed(cell_seed)
  replicate(n_replications, {
    sampled <- sample.int(length(participant_levels), length(participant_levels), replace = TRUE)
    tabulate(sampled, nbins = length(participant_levels))
  })
}

fit_bootstrap_grid <- function(frame, outcome_label, threshold_label,
                               C_values, weight_matrix) {
  y <- if (outcome_label == "Again") {
    as.integer(frame$next_first_ease == 1)
  } else {
    as.integer(frame$next_first_ease <= 2)
  }
  x <- design_matrix(frame)
  focal_index <- match("log_actual_gap", colnames(x))
  participant_levels <- sort(unique(as.character(frame$participant_code)))
  participant_index <- match(as.character(frame$participant_code), participant_levels)
  point_estimates <- numeric(length(C_values))
  point_beta <- NULL
  for (C_index in seq_along(C_values)) {
    point_beta <- ridge_logistic(
      x, y, C = C_values[C_index], initial_beta = point_beta
    )
    point_estimates[C_index] <- point_beta[focal_index]
  }
  fit_one_resample <- function(replication) {
    observation_weights <- weight_matrix[participant_index, replication]
    estimates <- numeric(length(C_values))
    beta <- NULL
    for (C_index in seq_along(C_values)) {
      beta <- ridge_logistic(
        x, y, case_weights = observation_weights, C = C_values[C_index],
        initial_beta = beta
      )
      estimates[C_index] <- beta[focal_index]
    }
    estimates
  }
  replication_ids <- seq_len(ncol(weight_matrix))
  if (.Platform$OS.type != "windows" && bootstrap_cores > 1L) {
    fitted_list <- parallel::mclapply(
      replication_ids, fit_one_resample,
      mc.cores = bootstrap_cores, mc.preschedule = TRUE
    )
  } else {
    fitted_list <- lapply(replication_ids, fit_one_resample)
  }
  estimate_matrix <- do.call(rbind, fitted_list)
  if (length(C_values) == 1L) {
    estimate_matrix <- matrix(estimate_matrix, ncol = 1L)
  }
  summaries <- vector("list", length(C_values))
  replications <- vector("list", length(C_values))

  for (C_index in seq_along(C_values)) {
    C_value <- C_values[C_index]
    estimates <- estimate_matrix[, C_index]
    successful <- is.finite(estimates)
    interval <- stats::quantile(
      estimates[successful], c(0.025, 0.975), names = FALSE, type = 7
    )
    summaries[[C_index]] <- data.frame(
      planned_interval_threshold_days = threshold_label,
      outcome_definition = outcome_label,
      regularization_C = C_value,
      n_transitions = nrow(frame),
      n_outcomes = sum(y),
      log_odds_coefficient = point_estimates[C_index],
      odds_ratio = exp(point_estimates[C_index]),
      ci_low = exp(interval[1L]),
      ci_high = exp(interval[2L]),
      n_bootstrap_requested = ncol(weight_matrix),
      n_bootstrap_successful = sum(successful),
      stringsAsFactors = FALSE
    )
    replications[[C_index]] <- data.frame(
      planned_interval_threshold_days = threshold_label,
      outcome_definition = outcome_label,
      regularization_C = C_value,
      replication = seq_len(ncol(weight_matrix)),
      log_odds_coefficient = estimates,
      odds_ratio = exp(estimates),
      successful = successful,
      stringsAsFactors = FALSE
    )
  }
  list(
    summary = do.call(rbind, summaries),
    replications = do.call(rbind, replications)
  )
}

# The public workflow begins with the minimized cross-day transition table.
# Raw exports, exact within-day timestamps, original identifiers, and the private
# event-to-transition de-identification step are intentionally not distributed.
public_data_audit <- brm_validate_public_data(project_dir)
transitions <- public_data_audit$data
public_data_validation <- public_data_audit$checks
required_columns <- public_data_audit$required_fields

# Sensitivity grid. The three positive thresholds select the same records because
# the only sub-day scheduled intervals in this export are 60, 330, and 600 s.
transitions$log_actual_gap <- log1p(transitions$next_gap_days)
transitions$log_planned_gap <- log1p(pmax(transitions$last_scheduled_interval_days, 0))
transitions$log_review_index <- log1p(transitions$daily_review_index)
transitions$current_failure <- as.integer(transitions$first_ease == 1)

thresholds <- c(0, 0.25, 0.5, 1)
C_values <- c(1, 10, 100)
outcomes <- c("Again", "Again or Hard")
summary_parts <- list()
replication_parts <- list()
part_index <- 0L

# Thresholds .25, .50, and 1.00 have identical inclusion. Common random numbers
# make this equivalence visible rather than adding irrelevant Monte Carlo noise.
unique_threshold_groups <- list(`0` = 0, positive = c(0.25, 0.5, 1))
# The same resampled participant multiplicities are used in every cell. These
# common random numbers make sensitivity contrasts more precise and ensure that
# the predeclared reviewer-sensitivity seed is the only bootstrap seed.
common_bootstrap_weights <- make_bootstrap_weights(
  transitions$participant_code, bootstrap_replications, seed
)
for (group_index in seq_along(unique_threshold_groups)) {
  threshold_group <- unique_threshold_groups[[group_index]]
  representative_threshold <- threshold_group[1L]
  frame <- transitions[
    transitions$last_scheduled_interval_days >= representative_threshold,
    , drop = FALSE
  ]
  participant_levels <- sort(unique(as.character(frame$participant_code)))
  for (outcome_index in seq_along(outcomes)) {
    fitted <- fit_bootstrap_grid(
      frame, outcomes[outcome_index], representative_threshold, C_values,
      common_bootstrap_weights
    )
    for (threshold in threshold_group) {
      part_index <- part_index + 1L
      threshold_summary <- fitted$summary
      threshold_replications <- fitted$replications
      threshold_summary$planned_interval_threshold_days <- threshold
      threshold_replications$planned_interval_threshold_days <- threshold
      summary_parts[[part_index]] <- threshold_summary
      replication_parts[[part_index]] <- threshold_replications
    }
  }
}
sensitivity_summary <- do.call(rbind, summary_parts)
sensitivity_replications <- do.call(rbind, replication_parts)
sensitivity_summary <- sensitivity_summary[order(
  sensitivity_summary$planned_interval_threshold_days,
  sensitivity_summary$outcome_definition,
  sensitivity_summary$regularization_C
), , drop = FALSE]
sensitivity_replications <- sensitivity_replications[order(
  sensitivity_replications$planned_interval_threshold_days,
  sensitivity_replications$outcome_definition,
  sensitivity_replications$regularization_C,
  sensitivity_replications$replication
), , drop = FALSE]
row.names(sensitivity_summary) <- NULL
row.names(sensitivity_replications) <- NULL
write.csv(sensitivity_summary, file.path(output_dir, "empirical_sensitivity_summary.csv"),
          row.names = FALSE)
write.csv(sensitivity_replications,
          file.path(output_dir, "empirical_sensitivity_bootstrap_replications.csv"),
          row.names = FALSE)

# Leave-one-participant-out analysis for the focal threshold and penalty.
focal_frame <- transitions[transitions$last_scheduled_interval_days >= 0.5, , drop = FALSE]
loo_rows <- list()
loo_index <- 0L
for (outcome_label in outcomes) {
  for (omitted in sort(unique(as.character(focal_frame$participant_code)))) {
    frame <- focal_frame[focal_frame$participant_code != omitted, , drop = FALSE]
    y <- if (outcome_label == "Again") {
      as.integer(frame$next_first_ease == 1)
    } else {
      as.integer(frame$next_first_ease <= 2)
    }
    x <- design_matrix(frame)
    beta <- ridge_logistic(x, y, C = 10)
    estimate <- beta[match("log_actual_gap", colnames(x))]
    loo_index <- loo_index + 1L
    loo_rows[[loo_index]] <- data.frame(
      planned_interval_threshold_days = 0.5,
      outcome_definition = outcome_label,
      regularization_C = 10,
      omitted_participant = omitted,
      n_transitions = nrow(frame),
      n_outcomes = sum(y),
      log_odds_coefficient = estimate,
      odds_ratio = exp(estimate),
      stringsAsFactors = FALSE
    )
  }
}
leave_one_out <- do.call(rbind, loo_rows)
write.csv(leave_one_out,
          file.path(output_dir, "empirical_leave_one_participant_out.csv"),
          row.names = FALSE)

# Crossed random-intercept benchmark. It is deliberately a benchmark rather than
# a replacement for the ridge sensitivity grid because only 12 learner clusters
# are available.
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required for the crossed random-intercept benchmark.")
}
mixed_frame <- focal_frame
mixed_frame$z_log_actual_gap <- standardize(mixed_frame$log_actual_gap)
mixed_frame$z_log_planned_gap <- standardize(mixed_frame$log_planned_gap)
mixed_frame$z_log_review_index <- standardize(mixed_frame$log_review_index)
mixed_frame$z_study_day <- standardize(mixed_frame$study_day)
mixed_frame$z_current_failure <- standardize(mixed_frame$current_failure)
mixed_frame$last_ease_factor <- factor(mixed_frame$last_ease)
mixed_frame$scheduler_encoding <- factor(mixed_frame$scheduler_encoding)
mixed_rows <- list()
for (outcome_index in seq_along(outcomes)) {
  outcome_label <- outcomes[outcome_index]
  mixed_frame$outcome <- if (outcome_label == "Again") {
    as.integer(mixed_frame$next_first_ease == 1)
  } else {
    as.integer(mixed_frame$next_first_ease <= 2)
  }
  fit <- lme4::glmer(
    outcome ~ z_log_actual_gap + z_log_planned_gap + z_log_review_index +
      z_study_day + z_current_failure + last_ease_factor + scheduler_encoding +
      (1 | participant_code) + (1 | card_code),
    data = mixed_frame,
    family = stats::binomial(),
    control = lme4::glmerControl(
      optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
      calc.derivs = TRUE
    )
  )
  coefficient_table <- summary(fit)$coefficients
  estimate <- coefficient_table["z_log_actual_gap", "Estimate"]
  standard_error <- coefficient_table["z_log_actual_gap", "Std. Error"]
  p_value <- coefficient_table["z_log_actual_gap", "Pr(>|z|)"]
  messages <- fit@optinfo$conv$lme4$messages
  converged <- is.null(messages) && identical(fit@optinfo$conv$opt, 0L)
  singular <- lme4::isSingular(fit, tol = 1e-4)
  mixed_rows[[outcome_index]] <- data.frame(
    planned_interval_threshold_days = 0.5,
    outcome_definition = outcome_label,
    n_transitions = nrow(mixed_frame),
    n_outcomes = sum(mixed_frame$outcome),
    log_odds_coefficient = estimate,
    standard_error = standard_error,
    odds_ratio = exp(estimate),
    ci_low = exp(estimate - stats::qnorm(0.975) * standard_error),
    ci_high = exp(estimate + stats::qnorm(0.975) * standard_error),
    p_value = p_value,
    converged = converged,
    singular = singular,
    optimizer = "bobyqa",
    stringsAsFactors = FALSE
  )
}
mixed_benchmark <- do.call(rbind, mixed_rows)
write.csv(mixed_benchmark,
          file.path(output_dir, "empirical_mixed_effects_benchmark.csv"),
          row.names = FALSE)

# Targeted, machine-readable validation.
validation <- list()
add_validation <- function(name, passed, observed, expected) {
  validation[[length(validation) + 1L]] <<- data.frame(
    check = name, passed = isTRUE(passed), observed = as.character(observed),
    expected = as.character(expected), stringsAsFactors = FALSE
  )
}
add_validation(
  "all public-data checks pass",
  all(public_data_validation$passed),
  paste(sum(public_data_validation$passed), "of", nrow(public_data_validation)),
  paste(nrow(public_data_validation), "of", nrow(public_data_validation))
)
add_validation("sensitivity grid has 24 cells", nrow(sensitivity_summary) == 24L,
               nrow(sensitivity_summary), 24L)
add_validation("all cells retain requested bootstraps",
               all(sensitivity_summary$n_bootstrap_successful == bootstrap_replications),
               min(sensitivity_summary$n_bootstrap_successful), bootstrap_replications)
add_validation("replication file has 24 times bootstrap rows",
               nrow(sensitivity_replications) == 24L * bootstrap_replications,
               nrow(sensitivity_replications), 24L * bootstrap_replications)
focal_again <- sensitivity_summary[
  sensitivity_summary$planned_interval_threshold_days == 0.5 &
    sensitivity_summary$outcome_definition == "Again" &
    sensitivity_summary$regularization_C == 10, , drop = FALSE
]
focal_hard <- sensitivity_summary[
  sensitivity_summary$planned_interval_threshold_days == 0.5 &
    sensitivity_summary$outcome_definition == "Again or Hard" &
    sensitivity_summary$regularization_C == 10, , drop = FALSE
]
add_validation("focal Again point estimate", abs(focal_again$odds_ratio - 0.758613) < 1e-4,
               sprintf("%.6f", focal_again$odds_ratio), "0.758613 +/- 0.0001")
add_validation("focal Again-or-Hard point estimate", abs(focal_hard$odds_ratio - 0.780638) < 1e-4,
               sprintf("%.6f", focal_hard$odds_ratio), "0.780638 +/- 0.0001")
add_validation("fixed-seed focal bootstrap intervals",
               abs(focal_again$ci_low - 0.348155) < 1e-4 &
                 abs(focal_again$ci_high - 2.547348) < 1e-4 &
                 abs(focal_hard$ci_low - 0.654587) < 1e-4 &
                 abs(focal_hard$ci_high - 2.205107) < 1e-4,
               sprintf("Again [%.6f, %.6f]; Again-or-Hard [%.6f, %.6f]",
                       focal_again$ci_low, focal_again$ci_high,
                       focal_hard$ci_low, focal_hard$ci_high),
               "Again [0.348155, 2.547348]; Again-or-Hard [0.654587, 2.205107]")
loo_again <- leave_one_out$odds_ratio[leave_one_out$outcome_definition == "Again"]
loo_hard <- leave_one_out$odds_ratio[leave_one_out$outcome_definition == "Again or Hard"]
add_validation("Again leave-one-out range and direction count",
               abs(min(loo_again) - 0.589392) < 1e-4 &
                 abs(max(loo_again) - 1.341560) < 1e-4 & sum(loo_again < 1) == 11L,
               sprintf("%.6f to %.6f; %d/12 below 1", min(loo_again), max(loo_again),
                       sum(loo_again < 1)), "0.589392 to 1.341560; 11/12 below 1")
add_validation("Again-or-Hard leave-one-out range and direction count",
               abs(min(loo_hard) - 0.720925) < 1e-4 &
                 abs(max(loo_hard) - 0.983165) < 1e-4 & sum(loo_hard < 1) == 12L,
               sprintf("%.6f to %.6f; %d/12 below 1", min(loo_hard), max(loo_hard),
                       sum(loo_hard < 1)), "0.720925 to 0.983165; 12/12 below 1")
add_validation("mixed models converge", all(mixed_benchmark$converged),
               paste(mixed_benchmark$converged, collapse = ", "), "TRUE, TRUE")
add_validation("mixed models are non-singular", !any(mixed_benchmark$singular),
               paste(mixed_benchmark$singular, collapse = ", "), "FALSE, FALSE")
add_validation("mixed Again benchmark reproduces",
               abs(mixed_benchmark$odds_ratio[1L] - 0.811383) < 1e-4 &
                 abs(mixed_benchmark$ci_low[1L] - 0.604634) < 1e-4 &
                 abs(mixed_benchmark$ci_high[1L] - 1.088826) < 1e-4 &
                 abs(mixed_benchmark$p_value[1L] - 0.163662) < 1e-4,
               sprintf("OR %.6f [%.6f, %.6f], p %.6f",
                       mixed_benchmark$odds_ratio[1L], mixed_benchmark$ci_low[1L],
                       mixed_benchmark$ci_high[1L], mixed_benchmark$p_value[1L]),
               "OR 0.811383 [0.604634, 1.088826], p 0.163662")
add_validation("mixed Again-or-Hard benchmark reproduces",
               abs(mixed_benchmark$odds_ratio[2L] - 0.801807) < 1e-4 &
                 abs(mixed_benchmark$ci_low[2L] - 0.648192) < 1e-4 &
                 abs(mixed_benchmark$ci_high[2L] - 0.991826) < 1e-4 &
                 abs(mixed_benchmark$p_value[2L] - 0.041790) < 1e-4,
               sprintf("OR %.6f [%.6f, %.6f], p %.6f",
                       mixed_benchmark$odds_ratio[2L], mixed_benchmark$ci_low[2L],
                       mixed_benchmark$ci_high[2L], mixed_benchmark$p_value[2L]),
               "OR 0.801807 [0.648192, 0.991826], p 0.041790")
validation_table <- do.call(rbind, validation)
write.csv(validation_table,
          file.path(validation_dir, "empirical_sensitivity_validation.csv"),
          row.names = FALSE)

manifest_lines <- c(
  "{",
  '  "analysis": "Empirical sensitivity audit from minimized public transitions",',
  sprintf('  "analysis_date": "%s",', format(Sys.Date(), "%Y-%m-%d")),
  sprintf('  "random_seed": %d,', seed),
  '  "input_file": "data/analysis_transitions_deidentified.csv",',
  sprintf('  "public_transition_rows": %d,', nrow(transitions)),
  sprintf('  "public_data_fields": %d,', length(required_columns)),
  sprintf('  "participant_codes": %d,',
          length(unique(transitions$participant_code))),
  sprintf('  "card_codes": %d,', length(unique(transitions$card_code))),
  sprintf('  "again_failures": %d,', sum(transitions$next_failure)),
  sprintf('  "public_data_checks_passed": %d,',
          sum(public_data_validation$passed)),
  sprintf('  "public_data_checks_total": %d,',
          nrow(public_data_validation)),
  sprintf('  "sensitivity_cells": %d,', nrow(sensitivity_summary)),
  sprintf('  "bootstrap_replications_per_cell": %d,', bootstrap_replications),
  sprintf('  "all_validation_checks_passed": %s,',
          tolower(as.character(all(validation_table$passed)))),
  '  "engine": "R base; lme4 for crossed random-intercept benchmarks",',
  '  "outputs": [',
  paste(sprintf('    "%s"', c(
    "empirical_sensitivity_summary.csv",
    "empirical_sensitivity_bootstrap_replications.csv",
    "empirical_leave_one_participant_out.csv",
    "empirical_mixed_effects_benchmark.csv"
  )), collapse = ",\n"),
  "  ]",
  "}"
)
writeLines(manifest_lines,
           file.path(output_dir, "empirical_sensitivity_manifest.json"),
           useBytes = TRUE)

report <- c(
  "# Empirical sensitivity audit",
  "",
  sprintf("Generated: %s", format(Sys.Date(), "%Y-%m-%d")),
  "",
  "## Public analysis input",
  "",
  sprintf(paste0(
    "The public workflow begins with **%s de-identified cross-day ",
    "transitions** and **%s analysis fields**. The table contains %s ",
    "pseudonymous participant codes and %s pseudonymous card codes. All ",
    "**%s/%s public-data ",
    "checks passed**. Raw exports, absolute dates, exact within-day timestamps, ",
    "original identifiers, linkage keys, and the private de-identification step ",
    "are outside this repository."
  ), format(nrow(transitions), big.mark = ","), length(required_columns),
  length(unique(transitions$participant_code)),
  length(unique(transitions$card_code)), sum(public_data_validation$passed),
  nrow(public_data_validation)),
  "",
  "## Sensitivity design",
  "",
  sprintf(paste0(
    "The audit evaluated 4 planned-interval thresholds x 2 outcome definitions ",
    "x 3 L2 penalties = **%s cells**. Every cell retained **%s successful ",
    "participant-cluster bootstrap estimates**. The .25-, .50-, and 1-day ",
    "thresholds selected the same 1,724 transitions because the only observed ",
    "sub-day plans were 60, 330, and 600 seconds."
  ), nrow(sensitivity_summary), bootstrap_replications),
  sprintf(paste0(
    "All cells used the same participant-resample multiplicities (common random ",
    "numbers) from the independently reserved sensitivity seed %s."
  ), seed),
  "",
  "## Focal results",
  "",
  sprintf(paste0(
    "At the .5-day threshold and C = 10, the Again outcome gave OR = %.3f, ",
    "bootstrap 95%% interval [%.3f, %.3f]. The Again-or-Hard outcome gave ",
    "OR = %.3f [%.3f, %.3f]."
  ), focal_again$odds_ratio, focal_again$ci_low, focal_again$ci_high,
  focal_hard$odds_ratio, focal_hard$ci_low, focal_hard$ci_high),
  sprintf(paste0(
    "The leave-one-participant-out range was %.3f-%.3f for Again (%s/12 below 1) ",
    "and %.3f-%.3f for Again-or-Hard (%s/12 below 1)."
  ), min(loo_again), max(loo_again), sum(loo_again < 1),
  min(loo_hard), max(loo_hard), sum(loo_hard < 1)),
  sprintf(paste0(
    "The crossed random-intercept benchmarks gave OR = %.3f [%.3f, %.3f], ",
    "p = %.3f for Again and OR = %.3f [%.3f, %.3f], p = %.3f for ",
    "Again-or-Hard. Both models converged and were non-singular."
  ), mixed_benchmark$odds_ratio[1L], mixed_benchmark$ci_low[1L],
  mixed_benchmark$ci_high[1L], mixed_benchmark$p_value[1L],
  mixed_benchmark$odds_ratio[2L], mixed_benchmark$ci_low[2L],
  mixed_benchmark$ci_high[2L], mixed_benchmark$p_value[2L]),
  "",
  "## Interpretation boundary",
  "",
  paste0(
    "These observational estimates are association diagnostics. The stable ",
    "below-one direction for the broader outcome is meaningful evidence of a ",
    "descriptive relation after the stated controls, while adaptive scheduling, ",
    "unlogged learner state, the 12-cluster sample, and outcome-definition ",
    "dependence prevent a causal spacing-effect interpretation."
  ),
  "",
  "## Validation",
  "",
  sprintf("Targeted validation status: **%s/%s checks passed**.",
          sum(validation_table$passed), nrow(validation_table))
)
writeLines(report,
           file.path(validation_dir, "empirical_sensitivity_audit_report.md"),
           useBytes = TRUE)

if (!all(validation_table$passed)) {
  stop("Empirical sensitivity validation failed: ",
       paste(validation_table$check[!validation_table$passed], collapse = "; "))
}
message(
  "R empirical sensitivity audit complete: ", nrow(transitions),
  " public transitions; ", nrow(sensitivity_summary),
  " sensitivity cells with ", bootstrap_replications, " bootstraps each."
)
