project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source(file.path(project_dir, "src", "R", "theme_brm.R"))
output_dir <- file.path(project_dir, "outputs")

as_numeric_safe <- function(x) {
  if (is.logical(x) || is.numeric(x)) return(as.numeric(x))
  normalized <- tolower(trimws(as.character(x)))
  boolean <- normalized %in% c("true", "false", "t", "f")
  if (all(boolean | is.na(normalized))) {
    return(ifelse(normalized %in% c("true", "t"), 1, 0))
  }
  suppressWarnings(as.numeric(normalized))
}

mean_safe <- function(x) mean(as_numeric_safe(x), na.rm = TRUE)
sd_safe <- function(x) {
  x <- as_numeric_safe(x)
  if (sum(is.finite(x)) < 2L) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}

split_summarise <- function(data, keys, summariser) {
  key_frame <- data[keys]
  key_text <- do.call(paste, c(lapply(key_frame, as.character), sep = "\r"))
  groups <- split(seq_len(nrow(data)), factor(key_text, levels = unique(key_text)))
  rows <- lapply(groups, function(index) {
    identifiers <- data[index[1L], keys, drop = FALSE]
    values <- summariser(data[index, , drop = FALSE])
    cbind(identifiers, as.data.frame(values, check.names = FALSE), stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

write_summary <- function(data, filename) {
  path <- file.path(output_dir, filename)
  write.csv(data, path, row.names = FALSE, na = "")
  message("R summary: ", filename, " [", nrow(data), " rows]")
}

recovery_metrics <- function(group, include_spearman = TRUE) {
  n <- nrow(group)
  values <- list(
    n_replications = n,
    mean_estimate = mean_safe(group$estimate),
    mean_truth = mean_safe(group$truth),
    bias = mean_safe(group$bias),
    absolute_bias = mean_safe(abs(group$bias)),
    rmse = sqrt(mean_safe(group$squared_error)),
    coverage_95 = mean_safe(group$covered),
    sign_error_rate = mean_safe(group$sign_error),
    convergence_rate = mean_safe(group$converged),
    monte_carlo_se_bias = sd_safe(group$bias) / sqrt(n),
    mean_failure_rate = mean_safe(group$failure_rate)
  )
  if (include_spearman) {
    values$mean_planned_actual_spearman <- mean_safe(group$planned_actual_spearman)
  }
  values
}

parameter <- brm_read_csv(file.path(output_dir, "parameter_recovery_replications.csv"))
parameter_summary <- split_summarise(
  parameter,
  c("n_participants", "adaptivity", "state_reliability", "estimator"),
  recovery_metrics
)
write_summary(parameter_summary, "parameter_recovery_summary.csv")

context <- brm_read_csv(file.path(output_dir, "context_violation_replications.csv"))
context_summary <- split_summarise(
  context,
  c("n_participants", "context_strength", "estimator"),
  recovery_metrics
)
write_summary(context_summary, "context_violation_summary.csv")

mrt <- brm_read_csv(file.path(output_dir, "micro_randomized_replications.csv"))
mrt_metrics <- function(group) {
  n <- nrow(group)
  list(
    n_replications = n,
    mean_estimate = mean_safe(group$estimate),
    mean_truth = mean_safe(group$truth),
    bias = mean_safe(group$bias),
    absolute_bias = mean_safe(abs(group$bias)),
    rmse = sqrt(mean_safe(group$squared_error)),
    coverage_95 = mean_safe(group$covered),
    power_positive = mean_safe(group$detected_positive),
    monte_carlo_se_bias = sd_safe(group$bias) / sqrt(n),
    mean_failure_rate = mean_safe(group$failure_rate)
  )
}
mrt_summary <- split_summarise(
  mrt,
  c("n_participants", "perturbation_log_days", "perturbation_compliance", "estimator"),
  mrt_metrics
)
write_summary(mrt_summary, "micro_randomized_summary.csv")

wcls_path <- file.path(output_dir, "micro_randomized_wcls_replications.csv")
if (file.exists(wcls_path)) {
  wcls <- brm_read_csv(wcls_path)
  wcls_metrics <- function(group) {
    n <- nrow(group)
    list(
      n_replications = n,
      mean_truth = mean_safe(group$truth),
      mean_estimate = mean_safe(group$estimate),
      bias = mean_safe(group$bias),
      absolute_bias = mean_safe(abs(group$bias)),
      rmse = sqrt(mean_safe(group$squared_error)),
      normal_coverage_95 = mean_safe(group$normal_covered),
      participant_t_coverage_95 = mean_safe(group$participant_t_covered),
      normal_power_positive = mean_safe(group$normal_detected_positive),
      participant_t_power_positive = mean_safe(group$participant_t_detected_positive),
      monte_carlo_se_bias = sd_safe(group$bias) / sqrt(n),
      mean_standard_error = mean_safe(group$standard_error),
      mean_failure_rate = mean_safe(group$failure_rate)
    )
  }
  wcls_summary <- split_summarise(
    wcls,
    c("n_participants", "perturbation_log_days",
      "perturbation_compliance", "estimator"),
    wcls_metrics
  )
  write_summary(wcls_summary, "micro_randomized_wcls_summary.csv")

  wcls_equivalence <- split_summarise(
    wcls,
    c("estimator"),
    function(group) {
      list(
        n_fits = nrow(group),
        maximum_absolute_estimate_difference =
          max(as_numeric_safe(group$absolute_estimate_difference), na.rm = TRUE),
        mean_absolute_estimate_difference =
          mean_safe(group$absolute_estimate_difference),
        maximum_absolute_standard_error_difference =
          max(as_numeric_safe(group$absolute_standard_error_difference), na.rm = TRUE),
        mean_absolute_standard_error_difference =
          mean_safe(group$absolute_standard_error_difference)
      )
    }
  )
  write_summary(
    wcls_equivalence,
    "micro_randomized_wcls_equivalence_summary.csv"
  )
}

context_proxy_path <- file.path(output_dir, "context_proxy_replications.csv")
if (file.exists(context_proxy_path)) {
  context_proxy <- brm_read_csv(context_proxy_path)
  context_proxy_summary <- split_summarise(
    context_proxy,
    c("proxy_r2_target", "estimator"),
    function(group) {
      n <- nrow(group)
      list(
        n_replications = n,
        mean_realized_proxy_r2 = mean_safe(group$proxy_r2_realized),
        maximum_absolute_proxy_r2_error = max(
          abs(as_numeric_safe(group$proxy_r2_realized) -
              as_numeric_safe(group$proxy_r2_target)),
          na.rm = TRUE
        ),
        mean_truth = mean_safe(group$truth),
        mean_estimate = mean_safe(group$estimate),
        bias = mean_safe(group$bias),
        absolute_bias = mean_safe(abs(group$bias)),
        rmse = sqrt(mean_safe(group$squared_error)),
        normal_coverage_95 = mean_safe(group$normal_covered),
        participant_t_coverage_95 = mean_safe(group$participant_t_covered),
        sign_error_rate = mean_safe(group$sign_error),
        convergence_rate = mean_safe(group$converged),
        monte_carlo_se_bias = sd_safe(group$bias) / sqrt(n),
        mean_failure_rate = mean_safe(group$failure_rate),
        mean_planned_actual_spearman = mean_safe(group$planned_actual_spearman)
      )
    }
  )
  write_summary(context_proxy_summary, "context_proxy_summary.csv")
}

if (exists("wcls", inherits = FALSE) &&
    exists("context_proxy", inherits = FALSE)) {
  extension_manifest <- data.frame(
    component = c("micro_randomized_wcls", "context_proxy"),
    n_replication_rows = c(nrow(wcls), nrow(context_proxy)),
    n_design_cells = c(
      nrow(unique(wcls[c(
        "n_participants", "perturbation_log_days",
        "perturbation_compliance", "estimator"
      )])),
      nrow(unique(context_proxy[c("proxy_r2_target", "estimator")]))
    ),
    replications_per_cell = c(250L, 160L),
    participant_grid = c("12|30|60|100", "30"),
    seed_stream = c(
      "SEED+90000000+N*100000+delta*100000+compliance*10+replication",
      "SEED+70000000+N*100000+replication; proxy noise=data_seed+700000000"
    ),
    stringsAsFactors = FALSE
  )
  write_summary(extension_manifest, "wcls_context_proxy_manifest.csv")

  target <- wcls[
    as_numeric_safe(wcls$n_participants) == 100 &
      abs(as_numeric_safe(wcls$perturbation_log_days) - 0.30) < 1e-12 &
      abs(as_numeric_safe(wcls$perturbation_compliance) - 1.00) < 1e-12 &
      as.character(wcls$estimator) == "covariate_adjusted_itt",
    , drop = FALSE
  ]
  target_truth <- mean_safe(target$truth)
  target_estimate <- mean_safe(target$estimate)
  target_t_power <- mean_safe(target$participant_t_detected_positive)
  wcls_cell_counts <- table(
    wcls$n_participants,
    wcls$perturbation_log_days,
    wcls$perturbation_compliance,
    wcls$estimator
  )
  proxy_cell_counts <- table(
    context_proxy$proxy_r2_target,
    context_proxy$estimator
  )
  proxy_curve <- context_proxy_summary[
    as.character(context_proxy_summary$estimator) == "plan_plus_proxy",
    , drop = FALSE
  ]
  proxy_curve <- proxy_curve[order(as_numeric_safe(proxy_curve$proxy_r2_target)), ]
  max_estimate_difference <- max(
    as_numeric_safe(wcls$absolute_estimate_difference), na.rm = TRUE
  )
  max_se_difference <- max(
    as_numeric_safe(wcls$absolute_standard_error_difference), na.rm = TRUE
  )
  max_proxy_r2_error <- max(
    abs(as_numeric_safe(context_proxy$proxy_r2_realized) -
        as_numeric_safe(context_proxy$proxy_r2_target)),
    na.rm = TRUE
  )
  extension_validation <- data.frame(
    check = c(
      "WCLS replication row count",
      "Context-proxy replication row count",
      "WCLS cell replication counts",
      "Context-proxy cell replication counts",
      "WCLS/LPM estimate equivalence",
      "WCLS/LPM standard-error equivalence",
      "Finite-sample proxy R2 construction",
      "Declared WCLS target truth",
      "Declared WCLS target estimate",
      "Declared WCLS target participant-t power",
      "Proxy bias decreases with proxy R2",
      "Proxy participant-t coverage is nondecreasing",
      "All extension fits converged",
      "Participant-t interval contains normal interval"
    ),
    observed = c(
      nrow(wcls),
      nrow(context_proxy),
      paste(range(as.numeric(wcls_cell_counts)), collapse = ".."),
      paste(range(as.numeric(proxy_cell_counts)), collapse = ".."),
      format(max_estimate_difference, digits = 17),
      format(max_se_difference, digits = 17),
      format(max_proxy_r2_error, digits = 17),
      format(target_truth, digits = 17),
      format(target_estimate, digits = 17),
      format(target_t_power, digits = 17),
      paste(format(as_numeric_safe(proxy_curve$bias), digits = 6), collapse = "|"),
      paste(format(as_numeric_safe(proxy_curve$participant_t_coverage_95),
                   digits = 6), collapse = "|"),
      mean(c(as_numeric_safe(wcls$converged),
             as_numeric_safe(context_proxy$converged))),
      mean(
        as_numeric_safe(wcls$participant_t_ci_low) <=
          as_numeric_safe(wcls$normal_ci_low) + 1e-15 &
          as_numeric_safe(wcls$participant_t_ci_high) >=
            as_numeric_safe(wcls$normal_ci_high) - 1e-15
      )
    ),
    passed = c(
      nrow(wcls) == 12000L,
      nrow(context_proxy) == 2400L,
      all(wcls_cell_counts == 250L),
      all(proxy_cell_counts == 160L),
      max_estimate_difference < 1e-12,
      max_se_difference < 1e-12,
      max_proxy_r2_error < 1e-10,
      abs(target_truth - 0.014127105635882507) < 1e-12,
      abs(target_estimate - 0.013839562567328042) < 1e-12,
      abs(target_t_power - 0.872) < 1e-12,
      all(diff(as_numeric_safe(proxy_curve$bias)) < 0),
      all(diff(as_numeric_safe(proxy_curve$participant_t_coverage_95)) >= 0),
      all(as_numeric_safe(wcls$converged) == 1) &&
        all(as_numeric_safe(context_proxy$converged) == 1),
      all(
        as_numeric_safe(wcls$participant_t_ci_low) <=
          as_numeric_safe(wcls$normal_ci_low) + 1e-15 &
          as_numeric_safe(wcls$participant_t_ci_high) >=
            as_numeric_safe(wcls$normal_ci_high) - 1e-15
      )
    ),
    stringsAsFactors = FALSE
  )
  write_summary(extension_validation, "wcls_context_proxy_validation.csv")
  if (!all(extension_validation$passed)) {
    stop("WCLS/context-proxy targeted validation failed")
  }
}

delay <- brm_read_csv(file.path(output_dir, "state_model_delay_recovery_replications.csv"))
delay_summary <- split_summarise(
  delay,
  c("learner_family", "scheduler_family", "n_stages", "estimator"),
  function(group) recovery_metrics(group, include_spearman = FALSE)
)
delay_summary <- delay_summary[c(
  "learner_family", "scheduler_family", "n_stages", "estimator",
  "n_replications", "mean_truth", "mean_estimate", "bias", "absolute_bias",
  "rmse", "coverage_95", "sign_error_rate", "convergence_rate",
  "monte_carlo_se_bias", "mean_failure_rate"
)]
write_summary(delay_summary, "state_model_delay_recovery_summary.csv")

diagnostics <- brm_read_csv(file.path(output_dir, "state_model_diagnostics_replications.csv"))
family_summary <- split_summarise(
  diagnostics,
  c("learner_family", "n_stages"),
  function(group) {
    family <- as.character(group$learner_family[1L])
    truth <- if (family == "RL") 0.25 else 0.20
    parameter_column <- if (family == "RL") {
      group$rl_parameter_within_family
    } else {
      group$bayesian_parameter_within_family
    }
    delta_nll <- if (family == "RL") {
      group$bayesian_nll - group$rl_nll
    } else {
      group$rl_nll - group$bayesian_nll
    }
    list(
      n_datasets = nrow(group),
      family_recovery_rate = mean(as.character(group$selected_family) == family),
      true_parameter = truth,
      mean_parameter_within_true_family = mean_safe(parameter_column),
      parameter_bias_within_true_family = mean_safe(parameter_column) - truth,
      parameter_exact_recovery_rate = mean(abs(as.numeric(parameter_column) - truth) < 1e-10),
      mean_delta_nll_wrong_minus_true = mean_safe(delta_nll)
    )
  }
)
write_summary(family_summary, "state_model_family_recovery_summary.csv")

state_rows <- list()
row_index <- 1L
state_keys <- c("learner_family", "scheduler_family", "n_stages")
state_key_text <- do.call(paste, c(lapply(diagnostics[state_keys], as.character), sep = "\r"))
state_groups <- split(seq_len(nrow(diagnostics)), factor(state_key_text, levels = unique(state_key_text)))
for (indices in state_groups) {
  group <- diagnostics[indices, , drop = FALSE]
  for (candidate in c("RL", "Bayesian")) {
    prefix <- tolower(candidate)
    state_rows[[row_index]] <- data.frame(
      learner_family = group$learner_family[1L],
      scheduler_family = group$scheduler_family[1L],
      n_stages = group$n_stages[1L],
      candidate_family = candidate,
      learner_family_match = candidate == group$learner_family[1L],
      scheduler_family_match = candidate == group$scheduler_family[1L],
      state_rmse = mean_safe(group[[paste0(prefix, "_state_rmse")]]),
      state_correlation = mean_safe(group[[paste0(prefix, "_state_correlation")]]),
      unsigned_pe_correlation = mean_safe(group[[paste0(prefix, "_unsigned_pe_correlation")]]),
      policy_rmse = mean_safe(group[[paste0(prefix, "_policy_rmse")]]),
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
state_summary <- do.call(rbind, state_rows)
write_summary(state_summary, "state_model_state_recovery_summary.csv")

decode_estimator <- function(estimator, learner_family, scheduler_family) {
  if (estimator == "naive") {
    return(list(design = "Naive", learner_match = NA, policy_match = NA))
  }
  if (estimator == "oracle") {
    return(list(design = "Oracle", learner_match = TRUE, policy_match = TRUE))
  }
  if (startsWith(estimator, "state_")) {
    candidate <- sub("^state_", "", estimator)
    return(list(design = "State proxy only", learner_match = candidate == learner_family,
                policy_match = NA))
  }
  if (startsWith(estimator, "logged_plan_")) {
    candidate <- sub("^logged_plan_", "", estimator)
    return(list(design = "Logged plan", learner_match = candidate == learner_family,
                policy_match = NA))
  }
  pieces <- strsplit(sub("^dual_", "", estimator), "_", fixed = TRUE)[[1L]]
  list(
    design = "Reconstructed policy",
    learner_match = pieces[1L] == learner_family,
    policy_match = pieces[2L] == scheduler_family
  )
}

decoded <- mapply(
  decode_estimator,
  delay_summary$estimator,
  delay_summary$learner_family,
  delay_summary$scheduler_family,
  SIMPLIFY = FALSE
)
layer_summary <- delay_summary
layer_summary$design <- vapply(decoded, `[[`, character(1L), "design")
layer_summary$learner_match <- vapply(decoded, function(x) {
  if (is.na(x$learner_match)) NA else as.logical(x$learner_match)
}, logical(1L))
layer_summary$policy_match <- vapply(decoded, function(x) {
  if (is.na(x$policy_match)) NA else as.logical(x$policy_match)
}, logical(1L))
write_summary(layer_summary, "state_model_layer_summary.csv")

observation_path <- file.path(
  output_dir, "state_model_observation_misspec_replications.csv"
)
if (file.exists(observation_path)) {
  observation <- brm_read_csv(observation_path)
  observation$selected_model <- as.character(observation$candidate_family) ==
    as.character(observation$selected_family)
  observation_summary <- split_summarise(
    observation,
    c(
      "learner_family", "scheduler_family", "n_stages",
      "observation_condition", "intercept_shift",
      "state_slope_multiplier", "candidate_family"
    ),
    function(group) {
      n <- nrow(group)
      selected <- as_numeric_safe(group$selected_model) > 0
      true_parameter <- mean_safe(group$true_parameter)
      list(
        n_replications = n,
        family_recovery_rate = mean_safe(group$family_recovered),
        candidate_selection_rate = mean_safe(group$selected_model),
        true_parameter = true_parameter,
        mean_best_parameter = mean_safe(group$best_parameter),
        parameter_bias = mean_safe(group$best_parameter) - true_parameter,
        parameter_exact_recovery_rate = mean(
          abs(as_numeric_safe(group$best_parameter) - true_parameter) < 1e-10,
          na.rm = TRUE
        ),
        mean_nll = mean_safe(group$nll),
        mean_state_rmse = mean_safe(group$state_rmse),
        mean_selected_state_rmse = if (any(selected)) {
          mean_safe(group$state_rmse[selected])
        } else {
          NA_real_
        },
        mean_state_correlation = mean_safe(group$state_correlation),
        mean_unsigned_pe_correlation = mean_safe(group$unsigned_pe_correlation),
        mean_failure_rate = mean_safe(group$failure_rate)
      )
    }
  )
  write_summary(
    observation_summary, "state_model_observation_misspec_summary.csv"
  )
}

calibration_path <- file.path(
  output_dir, "calibration_robustness_replications.csv"
)
if (file.exists(calibration_path)) {
  calibration <- brm_read_csv(calibration_path)
  calibration_summary <- split_summarise(
    calibration,
    c(
      "calibration_rank", "adaptivity", "natural_deviation_sd",
      "outcome_intercept", "estimator"
    ),
    function(group) recovery_metrics(group, include_spearman = TRUE)
  )
  write_summary(
    calibration_summary, "calibration_robustness_summary.csv"
  )
}

message("R completed all Monte Carlo summaries.")
