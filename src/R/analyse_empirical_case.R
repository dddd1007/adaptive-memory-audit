project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source(file.path(project_dir, "src", "R", "theme_brm.R"))
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "outputs")
seed <- 20260731L
bootstrap_replications <- as.integer(Sys.getenv("BRM_BOOTSTRAP_REPS", unset = "1200"))

standardize <- function(x) {
  x <- as.numeric(x)
  deviation <- stats::sd(x)
  if (!is.finite(deviation) || deviation <= 0) return(rep(0, length(x)))
  (x - mean(x)) / deviation
}

design_matrix <- function(frame, numeric, categorical) {
  columns <- list(`(Intercept)` = rep(1, nrow(frame)))
  for (name in numeric) columns[[name]] <- standardize(frame[[name]])
  for (name in categorical) {
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

ridge_logistic <- function(x, y, case_weights = NULL, C = 10, max_iter = 80L,
                           tolerance = 1e-8) {
  if (is.null(case_weights)) case_weights <- rep(1, length(y))
  keep <- case_weights > 0
  x <- x[keep, , drop = FALSE]
  y <- as.numeric(y[keep])
  case_weights <- as.numeric(case_weights[keep])
  if (length(unique(y)) < 2L) return(rep(NA_real_, ncol(x)))
  lambda <- 1 / C
  beta <- rep(0, ncol(x))
  penalty <- diag(lambda, ncol(x))
  for (iteration in seq_len(max_iter)) {
    eta <- drop(x %*% beta)
    probability <- plogis(pmax(-30, pmin(30, eta)))
    working_weight <- pmax(1e-8, probability * (1 - probability)) * case_weights
    working_response <- eta + (y - probability) / pmax(1e-8, probability * (1 - probability))
    weighted_x <- x * sqrt(working_weight)
    information <- crossprod(weighted_x) + penalty
    score_target <- crossprod(weighted_x, working_response * sqrt(working_weight))
    candidate <- tryCatch(
      drop(solve(information, score_target)),
      error = function(condition) drop(qr.solve(information, score_target, tol = 1e-10))
    )
    if (max(abs(candidate - beta)) < tolerance) {
      beta <- candidate
      break
    }
    beta <- candidate
  }
  beta
}

participant_bootstrap <- function(frame, x, y, focal_index, C, seed_offset) {
  participants <- as.character(frame$participant_code)
  unique_participants <- sort(unique(participants))
  set.seed(seed + seed_offset)
  estimates <- rep(NA_real_, bootstrap_replications)
  for (replication in seq_len(bootstrap_replications)) {
    selected <- sample(unique_participants, length(unique_participants), replace = TRUE)
    multiplicity <- table(selected)
    weights <- as.numeric(multiplicity[match(participants, names(multiplicity))])
    weights[is.na(weights)] <- 0
    beta <- ridge_logistic(x, y, case_weights = weights, C = C)
    estimates[replication] <- beta[focal_index]
  }
  estimates[is.finite(estimates)]
}

estimate_specification <- function(label, frame, numeric, categorical, C, seed_offset) {
  x <- design_matrix(frame, numeric, categorical)
  y <- as.numeric(frame$next_failure)
  focal_index <- match("log_actual_gap", colnames(x))
  beta <- ridge_logistic(x, y, C = C)
  bootstrap <- participant_bootstrap(frame, x, y, focal_index, C, seed_offset)
  interval <- stats::quantile(bootstrap, c(0.025, 0.975), names = FALSE, type = 7)
  data.frame(
    model = label,
    n_transitions = nrow(frame),
    n_failures = sum(y),
    log_odds_coefficient = beta[focal_index],
    odds_ratio = exp(beta[focal_index]),
    ci_low = exp(interval[1L]),
    ci_high = exp(interval[2L]),
    n_bootstrap_successful = length(bootstrap),
    regularization_C = C,
    stringsAsFactors = FALSE
  )
}

wilson_interval <- function(successes, trials, confidence = 0.95) {
  if (trials == 0L) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - confidence) / 2)
  p <- successes / trials
  denominator <- 1 + z^2 / trials
  center <- (p + z^2 / (2 * trials)) / denominator
  half <- z * sqrt(p * (1 - p) / trials + z^2 / (4 * trials^2)) / denominator
  c(center - half, center + half)
}

transitions <- brm_read_csv(
  file.path(data_dir, "analysis_transitions_deidentified.csv")
)
transitions <- transitions[
  is.finite(transitions$next_gap_days) & transitions$next_gap_days > 0,
  , drop = FALSE
]
transitions$next_failure <- as.integer(transitions$next_first_ease == 1)
transitions$log_actual_gap <- log1p(transitions$next_gap_days)
transitions$log_planned_gap <- log1p(pmax(transitions$last_scheduled_interval_days, 0))
transitions$log_review_index <- log1p(transitions$daily_review_index)
transitions$current_failure <- as.integer(transitions$first_ease == 1)

models <- list(
  estimate_specification(
    "Naive gap-only", transitions,
    c("log_actual_gap"), character(), 100, 1L
  ),
  estimate_specification(
    "Observed-history adjusted", transitions,
    c("log_actual_gap", "log_review_index", "study_day", "current_failure"),
    c("last_ease", "participant_code"), 100, 2L
  ),
  estimate_specification(
    "Crossed learner-item adjusted", transitions,
    c("log_actual_gap", "log_review_index", "study_day", "current_failure"),
    c("last_ease", "participant_code", "card_code", "scheduler_encoding"), 10, 3L
  )
)

eligible <- transitions[transitions$last_scheduled_interval_days >= 0.5, , drop = FALSE]
models[[4L]] <- estimate_specification(
  "Schedule-conditioned deviation", eligible,
  c("log_actual_gap", "log_planned_gap", "log_review_index", "study_day", "current_failure"),
  c("last_ease", "participant_code", "card_code", "scheduler_encoding"), 10, 4L
)
model_table <- do.call(rbind, models)
write.csv(model_table, file.path(output_dir, "empirical_model_comparison.csv"), row.names = FALSE)

breaks <- c(0, 1.25, 2.25, 4.25, 7.25, Inf)
labels <- c("~1 day", "~2 days", "3–4 days", "5–7 days", "≥8 days")
transitions$gap_bin <- cut(transitions$next_gap_days, breaks = breaks, labels = labels,
                           include.lowest = FALSE, right = TRUE)
curve_rows <- lapply(labels, function(label) {
  group <- transitions[as.character(transitions$gap_bin) == label, , drop = FALSE]
  n <- nrow(group)
  failures <- sum(group$next_failure)
  interval <- wilson_interval(failures, n)
  data.frame(
    gap_bin = label,
    n_transitions = n,
    n_failures = failures,
    failure_rate = failures / n,
    ci_low = interval[1L],
    ci_high = interval[2L],
    stringsAsFactors = FALSE
  )
})
curve <- do.call(rbind, curve_rows)
write.csv(curve, file.path(output_dir, "empirical_gap_failure_curve.csv"), row.names = FALSE)

message("R empirical analysis complete: ", nrow(transitions), " transitions, ",
        sum(transitions$next_failure), " failures, ", bootstrap_replications,
        " participant bootstraps per model.")
