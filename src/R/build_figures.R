project_dir <- Sys.getenv("BRM_PROJECT_DIR", unset = normalizePath(getwd(), mustWork = TRUE))
source(file.path(project_dir, "src", "R", "theme_brm.R"))
output_dir <- file.path(project_dir, "outputs")
figure_dir <- file.path(project_dir, "figures")
if (!dir.exists(figure_dir)) dir.create(figure_dir, recursive = TRUE)
legacy_stems <- c("figure_3_robustness_and_experiment", "figure_4_empirical_case")
legacy_files <- unlist(lapply(legacy_stems, function(stem) {
  file.path(figure_dir, paste0(stem, c(".png", ".svg")))
}))
legacy_files <- legacy_files[file.exists(legacy_files)]
if (length(legacy_files)) unlink(legacy_files)

read_output <- function(name) brm_read_csv(file.path(output_dir, name))

draw_node <- function(x, y, width, height, label, fill, border, cex = 0.78) {
  rect(x - width / 2, y - height / 2, x + width / 2, y + height / 2,
       col = fill, border = border, lwd = 1.5)
  text(x, y, label, cex = cex, col = BRM_PALETTE$ink, font = 2)
}

draw_figure_1 <- function() {
  brm_set_theme(mar = c(0.6, 0.6, 0.6, 0.6), oma = c(0.4, 0.5, 2.8, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  rect(0.025, 0.555, 0.975, 0.885, col = "#F2F6F9",
       border = BRM_PALETTE$blue_light, lwd = 1.5)
  text(0.05, 0.84, "Psychological state and adaptive observation policy", adj = 0,
       font = 2, cex = 0.92, col = BRM_PALETTE$blue)

  top_x <- c(0.10, 0.30, 0.50, 0.70, 0.90)
  top_width <- c(0.15, 0.15, 0.17, 0.16, 0.14)
  top_labels <- c(
    "Retrieval\nhistory",
    "Latent memory\nstate",
    "Recall expectation\nand prediction error",
    "Scheduler belief\nabout memory",
    "Planned\ninterval"
  )
  top_fill <- c("white", "white", "white", "#F7F4EC", "#F7F4EC")
  top_border <- c(rep(BRM_PALETTE$blue, 3), rep(BRM_PALETTE$gold, 2))
  for (index in seq_along(top_x)) {
    draw_node(top_x[index], 0.69, top_width[index], 0.13, top_labels[index],
              top_fill[index], top_border[index], cex = 0.70)
    if (index < length(top_x)) {
      arrows(top_x[index] + top_width[index] / 2 + 0.008, 0.69,
             top_x[index + 1] - top_width[index + 1] / 2 - 0.008, 0.69,
             length = 0.065, angle = 22, col = BRM_PALETTE$ink, lwd = 1.4)
    }
  }

  draw_node(0.75, 0.40, 0.17, 0.12, "Enacted\ninterval", "#F7F4EC",
            BRM_PALETTE$gold, cex = 0.73)
  draw_node(0.51, 0.40, 0.17, 0.12, "Retrieval outcome\nsuccess / failure", "white",
            BRM_PALETTE$blue, cex = 0.70)
  draw_node(0.30, 0.40, 0.15, 0.12, "Updated memory\nstate", "white",
            BRM_PALETTE$blue, cex = 0.72)
  draw_node(0.88, 0.22, 0.20, 0.10, "Attention, fatigue,\nand availability", "white",
            BRM_PALETTE$grey, cex = 0.67)

  arrows(0.90, 0.62, 0.78, 0.47, length = 0.065, angle = 22,
         col = BRM_PALETTE$gold, lwd = 1.5)
  arrows(0.655, 0.40, 0.605, 0.40, length = 0.065, angle = 22,
         col = BRM_PALETTE$ink, lwd = 1.4)
  arrows(0.42, 0.40, 0.385, 0.40, length = 0.065, angle = 22,
         col = BRM_PALETTE$blue, lwd = 1.4)
  arrows(0.30, 0.47, 0.30, 0.615, length = 0.065, angle = 22,
         col = BRM_PALETTE$blue, lwd = 1.4)
  arrows(0.57, 0.46, 0.655, 0.62, length = 0.060, angle = 22,
         col = BRM_PALETTE$gold, lwd = 1.3)
  text(0.625, 0.535, "logged response", cex = 0.61,
       col = BRM_PALETTE$muted, srt = 39)

  arrows(0.84, 0.27, 0.79, 0.34, length = 0.055, angle = 22,
         col = BRM_PALETTE$grey, lwd = 1.2, lty = 2)
  arrows(0.80, 0.22, 0.58, 0.35, length = 0.055, angle = 22,
         col = BRM_PALETTE$grey, lwd = 1.2, lty = 2)
  text(0.72, 0.255, "unlogged common causes", cex = 0.61,
       col = BRM_PALETTE$muted)

  audit_x <- c(0.14, 0.38, 0.62, 0.86)
  audit_labels <- c("Representation", "Updating", "Policy", "Identification")
  audit_border <- c(BRM_PALETTE$blue, BRM_PALETTE$blue,
                    BRM_PALETTE$gold, BRM_PALETTE$ink)
  for (index in seq_along(audit_x)) {
    rect(audit_x[index] - 0.10, 0.075, audit_x[index] + 0.10, 0.13,
         col = "white", border = audit_border[index], lwd = 1.3)
    text(audit_x[index], 0.102, audit_labels[index], cex = 0.66,
         col = BRM_PALETTE$ink, font = 2)
  }
  text(0.04, 0.155, "Audit targets", adj = 0, cex = 0.65,
       col = BRM_PALETTE$muted)
  brm_outer_title(
    "Psychological feedback in adaptive memory",
    "Retrieval measures and changes memory; the scheduler selects the delay on which memory is next observed"
  )
}

draw_figure_2 <- function() {
  recovery <- read_output("parameter_recovery_summary.csv")
  labels <- c(naive = "Naive", history = "History adjusted",
              schedule_conditioned = "Schedule conditioned", oracle = "Oracle")
  colors <- c(naive = BRM_PALETTE$gold, history = BRM_PALETTE$plum,
              schedule_conditioned = BRM_PALETTE$teal, oracle = BRM_PALETTE$blue)
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1))
  brm_set_theme(mar = c(5.0, 4.4, 4.2, 0.8), oma = c(0.5, 0.3, 2.8, 0.3))

  baseline <- recovery[recovery$n_participants == 12 &
                         abs(recovery$adaptivity - 0.75) < 1e-8 &
                         abs(recovery$state_reliability - 0.80) < 1e-8, ]
  baseline <- baseline[match(names(labels), baseline$estimator), ]
  x <- seq_along(labels)
  ylim <- range(c(baseline$mean_estimate - 1.96 * baseline$monte_carlo_se_bias,
                  baseline$mean_estimate + 1.96 * baseline$monte_carlo_se_bias,
                  baseline$mean_truth), finite = TRUE)
  ylim <- ylim + c(-0.05, 0.05)
  plot(x, baseline$mean_estimate, type = "n", xaxt = "n", xlab = "Estimator",
       ylab = "Standardized log-odds coefficient", ylim = ylim)
  brm_grid_y()
  abline(h = mean(baseline$mean_truth), lty = 2, col = BRM_PALETTE$ink, lwd = 1.3)
  axis(1, at = x, labels = unname(labels), tick = FALSE, cex.axis = 0.67)
  brm_errorbar(x, baseline$mean_estimate - 1.96 * baseline$monte_carlo_se_bias,
               baseline$mean_estimate + 1.96 * baseline$monte_carlo_se_bias,
               color = unname(colors[names(labels)]), width = 0.07)
  points(x, baseline$mean_estimate, pch = 21, bg = unname(colors[names(labels)]),
         col = "white", cex = 1.25, lwd = 1.2)
  brm_panel_title("A", "Empirical-like recovery", "N = 12; reliability = .80; adaptivity = .75")

  signs <- recovery[recovery$n_participants == 30 &
                      abs(recovery$state_reliability - 0.80) < 1e-8, ]
  plot(range(signs$adaptivity), c(0, 1), type = "n", xlab = "Scheduler adaptivity",
       ylab = "Sign-error rate")
  brm_grid_y(seq(0, 1, 0.2))
  for (estimator in names(labels)) {
    part <- signs[signs$estimator == estimator, ]
    part <- part[order(part$adaptivity), ]
    lines(part$adaptivity, part$sign_error_rate, col = colors[estimator],
          lwd = 2, type = "b", pch = 16, cex = 0.7)
  }
  brm_legend("topleft", legend = unname(labels), col = unname(colors[names(labels)]),
             lwd = 2, pch = 16)
  brm_panel_title("B", "Direction errors under adaptivity", "200 replications per condition")

  coverage <- recovery[abs(recovery$adaptivity - 0.75) < 1e-8 &
                         abs(recovery$state_reliability - 0.80) < 1e-8, ]
  plot(range(coverage$n_participants), c(0, 1), type = "n", xlab = "Participants",
       ylab = "95% interval coverage", log = "x", xaxt = "n")
  axis(1, at = c(12, 30, 100), labels = c("12", "30", "100"))
  brm_grid_y(seq(0, 1, 0.2))
  abline(h = 0.95, lty = 2, col = BRM_PALETTE$ink)
  for (estimator in names(labels)) {
    part <- coverage[coverage$estimator == estimator, ]
    part <- part[order(part$n_participants), ]
    lines(part$n_participants, part$coverage_95, col = colors[estimator],
          lwd = 2, type = "b", pch = 16, cex = 0.7)
  }
  brm_panel_title("C", "Coverage across sample size", "Adaptivity = .75; nominal coverage = .95")
  brm_outer_title("Parameter recovery under adaptive scheduling",
                  "Monte Carlo estimates from the calibrated closed-loop data-generating process")
}

draw_figure_3 <- function() {
  model <- read_output("state_model_family_recovery_summary.csv")
  state <- read_output("state_model_state_recovery_summary.csv")
  layer <- read_output("state_model_layer_summary.csv")
  layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1.15))
  brm_set_theme(mar = c(5.1, 4.4, 4.2, 0.8), oma = c(0.5, 0.3, 2.8, 0.3))

  plot(c(2, 8), c(0, 1), type = "n", xlab = "Transitions per learner-item pair",
       ylab = "Generating-family recovery", xaxt = "n")
  axis(1, at = c(2, 8), labels = c("2", "8"))
  brm_grid_y(seq(0, 1, 0.2))
  family_colors <- c(RL = BRM_PALETTE$blue, Bayesian = BRM_PALETTE$gold)
  for (family in c("RL", "Bayesian")) {
    part <- model[model$learner_family == family, ]
    part <- part[order(part$n_stages), ]
    lines(part$n_stages, part$family_recovery_rate, type = "b", pch = 21,
          bg = family_colors[family], col = family_colors[family], lwd = 2, cex = 1)
  }
  brm_legend("bottomright", legend = c("RL", "Bayesian"),
             col = family_colors, pch = 16, lwd = 2)
  brm_panel_title("A", "Model-family recovery", "200 datasets per learner family and history length")

  state_long <- state[state$n_stages == 8, ]
  matched <- brm_match_bool(state_long$learner_family_match, TRUE)
  state_values <- c(mean(state_long$state_rmse[matched]),
                    mean(state_long$state_rmse[!matched]))
  policy_matched <- brm_match_bool(state_long$scheduler_family_match, TRUE)
  policy_values <- c(mean(state_long$policy_rmse[policy_matched]),
                     mean(state_long$policy_rmse[!policy_matched]))
  values <- rbind(state_values, policy_values)
  bar_positions <- barplot(values, beside = TRUE, names.arg = c("Matched", "Mismatched"),
                           col = c(BRM_PALETTE$blue, BRM_PALETTE$gold), border = NA,
                           ylim = c(0, max(values) * 1.25), ylab = "RMSE")
  brm_grid_y()
  barplot(values, beside = TRUE, add = TRUE,
          names.arg = rep("", ncol(values)),
          col = c(BRM_PALETTE$blue, BRM_PALETTE$gold), border = NA, axes = FALSE)
  brm_legend("topleft", legend = c("Learner state", "Scheduler plan"),
             fill = c(BRM_PALETTE$blue, BRM_PALETTE$gold))
  brm_panel_title("B", "Latent-state and policy recovery", "Eight transitions; lower values indicate better recovery")

  long_layer <- layer[layer$n_stages == 8, ]
  design_order <- c("Naive", "State proxy only", "Reconstructed policy", "Logged plan", "Oracle")
  bias <- sapply(design_order, function(design) {
    mean(abs(long_layer$bias[long_layer$design == design]))
  })
  bar_cols <- c(BRM_PALETTE$gold, BRM_PALETTE$plum, BRM_PALETTE$teal,
                BRM_PALETTE$blue, BRM_PALETTE$blue_light)
  barplot(bias, horiz = TRUE, names.arg = design_order, las = 1,
          col = bar_cols, border = NA, xlab = "Absolute cell-level mean bias",
          xlim = c(0, max(bias) * 1.20), cex.names = 0.70)
  brm_grid_x()
  barplot(bias, horiz = TRUE, names.arg = rep("", length(bias)), add = TRUE,
          col = bar_cols, border = NA, axes = FALSE)
  brm_panel_title("C", "Delay-coefficient recovery", "Averaged over eight-transition design cells")
  brm_outer_title("Cognitive-state and scheduling-policy recovery",
                  "Recovering memory dynamics and identifying a spacing effect are separate validation targets")
}

draw_figure_4 <- function() {
  context <- read_output("context_violation_summary.csv")
  mrt <- read_output("micro_randomized_wcls_summary.csv")
  layout(matrix(1:3, nrow = 1))
  brm_set_theme(mar = c(5.0, 4.4, 4.2, 0.8), oma = c(0.5, 0.3, 2.8, 0.3))
  context <- context[context$n_participants == 30 &
                       context$estimator %in% c("schedule_conditioned", "oracle"), ]
  plot(range(context$context_strength), range(c(context$bias, 0)), type = "n",
       xlab = "Unlogged-context strength", ylab = "Mean coefficient bias")
  brm_grid_y()
  abline(h = 0, lty = 2, col = BRM_PALETTE$ink)
  ccols <- c(schedule_conditioned = BRM_PALETTE$gold, oracle = BRM_PALETTE$blue)
  for (estimator in names(ccols)) {
    part <- context[context$estimator == estimator, ]
    part <- part[order(part$context_strength), ]
    lines(part$context_strength, part$bias, type = "b", pch = 16,
          col = ccols[estimator], lwd = 2)
  }
  brm_legend("topleft", legend = c("Schedule conditioned", "Oracle"),
             col = unname(ccols), pch = 16, lwd = 2)
  brm_panel_title("A", "Sensitivity to unlogged context", "N = 30; 160 replications per condition")

  power <- mrt[mrt$estimator == "covariate_adjusted_itt" &
                 abs(mrt$perturbation_compliance - 1) < 1e-8, ]
  plot(range(power$n_participants), c(0, 1), type = "n", xlab = "Participants",
       ylab = "Power for a positive ITT effect", log = "x", xaxt = "n")
  axis(1, at = c(12, 30, 60, 100), labels = c("12", "30", "60", "100"))
  brm_grid_y(seq(0, 1, 0.2))
  perturbation_colors <- c(`0.15` = BRM_PALETTE$blue_light,
                           `0.3` = BRM_PALETTE$teal,
                           `0.45` = BRM_PALETTE$blue)
  for (delta in c(0.15, 0.30, 0.45)) {
    part <- power[abs(power$perturbation_log_days - delta) < 1e-8, ]
    part <- part[order(part$n_participants), ]
    lines(part$n_participants, part$participant_t_power_positive, type = "b", pch = 16,
          col = perturbation_colors[as.character(delta)], lwd = 2)
  }
  abline(h = 0.80, lty = 2, col = BRM_PALETTE$ink)
  brm_legend("bottomright", legend = c("±0.15", "±0.30", "±0.45"),
             col = unname(perturbation_colors), pch = 16, lwd = 2)
  brm_panel_title("B", "Micro-randomized WCLS power", "Participant-t; full compliance; 250 replications per cell")

  compliance <- mrt[mrt$estimator == "covariate_adjusted_itt" &
                      mrt$n_participants == 100, ]
  compliance$label <- ifelse(compliance$perturbation_compliance > 0.9,
                             "100%", "70%")
  values <- tapply(compliance$participant_t_power_positive,
                   list(compliance$perturbation_log_days, compliance$label), mean)
  values <- values[order(as.numeric(rownames(values))), c("70%", "100%"), drop = FALSE]
  barplot(t(values), beside = TRUE, names.arg = paste0("±", rownames(values)),
          col = c(BRM_PALETTE$gold_light, BRM_PALETTE$gold), border = NA,
          ylim = c(0, 1), ylab = "Power", xlab = "Perturbation in log days")
  brm_grid_y(seq(0, 1, 0.2))
  barplot(t(values), beside = TRUE, add = TRUE,
          names.arg = rep("", nrow(values)),
          col = c(BRM_PALETTE$gold_light, BRM_PALETTE$gold), border = NA, axes = FALSE)
  brm_legend("topleft", legend = c("70% compliance", "100% compliance"),
             fill = c(BRM_PALETTE$gold_light, BRM_PALETTE$gold))
  brm_panel_title("C", "Compliance and WCLS power", "Participant-t; N = 100")
  brm_outer_title("Robustness analysis and experimental validation",
                  "Policy logging supports local comparisons until unrecorded context drives schedule deviations")
}

draw_figure_5 <- function() {
  empirical <- read_output("empirical_model_comparison.csv")
  empirical_sensitivity <- read_output("empirical_sensitivity_summary.csv")
  curve <- read_output("empirical_gap_failure_curve.csv")
  calibration <- read_output("empirical_simulation_calibration.csv")
  layout(matrix(1:3, nrow = 1), widths = c(1.05, 1.05, 1))
  brm_set_theme(mar = c(5.2, 4.4, 4.2, 0.8), oma = c(0.5, 0.3, 2.8, 0.3))

  x <- seq_len(nrow(curve))
  plot(x, curve$failure_rate, type = "n", xaxt = "n", ylim = c(0, max(curve$ci_high) * 1.15),
       xlab = "Enacted-gap bin", ylab = "Observed failure rate")
  brm_grid_y()
  axis(1, at = x, labels = curve$gap_bin, tick = FALSE, cex.axis = 0.67)
  brm_errorbar(x, curve$ci_low, curve$ci_high, BRM_PALETTE$blue, width = 0.08)
  points(x, curve$failure_rate, pch = 21, bg = BRM_PALETTE$blue,
         col = "white", cex = 1.2)
  brm_panel_title("A", "Failure rates by enacted delay", "2,045 transitions; Wilson 95% intervals")

  labels <- c("Naive gap-only", "Observed-history adjusted",
              "Crossed learner-item adjusted", "Schedule-conditioned deviation")
  empirical <- empirical[match(labels, empirical$model), ]
  focal_sensitivity <- empirical_sensitivity[
    abs(empirical_sensitivity$planned_interval_threshold_days - 0.5) < 1e-10 &
      empirical_sensitivity$outcome_definition == "Again" &
      empirical_sensitivity$regularization_C == 10,
    , drop = FALSE
  ]
  if (nrow(focal_sensitivity) != 1L || anyNA(empirical$model)) {
    stop("Expected one .5-day/C=10/Again sensitivity row and all four empirical models.")
  }
  # Preserve the three primary-model intervals above and use the declared
  # fixed-seed sensitivity bootstrap only for the schedule-conditioned row.
  schedule_row <- match("Schedule-conditioned deviation", empirical$model)
  empirical$ci_low[schedule_row] <- focal_sensitivity$ci_low
  empirical$ci_high[schedule_row] <- focal_sensitivity$ci_high
  y <- rev(seq_len(nrow(empirical)))
  xlim <- range(c(empirical$ci_low, empirical$ci_high, 1), finite = TRUE)
  plot(empirical$odds_ratio, y, type = "n", yaxt = "n", log = "x", xlim = xlim,
       xlab = "Odds ratio per SD enacted log delay", ylab = "")
  brm_grid_x()
  abline(v = 1, lty = 2, col = BRM_PALETTE$ink)
  axis(2, at = y, labels = c("Naive", "History", "Learner-item", "Schedule deviation"),
       tick = FALSE, las = 1, cex.axis = 0.68)
  segments(empirical$ci_low, y, empirical$ci_high, y, col = BRM_PALETTE$teal, lwd = 2)
  points(empirical$odds_ratio, y, pch = 21, bg = BRM_PALETTE$teal,
         col = "white", cex = 1.15)
  brm_panel_title("B", "Estimator-dependent associations", "Participant-bootstrap 95% intervals")

  calibration$display <- calibration$calibration_target
  calibration$observed_scaled <- calibration$observed_value
  calibration$simulated_scaled <- calibration$simulated_mean
  count_row <- calibration$scale == "count"
  calibration$observed_scaled[count_row] <- calibration$observed_value[count_row] / 2400
  calibration$simulated_scaled[count_row] <- calibration$simulated_mean[count_row] / 2400
  ratio <- calibration$observed_scaled / calibration$simulated_scaled
  ycal <- rev(seq_len(nrow(calibration)))
  plot(ratio, ycal, type = "n", yaxt = "n", xlim = c(0.75, 1.25),
       xlab = "Observed / simulated signature", ylab = "")
  brm_grid_x(seq(0.8, 1.2, 0.1))
  abline(v = 1, lty = 2, col = BRM_PALETTE$ink)
  axis(2, at = ycal, labels = c("Transitions", "Failure rate", "Schedule correlation", "Naive OR"),
       tick = FALSE, las = 1, cex.axis = 0.66)
  segments(1, ycal, ratio, ycal, col = BRM_PALETTE$gold, lwd = 2)
  points(ratio, ycal, pch = 21, bg = BRM_PALETTE$gold,
         col = "white", cex = 1.15)
  brm_panel_title("C", "Empirical-simulation calibration", "Values near 1 indicate close calibration")
  brm_outer_title("Empirical case study and calibrated signatures",
                  "Association direction varies across specifications while aggregate signatures remain calibrated")
}

draw_supplementary_figure_s1 <- function() {
  brm_set_theme(mar = c(0.6, 0.6, 0.6, 0.6), oma = c(0.4, 0.5, 2.8, 0.5))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  labels <- c(
    "5,656 de-identified\nreview events",
    "5,334 events in the\n14-day primary window",
    "3,245 unique\nparticipant-card-days",
    "2,045 cross-day\ntransitions",
    "Next-day retrieval:\n199 Again outcomes"
  )
  x <- c(.10, .30, .50, .70, .90)
  fills <- c("white", "white", "#F2F6F9", "#F2F6F9", "#F7F4EC")
  borders <- c(BRM_PALETTE$blue, BRM_PALETTE$blue, BRM_PALETTE$teal,
               BRM_PALETTE$teal, BRM_PALETTE$gold)
  for (index in seq_along(x)) {
    draw_node(x[index], .66, .165, .15, labels[index], fills[index],
              borders[index], cex = .68)
    if (index < length(x)) {
      arrows(x[index] + .09, .66, x[index + 1] - .09, .66,
             length = .06, angle = 22, col = BRM_PALETTE$ink, lwd = 1.4)
    }
  }
  checks <- c(
    "Sort by participant, card, local day, and time",
    "First rated event defines the delayed-retrieval outcome",
    "Final event carries scheduler state and planned interval forward",
    "Join to the first rated event on the next observed study day"
  )
  y <- c(.39, .30, .21, .12)
  for (index in seq_along(checks)) {
    points(.11, y[index], pch = 21, bg = BRM_PALETTE$teal,
           col = "white", cex = 1.05)
    text(.14, y[index], checks[index], adj = 0, cex = .72,
         col = BRM_PALETTE$ink)
  }
  rect(.66, .10, .95, .42, col = "white", border = BRM_PALETTE$blue_light,
       lwd = 1.4)
  text(.685, .375, "Independent lineage audit", adj = 0, font = 2,
       cex = .78, col = BRM_PALETTE$blue)
  text(.685, .30, "23 / 23 checks passed", adj = 0, font = 2,
       cex = .82, col = BRM_PALETTE$ink)
  text(.685, .235, "Maximum gap error", adj = 0, cex = .68,
       col = BRM_PALETTE$muted)
  text(.685, .18, expression(3.55 %*% 10^-15 ~ days), adj = 0,
       cex = .76, col = BRM_PALETTE$ink)
  brm_outer_title(
    "Exact reconstruction of prospective memory transitions",
    "The prior day's final scheduler event is linked to the next observed day's first rated retrieval"
  )
}

draw_supplementary_figure_s2 <- function() {
  search <- read_output("calibration_search_candidates.csv")
  robust <- read_output("calibration_robustness_summary.csv")
  search <- search[order(search$rank), ]
  top <- search[search$rank <= 9, ]
  layout(matrix(1:3, nrow = 1), widths = c(1.15, 1, 1))
  brm_set_theme(mar = c(5.2, 4.5, 4.4, .8), oma = c(.5, .3, 2.8, .3))

  discrepancy <- rbind(
    (top$mean_failure_rate - top$target_failure_rate) / .015,
    (top$mean_scheduled_actual_spearman - top$target_scheduled_actual_spearman) / .05,
    (top$mean_naive_log_odds - top$target_naive_log_odds) / .10
  )
  matplot(top$rank, t(discrepancy), type = "b", pch = c(16, 17, 15),
          lty = 1, lwd = 1.8,
          col = c(BRM_PALETTE$blue, BRM_PALETTE$teal, BRM_PALETTE$gold),
          xlab = "Calibration rank", ylab = "Standardized target discrepancy")
  brm_grid_y()
  abline(h = 0, lty = 2, col = BRM_PALETTE$ink)
  brm_legend("topleft", legend = c("Failure", "Plan-enactment r", "Naive log OR"),
             col = c(BRM_PALETTE$blue, BRM_PALETTE$teal, BRM_PALETTE$gold),
             pch = c(16, 17, 15), lwd = 1.8)
  brm_panel_title("A", "Candidate-target agreement", "Forty common-random-number datasets per candidate")

  naive <- robust[robust$estimator == "naive", ]
  naive <- naive[order(naive$calibration_rank), ]
  plot(naive$calibration_rank, naive$sign_error_rate, type = "b", pch = 16,
       col = BRM_PALETTE$gold, lwd = 2, ylim = c(0, 1),
       xlab = "Calibration rank", ylab = "Naive sign-error rate")
  brm_grid_y(seq(0, 1, .2))
  brm_panel_title("B", "Sign reversal is robust", "Independent 100-replication runs for the five best candidates")

  selected <- robust[robust$estimator %in% c("schedule_conditioned", "oracle"), ]
  plot(range(selected$calibration_rank), range(selected$absolute_bias), type = "n",
       xlab = "Calibration rank", ylab = "Mean absolute coefficient error")
  brm_grid_y()
  estimator_colors <- c(schedule_conditioned = BRM_PALETTE$teal,
                        oracle = BRM_PALETTE$blue)
  for (estimator in names(estimator_colors)) {
    part <- selected[selected$estimator == estimator, ]
    part <- part[order(part$calibration_rank), ]
    lines(part$calibration_rank, part$absolute_bias, type = "b", pch = 16,
          col = estimator_colors[estimator], lwd = 2)
  }
  brm_legend("topleft", legend = c("Plan conditioned", "Oracle"),
             col = unname(estimator_colors), pch = 16, lwd = 2)
  brm_panel_title("C", "Estimator robustness", "Independent 100-replication stream")
  brm_outer_title(
    "Calibration search and post-selection robustness",
    "The empirical-like mechanism lies in a broader region with the same qualitative conclusions"
  )
}

draw_supplementary_figure_s3 <- function() {
  observation <- read_output("state_model_observation_misspec_summary.csv")
  proxy <- read_output("context_proxy_summary.csv")
  empirical <- read_output("empirical_sensitivity_summary.csv")
  layout(matrix(1:3, nrow = 1), widths = c(1.12, 1, 1.12))
  brm_set_theme(mar = c(6.7, 4.6, 4.5, .8), oma = c(.5, .3, 2.8, .3))

  observation <- observation[
    as.character(observation$candidate_family) == as.character(observation$learner_family), ]
  conditions <- c("correct", "intercept_minus_0.25", "intercept_plus_0.25",
                  "state_slope_0.80", "combined_plus_0.25_slope_0.80")
  condition_labels <- c("Correct", "b0 -.25", "b0 +.25", "Slope .80", "Combined")
  x <- seq_along(conditions)
  plot(c(1, length(conditions)), c(0, 1), type = "n", xaxt = "n",
       xlab = "", ylab = "Family-recovery rate")
  axis(1, at = x, labels = condition_labels, tick = FALSE, las = 2, cex.axis = .62)
  brm_grid_y(seq(0, 1, .2))
  family_colors <- c(RL = BRM_PALETTE$blue, Bayesian = BRM_PALETTE$gold)
  for (family in names(family_colors)) {
    values <- vapply(conditions, function(condition) {
      mean(observation$family_recovery_rate[
        observation$learner_family == family &
          observation$observation_condition == condition], na.rm = TRUE)
    }, numeric(1L))
    lines(x, values, type = "b", pch = 16, col = family_colors[family], lwd = 2)
  }
  brm_legend("bottomleft", legend = names(family_colors),
             col = unname(family_colors), pch = 16, lwd = 2)
  brm_panel_title("A", "Observation-model sensitivity", "Eight transitions; averaged across scheduler families")

  proxy_estimators <- c(plan_conditioned = "Plan only", plan_plus_proxy = "Plan + proxy",
                        oracle = "Oracle")
  proxy_colors <- c(plan_conditioned = BRM_PALETTE$gold,
                    plan_plus_proxy = BRM_PALETTE$teal,
                    oracle = BRM_PALETTE$blue)
  plot(range(proxy$proxy_r2_target), range(c(proxy$bias, 0)), type = "n",
       xlab = expression("Context-proxy " * R^2), ylab = "Mean coefficient bias")
  brm_grid_y()
  abline(h = 0, lty = 2, col = BRM_PALETTE$ink)
  for (estimator in names(proxy_estimators)) {
    part <- proxy[proxy$estimator == estimator, ]
    part <- part[order(part$proxy_r2_target), ]
    lines(part$proxy_r2_target, part$bias, type = "b", pch = 16,
          col = proxy_colors[estimator], lwd = 2)
  }
  brm_legend("topright", legend = unname(proxy_estimators),
             col = unname(proxy_colors), pch = 16, lwd = 2)
  brm_panel_title("B", "Partial context measurement", "N = 30; strong context; 160 paired replications")

  thresholds <- sort(unique(empirical$planned_interval_threshold_days))
  plot(range(thresholds) + c(-.05, .05), range(empirical$odds_ratio), type = "n",
       xlab = "Planned-interval threshold (days)", ylab = "Observed odds ratio",
       log = "y", xaxt = "n")
  axis(1, at = thresholds, labels = format(thresholds, trim = TRUE))
  brm_grid_y()
  abline(h = 1, lty = 2, col = BRM_PALETTE$ink)
  outcome_colors <- c("Again" = BRM_PALETTE$blue,
                      "Again or Hard" = BRM_PALETTE$gold)
  penalty_symbols <- c(`1` = 16, `10` = 17, `100` = 15)
  for (outcome in names(outcome_colors)) {
    for (penalty in c(1, 10, 100)) {
      part <- empirical[empirical$outcome_definition == outcome &
                          empirical$regularization_C == penalty, ]
      part <- part[order(part$planned_interval_threshold_days), ]
      offset <- if (outcome == "Again") -.012 else .012
      lines(part$planned_interval_threshold_days + offset, part$odds_ratio,
            col = outcome_colors[outcome], lwd = 1, lty = 3)
      points(part$planned_interval_threshold_days + offset, part$odds_ratio,
             pch = penalty_symbols[as.character(penalty)],
             col = outcome_colors[outcome], cex = .82)
    }
  }
  brm_legend("topright", legend = c("Again", "Again-or-Hard", "C_ridge = 1", "C_ridge = 10", "C_ridge = 100"),
             col = c(BRM_PALETTE$blue, BRM_PALETTE$gold,
                     rep(BRM_PALETTE$ink, 3)),
             pch = c(16, 16, 16, 17, 15))
  brm_panel_title("C", "Empirical specification grid", "Point estimates; 24 threshold-outcome-penalty cells")
  brm_outer_title(
    "Robustness boundaries requested in review",
    "Measurement error constrains model interpretation; context proxies reduce bias; empirical associations remain descriptive"
  )
}

brm_export(file.path(figure_dir, "figure_1_sequential_learning_program"), draw_figure_1,
           width = 10.4, height = 6.3)
brm_export(file.path(figure_dir, "figure_2_parameter_recovery"), draw_figure_2,
           width = 11.2, height = 5.8)
brm_export(file.path(figure_dir, "figure_3_state_model_misspecification"), draw_figure_3,
           width = 11.3, height = 5.8)
brm_export(file.path(figure_dir, "figure_4_robustness_and_experiment"), draw_figure_4,
           width = 11.2, height = 5.8)
brm_export(file.path(figure_dir, "figure_5_empirical_case"), draw_figure_5,
           width = 11.4, height = 5.8)
brm_export(file.path(figure_dir, "supplementary_figure_s1_event_reconstruction"),
           draw_supplementary_figure_s1, width = 11.2, height = 5.8)
brm_export(file.path(figure_dir, "supplementary_figure_s2_calibration_robustness"),
           draw_supplementary_figure_s2, width = 11.4, height = 5.8)
brm_export(file.path(figure_dir, "supplementary_figure_s3_reviewer_robustness"),
           draw_supplementary_figure_s3, width = 11.6, height = 6.0)

message("R generated 5 main and 3 supplementary publication figures in ", figure_dir)
