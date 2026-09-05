#!/usr/bin/env Rscript

# Independent summarization, figure production, and frozen-data validation.
# The data-generating processes and all Monte Carlo estimates originate in Julia.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]]) else normalizePath(getwd())
out <- file.path(root, "outputs_julia")
val_dir <- file.path(root, "validation")
dir.create(val_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(root, "figures_julia")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages(library(clubSandwich))

read_csv <- function(name) read.csv(file.path(out, name), stringsAsFactors = FALSE)

obs <- rbind(read_csv("observational_replications_spaced.csv"),
             read_csv("observational_replications_tutor.csv"))
mrt <- rbind(read_csv("mrt_effect_replications_spaced.csv"),
             read_csv("mrt_effect_replications_tutor.csv"))
mrt_null <- rbind(read_csv("mrt_null_replications_spaced.csv"),
                  read_csv("mrt_null_replications_tutor.csv"))
mrt_comp <- rbind(read_csv("mrt_compliance_replications_spaced.csv"),
                  read_csv("mrt_compliance_replications_tutor.csv"))

write.csv(obs, file.path(out, "observational_replications.csv"), row.names = FALSE)
write.csv(mrt, file.path(out, "mrt_effect_replications.csv"), row.names = FALSE)
write.csv(mrt_null, file.path(out, "mrt_null_replications.csv"), row.names = FALSE)
write.csv(mrt_comp, file.path(out, "mrt_compliance_replications.csv"), row.names = FALSE)

mc_summary <- function(df, groups) {
  key <- interaction(df[groups], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, key)
  rows <- lapply(pieces, function(x) {
    n <- nrow(x)
    est_sd <- sd(x$estimate)
    result <- c(as.list(x[1, groups, drop = FALSE]), list(
      n_replications = n,
      mean_truth = if ("truth" %in% names(x)) mean(x$truth) else NA_real_,
      mean_estimate = mean(x$estimate),
      bias = if ("truth" %in% names(x)) mean(x$estimate - x$truth) else NA_real_,
      empirical_sd = est_sd,
      mcse_mean = est_sd / sqrt(n),
      mean_model_se = if ("se" %in% names(x)) mean(x$se) else NA_real_,
      model_se_ratio = if ("se" %in% names(x)) mean(x$se) / est_sd else NA_real_,
      mean_observed_rate = if ("observed_rate" %in% names(x)) mean(x$observed_rate) else NA_real_,
      coverage = if ("covered" %in% names(x)) mean(as.logical(x$covered)) else NA_real_,
      coverage_mcse = if ("covered" %in% names(x)) {
        p <- mean(as.logical(x$covered)); sqrt(p * (1 - p) / n)
      } else NA_real_,
      sign_error = if ("sign_error" %in% names(x)) mean(as.logical(x$sign_error)) else NA_real_,
      rejection = if ("rejected" %in% names(x)) mean(as.logical(x$rejected)) else NA_real_,
      power = if ("detected" %in% names(x)) mean(as.logical(x$detected)) else NA_real_,
      mean_df = if ("df" %in% names(x)) mean(x$df) else NA_real_
    ))
    as.data.frame(result, check.names = FALSE)
  })
  ans <- do.call(rbind, rows)
  rownames(ans) <- NULL
  ans
}

obs_sum <- mc_summary(obs, c("architecture", "adaptivity", "estimator"))
mrt_sum <- mc_summary(mrt, c("architecture", "n_people"))
null_sum <- mc_summary(mrt_null, c("architecture", "n_people"))
comp_sum <- mc_summary(mrt_comp, c("architecture", "n_people", "compliance"))
context <- read_csv("context_sensitivity_replications.csv")
context_grouped <- context
context_grouped$proxy_r2[is.na(context_grouped$proxy_r2)] <- -1
context_sum <- mc_summary(context_grouped, c("context_strength", "proxy_r2"))
context_sum$proxy_r2[context_sum$proxy_r2 == -1] <- NA_real_
oracle <- read_csv("oracle_policy_replications.csv")
oracle_sum <- mc_summary(oracle, c("estimator"))
dynamic <- read_csv("dynamic_parameter_recovery_replications.csv")
dynamic_sum <- mc_summary(dynamic, c("architecture", "mechanism", "adaptivity",
                                    "context_strength", "observation_strength", "parameter"))
selection <- read_csv("observation_selection_stress_replications.csv")
selection_sum <- mc_summary(selection, c("architecture", "mechanism", "adaptivity",
                                        "context_strength", "observation_strength", "parameter"))
update <- read_csv("update_class_replications.csv")
update_key <- interaction(update[c("adaptivity", "true_family", "observation_model")], drop = TRUE)
update_sum <- do.call(rbind, lapply(split(update, update_key), function(x) {
  p <- mean(as.logical(x$recovered)); n <- nrow(x)
  data.frame(adaptivity=x$adaptivity[1], true_family=x$true_family[1],
             observation_model=x$observation_model[1], n_replications=n,
             recovery_rate=p, rate_mcse=sqrt(p*(1-p)/n))
}))
rownames(update_sum) <- NULL

write.csv(obs_sum, file.path(out, "observational_summary.csv"), row.names = FALSE)
write.csv(mrt_sum, file.path(out, "mrt_effect_summary.csv"), row.names = FALSE)
write.csv(null_sum, file.path(out, "mrt_null_summary.csv"), row.names = FALSE)
write.csv(comp_sum, file.path(out, "mrt_compliance_summary.csv"), row.names = FALSE)
write.csv(context_sum, file.path(out, "context_sensitivity_summary.csv"), row.names = FALSE)
write.csv(oracle_sum, file.path(out, "oracle_policy_summary.csv"), row.names = FALSE)
write.csv(dynamic_sum, file.path(out, "dynamic_parameter_recovery_summary.csv"), row.names = FALSE)
write.csv(selection_sum, file.path(out, "observation_selection_stress_summary.csv"), row.names = FALSE)
write.csv(update_sum, file.path(out, "update_class_summary.csv"), row.names = FALSE)

# Figures are produced from the Julia replication outputs.
draw_prediction <- function(device, filename) {
  device(file.path(fig_dir, filename), width=10.5, height=4.4)
  par(mfrow=c(1,2), mar=c(4.4,5.8,2.5,1.0), oma=c(0,0,0,0))
  for (arch in c("spaced","tutor")) {
    d <- obs[obs$architecture==arch & obs$estimator=="marginal",]
    aa <- sort(unique(d$adaptivity))
    est <- sapply(aa, function(a) mean(d$estimate[d$adaptivity==a]))
    pred <- sapply(aa, function(a) mean(d$prediction[d$adaptivity==a]))
    truth <- sapply(aa, function(a) mean(d$truth[d$adaptivity==a]))
    ylim <- range(c(est,pred,truth,0))
    plot(aa, truth, type="l", lty=2, lwd=2, col="#222222", ylim=ylim,
         xlab="Policy adaptivity", ylab="Standardized adverse-response coefficient",
         main=if (arch=="spaced") "Spaced retrieval" else "Adaptive tutor")
    abline(h=0,col="#999999",lwd=.8)
    lines(aa,pred,type="o",pch=16,lwd=2,col="#0072B2")
    lines(aa,est,type="o",pch=15,lwd=2,col="#D55E00")
    legend("bottomleft",c("Conditional truth","OVB prediction","Monte Carlo mean"),
           col=c("#222222","#0072B2","#D55E00"),lty=c(2,1,1),pch=c(NA,16,15),bty="n",cex=.82)
  }
  dev.off()
}
draw_prediction(function(...) png(..., units="in", res=400), "figure_analytic_prediction.png")
draw_prediction(pdf, "figure_analytic_prediction.pdf")

draw_recovery <- function(device, filename) {
  device(file.path(fig_dir, filename), width=10.5, height=4.4)
  par(mfrow=c(1,2), mar=c(4.4,5.8,2.5,1.0))
  palette <- c(marginal="#D55E00", recorded_belief="#009E73",
               reconstructed_history="#CC79A7", policy_log="#0072B2",
               oracle_state="#222222")
  labels <- c(marginal="Marginal", recorded_belief="Recorded belief",
              reconstructed_history="Analyst history", policy_log="Logged plan",
              oracle_state="Latent-state oracle")
  for (arch in c("spaced","tutor")) {
    s <- obs_sum[obs_sum$architecture==arch,]
    keep <- if (arch=="spaced") c("marginal","recorded_belief","policy_log","oracle_state") else
      c("marginal","recorded_belief","reconstructed_history","policy_log","oracle_state")
    ylim <- range(c(s$mean_truth,s$mean_estimate,0))
    plot(NA,xlim=range(s$adaptivity),ylim=ylim,xlab="Policy adaptivity",
         ylab="Standardized adverse-response coefficient",
         main=if(arch=="spaced") "Spaced retrieval" else "Adaptive tutor")
    abline(h=0,col="#999999",lwd=.8)
    truth <- aggregate(mean_truth~adaptivity,s,mean)
    lines(truth$adaptivity,truth$mean_truth,lty=2,lwd=2,col="#666666")
    for (nm in keep) {
      z <- s[s$estimator==nm,]; z <- z[order(z$adaptivity),]
      lines(z$adaptivity,z$mean_estimate,type="o",pch=16,lwd=1.8,col=palette[[nm]])
    }
    legend("topleft",c("Truth",labels[keep]),col=c("#666666",palette[keep]),
           lty=c(2,rep(1,length(keep))),pch=c(NA,rep(16,length(keep))),bty="n",cex=.72)
  }
  dev.off()
}
draw_recovery(function(...) png(...,units="in",res=400),"figure_estimator_recovery.png")
draw_recovery(pdf,"figure_estimator_recovery.pdf")

draw_boundary <- function(device, filename) {
  device(file.path(fig_dir,filename),width=10.5,height=4.4)
  par(mfrow=c(1,2),mar=c(4.4,5.8,2.5,1.0))
  cp <- context_sum[context_sum$context_strength==1 & !is.na(context_sum$proxy_r2),]
  plot(cp$proxy_r2,cp$bias,type="o",pch=16,lwd=2,col="#D55E00",
       xlab=expression("Context-proxy "*R^2),ylab="Bias of plan-adjusted coefficient",
       main="Unrecorded context")
  abline(h=0,col="#999999",lwd=.8)
  plot(NA,xlim=range(mrt_sum$n_people),ylim=c(.90,.98),xlab="Learner clusters",
       ylab="95% interval coverage",main="Randomized benchmark")
  abline(h=.95,lty=2,col="#666666")
  for (arch in c("spaced","tutor")) {
    z<-mrt_sum[mrt_sum$architecture==arch,]; z<-z[order(z$n_people),]
    lines(z$n_people,z$coverage,type="o",pch=if(arch=="spaced")16 else 15,lwd=2,
          col=if(arch=="spaced")"#0072B2" else "#009E73")
  }
  legend("bottomright",c("Spaced retrieval","Adaptive tutor"),
         col=c("#0072B2","#009E73"),pch=c(16,15),lty=1,bty="n")
  dev.off()
}
draw_boundary(function(...) png(...,units="in",res=400),"figure_assumption_boundary.png")
draw_boundary(pdf,"figure_assumption_boundary.pdf")

# Independent refits of frozen Julia-generated data.
coef_row <- function(fit, cluster, type, term, test="Satterthwaite") {
  ct <- coef_test(fit, vcov=type, cluster=cluster, test=test)
  rn <- rownames(ct)
  i <- match(term, rn)
  if (is.na(i)) i <- match(term, names(coef(fit)))
  data.frame(r_estimate=unname(coef(fit)[term]), r_se=ct$SE[i], r_df=ct$df_Satt[i])
}

validation_rows <- list()
v <- 1
for (arch in c("spaced","tutor")) {
  for (rep in 0:4) {
    d <- read.csv(file.path(val_dir,"frozen_data",sprintf("%s_primary_rep%d.csv",arch,rep)))
    d$action_z <- as.numeric(scale(d$actual)); d$plan_z <- as.numeric(scale(d$plan))
    d$belief_z <- as.numeric(scale(d$belief))
    if (arch=="spaced") {
      d$state_z <- as.numeric(scale(d$memory)); d$outcome <- d$failure
      formulas <- list(marginal=outcome~action_z,
                       recorded_belief=outcome~action_z+belief_z+previous,
                       policy_log=outcome~action_z+plan_z+previous,
                       oracle_state=outcome~action_z+state_z+previous)
    } else {
      d$recon_z <- as.numeric(scale(d$recon)); d$state_z <- as.numeric(scale(d$theta)); d$outcome <- 1-d$success
      formulas <- list(marginal=outcome~action_z,
                       recorded_belief=outcome~action_z+belief_z+previous,
                       reconstructed_history=outcome~action_z+recon_z+previous,
                       policy_log=outcome~action_z+plan_z+recon_z+previous,
                       oracle_state=outcome~action_z+state_z+previous)
    }
    for (nm in names(formulas)) {
      fit <- glm(formulas[[nm]], family=binomial(), data=d)
      rr <- coef_row(fit,d$person,"CR2","action_z")
      jr <- obs[obs$architecture==arch & obs$adaptivity==.70 & obs$replication==rep & obs$estimator==nm,]
      validation_rows[[v]] <- data.frame(architecture=arch,design="observational",replication=rep,
        estimator=nm,julia_estimate=jr$estimate,julia_se=jr$se,julia_df=jr$df,
        r_estimate=rr$r_estimate,r_se=rr$r_se,r_df=rr$r_df)
      v <- v+1
    }
  }
}

for (arch in c("spaced","tutor")) {
  for (rep in 0:2) {
    d <- read.csv(file.path(val_dir,"frozen_data",sprintf("%s_mrt_rep%d.csv",arch,rep)))
    d$centered <- d$assign-.5; d$plan_z <- as.numeric(scale(d$plan)); d$belief_z <- as.numeric(scale(d$belief))
    d$outcome <- if (arch=="spaced") d$failure else 1-d$success
    fit <- lm(outcome~centered+plan_z+belief_z+previous,data=d)
    rr <- coef_row(fit,d$person,"CR2","centered")
    jr <- mrt[mrt$architecture==arch & mrt$n_people==30 & mrt$replication==rep,]
    validation_rows[[v]] <- data.frame(architecture=arch,design="mrt",replication=rep,
      estimator="WCLS",julia_estimate=jr$estimate,julia_se=jr$se,julia_df=jr$df,
      r_estimate=rr$r_estimate,r_se=rr$r_se,r_df=rr$r_df)
    v <- v+1
  }
}

validation <- do.call(rbind,validation_rows)
validation$estimate_difference <- validation$julia_estimate-validation$r_estimate
validation$se_difference <- validation$julia_se-validation$r_se
validation$df_difference <- validation$julia_df-validation$r_df
write.csv(validation,file.path(val_dir,"julia_r_frozen_refit_comparison.csv"),row.names=FALSE)

checks <- data.frame(
  metric=c("maximum absolute coefficient difference","maximum absolute standard-error difference",
           "maximum absolute degrees-of-freedom difference",
           "Julia version","R version","clubSandwich version"),
  value=c(format(max(abs(validation$estimate_difference)),digits=12),
          format(max(abs(validation$se_difference)),digits=12),
          format(max(abs(validation$df_difference)),digits=12),
          sub("^julia_version=", "", as.character(readLines(file.path(out,"run_manifest.txt"),warn=FALSE)[2])),
          R.version.string,as.character(packageVersion("clubSandwich")))
)
write.csv(checks,file.path(val_dir,"validation_manifest.csv"),row.names=FALSE)
cat("R summaries and validation completed\n")
