args <- commandArgs(TRUE)
root <- args[1]
suppressPackageStartupMessages(library(clubSandwich))
ref <- read.csv(file.path(root,"validation/julia_r_frozen_refit_comparison.csv"))
rows <- lapply(seq_len(nrow(ref)),function(i) {
  z <- ref[i,]
  stem <- if(z$design=="observational") "primary" else "mrt"
  d <- read.csv(file.path(root,"validation/frozen_data",sprintf("%s_%s_rep%d.csv",z$architecture,stem,z$replication)))
  d$plan_z <- as.numeric(scale(d$plan)); d$belief_z <- as.numeric(scale(d$belief))
  d$outcome <- if(z$architecture=="spaced") d$failure else 1-d$success
  if(z$design=="observational") {
    d$action_z <- as.numeric(scale(d$actual))
    d$state_z <- as.numeric(scale(if(z$architecture=="spaced") d$memory else d$theta))
    if(z$architecture=="tutor") d$recon_z <- as.numeric(scale(d$recon))
    formulas <- list(marginal=outcome~action_z, recorded_belief=outcome~action_z+belief_z+previous,
      policy_log=if(z$architecture=="spaced") outcome~action_z+plan_z+previous else outcome~action_z+plan_z+recon_z+previous,
      oracle_state=outcome~action_z+state_z+previous,reconstructed_history=outcome~action_z+recon_z+previous)
    fit <- glm(formulas[[z$estimator]],family=binomial(),data=d,control=glm.control(epsilon=1e-12,maxit=100))
    term <- "action_z"
    converged <- fit$converged
  } else {
    d$centered <- d$assign-.5
    fit <- lm(outcome~centered+plan_z+belief_z+previous,data=d)
    term <- "centered"
    converged <- TRUE
  }
  original_fit <- fit
  if(z$design=="observational") {
    X <- model.matrix(fit)
    p <- fitted(fit)
    Z <- sqrt(p*(1-p))*X
    working_y <- as.vector(Z %*% coef(fit)) + (d$outcome-p)/sqrt(p*(1-p))
    fit <- lm(working_y~Z-1)
    names(fit$coefficients) <- colnames(Z)
    colnames(fit$qr$qr) <- colnames(Z)
  }
  ct <- coef_test(fit,vcov="CR2",cluster=d$person,test="Satterthwaite")
  k <- match(term,rownames(ct))
  data.frame(z[c("architecture","design","replication","estimator")],converged,
    julia_estimate=z$julia_estimate,julia_se=z$julia_se,julia_df=z$julia_df,
    r_estimate=unname(coef(original_fit)[term]),r_se=ct$SE[k],r_df=ct$df_Satt[k])
})
out <- do.call(rbind,rows)
stopifnot(all(out$converged))
write.csv(out,args[2],row.names=FALSE)
cat("Independent clubSandwich CR2 on Pearson linearization;",nrow(out),"refits.\n")
