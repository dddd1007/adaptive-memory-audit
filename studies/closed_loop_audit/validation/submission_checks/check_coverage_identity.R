# Independent arithmetic audit of existing summaries; no new Monte Carlo.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
dest <- normalizePath(args[1])
study <- normalizePath(args[2])
s <- read.csv(file.path(study,"outputs/confirmatory/summary.csv"))
x <- read.csv(file.path(study,"outputs/confirmatory/replication_results.csv"))
stopifnot(nrow(s)==12L, nrow(x)==12000L)
x$valid <- as.logical(x$valid); x$covered <- as.logical(x$covered)
for (i in seq_len(nrow(s))) {
  z <- x[x$lambda==s$lambda[i] & x$parameter==s$parameter[i],]
  stopifnot(nrow(z)==1000L, length(unique(z$replication))==1000L,
            sum(z$valid)==s$n_valid[i],
            abs(mean(z$valid & z$covered)-s$coverage_unconditional[i])<1e-12,
            abs(mean(z$covered[z$valid])-s$coverage_valid[i])<1e-12)
}
s$failure_rate <- s$n_failed/s$n_total
s$coverage_gap_pp <- 100*(s$coverage_valid-s$coverage_unconditional)
stopifnot(max(abs(s$coverage_unconditional-(1-s$failure_rate)*s$coverage_valid))<1e-12)
write.csv(s[c("lambda","parameter","n_total","n_valid","n_failed","failure_rate","coverage_valid","coverage_unconditional","coverage_gap_pp")],file.path(dest,"coverage_denominator_audit.csv"),row.names=FALSE)
cat("All 12 summary rows agree with stored replication outcomes and the denominator identity.\n")
cat("Coverage gap range in percentage points:", range(s$coverage_gap_pp),"\n")
print(aggregate(coverage_gap_pp~lambda,s,function(v)c(min=min(v),max=max(v))))

# Deterministic check of the new analytical identity on a finite state mixture.
# This is not a DGP calibration, Monte Carlo experiment or an effect estimate.
states <- c(-1.3,0.2,1.7); mass <- c(.2,.5,.3)
beta <- .7; gamma <- -.8; step <- 1e-5
identity_errors <- c()
for (intercept in c(-4,-1,1)) for (action in c(-.7,.4,1.1)) {
  q <- function(e) sum(mass*plogis(intercept+beta*e+gamma*states))
  r <- plogis(intercept+beta*action+gamma*states)
  marginal <- q(action)
  exact <- beta*(1-sum(mass*(r-marginal)^2)/(marginal*(1-marginal)))
  numerical <- (qlogis(q(action+step))-qlogis(q(action-step)))/(2*step)
  identity_errors <- c(identity_errors,abs(exact-numerical))
  stopifnot(exact>0, exact<beta)
}
stopifnot(max(identity_errors)<1e-7)
cat("Nine deterministic derivative identity checks passed; max difference:",max(identity_errors),"\n")
