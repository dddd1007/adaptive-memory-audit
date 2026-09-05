# Independent base-R implementation. Does not call Julia or reuse Julia derivatives.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
study <- normalizePath(args[1]); manifest <- read.csv(args[2]); output <- args[3]
numeric_toml <- function(path) {
  ans <- list(); section <- ""
  for(line in readLines(path)) {
    line <- trimws(sub("#.*$","",line))
    if(grepl("^\\[.*\\]$",line)) section <- sub("^\\[(.*)\\]$","\\1",line)
    if(grepl("^[A-Za-z_]+ *= *[-+0-9.eE]+$",line)) {
      part <- strsplit(line,"=",fixed=TRUE)[[1]]
      key <- paste(c(if(nzchar(section)) section,trimws(part[1])),collapse=".")
      ans[[key]] <- as.numeric(trimws(part[2]))
    }
  }
  ans
}
c <- numeric_toml(file.path(study,"config/design.toml"))
s <- numeric_toml(file.path(study,"config/execution.toml"))
gh <- function(q) {
  J <- matrix(0,q,q); J[cbind(1:(q-1),2:q)] <- sqrt(1:(q-1)); J <- J+t(J)
  e <- eigen(J,symmetric=TRUE); o <- order(e$values)
  list(nodes=e$values[o],weights=e$vectors[1,o]^2)
}
stopifnot(abs(sum(gh(21)$weights)-1)<1e-12,abs(sum(gh(21)$weights*gh(21)$nodes^2)-1)<1e-12)
read_data <- function(path) {
  d <- read.csv(path); n <- max(d$sequence); tt <- max(d$opportunity)
  mat <- function(col) { a<-matrix(NA_real_,n,tt); a[cbind(d$sequence,d$opportunity)]<-d[[col]];a }
  stopifnot(nrow(d)==n*tt,!anyDuplicated(d[c("sequence","opportunity")]))
  list(z0=d$z0[match(1:n,d$sequence)],a=mat("action"),r=mat("observed"),y=mat("y"),lambda=d$lambda[1],alpha=d$alpha[1])
}
softplus <- function(x) pmax(x,0)+log1p(exp(-abs(x)))
nll <- function(theta,d,q,aware=FALSE) {
  n <- length(d$z0); nq <- length(q$nodes)
  m <- outer(c[["estimation.posterior_mean_multiplier"]]*d$z0,rep(1,nq)) +
       outer(rep(1,n),sqrt(c[["estimation.posterior_variance"]])*q$nodes)
  ll <- matrix(rep(log(q$weights),each=n),n,nq); previous <- rep(0,n)
  phi <- plogis(theta[4])
  for(t in seq_len(ncol(d$a))) {
    rr <- d$r[,t]; yy <- ifelse(rr==1,d$y[,t],0)
    if(aware) {
      eta_r <- d$alpha+d$lambda*m
      ll <- ll+rr*eta_r-softplus(eta_r)
    }
    eta <- c[["response.intercept"]]+theta[1]*d$a[,t]+theta[2]*m+theta[3]*previous
    ll <- ll+rr*(yy*eta-softplus(eta))
    increment <- ifelse(rr==1,ifelse(yy==0,c[["learner.success_increment"]],c[["learner.failure_increment"]]),0)
    m <- pmin(pmax(phi*m+increment,-4.5),4.5)
    previous <- ifelse(rr==1,yy,previous)
  }
  mx <- apply(ll,1,max)
  -sum(mx+log(rowSums(exp(ll-mx))))
}
starts <- list(c(.1,-.4,.2,log(.65/.35)),c(.5,-.9,.6,log(.85/.15)),c(-.1,-.2,.5,log(.4/.6)))
rows <- list()
for(i in seq_len(nrow(manifest))) {
  z <- manifest[i,]; d <- read_data(file.path(study,z$dataset)); q <- gh(z$nodes)
  fn <- function(theta) nll(theta,d,q,as.logical(z$aware))
  theta_j <- as.numeric(z[paste0("theta",1:4)])
  objective_difference <- abs(fn(theta_j)-z$nll)
  fits <- lapply(starts,function(start) optim(start,fn,method="BFGS",control=list(reltol=1e-12,maxit=500,ndeps=rep(1e-5,4))))
  best <- fits[[which.min(vapply(fits,`[[`,numeric(1),"value"))]]
  theta <- best$par; H <- optimHess(theta,fn,control=list(ndeps=rep(1e-4,4)))
  eig <- eigen(H,symmetric=TRUE,only.values=TRUE)$values
  gradient <- sapply(1:4,function(j) { delta<-rep(0,4);delta[j]<-1e-5;(fn(theta+delta)-fn(theta-delta))/2e-5 })
  valid <- best$convergence==0 && min(eig)>0 && kappa(H,exact=TRUE)<1e10 && max(abs(gradient))<1e-3
  se <- rep(NA_real_,4); est <- c(theta[1:3],plogis(theta[4]))
  if(valid) { se<-sqrt(diag(solve(H))); se[4]<-se[4]*est[4]*(1-est[4]) }
  for(j in 1:4) {
    ed <- abs(est[j]-z[[paste0("estimate",j)]]); sd <- abs(se[j]-z[[paste0("se",j)]])
    pass <- valid && objective_difference<s$r_objective_atol &&
      ed<=s$r_parameter_atol+s$r_parameter_rtol*abs(z[[paste0("estimate",j)]]) &&
      sd<=s$r_se_atol+s$r_se_rtol*abs(z[[paste0("se",j)]])
    rows[[length(rows)+1]] <- data.frame(dataset=z$dataset,lambda=z$lambda,aware=z$aware,nodes=z$nodes,parameter=j,
      julia_estimate=z[[paste0("estimate",j)]],r_estimate=est[j],julia_se=z[[paste0("se",j)]],r_se=se[j],
      objective_difference=objective_difference,estimate_difference=ed,se_difference=sd,r_valid=valid,
      r_gradient_max=max(abs(gradient)),r_hessian_min=min(eig),pass=pass)
  }
  cat("Validated",z$dataset,"aware",z$aware,"objective difference",objective_difference,"valid",valid,"\n");flush.console()
}
out <- do.call(rbind,rows); write.csv(out,output,row.names=FALSE)
writeLines(capture.output(sessionInfo()),paste0(output,".session.txt"))
cat("Checks",nrow(out),"passed",sum(out$pass),"\n")
if(!all(out$pass)) quit(status=1)
