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
# Independent partitioned integration: construct breakpoints in each initial-state
# normal coordinate, then vectorize response paths across all quadrature points.
gl <- function(q) {
 J<-matrix(0,q,q);k<-1:(q-1);J[cbind(k,k+1)]<-k/sqrt(4*k*k-1);J<-J+t(J)
 e<-eigen(J,symmetric=TRUE);o<-order(e$values)
 list(nodes=e$values[o],weights=2*e$vectors[1,o]^2)
}
partition_all <- function(phi,d,q) {
 n<-length(d$z0);out<-vector("list",n)
 for(i in seq_len(n)) {
  a<-sqrt(c[["estimation.posterior_variance"]]);b<-c[["estimation.posterior_mean_multiplier"]]*d$z0[i]
  left<- -10;right<-10;cuts<-c(-10,-6,-3,0,3,6,10)
  for(t in seq_len(ncol(d$a)-1)) {
   inc<-if(d$r[i,t]==0) 0 else if(d$y[i,t]==0) c[["learner.success_increment"]] else c[["learner.failure_increment"]]
   a<-phi*a;b<-phi*b+inc
   if(a>.Machine$double.eps && left<right) {
    low<-(-4.5-b)/a;high<-(4.5-b)/a
    if(left<low && low<right) cuts<-c(cuts,low)
    if(left<high && high<right) cuts<-c(cuts,high)
    left<-max(left,low);right<-min(right,high)
   }
  }
  cuts<-sort(unique(cuts));half<-diff(cuts)/2;mid<-(head(cuts,-1)+tail(cuts,-1))/2
  z<-as.vector(outer(q$nodes,half)+rep(mid,each=length(q$nodes)))
  w<-rep(q$weights,length(half))*rep(half,each=length(q$nodes))*dnorm(z)
  out[[i]]<-list(m=c[["estimation.posterior_mean_multiplier"]]*d$z0[i]+sqrt(c[["estimation.posterior_variance"]])*z,w=w,id=rep(i,length(z)))
 }
 list(m=unlist(lapply(out,`[[`,"m")),w=unlist(lapply(out,`[[`,"w")),id=unlist(lapply(out,`[[`,"id")))
}
make_fn <- function(d,q,aware) {
 cache_phi<-NA_real_;cache<-NULL
 function(theta) {
  phi<-plogis(theta[4])
  if(is.na(cache_phi) || phi!=cache_phi) {cache<<-partition_all(phi,d,q);cache_phi<<-phi}
  id<-cache$id;m<-cache$m;ll<-log(cache$w);previous<-rep(0,length(d$z0))
  for(t in seq_len(ncol(d$a))) {
   rr<-d$r[id,t];yy<-ifelse(rr==1,d$y[id,t],0)
   if(aware) {eta_r<-d$alpha+d$lambda*m;ll<-ll+rr*eta_r-softplus(eta_r)}
   eta<-c[["response.intercept"]]+theta[1]*d$a[id,t]+theta[2]*m+theta[3]*previous[id]
   ll<-ll+rr*(yy*eta-softplus(eta))
   inc<-ifelse(rr==1,ifelse(yy==0,c[["learner.success_increment"]],c[["learner.failure_increment"]]),0)
   m<-pmin(pmax(phi*m+inc,-4.5),4.5)
   previous<-ifelse(d$r[,t]==1,d$y[,t],previous)
  }
  # Likelihoods for eight opportunities do not approach underflow in validation.
  -sum(log(rowsum(exp(ll),id,reorder=FALSE)[,1]))
 }
}
starts <- list(c(.1,-.4,.2,log(.65/.35)),c(.5,-.9,.6,log(.85/.15)),c(-.1,-.2,.5,log(.4/.6)))
rows <- list()
for(i in seq_len(nrow(manifest))) {
  z <- manifest[i,]; d <- read_data(file.path(study,z$dataset)); q <- gl(z$nodes)
  fn <- make_fn(d,q,as.logical(z$aware))
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
