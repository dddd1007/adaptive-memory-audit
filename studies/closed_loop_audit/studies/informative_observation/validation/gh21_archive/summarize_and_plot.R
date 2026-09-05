# Independently recompute summary quantities from the Julia replication table.
args<-commandArgs(TRUE);root<-normalizePath(args[1]);out<-file.path(root,"outputs/confirmatory")
x<-read.csv(file.path(out,"replication_results.csv"));j<-read.csv(file.path(out,"summary.csv"))
checks<-list()
for(i in seq_len(nrow(j))) {
 z<-j[i,]; a<-subset(x,lambda==z$lambda & parameter==z$parameter);g<-a[a$valid,]
 stopifnot(nrow(a)==1000,!anyDuplicated(a$replication),all(is.finite(g$estimate)),all(is.finite(g$se)))
 values<-c(mean_estimate=mean(g$estimate),bias=mean(g$estimate-g$truth),bias_mcse=sd(g$estimate)/sqrt(nrow(g)),empirical_sd=sd(g$estimate),mean_se=mean(g$se),se_sd_ratio=mean(g$se)/sd(g$estimate),coverage_valid=mean(g$covered),coverage_unconditional=mean(a$covered),mean_observed_rate=mean(a$observed_rate))
 for(metric in names(values)) checks[[length(checks)+1]]<-data.frame(lambda=z$lambda,parameter=z$parameter,metric=metric,difference=abs(values[[metric]]-z[[metric]]))
}
delta<-read.csv(file.path(out,"paired_contrasts.csv"))
for(i in seq_len(nrow(delta))) {
 z<-delta[i,];a<-subset(x,lambda==0 & parameter==z$parameter);b<-subset(x,lambda==z$lambda & parameter==z$parameter)
 pair<-merge(a,b,by="replication",suffixes=c("_0","_1"));v<-pair$valid_0&pair$valid_1
 de<-pair$estimate_1[v]-pair$estimate_0[v];dc<-as.numeric(pair$covered_1)-as.numeric(pair$covered_0)
 values<-c(bias_difference=mean(de),bias_difference_mcse=sd(de)/sqrt(length(de)),coverage_difference_unconditional=mean(dc),coverage_difference_unconditional_mcse=sd(dc)/sqrt(nrow(pair)),coverage_difference_paired_valid=mean(dc[v]),coverage_difference_paired_valid_mcse=sd(dc[v])/sqrt(sum(v)))
 for(metric in names(values)) checks[[length(checks)+1]]<-data.frame(lambda=z$lambda,parameter=z$parameter,metric=metric,difference=abs(values[[metric]]-z[[metric]]))
}
check<-do.call(rbind,checks);write.csv(check,file.path(root,"validation/summary_independent_checks.csv"),row.names=FALSE)
stopifnot(all(check$difference<1e-10))
dir.create(file.path(root,"figures"),showWarnings=FALSE)
draw<-function() {
 par(mfrow=c(2,2),mar=c(4,4.8,2.3,.8),family="serif")
 labels<-c(delay_coefficient="Delay coefficient",state_coefficient="State coefficient",previous_failure_coefficient="Previous failure",retention="Retention")
 for(p in names(labels)) {
  a<-j[j$parameter==p,];a<-a[order(a$lambda),]
  ylim<-range(c(0,a$bias-1.96*a$bias_mcse,a$bias+1.96*a$bias_mcse))
  plot(a$lambda,a$bias,type="b",pch=19,col="#17517A",ylim=ylim,xlab=expression(lambda),ylab="Bias in raw parameter units",main=labels[[p]],xaxt="n")
  axis(1,at=0:2);abline(h=0,lty=2,col="gray50")
  arrows(a$lambda,a$bias-1.96*a$bias_mcse,a$lambda,a$bias+1.96*a$bias_mcse,angle=90,code=3,length=.035,col="#17517A")
 }
}
pdf(file.path(root,"figures/informative_observation_bias.pdf"),width=7,height=5.5);draw();dev.off()
svg(file.path(root,"figures/informative_observation_bias.svg"),width=7,height=5.5);draw();dev.off()
png(file.path(root,"figures/informative_observation_bias.png"),width=7,height=5.5,units="in",res=300,type="cairo");draw();dev.off()
cat("Independent summary checks:",nrow(check),"maximum difference:",max(check$difference),"\n")
