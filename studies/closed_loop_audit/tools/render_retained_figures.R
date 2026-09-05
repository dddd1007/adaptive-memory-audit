args<-commandArgs(TRUE);root<-normalizePath(args[1]);out<-file.path(root,"outputs_julia");dest<-file.path(root,"figures")
obs<-rbind(read.csv(file.path(out,"observational_replications_spaced.csv")),read.csv(file.path(out,"observational_replications_tutor.csv")))
s<-read.csv(file.path(out,"observational_summary.csv"));cx<-read.csv(file.path(out,"context_sensitivity_summary.csv"));mrt<-read.csv(file.path(out,"mrt_effect_summary.csv"))
colours<-c(marginal="#B44514",recorded_belief="#1B7959",reconstructed_history="#A04B95",policy_log="#006A9A",oracle_state="#333333")
labels<-c(marginal="Marginal",recorded_belief="Recorded belief",reconstructed_history="Analyst history",policy_log="Logged plan",oracle_state="Latent-state oracle")
setup<-function(){layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),heights=c(4,1));par(family="serif",mar=c(4,4.4,2,.8),cex=.95)}
draw_prediction<-function(){
 setup()
 for(arch in c("spaced","tutor")) {
  d<-subset(obs,architecture==arch & estimator=="marginal");a<-sort(unique(d$adaptivity));est<-sapply(a,function(v) mean(d$estimate[d$adaptivity==v]));pred<-sapply(a,function(v) mean(d$prediction[d$adaptivity==v]));truth<-sapply(a,function(v) mean(d$truth[d$adaptivity==v]))
  plot(a,truth,type="l",lty=2,lwd=1.5,ylim=range(c(est,pred,truth,0)),xlab="Policy adaptivity",ylab="Standardized coefficient",main=if(arch=="spaced")"Spaced retrieval" else "Adaptive tutor",cex.main=1)
  abline(h=0,col="gray75");lines(a,pred,type="o",pch=16,col="#006A9A");lines(a,est,type="o",pch=15,col="#B44514")
 }
 par(mar=c(0,0,0,0));plot.new();legend("center",c("Conditional truth","OVB prediction","Monte Carlo mean"),col=c("#333333","#006A9A","#B44514"),lty=c(2,1,1),pch=c(NA,16,15),ncol=3,bty="n",cex=.95)
}
draw_recovery<-function(){
 setup()
 for(arch in c("spaced","tutor")) {
  d<-subset(s,architecture==arch);keep<-if(arch=="spaced")c("marginal","recorded_belief","policy_log","oracle_state") else names(colours)
  plot(NA,xlim=range(d$adaptivity),ylim=range(c(d$mean_truth,d$mean_estimate,0)),xlab="Policy adaptivity",ylab="Standardized coefficient",main=if(arch=="spaced")"Spaced retrieval" else "Adaptive tutor",cex.main=1)
  abline(h=0,col="gray75");truth<-aggregate(mean_truth~adaptivity,d,mean);lines(truth$adaptivity,truth$mean_truth,lty=2,lwd=1.5)
  for(nm in keep){z<-d[d$estimator==nm,];z<-z[order(z$adaptivity),];lines(z$adaptivity,z$mean_estimate,type="o",pch=match(nm,names(colours))+14,col=colours[[nm]])}
 }
 par(mar=c(0,0,0,0));plot.new();legend("center",c("Conditional truth",labels),col=c("black",colours),lty=c(2,rep(1,5)),pch=c(NA,15:19),ncol=3,bty="n",cex=.95)
}
draw_boundary<-function(){
 setup();z<-cx[cx$context_strength==1 & !is.na(cx$proxy_r2),];z<-z[order(z$proxy_r2),]
 plot(z$proxy_r2,z$bias,type="o",pch=16,col="#B44514",xlab=expression("Context-proxy "*R^2),ylab="Standardized coefficient bias",main="Unrecorded context",cex.main=1);abline(h=0,col="gray75")
 plot(NA,xlim=range(mrt$n_people),ylim=c(.90,.98),xlab="Learner clusters",ylab="95% interval coverage",main="Randomized benchmark",cex.main=1);abline(h=.95,lty=2,col="gray50")
 for(arch in c("spaced","tutor")){z<-mrt[mrt$architecture==arch,];z<-z[order(z$n_people),];lines(z$n_people,z$coverage,type="o",pch=if(arch=="spaced")16 else 15,col=if(arch=="spaced")"#006A9A" else "#1B7959")}
 par(mar=c(0,0,0,0));plot.new();legend("center",c("Spaced retrieval","Adaptive tutor","Nominal coverage"),col=c("#006A9A","#1B7959","gray50"),lty=c(1,1,2),pch=c(16,15,NA),ncol=3,bty="n",cex=.95)
}
for(nm in c("prediction","recovery","boundary")) {
 file<-c(prediction="figure_analytic_prediction_v5",recovery="figure_estimator_recovery_v5",boundary="figure_assumption_boundary_v5")[[nm]];draw<-get(paste0("draw_",nm))
 pdf(file.path(dest,paste0(file,".pdf")),width=7.2,height=4.7,pointsize=12);draw();dev.off()
 svg(file.path(dest,paste0(file,".svg")),width=7.2,height=4.7,pointsize=12);draw();dev.off()
 png(file.path(dest,paste0(file,".png")),width=7.2,height=4.7,units="in",res=300,type="cairo",pointsize=12);draw();dev.off()
}
