# Independent checks of the distributed summaries, refits and minimized public data.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
root<-normalizePath(args[1]);out<-normalizePath(args[2]);checks<-list()
add<-function(name,ok,detail){stopifnot(length(ok)==1L);checks[[length(checks)+1L]]<<-data.frame(check=name,passed=isTRUE(ok),detail=as.character(detail))}
read<-function(p)read.csv(file.path(root,p),stringsAsFactors=FALSE)
for(z in list(c('studies/informative_observation/validation/confirmatory_r_comparison.csv',36),c('studies/informative_observation/validation/smoke_r_comparison.csv',24),c('studies/informative_observation/validation/quadrature_comparison.csv',248),c('validation/pearson_cr2_comparison.csv',153))){
 x<-read(z[1]);add(z[1],nrow(x)==as.integer(z[2])&&all(as.logical(x$pass)),nrow(x))
}
x<-read('studies/informative_observation/validation/summary_independent_checks.csv')
add('archived_independent_summary',nrow(x)==156&&max(abs(x$difference))<1e-10,max(abs(x$difference)))
refit<-read.csv(file.path(out,'pearson_refits.csv'))
ok<-all(refit$converged)
for(metric in c('estimate','se','df')){
 delta<-abs(refit[[paste0('julia_',metric)]]-refit[[paste0('r_',metric)]])
 lim<-if(metric=='df')1e-5 else 1e-6
 add(paste0('current_Pearson_',metric),all(delta<=lim+1e-5*abs(refit[[paste0('julia_',metric)]])),max(delta))
}
add('current_refit_count_and_convergence',nrow(refit)==51&&ok,nrow(refit))
source(file.path(root,'empirical_workflow/src/R/theme_brm.R'))
source(file.path(root,'empirical_workflow/src/R/validate_public_data.R'))
v<-brm_validate_public_data(file.path(root,'empirical_workflow'),write_report=FALSE)
write.csv(v$checks,file.path(out,'public_data_validation.csv'),row.names=FALSE)
add('minimized_public_data_checks',all(v$checks$passed),nrow(v$checks))
z<-do.call(rbind,checks);write.csv(z,file.path(out,'public_evidence_checks.csv'),row.names=FALSE)
stopifnot(all(z$passed));cat(nrow(z),'public evidence checks passed.\n')
writeLines(capture.output(sessionInfo()),file.path(out,'R_session.txt'))
