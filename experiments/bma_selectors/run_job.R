args<-commandArgs(TRUE)
stopifnot(length(args)==1L)
job_id<-args[1L]
setwd("core")
source("production.R")
source("../ranking.R")
jobs<-read.csv("../jobs.csv",stringsAsFactors=FALSE)
job<-jobs[jobs$job==job_id,,drop=FALSE];stopifnot(nrow(job)==1L)
out<-file.path("../results",job_id);dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"complete.json")))quit(status=0)
if(job$kind=="simulation") {
  cases<-read.csv("config/cases.csv",stringsAsFactors=FALSE)
  d<-production_data(cases[job$case_index,,drop=FALSE],job$data_seed)
  seed<-job$data_seed+10110L
} else {
  a<-read.csv("../data/tecator.csv");train<-which(a$split=="train");test<-which(a$split=="test")
  channels<-as.integer(round(seq(1,100,length.out=16)))
  X<-as.matrix(a[,3L+channels]);xm<-colMeans(X[train,]);xs<-apply(X[train,],2,sd)
  ym<-mean(a$fat_percent[train]);ys<-sd(a$fat_percent[train])
  X<-sweep(sweep(X,2,xm),2,xs,"/");y<-(a$fat_percent-ym)/ys
  d<-list(X_train=X[train,],y_train=y[train],X_test=X[test,],y_test=y[test],
    groups=split(seq_len(16),rep(seq_len(4),each=4)),depth=4L,sigma=job$sigma,
    channels=channels,response_mean=ym,response_sd=ys,x_mean=xm,x_sd=xs,
    train_rows=train,test_rows=test)
  seed<-13001001L
}
atomic_rds(d,file.path(out,"data.rds"))
methods<-c("ars","mass","mass_cost")
order<-methods[((seq_along(methods)-1L+job$order_index)%%3)+1L]
answers<-list()
for(method in order) {
  path<-file.path(out,paste0(method,".rds"))
  condition<-list(method="adaptive_kl",epsilon=job$epsilon,budget=job$budget,
    seconds=job$seconds,algorithm_seed=seed,tree_order="forward")
  cat(format(Sys.time()),job_id,method,"starting\n");flush.console()
  ans<-if(file.exists(path))readRDS(path) else ranking_search(d,condition,method)
  if(!file.exists(path))atomic_rds(ans,path)
  # Preserve a numerical proposal-underflow failure as an outcome. The raw
  # checkpoint retains its original status and message for diagnosis.
  if(identical(ans$error,"Invalid screened configuration table."))ans$status<-"numerical-limit"
  if(ans$status %in% c("implementation-error","evidence_interval_conflict"))stop(ans$error)
  answers[[method]]<-ans
}
# No exhaustive reference exists in this process until every practical run ends.
load_production_audit()
source("../real_audit.R")
ref<-if(job$kind=="real")real_reference(d) else make_production_reference(d)
atomic_rds(ref,file.path(out,"reference.rds"))
pilot$epsilon<-job$epsilon
rows<-list()
for(method in methods) {
  ans<-answers[[method]];audit<-audit_pilot_search(ans,ref);looks<-audit_all_looks(ans,ref)
  atomic_rds(audit,file.path(out,paste0(method,"_audit.rds")))
  if(audit$summary$status=="violation" || (nrow(looks)&&!all(looks$covered)))stop("Audit failure")
  s<-ans$snapshot$state;v<-audit$summary
  rows[[method]]<-data.frame(job=job_id,kind=job$kind,case_index=job$case_index,
    replication=job$replication,sigma=job$sigma,epsilon=job$epsilon,method=method,
    status=ans$status,seconds=ans$total_seconds,calls=ans$total_score_calls,
    cost_lower=s$cost_lower%||%0,cost=s$cost_upper%||%Inf,exact_cost=v$exact_optimal_cost,
    bound=s$omitted_mass_bound%||%1,mass=v$retained_mass,feasible=v$mass_guarantee_verified,
    certified=v$verified_optimal,actual_optimal=v$actual_cost_optimal,
    models=v$selected_models,pip_error=v$maximum_pip_error,
    regions=sum(ans$nodes$selected),looks=nrow(ans$sampling_detail),
    audit=v$status,looks_covered=if(nrow(looks))all(looks$covered) else TRUE)
}
write.csv(do.call(rbind,rows),file.path(out,"outcomes.csv"),row.names=FALSE)
write_json(list(job=job_id,outcomes=do.call(rbind,rows),session=capture.output(sessionInfo())),
  file.path(out,"complete.json"))
cat(job_id,"completed\n")
