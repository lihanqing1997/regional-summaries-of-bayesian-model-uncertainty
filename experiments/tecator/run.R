args<-commandArgs(TRUE)
setwd("core")
source("production.R")
cfg$mass_conservation<-FALSE
cfg$envelope_optimizer<-"newton_span"
outdir<-"../validation/results/numerical_end_to_end"
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

a<-read.csv("../data/tecator.csv");train<-which(a$split=="train");test<-which(a$split=="test")
channels<-as.integer(round(seq(1,100,length.out=16)))
X<-as.matrix(a[,3L+channels]);xm<-colMeans(X[train,]);xs<-apply(X[train,],2,sd)
ym<-mean(a$fat_percent[train]);ys<-sd(a$fat_percent[train])
X<-sweep(sweep(X,2,xm),2,xs,"/");y<-(a$fat_percent-ym)/ys
datasets<-list()
for(sigma in c(.25,.5,1)) datasets[[paste0("tecator_sigma_",sigma)]]<-list(
 X_train=X[train,],y_train=y[train],X_test=X[test,],y_test=y[test],
 groups=split(seq_len(16),rep(seq_len(4),each=4)),depth=4L,sigma=sigma,
 channels=channels,response_mean=ym,response_sd=ys,x_mean=xm,x_sd=xs,
 train_rows=train,test_rows=test)

jobs<-do.call(rbind,lapply(names(datasets),function(name)data.frame(
 dataset=name,epsilon=if(grepl("tecator",name))c(.1,.05) else .1,
 stringsAsFactors=FALSE)))
write.csv(jobs,file.path(outdir,"planned_runs.csv"),row.names=FALSE)
for(i in seq_len(nrow(jobs))) {
 name<-paste0(jobs$dataset[i],"_eps_",jobs$epsilon[i]);d<-datasets[[jobs$dataset[i]]]
 path<-file.path(outdir,paste0(name,".rds"))
 if(file.exists(path))next
 condition<-list(method="adaptive_kl",epsilon=jobs$epsilon[i],budget=50000L,
   seconds=180,algorithm_seed=13001001L+i,tree_order="forward")
 cat(name,"starting\n");flush.console()
 ans<-production_search(d,condition)
 atomic_rds(list(version="method_revision_20261002",mass_conservation=FALSE,
   optimizer="newton_span",data=d,condition=condition,result=ans),path)
 cat(name,ans$status,"calls",ans$total_score_calls,"seconds",ans$total_seconds,"\n");flush.console()
 if(ans$status %in% c("implementation-error","evidence_interval_conflict"))stop(ans$error)
}

# Construct exact audit information only after every practical search is saved.
load_production_audit();source("../real_audit.R")
references<-lapply(names(datasets),function(name) {
 d<-datasets[[name]]
 ref<-if(grepl("tecator",name))real_reference(d) else make_production_reference(d)
 atomic_rds(ref,file.path(outdir,paste0(name,"_reference.rds")));ref
});names(references)<-names(datasets)
rows<-list()
for(i in seq_len(nrow(jobs))) {
 name<-paste0(jobs$dataset[i],"_eps_",jobs$epsilon[i])
 saved<-readRDS(file.path(outdir,paste0(name,".rds")));ans<-saved$result
 pilot$epsilon<-jobs$epsilon[i];ref<-references[[jobs$dataset[i]]]
 audit<-audit_pilot_search(ans,ref);looks<-audit_all_looks(ans,ref)
 stopifnot(audit$summary$status=="passed",!nrow(looks)||all(looks$covered))
 atomic_rds(list(audit=audit,looks=looks),file.path(outdir,paste0(name,"_audit.rds")))
 s<-ans$snapshot$state;v<-audit$summary
 rows[[i]]<-data.frame(dataset=jobs$dataset[i],epsilon=jobs$epsilon[i],
  status=ans$status,calls=ans$total_score_calls,seconds=ans$total_seconds,
  lower=s$cost_lower,cost=s$cost_upper,omitted_bound=s$omitted_mass_bound,
  audit_mass=v$retained_mass,exact_optimum=v$exact_optimal_cost,
  feasible=v$mass_guarantee_verified,certified=v$verified_optimal,
  audit=v$status,looks=nrow(looks),all_looks_covered=if(nrow(looks))all(looks$covered) else TRUE)
}
write.csv(do.call(rbind,rows),file.path(outdir,"outcomes.csv"),row.names=FALSE)
jsonlite::write_json(list(version="method_revision_20261002",runs=do.call(rbind,rows),
  session=capture.output(sessionInfo())),file.path(outdir,"complete.json"),
  pretty=TRUE,auto_unbox=TRUE,digits=16,na="null")
cat("End-to-end numerical validation passed.\n")
