source("production.R")
verify_production()
arg<-args_named();jobs<-read.csv("config/jobs.csv");job<-jobs[jobs$job==arg$job,,drop=FALSE]
stopifnot(nrow(job)==1L)
cases<-read.csv("config/cases.csv");case<-cases[cases$case_id==job$case_id,,drop=FALSE]
conditions<-read.csv("config/conditions.csv");conditions<-conditions[conditions$job==job$job,,drop=FALSE]
root<-file.path("results",job$job);dir.create(root,recursive=TRUE,showWarnings=FALSE)
rows<-list()
event<-function(stage,condition="") {
  write_json(list(job=job$job,stage=stage,condition=condition,updated=as.character(Sys.time())),
    "status/worker.json")
  cat(format(Sys.time()),stage,condition,"\n");flush.console()
}
save_rows<-function()write_json(list(job=job$job,case_id=job$case_id,replication=job$replication,
  outcomes=unname(rows)),file.path(root,"outcomes.json"))
fatal<-function(reason) {
  write_json(list(job=job$job,reason=reason,time=as.character(Sys.time())),"status/AUDIT_FAILURE.json")
  stop(reason)
}
tryCatch({
  dpath<-file.path(root,"data.rds")
  d<-if(file.exists(dpath))readRDS(dpath) else production_data(case,job$data_seed)
  if(!file.exists(dpath))atomic_rds(d,dpath)
  write_json(data_diagnostics(d),file.path(root,"data_diagnostics.json"))
  for(i in seq_len(nrow(conditions))) {
    cnd<-conditions[i,];path<-file.path(root,paste0(cnd$condition,".rds"))
    event("practical-search",cnd$condition)
    ans<-if(file.exists(path))readRDS(path) else tryCatch(production_search(d,cnd),
      error=function(e)list(status=error_status(e),error=conditionMessage(e)))
    if(!file.exists(path))atomic_rds(ans,path)
    rows[[cnd$condition]]<-search_row(cnd,ans);save_rows()
    if(!is.null(ans$total_score_calls)&&ans$total_score_calls>cnd$budget)fatal("Score budget exceeded")
    if(identical(ans$status,"evidence_interval_conflict"))fatal("Regional interval conflict")
    rm(ans);gc(FALSE)
  }
  event("audit-reference")
  load_production_audit()
  refpath<-file.path(root,"reference.rds")
  ref<-tryCatch(if(file.exists(refpath))readRDS(refpath) else make_production_reference(d),error=function(e)e)
  if(inherits(ref,"error")) {
    for(id in names(rows)){rows[[id]]$audit_status<-"reference-unavailable";rows[[id]]$audit_error<-conditionMessage(ref)}
    save_rows()
  } else {
    if(!file.exists(refpath))atomic_rds(ref,refpath)
    comparators<-production_comparators(ref,c(.01,.05,.1,.2))
    write.csv(comparators,file.path(root,"comparators.csv"),row.names=FALSE)
    for(i in seq_len(nrow(conditions))) {
      cnd<-conditions[i,];event("audit-search",cnd$condition)
      ans<-readRDS(file.path(root,paste0(cnd$condition,".rds")))
      pilot$epsilon<-cnd$epsilon
      checked<-tryCatch(limited(function() {
        audit<-audit_pilot_search(ans,ref);looks<-audit_all_looks(ans,ref)
        list(audit=audit,looks=looks)
      },180),error=function(e)e)
      if(inherits(checked,"error")) {
        rows[[cnd$condition]]$audit_status<-"audit-unavailable"
        rows[[cnd$condition]]$audit_error<-conditionMessage(checked)
      } else {
        atomic_rds(checked,file.path(root,paste0(cnd$condition,"_audit.rds")))
        a<-checked$audit$summary;r<-rows[[cnd$condition]]
        r$audit_status<-a$status;r$exact_cost<-a$exact_optimal_cost;r$retained_mass<-a$retained_mass
        r$verified_feasible<-a$mass_guarantee_verified;r$verified_optimal<-a$verified_optimal
        r$all_looks_covered<-if(nrow(checked$looks))all(checked$looks$covered) else TRUE
        r$maximum_pip_error<-a$maximum_pip_error;r$selected_models<-a$selected_models
        r$reference_seconds<-a$reference_enumeration_seconds
        rows[[cnd$condition]]<-r;save_rows()
        if(a$status=="violation"||!r$all_looks_covered||
           (isTRUE(a$certificate_available)&&!isTRUE(a$mass_guarantee_verified)))
          fatal("A reported certificate or regional interval failed its exact audit")
      }
      save_rows();rm(ans,checked);gc(FALSE)
    }
  }
  write_json(list(job=job$job,status="completed",searches=nrow(conditions),
    finished=as.character(Sys.time())),file.path(root,"complete.json"))
  event("dataset-completed")
},error=function(e) {
  save_rows()
  write_json(list(job=job$job,status="worker-failure",error=conditionMessage(e),
    finished=as.character(Sys.time())),file.path(root,"worker_failure.json"))
  stop(e)
})
