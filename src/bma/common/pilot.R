source("load_repair.R")
pilot<-jsonlite::fromJSON("config/pilot.json")
cfg$search_seconds_limit<-pilot$search_seconds
# Extract the original runner unchanged, without executing its old-case loop.
expressions<-parse("frozen_runner_source.R")
runner<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name("<-"))&&
  identical(x[[2]],as.name("run_one")),as.list(expressions))
stopifnot(length(runner)==1L);eval(runner[[1]],.GlobalEnv)
write_json<-function(x,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp<-paste0(path,".writing")
  jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,na="null",digits=16)
  if(file.exists(path))file.remove(path)
  if(!file.rename(tmp,path))stop("Cannot commit JSON: ",path)
}
verify_freeze<-function() {
  frozen<-jsonlite::fromJSON("status/freeze.json")
  if(unname(tools::md5sum("config/frozen_manifest.csv"))!=frozen$manifest_md5)
    stop("Pilot manifest changed after freezing")
  manifest<-read.csv("config/frozen_manifest.csv",stringsAsFactors=FALSE)
  got<-unname(tools::md5sum(manifest$path))
  if(anyNA(got)||!identical(tolower(got),tolower(manifest$md5)))stop("Frozen pilot file changed")
}
pilot_cases<-function() {
  data.frame(case_id=c("pairs_p12_b005","pairs_p16_b005","pairs_p20_b005",
    "pairs_p16_b020","triples_p12_b005","quartets_p12_b005",
    "large_quartets_p20_b005","large_quartets_p20_b020"),
    n=80L,p=c(12L,16L,20L,16L,12L,12L,20L,20L),
    sizes=c("2-2-2-2-2-2","2-2-2-2-2-2-2-2","2-2-2-2-2-2-2-2-2-2",
      "2-2-2-2-2-2-2-2","3-3-3-3","4-4-4","4-4-4-4-4","4-4-4-4-4"),
    between=c(.05,.05,.05,.2,.05,.05,.05,.2),design="random",active="1-2-3",
    signal="single",stress=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,TRUE,TRUE))
}
pilot_data<-function(case,seed) {
  d<-generate_data(case,seed)
  if(isTRUE(case$stress)) {
    p<-ncol(d$X_train);S<-matrix(case$between,p,p)
    for(g in d$groups)S[g,g]<-cfg$within_correlation
    diag(S)<-1
    old_signal<-as.numeric(d$X_train%*%d$beta)
    old_signal<-old_signal-mean(old_signal)
    target<-as.numeric(crossprod(d$beta,S%*%d$beta))
    beta<-numeric(p);beta[vapply(d$groups,`[`,integer(1),1L)]<-c(1.4,-1.12,.98,-.84,.7)
    beta<-beta*sqrt(target/as.numeric(crossprod(beta,S%*%beta)))
    signal<-as.numeric(d$X_train%*%beta);signal<-signal-mean(signal)
    d$y_train<-d$y_train-old_signal+signal
    d$beta<-beta;d$signal_variance<-as.numeric(crossprod(beta,S%*%beta))
  }
  d
}
region_cardinality<-function(region) {
  sizes<-lengths(region$groups)
  prod(vapply(seq_along(sizes),function(j)switch(region$states[j],
    free=2^sizes[j],zero=1,one=sizes[j],multi=2^sizes[j]-sizes[j]-1),numeric(1)))
}
# Observation only: preserve the evaluator's arguments, random draws and result.
instrument_search<-function(d,seed) {
  original<-ars_make_hybrid_evaluator;events<-list()
  wrapper<-function(log_score_fn,n_declared_nodes,cap=256L) {
    evaluator<-original(log_score_fn,n_declared_nodes,cap)
    function(node,remaining_score_budget=Inf) {
      start<-proc.time()[["elapsed"]]
      fit<-tryCatch(evaluator(node,remaining_score_budget),error=function(e)e)
      count<-region_cardinality(node$region)
      outcome<-if(inherits(fit,"error"))"error" else if(isTRUE(fit$preflight_failed))
        "preflight-blocked" else if(isTRUE(fit$budget_exhausted))"budget-blocked" else
        if(identical(fit$fit$method,"exact_regional_sum"))"exact" else "sampled"
      events[[length(events)+1L]]<<-data.frame(node_id=node$id,
        states=paste(node$region$states,collapse="/"),regional_models=count,
        exceeds_sum_cap=count>cap,rho=node$rho,outcome=outcome,
        score_calls=if(inherits(fit,"error"))NA_real_ else fit$score_evaluations%||%0,
        required_paths=if(inherits(fit,"error"))NA_real_ else fit$preflight$required_paths%||%NA_real_,
        relative_radius=if(inherits(fit,"error"))NA_real_ else fit$fit$relative_radius%||%NA_real_,
        seconds=proc.time()[["elapsed"]]-start)
      if(inherits(fit,"error"))stop(fit)
      fit
    }
  }
  assign("ars_make_hybrid_evaluator",wrapper,.GlobalEnv)
  on.exit(assign("ars_make_hybrid_evaluator",original,.GlobalEnv))
  result<-run_one(d,"hybrid_gaussian",seed)
  result$evidence_events<-if(length(events))do.call(rbind,events) else data.frame()
  # Save data-only checkpoints; closures can retain large unused score caches.
  result$trace<-result$search$trace
  if(!is.null(result$search))result$search$nodes<-NULL
  if(!is.null(result$snapshot))result$snapshot$nodes<-NULL
  result
}
