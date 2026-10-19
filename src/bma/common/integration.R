source("sampling.R")
original_hybrid_factory<-ars_make_hybrid_evaluator
make_revised_factory<-function(d,events) {
  function(log_score_fn,n_declared_nodes,cap=256L) {
    small<-original_hybrid_factory(log_score_fn,n_declared_nodes,cap)
    function(node,remaining_score_budget=Inf) {
      count<-region_cardinality(node$region)
      if(count<=cap)return(small(node,remaining_score_budget))
      sizes<-lengths(node$region$groups)
      options<-vapply(seq_along(sizes),function(j)switch(node$region$states[j],
        free=2^sizes[j],zero=1,one=sizes[j],multi=2^sizes[j]-1-sizes[j]),numeric(1))
      cells<-1
      for(k in options){if(cells*k>256)break;cells<-cells*k}
      n<-min(2500,remaining_score_budget)-cells
      if(n<100)return(list(budget_exhausted=TRUE,score_evaluations=0))
      bounded<-rrs_budgeted_log_score(log_score_fn,remaining_score_budget)
      proposal<-make_partitioned_envelope(d,bounded,node$region$states,256L)
      prep<-attr(bounded,"score_budget_usage")()
      stopifnot(prep==cells)
      fit<-sample_envelope_evidence(proposal,bounded,n,cfg$delta/n_declared_nodes)
      fit$score_evaluations<-attr(bounded,"score_budget_usage")()
      stopifnot(fit$score_evaluations==prep+n,fit$score_evaluations<=remaining_score_budget)
      fit$fit$preparation_calls<-prep;fit$fit$sampling_calls<-n
      fit$preflight<-list(required_paths=NA_real_,scheduled_paths=n)
      events$rows[[length(events$rows)+1L]]<-data.frame(node=node$id,states=paste(node$region$states,collapse="/"),
        regional_models=count,cells=cells,preparation_calls=prep,sampling_calls=n,
        relative_radius=fit$fit$relative_radius,positive_lower=is.finite(fit$log_lower),
        log_lower=fit$log_lower,log_upper=fit$log_upper,log_estimate=fit$log_evidence)
      fit
    }
  }
}
revised_search<-function(d,seed) {
  events<-new.env(parent=emptyenv());events$rows<-list()
  assign("ars_make_hybrid_evaluator",make_revised_factory(d,events),.GlobalEnv)
  on.exit(assign("ars_make_hybrid_evaluator",original_hybrid_factory,.GlobalEnv))
  ans<-instrument_search(d,seed)
  ans$sampling_detail<-if(length(events$rows))do.call(rbind,events$rows) else data.frame()
  ans
}
