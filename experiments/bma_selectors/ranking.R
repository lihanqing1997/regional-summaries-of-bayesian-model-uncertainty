# Source after core/production.R, with working directory set to core/.
Rcpp::sourceCpp("../ranked_prefix.cpp",cacheDir="build",showOutput=FALSE)
ranking_search <- function(d,condition,rule=c("ars","mass","mass_cost")) {
  rule<-match.arg(rule)
  if(rule=="ars")return(production_search(d,condition))
  saved_state<-ars_base_state;saved_advance<-advance_envelope
  estimates<-new.env(hash=TRUE,parent=emptyenv())
  on.exit({assign("ars_base_state",saved_state,.GlobalEnv)
    assign("advance_envelope",saved_advance,.GlobalEnv)},add=TRUE)
  assign("advance_envelope",function(record,d,node,score,remaining,alpha,grid=adaptive_grid) {
    fit<-saved_advance(record,d,node,score,remaining,alpha,grid)
    if(!isTRUE(fit$blocked))assign(node$id,fit$log_estimate,estimates)
    fit
  },.GlobalEnv)
  assign("ars_base_state",function(nodes,frontier_ids,epsilon) {
    # Estimates affect ordering only. Bounds alone decide feasibility.
    estimate<-vapply(nodes[frontier_ids],function(node) {
      if(exists(node$id,estimates,inherits=FALSE))return(get(node$id,estimates))
      z<-node$evidence_fit$log_evidence
      if(!is.null(z)&&is.finite(z))return(z)
      rrs_log_sum_exp(c(node$log_lower,node$log_upper))-log(2)
    },numeric(1))
    selector<-function(log_evidence_lower,costs,epsilon,
      log_evidence_upper=log_evidence_lower,mode=c("conservative","optimistic"),eligible=NULL) {
      mode<-match.arg(mode)
      if(mode=="optimistic")return(select_region_union_epsilon(log_evidence_lower,costs,
        epsilon,log_evidence_upper,mode,eligible))
      v<-rrs_validate_interval_frontier_inputs(log_evidence_lower,log_evidence_upper,costs,epsilon,eligible)
      ids<-which(v$eligible)
      priority<-estimate-if(rule=="mass_cost")log(costs) else 0
      ordered<-ids[order(-priority[ids],ids)]
      selected<-verified_ranked_prefix(v$log_lower,v$log_upper,as.integer(ordered),epsilon)
      interval<-if(length(selected))rrs_interval_omitted_mass(v$log_lower,v$log_upper,selected) else
        list(log_lower_selected=-Inf,log_upper_selected=-Inf,lower=1,upper=1)
      list(feasible=length(selected)>0,selected=selected,budget=if(length(selected))sum(costs[selected]) else Inf,
        solver_status="verified-ranked-prefix",precision_bits=256,
        log_evidence_lower=interval$log_lower_selected,log_evidence_upper=interval$log_upper_selected,
        omitted_mass_lower=interval$lower,omitted_mass_upper=interval$upper,epsilon=epsilon,mode=mode)
    }
    f<-saved_state;scope<-new.env(parent=environment(f))
    scope$select_region_union_epsilon<-selector;environment(f)<-scope
    f(nodes,frontier_ids,epsilon)
  },.GlobalEnv)
  production_search(d,condition)
}
