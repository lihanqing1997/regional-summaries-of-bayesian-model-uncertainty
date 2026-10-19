# Keep a previously verified certificate when a later partition has weaker
# bounds. The family of allowable terminal unions is fixed throughout a run.
ars_base_state<-rrs_tree_epsilon_search_state
ars_tracker<-NULL
ars_new_state_solver<-function() {
  previous_input<-NULL;previous_state<-NULL;best_incumbent<-NULL;best_lower<-0
  function(nodes,frontier_ids,epsilon) {
    if(!is.null(ars_tracker$deadline) && proc.time()[["elapsed"]]>ars_tracker$deadline)
      stop("Reached search time limit at an action boundary",call.=FALSE)
    t<-proc.time()[["elapsed"]]
    input<-lapply(nodes[frontier_ids],function(n)
      list(n$id,n$log_lower,n$log_upper,n$evaluated,n$terminal,n$cost_lower,n$region$cost))
    if(identical(input,previous_input)) {
      # Ancestor-only propagation can leave the frontier selection input intact.
      if(!is.null(ars_tracker)) {
        ars_tracker$nodes<-nodes;ars_tracker$frontier_ids<-frontier_ids
        ars_tracker$state<-previous_state
      }
      return(previous_state)
    }
    s<-ars_base_state(nodes,frontier_ids,epsilon)
    if(!is.null(best_incumbent) && best_incumbent$cost_upper<=s$cost_upper) {
      if(best_incumbent$cost_upper<s$cost_upper ||
         best_incumbent$omitted_mass_bound<s$omitted_mass_bound) {
        for(nm in c("incumbent","selected_ids","cost_upper","omitted_mass_bound","certificate"))
          s[[nm]]<-best_incumbent[[nm]]
        s$incumbent$selected<-match(s$selected_ids,frontier_ids)
        stopifnot(!anyNA(s$incumbent$selected))
      }
    }
    if(isTRUE(s$incumbent$feasible) && is.null(s$certificate))
      s$certificate<-list(frontier_ids=frontier_ids,
        log_lower=vapply(nodes[frontier_ids],`[[`,numeric(1),"log_lower"),
        log_upper=vapply(nodes[frontier_ids],`[[`,numeric(1),"log_upper"),
        selected_ids=s$selected_ids,omitted_mass_bound=s$omitted_mass_bound)
    if(isTRUE(s$incumbent$feasible))best_incumbent<<-s
    best_lower<<-max(best_lower,s$cost_lower)
    s$cost_lower<-best_lower
    s$bracket_consistent<-!is.finite(s$cost_upper)||best_lower<=s$cost_upper
    if(!s$bracket_consistent)stop("Cost bounds conflict; no certificate returned")
    s$cost_gap<-s$cost_upper-best_lower
    previous_input<<-input;previous_state<<-s
    # This snapshot is committed only after the two selection programs finish.
    if(!is.null(ars_tracker)) {
      ars_tracker$nodes<-nodes;ars_tracker$frontier_ids<-frontier_ids;ars_tracker$state<-s
      ars_tracker$selection_seconds<-ars_tracker$selection_seconds+proc.time()[["elapsed"]]-t
      ars_tracker$selection_calls<-ars_tracker$selection_calls+1L
    }
    s
  }
}
ars_base_search<-ars_add_mass_conservation(certified_tree_region_search)
certified_tree_region_search<-function(...,bound_node=NULL) {
  engine<-ars_base_search
  scope<-new.env(parent=environment(engine))
  scope$rrs_tree_epsilon_search_state<-ars_new_state_solver()
  if(!is.null(bound_node))scope$rrs_bound_tree_node<-bound_node
  environment(engine)<-scope
  engine(...)
}
select_certified_region_tree<-certified_tree_region_search
search_region_tree<-certified_tree_region_search

# Enumerate only the configurations in one requested region. This does not
# inspect full-posterior enumeration or a posterior-dependent reference table.
ars_regional_supports<-function(region,cap=256L) {
  tables<-region$configuration_tables
  if(is.null(tables))return(NULL)
  allowed<-lapply(seq_along(tables),function(j) {
    k<-rowSums(tables[[j]]$supports)
    switch(region$states[j],free=seq_along(k),zero=which(k==0),
           one=which(k==1),multi=which(k>=2))
  })
  count<-prod(lengths(allowed))
  if(count>cap)return(NULL)
  index<-as.matrix(expand.grid(allowed,KEEP.OUT.ATTRS=FALSE))
  supports<-matrix(0L,nrow=count,ncol=sum(lengths(region$groups)))
  for(j in seq_along(tables))supports[,region$groups[[j]]]<-tables[[j]]$supports[index[,j],,drop=FALSE]
  supports
}
ars_make_hybrid_evaluator<-function(log_score_fn,n_declared_nodes,cap=256L) {
  iid<-make_iid_bounded_evidence_evaluator(log_score_fn,n_declared_nodes,
    delta=cfg$delta,n_paths=cfg$minimum_paths,
    maximum_relative_radius=cfg$maximum_relative_radius,maximum_paths=cfg$maximum_paths)
  function(node,remaining_score_budget=Inf) {
    supports<-ars_regional_supports(node$region,cap)
    if(is.null(supports))return(iid(node,remaining_score_budget))
    count<-nrow(supports)
    pf<-regional_evidence_preflight(node$rho,cfg$maximum_relative_radius,
      cfg$delta/n_declared_nodes,remaining_paths=cfg$maximum_paths)
    # Direct summation is used only when it is cheaper than the scheduled
    # Monte Carlo calculation, or the range preflight rejects that calculation.
    if(isTRUE(pf$feasible) && count>max(cfg$minimum_paths,pf$required_paths))
      return(iid(node,remaining_score_budget))
    if(count>remaining_score_budget)return(list(budget_exhausted=TRUE,score_evaluations=0))
    bounded<-rrs_budgeted_log_score(log_score_fn,remaining_score_budget)
    logs<-rrs_evaluate_log_score(bounded,supports)
    z<-rrs_log_sum_exp(logs)
    used<-attr(bounded,"score_budget_usage")()
    stopifnot(used==count)
    fit<-list(method="exact_regional_sum",log_evidence=z,log_lower=z,log_upper=z,
              score_evaluations=used,regional_support_count=count,relative_radius=0)
    list(log_evidence=z,log_lower=z,log_upper=z,score_evaluations=used,fit=fit,
         preflight=list(required_paths=0,scheduled_paths=count,method="exact_regional_sum"))
  }
}

# A small internal region can be summed once and reused by all its descendants.
# The cache contains only supports actually evaluated under the search budget.
ars_make_cached_bounder<-function(score,budget,cap=256L,split_order) {
  cache<-new.env(hash=TRUE,parent=emptyenv())
  function(node,log_score_fn) {
    states<-node$region$states;ancestor<-states;record<-NULL
    repeat {
      key<-paste(ancestor,collapse="/")
      if(exists(key,cache,inherits=FALSE)){record<-get(key,cache);break}
      assigned<-split_order[ancestor[split_order]!="free"]
      if(!length(assigned))break
      ancestor[tail(assigned,1)]<-"free"
    }
    calls<-0L
    if(is.null(record)) {
      remaining<-budget-unname(score$diagnostics()[["calls"]])
      supports<-ars_regional_supports(node$region,min(cap,remaining))
      if(is.null(supports))return(rrs_bound_tree_node(node,log_score_fn))
      logs<-rrs_evaluate_log_score(rrs_budgeted_log_score(log_score_fn,remaining),supports)
      calls<-nrow(supports)
      counts<-vapply(node$region$groups,function(g)pmin(2L,rowSums(supports[,g,drop=FALSE])),
                     numeric(nrow(supports)))
      if(length(node$region$groups)==1)counts<-matrix(counts,ncol=1)
      record<-list(logs=logs,counts=counts)
      assign(paste(states,collapse="/"),record,cache)
    }
    inside<-rep(TRUE,length(record$logs))
    for(j in which(states!="free"))inside<-inside &
      record$counts[,j]==match(states[j],c("zero","one","multi"))-1L
    stopifnot(any(inside))
    z<-rrs_log_sum_exp(record$logs[inside])
    node$log_lower<-z;node$log_upper<-z;node$bound_score_evaluations<-calls
    node$anchor<-NULL
    node$evaluated<-isTRUE(node$terminal)
    node$evidence_fit<-list(method=if(calls>0)"exact_local_subtree" else "exact_cached_subtree",
      log_evidence=z,log_lower=z,log_upper=z,score_evaluations=calls,
      regional_support_count=sum(inside),relative_radius=0)
    node
  }
}
