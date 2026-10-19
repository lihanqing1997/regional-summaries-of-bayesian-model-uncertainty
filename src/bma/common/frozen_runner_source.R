run_one<-function(d,mode,seed) {
  ars_tracker<<-new.env(parent=emptyenv())
  ars_tracker$selection_seconds<-0;ars_tracker$selection_calls<-0L
  score<-NULL;prep_calls<-0;search<-NULL
  started<-proc.time()[["elapsed"]]
  ars_tracker$deadline<-started+cfg$search_seconds_limit
  ans<-tryCatch((function() {
    score<<-new_score(d)
    proposal<-make_screened_group_region_factory(score$score,ncol(d$X_train),d$groups)
    bound<-cb_make_practical_group_log_g_range_bound(d$X_train,d$y_train,
      sigma=d$sigma,tau2=cfg$tau2,inclusion_probability=cfg$inclusion_probability,
      groups=d$groups,configuration_tables=proposal$configuration_tables)
    tree<-new_group_state_tree(rep(cfg$inclusion_probability,ncol(d$X_train)),d$groups,
      split_order=seq_along(d$groups),max_depth=d$depth,rho=bound$rho,
      region_factory=proposal$region_factory,target_mutation=make_group_gibbs_target_mutation,
      state_cost=c(zero=1L,one=2L,multi=2L),base_cost=1L)
    prep_calls<<-unname(score$diagnostics()[["calls"]]);score$reset(clear_cache=FALSE)
    evaluator<-if(mode%in%c("hybrid","hybrid_cached","hybrid_gaussian"))ars_make_hybrid_evaluator(score$score,tree$n_declared_nodes) else
      make_iid_bounded_evidence_evaluator(score$score,tree$n_declared_nodes,
        delta=cfg$delta,n_paths=cfg$minimum_paths,maximum_relative_radius=cfg$maximum_relative_radius,
        maximum_paths=cfg$maximum_paths)
    set.seed(seed)
    search<<-certified_tree_region_search(tree,score$score,evaluator,epsilon=.1,
      cost_tolerance=cfg$cost_tolerance,max_score_evaluations=25000-prep_calls,
      max_actions=cfg$action_budget,
      bound_node=if(mode=="hybrid_cached")ars_make_cached_bounder(score,25000-prep_calls,
        cap=256L,split_order=seq_along(d$groups)) else if(mode=="hybrid_gaussian")
        ars_make_gaussian_bounder(d,score,25000-prep_calls) else NULL)
    list(status=search$status,error=NA_character_)
  })(),error=function(e)list(status=error_status(e),error=conditionMessage(e)))
  ans$total_seconds<-proc.time()[["elapsed"]]-started
  ans$mode<-mode;ans$epsilon<-.1;ans$preparation_calls<-prep_calls
  ans$score_diagnostics<-if(!is.null(score))score$diagnostics() else NULL
  ans$total_score_calls<-prep_calls+if(!is.null(score))unname(score$diagnostics()[["calls"]]) else 0
  ans$selection_seconds<-ars_tracker$selection_seconds;ans$selection_calls<-ars_tracker$selection_calls
  ans$search<-search
  if(!is.null(ars_tracker$state)) {
    ans$snapshot<-list(nodes=ars_tracker$nodes,frontier_ids=ars_tracker$frontier_ids,state=ars_tracker$state)
    ans$nodes<-node_table(list(nodes=ars_tracker$nodes,frontier_ids=ars_tracker$frontier_ids,
                             selected_ids=ars_tracker$state$selected_ids))
  }
  if(!is.null(search))stopifnot(ans$total_score_calls==prep_calls+search$total_score_evaluations)
  stopifnot(ans$total_score_calls<=25000)
  ans
}
