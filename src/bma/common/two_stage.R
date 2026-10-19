source("integration.R")
original_search_engine<-certified_tree_region_search
two_stage_search<-function(d,seed) {
  events<-new.env(parent=emptyenv());events$rows<-list()
  previous_delta<-cfg$delta
  # Split the joint error allowance between the initial search and refinement.
  cfg$delta<<-previous_delta/2
  on.exit({cfg$delta<<-previous_delta
    assign("certified_tree_region_search",original_search_engine,.GlobalEnv)},add=TRUE)
  engine<-function(tree,log_score_fn,evidence_evaluator,epsilon,cost_tolerance=0,
    max_score_evaluations=Inf,max_actions=Inf,bound_node=NULL) {
    base<-original_search_engine(tree,log_score_fn,evidence_evaluator,epsilon,cost_tolerance,
      max_score_evaluations,max_actions,bound_node=bound_node)
    if(base$certified || base$status=="evidence_interval_conflict")return(base)
    nodes<-base$nodes;frontier<-base$frontier_ids
    initial<-ars_tracker$state;solver<-ars_new_state_solver()
    environment(solver)$best_incumbent<-initial
    environment(solver)$best_lower<-initial$cost_lower
    state<-initial;used<-base$total_score_evaluations;additional<-0;attempted<-character()
    evaluator<-make_revised_factory(d,events)(log_score_fn,tree$n_declared_nodes)
    status<-base$status
    repeat {
      if(proc.time()[["elapsed"]]>ars_tracker$deadline){status<-"time-limit";break}
      ids<-frontier[vapply(nodes[frontier],function(n)isTRUE(n$terminal)&&isTRUE(n$preflight_blocked),logical(1))]
      ids<-setdiff(ids,attempted)
      if(!length(ids)){status<-"frontier_exhausted";break}
      priority<-vapply(nodes[ids],function(n)n$log_upper-log(n$region$cost),numeric(1))
      id<-ids[which.max(priority)];node<-nodes[[id]];attempted<-c(attempted,id)
      fit<-evaluator(node,max_score_evaluations-used)
      used<-used+fit$score_evaluations;additional<-additional+fit$score_evaluations
      if(isTRUE(fit$budget_exhausted)){status<-"budget_exhausted";break}
      if(isTRUE(fit$preflight_failed))next
      lo<-max(node$log_lower,fit$log_lower);hi<-min(node$log_upper,fit$log_upper)
      if(lo>hi+1e-10){status<-"evidence_interval_conflict";break}
      interval<-rrs_reconcile_log_interval(lo,hi)
      node$log_lower<-interval[["log_lower"]];node$log_upper<-interval[["log_upper"]]
      node$evidence_fit<-fit$fit;node$evidence_preflight<-fit$preflight
      node$preflight_blocked<-FALSE;node$evaluated<-TRUE;nodes[[id]]<-node
      state<-solver(nodes,frontier,epsilon)
      if(state$cost_upper>initial$cost_upper)stop("Refinement lost the valid initial solution")
      if(isTRUE(state$incumbent$feasible)&&state$omitted_mass_bound<=epsilon&&state$cost_gap<=cost_tolerance){
        status<-"certified";break
      }
    }
    base$nodes<-nodes;base$status<-status;base$registered_status<-rrs_registered_status(status)
    base$certified<-base$bounds_met<-identical(status,"certified")
    base$epsilon_feasible<-isTRUE(state$incumbent$feasible)&&state$omitted_mass_bound<=epsilon
    base$cost_optimality_met<-base$epsilon_feasible&&state$cost_gap<=cost_tolerance
    for(name in c("selected_ids","cost_lower","cost_upper","cost_gap","bracket_consistent","omitted_mass_bound"))base[[name]]<-state[[name]]
    base$total_score_evaluations<-used;base$evidence_score_evaluations<-base$evidence_score_evaluations+additional
    base$actions<-base$actions+length(attempted)
    base$refinement<-list(initial_cost=initial$cost_upper,final_cost=state$cost_upper,
      additional_calls=additional,attempted_regions=length(attempted),joint_delta=previous_delta)
    base
  }
  assign("certified_tree_region_search",engine,.GlobalEnv)
  ans<-instrument_search(d,seed)
  ans$sampling_detail<-if(length(events$rows))do.call(rbind,events$rows) else data.frame()
  ans
}
