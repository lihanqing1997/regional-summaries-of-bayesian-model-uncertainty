scaling_region_log_mass <- function(states,ref) {
  values <- strsplit(states,'/',fixed=TRUE)[[1]]
  depth <- sum(values!='free')
  if(depth)stopifnot(all(values[seq_len(depth)]!='free'))
  if(depth<length(values))stopifnot(all(values[(depth+1L):length(values)]=='free'))
  digits <- match(values[seq_len(depth)],c('zero','one','multi'))-1L
  stopifnot(!anyNA(digits))
  code <- if(depth)sum(digits*3^(seq_len(depth)-1L)) else 0
  ref$prefix_log_mass[[depth+1L]][code+1L]
}

scaling_audit_ars <- function(ans,ref,epsilon,optimum) {
  nodes <- ans$nodes;state <- ans$snapshot$state
  if(is.null(nodes)||!nrow(nodes)||is.null(state))return(list(
    audit_status='no-valid-frontier',feasible=FALSE,attained_optimal=FALSE,
    reported_feasible=FALSE,reported_certified=FALSE,certified_optimal=FALSE))
  exact <- vapply(nodes$states,scaling_region_log_mass,numeric(1),ref=ref)
  tolerance <- 1e-9*pmax(1,abs(exact))
  covered <- nodes$log_lower<=exact+tolerance & nodes$log_upper>=exact-tolerance
  selected <- which(nodes$selected)
  stopifnot(!length(selected)||all(!grepl('free',nodes$states[selected])))
  mass <- sum(exp(exact[selected]-ref$logZ))
  cost <- if(is.finite(state$cost_upper))state$cost_upper else NA_real_
  if(length(selected))stopifnot(sum(nodes$cost[selected])==cost)
  feasible <- length(selected)>0 && mass>=1-epsilon-1e-9
  attained <- feasible && !is.na(cost) && cost==optimum
  events <- ans$sampling_detail
  if(is.null(events)||!nrow(events)) {
    nlooks <- 0L;looks_covered <- TRUE
  } else {
    truth <- vapply(events$states,scaling_region_log_mass,numeric(1),ref=ref)
    tol <- 1e-9*pmax(1,abs(truth))
    nlooks <- nrow(events)
    looks_covered <- all(events$log_lower<=truth+tol & events$log_upper>=truth-tol)
  }
  reported_feasible <- isTRUE(state$incumbent$feasible) && state$omitted_mass_bound<=epsilon
  reported_certified <- identical(ans$status,'certified')
  list(audit_status='completed',cost=cost,lower=state$cost_lower,
       retained_mass=mass,feasible=feasible,attained_optimal=attained,
       reported_feasible=reported_feasible,reported_certified=reported_certified,
       certified_optimal=reported_certified && attained && state$cost_lower==cost,
       cost_bracket_covered=state$cost_lower<=optimum && state$cost_upper>=optimum,
       mass_bound_covered=1-mass<=state$omitted_mass_bound+1e-9,
       all_node_intervals_covered=all(covered),node_intervals=nrow(nodes),
       all_looks_covered=looks_covered,looks=nlooks,
       cost_ratio=cost/optimum,omitted_bound=state$omitted_mass_bound,
       preparation_seconds=ans$timing$proposal_preparation,
       sampling_interval_seconds=ans$timing$sampling_and_intervals,
       selection_seconds=ans$selection_seconds,cache_setup_seconds=ans$cache_setup_seconds,
       frontier_nodes=sum(nodes$frontier),
       terminal_frontier=sum(nodes$frontier & !grepl('free',nodes$states)))
}

scaling_audit_mc3 <- function(ans,ref,epsilon,optimum,method) {
  selected <- ans[[method]]
  if(is.null(selected))return(list(audit_status='no-summary',feasible=FALSE,
    attained_optimal=FALSE,reported_feasible=FALSE,reported_certified=FALSE,certified_optimal=FALSE))
  ids <- as.integer(selected$region_codes)+1L
  stopifnot(!anyDuplicated(ids),all(ids>=1 & ids<=length(ref$costs)),
    sum(ref$costs[ids])==selected$cost)
  mass <- sum(ref$posterior_mass[ids]);feasible <- mass>=1-epsilon-1e-9
  list(audit_status='completed',cost=selected$cost,retained_mass=mass,feasible=feasible,
       attained_optimal=feasible && selected$cost==optimum,
       reported_feasible=NA,reported_certified=FALSE,certified_optimal=FALSE,
       cost_ratio=selected$cost/optimum,estimated_mass=selected$estimated_mass,
       selection_seconds=ans$selection_seconds,cache_setup_seconds=ans$cache_setup_seconds,
       unique_models=ans$unique_models,draws=ans$draws,mean_acceptance=mean(ans$acceptance),
       empirical_selection=selected$status)
}

scaling_empty_row <- function(dataset,job,method,optimum) {
  list(dataset=dataset$dataset,case_index=dataset$case_index,family=dataset$family,
    p=dataset$p,group_size=dataset$group_size,between=dataset$between,
    replication=dataset$replication,budget=job$budget,method=method,job=job$job,
    status=NA_character_,calls=NA_real_,seconds=NA_real_,cost=NA_real_,lower=NA_real_,
    exact_cost=optimum,retained_mass=NA_real_,feasible=FALSE,attained_optimal=FALSE,
    reported_feasible=FALSE,reported_certified=FALSE,certified_optimal=FALSE,
    cost_bracket_covered=NA, mass_bound_covered=NA,all_node_intervals_covered=NA,
    node_intervals=0L,all_looks_covered=NA,looks=0L,cost_ratio=NA_real_,
    omitted_bound=NA_real_,estimated_mass=NA_real_,preparation_seconds=NA_real_,
    sampling_interval_seconds=NA_real_,selection_seconds=NA_real_,cache_setup_seconds=NA_real_,
    frontier_nodes=NA_integer_,terminal_frontier=NA_integer_,unique_models=NA_integer_,
    draws=NA_integer_,mean_acceptance=NA_real_,empirical_selection=NA_character_,
    error=NA_character_,audit_status=NA_character_)
}
