source("two_stage.R")

# A finite, predeclared grid permits adaptive interleaving and stopping by a
# union bound over node/look pairs. Never inspect a partial batch as a new look.
adaptive_grid <- c(256L, 512L, 1024L, 2048L, 4096L, 8192L)

envelope_cell_count <- function(region, cap=256L) {
  sizes <- lengths(region$groups)
  options <- vapply(seq_along(sizes), function(j) switch(region$states[j],
    free=2^sizes[j], zero=1, one=sizes[j], multi=2^sizes[j]-1-sizes[j]), numeric(1))
  cells <- 1L
  for(k in options) { if(cells*k>cap) break; cells <- cells*k }
  as.integer(cells)
}

log_interval_width <- function(lower, upper) {
  if(!is.finite(upper) || lower>=upper) return(-Inf)
  if(!is.finite(lower)) return(upper)
  upper + log(-expm1(lower-upper))
}

advance_envelope <- function(record, d, node, score, remaining, alpha,
                             grid=adaptive_grid) {
  look <- if(is.null(record)) 1L else record$look+1L
  if(look>length(grid)) return(list(blocked=TRUE, calls=0L))
  old_n <- if(is.null(record)) 0L else length(record$weights)
  prep_expected <- if(is.null(record)) envelope_cell_count(node$region) else 0L
  batch <- grid[look]-old_n
  if(prep_expected+batch>remaining) return(list(blocked=TRUE, calls=0L))
  guarded <- rrs_budgeted_log_score(score, remaining)
  if(is.null(record)) record <- list(
    proposal=make_partitioned_envelope(d, guarded, node$region$states, 256L),
    weights=numeric(), look=0L)
  prep <- attr(guarded,"score_budget_usage")()
  stopifnot(prep==prep_expected)
  x <- record$proposal$region$sample(batch)
  log_weights <- rrs_evaluate_log_score(guarded,x) -
    record$proposal$region$log_prob(x) - record$proposal$logU
  if(any(!is.finite(log_weights)) || any(log_weights>1e-10))
    stop("Invalid weight or violation of the pointwise envelope")
  record$weights <- c(record$weights,exp(pmin(0,log_weights)))
  record$look <- look
  ci <- bounded_mean_interval(record$weights, alpha)
  calls <- attr(guarded,"score_budget_usage")()
  stopifnot(length(record$weights)==grid[look],calls==prep+batch,calls<=remaining)
  list(blocked=FALSE,record=record,calls=calls,prep=prep,batch=batch,ci=ci,
    log_estimate=record$proposal$logU+log(ci$mean),
    log_lower=record$proposal$logU+log(ci$lower),
    log_upper=record$proposal$logU+log(ci$upper))
}

allocation_options <- function(ids,nodes,records,state,epsilon,remaining,
                               grid=adaptive_grid) {
  rows <- lapply(ids,function(id) {
    record <- records[[id]]
    look <- if(is.null(record)) 1L else record$look+1L
    if(look>length(grid)) return(NULL)
    n <- if(is.null(record)) 0L else length(record$weights)
    prep <- if(is.null(record)) envelope_cell_count(nodes[[id]]$region) else 0L
    calls <- prep+grid[look]-n
    if(calls>remaining) return(NULL)
    width <- log_interval_width(nodes[[id]]$log_lower,nodes[[id]]$log_upper)
    if(!is.finite(width)) return(NULL)
    # For candidate S, feasibility is sum(U_outside) <= epsilon/(1-epsilon)
    # * sum(L_inside). Width outside S has coefficient 1; inside has epsilon/
    # (1-epsilon). Use the current optimistic set as a cheaper candidate.
    inside <- id %in% state$optimistic_ids
    coefficient <- if(inside) epsilon/(1-epsilon) else 1
    # A scheduling heuristic, not a claimed prediction or confidence bound.
    # Untested regions get a discovery allowance; revisits use sqrt(n) scaling.
    shrink <- if(n==0L) .5 else 1-sqrt(n/grid[look])
    data.frame(id=id,look=look,expected_calls=calls,inside_candidate=inside,
      log_priority=width+log(coefficient)+log(shrink)-log(calls))
  })
  rows <- Filter(Negate(is.null),rows)
  if(!length(rows)) return(data.frame())
  do.call(rbind,rows)
}

adaptive_search <- function(d,seed) {
  events <- new.env(parent=emptyenv()); events$rows <- list()
  previous_delta <- cfg$delta
  cfg$delta <<- previous_delta/2
  on.exit({cfg$delta <<- previous_delta
    assign("certified_tree_region_search",original_search_engine,.GlobalEnv)},add=TRUE)
  engine <- function(tree,log_score_fn,evidence_evaluator,epsilon,cost_tolerance=0,
    max_score_evaluations=Inf,max_actions=Inf,bound_node=NULL) {
    base <- original_search_engine(tree,log_score_fn,evidence_evaluator,epsilon,
      cost_tolerance,max_score_evaluations,max_actions,bound_node=bound_node)
    if(base$certified || base$status=="evidence_interval_conflict") return(base)
    nodes <- base$nodes; frontier <- base$frontier_ids
    initial <- ars_tracker$state; state <- initial; solver <- ars_new_state_solver()
    if(isTRUE(initial$incumbent$feasible)) environment(solver)$best_incumbent <- initial
    environment(solver)$best_lower <- initial$cost_lower
    ids <- frontier[vapply(nodes[frontier],function(n)
      isTRUE(n$terminal)&&isTRUE(n$preflight_blocked),logical(1))]
    records <- list(); used <- base$total_score_evaluations; additional <- 0L
    alpha <- cfg$delta/(tree$n_declared_nodes*length(adaptive_grid))
    status <- base$status
    repeat {
      if(proc.time()[["elapsed"]]>ars_tracker$deadline) {status <- "time-limit";break}
      choices <- allocation_options(ids,nodes,records,state,epsilon,max_score_evaluations-used)
      if(!nrow(choices)) {
        unfinished <- any(vapply(ids,function(id)
          is.null(records[[id]]) || records[[id]]$look<length(adaptive_grid),logical(1)))
        status <- if(unfinished) "budget_exhausted" else "frontier_exhausted"
        break
      }
      choice <- choices[which.max(choices$log_priority),];id <- choice$id;node <- nodes[[id]]
      old_cost <- state$cost_upper
      fit <- advance_envelope(records[[id]],d,node,log_score_fn,
        max_score_evaluations-used,alpha)
      stopifnot(!fit$blocked,fit$calls==choice$expected_calls)
      used <- used+fit$calls;additional <- additional+fit$calls
      records[[id]] <- fit$record
      lo <- max(node$log_lower,fit$log_lower)
      hi <- min(node$log_upper,fit$log_upper,fit$record$proposal$logU)
      if(lo>hi+1e-10) {status <- "evidence_interval_conflict";break}
      interval <- rrs_reconcile_log_interval(lo,hi)
      node$log_lower <- interval[["log_lower"]];node$log_upper <- interval[["log_upper"]]
      node$evaluated <- TRUE;node$preflight_blocked <- FALSE
      node$evidence_fit <- list(method="staged_bounded_envelope_iid",
        n_paths=length(fit$record$weights),look=fit$record$look,alpha=alpha,
        relative_radius=fit$ci$relative_radius,logU=fit$record$proposal$logU)
      nodes[[id]] <- node
      # Log the completed look even if a later state update reaches its deadline.
      events$rows[[length(events$rows)+1L]] <- data.frame(node=id,
        states=paste(node$region$states,collapse="/"),look=fit$record$look,
        cumulative_draws=length(fit$record$weights),preparation_calls=fit$prep,
        sampling_calls=fit$batch,alpha=alpha,inside_candidate=choice$inside_candidate,
        log_priority=choice$log_priority,log_lower=node$log_lower,log_upper=node$log_upper,
        log_estimate=fit$log_estimate,relative_radius=fit$ci$relative_radius,
        positive_lower=is.finite(node$log_lower),cost_before=old_cost,cost_after=NA_real_)
      state <- solver(nodes,frontier,epsilon)
      stopifnot(state$cost_upper<=initial$cost_upper,state$cost_upper<=old_cost)
      events$rows[[length(events$rows)]]$cost_after <- state$cost_upper
      if(isTRUE(state$incumbent$feasible)&&state$omitted_mass_bound<=epsilon&&
         state$cost_gap<=cost_tolerance) {status <- "certified";break}
    }
    base$nodes <- nodes;base$status <- status
    base$registered_status <- if(status=="time-limit") status else rrs_registered_status(status)
    base$certified <- base$bounds_met <- identical(status,"certified")
    base$epsilon_feasible <- isTRUE(state$incumbent$feasible)&&state$omitted_mass_bound<=epsilon
    base$cost_optimality_met <- base$epsilon_feasible&&state$cost_gap<=cost_tolerance
    for(name in c("selected_ids","cost_lower","cost_upper","cost_gap","bracket_consistent","omitted_mass_bound"))
      base[[name]] <- state[[name]]
    base$total_score_evaluations <- used
    base$evidence_score_evaluations <- base$evidence_score_evaluations+additional
    base$actions <- base$actions+length(events$rows)
    base$refinement <- list(initial_cost=initial$cost_upper,final_cost=state$cost_upper,
      additional_calls=additional,regions=length(records),looks=length(events$rows),
      grid=adaptive_grid,alpha_per_look=alpha,joint_delta=previous_delta)
    base
  }
  assign("certified_tree_region_search",engine,.GlobalEnv)
  ans <- instrument_search(d,seed)
  ans$sampling_detail <- if(length(events$rows)) do.call(rbind,events$rows) else data.frame()
  ans
}
