source("adaptive.R")

# Parameterize only the old runner's literal global call cap in this isolated
# diagnostic. All search, proposal, CI and allocation functions are unchanged.
replace_budget_literal <- function(x) {
  if(is.numeric(x)&&length(x)==1L&&x==25000) return(quote(cfg$score_budget))
  if(is.call(x)) for(i in seq_along(x)) {
    value <- replace_budget_literal(x[[i]])
    # Assigning NULL through [[ would delete an argument rather than retain it.
    if(!is.null(value)) x[[i]] <- value
  }
  x
}
body(run_one) <- replace_budget_literal(body(run_one))

timed_diagnostic_search <- function(d,seed,budget=25000L,seconds=600) {
  old_budget <- cfg$score_budget;old_seconds <- cfg$search_seconds_limit
  saved <- list(prepare=make_partitioned_envelope,advance=advance_envelope,
    allocate=allocation_options,initial=original_search_engine,
    engine=certified_tree_region_search)
  timer <- new.env(parent=emptyenv())
  timer$prepare <- timer$advance <- timer$allocate <- timer$initial <- 0
  timer$initial_selection <- 0;timer$events <- list()
  cfg$score_budget <<- budget;cfg$search_seconds_limit <<- seconds
  on.exit({cfg$score_budget <<- old_budget;cfg$search_seconds_limit <<- old_seconds
    assign("make_partitioned_envelope",saved$prepare,.GlobalEnv)
    assign("advance_envelope",saved$advance,.GlobalEnv)
    assign("allocation_options",saved$allocate,.GlobalEnv)
    assign("original_search_engine",saved$initial,.GlobalEnv)
    assign("certified_tree_region_search",saved$engine,.GlobalEnv)},add=TRUE)
  assign("make_partitioned_envelope",function(...) {
    t <- proc.time()[["elapsed"]]
    on.exit(timer$prepare <- timer$prepare+proc.time()[["elapsed"]]-t)
    saved$prepare(...)
  },.GlobalEnv)
  assign("advance_envelope",function(record,d,node,score,remaining,alpha,grid=adaptive_grid) {
    t <- proc.time()[["elapsed"]];p <- timer$prepare
    result <- saved$advance(record,d,node,score,remaining,alpha,grid)
    elapsed <- proc.time()[["elapsed"]]-t;prep <- timer$prepare-p
    timer$advance <- timer$advance+elapsed
    if(!result$blocked) timer$events[[length(timer$events)+1L]] <- data.frame(
      node=node$id,look=result$record$look,preparation_seconds=prep,
      sampling_seconds=max(0,elapsed-prep),calls=result$calls)
    result
  },.GlobalEnv)
  assign("allocation_options",function(...) {
    t <- proc.time()[["elapsed"]]
    on.exit(timer$allocate <- timer$allocate+proc.time()[["elapsed"]]-t)
    saved$allocate(...)
  },.GlobalEnv)
  assign("original_search_engine",function(...) {
    t <- proc.time()[["elapsed"]];s <- ars_tracker$selection_seconds
    on.exit({timer$initial <- timer$initial+proc.time()[["elapsed"]]-t
      timer$initial_selection <- ars_tracker$selection_seconds-s})
    saved$initial(...)
  },.GlobalEnv)
  ans <- adaptive_search(d,seed)
  ans$timing <- list(total=ans$total_seconds,initial_phase=timer$initial,
    initial_selection=timer$initial_selection,
    proposal_preparation=timer$prepare,
    sampling_and_intervals=max(0,timer$advance-timer$prepare),
    selection=ans$selection_seconds,allocation=timer$allocate,
    other=ans$total_seconds-timer$advance-ans$selection_seconds-timer$allocate)
  ans$timing_events <- if(length(timer$events)) do.call(rbind,timer$events) else data.frame()
  ans$diagnostic_budget <- budget;ans$diagnostic_seconds_limit <- seconds
  ans
}
