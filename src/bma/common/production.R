source("refinement.R")
source("evaluate.R")

# Expose epsilon and tree order without changing the frozen search mathematics.
parameterize_runner <- function(expr) {
  if(!is.call(expr))return(expr)
  if(identical(expr[[1]],as.name("<-")) && identical(expr[[2]],quote(ans$epsilon)))
    expr[[3]] <- quote(cfg$run_epsilon)
  nm<-names(expr)
  if(!is.null(nm)) {
    if("epsilon"%in%nm)expr[[which(nm=="epsilon")]]<-quote(cfg$run_epsilon)
    if("split_order"%in%nm)expr[[which(nm=="split_order")]]<-quote(cfg$run_order)
  }
  for(i in seq_along(expr))if(is.call(expr[[i]]))expr[[i]]<-parameterize_runner(expr[[i]])
  expr
}
body(run_one)<-parameterize_runner(body(run_one))
cfg$run_epsilon<-.1
cfg$run_order<-integer()

verify_production <- function() {
  m<-read.csv("config/frozen_manifest.csv")
  stopifnot(identical(unname(tools::md5sum(m$path)),m$md5))
}
production_data <- function(case,seed) {
  case$stress<-as.logical(case$stress)
  pilot_data(case,seed)
}
production_search <- function(d,condition) {
  cfg$run_epsilon<<-condition$epsilon
  cfg$run_order<<-seq_along(d$groups)
  if(condition$tree_order=="reverse")cfg$run_order<<-rev(cfg$run_order)
  if(condition$tree_order=="shuffled"){set.seed(918271);cfg$run_order<<-sample(cfg$run_order)}
  saved<-allocation_options
  on.exit(assign("allocation_options",saved,.GlobalEnv),add=TRUE)
  if(condition$method=="equal_kl")assign("allocation_options",function(ids,nodes,records,state,epsilon,remaining,grid=adaptive_grid) {
    options<-saved(ids,nodes,records,state,epsilon,remaining,grid)
    if(nrow(options))options$log_priority<-vapply(options$id,function(id)
      if(is.null(records[[id]]))0 else -length(records[[id]]$weights),numeric(1))
    options
  },.GlobalEnv)
  ans<-refinement_search(d,condition$algorithm_seed,condition$budget,condition$seconds,
    if(condition$method=="adaptive_h")"original" else "kl")
  ans$condition<-as.list(condition)
  ans
}

value_or <- function(x,default=NA)if(is.null(x)||!length(x))default else x
search_row <- function(condition,ans) {
  s<-ans$snapshot$state
  c(as.list(condition),list(status=ans$status,error=value_or(ans$error,NA_character_),
    seconds_used=value_or(ans$total_seconds,NA_real_),calls=value_or(ans$total_score_calls,NA_real_),
    cost_lower=value_or(s$cost_lower,NA_real_),cost_upper=value_or(s$cost_upper,NA_real_),
    cost_gap=value_or(s$cost_gap,NA_real_),omitted_mass_bound=value_or(s$omitted_mass_bound,NA_real_),
    reported_feasible=isTRUE(s$incumbent$feasible),reported_certified=identical(ans$status,"certified"),
    sampled_looks=if(is.null(ans$sampling_detail))0L else nrow(ans$sampling_detail),
    audit_status="pending",exact_cost=NA_real_,retained_mass=NA_real_,
    verified_feasible=NA,verified_optimal=NA,all_looks_covered=NA))
}

# Audit/reference construction is deliberately sourced only after practical runs.
load_production_audit <- function() {
  source("audit/audit.R",local=.GlobalEnv)
  # Keep ranked-model description counts, but do not retain a million-row table.
  f<-make_pilot_reference
  code<-body(f);found<-0L
  for(i in seq_along(code))if(is.call(code[[i]]) &&
      identical(code[[i]][[1]],as.name("<-")) &&
      identical(code[[i]][[2]],quote(ans$support_log_scores))) {
    code[[i]]<-quote({
      ranked_probability<-sort(exp(ans$support_log_scores-ans$logZ),decreasing=TRUE)
      ans$ranked_model_cumulative<-cumsum(ranked_probability)
      ans$support_log_scores<-NULL
    });found<-found+1L
  }
  stopifnot(found==1L);body(f)<-code
  assign("make_production_reference",f,.GlobalEnv)
}
production_comparators <- function(ref,epsilons) {
  do.call(rbind,lapply(epsilons,function(eps) {
    rows<-list()
    for(method in c("ranked_regions","mass_per_cost_regions")) {
      ord<-if(method=="ranked_regions")order(-ref$posterior_mass) else order(-ref$posterior_mass/ref$costs)
      k<-which(cumsum(ref$posterior_mass[ord])>=1-eps)[1L];ids<-ord[seq_len(k)]
      rows[[method]]<-data.frame(epsilon=eps,method=method,cost=sum(ref$costs[ids]),
        entries=k,retained_mass=sum(ref$posterior_mass[ids]))
    }
    k<-which(ref$ranked_model_cumulative>=1-eps)[1L]
    rows$models<-data.frame(epsilon=eps,method="ranked_models",cost=NA_real_,entries=k,
      retained_mass=ref$ranked_model_cumulative[k])
    do.call(rbind,rows)
  }))
}
