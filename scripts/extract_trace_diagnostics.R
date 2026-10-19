# Read saved searches only; no search, simulation, or frozen-file modification.
args <- commandArgs(TRUE)
stopifnot(length(args)==3L)
components <- normalizePath(args[1],winslash='/')
selectors <- normalizePath(args[2],winslash='/')
dest <- args[3]
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
`%or%` <- function(x,y) if(is.null(x)||!length(x)) y else x
one <- function(path,study,job,method,case) {
  x <- readRDS(path)
  tr <- x$trace
  sd <- x$sampling_detail
  eps <- x$epsilon
  first <- which(is.finite(tr$cost_upper) & tr$omitted_mass_bound <= eps)
  first_calls <- first_cost <- first_action <- NA_real_
  phase <- if(is.null(tr)) 'unavailable' else 'none'
  if(length(first)) {
    z <- tr[first[1],]
    first_calls <- z$score_evaluations + (x$preparation_calls %or% 0)
    first_cost <- z$cost_upper; first_action <- z$action; phase <- 'base'
  } else if(!is.null(tr) && nrow(sd) && any(is.finite(sd$cost_after))) {
    ix <- which(is.finite(sd$cost_after))[1]
    extra <- x$search$refinement$additional_calls
    if(length(extra)==1) {
      first_calls <- x$total_score_calls-extra+cumsum(sd$preparation_calls+sd$sampling_calls)[ix]
    } else {
      warning('No complete refinement total: ',path)
    }
    first_cost <- sd$cost_after[ix]; phase <- 'refinement'
  }
  front <- x$nodes[x$nodes$frontier,,drop=FALSE]
  terminal <- !grepl('free',front$states,fixed=TRUE)
  state <- x$snapshot$state
  data.frame(study=study,job=job,case_index=case,method=method,
    action_history_available=!is.null(tr),first_feasible_phase=phase,first_feasible_calls=first_calls,
    first_feasible_cost=first_cost,first_feasible_action=first_action,
    final_cost=state$cost_upper,final_lower=state$cost_lower,
    final_calls=x$total_score_calls,final_seconds=x$total_seconds,
    frontier_nodes=nrow(front),exposed_terminals=sum(terminal),
    sampled_checkpoints=nrow(sd),sampled_regions=length(unique(sd$node)),
    envelope_seconds=x$timing$proposal_preparation %or% 0,
    sampling_seconds=x$timing$sampling_and_intervals %or% 0,
    selection_seconds=x$timing$selection %or% x$selection_seconds,
    initial_seconds=x$timing$initial_phase %or% NA_real_,
    stringsAsFactors=FALSE)
}
out <- list(); k<-0L
for(case in c(10:12,23:24)) for(rep in 1:100) {
  job<-sprintf('c%02d_r%03d',case,rep)
  file<-paste0(job,'_adaptive_kl_e10_b25000_s1_forward.rds')
  path<-file.path(components,'results',job,file)
  k<-k+1L;out[[k]]<-one(path,'primary',job,'ars',case)
}
for(case in c(7,8,23,24)) for(rep in 1:12) for(method in c('ars','mass','mass_cost')) {
  job<-sprintf('c%02d_r%03d',case,rep)
  path<-file.path(selectors,'results',job,paste0(method,'.rds'))
  k<-k+1L;out[[k]]<-one(path,'matched',job,method,case)
}
out<-do.call(rbind,out)
write.csv(out,file.path(dest,'saved_trace_diagnostics.csv'),row.names=FALSE)
cat('Read',nrow(out),'saved searches. First-feasible phase:\n')
print(with(out,table(study,first_feasible_phase)))
