index <- as.integer(commandArgs(TRUE)[1L])
setwd('core');source('production.R');source('../competitors/mc3.R')
stopifnot(file.exists('../scaling_freeze.json'))
jobs <- read.csv('../scaling_config/jobs.csv',stringsAsFactors=FALSE)
job <- jobs[index,];out <- file.path('../scaling_results',job$dataset)
path <- file.path(out,paste0(job$job,'.rds'))
if(!file.exists(path)) {
  d <- readRDS(file.path(out,'data.rds'))
  cfg$mass_conservation <- TRUE
  if(job$method=='mc3') {
    started <- proc.time()[['elapsed']];score <- new_score(d)
    setup_seconds <- proc.time()[['elapsed']]-started
    ans <- mc3_search(score$score,ncol(d$X_train),d$groups,job$budget,
      max(.001,job$seconds-setup_seconds),job$seed,
      epsilon=job$epsilon,inclusion_probability=cfg$inclusion_probability)
    ans$total_seconds <- proc.time()[['elapsed']]-started
    ans$cache_setup_seconds <- setup_seconds
    ans$score_diagnostics <- score$diagnostics()
    stopifnot(ans$calls==ans$score_diagnostics[['calls']],ans$calls<=job$budget)
  } else {
    cfg$envelope_optimizer <- if(job$method=='ars_bfgs')'bfgs' else 'newton_span'
    condition <- list(method='adaptive_kl',epsilon=job$epsilon,budget=job$budget,
      seconds=job$seconds,algorithm_seed=job$seed,tree_order='forward')
    ans <- production_search(d,condition)
  }
  atomic_rds(list(job=job,result=ans),path)
} else ans <- readRDS(path)$result
write_json(list(job=job$job,status=ans$status,seconds=ans$total_seconds,
  completed=as.character(Sys.time()),process=Sys.getpid()),file.path(out,paste0(job$job,'.json')))
cat(job$job,ans$status,round(ans$total_seconds,2),'seconds\n')
if(ans$status %in% c('implementation-error','evidence_interval_conflict'))quit(status=9L)
