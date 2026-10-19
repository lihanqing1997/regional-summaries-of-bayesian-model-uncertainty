index <- as.integer(commandArgs(TRUE)[1L])
setwd('core');source('production.R');source('../scaling_reference.R');source('../scaling_audit.R')
stopifnot(file.exists('../scaling_freeze.json'))
datasets <- read.csv('../scaling_config/datasets.csv',stringsAsFactors=FALSE)
dataset <- datasets[index,];out <- file.path('../scaling_results',dataset$dataset)
jobs <- read.csv('../scaling_config/jobs.csv',stringsAsFactors=FALSE)
jobs <- jobs[jobs$dataset==dataset$dataset,]
terminal <- vapply(jobs$job,function(name)file.exists(file.path(out,paste0(name,'.rds'))) ||
  file.exists(file.path(out,paste0(name,'.failure.json'))),logical(1))
stopifnot(all(terminal))
path <- file.path(out,'reference.rds')
if(!file.exists(path))atomic_rds(make_scaling_reference(readRDS(file.path(out,'data.rds'))),path)
ref <- readRDS(path)
opt <- ars_interval_knapsack(ref$leaf_log_mass-ref$logZ,ref$leaf_log_mass-ref$logZ,
  as.integer(ref$costs),rep(TRUE,length(ref$costs)),.1,TRUE)
stopifnot(opt$status=='optimal')
rows <- list()
for(i in seq_len(nrow(jobs))) {
  job <- jobs[i,];result_path <- file.path(out,paste0(job$job,'.rds'))
  methods <- if(job$method=='mc3')c('mc3_frequency','mc3_renormalized') else job$method
  for(method in methods) {
    row <- scaling_empty_row(dataset,job,method,opt$upper_cost)
    if(file.exists(result_path)) {
      ans <- readRDS(result_path)$result
      row$status <- ans$status;row$seconds <- ans$total_seconds;row$error <- ans$error %||% NA_character_
      if(job$method=='mc3') {
        row$calls <- ans$calls
        checked <- scaling_audit_mc3(ans,ref,job$epsilon,opt$upper_cost,sub('mc3_','',method))
      } else {
        row$calls <- ans$total_score_calls
        checked <- scaling_audit_ars(ans,ref,job$epsilon,opt$upper_cost)
      }
      for(name in names(checked))row[[name]] <- checked[[name]]
    } else {
      failed <- jsonlite::fromJSON(file.path(out,paste0(job$job,'.failure.json')))
      row$status <- failed$status;row$error <- failed$error;row$audit_status <- 'no-summary'
    }
    rows[[length(rows)+1L]] <- as.data.frame(row,stringsAsFactors=FALSE)
  }
}
results <- do.call(rbind,rows)
write.csv(results,file.path(out,'outcomes.csv'),row.names=FALSE,na='NA')
write_json(list(dataset=dataset$dataset,status='complete',rows=nrow(results),
  exact_cost=opt$upper_cost,reference_models=ref$enumerated_models,
  reference_seconds=ref$seconds,score_check_max_error=ref$score_check_max_error,
  interval_misses=sum(results$all_node_intervals_covered==FALSE,na.rm=TRUE),
  checkpoint_misses=sum(results$all_looks_covered==FALSE,na.rm=TRUE),
  false_certificates=sum(results$reported_certified & !results$certified_optimal,na.rm=TRUE),
  finished=as.character(Sys.time())),file.path(out,'audit_complete.json'))
cat(dataset$dataset,'audit complete; exact cost',opt$upper_cost,'\n')
