# Read and compare existing RDS records only; no model fitting or search.
args<-commandArgs(trailingOnly=FALSE)
script<-sub('^--file=','',args[grep('^--file=',args)])
root<-normalizePath(file.path(dirname(script),'..'),winslash='/')
folder<-file.path(root,'results/tecator')
rows<-read.csv(file.path(folder,'outcomes.csv'))
for(i in seq_len(nrow(rows))) {
 r<-rows[i,];name<-paste0(r$dataset,'_eps_',r$epsilon)
 x<-readRDS(file.path(folder,'records',paste0(name,'.rds')))
 a<-readRDS(file.path(folder,'records',paste0(name,'_audit.rds')))
 stopifnot(x$condition$algorithm_seed==13001001L+i,
   x$condition$budget==50000,x$condition$seconds==180,
   x$mass_conservation==FALSE,x$optimizer=='newton_span',
   x$result$total_score_calls==r$calls,
   abs(x$result$snapshot$state$cost_lower-r$lower)<1e-10,
   abs(x$result$snapshot$state$cost_upper-r$cost)<1e-10,
   abs(x$result$snapshot$state$omitted_mass_bound-r$omitted_bound)<1e-12,
   abs(a$audit$summary$retained_mass-r$audit_mass)<1e-12)
}
x<-readRDS(file.path(folder,'records/tecator_sigma_0.5_eps_0.1.rds'))
ref<-readRDS(file.path(folder,'records/tecator_sigma_0.5_reference.rds'))
nodes<-x$result$nodes[x$result$nodes$selected,,drop=FALSE]
states<-t(vapply(strsplit(nodes$states,'/',fixed=TRUE),function(s)match(s,c('zero','one','multi'))-1L,integer(4)))
patterns<-read.csv(file.path(folder,'real_patterns.csv'))
key<-function(x)apply(x,1,paste,collapse='/')
stopifnot(setequal(key(states),key(patterns[,1:4])))
ids<-apply(states,1,function(s)which(rowSums(sweep(ref$states,2,s,'!='))==0))
mass<-sum(ref$posterior_mass[ids])
stopifnot(abs(mass-sum(patterns$posterior_mass))<1e-12)
pip<-read.csv(file.path(folder,'real_channel_queries.csv'))
stopifnot(max(abs(pip$full-colSums(ref$leaf_pip_contribution)))<1e-12,
 max(abs(pip$retained-colSums(ref$leaf_pip_contribution[ids,,drop=FALSE])/mass))<1e-12)
cat('PASS: six Tecator search/audit records and the current primary union match the saved tables and channel probabilities.\n')
