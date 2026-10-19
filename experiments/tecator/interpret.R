# Opt-in interpretation calculation from saved Tecator records; no search is run.
setwd("core");source("production.R");load_production_audit()
saved<-readRDS("../validation/results/numerical_end_to_end/tecator_sigma_0.5_eps_0.1.rds")
d<-saved$data;a<-saved$result
ref<-readRDS("../validation/results/numerical_end_to_end/tecator_sigma_0.5_reference.rds")
nodes<-a$nodes[a$nodes$selected,,drop=FALSE]
states<-t(vapply(strsplit(nodes$states,"/",fixed=TRUE),function(x)
  match(x,c("zero","one","multi"))-1L,integer(4)))
ids<-apply(states,1,function(x)which(rowSums(sweep(ref$states,2,x,"!="))==0))
mass<-sum(ref$posterior_mass[ids])
patterns<-data.frame(group_A=states[,1],group_B=states[,2],group_C=states[,3],group_D=states[,4],
  cost=ref$costs[ids],models=ref$leaf_model_count[ids],posterior_mass=ref$posterior_mass[ids])
patterns<-patterns[order(-patterns$posterior_mass),]
write.csv(patterns,"../real_patterns.csv",row.names=FALSE)
groups<-data.frame(group=LETTERS[1:4],channels=vapply(d$groups,function(g)
  paste(d$channels[g],collapse=", "),character(1)),
  full_active=colSums(ref$states*0+ (ref$states>0)*ref$posterior_mass),
  retained_active=colSums((ref$states[ids,,drop=FALSE]>0)*ref$posterior_mass[ids])/mass)
write.csv(groups,"../real_group_queries.csv",row.names=FALSE)
pip<-data.frame(channel=d$channels,full=colSums(ref$leaf_pip_contribution),
  retained=colSums(ref$leaf_pip_contribution[ids,,drop=FALSE])/mass)
write.csv(pip,"../real_channel_queries.csv",row.names=FALSE)
# Post-search enumeration is used only to illustrate predictions and cost sensitivity.
native<-pilot_exact_reference(crossprod(d$X_train),as.numeric(crossprod(d$X_train,d$y_train)),
  sum(d$y_train^2),nrow(d$X_train),d$sigma,cfg$tau2,cfg$inclusion_probability,rep(0:3,each=4))
prob<-exp(native$support_log_scores-native$logZ)
G<-crossprod(d$X_train);xy<-crossprod(d$X_train,d$y_train)
full_beta<-retained_beta<-numeric(16)
for(mask in 0:65535) {
  support<-as.integer(intToBits(mask))[1:16];s<-which(support==1L)
  if(!length(s))next
  beta<-solve(G[s,s,drop=FALSE]+diag(1/cfg$tau2,length(s)),xy[s])
  full_beta[s]<-full_beta[s]+prob[mask+1L]*beta
  counts<-vapply(d$groups,function(g)min(2L,sum(support[g])),integer(1))
  inside<-any(rowSums(sweep(states,2,counts,"!="))==0)
  if(inside)retained_beta[s]<-retained_beta[s]+prob[mask+1L]*beta/mass
}
rmse<-function(beta)sqrt(mean((as.numeric(d$X_test%*%beta)-d$y_test)^2))*d$response_sd
write_json(list(selected_mass=mass,full_rmse=rmse(full_beta),retained_rmse=rmse(retained_beta),
  bound=a$snapshot$state$omitted_mass_bound,channels=d$channels,
  within_between_correlation=range(cor(d$X_train)[upper.tri(cor(d$X_train))])),"../real_interpretation.json")
print(patterns);print(groups)
