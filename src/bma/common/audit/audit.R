Rcpp::sourceCpp("audit/exact_reference.cpp",cacheDir="build",rebuild=FALSE,showOutput=FALSE)
make_pilot_reference<-function(d) {
  groups<-integer(ncol(d$X_train))
  for(j in seq_along(d$groups))groups[d$groups[[j]]]<-j-1L
  ans<-pilot_exact_reference(crossprod(d$X_train),as.numeric(crossprod(d$X_train,d$y_train)),
    sum(d$y_train^2),nrow(d$X_train),d$sigma,cfg$tau2,cfg$inclusion_probability,
    groups,seconds=pilot$audit_seconds)
  # Independent spot checks against the frozen R likelihood implementation.
  set.seed(98314);masks<-unique(c(0,length(ans$support_log_scores)-1,
    sample.int(length(ans$support_log_scores),min(500,length(ans$support_log_scores)))-1L))
  supports<-t(vapply(masks,function(i)as.integer(intToBits(i))[seq_len(ncol(d$X_train))],integer(ncol(d$X_train))))
  score<-new_score(d,FALSE)
  discrepancy<-max(abs(rrs_evaluate_log_score(score$score,supports)-ans$support_log_scores[masks+1L]))
  if(discrepancy>1e-8)stop("Native audit scores disagree with the frozen R model")
  ans$score_check_max_error<-discrepancy;ans$support_log_scores<-NULL
  ans$states<-as.matrix(expand.grid(rep(list(0:2),length(d$groups)),KEEP.OUT.ATTRS=FALSE))
  ans$costs<-1L+rowSums(1L+(ans$states>0))
  ans$posterior_mass<-exp(ans$leaf_log_mass-ans$logZ)
  ans
}
audit_pilot_search<-function(ans,ref) {
  start<-proc.time()[["elapsed"]]
  opt<-select_region_union_epsilon(ref$leaf_log_mass-ref$logZ,ref$costs,pilot$epsilon)
  if(opt$solver_status!="optimal")stop("Exact reference cost was not resolved")
  inside<-function(states) {
    keep<-rep(TRUE,nrow(ref$states))
    for(j in which(states!="free"))keep<-keep & ref$states[,j]==match(states[j],c("zero","one","multi"))-1L
    keep
  }
  selected<-rep(FALSE,nrow(ref$states));nodes<-ans$nodes
  checks<-list()
  if(!is.null(nodes))for(i in seq_len(nrow(nodes))) {
    mask<-inside(strsplit(nodes$states[i],"/",fixed=TRUE)[[1]])
    z<-rrs_log_sum_exp(ref$leaf_log_mass[mask]);tol<-1e-9*max(1,abs(z))
    checks[[i]]<-data.frame(node_id=nodes$id[i],exact_log_mass=z,
      lower_covered=nodes$log_lower[i]<=z+tol,upper_covered=z<=nodes$log_upper[i]+tol)
    if(nodes$selected[i])selected<-selected|mask
  }
  checks<-if(length(checks))do.call(rbind,checks) else data.frame()
  state<-ans$snapshot$state
  mass<-sum(ref$posterior_mass[selected]);global_pip<-colSums(ref$leaf_pip_contribution)
  has_certificate<-!is.null(state)&&is.finite(state$cost_upper)&&length(state$selected_ids)>0
  interval_coverage<-if(nrow(checks))all(checks$lower_covered & checks$upper_covered) else NA
  reported_bound<-state$omitted_mass_bound%||%1
  bounds_covered<-(1-mass)<=reported_bound+1e-9
  bracket_valid<-if(is.null(state))NA else state$cost_lower<=opt$budget&&state$cost_upper>=opt$budget
  verified<-has_certificate&&bounds_covered&&mass>=1-pilot$epsilon-1e-9
  if(has_certificate) {
    cert<-state$certificate
    cert_selected<-rep(FALSE,nrow(ref$states))
    for(id in cert$selected_ids)cert_selected<-cert_selected|
      inside(strsplit(nodes$states[match(id,nodes$id)],"/",fixed=TRUE)[[1]])
    stopifnot(identical(cert_selected,selected))
  }
  summary<-list(status=if(is.null(state))"no-valid-frontier" else if(isFALSE(interval_coverage)||isFALSE(bracket_valid)||
       (ans$status=="certified"&&!verified))"violation" else "passed",
    exact_optimal_cost=opt$budget,certificate_available=has_certificate,
    retained_mass=mass,actual_omitted_mass=1-mass,mass_guarantee_verified=verified,
    reported_bound_covered=bounds_covered,cost_bracket_valid=bracket_valid,
    interval_coverage=interval_coverage,selected_models=sum(ref$leaf_model_count[selected]),
    maximum_pip_error=if(mass>0)max(abs(colSums(ref$leaf_pip_contribution[selected,,drop=FALSE])/mass-global_pip)) else NA_real_,
    actual_cost_optimal=verified&&state$cost_upper==opt$budget,
    verified_optimal=ans$status=="certified"&&verified&&state$cost_upper==opt$budget&&state$cost_gap<1,
    reference_models=ref$enumerated_models,reference_enumeration_seconds=ref$seconds,
    reference_score_check_max_error=ref$score_check_max_error,
    audit_comparison_seconds=proc.time()[["elapsed"]]-start)
  list(summary=summary,node_checks=checks)
}
