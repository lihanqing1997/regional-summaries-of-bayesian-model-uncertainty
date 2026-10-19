# Isolated study harness. The released manuscript and software are never edited.
source(file.path("vendor", "region_restricted_smc.R"))
source(file.path("vendor", "exact_gaussian_helpers.R"))
source(file.path("vendor", "published_bma_comparators.R"))
source(file.path("vendor", "exact_bma.R"))
cfg <- jsonlite::fromJSON(file.path("config", "study.json"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

args_named <- function() {
  x <- commandArgs(TRUE)
  if (any(!grepl("^--[^=]+=", x))) stop("Use --name=value arguments")
  setNames(lapply(x, function(s) sub("^--[^=]+=", "", s)),
           sub("^--([^=]+)=.*", "\\1", x))
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  tmp <- paste0(path, ".writing")
  saveRDS(x, tmp, compress="gzip")
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) stop("Cannot commit checkpoint: ", path)
}
write_json <- function(x, path) {
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  jsonlite::write_json(x, path, auto_unbox=TRUE, pretty=TRUE, na="null")
}
limited <- function(fun, seconds) {
  setTimeLimit(cpu=Inf, elapsed=seconds, transient=TRUE)
  on.exit(setTimeLimit(cpu=Inf, elapsed=Inf, transient=FALSE))
  fun()
}
error_status <- function(e) {
  m <- conditionMessage(e)
  if (grepl("time limit", m)) "time-limit" else
  if (grepl("RESOURCE_LIMIT|allocate vector|cannot allocate", m)) "resource-limit" else
  if (grepl("NUMERICAL_LIMIT", m)) "numerical-limit" else
  "implementation-error"
}

# Exact optimization with an upper cost cap and dominance within equal-cost
# classes. An optimum no more costly than a verified greedy feasible set can
# use at most floor(cap/c) items of cost c, so only that many largest benefits
# of each cost are needed. No posterior mass is discarded from ratio checks.
original_cost_solver <- rrs_minimum_cost_log_benefit
capped_reference_solver <- function(log_benefit, costs, eligible,
                                       log_threshold, feasibility_check=NULL) {
  ids <- which(eligible)
  no <- list(feasible=FALSE, selected=integer(), budget=Inf)
  if (!length(ids)) return(no)
  if (is.infinite(log_threshold) && log_threshold < 0) {
    j <- ids[which.min(costs[ids])]
    return(list(feasible=TRUE,selected=j,budget=costs[j],log_benefit=log_benefit[j]))
  }
  check <- function(j) is.null(feasibility_check) || isTRUE(feasibility_check(j))
  if (rrs_log_sum_exp(log_benefit[ids]) < log_threshold || !check(ids)) return(no)
  greedy <- function(ordering) {
    vals <- log_benefit[ordering]
    shift <- max(vals)
    cumulative <- cumsum(exp(vals-shift))
    hit <- which(log(cumulative)+shift >= log_threshold)
    if (!length(hit)) return(ordering)
    k <- hit[1L]
    while (k < length(ordering) && !check(ordering[seq_len(k)])) k <- k+1L
    ordering[seq_len(k)]
  }
  a <- greedy(ids[order(-log_benefit[ids],ids)])
  b <- greedy(ids[order(-(log_benefit[ids]-log(costs[ids])),ids)])
  upper <- if (sum(costs[a]) <= sum(costs[b])) a else b
  cap <- as.integer(sum(costs[upper]))
  kept <- unlist(lapply(split(ids,costs[ids]),function(js) {
    js <- js[order(-log_benefit[js],js)]
    head(js,cap %/% costs[js[1L]])
  }),use.names=FALSE)
  kept <- sort(kept)
  bytes <- 4*length(kept)*(as.double(cap)+1)+64*(as.double(cap)+1)
  if (bytes > cfg$dp_matrix_mb*1024^2)
    stop(sprintf("RESOURCE_LIMIT: exact selection requires %.1f MB",bytes/1024^2))
  best <- rep(-Inf,cap+1L)
  take <- matrix(FALSE,length(kept),cap+1L)
  for (pos in seq_along(kept)) {
    j <- kept[pos]; cst <- as.integer(costs[j])
    prev <- best; src <- seq_len(cap+1L-cst); dst <- src+cst
    left <- prev[src]; right <- log_benefit[j]
    cand <- pmax(left,right)+log1p(exp(-abs(left-right)))
    cand[is.infinite(left) & left<0] <- right
    improve <- cand > prev[dst]
    improve[is.na(improve)] <- FALSE
    best[dst[improve]] <- cand[improve]
    take[pos,dst[improve]] <- TRUE
  }
  # Summing the same masses in different orders can differ at roundoff scale.
  # A relaxed shortlist is safe: every reconstructed set must still pass the
  # original omitted-mass ratio check before it can be returned.
  rounding <- 64*.Machine$double.eps*max(1,abs(log_threshold))
  for (bi in which(best >= log_threshold-rounding)) {
    k <- bi; selected <- integer()
    for (pos in rev(seq_along(kept))) if (take[pos,k]) {
      j <- kept[pos]; selected <- c(j,selected); k <- k-costs[j]
    }
    if (length(selected) && check(selected))
      return(list(feasible=TRUE,selected=selected,budget=sum(costs[selected]),
                  log_benefit=rrs_log_sum_exp(log_benefit[selected])))
  }
  # A verified feasible incumbent exists; returning it as the exact optimum
  # would be incorrect, so fail explicitly if numerical backtracking disagrees.
  stop("NUMERICAL_LIMIT: exact reference DP could not verify the reconstructed optimum")
}

# Preserve the production algorithm for the diagnosis. Guard the ORIGINAL
# solver's allocation rather than changing its decisions or tie-breaking.
rrs_minimum_cost_log_benefit <- function(log_benefit,costs,eligible,log_threshold,
                                       feasibility_check=NULL) {
  m<-sum(eligible);cap<-sum(costs[eligible])
  bytes<-4*m*(as.double(cap)+1)+64*(as.double(cap)+1)
  if(bytes>cfg$dp_matrix_mb*1024^2)
    stop(sprintf("RESOURCE_LIMIT: original search DP requires %.1f MB",bytes/1024^2))
  original_cost_solver(log_benefit,costs,eligible,log_threshold,feasibility_check)
}
exact_reference_optimum <- function(ref,epsilon) {
  previous<-rrs_minimum_cost_log_benefit
  assign("rrs_minimum_cost_log_benefit",capped_reference_solver,envir=.GlobalEnv)
  on.exit(assign("rrs_minimum_cost_log_benefit",previous,envir=.GlobalEnv))
  select_region_union_epsilon(log(ref$mass),ref$costs,epsilon)
}

groups_from_sizes <- function(sizes) split(seq_len(sum(sizes)),rep(seq_along(sizes),sizes))
case_registry <- function() {
  out <- list()
  add <- function(id,p,sizes,design,b,active,signal="single",orders="forward",family="main") {
    out[[length(out)+1L]] <<- data.frame(case_id=id,p=p,n=cfg$n,
      sizes=paste(sizes,collapse="-"),design=design,between=b,
      active=paste(active,collapse="-"),signal=signal,orders=orders,family=family)
  }
  for (p in cfg$dimensions) {
    add(paste0("main_p",p,"_orthogonal"),p,rep(2,p/2),"orthogonal",0,c(1,3,5))
    for (b in cfg$between_correlations) {
      label <- sprintf("b%03d",round(100*b))
      ord <- if (p==16 && b>0) "forward,reverse,shuffled" else "forward"
      add(paste0("main_p",p,"_",label),p,rep(2,p/2),"random",b,c(1,3,5),orders=ord)
    }
  }
  for (b in c(.05,.20)) {
    # Matched pair control is additional: active groups 1,2,3 in all four designs.
    for (kind in c("pairs","triples","quartets","mixed")) {
      sz <- switch(kind,pairs=rep(2,6),triples=rep(3,4),quartets=rep(4,3),mixed=c(2,3,3,4))
      add(sprintf("groups_%s_b%03d",kind,round(100*b)),12,sz,"random",b,c(1,2,3),family="group-size")
    }
    add(sprintf("double_p16_b%03d",round(100*b)),16,rep(2,8),"random",b,c(1,3,5),
        signal="double",orders="forward,reverse,shuffled",family="two-active")
  }
  do.call(rbind,out)
}
make_jobs <- function() {
  cases <- case_registry()
  jobs <- list(); k <- 0L
  for (r in seq_len(cfg$diagnostic_replications)) {
    k<-k+1L
    jobs[[k]] <- data.frame(job_id=sprintf("diagnostic_r%03d",r),phase="diagnostic",
      case_id="historical_perturbed",replication=r,case_index=0L,data_seed=202802000+100*r)
  }
  for (phase in c("pilot","production")) {
    reps <- if (phase=="pilot") cfg$pilot_replications else cfg$production_replications
    base <- if (phase=="pilot") 410000000 else 510000000
    for (r in seq_len(reps)) for (i in seq_len(nrow(cases))) {
      k<-k+1L
      jobs[[k]] <- data.frame(job_id=sprintf("%s_%s_r%03d",phase,cases$case_id[i],r),
        phase=phase,case_id=cases$case_id[i],replication=r,case_index=i,
        data_seed=base+100000*i+100*r)
    }
  }
  do.call(rbind,jobs)
}
search_registry <- function(jobs=make_jobs()) {
  cases<-case_registry();out<-list();k<-0L
  for(i in seq_len(nrow(jobs))) {
    job<-jobs[i,];diag<-job$phase=="diagnostic"
    ep<-if(diag).1 else cfg$epsilon
    orders<-if(diag)"forward" else strsplit(cases$orders[job$case_index],",",fixed=TRUE)[[1L]]
    budgets<-if(diag)cfg$diagnostic_budgets else cfg$score_budget
    seeds<-if(diag)seq_len(cfg$diagnostic_algorithm_seeds) else 1L
    plan<-expand.grid(epsilon=ep,tree_order=orders,budget=budgets,seed_index=seeds,
                     mode="analytical",stringsAsFactors=FALSE)
    if(diag) {
      oracle<-plan;oracle$mode<-"oracle_range"
      exact<-plan[1,,drop=FALSE];exact$mode<-"exact"
      plan<-rbind(plan,oracle,exact)
    }
    plan$search_name<-sprintf("search_%s_%d_s%d_e%03d_%s",plan$mode,plan$budget,
                             plan$seed_index,round(100*plan$epsilon),plan$tree_order)
    plan$algorithm_seed<-job$data_seed+10000+plan$seed_index*100+round(100*plan$epsilon)
    k<-k+1L;out[[k]]<-cbind(job[rep(1,nrow(plan)),],plan)
  }
  do.call(rbind,out)
}
historical_data <- function(seed) {
  d <- cb_block_orthogonal_case(p=12,n_train=30,n_test=500,rho=.95,
                               signal=1.4,sigma=1,seed=seed,heterogeneous=FALSE)
  set.seed(seed+1L)
  loadings <- rep(c(-1,1),6)
  f <- rnorm(30); ft <- rnorm(500)
  d$X_train <- scale(sqrt(.999)*d$X_train+sqrt(.001)*tcrossprod(f,loadings))
  d$X_test <- scale(sqrt(.999)*d$X_test+sqrt(.001)*tcrossprod(ft,loadings))
  d$y_train <- as.numeric(d$X_train%*%d$beta+rnorm(30))
  d$y_test <- as.numeric(d$X_test%*%d$beta+rnorm(500))
  ym <- mean(d$y_train); d$y_train<-d$y_train-ym; d$y_test<-d$y_test-ym
  d$depth<-5L; d$population_between<-NA_real_; d$signal_variance<-NA_real_
  d
}
generate_data <- function(case,seed) {
  sizes <- as.integer(strsplit(case$sizes,"-",fixed=TRUE)[[1L]])
  groups <- groups_from_sizes(sizes); p<-sum(sizes); n<-case$n; h<-length(groups)
  rho <- cfg$within_correlation; b<-case$between
  sigma_x <- matrix(b,p,p)
  for (g in groups) sigma_x[g,g] <- rho
  diag(sigma_x)<-1
  set.seed(seed)
  if (case$design=="orthogonal") {
    # Centered orthonormal basis gives the exact empirical correlation target.
    Q<-qr.Q(qr(scale(matrix(rnorm(n*p),n,p),center=TRUE,scale=FALSE)))
    X<-sqrt(n-1)*Q%*%chol(sigma_x)
  } else {
    global<-rnorm(n); latent<-matrix(rnorm(n*h),n,h); noise<-matrix(rnorm(n*p),n,p)
    X<-sqrt(b)*global+sqrt(1-rho)*noise
    for (j in seq_along(groups)) X[,groups[[j]]]<-X[,groups[[j]],drop=FALSE]+sqrt(rho-b)*latent[,j]
    X<-scale(X,center=TRUE,scale=TRUE)
  }
  active<-as.integer(strsplit(case$active,"-",fixed=TRUE)[[1L]])
  beta<-numeric(p); beta[vapply(groups[active],`[`,integer(1),1L)]<-c(1.4,-1.12,.98)
  baseline_beta<-beta
  if (case$signal=="double") {
    # Split the first positive effect over two positively correlated variables,
    # then match population Var(X beta) to the paired single-signal control.
    beta[groups[[active[1L]]][1:2]]<-c(.7,.7)
    beta<-beta*sqrt(as.numeric(crossprod(baseline_beta,sigma_x%*%baseline_beta))/
                    as.numeric(crossprod(beta,sigma_x%*%beta)))
  }
  y<-as.numeric(X%*%beta+rnorm(n,sd=cfg$sigma)); y<-y-mean(y)
  list(X_train=X,y_train=y,beta=beta,groups=groups,sigma=cfg$sigma,depth=h,
       population_between=b,signal_variance=as.numeric(crossprod(beta,sigma_x%*%beta)))
}
data_diagnostics <- function(d) {
  R<-cor(d$X_train); same<-outer(rep(seq_along(d$groups),lengths(d$groups)),
                              rep(seq_along(d$groups),lengths(d$groups)),"==")
  upper<-upper.tri(R)
  list(n=nrow(d$X_train),p=ncol(d$X_train),within_mean=mean(R[upper & same]),
       between_mean=mean(R[upper & !same]),between_max_abs=max(abs(R[upper & !same])),
       group_sizes=lengths(d$groups),beta=d$beta,signal_variance=d$signal_variance)
}
new_score <- function(d,memoize=TRUE) make_known_sigma_support_score(
  d$X_train,d$y_train,sigma=d$sigma,inclusion_probability=cfg$inclusion_probability,
  tau2=cfg$tau2,memoize=memoize)

# Enumeration is audit-only. Practical searches construct their own score and
# cannot read this object, even when a resumed audit already exists on disk.
enumerate_truth <- function(d) {
  supports<-enumerate_supports(ncol(d$X_train)); score<-new_score(d,FALSE)
  logs<-rrs_evaluate_log_score(score$score,supports)
  logZ<-cb_log_sum_exp(logs); posterior<-exp(logs-logZ)
  list(supports=supports,logs=logs,logZ=logZ,posterior=posterior,
       pip=as.numeric(crossprod(posterior,supports)))
}
make_reference <- function(d,truth) {
  h<-length(d$groups)
  counts<-vapply(d$groups,function(g) pmin(2L,rowSums(truth$supports[,g,drop=FALSE])),
                 numeric(nrow(truth$supports)))
  if (h==1) counts<-matrix(counts,ncol=1)
  processed<-seq_len(d$depth)
  codes<-as.integer(counts[,processed,drop=FALSE]%*%3^(seq_along(processed)-1))
  leaf_states<-as.matrix(expand.grid(rep(list(0:2),length(processed))))
  mass<-numeric(nrow(leaf_states)); agg<-rowsum(truth$posterior,codes,reorder=FALSE)
  mass[as.integer(rownames(agg))+1L]<-agg[,1L]
  costs<-1+rowSums(1L+(leaf_states>0))
  list(counts=counts,codes=codes,states=leaf_states,mass=mass,costs=costs)
}
region_membership <- function(states,ref) {
  selected<-rep(TRUE,nrow(ref$counts))
  for (j in which(states!="free")) selected<-selected &
    ref$counts[,j]==match(states[j],c("zero","one","multi"))-1L
  selected
}
make_oracle <- function(d,truth,proposal,ref) {
  cache<-new.env(hash=TRUE,parent=emptyenv())
  get<-function(states,depth=NULL) {
    key<-paste(states,collapse="/")
    if (exists(key,cache,inherits=FALSE)) return(base::get(key,cache))
    inside<-region_membership(states,ref)
    region<-proposal$region_factory(states,1,key)
    lg<-truth$logs[inside]-region$log_prob(truth$supports[inside,,drop=FALSE])
    z<-cb_log_sum_exp(truth$logs[inside])
    ans<-list(rho=max(lg)-min(lg),lower=min(lg),upper=max(lg),z=z)
    assign(key,ans,cache); ans
  }
  list(get=get,rho=function(states,depth=NULL)get(states)$rho,
       log_g=function(states,region=NULL){x<-get(states);c(x$lower,x$upper)},
       log_z=function(states,region=NULL)rep(get(states)$z,2))
}
node_table <- function(search) {
  do.call(rbind,lapply(search$nodes,function(nd) {
    pf<-nd$evidence_preflight
    data.frame(id=nd$id,states=paste(nd$region$states,collapse="/"),
      rho=nd$rho,log_lower=nd$log_lower,log_upper=nd$log_upper,
      cost=nd$region$cost,selected=nd$id%in%search$selected_ids,
      frontier=nd$id%in%search$frontier_ids,evaluated=nd$evaluated,
      blocked=isTRUE(nd$preflight_blocked),required_paths=pf$required_paths%||%NA_real_,
      scheduled_paths=pf$scheduled_paths%||%NA_real_,
      evidence_method=nd$evidence_fit$method%||%"none",stringsAsFactors=FALSE)
  }))
}
run_search <- function(d,epsilon,mode,budget,seed,order="forward",truth=NULL,ref=NULL) {
  started<-proc.time()[["elapsed"]]
  result<-list(epsilon=epsilon,mode=mode,budget=budget,algorithm_seed=seed,order=order,
    status="implementation-error",error=NA_character_,search=NULL,nodes=NULL)
  score<-NULL; prep<-NA_real_
  ans<-tryCatch(limited(function() {
    t<-proc.time()[["elapsed"]]; score<<-new_score(d)
    proposal<-make_screened_group_region_factory(score$score,ncol(d$X_train),d$groups)
    if (mode=="analytical") {
      bound<-cb_make_practical_group_log_g_range_bound(d$X_train,d$y_train,
        sigma=d$sigma,tau2=cfg$tau2,inclusion_probability=cfg$inclusion_probability,
        groups=d$groups,configuration_tables=proposal$configuration_tables)
      rho<-bound$rho; log_g<-NULL; log_z<-NULL
    } else {
      oracle<-make_oracle(d,truth,proposal,ref)
      rho<-oracle$rho; log_g<-oracle$log_g
      log_z<-if(mode=="exact") oracle$log_z else NULL
    }
    ord<-seq_along(d$groups)
    if (order=="reverse") ord<-rev(ord)
    if (order=="shuffled") {set.seed(918271);ord<-sample(ord)}
    tree<-new_group_state_tree(rep(cfg$inclusion_probability,ncol(d$X_train)),d$groups,
      split_order=ord,max_depth=d$depth,rho=rho,log_g_bounds=log_g,
      log_evidence_bounds=log_z,region_factory=proposal$region_factory,
      target_mutation=make_group_gibbs_target_mutation,
      state_cost=c(zero=1L,one=2L,multi=2L),base_cost=1L)
    prep<<-proc.time()[["elapsed"]]-t
    prep_calls<-unname(score$diagnostics()[["calls"]])
    score$reset(clear_cache=FALSE);set.seed(seed)
    remaining<-budget-prep_calls
    if(remaining<1) stop("Preparation consumes the evaluation budget")
    st<-proc.time()[["elapsed"]]
    search<-if(mode=="exact") select_minimum_cost_region_tree(tree,score$score,
      evidence_evaluator=function(...)stop("Exact node unexpectedly sampled"),
      epsilon=epsilon,cost_tolerance=cfg$cost_tolerance,
      max_score_evaluations=remaining,max_actions=cfg$action_budget) else
      fit_certified_region_selection(tree,score$score,epsilon,
        cost_tolerance=cfg$cost_tolerance,delta=cfg$delta,n_pilot=60L,
        n_paths=cfg$minimum_paths,maximum_relative_radius=cfg$maximum_relative_radius,
        maximum_paths=cfg$maximum_paths,evidence_method="iid",
        max_score_evaluations=remaining,max_actions=cfg$action_budget)
    result$status<-search$registered_status
    result$preparation_calls<-prep_calls
    result$search_seconds<-proc.time()[["elapsed"]]-st
    result$nodes<-node_table(search)
    result$trace<-search$trace
    result$search<-search[c("cost_lower","cost_upper","cost_gap","omitted_mass_bound",
      "epsilon_feasible","bounds_met","total_score_evaluations","bound_score_evaluations",
      "evidence_score_evaluations","budget_overshoot","actions","expansions")]
    result
  },cfg$search_seconds_limit),error=function(e) {
    result$status<-error_status(e);result$error<-conditionMessage(e);result
  })
  ans$preparation_seconds<-prep
  ans$total_seconds<-proc.time()[["elapsed"]]-started
  ans$score_diagnostics<-if(is.null(score)) NULL else score$diagnostics()
  ans
}

evaluate_mask <- function(mask,truth,epsilon) {
  mass<-sum(truth$posterior[mask]); out<-list(retained_mass=mass,omitted_mass=1-mass,
    target_met=mass>=1-epsilon-1e-10,model_count=sum(mask),maximum_pip_error=NA_real_)
  if(mass>0) out$maximum_pip_error<-max(abs(as.numeric(crossprod(
    truth$posterior[mask]/mass,truth$supports[mask,,drop=FALSE]))-truth$pip))
  out
}
audit_search <- function(ans,truth,ref,d=NULL) {
  if(is.null(ans$nodes)) return(list(audit_available=FALSE))
  nodes<-ans$nodes; mask<-rep(FALSE,nrow(truth$supports));covered<-logical(nrow(nodes))
  node_diagnostics<-NULL
  if(!is.null(d)) {
    sc<-new_score(d)
    prop<-make_screened_group_region_factory(sc$score,ncol(d$X_train),d$groups)
    root<-prop$region_factory(rep("free",length(d$groups)),1,"audit-root")
    logq<-root$log_prob(truth$supports)
    logg<-truth$logs-logq
    node_diagnostics<-nodes
    node_diagnostics$exact_rho<-NA_real_
    node_diagnostics$importance_cv_squared<-NA_real_
  }
  for(i in seq_len(nrow(nodes))) {
    inside<-region_membership(strsplit(nodes$states[i],"/",fixed=TRUE)[[1L]],ref)
    z<-cb_log_sum_exp(truth$logs[inside]);tol<-1e-9*max(1,abs(z))
    covered[i]<-nodes$log_lower[i]<=z+tol && z<=nodes$log_upper[i]+tol
    if(nodes$selected[i]) mask<-mask|inside
    if(!is.null(d)) {
      node_diagnostics$exact_rho[i]<-diff(range(logg[inside]))
      conditional_logq<-logq[inside]-cb_log_sum_exp(logq[inside])
      log_second_moment<-cb_log_sum_exp(2*truth$logs[inside]-conditional_logq)-2*z
      node_diagnostics$importance_cv_squared[i]<-expm1(max(0,log_second_moment))
    }
  }
  ev<-evaluate_mask(mask,truth,ans$epsilon)
  c(list(audit_available=TRUE,interval_coverage=all(covered),coverage_tolerance=1e-9,
    selected_regions=sum(nodes$selected),selected_cost=if(any(nodes$selected))
      sum(nodes$cost[nodes$selected]) else Inf,
    omitted_bound_covered=ev$omitted_mass<=ans$search$omitted_mass_bound+1e-10,
    node_diagnostics=node_diagnostics),ev)
}
reference_comparators <- function(d,truth,ref,epsilon) {
  results<-list()
  one<-function(name,fun) {
    start<-proc.time()[["elapsed"]]
    r<-tryCatch(limited(fun,cfg$reference_seconds_limit),error=function(e)
      list(status=error_status(e),error=conditionMessage(e)))
    r$method<-name;r$seconds<-proc.time()[["elapsed"]]-start
    r$status<-r$status%||%"evaluated";results[[name]]<<-r
  }
  one("exact_optimal_regions",function() {
    z<-exact_reference_optimum(ref,epsilon)
    c(list(cost=z$budget,entries=length(z$selected)),
      evaluate_mask((ref$codes+1L) %in% z$selected,truth,epsilon))
  })
  for(kind in c("ranked_regions","mass_per_cost_regions")) local({
    k<-kind
    one(k,function() {
      ord<-if(k=="ranked_regions") order(-ref$mass) else order(-ref$mass/ref$costs)
      count<-which(cumsum(ref$mass[ord])>=1-epsilon)[1L]
      ids<-ord[seq_len(count)]
      c(list(cost=sum(ref$costs[ids]),entries=length(ids)),
        evaluate_mask((ref$codes+1L) %in% ids,truth,epsilon))
    })
  })
  one("ranked_models",function() {
    x<-pbmc_hpp_credible_set(truth$supports,truth$posterior,epsilon)
    c(list(entries=sum(x$selected)),evaluate_mask(x$selected,truth,epsilon))
  })
  one("cartesian_supplied_groups_unscreened",function() {
    x<-pbmc_cartesian_set_for_partition(truth$supports,truth$posterior,d$groups,1-epsilon)
    c(list(entries=x$description_entries,factorized_mass=x$factorized_mass),
      evaluate_mask(x$selected,truth,epsilon))
  })
  one("cartesian_learned_partition_unscreened",function() {
    x<-pbmc_cartesian_credible_set(truth$supports,truth$posterior,epsilon,pip_threshold=0,M=2)
    c(list(entries=x$description_entries,factorized_mass=x$factorized_mass),
      evaluate_mask(x$selected,truth,epsilon))
  })
  results
}

run_job <- function(job) {
  root<-file.path("results",job$job_id);dir.create(root,recursive=TRUE,showWarnings=FALSE)
  event<-function(stage,detail="") {
    write_json(list(job_id=job$job_id,stage=stage,detail=detail,time=as.character(Sys.time())),
      file.path("status",paste0(job$job_id,".json")))
    message(format(Sys.time())," ",job$job_id," ",stage," ",detail)
  }
  checkpoint<-function(name,fun) {
    path<-file.path(root,paste0(name,".rds"))
    if(file.exists(path)) return(readRDS(path))
    event(name)
    out<-fun();atomic_rds(out,path);out
  }
  diagnostic<-job$phase=="diagnostic"
  case<-if(!diagnostic) case_registry()[job$case_index,,drop=FALSE] else NULL
  d<-checkpoint("data",function() if(diagnostic) historical_data(job$data_seed) else
    generate_data(case,job$data_seed))
  write_json(data_diagnostics(d),file.path(root,"data_description.json"))
  eps<-if(diagnostic) .1 else cfg$epsilon
  orders<-if(diagnostic) "forward" else strsplit(case$orders,",",fixed=TRUE)[[1L]]
  runs<-list()
  do_run<-function(mode,budget,seed_index,ep,ord,truth=NULL,ref=NULL) {
    name<-sprintf("search_%s_%d_s%d_e%03d_%s",mode,budget,seed_index,round(100*ep),ord)
    # Common random seed across modes/budgets; separate seeds across replicates.
    sd<-as.integer(job$data_seed+10000+seed_index*100+round(100*ep))
    out<-checkpoint(name,function() run_search(d,ep,mode,budget,sd,ord,truth,ref))
    runs[[name]]<<-out
  }
  budgets<-if(diagnostic) cfg$diagnostic_budgets else cfg$score_budget
  seeds<-if(diagnostic) seq_len(cfg$diagnostic_algorithm_seeds) else 1L
  # No exact posterior has been constructed or loaded before these searches.
  for(ep in eps) for(ord in orders) for(b in budgets) for(s in seeds)
    do_run("analytical",b,s,ep,ord)
  checkpoint("practical_complete",function()list(time=as.character(Sys.time()),
    searches=names(runs),all_returned=TRUE,exact_information_used=FALSE))
  truth<-checkpoint("truth",function() {
    t<-proc.time()[["elapsed"]]
    x<-limited(function()enumerate_truth(d),cfg$reference_seconds_limit)
    x$enumeration_seconds<-proc.time()[["elapsed"]]-t;x
  })
  ref<-make_reference(d,truth)
  if(diagnostic) {
    for(b in budgets) for(s in seeds) do_run("oracle_range",b,s,.1,"forward",truth,ref)
    do_run("exact",cfg$score_budget,1,.1,"forward",truth,ref)
  }
  for(name in names(runs)) local({
    nm<-name
    checkpoint(paste0("audit_",nm),function()limited(function()audit_search(runs[[nm]],truth,ref,
                                                                   if(diagnostic)d else NULL),
                                                    cfg$reference_seconds_limit))
  })
  for(ep in eps) checkpoint(sprintf("comparators_e%03d",round(100*ep)),
                            function()reference_comparators(d,truth,ref,ep))
  checkpoint("complete",function()list(time=as.character(Sys.time()),searches=length(runs),
    phase=job$phase,case_id=job$case_id,replication=job$replication))
  event("complete")
  invisible(TRUE)
}
