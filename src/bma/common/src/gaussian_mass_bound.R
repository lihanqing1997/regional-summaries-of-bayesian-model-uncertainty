# A proposal-independent bound on regional Gaussian evidence. For any support
# Gamma in the region, C_Gamma <= C_A where A is the union of allowed columns,
# so its quadratic fit is no greater than the fit using A. The determinant
# factorizes by groups. Replacing the earlier selected columns by all earlier
# allowed columns decreases each conditional determinant, producing an upper
# bound after summing the local prior-weighted inverse square roots.
ars_gaussian_mass_upper<-function(d,score,budget,tau2=cfg$tau2,pi=cfg$inclusion_probability) {
  X<-d$X_train;groups<-d$groups;p<-ncol(X);G<-crossprod(X)
  tables<-lapply(groups,function(g)as.matrix(expand.grid(rep(list(0:1),length(g)))))
  fit_cache<-new.env(hash=TRUE,parent=emptyenv())
  determinant_cache<-new.env(hash=TRUE,parent=emptyenv())
  chol_logdet<-function(A)2*sum(log(diag(chol(A))))
  function(states) {
    active_groups<-which(states!="zero")
    columns<-unlist(groups[active_groups],use.names=FALSE)
    key<-paste(active_groups,collapse="/")
    if(!nzchar(key))key<-"none"
    if(!exists(key,fit_cache,inherits=FALSE)) {
      if(unname(score$diagnostics()[["calls"]])>=budget)return(Inf)
      gamma<-integer(p);gamma[columns]<-1L
      value<-score$score(gamma)
      ldet<-if(length(columns))chol_logdet(diag(length(columns))+tau2*G[columns,columns,drop=FALSE]) else 0
      logprior<-length(columns)*log(pi)+(p-length(columns))*log1p(-pi)
      assign(key,value+.5*ldet-logprior,fit_cache)
    }
    fit<-get(key,fit_cache)
    local_sum<-function(order) {
      previous<-integer();total<-0
      for(j in order) {
        table<-tables[[j]];k<-rowSums(table)
        keep<-switch(states[j],free=rep(TRUE,length(k)),zero=k==0,one=k==1,multi=k>=2)
        dkey<-paste(j,paste(previous,collapse=","),sep="|")
        if(!exists(dkey,determinant_cache,inherits=FALSE)) {
          g<-groups[[j]];M<-tau2*G[g,g,drop=FALSE]
          if(length(previous))M<-M-tau2^2*G[g,previous,drop=FALSE]%*%
            solve(diag(length(previous))+tau2*G[previous,previous,drop=FALSE],G[previous,g,drop=FALSE])
          M<-(M+t(M))/2
          ld<-vapply(seq_len(nrow(table)),function(i) {
            s<-which(table[i,]==1L)
            if(length(s))chol_logdet(diag(length(s))+M[s,s,drop=FALSE]) else 0
          },numeric(1))
          assign(dkey,ld,determinant_cache)
        }
        ld<-get(dkey,determinant_cache)
        prior<-k*log(pi)+(length(groups[[j]])-k)*log1p(-pi)
        total<-total+rrs_log_sum_exp((prior-.5*ld)[keep])
        if(states[j]!="zero")previous<-c(previous,groups[[j]])
      }
      total
    }
    upper<-fit+min(local_sum(seq_along(groups)),local_sum(rev(seq_along(groups))))
    # As in the existing analytical bounds, pad floating-point matrix algebra.
    upper+1e-9*max(1,abs(upper),abs(fit))
  }
}
ars_make_gaussian_bounder<-function(d,score,budget) {
  upper<-ars_gaussian_mass_upper(d,score,budget)
  function(node,log_score_fn) {
    before<-unname(score$diagnostics()[["calls"]])
    node<-rrs_bound_tree_node(node,log_score_fn)
    u<-upper(node$region$states)
    if(node$log_lower>u+1e-8)stop("Gaussian evidence bound conflicts with the anchor lower bound")
    node$log_upper<-min(node$log_upper,u)
    node$bound_score_evaluations<-unname(score$diagnostics()[["calls"]])-before
    node
  }
}
