source("pilot.R")
make_product_region<-function(tables,allowed,groups,states) {
  normalizations<-lapply(seq_along(groups),function(j)
    rrs_stable_local_probabilities(rrs_table_log_weights(tables[[j]]),allowed[[j]]))
  probabilities<-lapply(normalizations,`[[`,"probabilities")
  draw<-function(n) {
    x<-matrix(0L,n,sum(lengths(groups)))
    for(j in seq_along(groups)) {
      ids<-sample.int(length(allowed[[j]]),n,TRUE,probabilities[[j]][allowed[[j]]])
      x[,groups[[j]]]<-tables[[j]]$supports[allowed[[j]][ids],,drop=FALSE]
    };x
  }
  density<-function(x) {
    x<-as.matrix(x);out<-numeric(nrow(x))
    if(ncol(x)!=sum(lengths(groups)) || any(!x %in% c(0L,1L)))
      stop("Invalid support passed to product proposal")
    for(j in seq_along(groups)) {
      code<-1L+as.integer(x[,groups[[j]],drop=FALSE]%*%2^(seq_along(groups[[j]])-1))
      out<-out+normalizations[[j]]$log_probabilities[code]
    };out
  }
  anchor<-function(){x<-matrix(0L,1,sum(lengths(groups)))
    for(j in seq_along(groups))x[,groups[[j]]]<-tables[[j]]$supports[which.max(probabilities[[j]]),];x}
  structure(list(sample=draw,log_prob=density,log_prob_inside=density,anchor=anchor,
    groups=groups,states=states,cost=1L,label="envelope",
    local_normalizations=normalizations,
    log_factor_ratio=sum(vapply(normalizations,`[[`,numeric(1),"log_factor_ratio"))),
    class="rrs_region")
}

finalize_base_envelope<-function(z,fit,allowed,columns,groups,states) {
  z$allowed<-allowed;z$columns<-columns
  z$region<-make_product_region(z$tables,allowed,groups,states)
  z$ideal_logU<-z$logU
  # Use the sampled proposal, including any positive tail safeguard.
  actual<-fit+z$region$log_factor_ratio
  z$logU<-actual+1e-9*max(1,abs(actual),abs(fit))
  z$log_factor_constant<-fit
  z
}
make_envelope_proposal<-function(d,score,states,constraints=NULL) {
  G<-crossprod(d$X_train);tau2<-cfg$tau2;pi<-cfg$inclusion_probability
  groups<-d$groups;p<-ncol(d$X_train)
  tables<-lapply(groups,function(g)as.matrix(expand.grid(rep(list(0:1),length(g)))))
  keep<-lapply(seq_along(groups),function(j){k<-rowSums(tables[[j]])
    if(!is.null(constraints))return(seq_along(k)%in%constraints[[j]])
    switch(states[j],free=rep(TRUE,length(k)),zero=k==0,one=k==1,multi=k>=2)})
  possible<-lapply(seq_along(groups),function(j)groups[[j]][colSums(tables[[j]][keep[[j]],,drop=FALSE])>0])
  allowed<-unlist(possible,use.names=FALSE)
  gamma<-integer(p);gamma[allowed]<-1L
  ld<-if(length(allowed))2*sum(log(diag(chol(diag(length(allowed))+tau2*G[allowed,allowed,drop=FALSE])))) else 0
  fit<-score(gamma)+.5*ld-length(allowed)*log(pi)-(p-length(allowed))*log1p(-pi)
  build<-function(order) {
    previous<-integer();out<-vector("list",length(groups));norm<-0
    for(j in order) {
      g<-groups[[j]];M<-tau2*G[g,g,drop=FALSE]
      if(length(previous))M<-M-tau2^2*G[g,previous,drop=FALSE]%*%
        solve(diag(length(previous))+tau2*G[previous,previous,drop=FALSE],G[previous,g,drop=FALSE])
      M<-(M+t(M))/2;table<-tables[[j]];k<-rowSums(table)
      logdet<-vapply(seq_len(nrow(table)),function(i){s<-which(table[i,]==1)
        if(length(s))2*sum(log(diag(chol(diag(length(s))+M[s,s,drop=FALSE])))) else 0},numeric(1))
      logs<-k*log(pi)+(length(g)-k)*log1p(-pi)-.5*logdet
      norm<-norm+rrs_log_sum_exp(logs[keep[[j]]])
      out[[j]]<-list(supports=table,log_weights=logs,
        probabilities=rrs_stable_local_probabilities(logs)$probabilities)
      previous<-c(previous,possible[[j]])
    }
    names(out)<-names(groups)
    list(tables=out,logU=fit+norm)
  }
  a<-build(seq_along(groups));b<-build(rev(seq_along(groups)))
  z<-if(a$logU<=b$logU)a else b
  finalize_base_envelope(z,fit,lapply(keep,which),allowed,groups,states)
}
make_variational_envelope<-function(d,base,initial_vector=NULL,span_basis=NULL,control=list()) {
  X<-d$X_train;y<-as.numeric(d$y_train);sigma2<-d$sigma^2;tau2<-cfg$tau2
  states<-base$region$states;groups<-d$groups;tables<-base$tables;columns<-base$columns
  keep<-base$allowed
  method<-control$method %||% cfg$envelope_optimizer %||% "newton_span"
  tolerance<-control$gap_tolerance %||% 1e-6
  maxit<-control$maxit %||% 60L
  if(is.null(initial_vector)) {
    v0<-y
    if(length(columns))v0<-y-tau2*X[,columns,drop=FALSE]%*%solve(
      diag(length(columns))+tau2*crossprod(X[,columns,drop=FALSE]),crossprod(X[,columns,drop=FALSE],y))
    v0<-as.numeric(v0)
  } else v0<-as.numeric(initial_vector)
  # A full QR basis (without numerical rank truncation) contains every allowed
  # column. The optimizer's orthogonal component is exactly that of y.
  if(is.null(span_basis)) {
    Q<-if(length(columns))qr.Q(qr(X[,columns,drop=FALSE],LAPACK=TRUE)) else matrix(0,nrow(X),0)
  } else Q<-span_basis
  r<-ncol(Q);Z<-crossprod(Q,X)
  yperp<-if(r)y-as.numeric(Q%*%crossprod(Q,y)) else y
  local_logs<-lapply(tables,rrs_table_log_weights)
  evaluate<-function(theta,need_hessian=TRUE) {
    v<-yperp+if(r)as.numeric(Q%*%theta) else rep(0,length(y))
    projection<-as.numeric(crossprod(X,v));total<-0;inclusion<-numeric(ncol(X));logs_out<-list()
    H<-if(need_hessian)diag(r)/sigma2 else NULL
    for(j in seq_along(groups)) {
      g<-groups[[j]];s<-tables[[j]]$supports
      logs<-local_logs[[j]]+as.numeric(s%*%(tau2*projection[g]^2/(2*sigma2)))
      norm<-rrs_log_sum_exp(logs[keep[[j]]])
      # No floored probabilities enter derivatives of the smooth ideal F.
      prob<-exp(logs[keep[[j]]]-norm);ss<-s[keep[[j]],,drop=FALSE]
      inc<-as.numeric(crossprod(prob,ss));inclusion[g]<-inc
      total<-total+norm;logs_out[[j]]<-logs
      if(need_hessian && r) {
        Zg<-Z[,g,drop=FALSE]
        covariance<-crossprod(ss,ss*prob)-tcrossprod(inc)
        B<-sweep(Zg,2,projection[g],"*")
        H<-H+(tau2/sigma2)*tcrossprod(sweep(Zg,2,sqrt(pmax(inc,0)),"*"))+
          (tau2/sigma2)^2*B%*%covariance%*%t(B)
      }
    }
    value<-sum(v^2-2*y*v)/(2*sigma2)+total
    gradient<-as.numeric((v-y+tau2*X%*%(inclusion*projection))/sigma2)
    list(value=value,gradient=gradient,projected_gradient=as.numeric(crossprod(Q,gradient)),
      hessian=if(need_hessian)(H+t(H))/2 else NULL,logs=logs_out,v=v,
      gap_bound=sigma2*sum(gradient^2)/2)
  }
  theta<-if(r)as.numeric(crossprod(Q,v0)) else numeric()
  initial<-evaluate(theta);result<-initial;iterations<-0L;backtracks<-0L
  trace<-data.frame(iteration=0L,value=result$value,gap=result$gap_bound)
  reason<-if(result$gap_bound<=tolerance)"gradient-gap" else "iteration-limit"
  if(method=="bfgs") {
    if(r) {
      fit<-optim(theta,function(t)evaluate(t,FALSE)$value,
        function(t)evaluate(t,FALSE)$projected_gradient,method="BFGS",
        control=list(maxit=100,reltol=1e-10))
      candidate<-evaluate(fit$par)
      if(candidate$value<=result$value){theta<-fit$par;result<-candidate}
      iterations<-unname(fit$counts[[1]])
      reason<-paste0("bfgs-",fit$convergence)
    }
  } else if(method=="newton_span") {
    for(iteration in seq_len(maxit)) {
      if(result$gap_bound<=tolerance || !r) {reason<-"gradient-gap";break}
      direction<-tryCatch(as.numeric(solve(result$hessian,-result$projected_gradient)),
        error=function(e)-sigma2*result$projected_gradient)
      slope<-sum(result$projected_gradient*direction)
      if(!is.finite(slope) || slope>=0) {
        direction<--sigma2*result$projected_gradient;slope<-sum(result$projected_gradient*direction)
      }
      accepted<-FALSE;step<-1
      for(line in 0:40) {
        trial<-evaluate(theta+step*direction,FALSE)
        if(is.finite(trial$value) && trial$value<=result$value+1e-4*step*slope) {
          accepted<-TRUE;break
        }
        step<-step/2;backtracks<-backtracks+1L
      }
      if(!accepted){reason<-"line-search-limit";break}
      theta<-theta+step*direction;result<-evaluate(theta);iterations<-iteration
      trace<-rbind(trace,data.frame(iteration=iteration,value=result$value,gap=result$gap_bound))
      if(result$gap_bound<=tolerance){reason<-"gradient-gap";break}
    }
  } else stop("Unknown envelope optimizer")
  improved<-lapply(seq_along(groups),function(j)list(supports=tables[[j]]$supports,
    log_weights=result$logs[[j]],
    probabilities=rrs_stable_local_probabilities(result$logs[[j]])$probabilities))
  names(improved)<-names(groups)
  region<-make_product_region(improved,base$allowed,groups,states)
  constant<--nrow(X)/2*log(2*pi*sigma2)+sum(result$v^2-2*y*result$v)/(2*sigma2)
  ideal_logU<--nrow(X)/2*log(2*pi*sigma2)+result$value
  actual_logU<-constant+region$log_factor_ratio
  diagnostics<-list(method=method,dimension=r,ambient_dimension=nrow(X),iterations=iterations,
    backtracks=backtracks,reason=reason,gradient_norm=sqrt(sum(result$gradient^2)),
    gap_bound=result$gap_bound,gap_tolerance=tolerance,
    gradient_target_met=result$gap_bound<=tolerance,initial_value=initial$value,
    final_value=result$value,ideal_logU=ideal_logU,
    final_vector=result$v,gradient=result$gradient,
    proposal_correction=actual_logU-ideal_logU,trace=trace)
  list(tables=improved,region=region,allowed=base$allowed,columns=columns,
    logU=actual_logU+1e-9*max(1,abs(actual_logU)),ideal_logU=ideal_logU,
    optimization_convergence=if(diagnostics$gradient_target_met)0L else 1L,
    optimization=diagnostics)
}
make_partitioned_envelope<-function(d,score,states,max_cells=256L) {
  initial<-lapply(seq_along(d$groups),function(j){m<-length(d$groups[[j]])
    x<-as.matrix(expand.grid(rep(list(0:1),m)));k<-rowSums(x)
    which(switch(states[j],free=rep(TRUE,length(k)),zero=k==0,one=k==1,multi=k>=2))})
  cells<-list(initial)
  for(j in seq_along(initial)) {
    if(length(cells)*length(initial[[j]])>max_cells)break
    cells<-unlist(lapply(cells,function(cell)lapply(initial[[j]],function(i){cell[[j]]<-i;cell})),recursive=FALSE)
  }
  proposals<-lapply(cells,function(cell)make_variational_envelope(d,make_envelope_proposal(d,score,states,cell)))
  logs<-vapply(proposals,`[[`,numeric(1),"logU");logU<-rrs_log_sum_exp(logs)
  probability<-pmax(exp(logs-logU),.Machine$double.xmin);probability<-probability/sum(probability)
  logU<-max(logs-log(probability))+1e-9*max(1,abs(logU))
  region<-proposals[[1]]$region
  region$sample<-function(n) {
    selected<-sample.int(length(cells),n,TRUE,probability);x<-matrix(0L,n,ncol(d$X_train))
    for(i in unique(selected)){rows<-which(selected==i);x[rows,]<-proposals[[i]]$region$sample(length(rows))};x
  }
  region$log_prob<-region$log_prob_inside<-function(x) {
    out<-rep(-Inf,nrow(as.matrix(x)))
    for(i in seq_along(cells)) {
      value<-log(probability[i])+proposals[[i]]$region$log_prob(x)
      # Cells are disjoint, so each support has exactly one finite density.
      inside<-is.finite(value);out[inside]<-value[inside]
    };out
  }
  diagnostics<-lapply(proposals,`[[`,"optimization")
  list(region=region,logU=logU,cells=length(cells),
    optimization=list(methods=unique(vapply(diagnostics,`[[`,character(1),"method")),
      max_gap_bound=max(vapply(diagnostics,`[[`,numeric(1),"gap_bound")),
      target_met=sum(vapply(diagnostics,`[[`,logical(1),"gradient_target_met")),
      cells=length(cells),iterations=sum(vapply(diagnostics,`[[`,numeric(1),"iterations")),
      max_dimension=max(vapply(diagnostics,`[[`,numeric(1),"dimension"))),
    cell_optimization=diagnostics)
}
make_conditional_mixture<-function(d,score,states,envelope,weight=.5) {
  groups<-d$groups;p<-ncol(d$X_train);anchor<-envelope$region$anchor()
  allowed<-lapply(seq_along(groups),function(j) {
    x<-envelope$tables[[j]]$supports;k<-rowSums(x)
    x[switch(states[j],free=rep(TRUE,length(k)),zero=k==0,one=k==1,multi=k>=2),,drop=FALSE]
  })
  candidates<-function(j,center) {
    out<-matrix(rep(center,each=nrow(allowed[[j]])),nrow=nrow(allowed[[j]]),ncol=p)
    out[,groups[[j]]]<-allowed[[j]];out
  }
  # Fixed deterministic preparation; inference draws are fresh afterwards.
  for(pass in 1:2)for(j in seq_along(groups)) {
    x<-candidates(j,anchor);logs<-rrs_evaluate_log_score(score,x)
    anchor<-x[which.max(logs),,drop=FALSE]
  }
  tables<-lapply(seq_along(groups),function(j) {
    x<-candidates(j,anchor);logs<-rrs_evaluate_log_score(score,x)
    table<-envelope$tables[[j]]$supports
    probs<-rep(.Machine$double.xmin,nrow(table))
    code<-1L+as.integer(allowed[[j]]%*%2^(seq_len(ncol(table))-1))
    probs[code]<-pmax(exp(logs-max(logs)),.Machine$double.xmin)
    list(supports=table,probabilities=probs)
  })
  names(tables)<-names(groups)
  adapted<-new_group_table_region(tables,states,groups,cost=1L)
  q<-envelope$region
  q$sample<-function(n) {
    take<-runif(n)<weight;out<-matrix(0L,n,p)
    if(any(take))out[take,]<-adapted$sample(sum(take))
    if(any(!take))out[!take,]<-envelope$region$sample(sum(!take))
    out
  }
  q$log_prob<-q$log_prob_inside<-function(x) {
    a<-log1p(-weight)+envelope$region$log_prob(x);b<-log(weight)+adapted$log_prob(x)
    m<-pmax(a,b);m+log(exp(a-m)+exp(b-m))
  }
  list(region=q,logU=envelope$logU-log1p(-weight),anchor=anchor)
}
bounded_mean_interval<-function(w,delta) {
  n<-length(w);stopifnot(n>=2,delta>0,delta<1,all(is.finite(w)),all(w>=0),all(w<=1))
  mean<-mean(w);variance<-var(w)
  # Two-sided Hoeffding gets delta/2. Each one-sided empirical Bernstein
  # bound gets delta/4. Their intersection therefore has failure <= delta.
  h<-sqrt(log(4/delta)/(2*n));L<-log(8/delta)
  eb<-sqrt(2*variance*L/n)+7*L/(3*(n-1));radius<-min(h,eb)
  list(mean=mean,lower=max(0,mean-radius),upper=min(1,mean+radius),
       radius=radius,hoeffding_radius=h,bernstein_radius=eb,
       relative_radius=if(mean>0)radius/mean else Inf,n=n)
}
sample_envelope_evidence<-function(proposal,score,n,delta) {
  x<-proposal$region$sample(n)
  logs<-rrs_evaluate_log_score(score,x)-proposal$region$log_prob(x)
  lw<-logs-proposal$logU
  if(any(!is.finite(lw)) || any(lw>0))stop("Importance weights exceed the proven envelope")
  w<-exp(lw);ci<-bounded_mean_interval(w,delta)
  list(log_evidence=proposal$logU+log(ci$mean),
    log_lower=proposal$logU+log(ci$lower),log_upper=proposal$logU+log(ci$upper),
    score_evaluations=n,fit=list(method="bounded_envelope_iid",relative_radius=ci$relative_radius,
      n_paths=n,logU=proposal$logU,mean_normalized_weight=ci$mean,
      hoeffding_radius=ci$hoeffding_radius,bernstein_radius=ci$bernstein_radius))
}
