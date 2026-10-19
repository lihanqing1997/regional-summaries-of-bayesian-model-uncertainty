source("pilot.R")
make_product_region<-function(tables,allowed,groups,states) {
  probabilities<-lapply(seq_along(groups),function(j){p<-numeric(nrow(tables[[j]]$supports))
    p[allowed[[j]]]<-tables[[j]]$probabilities[allowed[[j]]]
    p/sum(p)})
  draw<-function(n) {
    x<-matrix(0L,n,sum(lengths(groups)))
    for(j in seq_along(groups)) {
      ids<-sample.int(length(allowed[[j]]),n,TRUE,probabilities[[j]][allowed[[j]]])
      x[,groups[[j]]]<-tables[[j]]$supports[allowed[[j]][ids],,drop=FALSE]
    };x
  }
  density<-function(x) {
    x<-as.matrix(x);out<-numeric(nrow(x))
    for(j in seq_along(groups)) {
      code<-1L+as.integer(x[,groups[[j]],drop=FALSE]%*%2^(seq_along(groups[[j]])-1))
      out<-out+log(probabilities[[j]][code])
    };out
  }
  anchor<-function(){x<-matrix(0L,1,sum(lengths(groups)))
    for(j in seq_along(groups))x[,groups[[j]]]<-tables[[j]]$supports[which.max(probabilities[[j]]),];x}
  structure(list(sample=draw,log_prob=density,log_prob_inside=density,anchor=anchor,
    groups=groups,states=states,cost=1L,label="envelope"),class="rrs_region")
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
      out[[j]]<-list(supports=table,probabilities=exp(logs-max(logs)))
      previous<-c(previous,possible[[j]])
    }
    names(out)<-names(groups)
    list(tables=out,logU=fit+norm)
  }
  a<-build(seq_along(groups));b<-build(rev(seq_along(groups)))
  z<-if(a$logU<=b$logU)a else b
  # Same floating-point padding as the audited Gaussian evidence bound.
  z$logU<-z$logU+1e-9*max(1,abs(z$logU),abs(fit))
  z$allowed<-lapply(keep,which);z$columns<-allowed
  z$region<-make_product_region(z$tables,z$allowed,groups,states)
  z
}
make_variational_envelope<-function(d,base) {
  X<-d$X_train;y<-d$y_train;sigma2<-d$sigma^2;tau2<-cfg$tau2
  states<-base$region$states;groups<-d$groups
  # Base probabilities are proportional to prior(S)/sqrt(det_lower(S)).
  tables<-base$tables
  keep<-lapply(seq_along(groups),function(j)seq_len(nrow(tables[[j]]$supports))%in%base$allowed[[j]])
  # Recover the unshifted determinant/prior coefficients directly. The normalizer
  # constants lost by table scaling are recovered from the old fit and log U.
  columns<-base$columns
  v0<-y
  if(length(columns))v0<-y-tau2*X[,columns,drop=FALSE]%*%solve(
    diag(length(columns))+tau2*crossprod(X[,columns,drop=FALSE]),crossprod(X[,columns,drop=FALSE],y))
  v0<-as.numeric(v0)
  fitA<--nrow(X)/2*log(2*pi*sigma2)-sum(y*v0)/(2*sigma2)
  old_norm<-sum(vapply(seq_along(groups),function(j)
    rrs_log_sum_exp(log(tables[[j]]$probabilities)[keep[[j]]]),numeric(1)))
  offset<-base$logU-fitA-old_norm
  last_v<-NULL;last<-NULL
  evaluate<-function(v) {
    if(identical(v,last_v))return(last)
    projection<-as.numeric(crossprod(X,v));total<-0;inclusion<-numeric(ncol(X));logs_out<-list()
    for(j in seq_along(groups)) {
      s<-tables[[j]]$supports
      logs<-log(tables[[j]]$probabilities)+as.numeric(s%*%(tau2*projection[groups[[j]]]^2/(2*sigma2)))
      norm<-rrs_log_sum_exp(logs[keep[[j]]]);prob<-numeric(nrow(s));prob[keep[[j]]]<-exp(logs[keep[[j]]]-norm)
      inclusion[groups[[j]]]<-as.numeric(crossprod(prob,s));total<-total+norm;logs_out[[j]]<-logs
    }
    value<-sum(v^2-2*y*v)/(2*sigma2)+total
    grad<-as.numeric((v-y+tau2*X%*%(inclusion*projection))/sigma2)
    last_v<<-v;last<<-list(value=value,gradient=grad,logs=logs_out);last
  }
  fit<-optim(v0,function(v)evaluate(v)$value,function(v)evaluate(v)$gradient,
    method="BFGS",control=list(maxit=100,reltol=1e-10))
  v<-if(fit$value<=evaluate(v0)$value)fit$par else v0
  result<-evaluate(v)
  improved<-lapply(seq_along(groups),function(j)list(supports=tables[[j]]$supports,
    probabilities=pmax(exp(result$logs[[j]]-max(result$logs[[j]])),.Machine$double.xmin)))
  names(improved)<-names(groups)
  logU<--nrow(X)/2*log(2*pi*sigma2)+result$value+offset
  list(tables=improved,region=make_product_region(improved,base$allowed,groups,states),
    allowed=base$allowed,columns=columns,
    logU=logU+1e-9*max(1,abs(logU)),optimization_convergence=fit$convergence)
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
  list(region=region,logU=logU,cells=length(cells))
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
  if(any(lw>1e-10))stop("Importance weights exceed the proven envelope")
  w<-exp(pmin(0,lw));ci<-bounded_mean_interval(w,delta)
  list(log_evidence=proposal$logU+log(ci$mean),
    log_lower=proposal$logU+log(ci$lower),log_upper=proposal$logU+log(ci$upper),
    score_evaluations=n,fit=list(method="bounded_envelope_iid",relative_radius=ci$relative_radius,
      n_paths=n,logU=proposal$logU,mean_normalized_weight=ci$mean,
      hoeffding_radius=ci$hoeffding_radius,bernstein_radius=ci$bernstein_radius))
}
