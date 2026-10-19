source("diagnostic.R")

make_shared_preparation <- function(d) {
  context <- new.env(parent=emptyenv())
  context$G <- crossprod(d$X_train)
  context$tables <- lapply(d$groups,function(g)as.matrix(expand.grid(rep(list(0:1),length(g)))))
  context$values <- new.env(hash=TRUE,parent=emptyenv())
  context$hits <- context$misses <- c(logdet=0L,factors=0L,initial=0L)
  memo <- function(type,key,calculate) {
    name <- paste(type,key,sep="|")
    if(exists(name,context$values,inherits=FALSE)) {
      context$hits[type] <- context$hits[type]+1L
      return(get(name,context$values,inherits=FALSE))
    }
    context$misses[type] <- context$misses[type]+1L
    value <- calculate();assign(name,value,context$values);value
  }
  base <- function(d,score,states,constraints=NULL) {
    G <- context$G;tau2 <- cfg$tau2;pi <- cfg$inclusion_probability
    groups <- d$groups;p <- ncol(d$X_train);tables <- context$tables
    keep <- lapply(seq_along(groups),function(j) {
      k <- rowSums(tables[[j]])
      if(!is.null(constraints))return(seq_along(k)%in%constraints[[j]])
      switch(states[j],free=rep(TRUE,length(k)),zero=k==0,one=k==1,multi=k>=2)
    })
    possible <- lapply(seq_along(groups),function(j)
      groups[[j]][colSums(tables[[j]][keep[[j]],,drop=FALSE])>0])
    allowed <- unlist(possible,use.names=FALSE)
    gamma <- integer(p);gamma[allowed] <- 1L
    ld <- memo("logdet",paste(allowed,collapse=","),function()
      if(length(allowed))2*sum(log(diag(chol(diag(length(allowed))+tau2*G[allowed,allowed,drop=FALSE])))) else 0)
    # Keep the original score request even if its matrix calculation is cached.
    fit <- score(gamma)+.5*ld-length(allowed)*log(pi)-(p-length(allowed))*log1p(-pi)
    build <- function(order) {
      previous <- integer();out <- vector("list",length(groups));norm <- 0
      for(j in order) {
        logs <- memo("factors",paste(j,paste(previous,collapse=","),sep=":"),function() {
          g <- groups[[j]];M <- tau2*G[g,g,drop=FALSE]
          if(length(previous))M <- M-tau2^2*G[g,previous,drop=FALSE]%*%
            solve(diag(length(previous))+tau2*G[previous,previous,drop=FALSE],G[previous,g,drop=FALSE])
          M <- (M+t(M))/2;table <- tables[[j]];k <- rowSums(table)
          logdet <- vapply(seq_len(nrow(table)),function(i) {
            s <- which(table[i,]==1)
            if(length(s))2*sum(log(diag(chol(diag(length(s))+M[s,s,drop=FALSE])))) else 0
          },numeric(1))
          k*log(pi)+(length(groups[[j]])-k)*log1p(-pi)-.5*logdet
        })
        norm <- norm+rrs_log_sum_exp(logs[keep[[j]]])
        out[[j]] <- list(supports=tables[[j]],probabilities=exp(logs-max(logs)))
        previous <- c(previous,possible[[j]])
      }
      names(out) <- names(groups);list(tables=out,logU=fit+norm)
    }
    a <- build(seq_along(groups));b <- build(rev(seq_along(groups)))
    z <- if(a$logU<=b$logU)a else b
    z$logU <- z$logU+1e-9*max(1,abs(z$logU),abs(fit))
    z$allowed <- lapply(keep,which);z$columns <- allowed
    z$region <- make_product_region(z$tables,z$allowed,groups,states);z
  }
  initial <- function(columns) memo("initial",paste(columns,collapse=","),function() {
    X <- d$X_train;y <- d$y_train;tau2 <- cfg$tau2
    v0 <- y
    if(length(columns))v0 <- y-tau2*X[,columns,drop=FALSE]%*%solve(
      diag(length(columns))+tau2*crossprod(X[,columns,drop=FALSE]),crossprod(X[,columns,drop=FALSE],y))
    v0
  })
  # Preserve every optimizer operation. Substitute only its deterministic
  # starting-vector calculation with a memoized copy of that same expression.
  variational <- make_variational_envelope;code <- body(variational);replaced <- 0L
  for(i in seq_along(code)) if(is.call(code[[i]]) &&
    identical(code[[i]][[1]],as.name("if")) &&
    identical(code[[i]][[2]],quote(length(columns)))) {
    code[[i]] <- quote(v0 <- cached_initial_vector(columns));replaced <- replaced+1L
  }
  stopifnot(replaced==1L);body(variational) <- code
  scope <- new.env(parent=environment(variational));scope$cached_initial_vector <- initial
  environment(variational) <- scope
  list(base=base,variational=variational,stats=function()
    list(hits=context$hits,misses=context$misses,entries=length(ls(context$values))))
}

bernoulli_kl <- function(p,q) {
  if(p==0)return(-log1p(-q))
  if(p==1)return(-log(q))
  if(q==0||q==1)return(Inf)
  p*log(p/q)+(1-p)*log((1-p)/(1-q))
}
kl_mean_interval <- function(mean,n,delta) {
  stopifnot(mean>=0,mean<=1,n>=2,delta>0,delta<1)
  cutoff <- log(2/delta)/n
  if(mean==0) return(c(lower=0,upper=-expm1(-cutoff)))
  if(mean==1) return(c(lower=exp(-cutoff),upper=1))
  # Retain outer brackets, and pad outward, rather than return an inward root.
  a <- 0;b <- mean
  for(i in 1:60) {m <- (a+b)/2;if(bernoulli_kl(mean,m)>cutoff)a <- m else b <- m}
  lower <- max(0,a-1e-12)
  a <- mean;b <- 1
  for(i in 1:60) {m <- (a+b)/2;if(bernoulli_kl(mean,m)>cutoff)b <- m else a <- m}
  c(lower=lower,upper=min(1,b+1e-12))
}
original_mean_interval <- bounded_mean_interval
refined_mean_interval <- function(w,delta) {
  old <- original_mean_interval(w,delta/2)
  kl <- kl_mean_interval(mean(w),length(w),delta/2)
  lower <- max(old$lower,kl[["lower"]]);upper <- min(old$upper,kl[["upper"]])
  stopifnot(lower<=upper)
  old$lower <- lower;old$upper <- upper
  old$radius <- max(old$mean-lower,upper-old$mean)
  old$relative_radius <- if(old$mean>0)old$radius/old$mean else Inf
  old$kl_lower <- kl[["lower"]];old$kl_upper <- kl[["upper"]];old
}

refinement_search <- function(d,seed,budget=25000L,seconds=60,interval="original") {
  stopifnot(interval%in%c("original","kl"))
  old_base <- make_envelope_proposal;old_variational <- make_variational_envelope
  old_interval <- bounded_mean_interval
  on.exit({assign("make_envelope_proposal",old_base,.GlobalEnv)
    assign("make_variational_envelope",old_variational,.GlobalEnv)
    assign("bounded_mean_interval",old_interval,.GlobalEnv)},add=TRUE)
  # Include cache construction in the elapsed search time reported below.
  t <- proc.time()[["elapsed"]];shared <- make_shared_preparation(d)
  setup <- proc.time()[["elapsed"]]-t
  assign("make_envelope_proposal",shared$base,.GlobalEnv)
  assign("make_variational_envelope",shared$variational,.GlobalEnv)
  if(interval=="kl")assign("bounded_mean_interval",refined_mean_interval,.GlobalEnv)
  ans <- timed_diagnostic_search(d,seed,budget,max(.001,seconds-setup))
  ans$total_seconds <- ans$total_seconds+setup
  ans$timing$total <- ans$timing$total+setup
  ans$timing$other <- ans$timing$other+setup
  ans$cache_setup_seconds <- setup;ans$cache_stats <- shared$stats();ans$interval <- interval
  ans
}
