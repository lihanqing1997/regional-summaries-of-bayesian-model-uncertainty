make_scaling_data <- function(row) {
  p <- row$p;h <- p/row$group_size
  case <- data.frame(case_id=row$dataset,p=p,n=row$n,
    sizes=paste(rep(row$group_size,h),collapse='-'),design='random',
    between=row$between,active='1-2-3',signal='single',stress=FALSE)
  d <- production_data(case,row$data_seed)
  if(row$family=='stress') {
    covariance <- matrix(row$between,p,p)
    for(g in d$groups)covariance[g,g] <- cfg$within_correlation
    diag(covariance) <- 1
    old_signal <- as.numeric(d$X_train%*%d$beta);old_signal <- old_signal-mean(old_signal)
    target <- as.numeric(crossprod(d$beta,covariance%*%d$beta))
    beta <- numeric(p)
    beta[vapply(d$groups[1:5],`[`,integer(1),1L)] <- c(1.4,-1.12,.98,-.84,.70)
    beta <- beta*sqrt(target/as.numeric(crossprod(beta,covariance%*%beta)))
    signal <- as.numeric(d$X_train%*%beta);signal <- signal-mean(signal)
    d$y_train <- d$y_train-old_signal+signal;d$beta <- beta
    d$signal_variance <- as.numeric(crossprod(beta,covariance%*%beta))
    stopifnot(abs(d$signal_variance-target)<1e-10)
  }
  d
}
