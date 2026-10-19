# Source from core/ after production.R. This file is used only after searches.
Rcpp::sourceCpp('../native/scaling_reference.cpp',cacheDir='build',rebuild=FALSE,showOutput=FALSE)
make_scaling_reference <- function(d,seconds=900) {
  p <- ncol(d$X_train);h <- length(d$groups)
  group <- integer(p);for(j in seq_along(d$groups))group[d$groups[[j]]] <- j-1L
  set.seed(16005101L)
  masks <- sort(unique(c(0L,as.integer(2^p-1),sample.int(2^p,min(500,2^p))-1L)))
  ref <- scaling_exact_reference(crossprod(d$X_train),as.numeric(crossprod(d$X_train,d$y_train)),
    sum(d$y_train^2),nrow(d$X_train),d$sigma,cfg$tau2,cfg$inclusion_probability,
    group,masks,seconds)
  supports <- t(vapply(masks,function(i)as.integer(intToBits(i))[seq_len(p)],integer(p)))
  score <- new_score(d,FALSE)
  discrepancy <- max(abs(rrs_evaluate_log_score(score$score,supports)-ref$check_log_scores))
  stopifnot(discrepancy<1e-8,sum(ref$leaf_model_count)==2^p)
  ref$score_check_max_error <- discrepancy
  ref$states <- as.matrix(expand.grid(rep(list(0:2),h),KEEP.OUT.ATTRS=FALSE))
  ref$costs <- 1L+rowSums(1L+(ref$states>0))
  ref$posterior_mass <- exp(ref$leaf_log_mass-ref$logZ)
  stopifnot(abs(sum(ref$posterior_mass)-1)<1e-12)
  ref$prefix_log_mass <- scaling_prefix_log_masses(ref$leaf_log_mass,h)
  stopifnot(abs(ref$prefix_log_mass[[1L]]-ref$logZ)<1e-8)
  ref
}
