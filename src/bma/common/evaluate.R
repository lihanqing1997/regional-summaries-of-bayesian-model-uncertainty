# Run only after sourcing adaptive.R. References are loaded after all practical
# methods finish. Every inspected stochastic look is audited, not just the last.
audit_all_looks <- function(ans,reference) {
  events <- ans$sampling_detail
  if(is.null(events)||!nrow(events)) return(data.frame())
  checks <- lapply(seq_len(nrow(events)),function(i) {
    states <- strsplit(events$states[i],"/",fixed=TRUE)[[1]]
    keep <- rep(TRUE,nrow(reference$states))
    for(j in which(states!="free")) keep <- keep &
      reference$states[,j]==match(states[j],c("zero","one","multi"))-1L
    z <- rrs_log_sum_exp(reference$leaf_log_mass[keep]);tol <- 1e-9*max(1,abs(z))
    data.frame(node=events$node[i],look=events$look[i]%||%1L,exact_log_mass=z,
      covered=events$log_lower[i]<=z+tol && z<=events$log_upper[i]+tol)
  })
  do.call(rbind,checks)
}

result_row <- function(job,method,ans,audit,looks) {
  s <- ans$snapshot$state;a <- audit$summary;e <- ans$sampling_detail
  data.frame(job=job,method=method,status=ans$status,error=ans$error,
    seconds=ans$total_seconds,calls=ans$total_score_calls,cost_lower=s$cost_lower,
    cost_upper=s$cost_upper,exact_cost=a$exact_optimal_cost,
    retained_mass=a$retained_mass,verified_feasible=a$mass_guarantee_verified,
    certified_optimal=a$verified_optimal,audit=a$status,
    all_looks_covered=if(nrow(looks))all(looks$covered) else TRUE,
    sampled_regions=if(is.null(e)||!nrow(e))0L else length(unique(e$node)),
    looks=if(is.null(e))0L else nrow(e),
    preparation_calls=if(is.null(e))0L else sum(e$preparation_calls),
    sampling_calls=if(is.null(e))0L else sum(e$sampling_calls))
}
