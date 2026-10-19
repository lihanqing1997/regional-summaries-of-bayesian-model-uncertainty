source("base_study.R")
# Use the configured toolchain; optional override for nonstandard Windows installs.
if (.Platform$OS.type == "windows") {
  toolroot <- Sys.getenv("ARS_RTOOLS_HOME", Sys.getenv("RTOOLS44_HOME"))
  if (!nzchar(toolroot)) {
    candidates <- c("C:/rtools44", "C:/Program Files/R/rtools44")
    found <- candidates[dir.exists(candidates)]
    if (length(found)) toolroot <- found[1L]
  }
  if (nzchar(toolroot)) {
    Sys.setenv(RTOOLS44_HOME=toolroot)
    Sys.setenv(PATH=paste(file.path(toolroot,"usr/bin"),
      file.path(toolroot,"x86_64-w64-mingw32.static.posix/bin"),Sys.getenv("PATH"),sep=";"))
  }
}
Sys.setenv(PKG_LIBS="-lmpfr -lgmp")
Rcpp::sourceCpp("src/selection.cpp",cacheDir="build",rebuild=FALSE,showOutput=FALSE)
old_select_region_union_epsilon<-select_region_union_epsilon
select_region_union_epsilon<-function(log_evidence_lower,costs,epsilon,
  log_evidence_upper=log_evidence_lower,mode=c("conservative","optimistic"),eligible=NULL) {
  mode<-match.arg(mode)
  v<-rrs_validate_interval_frontier_inputs(log_evidence_lower,log_evidence_upper,costs,epsilon,eligible)
  answer<-ars_interval_knapsack(v$log_lower,v$log_upper,v$costs,v$eligible,v$epsilon,
                               mode=="conservative",memory_mb=cfg$dp_matrix_mb)
  interval<-if(length(answer$selected))
    rrs_interval_omitted_mass(v$log_lower,v$log_upper,answer$selected) else
    list(log_lower_selected=-Inf,log_upper_selected=-Inf,lower=1,upper=1)
  list(feasible=answer$feasible,selected=answer$selected,
    budget=if(mode=="conservative")answer$upper_cost else answer$lower_cost,
    solver_lower_cost=answer$lower_cost,solver_upper_cost=answer$upper_cost,
    solver_status=answer$status,precision_bits=answer$bits,
    log_evidence_lower=interval$log_lower_selected,log_evidence_upper=interval$log_upper_selected,
    omitted_mass_lower=interval$lower,omitted_mass_upper=interval$upper,
    epsilon=epsilon,mode=mode,frontier=NULL)
}
source("src/search_repair.R")
source("src/gaussian_mass_bound.R")
