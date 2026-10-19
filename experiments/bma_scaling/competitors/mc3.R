# Independent model-space MCMC and empirical regional-summary selection.
# No ARS tree, bound, envelope, confidence sequence, or allocator is used.
# The calling harness supplies the common posterior score, so targets match.

mc3_select_regions <- function(weights,costs,target,memory_mb=192) {
  stopifnot(length(weights)==length(costs),all(is.finite(weights)),all(weights>=0),
            all(costs>=1),all(costs==as.integer(costs)),target>0,target<=sum(weights))
  ids <- which(weights>0)
  candidates <- lapply(list(order(-weights[ids],ids),order(-weights[ids]/costs[ids],ids)),function(o) {
    ranked <- ids[o];ranked[seq_len(which(cumsum(weights[ranked])>=target)[1L])]
  })
  incumbent <- candidates[[which.min(vapply(candidates,function(z)sum(costs[z]),numeric(1)))]]
  cap <- as.integer(sum(costs[incumbent]))
  # A set of cost <= cap needs at most floor(cap/c) items of any cost c.
  kept <- sort(unlist(lapply(split(ids,costs[ids]),function(z)
    head(z[order(-weights[z],z)],cap%/%costs[z[1L]])),use.names=FALSE))
  if(4*length(kept)*(cap+1)>memory_mb*1024^2)
    return(list(selected=incumbent,cost=cap,status='greedy-memory-fallback'))
  best <- rep(-Inf,cap+1L);best[1L] <- 0
  decisions <- matrix(FALSE,length(kept),cap+1L)
  for(i in seq_along(kept)) {
    j <- kept[i];cost <- costs[j];source <- seq_len(cap+1L-cost);destination <- source+cost
    proposed <- best[source]+weights[j];improve <- proposed>best[destination]
    best[destination[improve]] <- proposed[improve]
    decisions[i,destination[improve]] <- TRUE
  }
  possible <- which(best>=target)
  if(!length(possible))return(list(selected=incumbent,cost=cap,status='greedy-roundoff-fallback'))
  budget <- possible[1L];selected <- integer()
  for(i in rev(seq_along(kept)))if(decisions[i,budget]) {
    selected <- c(kept[i],selected);budget <- budget-costs[kept[i]]
  }
  if(sum(weights[selected])<target)
    return(list(selected=incumbent,cost=cap,status='greedy-roundoff-fallback'))
  list(selected=selected,cost=sum(costs[selected]),status='empirical-optimum')
}

mc3_region_code <- function(gamma,groups) {
  states <- vapply(groups,function(g)min(2L,sum(gamma[g])),numeric(1))
  sum(states*3^(seq_along(groups)-1L))
}
mc3_region_cost <- function(code,h) {
  1L+sum(1L+(floor(code/3^(seq_len(h)-1L))%%3>0))
}

mc3_search <- function(score,p,groups,budget,seconds,seed,epsilon=.1,
                      burn_fraction=.25,chains=4L,inclusion_probability=.1) {
  stopifnot(p>=2L,p<=30L,budget>=chains*8L,chains==4L,seconds>0)
  started <- proc.time()[['elapsed']];deadline <- started+seconds
  set.seed(seed)
  current <- cbind(rep(0L,p),rep(1L,p),rbinom(p,1,inclusion_probability),rbinom(p,1,.5))
  values <- vapply(seq_len(chains),function(j)score(current[,j]),numeric(1))
  calls <- chains;iterations <- as.integer(budget-chains)
  codes <- numeric(iterations);chain_ids <- integer(iterations);model_size <- integer(iterations)
  accepted <- integer(chains);proposed <- integer(chains);chain_steps <- integer(chains)
  chain_kept <- integer(chains)
  burn <- floor((iterations%/%chains)*burn_fraction)
  visited <- new.env(hash=TRUE,parent=emptyenv())
  save_state <- function(gamma,value) {
    key <- paste0(gamma,collapse='')
    if(!exists(key,visited,inherits=FALSE))
      assign(key,c(code=mc3_region_code(gamma,groups),log_score=value),visited)
  }
  for(j in seq_len(chains))save_state(current[,j],values[j])
  count <- 0L;stopped <- FALSE
  for(t in seq_len(iterations)) {
    if(t%%256L==1L && proc.time()[['elapsed']]>deadline){stopped<-TRUE;break}
    chain <- 1L+(t-1L)%%chains;candidate <- current[,chain]
    if(runif(1)<.5) {
      j <- sample.int(p,1L);candidate[j] <- 1L-candidate[j]
    } else {
      active <- which(candidate==1L);inactive <- which(candidate==0L)
      if(length(active)&&length(inactive)) {
        a <- active[sample.int(length(active),1L)];b <- inactive[sample.int(length(inactive),1L)]
        candidate[c(a,b)] <- candidate[c(b,a)]
      }
      # At empty/full supports the swap branch is a self-transition.
    }
    value <- score(candidate);calls <- calls+1L;proposed[chain] <- proposed[chain]+1L
    # Toggle probabilities are symmetric; swaps preserve model size and the
    # number of available swap pairs, so the reverse proposal has equal mass.
    if(log(runif(1))<min(0,value-values[chain])) {
      current[,chain] <- candidate;values[chain] <- value;accepted[chain] <- accepted[chain]+1L
    }
    save_state(current[,chain],values[chain])
    chain_steps[chain] <- chain_steps[chain]+1L
    if(chain_steps[chain]>burn) {
      count <- count+1L;codes[count] <- mc3_region_code(current[,chain],groups)
      chain_ids[count] <- chain;model_size[count] <- sum(current[,chain]);chain_kept[chain] <- chain_kept[chain]+1L
    }
  }
  sampling_seconds <- proc.time()[['elapsed']]-started
  if(count==0L)return(list(status='no-post-burn-in-samples',calls=calls,
    total_seconds=sampling_seconds,certified=FALSE))
  codes <- codes[seq_len(count)];chain_ids <- chain_ids[seq_len(count)];model_size <- model_size[seq_len(count)]
  tab <- table(codes);region_codes <- as.numeric(names(tab));counts <- as.numeric(tab)
  costs <- vapply(region_codes,mc3_region_cost,numeric(1),h=length(groups))
  # Integer frequencies give an exact empirical count threshold.
  frequency <- mc3_select_regions(counts,costs,ceiling((1-epsilon)*sum(counts)))
  frequency$region_codes <- region_codes[frequency$selected]
  frequency$estimated_mass <- sum(counts[frequency$selected])/sum(counts)
  unique_models <- do.call(rbind,mget(ls(visited,all.names=TRUE),visited,inherits=FALSE))
  shifted <- exp(unique_models[,'log_score']-max(unique_models[,'log_score']))
  regional <- tapply(shifted,unique_models[,'code'],sum)
  region_codes_rm <- as.numeric(names(regional));rm_weights <- as.numeric(regional)/sum(regional)
  rm_costs <- vapply(region_codes_rm,mc3_region_cost,numeric(1),h=length(groups))
  renormalized <- mc3_select_regions(rm_weights,rm_costs,1-epsilon)
  renormalized$region_codes <- region_codes_rm[renormalized$selected]
  renormalized$estimated_mass <- sum(rm_weights[renormalized$selected])
  list(status=if(stopped)'time-limit' else 'budget-complete',calls=calls,
       total_seconds=proc.time()[['elapsed']]-started,sampling_seconds=sampling_seconds,
       selection_seconds=proc.time()[['elapsed']]-started-sampling_seconds,
       certified=FALSE,frequency=frequency,renormalized=renormalized,
       draws=count,chain_draws=chain_kept,chain_steps=chain_steps,
       acceptance=accepted/pmax(1,proposed),burn_per_chain=burn,
       region_codes=codes,chain_ids=chain_ids,model_size=model_size,
       unique_models=nrow(unique_models),region_frequency_codes=region_codes,
       region_frequency=counts/sum(counts),renormalized_codes=region_codes_rm,
       renormalized_mass=rm_weights,seed=seed)
}
