# Published Bayesian-model-set comparators for the exact-reference BMA study.
#
# Occam's Window follows Madigan and Raftery (1994), equations (5)--(7):
# first retain models within a factor c of the posterior mode, then remove a
# model if a proper submodel remaining in that window has larger posterior
# probability.  The default c = 20 is the value used in their examples.
#
# Cartesian credible sets follow Griffin (2025), Algorithms 1--2, using exact
# posterior marginals rather than their Monte Carlo estimates.  Low-PIP
# variables are screened at 0.04, blocks are merged greedily by the largest
# reduction in KL divergence (their mutual information), and the published
# M = 2 ease-of-understanding penalty is used.  Coverage under the factorized
# approximation and coverage under the original posterior are both returned;
# only the latter is a fair exact-reference coverage audit.

pbmc_validate_posterior <- function(supports, posterior) {
  supports <- as.matrix(supports)
  storage.mode(supports) <- "integer"
  posterior <- as.numeric(posterior)
  if (nrow(supports) < 1L || ncol(supports) < 1L ||
      length(posterior) != nrow(supports) ||
      any(!supports %in% c(0L, 1L)) ||
      any(!is.finite(posterior)) || any(posterior < 0) ||
      !is.finite(sum(posterior)) || sum(posterior) <= 0) {
    stop("Invalid support posterior.", call. = FALSE)
  }
  keys <- apply(supports, 1L, paste0, collapse = "")
  if (anyDuplicated(keys)) {
    stop("Support rows must be distinct.", call. = FALSE)
  }
  list(supports = supports, posterior = posterior / sum(posterior), keys = keys)
}

pbmc_enumerate_supports <- function(p) {
  p <- as.integer(p)
  if (length(p) != 1L || is.na(p) || p < 0L || p > 30L) {
    stop("p must be one integer between zero and 30.", call. = FALSE)
  }
  if (p == 0L) return(matrix(integer(0), nrow = 1L, ncol = 0L))
  ids <- 0:(2^p - 1)
  out <- sapply(seq_len(p), function(j) {
    as.integer(bitwAnd(ids, bitwShiftL(1L, j - 1L)) != 0L)
  })
  matrix(out, ncol = p)
}

pbmc_support_masks <- function(supports) {
  p <- ncol(supports)
  if (p > 30L) stop("Bit-mask helpers support at most 30 variables.", call. = FALSE)
  as.integer(supports %*% 2^(seq_len(p) - 1L))
}

pbmc_evaluate_set <- function(supports, posterior, selected, epsilon) {
  checked <- pbmc_validate_posterior(supports, posterior)
  supports <- checked$supports
  posterior <- checked$posterior
  selected <- as.logical(selected)
  epsilon <- as.numeric(epsilon)
  if (length(selected) != nrow(supports) || anyNA(selected) ||
      length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 1) {
    stop("Invalid selected set or omitted-mass tolerance.", call. = FALSE)
  }
  retained_mass <- sum(posterior[selected])
  full_pip <- colSums(supports * posterior)
  if (retained_mass > 0) {
    restricted <- posterior * selected / retained_mass
    pip_error <- abs(colSums(supports * restricted) - full_pip)
  } else {
    pip_error <- rep(NA_real_, ncol(supports))
  }
  list(
    target_mass = 1 - epsilon,
    retained_mass = retained_mass,
    omitted_mass = 1 - retained_mass,
    feasible = retained_mass + 1e-12 >= 1 - epsilon,
    full_model_count = sum(selected),
    full_space_fraction = mean(selected),
    total_variation_to_restricted = 1 - retained_mass,
    reverse_kullback_leibler = if (retained_mass > 0) -log(retained_mass) else Inf,
    maximum_pip_error = if (retained_mass > 0) max(pip_error) else NA_real_,
    mean_pip_error = if (retained_mass > 0) mean(pip_error) else NA_real_
  )
}

pbmc_hpp_credible_set <- function(supports, posterior, epsilon) {
  checked <- pbmc_validate_posterior(supports, posterior)
  epsilon <- as.numeric(epsilon)
  if (length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 1) {
    stop("epsilon must lie in [0,1).", call. = FALSE)
  }
  order_index <- order(-checked$posterior, checked$keys)
  number <- which(cumsum(checked$posterior[order_index]) >= 1 - epsilon)[1L]
  selected <- rep(FALSE, nrow(checked$supports))
  selected[order_index[seq_len(number)]] <- TRUE
  list(
    method = "highest-posterior-probability set",
    selected = selected,
    selected_indices = which(selected),
    description_entries = sum(selected),
    evaluation = pbmc_evaluate_set(
      checked$supports, checked$posterior, selected, epsilon
    )
  )
}

pbmc_occams_window <- function(supports, posterior, c = 20) {
  checked <- pbmc_validate_posterior(supports, posterior)
  supports <- checked$supports
  posterior <- checked$posterior
  c <- as.numeric(c)
  if (length(c) != 1L || !is.finite(c) || c < 1) {
    stop("The Occam's Window factor c must be finite and at least one.", call. = FALSE)
  }
  if (ncol(supports) > 30L) {
    stop("Occam's Window currently supports at most 30 variables.", call. = FALSE)
  }
  in_window <- posterior * c + .Machine$double.eps >= max(posterior)
  masks <- pbmc_support_masks(supports)
  table_size <- 2^ncol(supports)
  best_submodel <- rep(-Inf, table_size)
  best_submodel[masks[in_window] + 1L] <- posterior[in_window]
  if (ncol(supports) > 0L) {
    for (variable in seq_len(ncol(supports))) {
      bit <- bitwShiftL(1L, variable - 1L)
      containing <- which(bitwAnd(0:(table_size - 1L), bit) != 0L)
      without <- bitwXor(containing - 1L, bit) + 1L
      best_submodel[containing] <- pmax(
        best_submodel[containing], best_submodel[without]
      )
    }
  }
  has_better_submodel <- rep(FALSE, nrow(supports))
  candidates <- which(in_window)
  has_better_submodel[candidates] <-
    best_submodel[masks[candidates] + 1L] > posterior[candidates]
  selected <- in_window & !has_better_submodel
  list(
    method = "Occam's Window",
    selected = selected,
    selected_indices = which(selected),
    c = c,
    initially_in_window = sum(in_window),
    removed_by_submodel_rule = sum(in_window & has_better_submodel),
    description_entries = sum(selected),
    source_definition = "Madigan and Raftery (1994), equations (5)--(7)"
  )
}

pbmc_validate_partition <- function(blocks, variables) {
  blocks <- lapply(blocks, as.integer)
  variables <- as.integer(variables)
  if (length(blocks) < 1L || any(lengths(blocks) < 1L) ||
      !identical(sort(unlist(blocks, use.names = FALSE)), sort(variables)) ||
      anyDuplicated(unlist(blocks, use.names = FALSE))) {
    stop("Blocks must form a partition of the supplied variables.", call. = FALSE)
  }
  blocks
}

pbmc_block_marginal <- function(supports, posterior, block) {
  local <- supports[, block, drop = FALSE]
  code <- as.integer(local %*% 2^(seq_along(block) - 1L))
  probabilities <- numeric(2^length(block))
  collapsed <- rowsum(posterior, code, reorder = FALSE)
  probabilities[as.integer(rownames(collapsed)) + 1L] <- collapsed[, 1L]
  list(
    states = pbmc_enumerate_supports(length(block)),
    probabilities = probabilities,
    code = code
  )
}

pbmc_cartesian_set_for_partition <- function(supports,
                                             posterior,
                                             blocks,
                                             lambda) {
  checked <- pbmc_validate_posterior(supports, posterior)
  supports <- checked$supports
  posterior <- checked$posterior
  lambda <- as.numeric(lambda)
  blocks <- pbmc_validate_partition(blocks, seq_len(ncol(supports)))
  if (length(lambda) != 1L || !is.finite(lambda) || lambda <= 0 || lambda > 1) {
    stop("lambda must lie in (0,1].", call. = FALSE)
  }
  marginals <- lapply(blocks, function(block) {
    pbmc_block_marginal(supports, posterior, block)
  })
  retained <- lapply(marginals, function(x) seq_along(x$probabilities))
  retained_mass <- rep(1, length(blocks))
  repeat {
    removable <- lengths(retained) > 1L
    if (!any(removable)) break
    scores <- rep(Inf, length(blocks))
    for (block_index in which(removable)) {
      probabilities <- marginals[[block_index]]$probabilities[retained[[block_index]]]
      scores[block_index] <-
        min(probabilities) / sum(probabilities) * length(probabilities)
    }
    if (!any(is.finite(scores))) break
    chosen_block <- which.min(scores)
    probabilities <- marginals[[chosen_block]]$probabilities[retained[[chosen_block]]]
    remove_position <- which.min(probabilities)
    proposed_retained <- retained
    proposed_retained[[chosen_block]] <-
      proposed_retained[[chosen_block]][-remove_position]
    proposed_mass <- retained_mass
    proposed_mass[chosen_block] <- sum(
      marginals[[chosen_block]]$probabilities[proposed_retained[[chosen_block]]]
    )
    if (prod(proposed_mass) + 1e-14 < lambda) break
    retained <- proposed_retained
    retained_mass <- proposed_mass
  }
  selected <- rep(TRUE, nrow(supports))
  block_sets <- vector("list", length(blocks))
  for (block_index in seq_along(blocks)) {
    selected_codes <- retained[[block_index]] - 1L
    selected <- selected & marginals[[block_index]]$code %in% selected_codes
    block_sets[[block_index]] <- list(
      variables = blocks[[block_index]],
      states = marginals[[block_index]]$states[retained[[block_index]], , drop = FALSE],
      state_probabilities =
        marginals[[block_index]]$probabilities[retained[[block_index]]]
    )
  }
  state_counts <- lengths(retained)
  list(
    selected = selected,
    block_sets = block_sets,
    factorized_mass = prod(retained_mass),
    exact_mass = sum(posterior[selected]),
    cartesian_cardinality = prod(state_counts),
    description_entries = sum(state_counts),
    state_counts = state_counts
  )
}

pbmc_mutual_information <- function(supports, posterior, first, second) {
  first_marginal <- pbmc_block_marginal(supports, posterior, first)
  second_marginal <- pbmc_block_marginal(supports, posterior, second)
  joint <- c(first, second)
  joint_marginal <- pbmc_block_marginal(supports, posterior, joint)
  first_code <- joint_marginal$states[, seq_along(first), drop = FALSE] %*%
    2^(seq_along(first) - 1L)
  second_columns <- length(first) + seq_along(second)
  second_code <- joint_marginal$states[, second_columns, drop = FALSE] %*%
    2^(seq_along(second) - 1L)
  joint_probability <- joint_marginal$probabilities
  positive <- joint_probability > 0
  sum(joint_probability[positive] * log(
    joint_probability[positive] /
      (first_marginal$probabilities[first_code[positive] + 1L] *
         second_marginal$probabilities[second_code[positive] + 1L])
  ))
}

pbmc_merge_highest_information_blocks <- function(supports, posterior, blocks) {
  if (length(blocks) < 2L) return(blocks)
  pairs <- utils::combn(seq_along(blocks), 2L)
  information <- apply(pairs, 2L, function(pair) {
    pbmc_mutual_information(
      supports, posterior, blocks[[pair[1L]]], blocks[[pair[2L]]]
    )
  })
  chosen <- pairs[, which.max(information)]
  merged <- sort(c(blocks[[chosen[1L]]], blocks[[chosen[2L]]]))
  out <- blocks[-chosen]
  append(out, list(merged), after = chosen[1L] - 1L)
}

pbmc_cartesian_credible_set <- function(supports,
                                        posterior,
                                        epsilon,
                                        pip_threshold = 0.04,
                                        M = 2) {
  checked <- pbmc_validate_posterior(supports, posterior)
  supports <- checked$supports
  posterior <- checked$posterior
  epsilon <- as.numeric(epsilon)
  pip_threshold <- as.numeric(pip_threshold)
  M <- as.numeric(M)
  if (length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 1 ||
      length(pip_threshold) != 1L || !is.finite(pip_threshold) ||
      pip_threshold < 0 || pip_threshold >= 1 ||
      length(M) != 1L || !is.finite(M) || M <= 0) {
    stop("Invalid Cartesian credible-set tuning values.", call. = FALSE)
  }
  lambda <- 1 - epsilon
  pip <- colSums(supports * posterior)
  screened_variables <- which(pip > pip_threshold)
  excluded_variables <- setdiff(seq_len(ncol(supports)), screened_variables)
  screening_membership <- if (length(excluded_variables) == 0L) {
    rep(TRUE, nrow(supports))
  } else {
    rowSums(supports[, excluded_variables, drop = FALSE]) == 0L
  }
  screening_mass <- sum(posterior[screening_membership])
  if (length(screened_variables) == 0L) {
    selected <- screening_membership
    return(list(
      method = "Cartesian credible set",
      selected = selected,
      selected_indices = which(selected),
      screened_variables = integer(0),
      excluded_variables = excluded_variables,
      screening_mass = screening_mass,
      blocks = list(),
      block_sets = list(),
      factorized_mass = 1,
      exact_mass = sum(posterior[selected]),
      cartesian_cardinality = 1,
      description_entries = 1,
      published_criterion = 0,
      pip_threshold = pip_threshold,
      M = M,
      source_definition = "Griffin (2025), Algorithms 1--2"
    ))
  }
  reduced_supports <- supports[, screened_variables, drop = FALSE]
  blocks <- lapply(seq_along(screened_variables), function(index) index)
  candidates <- vector("list", length(blocks))
  for (stage in seq_along(candidates)) {
    candidate <- pbmc_cartesian_set_for_partition(
      reduced_supports, posterior, blocks, lambda
    )
    block_sizes <- lengths(blocks)
    candidate$criterion <- sum(log(candidate$state_counts)) +
      length(blocks) * log(M) + sum(lgamma(block_sizes))
    candidate$blocks <- lapply(blocks, function(block) screened_variables[block])
    for (block_index in seq_along(candidate$block_sets)) {
      candidate$block_sets[[block_index]]$variables <-
        screened_variables[candidate$block_sets[[block_index]]$variables]
    }
    candidate$selected <- candidate$selected & screening_membership
    candidate$exact_mass <- sum(posterior[candidate$selected])
    candidates[[stage]] <- candidate
    if (length(blocks) > 1L) {
      blocks <- pbmc_merge_highest_information_blocks(
        reduced_supports, posterior, blocks
      )
    }
  }
  chosen_index <- which.min(vapply(candidates, `[[`, numeric(1), "criterion"))
  chosen <- candidates[[chosen_index]]
  list(
    method = "Cartesian credible set",
    selected = chosen$selected,
    selected_indices = which(chosen$selected),
    screened_variables = screened_variables,
    excluded_variables = excluded_variables,
    screening_mass = screening_mass,
    blocks = chosen$blocks,
    block_sets = chosen$block_sets,
    factorized_mass = chosen$factorized_mass,
    exact_mass = chosen$exact_mass,
    cartesian_cardinality = chosen$cartesian_cardinality,
    description_entries = chosen$description_entries,
    published_criterion = chosen$criterion,
    partition_stage = chosen_index,
    pip_threshold = pip_threshold,
    M = M,
    source_definition = "Griffin (2025), Algorithms 1--2",
    criterion_convention =
      "author R code: sum(log block-set sizes) plus M=2 partition penalty"
  )
}

pbmc_comparison_table <- function(supports,
                                  posterior,
                                  epsilon = c(0.01, 0.05, 0.10, 0.20),
                                  occam_c = 20,
                                  ccs_pip_threshold = 0.04,
                                  ccs_M = 2) {
  checked <- pbmc_validate_posterior(supports, posterior)
  epsilon <- as.numeric(epsilon)
  if (length(epsilon) < 1L || any(!is.finite(epsilon)) ||
      any(epsilon < 0 | epsilon >= 1)) {
    stop("Every epsilon must lie in [0,1).", call. = FALSE)
  }
  occam <- pbmc_occams_window(
    checked$supports, checked$posterior, c = occam_c
  )
  rows <- vector("list", 3L * length(epsilon))
  cursor <- 0L
  for (value in epsilon) {
    hpp <- pbmc_hpp_credible_set(checked$supports, checked$posterior, value)
    ccs <- pbmc_cartesian_credible_set(
      checked$supports,
      checked$posterior,
      epsilon = value,
      pip_threshold = ccs_pip_threshold,
      M = ccs_M
    )
    methods <- list(hpp = hpp, occam = occam, ccs = ccs)
    labels <- c(
      hpp = "highest-posterior-probability set",
      occam = "Occam's Window",
      ccs = "Cartesian credible set"
    )
    native_cost <- c(
      hpp = "number of listed full models",
      occam = "number of listed full models",
      ccs = "number of listed block submodels"
    )
    for (name in names(methods)) {
      cursor <- cursor + 1L
      method <- methods[[name]]
      audit <- if (!is.null(method$evaluation)) {
        method$evaluation
      } else {
        pbmc_evaluate_set(
          checked$supports, checked$posterior, method$selected, value
        )
      }
      rows[[cursor]] <- data.frame(
        method = unname(labels[name]),
        epsilon = value,
        target_mass = audit$target_mass,
        retained_mass = audit$retained_mass,
        omitted_mass = audit$omitted_mass,
        target_met_under_exact_posterior = audit$feasible,
        full_model_count = audit$full_model_count,
        native_description_entries = method$description_entries,
        native_cost_definition = unname(native_cost[name]),
        maximum_pip_error = audit$maximum_pip_error,
        reverse_kullback_leibler = audit$reverse_kullback_leibler,
        factorized_target_mass = if (name == "ccs") {
          method$factorized_mass
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
