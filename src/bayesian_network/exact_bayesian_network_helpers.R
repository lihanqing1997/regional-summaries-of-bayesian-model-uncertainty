bn_pair_table <- function(n_nodes) {
  n_nodes <- as.integer(n_nodes)
  if (length(n_nodes) != 1L || !is.finite(n_nodes) || n_nodes < 2L) {
    stop("n_nodes must be one integer of at least two.", call. = FALSE)
  }
  pairs <- t(utils::combn(seq_len(n_nodes), 2L))
  colnames(pairs) <- c("lower", "upper")
  pairs
}

bn_state_ids <- function(states) {
  states <- as.matrix(states)
  if (ncol(states) < 1L || any(!states %in% 0:2)) {
    stop("Bayesian-network states must be ternary matrices.", call. = FALSE)
  }
  as.integer(1 + states %*% 3^(seq_len(ncol(states)) - 1L))
}

bn_state_to_adjacency <- function(state, pairs, n_nodes) {
  state <- as.integer(state)
  pairs <- as.matrix(pairs)
  if (length(state) != nrow(pairs) || any(!state %in% 0:2)) {
    stop("The ternary state does not match the node-pair table.", call. = FALSE)
  }
  adjacency <- matrix(FALSE, nrow = n_nodes, ncol = n_nodes)
  forward <- which(state == 1L)
  reverse <- which(state == 2L)
  if (length(forward) > 0L) {
    adjacency[pairs[forward, , drop = FALSE]] <- TRUE
  }
  if (length(reverse) > 0L) {
    adjacency[pairs[reverse, 2:1, drop = FALSE]] <- TRUE
  }
  adjacency
}

bn_topological_order <- function(adjacency) {
  adjacency <- as.matrix(adjacency) != 0
  if (nrow(adjacency) != ncol(adjacency) || any(diag(adjacency))) {
    stop("adjacency must be a square matrix with an empty diagonal.",
         call. = FALSE)
  }
  n_nodes <- nrow(adjacency)
  indegree <- colSums(adjacency)
  available <- which(indegree == 0L)
  order <- integer(0)
  while (length(available) > 0L) {
    node <- available[1L]
    available <- available[-1L]
    order <- c(order, node)
    children <- which(adjacency[node, ])
    if (length(children) > 0L) {
      indegree[children] <- indegree[children] - 1L
      newly_available <- children[indegree[children] == 0L]
      available <- sort(unique(c(available, newly_available)))
    }
  }
  if (length(order) == n_nodes) order else integer(0)
}

bn_is_acyclic_state <- function(state, pairs, n_nodes) {
  length(bn_topological_order(
    bn_state_to_adjacency(state, pairs, n_nodes)
  )) == n_nodes
}

bn_enumerate_dags <- function(n_nodes) {
  pairs <- bn_pair_table(n_nodes)
  n_pairs <- nrow(pairs)
  ids <- 0:(3^n_pairs - 1L)
  states <- vapply(seq_len(n_pairs), function(index) {
    as.integer((ids %/% 3^(index - 1L)) %% 3L)
  }, integer(length(ids)))
  if (is.null(dim(states))) states <- matrix(states, ncol = 1L)
  acyclic <- vapply(seq_len(nrow(states)), function(index) {
    bn_is_acyclic_state(states[index, ], pairs, n_nodes)
  }, logical(1))
  states <- states[acyclic, , drop = FALSE]
  colnames(states) <- apply(pairs, 1L, function(pair) {
    paste0(pair[1L], "-", pair[2L])
  })
  list(
    n_nodes = as.integer(n_nodes),
    pairs = pairs,
    states = states,
    state_ids = bn_state_ids(states)
  )
}

bn_true_five_node_state <- function(pairs = bn_pair_table(5L)) {
  directed_edges <- matrix(c(
    1L, 3L,
    2L, 3L,
    2L, 4L,
    3L, 5L,
    4L, 5L
  ), ncol = 2L, byrow = TRUE)
  state <- integer(nrow(pairs))
  for (index in seq_len(nrow(directed_edges))) {
    from <- directed_edges[index, 1L]
    to <- directed_edges[index, 2L]
    pair_index <- which(
      pairs[, 1L] == min(from, to) & pairs[, 2L] == max(from, to)
    )
    state[pair_index] <- if (from < to) 1L else 2L
  }
  state
}

bn_simulate_binary <- function(n,
                               state,
                               pairs,
                               signal = 1,
                               seed = NULL,
                               intercept = 0) {
  n <- as.integer(n)
  n_nodes <- max(pairs)
  if (length(n) != 1L || n < 1L || !is.finite(signal) || signal < 0) {
    stop("Invalid Bayesian-network simulation inputs.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))
  adjacency <- bn_state_to_adjacency(state, pairs, n_nodes)
  order <- bn_topological_order(adjacency)
  if (length(order) != n_nodes) {
    stop("The data-generating graph must be acyclic.", call. = FALSE)
  }
  intercept <- rep_len(as.numeric(intercept), n_nodes)
  data <- matrix(0L, nrow = n, ncol = n_nodes)
  for (node in order) {
    parents <- which(adjacency[, node])
    linear_predictor <- rep(intercept[node], n)
    if (length(parents) > 0L) {
      centered_parents <- 2 * data[, parents, drop = FALSE] - 1
      linear_predictor <- linear_predictor +
        signal * rowSums(centered_parents) / sqrt(length(parents))
    }
    data[, node] <- stats::rbinom(n, size = 1L,
                                  prob = stats::plogis(linear_predictor))
  }
  colnames(data) <- paste0("X", seq_len(n_nodes))
  data
}

bn_bdeu_local_scores <- function(data, equivalent_sample_size = 1) {
  data <- as.matrix(data)
  storage.mode(data) <- "integer"
  n_nodes <- ncol(data)
  if (nrow(data) < 1L || n_nodes < 2L || any(!data %in% 0:1) ||
      length(equivalent_sample_size) != 1L ||
      !is.finite(equivalent_sample_size) || equivalent_sample_size <= 0) {
    stop("BDeu scoring needs a nonempty binary data matrix and positive ESS.",
         call. = FALSE)
  }
  scores <- matrix(NA_real_, nrow = n_nodes, ncol = 2^n_nodes)
  for (node in seq_len(n_nodes)) {
    for (mask in 0:(2^n_nodes - 1L)) {
      if (bitwAnd(mask, bitwShiftL(1L, node - 1L)) != 0L) next
      parents <- which(vapply(seq_len(n_nodes), function(parent) {
        bitwAnd(mask, bitwShiftL(1L, parent - 1L)) != 0L
      }, logical(1)))
      n_configurations <- 2^length(parents)
      parent_configuration <- if (length(parents) == 0L) {
        integer(nrow(data))
      } else {
        as.integer(data[, parents, drop = FALSE] %*%
                     2^(seq_along(parents) - 1L))
      }
      alpha_parent <- equivalent_sample_size / n_configurations
      alpha_cell <- equivalent_sample_size / (2 * n_configurations)
      local <- 0
      for (configuration in 0:(n_configurations - 1L)) {
        rows <- parent_configuration == configuration
        counts <- tabulate(data[rows, node] + 1L, nbins = 2L)
        total <- sum(counts)
        local <- local + lgamma(alpha_parent) -
          lgamma(alpha_parent + total) +
          sum(lgamma(alpha_cell + counts) - lgamma(alpha_cell))
      }
      scores[node, mask + 1L] <- local
    }
  }
  scores
}

bn_parent_masks <- function(states, pairs, n_nodes) {
  states <- as.matrix(states)
  masks <- matrix(0L, nrow = nrow(states), ncol = n_nodes)
  for (pair_index in seq_len(nrow(pairs))) {
    lower <- pairs[pair_index, 1L]
    upper <- pairs[pair_index, 2L]
    forward <- states[, pair_index] == 1L
    reverse <- states[, pair_index] == 2L
    masks[forward, upper] <- masks[forward, upper] +
      bitwShiftL(1L, lower - 1L)
    masks[reverse, lower] <- masks[reverse, lower] +
      bitwShiftL(1L, upper - 1L)
  }
  masks
}

bn_score_dags <- function(states,
                          pairs,
                          local_scores,
                          edge_probability = 0.2) {
  states <- as.matrix(states)
  n_nodes <- nrow(local_scores)
  n_pairs <- nrow(pairs)
  if (ncol(states) != n_pairs ||
      ncol(local_scores) != 2^n_nodes ||
      length(edge_probability) != 1L ||
      !is.finite(edge_probability) ||
      edge_probability <= 0 || edge_probability >= 1) {
    stop("Invalid DAG scoring inputs.", call. = FALSE)
  }
  parent_masks <- bn_parent_masks(states, pairs, n_nodes)
  log_likelihood <- numeric(nrow(states))
  for (node in seq_len(n_nodes)) {
    log_likelihood <- log_likelihood +
      local_scores[cbind(node, parent_masks[, node] + 1L)]
  }
  edge_count <- rowSums(states != 0L)
  log_prior <- edge_count * log(edge_probability / 2) +
    (n_pairs - edge_count) * log1p(-edge_probability)
  list(
    log_likelihood = log_likelihood,
    log_prior = log_prior,
    log_score = log_likelihood + log_prior,
    edge_count = edge_count,
    parent_masks = parent_masks
  )
}

bn_make_precomputed_state_score <- function(log_scores, n_pairs) {
  log_scores <- as.numeric(log_scores)
  lookup <- rep(NA_real_, 3^n_pairs)
  seen <- rep(FALSE, length(lookup))
  calls <- 0L
  evaluations <- 0L
  matrix_score <- function(states) {
    states <- as.matrix(states)
    if (ncol(states) != n_pairs) {
      stop("The score oracle received the wrong state dimension.",
           call. = FALSE)
    }
    ids <- bn_state_ids(states)
    values <- lookup[ids]
    if (anyNA(values)) {
      stop("The score oracle received a cyclic or undeclared state.",
           call. = FALSE)
    }
    calls <<- calls + length(ids)
    new_ids <- unique(ids[!seen[ids]])
    if (length(new_ids) > 0L) {
      seen[new_ids] <<- TRUE
      evaluations <<- evaluations + length(new_ids)
    }
    values
  }
  score <- function(state) {
    matrix_score(matrix(state, nrow = 1L))[1L]
  }
  attr(score, "matrix_score") <- matrix_score
  initialize <- function(states) {
    ids <- bn_state_ids(states)
    if (length(ids) != length(log_scores) || anyDuplicated(ids)) {
      stop("Precomputed scores do not match the declared DAG states.",
           call. = FALSE)
    }
    lookup[ids] <<- log_scores
    invisible(NULL)
  }
  reset <- function() {
    seen[] <<- FALSE
    calls <<- 0L
    evaluations <<- 0L
    invisible(NULL)
  }
  diagnostics <- function() c(calls = calls, evaluations = evaluations)
  structure(
    list(score = score, initialize = initialize, reset = reset,
         diagnostics = diagnostics),
    class = "bn_precomputed_state_score"
  )
}

bn_skeleton_keys <- function(states) {
  states <- as.matrix(states)
  apply(1L * (states != 0L), 1L, paste0, collapse = "")
}

bn_skeleton_reference <- function(states,
                                  log_scores,
                                  base_cost = 1L,
                                  edge_cost = 1L) {
  log_scores <- as.numeric(log_scores)
  log_normalizer <- rrs_log_sum_exp(log_scores)
  posterior <- exp(log_scores - log_normalizer)
  keys <- bn_skeleton_keys(states)
  masses <- rowsum(posterior, keys, reorder = TRUE)[, 1L]
  patterns <- do.call(rbind, strsplit(names(masses), "", fixed = TRUE))
  storage.mode(patterns) <- "integer"
  costs <- as.integer(base_cost + edge_cost * rowSums(patterns))
  list(
    keys = names(masses),
    mass = as.numeric(masses),
    log_mass = log(as.numeric(masses)),
    costs = costs,
    patterns = patterns,
    state_keys = keys,
    posterior = posterior,
    log_normalizer = log_normalizer
  )
}

bn_ranked_cover <- function(mass, costs, epsilon) {
  mass <- as.numeric(mass)
  costs <- as.numeric(costs)
  if (length(mass) != length(costs) || any(mass < 0) ||
      abs(sum(mass) - 1) > 1e-8 || epsilon < 0 || epsilon >= 1) {
    stop("Invalid ranked-cover inputs.", call. = FALSE)
  }
  order <- order(mass, decreasing = TRUE)
  count <- which(cumsum(mass[order]) >= 1 - epsilon)[1L]
  selected <- order[seq_len(count)]
  list(
    selected = selected,
    count = count,
    cost = sum(costs[selected]),
    retained_mass = sum(mass[selected])
  )
}

new_skeleton_region_tree <- function(states,
                                     log_scores,
                                     log_prior,
                                     split_order = seq_len(ncol(states)),
                                     base_cost = 1L,
                                     edge_cost = 1L,
                                     exact_evidence_bounds = TRUE) {
  states <- as.matrix(states)
  log_scores <- as.numeric(log_scores)
  log_prior <- as.numeric(log_prior)
  n_pairs <- ncol(states)
  split_order <- as.integer(split_order)
  if (nrow(states) < 1L || length(log_scores) != nrow(states) ||
      length(log_prior) != nrow(states) || any(!is.finite(log_scores)) ||
      any(!is.finite(log_prior)) || length(split_order) != n_pairs ||
      !setequal(split_order, seq_len(n_pairs)) ||
      any(c(base_cost, edge_cost) < 1) ||
      any(abs(c(base_cost, edge_cost) - round(c(base_cost, edge_cost))) >
          sqrt(.Machine$double.eps))) {
    stop("Invalid skeleton-region tree inputs.", call. = FALSE)
  }
  base_cost <- as.integer(round(base_cost))
  edge_cost <- as.integer(round(edge_cost))
  skeletons <- states != 0L
  cache <- new.env(hash = TRUE, parent = emptyenv())

  summary_for <- function(constraints) {
    key <- paste(constraints, collapse = "")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    fixed <- which(constraints >= 0L)
    inside <- if (length(fixed) == 0L) {
      rep(TRUE, nrow(states))
    } else {
      rowSums(sweep(
        skeletons[, fixed, drop = FALSE],
        2L,
        constraints[fixed] == 1L,
        `!=`
      )) == 0L
    }
    indices <- which(inside)
    if (length(indices) == 0L) {
      stop("A declared partial skeleton contains no DAG.", call. = FALSE)
    }
    log_proposal_normalizer <- rrs_log_sum_exp(log_prior[indices])
    log_proposal <- log_prior[indices] - log_proposal_normalizer
    log_g <- log_scores[indices] - log_proposal
    value <- list(
      indices = indices,
      proposal_probability = exp(log_proposal),
      log_g_bounds = range(log_g),
      rho = diff(range(log_g)),
      log_evidence = rrs_log_sum_exp(log_scores[indices])
    )
    assign(key, value, envir = cache)
    value
  }

  node_from_constraints <- function(constraints, depth, id) {
    summary <- summary_for(constraints)
    terminal <- depth >= n_pairs
    cost_lower <- base_cost + edge_cost * sum(constraints == 1L)
    region <- new_enumerated_state_region(
      states = states[summary$indices, , drop = FALSE],
      proposal_prob = summary$proposal_probability,
      cost = cost_lower,
      label = paste0("depth-", depth, ":", paste(constraints, collapse = ""))
    )
    region$constraints <- constraints
    anchor_state <- states[summary$indices[which.max(
      summary$proposal_probability
    )], , drop = FALSE]
    region$anchor <- local({
      fixed_anchor <- anchor_state
      function() fixed_anchor
    })
    region$log_g_bounds <- local({
      bounds <- summary$log_g_bounds
      function() bounds
    })
    if (isTRUE(exact_evidence_bounds)) {
      region$log_evidence_bounds <- local({
        value <- summary$log_evidence
        function() c(value, value)
      })
    }
    new_region_tree_node(
      id = id,
      region = region,
      children = character(0),
      rho = summary$rho,
      cost_lower = cost_lower,
      terminal = terminal
    )
  }

  expand_node <- function(node) {
    constraints <- node$region$constraints
    depth <- sum(constraints >= 0L)
    if (depth >= n_pairs) {
      stop("A complete skeleton node cannot be expanded.", call. = FALSE)
    }
    pair_index <- split_order[depth + 1L]
    lapply(0:1, function(value) {
      child_constraints <- constraints
      child_constraints[pair_index] <- value
      child_id <- paste0(node$id, "/e", pair_index, "-", value)
      node_from_constraints(child_constraints, depth + 1L, child_id)
    })
  }

  root <- node_from_constraints(rep(-1L, n_pairs), 0L, "root")
  tree <- new_lazy_region_tree(
    root_node = root,
    n_declared_nodes = 2^(n_pairs + 1L) - 1L,
    expand_node = expand_node
  )
  tree$summary_for <- summary_for
  tree$split_order <- split_order
  tree
}

bn_run_exact_skeleton_dataset <- function(protocol_id,
                                          run_mode,
                                          registered_production,
                                          dag_space,
                                          true_state,
                                          sample_size,
                                          signal,
                                          replication,
                                          epsilon_grid,
                                          tree_orders,
                                          equivalent_sample_size,
                                          edge_probability,
                                          data_seed,
                                          algorithm_seed_for,
                                          score_budget = 10000,
                                          action_budget = 10000) {
  data <- bn_simulate_binary(
    n = sample_size,
    state = true_state,
    pairs = dag_space$pairs,
    signal = signal,
    seed = data_seed
  )
  local_scores <- bn_bdeu_local_scores(
    data,
    equivalent_sample_size = equivalent_sample_size
  )
  scored <- bn_score_dags(
    states = dag_space$states,
    pairs = dag_space$pairs,
    local_scores = local_scores,
    edge_probability = edge_probability
  )
  reference <- bn_skeleton_reference(
    dag_space$states,
    scored$log_score,
    base_cost = 1L,
    edge_cost = 1L
  )
  true_skeleton_key <- bn_skeleton_keys(matrix(true_state, nrow = 1L))
  true_skeleton_mass <- reference$mass[
    match(true_skeleton_key, reference$keys)
  ]
  full_skeleton_edge_probability <- colSums(
    (dag_space$states != 0L) * reference$posterior
  )
  state_cost <- 1L + scored$edge_count
  score <- bn_make_precomputed_state_score(
    scored$log_score,
    n_pairs = nrow(dag_space$pairs)
  )
  score$initialize(dag_space$states)

  split_order_for <- function(tree_order) {
    if (tree_order == "lexicographic") {
      seq_len(nrow(dag_space$pairs))
    } else if (tree_order == "reverse") {
      rev(seq_len(nrow(dag_space$pairs)))
    } else {
      stop("Unknown tree order.", call. = FALSE)
    }
  }

  rows <- list()
  row_index <- 0L
  for (epsilon in epsilon_grid) {
    ranked_dag <- bn_ranked_cover(
      reference$posterior,
      state_cost,
      epsilon
    )
    ranked_skeleton <- bn_ranked_cover(
      reference$mass,
      reference$costs,
      epsilon
    )
    exact_optimum <- select_region_union_epsilon(
      log_evidence_lower = reference$log_mass,
      log_evidence_upper = reference$log_mass,
      costs = reference$costs,
      epsilon = epsilon,
      mode = "conservative"
    )
    for (tree_order in tree_orders) {
      algorithm_seed <- algorithm_seed_for(data_seed, tree_order, epsilon)
      set.seed(algorithm_seed)
      tree <- new_skeleton_region_tree(
        states = dag_space$states,
        log_scores = scored$log_score,
        log_prior = scored$log_prior,
        split_order = split_order_for(tree_order),
        base_cost = 1L,
        edge_cost = 1L,
        exact_evidence_bounds = TRUE
      )
      root_rho <- tree$nodes[[tree$root_id]]$rho
      score$reset()
      elapsed <- system.time({
        search <- select_minimum_cost_region_tree(
          tree = tree,
          log_score_fn = score$score,
          evidence_evaluator = function(node, remaining_score_budget = Inf) {
            stop("An exact terminal was unexpectedly reevaluated.",
                 call. = FALSE)
          },
          epsilon = epsilon,
          cost_tolerance = 0,
          max_score_evaluations = score_budget,
          max_actions = action_budget
        )
      })[["elapsed"]]
      selected_nodes <- search$nodes[search$selected_ids]
      selected_keys <- vapply(selected_nodes, function(node) {
        paste(node$region$constraints, collapse = "")
      }, character(1))
      selected_reference_rows <- match(selected_keys, reference$keys)
      selected_cost <- if (length(selected_nodes) == 0L) {
        Inf
      } else {
        sum(vapply(selected_nodes, function(node) {
          node$region$cost
        }, numeric(1)))
      }
      retained_mass <- if (length(selected_reference_rows) == 0L) {
        0
      } else {
        sum(reference$mass[selected_reference_rows])
      }
      selected_states <- reference$state_keys %in% selected_keys
      maximum_skeleton_edge_error <- if (retained_mass > 0) {
        restricted_posterior <-
          reference$posterior * selected_states / retained_mass
        max(abs(
          colSums((dag_space$states != 0L) * restricted_posterior) -
            full_skeleton_edge_probability
        ))
      } else {
        NA_real_
      }
      diagnostics <- score$diagnostics()
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        protocol_id = protocol_id,
        run_mode = run_mode,
        registered_production = registered_production,
        n_nodes = dag_space$n_nodes,
        n_dags = nrow(dag_space$states),
        sample_size = sample_size,
        signal = signal,
        replication = replication,
        epsilon = epsilon,
        tree_order = tree_order,
        interval_mode = "exact",
        score_oracle = "precomputed_exact_DAG_table",
        evidence_method = "exact",
        data_seed = data_seed,
        algorithm_seed = algorithm_seed,
        score_budget = score_budget,
        action_budget = action_budget,
        raw_status = search$status,
        status = search$registered_status,
        exact_optimal_cost = exact_optimum$budget,
        cost_lower = search$cost_lower,
        cost_upper = search$cost_upper,
        cost_gap = search$cost_gap,
        selected_cost = selected_cost,
        exact_cost_optimal = is.finite(selected_cost) &&
          selected_cost == exact_optimum$budget,
        cost_sandwich_covered =
          search$cost_lower <= exact_optimum$budget &&
          exact_optimum$budget <= search$cost_upper,
        selected_skeletons = length(selected_keys),
        retained_mass = retained_mass,
        true_omitted_mass = 1 - retained_mass,
        epsilon_feasible = 1 - retained_mass <= epsilon + 1e-12,
        maximum_skeleton_edge_error = maximum_skeleton_edge_error,
        ranked_dag_count = ranked_dag$count,
        ranked_dag_cost = ranked_dag$cost,
        ranked_dag_retained_mass = ranked_dag$retained_mass,
        ranked_skeleton_count = ranked_skeleton$count,
        ranked_skeleton_cost = ranked_skeleton$cost,
        ranked_skeleton_retained_mass = ranked_skeleton$retained_mass,
        cost_ratio_vs_ranked_dags = selected_cost / ranked_dag$cost,
        cost_ratio_vs_ranked_skeletons =
          selected_cost / ranked_skeleton$cost,
        true_skeleton_posterior_mass = true_skeleton_mass,
        effective_dag_count = 1 / sum(reference$posterior^2),
        effective_skeleton_count = 1 / sum(reference$mass^2),
        root_log_weight_range = root_rho,
        search_actions = search$actions,
        search_expansions = search$expansions,
        frontier_size = length(search$frontier_ids),
        raw_score_requests = unname(diagnostics[["calls"]]),
        distinct_score_evaluations =
          unname(diagnostics[["evaluations"]]),
        elapsed_seconds = elapsed,
        error_message = NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}
