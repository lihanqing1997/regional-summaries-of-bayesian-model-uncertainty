if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("The Bayesian-network summary requires the yaml package.", call. = FALSE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) > 3L) {
  stop("Use at most: input.csv summary.csv audit.csv", call. = FALSE)
}
input_path <- if (length(arguments) >= 1L) arguments[[1L]] else file.path(
  "sim", "output", "bayesian_network", "bayesian_network_exact.csv"
)
summary_path <- if (length(arguments) >= 2L) arguments[[2L]] else file.path(
  "sim", "output", "bayesian_network", "bayesian_network_summary.csv"
)
audit_path <- if (length(arguments) >= 3L) arguments[[3L]] else file.path(
  "sim", "output", "bayesian_network", "bayesian_network_audit.csv"
)
protocol <- yaml::read_yaml(file.path(
  "sim", "config", "bayesian_network_protocol.yml"
))
data <- read.csv(input_path, stringsAsFactors = FALSE)

group_key <- interaction(
  data$sample_size,
  data$signal,
  data$epsilon,
  data$tree_order,
  drop = TRUE,
  lex.order = TRUE
)
summary_rows <- lapply(split(data, group_key), function(cell) {
  data.frame(
    sample_size = cell$sample_size[1L],
    signal = cell$signal[1L],
    epsilon = cell$epsilon[1L],
    tree_order = cell$tree_order[1L],
    registered_runs = nrow(cell),
    bounds_met_rate = mean(cell$status == "bounds-met"),
    exact_cost_optimal_rate = mean(cell$exact_cost_optimal),
    epsilon_feasible_rate = mean(cell$epsilon_feasible),
    median_selected_cost = stats::median(cell$selected_cost),
    median_ranked_dag_cost = stats::median(cell$ranked_dag_cost),
    mean_cost_ratio_vs_ranked_dags = mean(cell$cost_ratio_vs_ranked_dags),
    median_cost_ratio_vs_ranked_dags =
      stats::median(cell$cost_ratio_vs_ranked_dags),
    q10_cost_ratio_vs_ranked_dags = stats::quantile(
      cell$cost_ratio_vs_ranked_dags, 0.10, names = FALSE
    ),
    q90_cost_ratio_vs_ranked_dags = stats::quantile(
      cell$cost_ratio_vs_ranked_dags, 0.90, names = FALSE
    ),
    mean_cost_ratio_vs_ranked_skeletons =
      mean(cell$cost_ratio_vs_ranked_skeletons),
    mean_true_omitted_mass = mean(cell$true_omitted_mass),
    mean_maximum_skeleton_edge_error =
      mean(cell$maximum_skeleton_edge_error),
    q90_maximum_skeleton_edge_error = stats::quantile(
      cell$maximum_skeleton_edge_error, 0.90, names = FALSE
    ),
    mean_selected_skeletons = mean(cell$selected_skeletons),
    mean_search_actions = mean(cell$search_actions),
    mean_raw_score_requests = mean(cell$raw_score_requests),
    mean_elapsed_seconds = mean(cell$elapsed_seconds),
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
summary <- summary[order(
  summary$sample_size,
  summary$signal,
  summary$epsilon,
  summary$tree_order
), , drop = FALSE]

expected_rows <-
  length(unlist(protocol$data_model$sample_sizes)) *
  length(unlist(protocol$data_model$signal_levels)) *
  as.integer(protocol$data_model$production_replications_per_cell) *
  length(unlist(protocol$selection$epsilon)) *
  length(unlist(protocol$selection$tree_orders))
key <- paste(
  data$sample_size, data$signal, data$replication,
  data$epsilon, data$tree_order, sep = "/"
)
matched_key <- paste(
  data$sample_size, data$signal, data$replication,
  data$epsilon, sep = "/"
)
tree_cost_agreement <- all(vapply(
  split(data$selected_cost, matched_key),
  function(values) length(unique(values)) == 1L,
  logical(1)
))
cell_medians <- aggregate(
  cost_ratio_vs_ranked_dags ~ sample_size + signal + epsilon + tree_order,
  data,
  stats::median
)
checks <- data.frame(
  check = c(
    "expected_registered_rows",
    "unique_registered_keys",
    "all_bounds_met",
    "all_exact_cost_optimal",
    "all_cost_sandwich_covered",
    "all_epsilon_feasible",
    "tree_order_cost_agreement",
    "cell_median_cost_ratio_below_half",
    "no_implementation_error"
  ),
  passed = c(
    nrow(data) == expected_rows,
    !anyDuplicated(key),
    all(data$status == "bounds-met"),
    all(data$exact_cost_optimal),
    all(data$cost_sandwich_covered),
    all(data$epsilon_feasible),
    tree_cost_agreement,
    all(cell_medians$cost_ratio_vs_ranked_dags < 0.50),
    !any(data$status == "implementation-error")
  ),
  stringsAsFactors = FALSE
)
checks$value <- c(
  nrow(data),
  length(unique(key)),
  mean(data$status == "bounds-met"),
  mean(data$exact_cost_optimal),
  mean(data$cost_sandwich_covered),
  mean(data$epsilon_feasible),
  as.numeric(tree_cost_agreement),
  max(cell_medians$cost_ratio_vs_ranked_dags),
  sum(data$status == "implementation-error")
)
checks$required <- c(
  expected_rows, expected_rows, 1, 1, 1, 1, 1, 0.50, 0
)

dir.create(dirname(summary_path), recursive = TRUE, showWarnings = FALSE)
write.csv(summary, summary_path, row.names = FALSE, na = "")
write.csv(checks, audit_path, row.names = FALSE, na = "")
if (!all(checks$passed)) {
  print(checks[!checks$passed, , drop = FALSE])
  stop("The Bayesian-network production audit failed.", call. = FALSE)
}
cat(sprintf(
  "Bayesian-network audit passed: %d rows and %d checks.\n",
  nrow(data), nrow(checks)
))
