source(file.path("sim", "src", "region_restricted_smc.R"))
source(file.path("sim", "src", "exact_bayesian_network_helpers.R"))

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("The Bayesian-network study requires the yaml package.", call. = FALSE)
}

parse_arguments <- function() {
  raw <- commandArgs(trailingOnly = TRUE)
  out <- list()
  for (argument in raw) {
    if (!startsWith(argument, "--") || !grepl("=", argument, fixed = TRUE)) {
      stop("Arguments must have the form --name=value.", call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", argument), "=", fixed = TRUE)[[1L]]
    out[[pieces[1L]]] <- paste(pieces[-1L], collapse = "=")
  }
  out
}

split_values <- function(value) {
  if (is.null(value) || !nzchar(value)) NULL else
    strsplit(value, ",", fixed = TRUE)[[1L]]
}

arguments <- parse_arguments()
allowed_arguments <- c(
  "sample-size", "signal", "epsilon", "tree-order", "replications",
  "replication-ids", "output"
)
unknown_arguments <- setdiff(names(arguments), allowed_arguments)
if (length(unknown_arguments) > 0L) {
  stop(paste("Unknown arguments:", paste(unknown_arguments, collapse = ", ")),
       call. = FALSE)
}

protocol <- yaml::read_yaml(file.path(
  "sim", "config", "bayesian_network_protocol.yml"
))
if (!identical(protocol$status, "frozen_before_production")) {
  stop("The Bayesian-network protocol must be frozen before production.",
       call. = FALSE)
}

registered_sample_sizes <- as.integer(unlist(protocol$data_model$sample_sizes))
registered_signals <- as.numeric(unlist(protocol$data_model$signal_levels))
registered_epsilon <- as.numeric(unlist(protocol$selection$epsilon))
registered_tree_orders <- as.character(unlist(protocol$selection$tree_orders))
sample_sizes <- if (is.null(arguments$`sample-size`)) {
  registered_sample_sizes
} else {
  as.integer(split_values(arguments$`sample-size`))
}
signals <- if (is.null(arguments$signal)) {
  registered_signals
} else {
  as.numeric(split_values(arguments$signal))
}
epsilon_grid <- if (is.null(arguments$epsilon)) {
  registered_epsilon
} else {
  as.numeric(split_values(arguments$epsilon))
}
tree_orders <- if (is.null(arguments$`tree-order`)) {
  registered_tree_orders
} else {
  split_values(arguments$`tree-order`)
}
if (any(!sample_sizes %in% registered_sample_sizes) ||
    any(!signals %in% registered_signals) ||
    any(!epsilon_grid %in% registered_epsilon) ||
    any(!tree_orders %in% registered_tree_orders)) {
  stop("A requested cell is outside the frozen protocol.", call. = FALSE)
}

maximum_replications <- as.integer(
  protocol$data_model$production_replications_per_cell
)
if (!is.null(arguments$replications) && !is.null(arguments$`replication-ids`)) {
  stop("Use only one of --replications or --replication-ids.",
       call. = FALSE)
}
replications <- if (!is.null(arguments$`replication-ids`)) {
  as.integer(split_values(arguments$`replication-ids`))
} else if (!is.null(arguments$replications)) {
  seq_len(as.integer(arguments$replications))
} else {
  seq_len(maximum_replications)
}
if (any(!is.finite(replications)) || any(replications < 1L) ||
    any(replications > maximum_replications) || anyDuplicated(replications)) {
  stop("Invalid production replication request.", call. = FALSE)
}

default_output <- file.path(
  "sim", "output", "bayesian_network", "bayesian_network_exact.csv"
)
output_path <- if (is.null(arguments$output)) default_output else arguments$output
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

n_nodes <- as.integer(protocol$state_space$variables)
dag_space <- bn_enumerate_dags(n_nodes)
true_state <- bn_true_five_node_state(dag_space$pairs)
equivalent_sample_size <- as.numeric(
  protocol$posterior$equivalent_sample_size
)
edge_probability <- as.numeric(protocol$posterior$edge_probability)
score_budget <- as.numeric(protocol$production_settings$score_budget)
action_budget <- as.numeric(protocol$production_settings$action_budget)

data_seed_for <- function(sample_size, signal, replication) {
  as.integer(protocol$seed_registry$production_base) +
    100000L * match(sample_size, registered_sample_sizes) +
    10000L * match(signal, registered_signals) +
    as.integer(replication)
}

algorithm_seed_for <- function(data_seed, tree_order, epsilon) {
  data_seed + 100L * match(tree_order, registered_tree_orders) +
    10L * match(epsilon, registered_epsilon)
}

bind_output_rows <- function(rows) {
  all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))
  rows <- lapply(rows, function(row) {
    missing <- setdiff(all_names, names(row))
    for (name in missing) row[[name]] <- NA
    row[all_names]
  })
  do.call(rbind, rows)
}

error_rows <- function(sample_size, signal, replication, error) {
  rows <- expand.grid(
    epsilon = epsilon_grid,
    tree_order = tree_orders,
    stringsAsFactors = FALSE
  )
  rows$protocol_id <- protocol$protocol_id
  rows$run_mode <- "production"
  rows$registered_production <- TRUE
  rows$n_nodes <- n_nodes
  rows$n_dags <- nrow(dag_space$states)
  rows$sample_size <- sample_size
  rows$signal <- signal
  rows$replication <- replication
  rows$interval_mode <- "exact"
  rows$score_oracle <- "precomputed_exact_DAG_table"
  rows$evidence_method <- "exact"
  rows$data_seed <- data_seed_for(sample_size, signal, replication)
  rows$algorithm_seed <- mapply(
    algorithm_seed_for,
    data_seed = rows$data_seed,
    tree_order = rows$tree_order,
    epsilon = rows$epsilon
  )
  rows$score_budget <- score_budget
  rows$action_budget <- action_budget
  rows$raw_status <- "implementation_error"
  rows$status <- "implementation-error"
  rows$error_message <- conditionMessage(error)
  rows
}

rows <- list()
index <- 0L
for (sample_size in sample_sizes) {
  for (signal in signals) {
    for (replication in replications) {
      data_seed <- data_seed_for(sample_size, signal, replication)
      dataset_rows <- tryCatch(
        bn_run_exact_skeleton_dataset(
          protocol_id = protocol$protocol_id,
          run_mode = "production",
          registered_production = TRUE,
          dag_space = dag_space,
          true_state = true_state,
          sample_size = sample_size,
          signal = signal,
          replication = replication,
          epsilon_grid = epsilon_grid,
          tree_orders = tree_orders,
          equivalent_sample_size = equivalent_sample_size,
          edge_probability = edge_probability,
          data_seed = data_seed,
          algorithm_seed_for = algorithm_seed_for,
          score_budget = score_budget,
          action_budget = action_budget
        ),
        error = function(error) error_rows(
          sample_size, signal, replication, error
        )
      )
      for (row in seq_len(nrow(dataset_rows))) {
        index <- index + 1L
        rows[[index]] <- dataset_rows[row, , drop = FALSE]
      }
      output <- bind_output_rows(rows)
      write.csv(output, paste0(output_path, ".partial"), row.names = FALSE,
                na = "")
      message(sprintf(
        "sample_size=%d signal=%.1f replication=%d rows=%d statuses=%s",
        sample_size, signal, replication, nrow(output),
        paste(unique(dataset_rows$status), collapse = ",")
      ))
    }
  }
}

output <- bind_output_rows(rows)
write.csv(output, output_path, row.names = FALSE, na = "")
partial_path <- paste0(output_path, ".partial")
if (file.exists(partial_path)) unlink(partial_path)
cat(sprintf(
  "Wrote %d registered Bayesian-network rows to %s.\n",
  nrow(output), normalizePath(output_path, winslash = "/", mustWork = FALSE)
))
