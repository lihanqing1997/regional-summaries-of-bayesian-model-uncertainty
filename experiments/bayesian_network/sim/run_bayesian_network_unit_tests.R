source(file.path("sim", "src", "region_restricted_smc.R"))
source(file.path("sim", "src", "exact_bayesian_network_helpers.R"))

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

dag3 <- bn_enumerate_dags(3L)
dag4 <- bn_enumerate_dags(4L)
dag5 <- bn_enumerate_dags(5L)
expect_true(nrow(dag3$states) == 25L, "Three nodes should have 25 DAGs.")
expect_true(nrow(dag4$states) == 543L, "Four nodes should have 543 DAGs.")
expect_true(nrow(dag5$states) == 29281L, "Five nodes should have 29,281 DAGs.")
expect_true(
  all(vapply(seq_len(nrow(dag4$states)), function(index) {
    bn_is_acyclic_state(dag4$states[index, ], dag4$pairs, 4L)
  }, logical(1))),
  "The enumerator returned a cyclic graph."
)

set.seed(4201)
equivalence_data <- matrix(
  stats::rbinom(600L, size = 1L, prob = 0.5),
  ncol = 3L
)
local_scores <- bn_bdeu_local_scores(equivalence_data, equivalent_sample_size = 1)
chain <- matrix(c(1L, 0L, 1L), nrow = 1L)
fork <- matrix(c(2L, 0L, 1L), nrow = 1L)
equivalent_scores <- bn_score_dags(
  rbind(chain, fork),
  dag3$pairs,
  local_scores,
  edge_probability = 0.2
)$log_score
expect_true(
  abs(diff(equivalent_scores)) < 1e-10,
  "BDeu did not assign equal scores to equivalent chain and fork DAGs."
)

true_state <- c(1L, 0L, 1L)
data <- bn_simulate_binary(
  n = 120L,
  state = true_state,
  pairs = dag3$pairs,
  signal = 1,
  seed = 901L
)
local_scores <- bn_bdeu_local_scores(data, equivalent_sample_size = 1)
scored <- bn_score_dags(
  dag3$states,
  dag3$pairs,
  local_scores,
  edge_probability = 0.2
)
score <- bn_make_precomputed_state_score(scored$log_score, nrow(dag3$pairs))
score$initialize(dag3$states)
reference <- bn_skeleton_reference(dag3$states, scored$log_score)
epsilon <- 0.10
optimum <- select_region_union_epsilon(
  log_evidence_lower = reference$log_mass,
  log_evidence_upper = reference$log_mass,
  costs = reference$costs,
  epsilon = epsilon,
  mode = "conservative"
)
tree <- new_skeleton_region_tree(
  states = dag3$states,
  log_scores = scored$log_score,
  log_prior = scored$log_prior,
  exact_evidence_bounds = TRUE
)
score$reset()
fit <- select_minimum_cost_region_tree(
  tree = tree,
  log_score_fn = score$score,
  evidence_evaluator = function(node, remaining_score_budget = Inf) {
    stop("An exact terminal was unexpectedly reevaluated.", call. = FALSE)
  },
  epsilon = epsilon,
  cost_tolerance = 0,
  max_score_evaluations = 10000,
  max_actions = 10000
)
selected_cost <- sum(vapply(
  fit$nodes[fit$selected_ids],
  function(node) node$region$cost,
  numeric(1)
))
expect_true(fit$bounds_met, "The exact skeleton search did not close.")
expect_true(
  selected_cost == optimum$budget,
  "The exact skeleton search did not recover the reference optimum."
)
expect_true(
  identical(unname(score$diagnostics()[["calls"]]),
            unname(fit$total_score_evaluations)),
  "Score-request accounting disagrees with the search."
)

cat("Bayesian-network helper and integration tests passed.\n")
