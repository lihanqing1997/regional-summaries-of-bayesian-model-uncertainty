# Paper-to-code map

Target: **Regional Summaries of Bayesian Model Uncertainty**, current TMLR manuscript. Paths are relative to this package. The manuscript and its anonymous supplement are separate submission artifacts.

| Paper item | Saved evidence | Code / settings |
|---|---|---|
| Figure 2; Table 4; Table 7(a,b); main and stress comparisons at both budgets | `results/bma_scaling/outcomes.csv` (1,200 rows), saved pooled summaries and paired intervals | `experiments/bma_scaling/`; `src/bma/common/` + `src/bma/newton/`; `scripts/plot_scaling.py` |
| Figure 3; Tables 8–9; interval, allocation, budget and tree-order comparisons | `results/bma_components/outcomes.csv` (10,400 rows), `figure3_data.csv` (200 paired outcomes), `paired_comparisons.csv` | `experiments/bma_components/config/`; `src/bma/common/` + `src/bma/staged/`; `scripts/plot_components.R` |
| Table 10 and search-diagnostic discussion | The 500 primary rows in `results/bma_components/saved_trace_diagnostics.csv`, joined to the 90% primary outcomes | `scripts/build_tables.py`; original extraction logic in `scripts/extract_trace_diagnostics.R` |
| Exact-mass ranked-region comparison in the component-study appendix | `results/bma_components/exact_reference_comparators.csv`: all 2,400 datasets at the 90% target | `src/bma/common/production.R` (`production_comparators`), independent audit enumerator |
| Matched-selector comparison in the appendix | `results/bma_selectors/outcomes.csv`: 144 outcomes, three methods on 48 datasets; matched trace rows in `saved_trace_diagnostics.csv` | `experiments/bma_selectors/`, staged BMA implementation |
| Figure 4; network-simulation results and tree-order comparison | `results/bayesian_network/bayesian_network_exact.csv`: 1,600 outcomes on 200 datasets | `experiments/bayesian_network/`, `src/bayesian_network/`, `scripts/plot_network.py` |
| Table 5; Tecator interpretation and prediction statements | `results/tecator/real_patterns.csv`, group/channel queries, `real_interpretation.json`, current saved primary record | `experiments/tecator/interpret.R`; `data/tecator/` |
| Table 11 and Tecator sensitivities | `results/tecator/outcomes.csv` (six rows); all six search records/audits and three numerical references in `records/` | `experiments/tecator/run.R`, Newton BMA implementation with mass propagation disabled |
| Figure 5; Table 12; Sachs interpretation and ranked-DAG comparison | `results/sachs/`, with all three specifications, saved regional sums, graph selections and comparisons | `src/sachs/network.py`; `experiments/sachs/protocol.json`; `scripts/plot_sachs.py` |

`scripts/build_tables.py` rebuilds the empirical table values and `scripts/rebuild.py` produces all four empirical figures. Tables 1–3 and 6 are conceptual comparisons, stipulated worked examples and notation. Figure 1 is a conceptual tree drawn directly in the manuscript's LaTeX source. They have no simulation inputs and are not duplicated in this code package.

The published intervals in Table 9 and the main BMA text come from the saved original bootstrap analyses. The ordinary rebuild does not draw new bootstrap samples. The optional `recompute_component_intervals.py` and `recompute_scaling_intervals.py` preserve the original resampling formulas and seeds.

Interpretation records for Tecator were computed for the same retained union as the current primary search. They supply its audited mass, group/channel probabilities and test predictions. The current search's omitted-mass bound and complexity bracket come from `outcomes.csv` and the current RDS record, not an earlier interpretation checkpoint.
