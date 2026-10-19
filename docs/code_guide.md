# Code guide

The release keeps the scientific implementations used for the reported experiments. It reorganizes files and supplies portable launchers; it does not merge historically different algorithms into a single revised implementation.

## BMA

`src/bma/common/` holds files that are identical across the staged and Newton versions. `src/bma/staged/` and `src/bma/newton/` supply the files that differ. A prepared workspace overlays exactly one version on the common files. The network simulation also reuses the staged `vendor/region_restricted_smc.R`, stored only once.

The practical entry point is `production_search(d, condition)` in `common/production.R`. The data object contains training predictors and response, predictor groups, tree depth and noise scale; real-data objects also retain standardization information. The condition supplies tolerance, score/time ceilings, algorithm seed and group order. Study-specific workers construct these objects from the supplied settings and seeds.

- `common/frozen_runner_source.R` contains the original `run_one` function, extracted without changing its expression. The obsolete top-level historical launcher has been removed.
- `base_study.R` and the version's `vendor/` files define the posterior scorer, data construction and regional objects.
- `src/gaussian_mass_bound.R` and `sampling.R` implement bounds, proposals and regional evidence calculations.
- `integration.R`, `two_stage.R`, `adaptive.R` and `refinement.R` organize exploration and refinement.
- `common/src/selection.cpp` implements the interval selector with directed rounding. The Newton version additionally contains `src/mass_conservation.cpp` and its R wrapper.
- `common/audit/` and the scaling study's streaming enumerator supply independent post-search numerical references.

Some internal helper files retain inherited names such as `pilot.R` and `load_repair.R`. They are dependencies of the production code, not additional pilot studies. Their names are preserved to keep source relationships and provenance easy to check. Development launchers, logs, archived releases and unrelated tests are not included.

## Bayesian networks

The five-node simulation combines `src/bayesian_network/exact_bayesian_network_helpers.R` with the staged common tree-search library. Its protocol explicitly uses complete numerical enumeration rather than stochastic mass refinement.

`src/sachs/network.py` is the separate real-data implementation: local BDeu scores, longest-path-layer DAG sums, largest-mass regional expansion and HiGHS selection. The small `prepare.py` alongside it supplies only local path/JSON helpers. The Sachs preprocessing script retains the original supplied-discretization transformation and reads the packaged raw inputs.

## Presentation and analysis

`build_tables.py` uses original saved outcomes to regenerate descriptive values. Plotting scripts are adapted from the final manuscript plotting code by changing paths and retained-row counts. Original saved bootstrap intervals are supplied; optional resampling scripts preserve the original analysis seeds and order of the retained comparisons.

The original submission supplement and historical repository remain separate. `source_provenance.json` records copied source hashes and deterministic filters. `MANIFEST.json` records every file in this clean release, excluding the manifest itself and generated `build/` outputs.
