# Reproduction details

## Saved-result workflow

`verify_package.py` checks file hashes, registry coverage, required counts, retained comparison families, figure inputs and recorded scientific outcomes. `rebuild.py` recomputes descriptive table entries and renders plots from saved data. It does not source the search implementations or call experiment launchers. It uses the existing confidence intervals and checks their point estimates against the saved paired observations.

The four files in `figures/` are byte-identical to the current manuscript assets. Rebuild outputs go elsewhere. Figure 3 uses R graphics followed by Ghostscript font embedding. Python figures use embedded TrueType fonts. Cross-platform font rendering and PDF metadata need not be byte-identical.

## Full experiment dependencies and preparation

- Python: install `requirements-experiments.txt` for the Sachs numerical study, in addition to `requirements.txt` for presentation and analysis.
- R: `Rscript requirements-r.R` installs Rcpp, jsonlite and yaml. Rmpfr is not required; the C++ selectors use MPFR directly.
- Native compilation: a C++17 compiler, MPFR headers/libraries and GMP headers/libraries must be visible to R's toolchain. Configure `PKG_CPPFLAGS` and library search paths when these are outside standard locations. On Windows, the loader accepts `ARS_RTOOLS_HOME` or `RTOOLS44_HOME` and otherwise checks the original Rtools locations.

`run_experiment.py` copies common code and the appropriate version overlay into a fresh working directory, installs the unchanged study settings and writes a new integrity manifest for that prepared workspace. This manifest records the clean release; it does not claim that the reorganized files are the original preregistration. Original identifiers, seeds and source hashes remain available in the release provenance.

The default launcher only prepares files. `--execute` starts scientific computation. The launcher does not copy saved results into a new experiment, so reruns start from fresh workspaces. It runs jobs sequentially. The original component study used one to three concurrent workers, and its 60-second limits can therefore produce different outcomes under a serial rerun. Seed preservation does not imply identical timing-dependent results.

| Study | Registered scope retained here | Numerical implementation |
|---|---|---|
| `bma_scaling` | 200 datasets, two budgets, 800 practical searches, 1,200 estimator outcomes | Newton envelopes, mass propagation enabled; independent four-chain MC3 |
| `bma_components` | 2,400 datasets, 9,600 primary plus 800 reported component/order/budget outcomes | Staged implementation, BFGS envelopes, mass propagation disabled |
| `bma_selectors` | 48 datasets, three selectors, 144 outcomes | Staged implementation, shared traversal/proposals/refinement |
| `bayesian_network` | 200 datasets, four tolerances, two tree orders, 1,600 outcomes | Five-node exhaustive numerical masses and common tree selection |
| `tecator` | Three noise levels, two targets, six searches | Newton envelopes, mass propagation disabled |
| `sachs` | Primary plus two prior sensitivities; 90/95/99% targets | Separate longest-path-layer numerical sums and largest-mass search |

Practical BMA searches finish before exhaustive audit references are constructed. The network simulation deliberately uses complete numerical masses before tree selection. Sachs uses its separate numerical algorithm; it is not a stochastic validation of Algorithm 1.

## Optional statistical and interpretation computations

The saved-result workflow is sufficient to recover the presented evidence. For an explicitly requested reanalysis, the optional interval scripts accept a new output directory and recompute the original bootstrap calculations. They are never invoked by `rebuild.py`.

`experiments/tecator/interpret.R` contains the original group, channel and prediction calculations adapted to the current saved record. It additionally enumerates supports for predictions, so it is not part of the ordinary saved-result checks. To run it explicitly, prepare a Tecator workspace, copy `results/tecator/records/` to that workspace's `validation/results/numerical_end_to_end/`, then run `Rscript interpret.R` from the workspace. No search is performed by this interpretation script.

`Rscript scripts/check_saved_tecator.R` checks the supplied six search records, audits and primary retained union without performing enumeration. `extract_trace_diagnostics.R` accepts three arguments: a completed component workspace, a completed selector workspace and an output directory. It extracts diagnostics from those full RDS archives; the resulting 644-row export is already supplied here.

## Recorded environment and limits of this release check

The main BMA results were recorded on Windows 11 build 26300, R 4.4.2 and Python 3.12.14, with one practical worker. Processor and RAM information was not recorded. Search timing includes preparation and selection but excludes data loading, interpreter startup, compilation and auditing. Equal score ceilings are not equal elapsed-time budgets. Sachs reused saved regional caches and does not provide a cold-start timing comparison. See `recorded_environment.json`.

This release was checked by rebuilding saved tables/figures, reading the Tecator saved records, parsing sources and preparing every study in a clean directory. Full simulations and scientific search reruns were not performed. The package includes complete retained outcome tables and seeds, plus the compact records needed for the reported diagnostics and real-data interpretations. It does not include every large intermediate RDS file from the simulation archive.
