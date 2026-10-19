# Regional Summaries of Bayesian Model Uncertainty

Code, data and saved results for the current TMLR manuscript, including its appendices. Start with the saved-result workflow below. It does not fit models, run searches, generate datasets or recompute bootstrap intervals.

## Rebuild the reported results

Use Python 3.12 or later, R, and Ghostscript. Install the Python plotting dependencies with `python -m pip install -r requirements.txt`. R's standard graphics package draws Figure 3; no additional R package is needed for this rebuild. Make `Rscript` and `gs` available on your executable path, or pass their locations with `--rscript` and `--gs`.

```sh
python scripts/verify_package.py
python scripts/rebuild.py
```

The rebuild writes four empirical figures and Tables 4, 5 and 7–12 to `build/`, alongside their numerical values. It checks the paired point estimates against the original analysis and reads the saved bootstrap endpoints. The supplied `figures/` PDFs are the exact manuscript assets. Rendering can differ across graphics-library versions; the underlying numerical values are checked separately.

`requirements-rebuild-tested.txt` records the Python versions used for the successful release check. The check environment is recorded separately from the original experiment environment in `docs/`.

## Contents

| Directory | Contents |
|---|---|
| `src/` | BMA implementations and the two Bayesian-network implementations |
| `experiments/` | Data-generation code, registered settings, seeds and study workers |
| `data/` | Tecator and Sachs inputs, with source information |
| `results/` | Saved outcomes, reported comparisons and supporting audit records |
| `scripts/` | Saved-result rebuild, integrity checks and explicit experiment launchers |
| `figures/` | The four external PDFs used in the manuscript |
| `docs/` | Paper-to-code map, reproduction details, provenance and environment |

The retained results comprise 1,200 main/stress BMA estimator outcomes, 10,400 component-study outcomes, 144 matched-selector outcomes, 1,600 network-simulation outcomes, all six Tecator conditions, and all three Sachs posterior specifications. The component study includes 9,600 primary searches on 2,400 independent datasets. Repeated targets and variants do not constitute new independent datasets. The two MC3 estimators share each sampling run.

Every attempt belonging to a retained comparison remains included, including failures and nonreturns. The unreported optimizer ablation, sampler-repeat family, larger-space extension, obsolete real-data analysis and historical development runs are excluded. Selection rules and original input hashes are recorded in `docs/source_provenance.json`.

## Rerun an experiment explicitly

Full experiments are optional and can be expensive. They need the additional dependencies in `requirements-experiments.txt`, the R packages listed in `requirements-r.R`, and a C++17 toolchain with MPFR and GMP development libraries. See [reproduction details](docs/reproduction.md).

```sh
# Prepare a new workspace without running any scientific computation.
python scripts/run_experiment.py bma_scaling --output ../scaling-run

# Prepare another new workspace and explicitly start the full experiment.
python scripts/run_experiment.py bma_scaling --output ../scaling-full --execute
```

Available study names are `bma_scaling`, `bma_components`, `bma_selectors`, `bayesian_network`, `tecator` and `sachs`. The launcher requires a new directory outside this package, so it cannot overwrite the supplied results. Time-limited searches may produce different stopping outcomes on a different machine. The packaged results remain the evidence used in the paper.

See [the paper map](docs/paper_map.md) for exact figure/table inputs and [the code guide](docs/code_guide.md) for the implementation entry points. Saved records and the clean-workspace preparation were checked for this release; the full experiments were not rerun.

## Citation and licenses

`CITATION.cff` identifies the manuscript authors and this local code-package version. No publication DOI or GitHub release is claimed. The existing repository's MIT notice is preserved verbatim in `LICENSE`; inherited notices are retained under `docs/licenses/`. The original Tecator redistribution note is included with its data. Data-source attribution and terms remain separate from the software license; see [data sources](data/README.md).

This is a named package prepared for the public repository. The anonymous TMLR submission package is maintained separately.
