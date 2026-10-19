# Data sources

## Tecator

`tecator/tecator_source.txt` is the original StatLib file, including the Tecator permission and attribution note. The measurements were recorded with the Tecator Infratec Food and Feed Analyzer. The original notice permits redistribution with that note attached.

Source: https://lib.stat.cmu.edu/datasets/tecator

`tecator.csv` preserves the supplied 240 observations and split labels. The paper uses the first 172 observations for training and the next 43 for testing, excluding the 25 extrapolation observations. It selects 16 channels, standardizes using training statistics and groups adjacent selected channels into four groups. These steps are in `experiments/tecator/run.R`; the original 100-channel input is retained so preprocessing is reproducible.

## Sachs

`sachs/raw/` contains the 853-observation data and supplied three-level discretization from the bnlearn book companion. `sachs/processed/` retains the exact processed arrays used in the paper. The source URLs, input hashes, variable names and processing notes are recorded in `sachs/source_manifest.json`.

The paper cites Sachs et al. (2005). No new data license is asserted by this package; the software MIT license does not replace the source data's attribution or terms. `experiments/sachs/prepare_data.py` contains the original processing steps restricted to these data and uses the supplied raw copies, without downloading or refitting discretization boundaries.

Synthetic datasets are specified by the complete retained seed registries and generation code under `experiments/`. Their large intermediate matrices are regenerated only when a full experiment is explicitly executed.
