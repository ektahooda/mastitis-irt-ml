# Mastitis detection in Murrah buffaloes from infrared thermography

R pipeline accompanying the manuscript:

> **[Non-invasive mastitis screening in Murrah buffaloes using quarter-level infrared thermography and machine learning]**
> [Ekta Hooda1#, Sunesh Balhara1#, Mehar Singh Khatkar2, Ashok Boora1, Sarita Yadav1, Manish Tiwari1, Sanjay Choudhary1, Savita Nandal1, SK Phulia1, Mustafa Hasan Jan1, FC Tuteja1 and Ashok Kumar Balhara1*]

This repository contains the modelling code used to train and evaluate machine-learning classifiers for clinical mastitis (CM), subclinical mastitis (SCM), and any-mastitis detection in Murrah buffaloes using udder-quarter infrared thermography (IRT) features.

Figure-production code for the manuscript is **not** included here — this repository is the models, feature engineering, and quantitative evaluation only. The reported metrics can be reproduced from the CSV files the pipeline writes.

---

## What the pipeline does

Starting from an Excel sheet of per-animal thermal readings (four udder quarters + four teat quarters, mean and max temperatures) plus animal metadata (parity, DIM, age) and environmental covariates (ambient temperature, relative humidity), `R/pipeline.R`:

1. **Engineers features** — inter-quarter differentials, standard deviations, ranges, and asymmetry indices from the raw thermal columns.
2. **Splits animal-aware** — no animal appears in both train and test sets (following Bobbo et al., 2023, *J Dairy Sci*).
3. **Trains three feature-set variants** — `ANIMAL` (metadata only), `THERMALENV` (thermal + environment), `SIMPLE` (headline set: thermal + environment + metadata + engineered asymmetry).
4. **Fits Random Forest and XGBoost separately** — no ensembling; each model is reported on its own.
5. **Selects features per (task × feature-set)** with Boruta.
6. **Calibrates probabilities** — isotonic recalibration on cross-validation out-of-fold predictions for both RF and XGB.
7. **Bakes decision thresholds into the model** — Youden's-J thresholds derived on calibrated OOF scores, packaged with each classifier as a self-contained object with `predict_prob()` and `predict_class()` methods.
8. **Evaluates four tasks** — `healthy vs CM`, `healthy vs SCM`, `healthy vs any mastitis`, and a three-class `healthy / CM / SCM` formulation on the held-out test set.
9. **Writes summary CSVs** with 2000-iteration bootstrap CIs.

The pipeline uses checkpointing: intermediate Boruta selections and trained models are cached to `results/models/`, so re-runs skip completed steps.

---

## Repository layout

```
mastitis-irt-ml/
├── R/
│   └── pipeline.R                # the modelling pipeline
├── data/
│   └── README.md                 # data dictionary (raw data not shared)
├── results/                      # created on first run; git-ignored
│   ├── models/                   # trained classifier bundles (.rds)
│   └── metrics/                  # binary_metrics.csv, threeclass_metrics.csv
├── LICENSE                       # MIT
├── CITATION.cff                  # how to cite this code
├── .gitignore
└── README.md                     # this file
```

---

## Getting started

### Requirements

- **R** ≥ 4.2
- **RStudio** (recommended)
- **xgboost 1.7.x** — the pipeline uses `caret`'s `xgbTree` method, which is **not compatible with xgboost ≥ 3.0**. Install the 1.7.x version once:

  ```r
  # any OS
  remotes::install_version("xgboost", version = "1.7.9.1")
  ```

  The pipeline checks this at startup and stops with a clear message if a 3.x version is loaded.

The remaining packages install automatically on first run:

```
readxl, dplyr, lubridate, caret, randomForest, xgboost, pROC, ggplot2,
tibble, doParallel, tidyr, foreach, matrixStats, MLmetrics, reshape2,
Boruta, here
```

### Running the pipeline

1. Clone the repo:
   ```bash
   git clone https://github.com/ektahooda/mastitis-irt-ml.git
   cd mastitis-irt-ml
   ```
2. Place your dataset at `data/mastitis_data.xlsx` (see [`data/README.md`](data/README.md) for the expected schema).
3. Open the project in RStudio and source `R/pipeline.R`.

That's it — the pipeline resolves all paths relative to the repo root using the `here` package, so no manual path editing is needed. Runtime is on the order of tens of minutes on a laptop, depending on cores available.

### Outputs

After a successful run, `results/` contains:

- **`models/*.rds`** — one bundle per (task × feature-set). Each bundle holds `rf_classifier` and `xgb_classifier` objects. Use them on new data like this:

  ```r
  bundle <- readRDS("results/models/healthy_vs_CM_SIMPLE_bundle.rds")
  probs  <- bundle$xgb_classifier$predict_prob(new_data)
  preds  <- bundle$xgb_classifier$predict_class(new_data)
  ```

  The threshold used inside `predict_class()` is the model's own attribute (derived from calibrated CV OOF predictions in training) — no external threshold argument to pass.

- **`metrics/binary_metrics.csv`** — F1, MCC, AUC, sensitivity, specificity, confusion counts for every (task × feature-set × algorithm) combination on the held-out test set. This is the source for Table 2 in the manuscript.

- **`metrics/threeclass_metrics.csv`** — per-class AUC / sensitivity / specificity / F1 plus overall accuracy and macro-AUC for the three-class model.

---

## Data availability

The raw dataset is **not included** in this repository. It is available from the corresponding author on reasonable request, subject to institutional data-sharing policies. See [`data/README.md`](data/README.md) for the expected column schema so that the pipeline can be run on comparably structured data.

---

## Reproducibility notes

- `set.seed(2024)` is set at the top of the pipeline; individual modelling steps also seed inside. The exact numerical outputs depend on the package versions listed above — different `xgboost` or `randomForest` versions may shift metrics in the third decimal.
- Splits are animal-aware and stratified by outcome. Bootstrap CIs use 2000 iterations.
- Calibration is done on out-of-fold predictions inside the training set only — the held-out test set is never touched during model or threshold selection.

---

## Citation

If you use this code, please cite the manuscript above and the software release (see [`CITATION.cff`](CITATION.cff)).

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

[Ashok Kumar Balhara] — [balharaak@gmail.com] [ICAR-CIRB]
