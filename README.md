# Mastitis detection in Murrah buffaloes from infrared thermography

R pipeline accompanying the manuscript:

> **[Non-invasive mastitis screening in Murrah buffaloes using quarter-level infrared thermography and machine learning]**
> [Ekta Hooda1#, Sunesh Balhara1#, Mehar Singh Khatkar2, Ashok Boora1, Sarita Yadav1, Manish Tiwari1, Sanjay Choudhary1, Savita Nandal1, SK Phulia1, Mustafa Hasan Jan1, FC Tuteja1 and Ashok Kumar Balhara1*]


This repository contains the full analysis code used to build and evaluate machine-learning classifiers for clinical mastitis (CM), subclinical mastitis (SCM), and any-mastitis detection in Murrah buffaloes using udder-quarter infrared thermography (IRT) features.

---

## What the pipeline does

Starting from an Excel sheet of per-animal thermal readings (four udder quarters + four teat quarters, mean and max temperatures) plus animal metadata (parity, DIM, age) and environmental covariates (ambient temperature, relative humidity), the pipeline:

1. **Engineers features** — inter-quarter differentials, standard deviations, ranges, and asymmetry indices from the raw thermal columns.
2. **Splits animal-aware** — no animal appears in both train and test sets.
3. **Trains three feature-set variants** — `ANIMAL` (metadata only), `THERMALENV` (thermal + environment), `SIMPLE` (headline set: thermal + environment + metadata + engineered asymmetry).
4. **Fits Random Forest and XGBoost separately** — no ensembling; each model reported on its own.
5. **Calibrates probabilities** — isotonic recalibration on cross-validation out-of-fold predictions for both RF and XGB.
6. **Bakes decision thresholds into the model** — Youden's-J thresholds derived on calibrated OOF scores, packaged with each classifier as a self-contained object.
7. **Evaluates four tasks** — `healthy vs CM`, `healthy vs SCM`, `healthy vs any-mastitis`, and a three-class `healthy / CM / SCM` formulation.
8. **Runs feature-importance analysis** — Boruta, permutation importance, and SHAP (via `fastshap`).


---

## Repository layout

```
mastitis-irt-ml/
├── R/
│   └── pipeline_final_12june.R   # the analysis pipeline
├── data/
│   └── README.md                 # data dictionary (raw data not shared)
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
- **xgboost 1.7.x** — the pipeline uses `caret`'s `xgbTree` method, which is **not compatible with xgboost ≥ 3.0**. Install the 1.7.x binary before running:

  ```r
  # Windows
  install.packages("xgboost_1.7.9.1.zip", repos = NULL, type = "binary")
  # macOS / Linux
  remotes::install_version("xgboost", version = "1.7.9.1")
  ```

The remaining packages install automatically on first run:

```
readxl, dplyr, lubridate, caret, randomForest, xgboost, pROC, ggplot2,
tibble, doParallel, tidyr, foreach, matrixStats, gridExtra, MLmetrics,
reshape2, Boruta, fastshap
```

### Running the pipeline

1. Clone this repo:
   ```bash
   git clone https://github.com/ektahooda/mastitis-irt-ml.git
   cd mastitis-irt-ml
   ```
2. Open the project in RStudio.
3. Edit the two paths at the top of `R/pipeline_final_12june.R`:
   ```r
   data_path      <- "path/to/your/data.xlsx"
   checkpoint_dir <- "path/to/output/folder"
   ```
4. Source the script. Runtime is on the order of tens of minutes on a laptop, depending on cores available for the parallel backend.

Outputs land in `checkpoint_dir/`:
- `manuscript_main/` — main-text figures and CSV tables
- `manuscript_supp/` — supplementary figures
- `*.rds` — trained model bundles (each includes the caret model, isotonic calibrator, threshold, and feature list)

---

## Data availability

The raw dataset is **not included** in this repository. It is available from the corresponding author on reasonable request, subject to institutional data-sharing policies. See `data/README.md` for the expected column schema so that the pipeline can be run on comparably structured data.

---

## Reproducibility notes

- The pipeline uses `set.seed()` inside each modelling step. The exact numerical outputs depend on the versions listed above; different xgboost or randomForest versions may shift metrics in the third decimal.
- Splits are animal-aware and stratified by outcome. Bootstrap CIs use 2000 iterations.
- Calibration is done on out-of-fold predictions inside the training set only — the held-out test set is never touched during model or threshold selection.

---

## Citation

If you use this code, please cite the manuscript above and the software release (see `CITATION.cff`).

## License

MIT — see [`LICENSE`](LICENSE).

## Contact

[Ashok Kumar Balhara] — [balharaak@gmail.com] [ICAR-CIRB]
[Institution]
