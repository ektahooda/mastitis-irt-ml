# ================================================================
# Mastitis-IRT classification pipeline
# ================================================================

# This script trains and evaluates Random Forest and XGBoost
# classifiers for four mastitis-detection tasks in Murrah buffaloes
# from udder-quarter infrared thermography features:
#
#   1. healthy vs CM       (binary)
#   2. healthy vs SCM      (binary)
#   3. healthy vs any      (binary)
#   4. healthy / CM / SCM  (three-class)
#
# Key design points:
#
#   - Animal-aware train/test split (no leakage across animals)
#   - Isotonic probability calibration on CV out-of-fold predictions
#   - Decision thresholds (Youden's J) baked into each classifier
#     object as a model attribute
#   - Boruta feature selection per (task x feature-set)
#   - RF and XGB reported separately (no ensembling)
#
# Outputs written to results_dir/ :
#
#   models/    trained-classifier bundles (.rds), one per task
#              each bundle contains rf_classifier and xgb_classifier
#              objects with predict_prob() and predict_class() methods
#   metrics/   binary_metrics.csv and threeclass_metrics.csv
#
# Figure-production code for the manuscript is NOT included here;
# see the manuscript itself and the metric CSVs for the reported
# numbers.
# ================================================================


# ---------------------------------------------------------------
# SECTION 0 — SETUP
# ---------------------------------------------------------------

# NOTE: this pipeline uses caret's xgbTree method, which is not
# compatible with xgboost >= 3.0. Install xgboost 1.7.x once before
# running:
#
#   remotes::install_version("xgboost", version = "1.7.9.1")
#
# (or the equivalent Windows/Mac binary).

packages <- c("readxl", "dplyr", "lubridate", "caret", "randomForest",
              "xgboost", "pROC", "ggplot2", "tibble",
              "doParallel", "tidyr", "foreach", "matrixStats",
              "MLmetrics", "reshape2", "Boruta", "here")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

xgb_ver <- as.character(packageVersion("xgboost"))
cat("xgboost version:", xgb_ver, "\n")
if (as.integer(strsplit(xgb_ver, "\\.")[[1]][1]) >= 3) {
  stop("xgboost >= 3.0 detected. caret::xgbTree requires 1.7.x.\n",
       "  Install with: remotes::install_version('xgboost', '1.7.9.1')")
}

# ---------------- PATHS ----------------
# Uses `here` to resolve paths relative to the repository root, so
# the pipeline works on any machine without editing. Put your input
# .xlsx in data/ (see data/README.md for the expected schema).

data_path   <- here::here("data", "mastitis_data.xlsx")
results_dir <- here::here("results")
models_dir  <- file.path(results_dir, "models")
metrics_dir <- file.path(results_dir, "metrics")

dir.create(models_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(metrics_dir, showWarnings = FALSE, recursive = TRUE)

# For backward compatibility with the original section names below,
# keep `checkpoint_dir` and `fig_main` as aliases.
checkpoint_dir <- models_dir
fig_main       <- metrics_dir

if (!file.exists(data_path)) {
  stop("Input file not found: ", data_path,
       "\n  Place your dataset at data/mastitis_data.xlsx",
       "\n  or edit `data_path` above.")
}

set.seed(2024)


# SECTION 1 — DATA LOADING & FEATURE ENGINEERING
# Only SIMPLE features are computed (no _adj, no z-scores,
# no DIM splines, no interaction terms).
# ================================================================

engineer_features <- function(df) {
  df %>%
    mutate(
      DOB = as.Date(DOB), DOR = as.Date(DOR), DOC = as.Date(DOC),
      Age_years  = as.numeric(DOR - DOB) / 365.25,
      DIM        = as.numeric(DOR - DOC),
      
      # 3-category clinical lactation stage convention
      Lact_stage = case_when(DIM <= 90  ~ 1L,   # early
                             DIM <= 180 ~ 2L,   # mid
                             TRUE       ~ 3L),  # late
      
      # Teat-udder differentials (inflammation proxy)
      Diff_T1 = U1 - T1, Diff_T2 = U2 - T2,
      Diff_T3 = U3 - T3, Diff_T4 = U4 - T4,
      Diff_T1max = U1max - T1max, Diff_T2max = U2max - T2max,
      Diff_T3max = U3max - T3max, Diff_T4max = U4max - T4max,
      
      # Quarter averages
      Avg_T    = (T1 + T2 + T3 + T4) / 4,
      Avg_U    = (U1 + U2 + U3 + U4) / 4,
      Avg_Tmax = (T1max + T2max + T3max + T4max) / 4,
      Avg_Umax = (U1max + U2max + U3max + U4max) / 4,
      
      # Inter-quarter range
      T_range    = pmax(T1,T2,T3,T4)    - pmin(T1,T2,T3,T4),
      U_range    = pmax(U1,U2,U3,U4)    - pmin(U1,U2,U3,U4),
      Tmax_range = pmax(T1max,T2max,T3max,T4max) -
        pmin(T1max,T2max,T3max,T4max),
      Umax_range = pmax(U1max,U2max,U3max,U4max) -
        pmin(U1max,U2max,U3max,U4max),
      
      # Inter-quarter SD
      T_sd    = rowSds(as.matrix(cbind(T1,T2,T3,T4)),             na.rm=TRUE),
      U_sd    = rowSds(as.matrix(cbind(U1,U2,U3,U4)),             na.rm=TRUE),
      Tmax_sd = rowSds(as.matrix(cbind(T1max,T2max,T3max,T4max)), na.rm=TRUE),
      Umax_sd = rowSds(as.matrix(cbind(U1max,U2max,U3max,U4max)), na.rm=TRUE),
      
      # Max single-quarter differential
      Max_diff     = pmax(Diff_T1,    Diff_T2,    Diff_T3,    Diff_T4),
      Max_diff_max = pmax(Diff_T1max, Diff_T2max, Diff_T3max, Diff_T4max),
      
      # Labels
      Label_CM  = factor(case_when(
        Status == "CM"      ~ "disease",
        Status == "healthy" ~ "healthy",
        TRUE                ~ NA_character_),
        levels = c("healthy","disease")),
      Label_SCM = factor(case_when(
        Status == "SCM"     ~ "disease",
        Status == "healthy" ~ "healthy",
        TRUE                ~ NA_character_),
        levels = c("healthy","disease")),
      Label_Any = factor(ifelse(Status %in% c("CM","SCM"),
                                "disease","healthy"),
                         levels = c("healthy","disease")),
      Label_3class = factor(case_when(
        Status == "CM"      ~ "CM",
        Status == "SCM"     ~ "SCM",
        Status == "healthy" ~ "healthy",
        TRUE                ~ NA_character_),
        levels = c("healthy","CM","SCM"))
    )
}

cat("\n", strrep("=",60), "\n  LOADING DATA\n", strrep("=",60), "\n")

df_raw <- read_excel(data_path)
cat("Raw dimensions:", nrow(df_raw), "x", ncol(df_raw), "\n")
cat("Status values:\n"); print(table(df_raw$Status))

junk_cols <- grep("^\\.\\.\\.\\d+$", colnames(df_raw), value = TRUE)
if (length(junk_cols) > 0) {
  df_raw <- df_raw %>% select(-all_of(junk_cols))
}

numeric_cols <- c("T1","T2","T3","T4","U1","U2","U3","U4",
                  "T1max","T2max","T3max","T4max",
                  "U1max","U2max","U3max","U4max",
                  "AT","RH","Parity")
cols_present <- intersect(numeric_cols, colnames(df_raw))
df_raw[cols_present] <- lapply(df_raw[cols_present], as.numeric)

required_cols <- c("animal_no","Status","DOB","DOR","DOC",
                   "T1","T2","T3","T4","U1","U2","U3","U4",
                   "T1max","T2max","T3max","T4max",
                   "U1max","U2max","U3max","U4max",
                   "AT","RH","Parity")
required_present <- intersect(required_cols, colnames(df_raw))

df <- df_raw %>%
  engineer_features() %>%
  filter(if_all(all_of(required_present), ~ !is.na(.)))

cat("Rows after NA filter:", nrow(df), "\n")
cat("Label distributions:\n")
cat("  CM     :"); print(table(df$Label_CM))
cat("  SCM    :"); print(table(df$Label_SCM))
cat("  Any    :"); print(table(df$Label_Any))
cat("  3-class:"); print(table(df$Label_3class))
cat("Lactation stage counts (1=early, 2=mid, 3=late):\n")
print(table(df$Lact_stage))


# ================================================================
# SECTION 2 — FEATURE SET DEFINITIONS
# ================================================================

features_animal_only <- c("Parity", "Age_years", "DIM", "Lact_stage")

features_thermal_env <- c(
  "T1","T2","T3","T4","U1","U2","U3","U4",
  "T1max","T2max","T3max","T4max",
  "U1max","U2max","U3max","U4max",
  "AT","RH")

features_simple <- c(
  "T1","T2","T3","T4","U1","U2","U3","U4",
  "T1max","T2max","T3max","T4max",
  "U1max","U2max","U3max","U4max",
  "Avg_T","Avg_U","Avg_Tmax","Avg_Umax",
  "Diff_T1","Diff_T2","Diff_T3","Diff_T4",
  "Diff_T1max","Diff_T2max","Diff_T3max","Diff_T4max",
  "T_range","U_range","Tmax_range","Umax_range",
  "T_sd","U_sd","Tmax_sd","Umax_sd",
  "Max_diff","Max_diff_max",
  "AT","RH",
  "Parity","Age_years","DIM","Lact_stage")

cat("\nFeature set sizes:\n")
cat("  ANIMAL    :", length(features_animal_only), "features\n")
cat("  THERMALENV:", length(features_thermal_env), "features\n")
cat("  SIMPLE    :", length(features_simple),      "features\n")


# ================================================================
# SECTION 3 — HELPER FUNCTIONS
# ================================================================

save_ckpt <- function(obj, label) {
  path <- file.path(checkpoint_dir, paste0(label, ".rds"))
  saveRDS(obj, path); cat("  [SAVED]", path, "\n")
}
load_ckpt <- function(label) {
  path <- file.path(checkpoint_dir, paste0(label, ".rds"))
  if (file.exists(path)) { cat("  [LOADED]", path, "\n"); return(readRDS(path)) }
  NULL
}

calc_mcc <- function(TP, TN, FP, FN) {
  denom <- sqrt(as.numeric(TP+FP) * as.numeric(TP+FN) *
                  as.numeric(TN+FP) * as.numeric(TN+FN))
  if (denom == 0) return(NA_real_)
  (TP*TN - FP*FN) / denom
}

calc_auprc <- function(probs, labels) {
  y_bin <- as.integer(labels == "disease")
  ord   <- order(probs, decreasing = TRUE)
  y_ord <- y_bin[ord]
  n_pos <- sum(y_bin)
  if (n_pos == 0) return(NA_real_)
  tp <- cumsum(y_ord); fp <- cumsum(1 - y_ord)
  rec <- tp / n_pos
  pre <- tp / (tp + fp); pre[is.nan(pre)] <- 0
  n <- length(rec)
  sum((rec[2:n] - rec[1:(n-1)]) * (pre[2:n] + pre[1:(n-1)]) / 2)
}

calc_ece <- function(probs, labels, n_bins = 10) {
  y_bin  <- as.integer(labels == "disease")
  breaks <- seq(0, 1, length.out = n_bins + 1)
  bins   <- cut(probs, breaks, include.lowest = TRUE)
  n      <- length(probs); ece <- 0
  for (b in levels(bins)) {
    idx <- which(bins == b)
    if (length(idx) == 0) next
    ece <- ece + length(idx)/n *
      abs(mean(probs[idx]) - mean(y_bin[idx]))
  }
  ece
}

compute_youden_threshold <- function(obs, probs) {
  thresholds <- seq(0.02, 0.98, by = 0.01)
  youden_v   <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    preds <- factor(ifelse(probs >= thresholds[i], "disease", "healthy"),
                    levels = c("healthy", "disease"))
    cm <- suppressWarnings(confusionMatrix(preds, obs, positive = "disease"))
    s <- cm$byClass["Sensitivity"]; sp <- cm$byClass["Specificity"]
    youden_v[i] <- if (!is.na(s) && !is.na(sp)) s + sp - 1 else -Inf
  }
  thresholds[which.max(youden_v)]
}
# Specificity-prioritised threshold — used for CM where
# Youden's high-sensitivity operating point gives too many FPs.
# Picks the highest-spec operating point on OOF data subject
# to a minimum sensitivity floor.
compute_spec_priority_threshold <- function(obs, probs,
                                            min_sens = 0.85) {
  thresholds <- seq(0.02, 0.98, by = 0.01)
  sens_v <- spec_v <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    preds <- factor(ifelse(probs >= thresholds[i],
                           "disease", "healthy"),
                    levels = c("healthy", "disease"))
    cm <- suppressWarnings(
      confusionMatrix(preds, obs, positive = "disease"))
    sens_v[i] <- cm$byClass["Sensitivity"]
    spec_v[i] <- cm$byClass["Specificity"]
  }
  ok <- which(sens_v >= min_sens & !is.na(spec_v))
  if (length(ok) == 0) {
    return(thresholds[which.max(sens_v + spec_v - 1)])
  }
  thresholds[ok[which.max(spec_v[ok])]]
}

# ================================================================
# make_classifier() — packages a trained algorithm together with
# its calibrator, threshold, and feature list into one object
# with predict_prob() and predict_class() methods.
#
# The threshold is a model attribute (derived from CV OOF, training
# data only), not a post-hoc choice. This makes deployment trivial:
# pass new data to predict_class() and get a class prediction back.
# ================================================================

make_classifier <- function(model, calibrator, threshold, features,
                            algorithm_name) {
  cls <- list(
    algorithm  = algorithm_name,
    model      = model,
    calibrator = calibrator,
    threshold  = threshold,
    features   = features
  )
  cls$predict_prob <- function(X_new) {
    X_subset <- as.data.frame(X_new[, features, drop = FALSE])
    raw      <- predict(model, X_subset, type = "prob")[, "disease"]
    pmin(pmax(calibrator(raw), 0), 1)
  }
  cls$predict_class <- function(X_new) {
    probs <- cls$predict_prob(X_new)
    factor(ifelse(probs >= threshold, "disease", "healthy"),
           levels = c("healthy", "disease"))
  }
  class(cls) <- c("mastitis_classifier", "list")
  cls
}


# ================================================================
# SECTION 4 — BINARY CLASSIFIER FUNCTION
# Trains RF and XGB on a given (task, feature_set), calibrates
# both, finds Youden thresholds on calibrated OOF, returns
# self-contained classifier objects.
# ================================================================

run_classifier <- function(df, label_col, task_name,
                           feature_cols, variant_tag) {
  
  cat("\n", strrep("=", 60), "\n")
  cat(" TASK:", task_name, "| VARIANT:", variant_tag, "\n")
  cat(" Started:", format(Sys.time(), "%H:%M:%S"), "\n")
  cat(strrep("=", 60), "\n")
  
  task_tag <- paste0(gsub(" ", "_", task_name), "_", variant_tag)
  
  task_df <- df %>%
    mutate(outcome = .data[[label_col]]) %>%
    filter(!is.na(outcome)) %>%
    select(animal_no, outcome, all_of(feature_cols))
  
  cat("Class distribution:\n"); print(table(task_df$outcome))
  minority_ratio <- round(
    sum(task_df$outcome == "healthy") /
      sum(task_df$outcome == "disease"), 2)
  cat("Imbalance ratio:", minority_ratio, ":1\n")
  class_wts <- c(healthy = 1, disease = minority_ratio)
  
  # -------- Animal-aware 75/25 split --------
  set.seed(42)
  uniq      <- task_df %>% distinct(animal_no, outcome)
  tr_idx    <- createDataPartition(uniq$outcome, p = 0.75, list = FALSE)
  train_ids <- uniq$animal_no[tr_idx]
  test_ids  <- uniq$animal_no[-tr_idx]
  cat("Train animals:", length(train_ids),
      "| Test animals:", length(test_ids), "\n")
  
  # -------- Boruta --------
  boruta_ckpt <- file.path(checkpoint_dir,
                           paste0(task_tag, "_Boruta.rds"))
  if (file.exists(boruta_ckpt)) {
    cat("[Boruta] Loaded from checkpoint\n")
    boruta_result <- readRDS(boruta_ckpt)
  } else {
    cat("[Boruta] Running on training data...\n")
    train_boruta <- task_df %>%
      filter(animal_no %in% train_ids) %>%
      select(outcome, all_of(feature_cols)) %>%
      mutate(outcome = as.factor(outcome))
    set.seed(42)
    boruta_result <- Boruta(outcome ~ ., data = train_boruta,
                            maxRuns = 200, pValue = 0.01, doTrace = 0)
    boruta_result <- TentativeRoughFix(boruta_result)
    saveRDS(boruta_result, boruta_ckpt)
  }
  
  boruta_dec <- data.frame(
    Feature  = names(boruta_result$finalDecision),
    Decision = as.character(boruta_result$finalDecision))
  confirmed <- boruta_dec$Feature[boruta_dec$Decision == "Confirmed"]
  cat(sprintf("Boruta confirmed: %d of %d\n",
              length(confirmed), length(feature_cols)))
  if (length(confirmed) < 2) {
    cat("[WARN] <2 confirmed — using all candidates\n")
    confirmed <- feature_cols
  }
  all_feature_cols <- confirmed
  
  # -------- Animal-aware CV folds --------
  set.seed(42)
  train_animals_df <- task_df %>%
    filter(animal_no %in% train_ids) %>%
    mutate(.row_idx = row_number())
  
  unique_train_animals <- train_animals_df %>%
    distinct(animal_no, outcome) %>%
    group_by(outcome) %>%
    mutate(fold = sample(rep(1:5, length.out = n()))) %>%
    ungroup() %>%
    distinct(animal_no, .keep_all = TRUE) %>%
    select(animal_no, fold)
  
  train_with_fold <- train_animals_df %>%
    left_join(unique_train_animals, by = "animal_no",
              relationship = "many-to-one")
  
  fold_index <- lapply(1:5,
                       function(k) which(train_with_fold$fold != k))
  
  ctrl <- trainControl(
    method          = "cv", number = 5,
    index           = fold_index,
    classProbs      = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    verboseIter     = TRUE,
    allowParallel   = TRUE)
  
  train_df <- task_df %>%
    filter(animal_no %in% train_ids) %>%
    select(outcome, all_of(all_feature_cols))
  test_df  <- task_df %>%
    filter(animal_no %in% test_ids) %>%
    select(outcome, all_of(all_feature_cols))
  
  X_train <- as.data.frame(train_df %>% select(all_of(all_feature_cols)))
  y_train <- factor(train_df$outcome, levels = c("healthy","disease"))
  X_test  <- as.data.frame(test_df  %>% select(all_of(all_feature_cols)))
  y_test  <- factor(test_df$outcome,  levels = c("healthy","disease"))
  
  cat("Train rows:", nrow(X_train), "| Test rows:", nrow(X_test), "\n")
  
  n_cores <- max(1, detectCores() - 1)
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # -------- MODEL A: Random Forest --------
  rf_path <- file.path(checkpoint_dir, paste0(task_tag, "_RF.rds"))
  if (file.exists(rf_path)) {
    cat("[RF] Loaded from checkpoint\n")
    rf_model <- readRDS(rf_path)
  } else {
    cat("\n[RF] Training... ", format(Sys.time(),"%H:%M:%S"), "\n")
    rf_grid <- expand.grid(mtry = unique(pmax(1, c(
      floor(sqrt(ncol(X_train))),
      floor(ncol(X_train) / 4),
      floor(ncol(X_train) / 3),
      floor(ncol(X_train) / 2)))))
    set.seed(42)
    rf_model <- train(
      x = X_train, y = y_train, method = "rf", metric = "ROC",
      trControl = ctrl, tuneGrid = rf_grid,
      classwt = class_wts, ntree = 1000, importance = TRUE)
    cat("[RF] Done. Best mtry:", rf_model$bestTune$mtry, "\n")
    saveRDS(rf_model, rf_path); cat("  [SAVED]", rf_path, "\n")
  }
  
  # -------- MODEL B: XGBoost --------
  xgb_path <- file.path(checkpoint_dir, paste0(task_tag, "_XGB.rds"))
  if (file.exists(xgb_path)) {
    cat("[XGB] Loaded from checkpoint\n")
    xgb_model <- readRDS(xgb_path)
  } else {
    cat("\n[XGB] Training... ", format(Sys.time(),"%H:%M:%S"), "\n")
    ctrl_xgb <- trainControl(
      method = "cv", number = 5, index = fold_index,
      classProbs = TRUE, summaryFunction = twoClassSummary,
      savePredictions = "final", verboseIter = TRUE,
      allowParallel = FALSE)
    xgb_grid <- expand.grid(
      nrounds          = c(100, 300, 500),
      max_depth        = c(3, 5, 7),
      eta              = c(0.01, 0.05, 0.1),
      gamma            = c(0, 0.1),
      colsample_bytree = c(0.6, 0.8),
      min_child_weight = c(1, 3),
      subsample        = c(0.7, 0.9))
    cat("[XGB] Grid size:", nrow(xgb_grid), "\n")
    sample_weights <- ifelse(y_train == "disease", minority_ratio, 1)
    set.seed(42)
    xgb_model <- train(
      x = X_train, y = y_train, method = "xgbTree", metric = "ROC",
      trControl = ctrl_xgb, tuneGrid = xgb_grid,
      weights = sample_weights, verbose = 0)
    cat("[XGB] Done. Best params:\n"); print(xgb_model$bestTune)
    saveRDS(xgb_model, xgb_path); cat("  [SAVED]", xgb_path, "\n")
  }
  
  stopCluster(cl); registerDoSEQ()
  
  # -------- Recover OOF predictions --------
  rf_cv <- rf_model$pred %>% filter(mtry == rf_model$bestTune$mtry)
  xgb_cv <- xgb_model$pred %>%
    filter(nrounds          == xgb_model$bestTune$nrounds,
           max_depth        == xgb_model$bestTune$max_depth,
           eta              == xgb_model$bestTune$eta,
           gamma            == xgb_model$bestTune$gamma,
           colsample_bytree == xgb_model$bestTune$colsample_bytree,
           min_child_weight == xgb_model$bestTune$min_child_weight,
           subsample        == xgb_model$bestTune$subsample)
  
  # ----------------------------------------------------------------
  # CALIBRATION + THRESHOLDS — both derived from CV OOF only
  # ----------------------------------------------------------------
  iso_fit_rf  <- isoreg(rf_cv$disease,
                        as.integer(rf_cv$obs  == "disease"))
  iso_fit_xgb <- isoreg(xgb_cv$disease,
                        as.integer(xgb_cv$obs == "disease"))
  iso_fn      <- as.stepfun(iso_fit_rf)     # RF calibrator
  xgb_iso_fn  <- as.stepfun(iso_fit_xgb)    # XGB calibrator
  
  # Calibrate OOF preds and find Youden's J on calibrated scale
  rf_cv$disease_cal  <- pmin(pmax(iso_fn(rf_cv$disease),       0), 1)
  xgb_cv$disease_cal <- pmin(pmax(xgb_iso_fn(xgb_cv$disease),  0), 1)
  
  # Find the threshold that maximises MCC on calibrated OOF.
  # MCC balances TP/TN/FP/FN explicitly and tends to land at
  # higher thresholds than Youden, reducing false alarms.
  compute_mcc_threshold <- function(obs, probs) {
    thresholds <- seq(0.02, 0.98, by = 0.01)
    mcc_v <- numeric(length(thresholds))
    for (i in seq_along(thresholds)) {
      preds <- factor(ifelse(probs >= thresholds[i],
                             "disease", "healthy"),
                      levels = c("healthy", "disease"))
      cm <- suppressWarnings(
        confusionMatrix(preds, obs, positive = "disease"))
      TP <- cm$table["disease","disease"]; TN <- cm$table["healthy","healthy"]
      FP <- cm$table["disease","healthy"]; FN <- cm$table["healthy","disease"]
      mcc_v[i] <- calc_mcc(TP, TN, FP, FN)
    }
    mcc_v[is.na(mcc_v)] <- -Inf
    thresholds[which.max(mcc_v)]
  }
  if (grepl("CM", task_name) && !grepl("SCM", task_name)) {
    thr_rf  <- compute_mcc_threshold(rf_cv$obs,  rf_cv$disease_cal)
    thr_xgb <- compute_mcc_threshold(xgb_cv$obs, xgb_cv$disease_cal)
    cat(sprintf("Thresholds [CM | max-MCC on OOF] RF=%.2f XGB=%.2f\n",
                thr_rf, thr_xgb))
  } else {
    thr_rf  <- compute_youden_threshold(rf_cv$obs,  rf_cv$disease_cal)
    thr_xgb <- compute_youden_threshold(xgb_cv$obs, xgb_cv$disease_cal)
    cat(sprintf("Thresholds [Youden, calibrated] RF=%.2f XGB=%.2f\n",
                thr_rf, thr_xgb))
  }
  # ----------------------------------------------------------------
  # BUILD CLASSIFIER OBJECTS — model + calibrator + threshold
  # ----------------------------------------------------------------
  rf_classifier  <- make_classifier(
    model = rf_model, calibrator = iso_fn, threshold = thr_rf,
    features = all_feature_cols, algorithm_name = "RandomForest")
  xgb_classifier <- make_classifier(
    model = xgb_model, calibrator = xgb_iso_fn, threshold = thr_xgb,
    features = all_feature_cols, algorithm_name = "XGBoost")
  
  # Test-set predictions via the classifier objects
  rf_probs_raw  <- predict(rf_model,  X_test, type = "prob")[, "disease"]
  xgb_probs_raw <- predict(xgb_model, X_test, type = "prob")[, "disease"]
  rf_probs  <- rf_classifier$predict_prob(X_test)
  xgb_probs <- xgb_classifier$predict_prob(X_test)
  
  # ----------------------------------------------------------------
  # EVALUATE
  # ----------------------------------------------------------------
  evaluate <- function(probs, thr, model_name) {
    preds <- factor(ifelse(probs >= thr, "disease", "healthy"),
                    levels = c("healthy", "disease"))
    cm   <- confusionMatrix(preds, y_test, positive = "disease")
    rco  <- roc(y_test, probs, quiet = TRUE,
                levels = c("healthy", "disease"))
    TP <- cm$table["disease","disease"]; TN <- cm$table["healthy","healthy"]
    FP <- cm$table["disease","healthy"]; FN <- cm$table["healthy","disease"]
    y_bin <- as.integer(y_test == "disease")
    list(
      cm = cm, roc = rco, probs = probs, threshold = thr,
      AUC = as.numeric(auc(rco)),
      AUPRC = calc_auprc(probs, y_test),
      Sens = as.numeric(cm$byClass["Sensitivity"]),
      Spec = as.numeric(cm$byClass["Specificity"]),
      Prec = as.numeric(cm$byClass["Precision"]),
      NPV  = as.numeric(cm$byClass["Neg Pred Value"]),
      F1   = as.numeric(cm$byClass["F1"]),
      Kappa= as.numeric(cm$overall["Kappa"]),
      MCC  = calc_mcc(TP, TN, FP, FN),
      BalAcc=as.numeric(cm$byClass["Balanced Accuracy"]),
      Acc  = as.numeric(cm$overall["Accuracy"]),
      Brier= mean((probs - y_bin)^2),
      ECE  = calc_ece(probs, y_test),
      TP = TP, FN = FN, FP = FP, TN = TN,
      name = model_name)
  }
  
  res_rf  <- evaluate(rf_probs,  thr_rf,  "Random Forest")
  res_xgb <- evaluate(xgb_probs, thr_xgb, "XGBoost")
  
  cat(sprintf("\n[%s | RF ] AUC=%.3f Sens=%.3f Spec=%.3f MCC=%.3f F1=%.3f\n",
              variant_tag, res_rf$AUC, res_rf$Sens, res_rf$Spec,
              res_rf$MCC, res_rf$F1))
  cat(sprintf("[%s | XGB] AUC=%.3f Sens=%.3f Spec=%.3f MCC=%.3f F1=%.3f\n",
              variant_tag, res_xgb$AUC, res_xgb$Sens, res_xgb$Spec,
              res_xgb$MCC, res_xgb$F1))
  
  # ----------------------------------------------------------------
  # BUNDLE — classifier objects first, supporting data after
  # ----------------------------------------------------------------
  bundle <- list(
    # === PRIMARY: self-contained classifier objects =============
    rf_classifier  = rf_classifier,
    xgb_classifier = xgb_classifier,
    
    # === metadata ===============================================
    task = task_name, variant = variant_tag,
    confirmed_features = all_feature_cols,
    n_features = length(all_feature_cols),
    minority_ratio = minority_ratio,
    
    # === supporting data ========================================
    train_df = train_df, test_df = test_df,
    X_train = X_train, y_train = y_train,
    X_test  = X_test,  y_test  = y_test,
    boruta  = boruta_result,
    rf_cv = rf_cv, xgb_cv = xgb_cv,
    
    # === test predictions (raw & calibrated) ====================
    rf_probs_raw = rf_probs_raw, rf_probs = rf_probs,
    xgb_probs_raw = xgb_probs_raw, xgb_probs = xgb_probs,
    
    # === evaluation results =====================================
    res_rf = res_rf, res_xgb = res_xgb,
    
    # === components (also reachable via $rf_classifier$model etc.)
    rf_model = rf_model, xgb_model = xgb_model,
    iso_fn = iso_fn, xgb_iso_fn = xgb_iso_fn,
    thr_rf = thr_rf, thr_xgb = thr_xgb)
  
  bundle_path <- file.path(checkpoint_dir,
                           paste0(task_tag, "_bundle.rds"))
  saveRDS(bundle, bundle_path)
  cat("  [SAVED]", bundle_path, "\n")
  
  bundle
}


# ================================================================
# SECTION 5 — RUN ALL BINARY MODELS
# 3 variants × 3 tasks × 2 algorithms = 18 models total
# ================================================================

TASKS <- list(
  list(label = "Label_CM",  name = "healthy vs CM",
       short = "Clinical Mastitis (CM)"),
  list(label = "Label_SCM", name = "healthy vs SCM",
       short = "Subclinical Mastitis (SCM)"),
  list(label = "Label_Any", name = "healthy vs Any Disease",
       short = "Any Disease"))

VARIANTS <- list(
  list(tag = "ANIMAL",     feats = "features_animal_only",
       label = "Animal-only"),
  list(tag = "THERMALENV", feats = "features_thermal_env",
       label = "Thermal+Env"),
  list(tag = "SIMPLE",     feats = "features_simple",
       label = "Simple (final)"))

cat("\n", strrep("#", 70), "\n")
cat("  BINARY MODEL TRAINING\n")
cat("  Started:", format(Sys.time(), "%H:%M:%S"), "\n")
cat(strrep("#", 70), "\n")

bundles <- list()
for (vr in VARIANTS) {
  for (tsk in TASKS) {
    key <- paste(tsk$name, vr$tag, sep = " | ")
    bundles[[key]] <- run_classifier(
      df, tsk$label, tsk$name,
      feature_cols = get(vr$feats),
      variant_tag  = vr$tag)
  }
}


# ================================================================
# SECTION 6 — THREE-CLASS MODEL (healthy / CM / SCM)
# SIMPLE feature set, RF + XGB (multiclass uses argmax;
# no threshold or calibration needed)
# ================================================================

cat("\n", strrep("#", 70), "\n")
cat("  THREE-CLASS MODEL (healthy / CM / SCM)\n")
cat("  Started:", format(Sys.time(), "%H:%M:%S"), "\n")
cat(strrep("#", 70), "\n")

df3 <- df %>% filter(!is.na(Label_3class))
cat("Three-class distribution:\n"); print(table(df3$Label_3class))

set.seed(42)
uniq3   <- df3 %>% distinct(animal_no, Label_3class)
tr_idx3 <- createDataPartition(uniq3$Label_3class, p = 0.75, list = FALSE)
train_ids3 <- uniq3$animal_no[tr_idx3]
test_ids3  <- uniq3$animal_no[-tr_idx3]

train_df3 <- df3 %>% filter(animal_no %in% train_ids3)
test_df3  <- df3 %>% filter(animal_no %in% test_ids3)
cat("Train animals:", length(train_ids3),
    "| Test animals:", length(test_ids3), "\n")
cat("Train rows:", nrow(train_df3), "| Test rows:", nrow(test_df3), "\n")

boruta3_path <- file.path(checkpoint_dir, "ThreeClass_Boruta.rds")
if (file.exists(boruta3_path)) {
  cat("[Boruta 3-class] Loaded from checkpoint\n")
  boruta3 <- readRDS(boruta3_path)
} else {
  cat("[Boruta 3-class] Running...\n")
  set.seed(42)
  boruta3 <- Boruta(
    x = as.data.frame(train_df3 %>% select(all_of(features_simple))),
    y = train_df3$Label_3class,
    maxRuns = 200, pValue = 0.01, doTrace = 0)
  boruta3 <- TentativeRoughFix(boruta3)
  saveRDS(boruta3, boruta3_path)
}
confirmed3 <- names(boruta3$finalDecision)[
  boruta3$finalDecision == "Confirmed"]
if (length(confirmed3) < 2) confirmed3 <- features_simple
cat("Boruta confirmed:", length(confirmed3), "of",
    length(features_simple), "\n")

set.seed(42)
train_an3 <- train_df3 %>% mutate(.row_idx = row_number())
uniq_train_an3 <- train_an3 %>%
  distinct(animal_no, Label_3class) %>%
  group_by(Label_3class) %>%
  mutate(fold = sample(rep(1:5, length.out = n()))) %>%
  ungroup() %>%
  distinct(animal_no, .keep_all = TRUE) %>%
  select(animal_no, fold)
train_with_fold3 <- train_an3 %>%
  left_join(uniq_train_an3, by = "animal_no",
            relationship = "many-to-one")
fold_index3 <- lapply(1:5, function(k) which(train_with_fold3$fold != k))

X_train3 <- as.data.frame(train_df3 %>% select(all_of(confirmed3)))
y_train3 <- train_df3$Label_3class
X_test3  <- as.data.frame(test_df3  %>% select(all_of(confirmed3)))
y_test3  <- test_df3$Label_3class

freq3 <- table(y_train3)
class_wts3 <- as.numeric(max(freq3) / freq3); names(class_wts3) <- names(freq3)
sw3 <- class_wts3[as.character(y_train3)]

ctrl3 <- trainControl(
  method = "cv", number = 5, index = fold_index3,
  classProbs = TRUE, summaryFunction = multiClassSummary,
  savePredictions = "final", verboseIter = TRUE,
  allowParallel = TRUE)
ctrl3_seq <- ctrl3; ctrl3_seq$allowParallel <- FALSE

rf3_path <- file.path(checkpoint_dir, "ThreeClass_RF.rds")
if (file.exists(rf3_path)) {
  cat("[RF 3-class] Loaded from checkpoint\n")
  rf3 <- readRDS(rf3_path)
} else {
  cat("\n[RF 3-class] Training...\n")
  cl3 <- makeCluster(max(1, detectCores() - 1))
  registerDoParallel(cl3)
  rf_grid3 <- expand.grid(mtry = unique(pmax(1, c(
    floor(sqrt(ncol(X_train3))),
    floor(ncol(X_train3) / 4),
    floor(ncol(X_train3) / 3),
    floor(ncol(X_train3) / 2)))))
  set.seed(42)
  rf3 <- train(
    x = X_train3, y = y_train3, method = "rf",
    metric = "Mean_Balanced_Accuracy",
    trControl = ctrl3, tuneGrid = rf_grid3,
    classwt = class_wts3, ntree = 1000, importance = TRUE)
  stopCluster(cl3); registerDoSEQ()
  saveRDS(rf3, rf3_path)
}

xgb3_path <- file.path(checkpoint_dir, "ThreeClass_XGB.rds")
if (file.exists(xgb3_path)) {
  cat("[XGB 3-class] Loaded from checkpoint\n")
  xgb3 <- readRDS(xgb3_path)
} else {
  cat("\n[XGB 3-class] Training...\n")
  xgb_grid3 <- expand.grid(
    nrounds          = c(100, 300, 500),
    max_depth        = c(3, 5, 7),
    eta              = c(0.01, 0.05, 0.1),
    gamma            = c(0, 0.1),
    colsample_bytree = c(0.6, 0.8),
    min_child_weight = c(1, 3),
    subsample        = c(0.7, 0.9))
  set.seed(42)
  xgb3 <- train(
    x = X_train3, y = y_train3, method = "xgbTree",
    metric = "Mean_Balanced_Accuracy",
    trControl = ctrl3_seq, tuneGrid = xgb_grid3,
    weights = sw3, verbose = 0)
  saveRDS(xgb3, xgb3_path)
}

rf3_pred  <- predict(rf3, X_test3)
xgb3_pred <- predict(xgb3, X_test3)
rf3_cm    <- confusionMatrix(rf3_pred,  y_test3)
xgb3_cm   <- confusionMatrix(xgb3_pred, y_test3)

rf3_probs  <- predict(rf3,  X_test3, type = "prob")
xgb3_probs <- predict(xgb3, X_test3, type = "prob")

rf3_auc_ovr <- sapply(levels(y_test3), function(cl)
  as.numeric(auc(roc(as.integer(y_test3 == cl),
                     rf3_probs[, cl], quiet = TRUE))))
xgb3_auc_ovr <- sapply(levels(y_test3), function(cl)
  as.numeric(auc(roc(as.integer(y_test3 == cl),
                     xgb3_probs[, cl], quiet = TRUE))))

cat("\n[RF 3-class] CM:\n");  print(rf3_cm$table)
cat("Per-class AUC:");  print(round(rf3_auc_ovr, 3))
cat("[XGB 3-class] CM:\n"); print(xgb3_cm$table)
cat("Per-class AUC:"); print(round(xgb3_auc_ovr, 3))

bundle3 <- list(
  rf = rf3, xgb = xgb3,
  X_train = X_train3, y_train = y_train3,
  X_test = X_test3, y_test = y_test3,
  confirmed_features = confirmed3, boruta = boruta3,
  rf_pred = rf3_pred, xgb_pred = xgb3_pred,
  rf_probs = rf3_probs, xgb_probs = xgb3_probs,
  rf_cm = rf3_cm, xgb_cm = xgb3_cm,
  rf_auc_ovr = rf3_auc_ovr, xgb_auc_ovr = xgb3_auc_ovr)
saveRDS(bundle3, file.path(checkpoint_dir, "ThreeClass_bundle.rds"))


# ================================================================
# SECTION 14 — SUMMARY TABLES
# ================================================================

cat("\n[Summary tables]\n")

bin_rows <- lapply(names(bundles), function(key) {
  b <- bundles[[key]]
  parts <- strsplit(key, " \\| ")[[1]]
  task_nm <- parts[1]; var_tag <- parts[2]
  do.call(rbind, lapply(list(b$res_rf, b$res_xgb), function(r)
    data.frame(
      Task = task_nm, Feature_Set = var_tag, Model = r$name,
      N_Features = b$n_features, Threshold = round(r$threshold, 3),
      AUC = round(r$AUC, 3), AUPRC = round(r$AUPRC, 3),
      Sensitivity = round(r$Sens, 3),
      Specificity = round(r$Spec, 3),
      Precision   = round(r$Prec, 3), NPV = round(r$NPV, 3),
      F1 = round(r$F1, 3), Kappa = round(r$Kappa, 3),
      MCC = round(r$MCC, 3), Bal_Acc = round(r$BalAcc, 3),
      Accuracy = round(r$Acc, 3), Brier = round(r$Brier, 4),
      ECE = round(r$ECE, 4),
      TP_caught = r$TP, FN_missed = r$FN,
      FP_alarm = r$FP, TN_cleared = r$TN, row.names = NULL)))
})
binary_summary <- do.call(rbind, bin_rows)
write.csv(binary_summary, file.path(fig_main, "binary_metrics.csv"),
          row.names = FALSE)

build_3c_row <- function(cm, probs, y, auc_ovr, nm, nf) data.frame(
  Task = "healthy vs CM vs SCM (3-class)",
  Feature_Set = "SIMPLE", Model = nm, N_Features = nf,
  Accuracy = round(cm$overall["Accuracy"], 3),
  Kappa = round(cm$overall["Kappa"], 3),
  Macro_AUC = round(mean(auc_ovr), 3),
  AUC_healthy = round(auc_ovr["healthy"], 3),
  AUC_CM = round(auc_ovr["CM"], 3),
  AUC_SCM = round(auc_ovr["SCM"], 3),
  Sens_healthy = round(cm$byClass["Class: healthy","Sensitivity"], 3),
  Sens_CM = round(cm$byClass["Class: CM","Sensitivity"], 3),
  Sens_SCM = round(cm$byClass["Class: SCM","Sensitivity"], 3),
  Spec_healthy = round(cm$byClass["Class: healthy","Specificity"], 3),
  Spec_CM = round(cm$byClass["Class: CM","Specificity"], 3),
  Spec_SCM = round(cm$byClass["Class: SCM","Specificity"], 3),
  F1_healthy = round(cm$byClass["Class: healthy","F1"], 3),
  F1_CM = round(cm$byClass["Class: CM","F1"], 3),
  F1_SCM = round(cm$byClass["Class: SCM","F1"], 3),
  row.names = NULL)

threeclass_summary <- rbind(
  build_3c_row(rf3_cm,  rf3_probs,  y_test3, rf3_auc_ovr,
               "Random Forest", length(confirmed3)),
  build_3c_row(xgb3_cm, xgb3_probs, y_test3, xgb3_auc_ovr,
               "XGBoost",       length(confirmed3))
)
write.csv(threeclass_summary,
          file.path(fig_main, "threeclass_metrics.csv"),
          row.names = FALSE)




# ================================================================
# DONE
# ================================================================

cat("\n", strrep("#", 70), "\n  PIPELINE COMPLETE\n",
    "  Finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
    "  Trained models : ", models_dir,  "\n",
    "  Metric tables  : ", metrics_dir, "\n",
    strrep("#", 70), "\n", sep = "")

cat("\n  Headline (SIMPLE + RF) results:\n")
print(binary_summary %>%
        filter(Feature_Set == "SIMPLE", Model == "Random Forest") %>%
        select(Task, AUC, Sensitivity, Specificity, F1, MCC,
               TP_caught, FN_missed, FP_alarm, TN_cleared),
      row.names = FALSE)

cat("\n  Three-class results:\n")
print(threeclass_summary %>%
        select(Model, Accuracy, Kappa, Macro_AUC,
               AUC_healthy, AUC_CM, AUC_SCM),
      row.names = FALSE)

cat("\n  Each bundle in", models_dir, "\n",
    "  contains rf_classifier and xgb_classifier objects with\n",
    "  embedded thresholds. Load with:\n\n",
    "    bundle <- readRDS('results/models/healthy_vs_CM_SIMPLE_bundle.rds')\n",
    "    preds  <- bundle$xgb_classifier$predict_class(new_data)\n\n",
    sep = "")
