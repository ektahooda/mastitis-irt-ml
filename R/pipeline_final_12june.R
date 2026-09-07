# ================================================================
# Cattle Disease Detection — Complete Pipeline (FINAL v2)
# ================================================================
# Murrah buffalo IRT → mastitis classification
#
# Headline model: SIMPLE feature set + Random Forest
# Variants compared: ANIMAL, THERMALENV, SIMPLE
# Algorithms     : Random Forest, XGBoost (NOT ensembled)
# Tasks          : healthy_vs_CM, healthy_vs_SCM, healthy_vs_Any
#                  + three-class (healthy / CM / SCM)
#
# KEY DESIGN POINTS (changed from v1):
#
# 1. THRESHOLD AS PART OF THE MODEL
#    Each trained algorithm is packaged as a self-contained
#    classifier object: (caret_model + isotonic_calibrator +
#    threshold + feature_list + predict methods). The threshold
#    is a model attribute derived from training data only — not
#    an external post-hoc choice. See make_classifier().
#
# 2. CALIBRATION — BOTH MODELS
#    RF and XGB are BOTH isotonically recalibrated on CV out-of-
#    fold predictions. Decision thresholds (Youden's J) are
#    computed on the CALIBRATED OOF scale so they apply
#    consistently to calibrated test probabilities.
#
# 3. LACTATION STAGE — 3 CATEGORIES
#    Early (DIM ≤ 90), Mid (91–180), Late (>180).
#
# 4. ANIMAL-AWARE SPLITTING THROUGHOUT
#    No records from one animal in both train and test
#    (Bobbo et al. 2023, J Dairy Sci).
#
# 5. NO ENSEMBLE
#    RF and XGB are evaluated and reported separately.
#
# Outputs:
#   <checkpoint_dir>/manuscript_main/   main text figures + tables
#   <checkpoint_dir>/manuscript_supp/   supplementary figures
#   <checkpoint_dir>/*.rds              model checkpoints
# ================================================================


# ================================================================
# SECTION 0 — SETUP
# ================================================================

# Ensure xgboost 1.7.x is installed (caret xgbTree incompatible with 3.x).
# Run this ONCE then comment out:
# install.packages("C:/Users/Admin/Downloads/xgboost_1.7.9.1.zip",
#                  repos = NULL, type = "binary")

packages <- c("readxl","dplyr","lubridate","caret","randomForest",
              "xgboost","pROC","ggplot2","tibble",
              "doParallel","tidyr","foreach","matrixStats",
              "gridExtra","MLmetrics","reshape2")
installed <- packages %in% rownames(installed.packages())
if (any(!installed)) install.packages(packages[!installed])
invisible(lapply(packages, library, character.only = TRUE))

if (!"Boruta" %in% rownames(installed.packages())) {
  install.packages("Boruta", type = "binary")
}
library(Boruta)

if (!"fastshap" %in% rownames(installed.packages())) {
  install.packages("fastshap")
}
library(fastshap)

cat("xgboost version:", as.character(packageVersion("xgboost")),
    "(must be < 3.0)\n")

# ---------------- USER PATHS — EDIT IF NEEDED ----------------
data_path      <- "D://CNN//complete data 22sept.xlsx"
checkpoint_dir <- "D://CNN//checkpoints_final"
dir.create(checkpoint_dir, showWarnings = FALSE, recursive = TRUE)

fig_main <- file.path(checkpoint_dir, "manuscript_main")
fig_supp <- file.path(checkpoint_dir, "manuscript_supp")
dir.create(fig_main, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_supp, showWarnings = FALSE, recursive = TRUE)

cat("Checkpoints :", checkpoint_dir, "\n")
cat("Main figs   :", fig_main, "\n")
cat("Supp figs   :", fig_supp, "\n")


# ================================================================
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
# SECTION 7 — PLOTTING HELPERS
# ================================================================

plot_cm_binary <- function(cm_obj, title, subtitle = "",
                           colour_high = "#08519c") {
  df_cm <- as.data.frame(cm_obj$table)
  names(df_cm)[names(df_cm) == "Freq"] <- "Count"
  pred_low <- tolower(as.character(df_cm$Prediction))
  ref_low  <- tolower(as.character(df_cm$Reference))
  df_cm$Predicted <- factor(
    ifelse(pred_low == "disease", "Disease", "Healthy"),
    levels = c("Healthy", "Disease"))
  df_cm$Actual <- factor(
    ifelse(ref_low == "disease", "Disease", "Healthy"),
    levels = c("Disease", "Healthy"))
  df_cm$Label <- mapply(function(p, r) {
    if (p == "healthy" && r == "healthy") "TN"
    else if (p == "disease" && r == "healthy") "FP"
    else if (p == "healthy" && r == "disease") "FN"
    else "TP"
  }, pred_low, ref_low)
  df_cm$Pct <- 100 * df_cm$Count / sum(df_cm$Count)
  ggplot(df_cm, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%s\n%d\n(%.1f%%)",
                                  Label, Count, Pct)),
              size = 5.5, fontface = "bold", color = "white",
              lineheight = 0.95) +
    scale_fill_gradient(low = "#9ecae1", high = colour_high) +
    labs(title = title, subtitle = subtitle,
         x = "Predicted class", y = "Actual class") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none",
          panel.grid     = element_blank(),
          axis.text      = element_text(face = "bold", size = 12),
          plot.title     = element_text(face = "bold", size = 14),
          plot.subtitle  = element_text(size = 11, colour = "grey30"))
}

plot_cm_3class <- function(cm_obj, title, subtitle = "",
                           colour_high = "#08519c") {
  df_cm <- as.data.frame(cm_obj$table)
  names(df_cm) <- c("Predicted","Actual","Count")
  df_cm$Pct <- 100 * df_cm$Count / sum(df_cm$Count)
  df_cm$Predicted <- factor(df_cm$Predicted,
                            levels = rev(levels(df_cm$Predicted)))
  ggplot(df_cm, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Pct)),
              size = 5, fontface = "bold", color = "white",
              lineheight = 0.95) +
    scale_fill_gradient(low = "#deebf7", high = colour_high) +
    labs(title = title, subtitle = subtitle,
         x = "Actual class", y = "Predicted class") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          panel.grid     = element_blank(),
          axis.text      = element_text(face = "bold", size = 12),
          plot.title     = element_text(face = "bold", size = 14),
          plot.subtitle  = element_text(size = 11, colour = "grey30"))
}

plot_roc_comparison <- function(roc_list, title) {
  dfs <- lapply(names(roc_list), function(nm) {
    r <- roc_list[[nm]]
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities,
               Model = sprintf("%s (AUC=%.3f)", nm, as.numeric(auc(r))))
  })
  rdf <- do.call(rbind, dfs)
  ggplot(rdf, aes(x = FPR, y = TPR, colour = Model)) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey60") +
    labs(title = title, x = "False Positive Rate (1 - Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    theme_minimal(base_size = 13) +
    theme(legend.position = c(0.65, 0.25),
          legend.background = element_rect(fill = "white", colour = NA),
          legend.title = element_blank(),
          plot.title = element_text(face = "bold", size = 14))
}

plot_calibration <- function(probs, y_test, title, n_bins = 10) {
  y_bin <- as.integer(y_test == "disease")
  brks  <- seq(0, 1, length.out = n_bins + 1)
  bins  <- cut(probs, brks, include.lowest = TRUE)
  cal <- data.frame(probs = probs, y = y_bin, bin = bins) %>%
    group_by(bin) %>%
    summarise(mean_pred = mean(probs), mean_obs = mean(y),
              n = n(), .groups = "drop") %>%
    filter(n > 0)
  brier <- mean((probs - y_bin)^2)
  ece   <- calc_ece(probs, y_test, n_bins)
  ggplot(cal, aes(x = mean_pred, y = mean_obs)) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey60") +
    geom_point(aes(size = n), colour = "#08519c", alpha = 0.85) +
    geom_line(colour = "#08519c") +
    scale_size_continuous(range = c(2, 8)) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title,
         subtitle = sprintf("Brier=%.4f | ECE=%.4f", brier, ece),
         x = "Mean predicted probability",
         y = "Observed event rate", size = "Bin size") +
    theme_minimal(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 11, colour = "grey30"))
}

plot_boruta_importance <- function(boruta_obj, title, top_n = 30) {
  imp <- attStats(boruta_obj)
  imp$Feature  <- rownames(imp)
  imp <- imp[!grepl("^shadow", imp$Feature), ]
  imp$Decision <- factor(imp$decision,
                         levels = c("Confirmed","Tentative","Rejected"))
  imp <- imp[order(-imp$meanImp), ]
  if (nrow(imp) > top_n) imp <- imp[1:top_n, ]
  imp$Feature <- factor(imp$Feature, levels = rev(imp$Feature))
  ggplot(imp, aes(x = meanImp, y = Feature, fill = Decision)) +
    geom_col() +
    geom_errorbarh(aes(xmin = minImp, xmax = maxImp), height = 0.25,
                   colour = "grey40") +
    scale_fill_manual(values = c(Confirmed = "#08519c",
                                 Tentative = "#fdae6b",
                                 Rejected  = "grey70")) +
    labs(title = title, x = "Boruta importance", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

plot_perm_importance <- function(perm_df, title, top_n = 25) {
  perm_df <- perm_df[order(-perm_df$Importance), ]
  if (nrow(perm_df) > top_n) perm_df <- perm_df[1:top_n, ]
  perm_df$Feature <- factor(perm_df$Feature, levels = rev(perm_df$Feature))
  ggplot(perm_df, aes(x = Importance, y = Feature)) +
    geom_col(fill = "#08519c") +
    geom_errorbarh(aes(xmin = Importance - SD, xmax = Importance + SD),
                   height = 0.25, colour = "grey40") +
    labs(title = title, x = "Mean AUC drop (permutation importance)",
         y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

plot_shap_summary <- function(shap_mat, X_test, title, top_n = 15) {
  shap_long <- as.data.frame(shap_mat) %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(-row_id, names_to = "Feature", values_to = "SHAP")
  x_long <- as.data.frame(X_test) %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(-row_id, names_to = "Feature", values_to = "Value")
  shap_df <- shap_long %>% left_join(x_long, by = c("row_id", "Feature")) %>%
    group_by(Feature) %>%
    mutate(Value_norm = (Value - min(Value, na.rm = TRUE)) /
             (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE) + 1e-9)) %>%
    ungroup()
  feat_order <- shap_df %>%
    group_by(Feature) %>%
    summarise(MAE = mean(abs(SHAP), na.rm = TRUE)) %>%
    arrange(-MAE) %>% head(top_n) %>% pull(Feature)
  shap_df <- shap_df %>% filter(Feature %in% feat_order)
  shap_df$Feature <- factor(shap_df$Feature, levels = rev(feat_order))
  ggplot(shap_df, aes(x = SHAP, y = Feature, colour = Value_norm)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_jitter(height = 0.25, alpha = 0.6, size = 1.4) +
    scale_colour_gradient(low = "#3182bd", high = "#e6550d",
                          name = "Feature\nvalue",
                          breaks = c(0, 1), labels = c("Low", "High")) +
    labs(title = title,
         x = "SHAP value (impact on disease probability)", y = NULL) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}


# ================================================================
# SECTION 8 — CONFUSION MATRIX FIGURES
# ================================================================

cat("\n", strrep("#", 70), "\n  CONFUSION MATRIX FIGURES\n",
    strrep("#", 70), "\n")

# MAIN: SIMPLE + RF for each task
main_cms <- lapply(TASKS, function(tsk) {
  b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
  r <- b$res_rf
  st <- sprintf("AUC=%.3f | Sens=%.3f | Spec=%.3f | MCC=%.3f",
                r$AUC, r$Sens, r$Spec, r$MCC)
  plot_cm_binary(r$cm, title = tsk$short, subtitle = st)
})
g_main <- gridExtra::arrangeGrob(grobs = main_cms, ncol = 3,
                                 top = "Random Forest")
ggsave(file.path(fig_main, "Fig_CM_SIMPLE_RF.png"),
       g_main, width = 16, height = 5.5, dpi = 300, bg = "white")
cat("  [SAVED] Fig_CM_SIMPLE_RF.png\n")

# SUPP: every variant per task
for (tsk in TASKS) {
  panels <- list()
  for (vr in VARIANTS) {
    b <- bundles[[paste(tsk$name, vr$tag, sep = " | ")]]
    for (mt in c("RF","XGB")) {
      r <- if (mt == "RF") b$res_rf else b$res_xgb
      mtxt <- if (mt == "RF") "Random Forest" else "XGBoost"
      st <- sprintf("AUC=%.3f | Sens=%.3f | Spec=%.3f | MCC=%.3f",
                    r$AUC, r$Sens, r$Spec, r$MCC)
      panels[[length(panels) + 1]] <- plot_cm_binary(
        r$cm, title = paste0(vr$label, " — ", mtxt), subtitle = st)
    }
  }
  g <- gridExtra::arrangeGrob(grobs = panels, ncol = 2,
                              top = paste0("Confusion matrices — ", tsk$short))
  ggsave(file.path(fig_supp,
                   paste0("FigS_CM_", gsub(" ","_",tsk$name), ".png")),
         g, width = 11, height = 14, dpi = 300, bg = "white")
}

# THREE-CLASS — main + supp
best3_is_rf <- rf3_cm$overall["Kappa"] >= xgb3_cm$overall["Kappa"]
best3_cm <- if (best3_is_rf) rf3_cm else xgb3_cm
best3_nm <- if (best3_is_rf) "Random Forest" else "XGBoost"
best3_au <- if (best3_is_rf) rf3_auc_ovr else xgb3_auc_ovr

p3_main <- plot_cm_3class(best3_cm,
                          title = paste0("Three-class — ", best3_nm),
                          subtitle = sprintf("Accuracy=%.3f | Kappa=%.3f | Macro AUC=%.3f",
                                             best3_cm$overall["Accuracy"],
                                             best3_cm$overall["Kappa"], mean(best3_au)))
ggsave(file.path(fig_main, "Fig_CM_ThreeClass_best.png"),
       p3_main, width = 6.5, height = 6, dpi = 300, bg = "white")

p3_rf  <- plot_cm_3class(rf3_cm,  "Three-class — Random Forest",
                         sprintf("Accuracy=%.3f | Kappa=%.3f | Macro AUC=%.3f",
                                 rf3_cm$overall["Accuracy"], rf3_cm$overall["Kappa"],
                                 mean(rf3_auc_ovr)))
p3_xgb <- plot_cm_3class(xgb3_cm, "Three-class — XGBoost",
                         sprintf("Accuracy=%.3f | Kappa=%.3f | Macro AUC=%.3f",
                                 xgb3_cm$overall["Accuracy"], xgb3_cm$overall["Kappa"],
                                 mean(xgb3_auc_ovr)))
g3 <- gridExtra::arrangeGrob(p3_rf, p3_xgb, ncol = 2)
ggsave(file.path(fig_supp, "FigS_CM_ThreeClass.png"),
       g3, width = 12, height = 6, dpi = 300, bg = "white")


# ================================================================
# SECTION 9 — ROC CURVES
# ================================================================

cat("\n[ROC curves]\n")
for (tsk in TASKS) {
  roc_all <- list()
  for (vr in VARIANTS) {
    b <- bundles[[paste(tsk$name, vr$tag, sep = " | ")]]
    roc_all[[paste(vr$label, "| RF")]]  <- b$res_rf$roc
    roc_all[[paste(vr$label, "| XGB")]] <- b$res_xgb$roc
  }
  p_all <- plot_roc_comparison(roc_all,
                               paste0("ROC curves — ", tsk$short))
  ggsave(file.path(fig_supp,
                   paste0("FigS_ROC_", gsub(" ","_",tsk$name), ".png")),
         p_all, width = 7.5, height = 6, dpi = 300, bg = "white")
  roc_main <- roc_all[grepl("Simple", names(roc_all))]
  p_main <- plot_roc_comparison(roc_main,
                                paste0("ROC — ", tsk$short, "Full Model"))
  ggsave(file.path(fig_main,
                   paste0("Fig_ROC_", gsub(" ","_",tsk$name), ".png")),
         p_main, width = 6, height = 5.5, dpi = 300, bg = "white")
}


# ================================================================
# SECTION 10 — CALIBRATION PLOTS (4 panels per task)
# ================================================================

cat("\n[Calibration plots]\n")
for (tsk in TASKS) {
  b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
  p_rf_raw  <- plot_calibration(b$rf_probs_raw,  b$y_test, "RF (uncalibrated)")
  p_rf_cal  <- plot_calibration(b$rf_probs,      b$y_test, "RF (isotonic-calibrated)")
  p_xgb_raw <- plot_calibration(b$xgb_probs_raw, b$y_test, "XGB (uncalibrated)")
  p_xgb_cal <- plot_calibration(b$xgb_probs,     b$y_test, "XGB (isotonic-calibrated)")
  g_cal <- gridExtra::arrangeGrob(p_rf_raw, p_rf_cal, p_xgb_raw, p_xgb_cal,
                                  ncol = 2, top = paste0("Calibration — ", tsk$short))
  ggsave(file.path(fig_supp,
                   paste0("FigS_Cal_", gsub(" ","_",tsk$name), ".png")),
         g_cal, width = 11, height = 9, dpi = 300, bg = "white")
}


# ================================================================
# SECTION 11 — BORUTA IMPORTANCE PLOTS (SIMPLE)
# ================================================================

cat("\n[Boruta importance plots]\n")
for (tsk in TASKS) {
  b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
  p <- plot_boruta_importance(b$boruta,
                              title = paste0("Boruta importance — ", tsk$short, " (SIMPLE)"))
  ggsave(file.path(fig_supp,
                   paste0("FigS_Boruta_", gsub(" ","_",tsk$name), ".png")),
         p, width = 8, height = 9, dpi = 300, bg = "white")
}
p3b <- plot_boruta_importance(boruta3,
                              title = "Boruta importance — Three-class model (SIMPLE)")
ggsave(file.path(fig_supp, "FigS_Boruta_ThreeClass.png"),
       p3b, width = 8, height = 9, dpi = 300, bg = "white")


# ================================================================
# SECTION 12 — PERMUTATION IMPORTANCE (SIMPLE + RF)
# Uses the classifier object's calibrator automatically
# ================================================================

cat("\n[Permutation importance — SIMPLE+RF]\n")

permute_importance <- function(classifier, X_test, y_test,
                               n_reps = 20, seed = 42) {
  feat_names <- classifier$features
  set.seed(seed)
  base_probs <- classifier$predict_prob(X_test)
  base_auc <- as.numeric(auc(roc(y_test, base_probs, quiet = TRUE,
                                 levels = c("healthy","disease"))))
  out <- data.frame(Feature = feat_names,
                    Importance = NA_real_, SD = NA_real_)
  for (i in seq_along(feat_names)) {
    aucs <- numeric(n_reps)
    for (r in 1:n_reps) {
      X_perm <- X_test
      X_perm[[feat_names[i]]] <- sample(X_perm[[feat_names[i]]])
      probs <- classifier$predict_prob(X_perm)
      aucs[r] <- as.numeric(auc(roc(y_test, probs, quiet = TRUE,
                                    levels = c("healthy","disease"))))
    }
    out$Importance[i] <- base_auc - mean(aucs)
    out$SD[i]         <- sd(aucs)
  }
  out
}

for (tsk in TASKS) {
  ckpt_path <- file.path(checkpoint_dir,
                         paste0(gsub(" ","_",tsk$name),
                                "_SIMPLE_PermImp.rds"))
  if (file.exists(ckpt_path)) {
    perm_df <- readRDS(ckpt_path)
  } else {
    cat("  ", tsk$name, "... ")
    b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
    perm_df <- permute_importance(b$rf_classifier, b$X_test, b$y_test,
                                  n_reps = 20)
    saveRDS(perm_df, ckpt_path)
    cat("done\n")
  }
  p <- plot_perm_importance(perm_df,
                            paste0("Permutation importance — ", tsk$short, " | SIMPLE+RF"))
  ggsave(file.path(fig_main,
                   paste0("Fig_PermImp_", gsub(" ","_",tsk$name), ".png")),
         p, width = 8, height = 8, dpi = 300, bg = "white")
}


# ================================================================
# SECTION 13 — SHAP (SIMPLE + RF via fastshap)
# ================================================================

cat("\n[SHAP values — SIMPLE+RF via fastshap]\n")

predict_rf_prob <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, "disease"]
}

for (tsk in TASKS) {
  shap_path <- file.path(checkpoint_dir,
                         paste0(gsub(" ","_",tsk$name),
                                "_SIMPLE_SHAP.rds"))
  if (file.exists(shap_path)) {
    shap_mat <- readRDS(shap_path)
  } else {
    cat("  Computing SHAP for", tsk$name, "...\n")
    b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
    set.seed(42)
    shap_mat <- fastshap::explain(
      object = b$rf_classifier$model$finalModel,
      X = b$X_train, newdata = b$X_test,
      pred_wrapper = predict_rf_prob, nsim = 30)
    saveRDS(shap_mat, shap_path)
  }
  b <- bundles[[paste(tsk$name, "SIMPLE", sep = " | ")]]
  p <- plot_shap_summary(shap_mat, b$X_test,
                         paste0("SHAP summary — ", tsk$short, " | SIMPLE+RF"))
  ggsave(file.path(fig_main,
                   paste0("Fig_SHAP_", gsub(" ","_",tsk$name), ".png")),
         p, width = 9, height = 7, dpi = 300, bg = "white")
}


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
    "  Finished:", format(Sys.time(), "%H:%M:%S"), "\n",
    "  Main figures :", fig_main, "\n",
    "  Supp figures :", fig_supp, "\n",
    strrep("#", 70), "\n")

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

cat("\n  NOTE: Each task's bundle contains rf_classifier and\n")
cat("  xgb_classifier objects with embedded thresholds. For\n")
cat("  prediction on new data:\n")
cat("    b <- readRDS('healthy_vs_CM_SIMPLE_bundle.rds')\n")
cat("    predicted_classes <- b$rf_classifier$predict_class(new_data)\n")

library(Boruta)

tasks <- list(
  list(name = "healthy_vs_CM",          label = "CM"),
  list(name = "healthy_vs_SCM",         label = "SCM"),
  list(name = "healthy_vs_Any_Disease", label = "Any Disease")
)

boruta_summary <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_Boruta.rds")))
  dec <- data.frame(Feature  = names(b$finalDecision),
                    Decision = as.character(b$finalDecision))
  confirmed <- dec$Feature[dec$Decision == "Confirmed"]
  rejected  <- dec$Feature[dec$Decision == "Rejected"]
  cat(sprintf("\n=== %s ===  Confirmed=%d  Rejected=%d\n",
              tsk$label, length(confirmed), length(rejected)))
  cat("Confirmed:\n  ", paste(confirmed, collapse = ", "), "\n")
  cat("Rejected :\n  ", paste(rejected,  collapse = ", "), "\n")
  boruta_summary[[tsk$label]] <- list(confirmed = confirmed,
                                      rejected  = rejected)
}

# Features rejected across ALL three tasks
all_rej <- Reduce(intersect,
                  lapply(boruta_summary, function(x) x$rejected))
cat("\nRejected by Boruta in ALL three tasks:\n  ",
    paste(all_rej, collapse = ", "), "\n")

# 3-class Boruta
b3 <- readRDS(file.path(checkpoint_dir, "ThreeClass_Boruta.rds"))
dec3 <- data.frame(Feature = names(b3$finalDecision),
                   Decision = as.character(b3$finalDecision))
cat("\n=== Three-class ===  Confirmed=", 
    sum(dec3$Decision == "Confirmed"),
    "  Rejected=", sum(dec3$Decision == "Rejected"), "\n", sep = "")

cat("\n=== PERMUTATION IMPORTANCE — top 10 per task (SIMPLE+RF) ===\n")
for (tsk in tasks) {
  perm_df <- readRDS(file.path(checkpoint_dir,
                               paste0(tsk$name, "_SIMPLE_PermImp.rds")))
  top10 <- perm_df[order(-perm_df$Importance), ][1:10, ]
  top10$Importance <- round(top10$Importance, 4)
  top10$SD         <- round(top10$SD, 4)
  cat(sprintf("\n--- %s ---\n", tsk$label))
  print(top10, row.names = FALSE)
}
cat("\n=== SHAP — top 10 features by mean |SHAP| per task ===\n")
for (tsk in tasks) {
  shap_mat <- readRDS(file.path(checkpoint_dir,
                                paste0(tsk$name, "_SIMPLE_SHAP.rds")))
  s_summary <- data.frame(
    Feature       = colnames(shap_mat),
    Mean_abs_SHAP = colMeans(abs(shap_mat)))
  s_summary <- s_summary[order(-s_summary$Mean_abs_SHAP), ][1:10, ]
  s_summary$Mean_abs_SHAP <- round(s_summary$Mean_abs_SHAP, 4)
  cat(sprintf("\n--- %s ---\n", tsk$label))
  print(s_summary, row.names = FALSE)
}
library(pROC); library(caret)

bootstrap_metrics <- function(y_test, probs, thr,
                              n_boot = 2000, seed = 42) {
  set.seed(seed)
  n <- length(y_test)
  out <- vector("list", n_boot)
  for (i in 1:n_boot) {
    idx <- sample.int(n, n, replace = TRUE)
    y_b <- y_test[idx]; p_b <- probs[idx]
    if (length(unique(y_b)) < 2) next
    preds <- factor(ifelse(p_b >= thr, "disease", "healthy"),
                    levels = c("healthy", "disease"))
    cm <- suppressWarnings(
      confusionMatrix(preds, y_b, positive = "disease"))
    rco <- suppressMessages(
      roc(y_b, p_b, quiet = TRUE, levels = c("healthy","disease")))
    TP <- cm$table["disease","disease"]
    TN <- cm$table["healthy","healthy"]
    FP <- cm$table["disease","healthy"]
    FN <- cm$table["healthy","disease"]
    out[[i]] <- c(
      AUC      = as.numeric(auc(rco)),
      Sens     = as.numeric(cm$byClass["Sensitivity"]),
      Spec     = as.numeric(cm$byClass["Specificity"]),
      Prec     = as.numeric(cm$byClass["Precision"]),
      F1       = as.numeric(cm$byClass["F1"]),
      Accuracy = as.numeric(cm$overall["Accuracy"]),
      MCC      = calc_mcc(TP, TN, FP, FN))
  }
  m <- do.call(rbind, out[!sapply(out, is.null)])
  data.frame(
    Metric    = colnames(m),
    Estimate  = round(apply(m, 2, mean,     na.rm = TRUE), 3),
    Lower2.5  = round(apply(m, 2, quantile, 0.025, na.rm = TRUE), 3),
    Upper97.5 = round(apply(m, 2, quantile, 0.975, na.rm = TRUE), 3))
}

# Both RF and XGB for each task
ci_all <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_bundle.rds")))
  for (mtype in c("RF", "XGB")) {
    probs <- if (mtype == "RF") b$rf_probs else b$xgb_probs
    thr   <- if (mtype == "RF") b$thr_rf  else b$thr_xgb
    cat(sprintf("\n=== %s | %s | threshold=%.2f ===\n",
                tsk$label, mtype, thr))
    ci <- bootstrap_metrics(b$y_test, probs, thr, n_boot = 2000)
    print(ci, row.names = FALSE)
    ci$Task  <- tsk$label
    ci$Model <- mtype
    ci_all[[paste(tsk$label, mtype, sep = "_")]] <- ci
  }
}

all_ci_df <- do.call(rbind, ci_all)
write.csv(all_ci_df, file.path(fig_main, "bootstrap_CIs_binary.csv"),
          row.names = FALSE)
cat("\nSaved: bootstrap_CIs_binary.csv\n")

b3 <- readRDS(file.path(checkpoint_dir, "ThreeClass_bundle.rds"))

bootstrap_3class <- function(y_test, pred_classes, probs,
                             n_boot = 2000, seed = 42) {
  set.seed(seed)
  n <- length(y_test)
  classes <- levels(y_test)
  out <- vector("list", n_boot)
  for (i in 1:n_boot) {
    idx <- sample.int(n, n, replace = TRUE)
    y_b <- y_test[idx]; pc_b <- pred_classes[idx]; pp_b <- probs[idx, ]
    if (length(unique(y_b)) < 3) next
    cm <- suppressWarnings(confusionMatrix(pc_b, y_b))
    per_class_auc <- sapply(classes, function(cl) {
      as.numeric(auc(roc(as.integer(y_b == cl), pp_b[, cl],
                         quiet = TRUE)))
    })
    out[[i]] <- c(
      Accuracy    = as.numeric(cm$overall["Accuracy"]),
      Kappa       = as.numeric(cm$overall["Kappa"]),
      Macro_AUC   = mean(per_class_auc),
      AUC_healthy = per_class_auc["healthy"],
      AUC_CM      = per_class_auc["CM"],
      AUC_SCM     = per_class_auc["SCM"],
      F1_healthy  = cm$byClass["Class: healthy", "F1"],
      F1_CM       = cm$byClass["Class: CM",      "F1"],
      F1_SCM      = cm$byClass["Class: SCM",     "F1"])
  }
  m <- do.call(rbind, out[!sapply(out, is.null)])
  data.frame(
    Metric    = colnames(m),
    Estimate  = round(apply(m, 2, mean,     na.rm = TRUE), 3),
    Lower2.5  = round(apply(m, 2, quantile, 0.025, na.rm = TRUE), 3),
    Upper97.5 = round(apply(m, 2, quantile, 0.975, na.rm = TRUE), 3))
}

cat("\n=== Three-class | Random Forest ===\n")
ci_rf3 <- bootstrap_3class(b3$y_test, b3$rf_pred, b3$rf_probs)
print(ci_rf3, row.names = FALSE)

cat("\n=== Three-class | XGBoost ===\n")
ci_xgb3 <- bootstrap_3class(b3$y_test, b3$xgb_pred, b3$xgb_probs)
print(ci_xgb3, row.names = FALSE)

ci_rf3$Model  <- "RandomForest"
ci_xgb3$Model <- "XGBoost"
write.csv(rbind(ci_rf3, ci_xgb3),
          file.path(fig_main, "bootstrap_CIs_threeclass.csv"),
          row.names = FALSE)

library(ggplot2); library(dplyr); library(tidyr); library(Boruta)

tasks <- list(
  list(name = "healthy_vs_CM",          label = "CM"),
  list(name = "healthy_vs_SCM",         label = "SCM"),
  list(name = "healthy_vs_Any_Disease", label = "Any")
)

# Helper: rank → score in [0,1] (top rank = 1, rank > 20 floored at 0)
rank_to_score <- function(rank, max_meaningful = 20) {
  pmax(0, 1 - (rank - 1) / max_meaningful)
}

rows <- list()
for (tsk in tasks) {
  # Boruta: 1 for Confirmed, 0 for Rejected, 0.5 for Tentative
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_Boruta.rds")))
  bor <- data.frame(
    Feature = names(b$finalDecision),
    Score   = case_when(
      as.character(b$finalDecision) == "Confirmed" ~ 1,
      as.character(b$finalDecision) == "Tentative" ~ 0.5,
      TRUE ~ 0),
    Method = "Boruta", Task = tsk$label,
    stringsAsFactors = FALSE)
  
  # Permutation: rank → score
  perm <- readRDS(file.path(checkpoint_dir,
                            paste0(tsk$name, "_SIMPLE_PermImp.rds")))
  perm <- perm[order(-perm$Importance), ]
  perm$Rank <- seq_len(nrow(perm))
  perm <- data.frame(
    Feature = perm$Feature, Score = rank_to_score(perm$Rank),
    Method = "Permutation", Task = tsk$label,
    stringsAsFactors = FALSE)
  
  # SHAP: rank by mean|SHAP| → score
  shap_mat <- readRDS(file.path(checkpoint_dir,
                                paste0(tsk$name, "_SIMPLE_SHAP.rds")))
  shap <- data.frame(
    Feature = colnames(shap_mat),
    MeanAbs = colMeans(abs(shap_mat)))
  shap <- shap[order(-shap$MeanAbs), ]
  shap$Rank <- seq_len(nrow(shap))
  shap <- data.frame(
    Feature = shap$Feature, Score = rank_to_score(shap$Rank),
    Method = "SHAP", Task = tsk$label,
    stringsAsFactors = FALSE)
  
  rows <- c(rows, list(bor, perm, shap))
}

all_long <- bind_rows(rows)

# Overall consistency: sum of Score across all 9 cells per feature
feat_order <- all_long %>%
  group_by(Feature) %>%
  summarise(TotalScore = sum(Score, na.rm = TRUE), .groups = "drop") %>%
  arrange(-TotalScore)

top_n <- 20
top_features <- feat_order$Feature[seq_len(min(top_n, nrow(feat_order)))]

# Build plotting frame
# Fill missing cells with score=0 — features rejected by Boruta
# for a given task were not evaluated in permutation/SHAP, but
# their effective importance is near zero by construction.
all_long <- expand.grid(
  Feature = unique(all_long$Feature),
  Method  = c("Boruta", "Permutation", "SHAP"),
  Task    = c("CM", "SCM", "Any"),
  stringsAsFactors = FALSE) %>%
  left_join(all_long, by = c("Feature", "Method", "Task")) %>%
  mutate(Score = ifelse(is.na(Score), 0, Score))
plot_df <- all_long %>%
  filter(Feature %in% top_features) %>%
  mutate(
    Feature = factor(Feature, levels = rev(top_features)),
    Method  = factor(Method,  levels = c("Boruta", "Permutation", "SHAP")),
    Task    = factor(Task,    levels = c("CM", "SCM", "Any")),
    Column  = paste(Method, Task, sep = "\n"))

col_levels <- c(
  paste("Boruta",      c("CM","SCM","Any"), sep = "\n"),
  paste("Permutation", c("CM","SCM","Any"), sep = "\n"),
  paste("SHAP",        c("CM","SCM","Any"), sep = "\n"))
plot_df$Column <- factor(plot_df$Column, levels = col_levels)

p <- ggplot(plot_df, aes(x = Column, y = Feature, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_vline(xintercept = c(3.5, 6.5),
             color = "black", linewidth = 0.8) +
  scale_fill_gradient2(
    low      = "#d73027", mid    = "#fee090", high    = "#4575b4",
    midpoint = 0.5,
    breaks   = c(0, 0.5, 1),
    labels   = c("Low", "Mid", "High"),
    name     = NULL,
    limits   = c(0, 1)) +
  annotate("text", x = 2, y = top_n + 1.0,
           label = "Boruta",      fontface = "bold", size = 4.5) +
  annotate("text", x = 5, y = top_n + 1.0,
           label = "Permutation", fontface = "bold", size = 4.5) +
  annotate("text", x = 8, y = top_n + 1.0,
           label = "SHAP",        fontface = "bold", size = 4.5) +
  labs(x = NULL, y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x       = element_text(size = 9),
    axis.text.y       = element_text(size = 10),
    panel.grid        = element_blank(),
    plot.margin       = margin(t = 25, r = 10, b = 10, l = 10),
    legend.position   = "right")

print(p)

ggsave(file.path(fig_supp, "Fig_FeatureConsistency_Heatmap.png"),
       p, width = 10, height = 9, dpi = 300, bg = "white")

cat("\n[SAVED] Fig_FeatureConsistency_Heatmap.png\n")
cat("\nTop-", top_n, " features by cross-method consistency:\n", sep = "")
print(feat_order[1:top_n, ], row.names = FALSE)

library(dplyr)
library(tidyr)
library(Boruta)

tasks <- list(
  list(name = "healthy_vs_CM",          label = "CM"),
  list(name = "healthy_vs_SCM",         label = "SCM"),
  list(name = "healthy_vs_Any_Disease", label = "Any")
)

all_rankings <- list()

for (tsk in tasks) {
  # ---- BORUTA — full ranking by mean importance ----
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_Boruta.rds")))
  bor <- attStats(b)
  bor$Feature <- rownames(bor)
  bor <- bor[!grepl("^shadow", bor$Feature), ]
  bor <- bor[order(-bor$meanImp), ]
  bor$Boruta_rank <- seq_len(nrow(bor))
  bor <- bor[, c("Feature", "meanImp", "decision", "Boruta_rank")]
  names(bor) <- c("Feature", "Boruta_meanImp", "Boruta_decision", "Boruta_rank")
  
  # ---- PERMUTATION — full ranking ----
  perm <- readRDS(file.path(checkpoint_dir,
                            paste0(tsk$name, "_SIMPLE_PermImp.rds")))
  perm <- perm[order(-perm$Importance), ]
  perm$Perm_rank <- seq_len(nrow(perm))
  perm <- perm[, c("Feature", "Importance", "Perm_rank")]
  names(perm) <- c("Feature", "Perm_importance", "Perm_rank")
  
  # ---- SHAP — full ranking by mean|SHAP| ----
  shap_mat <- readRDS(file.path(checkpoint_dir,
                                paste0(tsk$name, "_SIMPLE_SHAP.rds")))
  shap <- data.frame(
    Feature      = colnames(shap_mat),
    SHAP_meanAbs = colMeans(abs(shap_mat)))
  shap <- shap[order(-shap$SHAP_meanAbs), ]
  shap$SHAP_rank <- seq_len(nrow(shap))
  
  # Outer join — keep all features even if rejected/unmeasured
  merged <- bor %>%
    full_join(perm, by = "Feature") %>%
    full_join(shap, by = "Feature")
  merged$Task <- tsk$label
  all_rankings[[tsk$label]] <- merged
  
  # ---- Print per-task summary ----
  cat("\n=================================\n  TASK:", tsk$label, "\n=================================\n")
  
  cat("\n[BORUTA] top 10 by meanImp:\n")
  print(merged %>% arrange(Boruta_rank) %>% head(10) %>%
          select(Feature, Boruta_meanImp, Boruta_decision),
        row.names = FALSE)
  
  cat("\n[PERMUTATION] full ranking:\n")
  print(merged %>% filter(!is.na(Perm_rank)) %>% arrange(Perm_rank) %>%
          select(Feature, Perm_importance, Perm_rank),
        row.names = FALSE)
  
  cat("\n[SHAP] full ranking:\n")
  print(merged %>% filter(!is.na(SHAP_rank)) %>% arrange(SHAP_rank) %>%
          select(Feature, SHAP_meanAbs, SHAP_rank),
        row.names = FALSE)
}

all_long <- bind_rows(all_rankings)


# ================================================================
# CROSS-METHOD CONSISTENCY SCORING
# Score each feature 0-9: confirmed by Boruta (3 tasks)
#                       + in top 10 by Permutation (3 tasks)
#                       + in top 10 by SHAP (3 tasks)
# ================================================================
consistency <- all_long %>%
  group_by(Feature) %>%
  summarise(
    Boruta_confirmed_tasks = sum(Boruta_decision == "Confirmed", na.rm = TRUE),
    Boruta_rejected_tasks  = sum(Boruta_decision == "Rejected",  na.rm = TRUE),
    Perm_top10_tasks       = sum(!is.na(Perm_rank) & Perm_rank <= 10),
    SHAP_top10_tasks       = sum(!is.na(SHAP_rank) & SHAP_rank <= 10),
    .groups = "drop") %>%
  mutate(ConsistencyScore = Boruta_confirmed_tasks +
           Perm_top10_tasks + SHAP_top10_tasks) %>%
  arrange(-ConsistencyScore)

cat("\n=================================\n  CROSS-METHOD CONSISTENCY\n=================================\n")
cat("\nScore = #Boruta-confirmed tasks + #Permutation top-10 tasks + #SHAP top-10 tasks\n")
cat("Maximum = 9 (confirmed by all methods in all 3 tasks)\n\n")
print(as.data.frame(consistency), row.names = FALSE)

# Category labels
fully_consistent <- consistency$Feature[consistency$ConsistencyScore == 9]
mostly_consistent <- consistency$Feature[consistency$ConsistencyScore %in% 6:8]
rejected_everywhere <- consistency$Feature[
  consistency$Boruta_rejected_tasks == 3 &
    consistency$Perm_top10_tasks      == 0 &
    consistency$SHAP_top10_tasks      == 0]

cat("\nFully consistent (score = 9):\n  ",
    paste(fully_consistent, collapse = ", "), "\n")
cat("\nMostly consistent (score 6-8):\n  ",
    paste(mostly_consistent, collapse = ", "), "\n")
cat("\nRejected by Boruta in all tasks AND not top-10 in either Perm or SHAP:\n  ",
    paste(rejected_everywhere, collapse = ", "), "\n")


# ================================================================
# DIVERGENCES — Permutation vs SHAP
# Find features where ranks differ substantially per task
# ================================================================
cat("\n=================================\n  PERMUTATION vs SHAP DIVERGENCES\n=================================\n")
for (tsk in tasks) {
  m <- all_rankings[[tsk$label]]
  m <- m %>% filter(!is.na(Perm_rank) & !is.na(SHAP_rank))
  m$Rank_diff <- m$SHAP_rank - m$Perm_rank   # positive = Perm ranks higher
  
  cat("\n--- ", tsk$label,
      " (positive = Perm ranks feature higher than SHAP) ---\n", sep = "")
  
  cat("\nLargest positive divergences (Perm thinks it matters, SHAP disagrees):\n")
  print(m %>% arrange(-Rank_diff) %>% head(5) %>%
          select(Feature, Perm_rank, SHAP_rank, Rank_diff),
        row.names = FALSE)
  
  cat("\nLargest negative divergences (SHAP thinks it matters, Perm disagrees):\n")
  print(m %>% arrange(Rank_diff) %>% head(5) %>%
          select(Feature, Perm_rank, SHAP_rank, Rank_diff),
        row.names = FALSE)
}

# Save for reference / supplementary
write.csv(all_long,
          file.path(fig_main, "feature_importance_full.csv"),
          row.names = FALSE)
write.csv(consistency,
          file.path(fig_main, "feature_consistency_scores.csv"),
          row.names = FALSE)
cat("\n[SAVED] feature_importance_full.csv, feature_consistency_scores.csv\n")


#################################
# ================================================================
# Model Validation and Comparison — generates all numbers
# needed for the manuscript paragraph.
# Requires: pROC, caret, and the helper calc_ece() from the
# pipeline. Reads saved bundles from checkpoint_dir.
# ================================================================
library(pROC); library(caret); library(dplyr)

tasks <- list(
  list(name = "healthy_vs_CM",          label = "CM"),
  list(name = "healthy_vs_SCM",         label = "SCM"),
  list(name = "healthy_vs_Any_Disease", label = "Any Disease"))

# Best model per task (matches your headline)
best_model <- c(CM = "XGB", SCM = "RF", `Any Disease` = "RF")

# ================================================================
# 1. DeLong's test (RF vs XGB AUC comparison) per task
# ================================================================
cat("\n=================================\n  1. DeLong's test (RF vs XGB AUC)\n=================================\n")

delong_results <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_bundle.rds")))
  roc_rf  <- roc(b$y_test, b$rf_probs,  quiet = TRUE,
                 levels = c("healthy","disease"))
  roc_xgb <- roc(b$y_test, b$xgb_probs, quiet = TRUE,
                 levels = c("healthy","disease"))
  d <- roc.test(roc_rf, roc_xgb, method = "delong")
  cat(sprintf("\n[%s] AUC_RF=%.3f  AUC_XGB=%.3f  Z=%.3f  p=%.3f\n",
              tsk$label, as.numeric(auc(roc_rf)),
              as.numeric(auc(roc_xgb)),
              d$statistic, d$p.value))
  delong_results[[tsk$label]] <- data.frame(
    Task     = tsk$label,
    AUC_RF   = round(as.numeric(auc(roc_rf)),  3),
    AUC_XGB  = round(as.numeric(auc(roc_xgb)), 3),
    Z        = round(as.numeric(d$statistic), 3),
    p_value  = signif(d$p.value, 3))
}

# ================================================================
# 2. McNemar's test (RF vs XGB error patterns) per task
# ================================================================
cat("\n=================================\n  2. McNemar's test (RF vs XGB errors)\n=================================\n")

mcnemar_results <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_bundle.rds")))
  rf_pred  <- factor(ifelse(b$rf_probs  >= b$thr_rf,  "disease", "healthy"),
                     levels = c("healthy","disease"))
  xgb_pred <- factor(ifelse(b$xgb_probs >= b$thr_xgb, "disease", "healthy"),
                     levels = c("healthy","disease"))
  rf_correct  <- as.integer(rf_pred  == b$y_test)
  xgb_correct <- as.integer(xgb_pred == b$y_test)
  tab <- table(RF = rf_correct, XGB = xgb_correct)
  # Need a 2x2 table
  if (all(dim(tab) == c(2,2))) {
    m <- mcnemar.test(tab, correct = TRUE)
    cat(sprintf("\n[%s] McNemar's chi-squared=%.3f  df=%d  p=%.3f\n",
                tsk$label, as.numeric(m$statistic),
                as.numeric(m$parameter), m$p.value))
    cat(" Discordance: RF-correct/XGB-wrong=", tab["1","0"], 
        " | RF-wrong/XGB-correct=", tab["0","1"], "\n", sep = "")
    mcnemar_results[[tsk$label]] <- data.frame(
      Task         = tsk$label,
      chi_squared  = round(as.numeric(m$statistic), 3),
      df           = as.numeric(m$parameter),
      p_value      = signif(m$p.value, 3),
      RF_only_correct  = tab["1","0"],
      XGB_only_correct = tab["0","1"])
  } else {
    cat(sprintf("\n[%s] McNemar: insufficient discordance to compute\n",
                tsk$label))
  }
}

# ================================================================
# 3. Permutation test (model significance vs label-shuffled null)
#    For each task, shuffle test labels n_perm times and compute
#    chance-level AUC distribution; observed AUC's p-value is
#    the fraction of shuffled AUCs >= observed.
# ================================================================
cat("\n=================================\n  3. Permutation test (model vs chance, ",
    "best model per task)\n=================================\n", sep = "")

n_perm <- 500
perm_results <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_bundle.rds")))
  bm <- best_model[tsk$label]
  probs <- if (bm == "RF") b$rf_probs else b$xgb_probs
  
  obs_auc <- as.numeric(auc(roc(b$y_test, probs, quiet = TRUE,
                                levels = c("healthy","disease"))))
  set.seed(42)
  null_auc <- numeric(n_perm)
  for (i in 1:n_perm) {
    y_shuf <- sample(b$y_test)
    null_auc[i] <- as.numeric(auc(roc(y_shuf, probs, quiet = TRUE,
                                      levels = c("healthy","disease"))))
  }
  p_perm <- (sum(null_auc >= obs_auc) + 1) / (n_perm + 1)
  cat(sprintf("\n[%s | %s] Observed AUC=%.3f  null mean=%.3f  null max=%.3f  p<%s\n",
              tsk$label, bm, obs_auc, mean(null_auc), max(null_auc),
              ifelse(p_perm < 1/n_perm, sprintf("%.4f", 1/n_perm),
                     sprintf("%.4f", p_perm))))
  perm_results[[tsk$label]] <- data.frame(
    Task        = tsk$label,
    Model       = bm,
    Observed_AUC = round(obs_auc, 3),
    Null_mean   = round(mean(null_auc), 3),
    Null_max    = round(max(null_auc), 3),
    p_value     = p_perm)
}

# ================================================================
# 4. Calibration: Brier and ECE — raw and calibrated
#    For BOTH RF and XGB (since you calibrated both)
# ================================================================
cat("\n=================================\n  4. Calibration — Brier and ECE (raw vs calibrated)\n=================================\n")

calib_results <- list()
for (tsk in tasks) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(tsk$name, "_SIMPLE_bundle.rds")))
  y_bin <- as.integer(b$y_test == "disease")
  
  brier_rf_raw  <- mean((b$rf_probs_raw  - y_bin)^2)
  brier_rf_cal  <- mean((b$rf_probs      - y_bin)^2)
  brier_xgb_raw <- mean((b$xgb_probs_raw - y_bin)^2)
  brier_xgb_cal <- mean((b$xgb_probs     - y_bin)^2)
  
  ece_rf_raw  <- calc_ece(b$rf_probs_raw,  b$y_test)
  ece_rf_cal  <- calc_ece(b$rf_probs,      b$y_test)
  ece_xgb_raw <- calc_ece(b$xgb_probs_raw, b$y_test)
  ece_xgb_cal <- calc_ece(b$xgb_probs,     b$y_test)
  
  cat(sprintf("\n[%s]\n", tsk$label))
  cat(sprintf("  RF  Brier: raw=%.4f calibrated=%.4f  |  ECE: raw=%.4f calibrated=%.4f\n",
              brier_rf_raw, brier_rf_cal, ece_rf_raw, ece_rf_cal))
  cat(sprintf("  XGB Brier: raw=%.4f calibrated=%.4f  |  ECE: raw=%.4f calibrated=%.4f\n",
              brier_xgb_raw, brier_xgb_cal, ece_xgb_raw, ece_xgb_cal))
  
  calib_results[[tsk$label]] <- data.frame(
    Task = tsk$label,
    Brier_RF_raw  = round(brier_rf_raw,  4),
    Brier_RF_cal  = round(brier_rf_cal,  4),
    Brier_XGB_raw = round(brier_xgb_raw, 4),
    Brier_XGB_cal = round(brier_xgb_cal, 4),
    ECE_RF_raw    = round(ece_rf_raw,    4),
    ECE_RF_cal    = round(ece_rf_cal,    4),
    ECE_XGB_raw   = round(ece_xgb_raw,   4),
    ECE_XGB_cal   = round(ece_xgb_cal,   4))
}

# ================================================================
# 5. Three-class permutation test (against macro AUC null)
# ================================================================
cat("\n=================================\n  5. Three-class permutation test\n=================================\n")

b3 <- readRDS(file.path(checkpoint_dir, "ThreeClass_bundle.rds"))
classes <- levels(b3$y_test)

macro_auc_obs <- mean(b3$rf_auc_ovr)

set.seed(42)
null_macro <- numeric(n_perm)
for (i in 1:n_perm) {
  y_shuf <- sample(b3$y_test)
  per_class <- sapply(classes, function(cl) {
    as.numeric(auc(roc(as.integer(y_shuf == cl),
                       b3$rf_probs[, cl], quiet = TRUE)))
  })
  null_macro[i] <- mean(per_class)
}
p_perm_3c <- (sum(null_macro >= macro_auc_obs) + 1) / (n_perm + 1)
cat(sprintf("\n[3-class | RF] Observed macro AUC=%.3f  null mean=%.3f  null max=%.3f  p<%s\n",
            macro_auc_obs, mean(null_macro), max(null_macro),
            ifelse(p_perm_3c < 1/n_perm, sprintf("%.4f", 1/n_perm),
                   sprintf("%.4f", p_perm_3c))))

# ================================================================
# Combine and save
# ================================================================
write.csv(do.call(rbind, delong_results),
          file.path(fig_main, "delong_test.csv"), row.names = FALSE)
write.csv(do.call(rbind, mcnemar_results),
          file.path(fig_main, "mcnemar_test.csv"), row.names = FALSE)
write.csv(do.call(rbind, perm_results),
          file.path(fig_main, "permutation_test.csv"), row.names = FALSE)
write.csv(do.call(rbind, calib_results),
          file.path(fig_main, "calibration_summary.csv"), row.names = FALSE)

cat("\n=================================\n  Saved:\n")
cat("  delong_test.csv         — DeLong AUC comparison\n")
cat("  mcnemar_test.csv        — McNemar error-pattern test\n")
cat("  permutation_test.csv    — Permutation significance test\n")
cat("  calibration_summary.csv — Brier and ECE, raw and calibrated\n")
cat("=================================\n")





library(ggplot2); library(gridExtra); library(dplyr)

# ----------------------------------------------------------------
# Plot helpers (inline — works even if pipeline not in session)
# ----------------------------------------------------------------
plot_cm_binary <- function(cm_obj, title, subtitle = "",
                           colour_high = "#08519c") {
  df_cm <- as.data.frame(cm_obj$table)
  names(df_cm)[names(df_cm) == "Freq"] <- "Count"
  pred_low <- tolower(as.character(df_cm$Prediction))
  ref_low  <- tolower(as.character(df_cm$Reference))
  df_cm$Predicted <- factor(
    ifelse(pred_low == "disease", "Disease", "Healthy"),
    levels = c("Healthy", "Disease"))
  df_cm$Actual <- factor(
    ifelse(ref_low == "disease", "Disease", "Healthy"),
    levels = c("Disease", "Healthy"))
  df_cm$Label <- mapply(function(p, r) {
    if (p == "healthy" && r == "healthy") "TN"
    else if (p == "disease" && r == "healthy") "FP"
    else if (p == "healthy" && r == "disease") "FN"
    else "TP"
  }, pred_low, ref_low)
  df_cm$Pct <- 100 * df_cm$Count / sum(df_cm$Count)
  ggplot(df_cm, aes(x = Predicted, y = Actual, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%s\n%d\n(%.1f%%)",
                                  Label, Count, Pct)),
              size = 5, fontface = "bold", color = "white",
              lineheight = 0.95) +
    scale_fill_gradient(low = "#9ecae1", high = colour_high) +
    labs(title = title, subtitle = subtitle,
         x = "Predicted class", y = "Actual class") +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.text       = element_text(face = "bold", size = 11),
      axis.title      = element_text(face = "bold", size = 12),
      
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5   # center title
      ),
      
      plot.subtitle = element_text(
        size = 10,
        colour = "grey30",
        hjust = 0.5   # center subtitle
      )
    )
}

plot_cm_3class <- function(cm_obj, title, subtitle = "",
                           colour_high = "#08519c") {
  df_cm <- as.data.frame(cm_obj$table)
  names(df_cm) <- c("Predicted","Actual","Count")
  df_cm$Pct <- 100 * df_cm$Count / sum(df_cm$Count)
  df_cm$Predicted <- factor(df_cm$Predicted,
                            levels = rev(levels(df_cm$Predicted)))
  ggplot(df_cm, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white", linewidth = 1.5) +
    geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Pct)),
              size = 4.5, fontface = "bold", color = "white",
              lineheight = 0.95) +
    scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
    labs(title = title, subtitle = subtitle,
         x = "Actual class", y = "Predicted class") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none",
          panel.grid     = element_blank(),
          axis.text      = element_text(face = "bold", size = 11),
          plot.title     = element_text(face = "bold", size = 13),
          plot.subtitle  = element_text(size = 10, colour = "grey30"))
}

# ================================================================
# BUILD PANELS — best model per task
# ================================================================

# Binary: CM=XGB, SCM=RF, Any=RF
best_binary <- list(
  list(task  = "healthy_vs_CM",
       short = "Clinical Mastitis", model = "XGB"),
  list(task  = "healthy_vs_SCM",
       short = "Subclinical Mastitis", model = "RF"),
  list(task  = "healthy_vs_Any_Disease",
       short = "Any Mastitis", model = "RF"))

binary_panels <- lapply(best_binary, function(bm) {
  
  b <- readRDS(file.path(
    checkpoint_dir,
    paste0(bm$task, "_SIMPLE_bundle.rds")
  ))
  
  r <- if (bm$model == "RF") b$res_rf else b$res_xgb
  
  plot_cm_binary(
    r$cm,
    title = bm$short,
    subtitle = NULL
  )
})

# Three-class: RF and XGB
b3 <- readRDS(file.path(checkpoint_dir, "ThreeClass_bundle.rds"))

threeclass_panels <- list(
  
  plot_cm_3class(
    b3$rf_cm,
    title = "Random Forest",
    subtitle = NULL
  ),
  
  plot_cm_3class(
    b3$xgb_cm,
    title = "XGBoost",
    subtitle = NULL
  )
)
  plot_cm_3class(b3$xgb_cm,
                 title    = "Three-class — XGBoost",
                 subtitle = sprintf("Accuracy=%.3f | Kappa=%.3f | Macro AUC=%.3f",
                                    b3$xgb_cm$overall["Accuracy"],
                                    b3$xgb_cm$overall["Kappa"],
                                    mean(b3$xgb_auc_ovr))))
                                    mean(b3$xgb_auc_ovr))))

# ================================================================
# COMBINE — top row 3 binary, bottom row 2 three-class
# Three-class panels naturally wider since they have a 3×3 grid
# ================================================================

top_row    <- arrangeGrob(grobs = binary_panels,     ncol = 3,
                          top = "Binary tasks — best model per task")
bottom_row <- arrangeGrob(grobs = threeclass_panels, ncol = 2,
                          top = "Three-class model (healthy / CM / SCM)")

combined <- arrangeGrob(top_row, bottom_row,
                        ncol = 1, heights = c(1, 1.1))

out_path <- file.path(fig_main, "Fig_CM_BestModels_Combined.png")
ggsave(out_path, combined,
       width = 15, height = 11, dpi = 300, bg = "white")

cat("[SAVED]", out_path, "\n")
ggsave(file.path(fig_main, "Fig_CM_BestModels_Binary.png"),
       top_row,    width = 16, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(fig_main, "Fig_CM_BestModels_ThreeClass.png"),
       bottom_row, width = 12, height = 6,   dpi = 300, bg = "white")



library(ggplot2); library(pROC); library(dplyr)

# Best model per task (matches the confusion matrix panel)
best_binary <- list(
  list(task  = "healthy_vs_CM",
       short = "Clinical Mastitis (CM)",         model = "XGB",
       model_full = "XGBoost"),
  list(task  = "healthy_vs_SCM",
       short = "Subclinical Mastitis (SCM)",    model = "RF",
       model_full = "Random Forest"),
  list(task  = "healthy_vs_Any_Disease",
       short = "Any Disease (CM + SCM)",        model = "RF",
       model_full = "Random Forest"))

# Collect ROC curve data
dfs <- lapply(best_binary, function(bm) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(bm$task, "_SIMPLE_bundle.rds")))
  r <- if (bm$model == "RF") b$res_rf else b$res_xgb
  data.frame(
    FPR   = 1 - r$roc$specificities,
    TPR   = r$roc$sensitivities,
    Model = sprintf("%s — %s (AUC = %.3f)",
                    bm$short, bm$model_full, r$AUC))
})
rdf <- do.call(rbind, dfs)

# Plot
p <- ggplot(rdf, aes(x = FPR, y = TPR, colour = Model)) +
  geom_line(linewidth = 1.3) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = c("#1b9e77", "#d95f02", "#7570b3")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title    = "ROC curves — best model per task",
       subtitle = "Binary classification on held-out test set",
       x = "False Positive Rate (1 − Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position    = c(0.62, 0.18),
    legend.background  = element_rect(fill = "white", colour = "grey80"),
    legend.title       = element_blank(),
    legend.text        = element_text(size = 10),
    plot.title         = element_text(face = "bold", size = 14),
    plot.subtitle      = element_text(size = 11, colour = "grey30"),
    panel.grid.minor   = element_blank())

out_path <- file.path(fig_main, "Fig_ROC_BestModels_Combined.png")
ggsave(out_path, p, width = 8, height = 7, dpi = 300, bg = "white")

cat("[SAVED]", out_path, "\n")
library(gridExtra)

panels <- lapply(best_binary, function(bm) {
  b <- readRDS(file.path(checkpoint_dir,
                         paste0(bm$task, "_SIMPLE_bundle.rds")))
  best_lbl <- if (bm$model == "RF") "RF (best)" else "XGB (best)"
  other_lbl <- if (bm$model == "RF") "XGB"        else "RF"
  
  rdf_task <- bind_rows(
    data.frame(FPR = 1 - b$res_rf$roc$specificities,
               TPR = b$res_rf$roc$sensitivities,
               Model = sprintf("Random Forest (AUC=%.3f)%s",
                               b$res_rf$AUC,
                               if (bm$model == "RF") " — primary" else "")),
    data.frame(FPR = 1 - b$res_xgb$roc$specificities,
               TPR = b$res_xgb$roc$sensitivities,
               Model = sprintf("XGBoost (AUC=%.3f)%s",
                               b$res_xgb$AUC,
                               if (bm$model == "XGB") " — primary" else "")))
  
  ggplot(rdf_task, aes(x = FPR, y = TPR, colour = Model)) +
    geom_line(linewidth = 1.1) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "dashed", colour = "grey60") +
    scale_colour_manual(values = c("#1f78b4", "#e31a1c")) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = bm$short,
         x = "1 − Specificity", y = "Sensitivity") +
    theme_minimal(base_size = 12) +
    theme(legend.position = c(0.62, 0.18),
          legend.title    = element_blank(),
          legend.text     = element_text(size = 9),
          legend.background = element_rect(fill = "white", colour = "grey80"),
          plot.title      = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank())
})

g <- arrangeGrob(grobs = panels, ncol = 3,
                 top = "ROC curves — Random Forest vs XGBoost per task")
ggsave(file.path(fig_main, "Fig_ROC_Panel_AllModels.png"),
       g, width = 16, height = 5.5, dpi = 300, bg = "white")

# ================================================================
# Animal-averaged SHAP importance across the three binary tasks
# Replicates the original figure style with SIMPLE features
# ================================================================

library(dplyr); library(tidyr); library(ggplot2)

# ----------------------------------------------------------------
# Step 1 — Inspect SHAP bundle structure (run once to confirm fields)
# ----------------------------------------------------------------
shap_probe <- readRDS(file.path(checkpoint_dir,
                                "healthy_vs_CM_SIMPLE_SHAP.rds"))
str(shap_probe, max.level = 1)
# Expected fields: rf_shap (matrix), xgb_shap (matrix),
#                   animal_no_test (vector) — adjust below if different

# ----------------------------------------------------------------
# ----------------------------------------------------------------
# Step 2 (REPLACEMENT) — Load and summarise SHAP per task
# ----------------------------------------------------------------
# Step 2 (FINAL) — Load SHAP, pull animal IDs from bundle$test_df
# ----------------------------------------------------------------
load_shap_summary <- function(task_name, variant = "SIMPLE") {
  
  shap_path   <- file.path(checkpoint_dir,
                           paste0(task_name, "_", variant, "_SHAP.rds"))
  bundle_path <- file.path(checkpoint_dir,
                           paste0(task_name, "_", variant, "_bundle.rds"))
  
  shap_obj <- readRDS(shap_path)
  shap_mat <- as.matrix(shap_obj)
  
  bundle     <- readRDS(bundle_path)
  animal_ids <- bundle$test_df$animal_no
  
  # Sanity check
  if (length(animal_ids) != nrow(shap_mat)) {
    cat("[", task_name, "] WARNING: animal IDs (", length(animal_ids),
        ") and SHAP rows (", nrow(shap_mat), ") differ — falling back to record-level\n")
    feature_importance <- colMeans(abs(shap_mat))
  } else {
    df <- as.data.frame(abs(shap_mat)) %>%
      mutate(animal_no = animal_ids) %>%
      group_by(animal_no) %>%
      summarise(across(everything(), mean), .groups = "drop") %>%
      select(-animal_no)
    feature_importance <- colMeans(df)
    cat("[", task_name, "] Animal-averaged across",
        length(unique(animal_ids)), "animals\n")
  }
  
  data.frame(
    Feature       = colnames(shap_mat),
    Mean_abs_SHAP = as.numeric(feature_importance),
    Task          = task_name,
    stringsAsFactors = FALSE)
}

# ----------------------------------------------------------------
# Step 3 — Apply loader to all three tasks (this is the lapply)
# ----------------------------------------------------------------
tasks <- list(
  list(task = "healthy_vs_CM",          label = "CM"),
  list(task = "healthy_vs_SCM",         label = "SCM"),
  list(task = "healthy_vs_Any_Disease", label = "Any"))

shap_all <- bind_rows(lapply(tasks, function(t) {
  df <- load_shap_summary(t$task)
  df$Task <- t$label
  df
}))

cat("\nFeatures found:", length(unique(shap_all$Feature)), "\n")
cat("Rows in shap_all:", nrow(shap_all), "(should be ~3 × features)\n")

# ----------------------------------------------------------------
# Step 4 — Plot
# ----------------------------------------------------------------
p <- ggplot(plot_data,
            aes(x = Mean_abs_SHAP, y = Feature, fill = Task)) +
  geom_col(position = position_dodge(width = 0.8),
           width    = 0.75,
           colour   = NA) +
  scale_fill_manual(values = c("Any" = "#4daf4a",
                               "CM"  = "#e41a1c",
                               "SCM" = "#377eb8")) +
  labs(title    = "Animal-averaged SHAP importance across tasks",
       subtitle = "Mean |SHAP value| per animal — fastshap, RF",
       x = "Mean |SHAP value|",
       y = NULL,
       fill = "Task") +
  theme_minimal(base_size = 12) +
  theme(legend.position    = "bottom",
        plot.title         = element_text(face = "bold", size = 14),
        plot.subtitle      = element_text(colour = "grey30", size = 11),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 10),
        axis.title.x       = element_text(size = 11))

# ----------------------------------------------------------------
# Step 5 — Save (PNG for viewing + TIFF for publication)
# ----------------------------------------------------------------
out_png  <- file.path(fig_supp, "SuppFig_SHAP_animal_averaged.png")
out_tiff <- file.path(fig_supp, "SuppFig_SHAP_animal_averaged.tiff")

ggsave(out_png,  p, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(out_tiff, p, width = 9, height = 7, dpi = 300, bg = "white",
       compression = "lzw")

cat("[SAVED]", out_png, "\n")
cat("[SAVED]", out_tiff, "\n")

b <- readRDS(file.path(checkpoint_dir, "healthy_vs_CM_SIMPLE_bundle.rds"))
colnames(b$test_df)
