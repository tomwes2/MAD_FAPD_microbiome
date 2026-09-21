#!/usr/bin/env Rscript
# Subject-grouped RF cross-validation + label-permutation test.
#
# Why: OOB accuracy can be optimistic with paired fecal+cecal samples from the
# same adolescent. This script holds out whole subjects (LOSO) and tests whether
# observed accuracy exceeds a label-permutation null under the same CV scheme.
#
# Design choices (strength over speed, still practical):
#   - Leave-one-subject-out CV (n = 25 folds)
#   - ntree = 2000 (matches primary RF)
#   - 999 label permutations
#   - Feature stability: top-importance recurrence across LOSO folds
#
# Outputs: analysis_output/random_forest/rf_subject_grouped_cv_*
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath(".")
source(file.path(script_dir, "R", "paths.R"))
paths <- mad_init_paths()
data_dir <- paths$data_dir
output_root <- paths$output_root
setwd(data_dir)

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(randomForest)
})
source(file.path(data_dir, "R", "microbiome_helpers.R"))
source(file.path(data_dir, "R", "rf_importance_helpers.R"))

n_top_taxa <- as.integer(Sys.getenv("RF_CV_TOP_TAXA", "100"))
ntree <- as.integer(Sys.getenv("RF_CV_NTREE", "2000"))
n_perm <- as.integer(Sys.getenv("RF_CV_NPERM", "999"))
seed <- as.integer(Sys.getenv("RF_CV_SEED", "42"))
n_cores <- as.integer(Sys.getenv("RF_NCORES", "1"))
top_k_stability <- as.integer(Sys.getenv("RF_CV_TOP_K", "20"))

message("RF subject-grouped CV settings:")
message("  n_top_taxa=", n_top_taxa, " ntree=", ntree, " n_perm=", n_perm,
        " seed=", seed, " n_cores=", n_cores)

mad_objs <- load_mad_phyloseq(data_dir)
fvc_ps <- mad_objs$fvc_ps
rf_data <- prepare_rf_matrix(fvc_ps, "sample_type", n_top_taxa = n_top_taxa, add_random_feature = TRUE)
if (is.null(rf_data)) stop("Could not build RF matrix.")

meta <- data.frame(sample_data(fvc_ps), check.names = FALSE)
meta <- meta[rownames(rf_data$df), , drop = FALSE]
subjects <- as.character(meta$ssn)
df <- rf_data$df
stopifnot(nrow(df) == length(subjects))
stopifnot(all(table(subjects) == 2L))

uniq_subjects <- sort(unique(subjects))
message("Samples=", nrow(df), " subjects=", length(uniq_subjects),
        " features=", ncol(df) - 1L)

#' Leave-one-subject-out CV.
#' Returns list with overall accuracy, per-subject fold accuracies, predictions,
#' and optional top-feature ranks per fold (when return_importance=TRUE).
loso_cv <- function(df, subjects, ntree = 2000, seed = 42,
                    return_importance = FALSE, top_k = 20) {
  uniq <- sort(unique(subjects))
  preds <- rep(NA_character_, nrow(df))
  names(preds) <- rownames(df)
  truth <- as.character(df$y)
  fold_acc <- numeric(length(uniq))
  names(fold_acc) <- uniq
  top_hits <- character(0)

  for (i in seq_along(uniq)) {
    sid <- uniq[[i]]
    test_idx <- which(subjects == sid)
    train_idx <- which(subjects != sid)
    set.seed(seed + i)
    rf <- randomForest(
      y ~ .,
      data = df[train_idx, , drop = FALSE],
      ntree = ntree,
      importance = return_importance
    )
    pr <- as.character(predict(rf, df[test_idx, , drop = FALSE]))
    preds[test_idx] <- pr
    fold_acc[[sid]] <- mean(pr == truth[test_idx])

    if (return_importance) {
      imp <- as.data.frame(importance(rf, type = 1))
      imp$feature <- rownames(imp)
      imp <- imp[order(-imp$MeanDecreaseAccuracy), , drop = FALSE]
      top_hits <- c(top_hits, head(imp$feature, top_k))
    }
  }

  overall <- mean(preds == truth)
  bal <- {
    lev <- levels(df$y)
    mean(vapply(lev, function(L) {
      idx <- truth == L
      if (!any(idx)) return(NA_real_)
      mean(preds[idx] == truth[idx])
    }, numeric(1)))
  }

  out <- list(
    overall_accuracy = overall,
    balanced_accuracy = bal,
    fold_accuracy = fold_acc,
    fold_accuracy_mean = mean(fold_acc),
    fold_accuracy_sd = stats::sd(fold_acc),
    predictions = preds,
    truth = truth,
    confusion = table(truth = truth, pred = preds)
  )
  if (return_importance) {
    tab <- sort(table(top_hits), decreasing = TRUE)
    out$feature_stability <- data.frame(
      feature_safe = names(tab),
      n_folds_in_top = as.integer(tab),
      fraction_folds = as.numeric(tab) / length(uniq),
      stringsAsFactors = FALSE
    )
    out$feature_stability$taxon <- unname(rf_data$name_map[out$feature_stability$feature_safe])
  }
  out
}

# --- Observed LOSO CV ---
message("\n=== Observed leave-one-subject-out CV ===")
t0 <- Sys.time()
obs <- loso_cv(
  df, subjects,
  ntree = ntree,
  seed = seed,
  return_importance = TRUE,
  top_k = top_k_stability
)
elapsed_obs <- as.numeric(Sys.time() - t0, units = "secs")
message(sprintf(
  "Observed overall accuracy: %.1f%% | fold mean±SD: %.1f%% ± %.1f%% | balanced: %.1f%% | %.1f sec",
  100 * obs$overall_accuracy,
  100 * obs$fold_accuracy_mean,
  100 * obs$fold_accuracy_sd,
  100 * obs$balanced_accuracy,
  elapsed_obs
))
print(obs$confusion)

# Full-data OOB for side-by-side reporting (not the generalization claim)
set.seed(seed)
rf_full <- randomForest(y ~ ., data = df, ntree = ntree, importance = TRUE)
oob_acc <- 1 - rf_full$err.rate[nrow(rf_full$err.rate), "OOB"]
mtry_used <- rf_full$mtry
nodesize_used <- formals(randomForest)$nodesize
if (is.null(nodesize_used)) nodesize_used <- 1

message(sprintf("Full-data OOB accuracy (diagnostic only): %.1f%%", 100 * oob_acc))
message("mtry=", mtry_used, " (classification default floor(sqrt(p))); nodesize default=", nodesize_used)

# --- Label permutation null under identical LOSO CV ---
message("\n=== Label permutation test (n=", n_perm, ") ===")
message("Estimated wall time ~ ", round(elapsed_obs * n_perm / 60, 1), " min (serial); less if RF_NCORES>1")

run_one_perm <- function(b) {
  set.seed(seed + 10000L + b)
  df_b <- df
  df_b$y <- sample(df$y) # shuffle site labels across samples
  # Keep factor levels
  df_b$y <- factor(df_b$y, levels = levels(df$y))
  res <- loso_cv(
    df_b, subjects,
    ntree = ntree,
    seed = seed + 20000L + b,
    return_importance = FALSE
  )
  c(
    overall_accuracy = res$overall_accuracy,
    fold_accuracy_mean = res$fold_accuracy_mean,
    balanced_accuracy = res$balanced_accuracy
  )
}

perm_mat <- NULL
if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
  message("Running permutations with parallel::mclapply(mc.cores=", n_cores, ")")
  perm_list <- parallel::mclapply(seq_len(n_perm), run_one_perm, mc.cores = n_cores)
  # Filter failed
  ok <- vapply(perm_list, function(x) is.numeric(x) && length(x) == 3L, logical(1))
  if (!all(ok)) {
    warning("Some parallel jobs failed; falling back to serial for failures.")
    for (i in which(!ok)) perm_list[[i]] <- run_one_perm(i)
  }
  perm_mat <- do.call(rbind, perm_list)
} else {
  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = 3L)
  colnames(perm_mat) <- c("overall_accuracy", "fold_accuracy_mean", "balanced_accuracy")
  t_perm0 <- Sys.time()
  for (b in seq_len(n_perm)) {
    perm_mat[b, ] <- run_one_perm(b)
    if (b %% 50L == 0L || b == n_perm) {
      elapsed <- as.numeric(Sys.time() - t_perm0, units = "secs")
      eta <- elapsed / b * (n_perm - b)
      message(sprintf(
        "  perm %d/%d | null mean overall=%.3f | elapsed=%.1f min | ETA=%.1f min",
        b, n_perm, mean(perm_mat[seq_len(b), "overall_accuracy"]), elapsed / 60, eta / 60
      ))
    }
  }
}

perm_overall <- perm_mat[, "overall_accuracy"]
# Empirical one-sided p: probability null >= observed
p_overall <- (1 + sum(perm_overall >= obs$overall_accuracy)) / (1 + n_perm)
p_foldmean <- (1 + sum(perm_mat[, "fold_accuracy_mean"] >= obs$fold_accuracy_mean)) / (1 + n_perm)
p_balanced <- (1 + sum(perm_mat[, "balanced_accuracy"] >= obs$balanced_accuracy)) / (1 + n_perm)

message(sprintf(
  "Permutation p (overall accuracy): %.4g | null mean±SD: %.1f%% ± %.1f%%",
  p_overall,
  100 * mean(perm_overall),
  100 * stats::sd(perm_overall)
))

# --- Write outputs ---
out_dir <- file.path(output_root, "analysis_output", "random_forest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summary_df <- data.frame(
  metric = c(
    "loso_overall_accuracy",
    "loso_fold_accuracy_mean",
    "loso_fold_accuracy_sd",
    "loso_balanced_accuracy",
    "full_model_oob_accuracy",
    "permutation_p_overall",
    "permutation_p_fold_mean",
    "permutation_p_balanced",
    "null_overall_mean",
    "null_overall_sd",
    "n_subjects",
    "n_samples",
    "n_features_plus_noise",
    "ntree",
    "n_permutations",
    "mtry",
    "nodesize_default",
    "seed"
  ),
  value = c(
    obs$overall_accuracy,
    obs$fold_accuracy_mean,
    obs$fold_accuracy_sd,
    obs$balanced_accuracy,
    oob_acc,
    p_overall,
    p_foldmean,
    p_balanced,
    mean(perm_overall),
    stats::sd(perm_overall),
    length(uniq_subjects),
    nrow(df),
    ncol(df) - 1L,
    ntree,
    n_perm,
    mtry_used,
    as.numeric(nodesize_used),
    seed
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(out_dir, "rf_subject_grouped_cv_summary.csv"), row.names = FALSE)

fold_df <- data.frame(
  subject = names(obs$fold_accuracy),
  fold_accuracy = as.numeric(obs$fold_accuracy),
  stringsAsFactors = FALSE
)
write.csv(fold_df, file.path(out_dir, "rf_subject_grouped_cv_folds.csv"), row.names = FALSE)

pred_df <- data.frame(
  sampleID = names(obs$predictions),
  subject = subjects,
  truth = obs$truth,
  predicted = unname(obs$predictions),
  correct = obs$predictions == obs$truth,
  stringsAsFactors = FALSE
)
write.csv(pred_df, file.path(out_dir, "rf_subject_grouped_cv_predictions.csv"), row.names = FALSE)

write.csv(obs$feature_stability, file.path(out_dir, "rf_subject_grouped_cv_feature_stability.csv"), row.names = FALSE)

perm_df <- data.frame(
  permutation = seq_len(n_perm),
  overall_accuracy = perm_mat[, "overall_accuracy"],
  fold_accuracy_mean = perm_mat[, "fold_accuracy_mean"],
  balanced_accuracy = perm_mat[, "balanced_accuracy"],
  stringsAsFactors = FALSE
)
write.csv(perm_df, file.path(out_dir, "rf_subject_grouped_cv_permutation_null.csv"), row.names = FALSE)

# Simple histogram of null vs observed
png(file.path(out_dir, "rf_subject_grouped_cv_permutation_hist.png"),
    width = 8, height = 5, units = "in", res = 150)
op <- par(mar = c(5, 4, 3, 1))
hist(
  perm_overall,
  breaks = 30,
  main = "Subject-grouped RF accuracy: permutation null",
  xlab = "Leave-one-subject-out overall accuracy",
  col = "grey85",
  border = "white",
  xlim = range(c(perm_overall, obs$overall_accuracy, 0.4, 1))
)
abline(v = obs$overall_accuracy, col = "firebrick", lwd = 2)
abline(v = mean(perm_overall), col = "grey40", lwd = 1, lty = 2)
legend(
  "topright",
  legend = c(
    sprintf("Observed = %.1f%%", 100 * obs$overall_accuracy),
    sprintf("Null mean = %.1f%%", 100 * mean(perm_overall)),
    sprintf("Permutation p = %.4g", p_overall)
  ),
  bty = "n"
)
par(op)
dev.off()

# Manuscript-ready text snippet
snippet <- sprintf(
  paste0(
    "Subject-grouped leave-one-subject-out cross-validation (n = %d subjects; ntree = %d) ",
    "achieved %.1f%% overall accuracy (per-subject fold mean ± SD: %.1f%% ± %.1f%%; ",
    "balanced accuracy %.1f%%). Under %d random shuffles of GI-site labels with the same ",
    "CV scheme, mean null accuracy was %.1f%% ± %.1f%% (permutation p = %.4g). ",
    "Full-model out-of-bag accuracy was %.1f%% and is reported only as a within-forest ",
    "diagnostic, not as independent generalization performance. Hyperparameters: ",
    "mtry = %d (package default floor(sqrt(p)) for classification); nodesize = %s (package default)."
  ),
  length(uniq_subjects), ntree,
  100 * obs$overall_accuracy,
  100 * obs$fold_accuracy_mean, 100 * obs$fold_accuracy_sd,
  100 * obs$balanced_accuracy,
  n_perm,
  100 * mean(perm_overall), 100 * stats::sd(perm_overall), p_overall,
  100 * oob_acc,
  mtry_used, as.character(nodesize_used)
)
writeLines(snippet, file.path(out_dir, "rf_subject_grouped_cv_manuscript_snippet.txt"))
message("\nManuscript snippet:\n", snippet)

message("\nWrote outputs to ", out_dir)
message("Done.")
