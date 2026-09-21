#!/usr/bin/env Rscript
# Standalone robust RF importance analysis (permutation, drop-column, collinearity).
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
  library(ggplot2)
  library(ComplexHeatmap)
  library(viridisLite)
})

source(file.path(data_dir, "R", "microbiome_helpers.R"))
source(file.path(data_dir, "R", "rf_importance_helpers.R"))

mad_objs <- load_mad_phyloseq(data_dir)
fvc_ps <- mad_objs$fvc_ps
rf_dir <- file.path(output_root, "analysis_output", "random_forest")
dir.create(rf_dir, recursive = TRUE, showWarnings = FALSE)

n_cores <- as.integer(Sys.getenv("RF_NCORES", "1"))
message("Running RF deep analysis (n_cores=", n_cores, ")...")

rf_deep <- run_rf_deep_analysis(
  ps = fvc_ps,
  response_var = "sample_type",
  label = "pain_fecal_vs_cecal",
  out_dir = rf_dir,
  n_top_taxa = 100,
  drop_col_top_n = 30,
  cor_top_n = 40,
  ntree_full = 2000,
  ntree_fast = 500,
  perm_repeats = 10,
  cor_threshold = 0.7,
  n_cores = n_cores
)

if (is.null(rf_deep)) {
  stop("RF deep analysis returned NULL — check sample_type levels in pain cohort.")
}

message("OOB accuracy: ", round(100 * rf_deep$summary$oob_accuracy, 1), "%")
message("Random feature rank: ", rf_deep$summary$random_feature_rank, " / ", nrow(rf_deep$perm_se))
message("Outputs written to: ", rf_dir)
