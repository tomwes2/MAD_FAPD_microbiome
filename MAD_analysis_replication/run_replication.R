#!/usr/bin/env Rscript
# Paper 1 replication: paired fecal vs cecal (pain cohort only).
# Inputs: MAD_DATA_DIR (default: this kit). Outputs: MAD_OUTPUT_DIR.
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath(".")
source(file.path(script_dir, "R", "paths.R"))
paths <- mad_init_paths()
data_dir <- paths$data_dir
output_root <- paths$output_root
setwd(data_dir)
message("data_dir:    ", data_dir)
message("output_root: ", output_root)

run_cmd <- function(script) {
  message("\n=== ", script, " ===")
  status <- system2("Rscript", script, stdout = "", stderr = "")
  if (!identical(status, 0L)) stop(script, " failed with exit code ", status)
}

# 1) Taxonomy: pain cohort stool vs cecum
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  stop("Install rmarkdown to render microbiome_pain_fecal_vs_cecal.Rmd")
}
reports_dir <- file.path(output_root, "reports")
dir.create(reports_dir, recursive = TRUE, showWarnings = FALSE)
message("\n=== microbiome_pain_fecal_vs_cecal.Rmd ===")
rmarkdown::render(
  "microbiome_pain_fecal_vs_cecal.Rmd",
  output_dir = reports_dir,
  intermediates_dir = file.path(reports_dir, "microbiome_pain_fecal_vs_cecal_tmp"),
  knit_root_dir = data_dir,
  quiet = FALSE
)

# Combine PERMANOVA tables if present
analysis_dir <- file.path(output_root, "analysis_output")
perm_parts <- file.path(analysis_dir, c("permanova_pain_fecal_vs_cecal.csv", "permanova_fecal.csv"))
perm_parts <- perm_parts[file.exists(perm_parts)]
if (length(perm_parts) && requireNamespace("dplyr", quietly = TRUE)) {
  combined <- dplyr::bind_rows(lapply(perm_parts, read.csv, stringsAsFactors = FALSE))
  write.csv(combined, file.path(analysis_dir, "permanova_results.csv"), row.names = FALSE)
}

# 2) Pathways, RF, network
run_cmd("run_functional_pathway_analysis.R")
run_cmd("annotate_direction_outputs.R")
run_cmd("run_rf_importance_analysis.R")
run_cmd("run_rf_subject_grouped_cv.R")
run_cmd("run_pathway_network_analysis.R")

message("\nDone.")
message("Taxonomy:     ", file.path(output_root, "analysis_output/"))
message("Pathways:     ", file.path(output_root, "functional_analysis_output/"))
message("Network:      ", file.path(output_root, "functional_analysis_output/pathway_network/"))
message("RF outputs:   ", file.path(output_root, "analysis_output/random_forest/"))
