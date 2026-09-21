#!/usr/bin/env Rscript
# Post-process existing analysis CSV/TSV files to add direction columns.
# Run after MaAsLin/Wilcoxon (no need to re-fit models):
#   Rscript annotate_direction_outputs.R

suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath(".")
source(file.path(script_dir, "R", "paths.R"))
paths <- mad_init_paths()
data_dir <- paths$data_dir
output_root <- paths$output_root
setwd(data_dir)
source(file.path(data_dir, "R", "direction_helpers.R"))

annotate_csv <- function(path, fun, ...) {
  if (!file.exists(path)) return(invisible(NULL))
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  df <- fun(df, ...)
  write.csv(df, path, row.names = FALSE)
  message("Annotated: ", path)
}

annotate_tsv <- function(path, fun, ...) {
  if (!file.exists(path)) return(invisible(NULL))
  df <- read.delim(path, check.names = FALSE)
  df <- fun(df, ...)
  write.table(df, path, sep = "\t", row.names = FALSE, quote = FALSE)
  message("Annotated: ", path)
}

message("=== Functional (pathway) outputs ===")
func <- file.path(output_root, "functional_analysis_output")
annotate_csv(file.path(func, "paired_wilcox_pain_gi_pathways.csv"), annotate_paired_wilcox_direction)
strip_direction_cols <- function(df) {
  keep <- !grepl("^(maaslin_|wilcox_|enriched_in|direction_)", names(df))
  df[, keep, drop = FALSE]
}

annotate_csv(file.path(func, "consensus_pain_gi_site_pathways.csv"), function(df) {
  df <- strip_direction_cols(df)
  maas <- read.delim(
    file.path(func, "maaslin2/pain_fecal_vs_cecal/all_results.tsv"),
    check.names = FALSE
  ) %>% annotate_maaslin_direction()
  wilcox <- read.csv(file.path(func, "paired_wilcox_pain_gi_pathways.csv"), stringsAsFactors = FALSE)
  add_consensus_direction(df, maas, wilcox)
})

for (sub in c("pain_fecal_vs_cecal", "pain_fecal_vs_cecal_paired_ssn", "fecal_metadata_associations")) {
  annotate_tsv(file.path(func, "maaslin2", sub, "all_results.tsv"), annotate_maaslin_direction)
  annotate_tsv(file.path(func, "maaslin2", sub, "significant_results.tsv"), annotate_maaslin_direction)
}
annotate_csv(file.path(func, "fecal_pathway_spearman_significant.csv"), annotate_spearman_direction)
annotate_csv(file.path(func, "top_pathways_fecal_vs_cecal.csv"), annotate_maaslin_direction)

message("=== Taxonomic outputs ===")
tax <- file.path(output_root, "analysis_output")
annotate_csv(file.path(tax, "paired_wilcox_pain_gi.csv"), annotate_paired_wilcox_direction)
for (sub in c("pain_fecal_vs_cecal", "pain_fecal_vs_cecal_paired_ssn",
              "fecal_metadata_associations", "fecal_pain_vs_no_pain")) {
  annotate_tsv(file.path(tax, "maaslin2", sub, "all_results.tsv"), annotate_maaslin_direction)
  annotate_tsv(file.path(tax, "maaslin2", sub, "significant_results.tsv"), annotate_maaslin_direction)
}
annotate_csv(file.path(tax, "fecal_spearman_significant.csv"), annotate_spearman_direction)
annotate_csv(file.path(tax, "fecal_maaslin_significant_by_metadata.csv"), annotate_maaslin_direction)

cons_path <- file.path(func, "consensus_pain_gi_site_pathways.csv")
cons_tax  <- file.path(tax, "consensus/consensus_pain_gi_site.csv")
if (file.exists(cons_tax)) {
  df <- strip_direction_cols(read.csv(cons_tax, stringsAsFactors = FALSE))
  maas <- read.delim(file.path(tax, "maaslin2/pain_fecal_vs_cecal/all_results.tsv"), check.names = FALSE) %>%
    annotate_maaslin_direction()
  wilcox <- read.csv(file.path(tax, "paired_wilcox_pain_gi.csv"), stringsAsFactors = FALSE)
  write.csv(add_consensus_direction(df, maas, wilcox), cons_tax, row.names = FALSE)
  message("Annotated: ", cons_tax)
}

message("Done. See R/direction_helpers.R for interpretation rules.")
