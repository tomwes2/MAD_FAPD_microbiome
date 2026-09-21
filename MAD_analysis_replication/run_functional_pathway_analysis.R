#!/usr/bin/env Rscript
# HUMAnN pathway analysis — mirrors microbiome_MAD_full_analysis cohorts & methods.
#
# Cohorts:
#   fvc = pain subjects with fecal + cecal samples (GI site comparison)
#   fps = all fecal samples (metadata associations; pop not in primary model)
#
# Direction (fecal vs cecal):
#   MaAsLin2 reference = "sample_type,Fecal" (stool is baseline).
#   Positive coef on Cecal rows => higher in CECUM; negative => higher in FECAL.
#   Paired Wilcoxon median_diff = median(cecal - fecal) within subject.
#
# After running, use annotate_direction_outputs.R if you need to refresh direction
# columns on saved tables without re-fitting MaAsLin2.
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
  library(dplyr)
  library(purrr)
  library(tibble)
  library(vegan)
  library(Maaslin2)
})
source(file.path(data_dir, "R", "direction_helpers.R"))

output_dir  <- file.path(output_root, "functional_analysis_output")
maaslin_dir <- file.path(output_dir, "maaslin2")
dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)

pathway_id <- function(feature) sub("\\|.*", "", feature)

load_sample_id_map <- function(csv_path = file.path(data_dir, "sample_ids.csv")) {
  if (!file.exists(csv_path)) return(NULL)
  m <- read.csv(csv_path, stringsAsFactors = FALSE)
  key <- sub("_.*", "", m$file_header)
  stats::setNames(m$sample_id, key)
}

humann_col_to_sample <- function(col, id_lookup = NULL) {
  x <- sub("_clean_Abundance$", "", col)
  key <- sub("_.*", "", x)
  if (!is.null(id_lookup) && key %in% names(id_lookup)) return(id_lookup[[key]])
  key
}

load_pathway_matrix <- function(tsv_path, id_lookup = NULL) {
  raw <- read.delim(tsv_path, check.names = FALSE, quote = "")
  feat <- raw[[1]]
  mat  <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- vapply(colnames(mat), humann_col_to_sample, character(1), id_lookup = id_lookup)
  keep_feat <- !grepl("^UNMAPPED$|^UNINTEGRATED", feat)
  mat <- mat[keep_feat, , drop = FALSE]
  feat <- feat[keep_feat]
  pid <- pathway_id(feat)
  split_idx <- split(seq_along(pid), pid)
  collapsed <- t(vapply(split_idx, function(ii) colSums(mat[ii, , drop = FALSE]), numeric(ncol(mat))))
  rownames(collapsed) <- names(split_idx)
  collapsed
}

`%||%` <- function(x, y) if (is.null(x)) y else x

prepare_maaslin_metadata <- function(meta, sample_ids, fixed_effects,
                                    extra_cols = c("ssn", "household_id")) {
  meta <- meta[sample_ids, , drop = FALSE]
  cols <- intersect(c(fixed_effects, extra_cols), colnames(meta))
  meta <- meta[, cols, drop = FALSE]
  rn <- rownames(meta)
  meta <- as.data.frame(lapply(meta, function(x) {
    if (is.data.frame(x)) as.vector(x[[1]]) else as.vector(x)
  }), stringsAsFactors = FALSE)
  rownames(meta) <- rn
  cat_vars <- c(
    "sample_type", "pop", "sex", "obese", "household_shared",
    "phq_difficulty", "gad_difficulty", "calprotectin_high", "ssn"
  )
  for (nm in intersect(cat_vars, names(meta))) meta[[nm]] <- factor(meta[[nm]])
  keep <- vapply(meta, function(x) length(unique(na.omit(x))) > 1, logical(1))
  meta[, keep, drop = FALSE]
}

run_maaslin_mat <- function(mat, meta, fixed_effects, output_subdir,
                            reference = NULL, random_effects = NULL) {
  out <- file.path(maaslin_dir, output_subdir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dat <- as.data.frame(t(mat), check.names = FALSE)
  meta <- prepare_maaslin_metadata(meta, rownames(dat), fixed_effects)
  meta <- meta[rownames(dat), , drop = FALSE]
  fx <- intersect(fixed_effects, colnames(meta))
  if (length(fx) == 0) stop("No valid fixed effects for ", output_subdir)
  args <- list(
    input_data = dat,
    input_metadata = meta,
    output = out,
    min_abundance = 0,
    min_prevalence = 0.1,
    max_significance = 0.25,
    normalization = "NONE",
    transform = "NONE",
    analysis_method = "LM",
    fixed_effects = fx,
    plot_scatter = FALSE
  )
  if (!is.null(reference)) args$reference <- reference
  re <- intersect(random_effects %||% character(0), colnames(meta))
  if (length(re) > 0) args$random_effects <- re
  message("MaAsLin2: ", output_subdir, " (", nrow(dat), " samples, ", ncol(dat), " pathways)")
  do.call(Maaslin2::Maaslin2, args)
}

run_permanova <- function(mat, meta, terms, label) {
  dist_mat <- vegdist(t(mat), method = "bray")
  bind_rows(lapply(terms, function(v) {
    f <- as.formula(paste("dist_mat ~", v))
    fit <- adonis2(f, data = meta, permutations = 9999)
    tibble(cohort = label, metadata = v, R2 = fit$R2[1], F = fit$F[1], p_value = fit$`Pr(>F)`[1])
  }))
}

paired_wilcox_screen <- function(mat, meta, group_var = "sample_type",
                                 p_thresh = 0.25, min_prev = 0.1) {
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  meta_df <- data.frame(
    ssn = meta$ssn,
    gi_site = as.character(meta[[group_var]]),
    row.names = rownames(meta),
    stringsAsFactors = FALSE
  )
  paired_ssn <- meta_df %>%
    group_by(ssn, gi_site) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(ssn) %>%
    summarise(levels = n(), .groups = "drop") %>%
    filter(levels == 2) %>%
    pull(ssn)
  lv <- levels(factor(meta[[group_var]]))
  purrr::map_dfr(taxa, function(pathway) {
    d <- meta_df %>% filter(ssn %in% paired_ssn)
    d$abund <- as.numeric(mat[pathway, rownames(d)])
    g1 <- d %>% filter(gi_site == lv[1]) %>% arrange(ssn)
    g2 <- d %>% filter(gi_site == lv[2]) %>% arrange(ssn)
    if (nrow(g1) != nrow(g2) || nrow(g1) == 0) return(NULL)
    tt <- wilcox.test(g1$abund, g2$abund, paired = TRUE)
    tibble(feature = pathway, comparison = paste(lv, collapse = "_vs_"),
           pval = tt$p.value, median_diff = median(g2$abund - g1$abund))
  }) %>%
    mutate(qval = p.adjust(pval, "BH")) %>%
    filter(qval < p_thresh) %>%
    arrange(qval)
}

lefse_screen <- function(mat, meta, group_var, p_thresh = 0.25, min_prev = 0.1) {
  grp <- factor(meta[[group_var]])
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  purrr::map_dfr(taxa, function(pathway) {
    x <- as.numeric(mat[pathway, rownames(meta)])
    kw <- kruskal.test(x ~ grp)
    tibble(feature = pathway, kruskal_p = kw$p.value)
  }) %>%
    mutate(comparison = group_var, qval = p.adjust(kruskal_p, "BH")) %>%
    filter(qval < p_thresh) %>%
    arrange(qval)
}

spearman_screen <- function(mat, meta, meta_var, p_thresh = 0.25, min_prev = 0.1) {
  x <- as.numeric(meta[[meta_var]])
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  purrr::map_dfr(taxa, function(pathway) {
    y <- as.numeric(mat[pathway, rownames(meta)])
    ct <- cor.test(y, x, method = "spearman", exact = FALSE)
    tibble(feature = pathway, rho = unname(ct$estimate), pval = ct$p.value)
  }) %>%
    mutate(metadata = meta_var, qval = p.adjust(pval, "BH")) %>%
    filter(qval < p_thresh) %>%
    arrange(qval)
}

# normalize_feature() lives in R/direction_helpers.R (sourced above)

build_consensus <- function(cohort_label, comparison, method_hits, min_methods = 2) {
  base <- bind_rows(lapply(names(method_hits), function(m) {
    method_hits[[m]] %>%
      transmute(feature, feature_norm = normalize_feature(feature), method = m, comparison)
  }))
  if (nrow(base) == 0) return(tibble())
  display <- base %>%
    group_by(feature_norm) %>%
    summarise(feature = feature[which.max(nchar(feature))][1], .groups = "drop")
  base %>%
    group_by(feature_norm, comparison) %>%
    summarise(n_methods = n_distinct(method), methods = paste(sort(unique(method)), collapse = "+"), .groups = "drop") %>%
    filter(n_methods >= min_methods) %>%
    left_join(display, by = "feature_norm") %>%
    select(feature, feature_norm, comparison, n_methods, methods) %>%
    mutate(cohort = cohort_label) %>%
    arrange(desc(n_methods), feature)
}

read_maaslin_hits <- function(out_subdir, metadata_terms, q_thresh = 0.25) {
  f <- file.path(maaslin_dir, out_subdir, "significant_results.tsv")
  if (!file.exists(f)) f <- file.path(maaslin_dir, out_subdir, "all_results.tsv")
  if (!file.exists(f)) return(tibble())
  read.delim(f, check.names = FALSE) %>%
    as_tibble() %>%
    filter(metadata %in% metadata_terms, qval < q_thresh) %>%
    transmute(feature, metadata, qval, coef, comparison = metadata)
}

# --- Load ---
message("Loading HUMAnN pathways...")
id_lookup <- load_sample_id_map()
if (!is.null(id_lookup)) message("Using sample_ids.csv (", length(id_lookup), " mappings)")
path_mat <- load_pathway_matrix(
  file.path(data_dir, "combined_pathabundance_relab.tsv"),
  id_lookup = id_lookup
)
mad_samdat_tab <- read.csv(file.path(data_dir, "subject_data.csv"), check.names = FALSE)
rownames(mad_samdat_tab) <- mad_samdat_tab$sampleID
common <- intersect(colnames(path_mat), rownames(mad_samdat_tab))
path_mat <- path_mat[, common, drop = FALSE]
mad_samdat_tab <- mad_samdat_tab[common, , drop = FALSE]
message("Pathways: ", nrow(path_mat), " | Samples: ", ncol(path_mat))

mad_samdat_tab$obese <- factor(ifelse(mad_samdat_tab$bmi > 30, "obese", "non_obese"), levels = c("non_obese", "obese"))
mad_samdat_tab$calprotectin_binary <- factor(
  ifelse(mad_samdat_tab$calprotectin_high %in% c(1, "1", TRUE), "high", "normal"),
  levels = c("normal", "high")
)
mad_samdat_tab$pop <- factor(mad_samdat_tab$pop, levels = c("no_pain", "pain"))
# Factor order matters for Wilcoxon: level 1 = Fecal, level 2 = Cecal (see median_diff sign).
mad_samdat_tab$sample_type <- factor(mad_samdat_tab$sample_type, levels = c("Fecal", "Cecal"))
mad_samdat_tab$sex <- factor(mad_samdat_tab$sex)
mad_samdat_tab$phq_difficulty <- factor(
  mad_samdat_tab$phq_difficulty, levels = 0:3,
  labels = c("PHQ: Not difficult", "PHQ: Somewhat difficult", "PHQ: Very difficult", "PHQ: Extremely difficult")
)
mad_samdat_tab$gad_difficulty <- factor(
  mad_samdat_tab$gad_difficulty, levels = 0:3,
  labels = c("GAD: Not difficult", "GAD: Somewhat difficult", "GAD: Very difficult", "GAD: Extremely difficult")
)
mad_samdat_tab$household_id <- ifelse(
  mad_samdat_tab$household == 0, paste0("unique_", mad_samdat_tab$ssn), paste0("hh_", mad_samdat_tab$household)
)
mad_samdat_tab$household_shared <- factor(
  ifelse(mad_samdat_tab$household == 0, "no_shared", "shared"), levels = c("no_shared", "shared")
)

fecal_predictors <- c(
  "calprotectin", "calprotectin_high", "age", "sex", "bmi",
  "phq_total", "phq_difficulty", "gad_total", "gad_difficulty", "household_shared"
)

fvc_samples <- rownames(mad_samdat_tab)[mad_samdat_tab$pop == "pain"]
fps_samples <- rownames(mad_samdat_tab)[mad_samdat_tab$sample_type == "Fecal"]
fvc_mat  <- path_mat[, fvc_samples, drop = FALSE]
fvc_meta <- mad_samdat_tab[fvc_samples, , drop = FALSE]
fps_mat  <- path_mat[, fps_samples, drop = FALSE]
fps_meta <- mad_samdat_tab[fps_samples, , drop = FALSE]

write.csv(
  tibble(cohort = c("Pain fecal+cecal", "Fecal only"), n_samples = c(ncol(fvc_mat), ncol(fps_mat))),
  file.path(output_dir, "cohort_summary.csv"), row.names = FALSE
)

# --- PERMANOVA ---
message("PERMANOVA...")
perm_all <- bind_rows(
  run_permanova(fvc_mat, fvc_meta, c("sample_type", "calprotectin_high", "phq_total", "gad_total"), "Pain fecal vs cecal"),
  run_permanova(fps_mat, fps_meta, fecal_predictors, "Fecal metadata")
)
write.csv(perm_all, file.path(output_dir, "pathway_permanova_results.csv"), row.names = FALSE)

# --- Pain cohort MaAsLin ---
maaslin_fvc <- run_maaslin_mat(
  fvc_mat, fvc_meta,
  fixed_effects = c(
    "sample_type", "calprotectin", "calprotectin_high",
    "age", "sex", "bmi", "obese",
    "phq_total", "phq_difficulty", "gad_total", "gad_difficulty", "household_shared"
  ),
  output_subdir = "pain_fecal_vs_cecal",
  reference = "sample_type,Fecal",
  random_effects = "household_id"
)

out_paired <- file.path(maaslin_dir, "pain_fecal_vs_cecal_paired_ssn")
dir.create(out_paired, recursive = TRUE, showWarnings = FALSE)
dat <- t(fvc_mat)
meta_p <- prepare_maaslin_metadata(
  fvc_meta, rownames(dat),
  c("sample_type", "calprotectin", "calprotectin_high", "age", "sex", "bmi", "obese", "phq_total", "gad_total")
)
message("MaAsLin2: pain_fecal_vs_cecal_paired_ssn")
Maaslin2(
  input_data = dat,
  input_metadata = meta_p,
  output = out_paired,
  min_abundance = 0, min_prevalence = 0.1, max_significance = 0.25,
  normalization = "NONE", transform = "NONE",
  analysis_method = "LM",
  fixed_effects = c("sample_type"),
  random_effects = c("ssn"),
  reference = "sample_type,Fecal",
  plot_scatter = FALSE
)

paired_fvc <- paired_wilcox_screen(fvc_mat, fvc_meta, "sample_type") %>%
  annotate_paired_wilcox_direction()
lefse_fvc  <- lefse_screen(fvc_mat, fvc_meta, "sample_type")
write.csv(paired_fvc, file.path(output_dir, "paired_wilcox_pain_gi_pathways.csv"), row.names = FALSE)
write.csv(lefse_fvc,  file.path(output_dir, "lefse_style_pain_gi_pathways.csv"), row.names = FALSE)

# --- Fecal metadata MaAsLin ---
maaslin_fecal_meta <- run_maaslin_mat(
  fps_mat, fps_meta,
  fixed_effects = fecal_predictors,
  output_subdir = "fecal_metadata_associations",
  random_effects = "household_id"
)
maaslin_fecal_sig <- maaslin_fecal_meta$results %>%
  as_tibble() %>%
  filter(qval < 0.25, metadata %in% fecal_predictors) %>%
  annotate_maaslin_direction() %>%
  arrange(qval)
write.csv(maaslin_fecal_sig, file.path(output_dir, "fecal_pathway_maaslin_significant_by_metadata.csv"), row.names = FALSE)

fecal_continuous <- c("calprotectin", "age", "bmi", "phq_total", "gad_total")
fecal_categorical <- c("calprotectin_high", "sex", "phq_difficulty", "gad_difficulty", "household_shared")
spearman_hits <- bind_rows(lapply(fecal_continuous, function(v) spearman_screen(fps_mat, fps_meta, v))) %>%
  annotate_spearman_direction()
lefse_fecal_meta <- bind_rows(lapply(fecal_categorical, function(v) {
  lefse_screen(fps_mat, fps_meta, v) %>% mutate(metadata = v)
}))
write.csv(spearman_hits, file.path(output_dir, "fecal_pathway_spearman_significant.csv"), row.names = FALSE)
write.csv(lefse_fecal_meta, file.path(output_dir, "fecal_pathway_lefse_style_significant.csv"), row.names = FALSE)

# --- Consensus (with direction from MaAsLin + paired Wilcoxon) ---
q_thresh <- 0.25
maaslin_fvc_all <- read.delim(
  file.path(maaslin_dir, "pain_fecal_vs_cecal/all_results.tsv"),
  check.names = FALSE
) %>% annotate_maaslin_direction()

consensus_fvc <- build_consensus(
  "pain_pathways", "sample_type",
  list(
    maaslin2 = read_maaslin_hits("pain_fecal_vs_cecal", "sample_type", q_thresh) %>% distinct(feature, .keep_all = TRUE),
    maaslin2_paired = read_maaslin_hits("pain_fecal_vs_cecal_paired_ssn", "sample_type", q_thresh) %>% distinct(feature, .keep_all = TRUE),
    paired_wilcox = paired_fvc %>% distinct(feature, .keep_all = TRUE),
    lefse = lefse_fvc %>% distinct(feature, .keep_all = TRUE)
  )
) %>%
  add_consensus_direction(maaslin_fvc_all, paired_fvc)
write.csv(consensus_fvc, file.path(output_dir, "consensus_pain_gi_site_pathways.csv"), row.names = FALSE)

# Annotate full MaAsLin tables on disk for downstream use
write_maaslin_with_direction(file.path(maaslin_dir, "pain_fecal_vs_cecal/all_results.tsv"))
write_maaslin_with_direction(file.path(maaslin_dir, "pain_fecal_vs_cecal/significant_results.tsv"))
write_maaslin_with_direction(file.path(maaslin_dir, "pain_fecal_vs_cecal_paired_ssn/all_results.tsv"))
write_maaslin_with_direction(file.path(maaslin_dir, "fecal_metadata_associations/all_results.tsv"))

maaslin_fvc_hits <- read_maaslin_hits("pain_fecal_vs_cecal", "sample_type", q_thresh) %>%
  left_join(
    maaslin_fvc_all %>% filter(metadata == "sample_type", value == "Cecal") %>%
      select(feature, enriched_in, direction_label, coef_interpretation),
    by = "feature"
  )
write.csv(maaslin_fvc_hits %>% arrange(qval) %>% slice_head(n = 30),
          file.path(output_dir, "top_pathways_fecal_vs_cecal.csv"), row.names = FALSE)
write.csv(maaslin_fecal_sig %>% group_by(metadata) %>% slice_min(qval, n = 10) %>% ungroup(),
          file.path(output_dir, "top_pathways_fecal_metadata.csv"), row.names = FALSE)

# Metadata-specific MaAsLin exports
for (term in unique(maaslin_fvc$results$metadata)) {
  sub <- maaslin_fvc$results %>% as_tibble() %>% filter(metadata == term, qval < 0.25) %>% arrange(qval)
  if (nrow(sub) > 0) {
    write.csv(sub, file.path(output_dir, paste0("pain_cohort_pathways_by_", term, ".csv")), row.names = FALSE)
  }
}

message("Done. Outputs in ", output_dir)
