#!/usr/bin/env Rscript
# Build manuscript tables, figures, supplementary analyses, and Results documents.
#
# Prerequisite: functional + taxonomic analysis outputs (see run_functional_pathway_analysis.R).
# Direction columns: sourced from R/direction_helpers.R (MaAsLin Fecal reference; Wilcoxon
# median_diff = cecal - fecal). Run annotate_direction_outputs.R to refresh direction on CSVs.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
})
has_rf <- requireNamespace("randomForest", quietly = TRUE)
if (has_rf) suppressPackageStartupMessages(library(randomForest))
args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg)) else normalizePath(".")
source(file.path(script_dir, "R", "paths.R"))
paths <- mad_init_paths()
data_dir <- paths$data_dir
output_root <- paths$output_root
setwd(data_dir)
source(file.path(data_dir, "R", "manuscript_utils.R"))
source(file.path(data_dir, "R", "direction_helpers.R"))

analysis_dir <- file.path(output_root, "analysis_output")
func_dir <- file.path(output_root, "functional_analysis_output")
out_root <- file.path(output_root, "manuscript_outputs")
p1_dir   <- file.path(out_root, "paper1_fecal_vs_cecal")
p2_dir   <- file.path(out_root, "paper2_fecal_metadata")
shared   <- file.path(out_root, "shared")
for (d in c(p1_dir, p2_dir, shared, file.path(p1_dir, "figures"), file.path(p1_dir, "tables"),
             file.path(p2_dir, "figures"), file.path(p2_dir, "tables"), file.path(shared, "figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("=== Loading data ===")
meta <- read.csv(file.path(data_dir, "subject_data.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sampleID
meta <- prepare_meta(meta)

meta_with_id <- function(m) {
  out <- data.frame(m, sample_id_row = rownames(m), stringsAsFactors = FALSE, check.names = FALSE)
  tibble::as_tibble(out)
}

id_lookup <- load_sample_id_map(file.path(data_dir, "sample_ids.csv"))
path_mat <- load_pathway_matrix(file.path(data_dir, "combined_pathabundance_relab.tsv"), id_lookup)
common <- intersect(colnames(path_mat), rownames(meta))
path_mat <- path_mat[, common, drop = FALSE]
meta <- meta[common, , drop = FALSE]

fvc_ids <- rownames(meta)[meta$pop == "pain"]
fps_ids <- rownames(meta)[meta$sample_type == "Fecal"]
fvc_mat <- path_mat[, fvc_ids, drop = FALSE]
fvc_meta <- meta[fvc_ids, , drop = FALSE]
fps_mat <- path_mat[, fps_ids, drop = FALSE]
fps_meta <- meta[fps_ids, , drop = FALSE]

# Sample flow
humann_cols <- strsplit(readLines(file.path(data_dir, "combined_pathabundance_relab.tsv"), n = 1), "\t")[[1]][-1]
flow <- tibble(
  step = c(
    "Metadata rows (subject_data.csv)",
    "HUMAnN columns (combined_pathabundance_relab.tsv)",
    "Mapped to metadata (sample_ids.csv)",
    "Pain cohort (fecal + cecal)",
    "Fecal-only cohort",
    "Paired pain subjects (both GI sites)",
    "Missing HUMAnN (MAD059C only)"
  ),
  n = c(
    nrow(read.csv(file.path(data_dir, "subject_data.csv"))),
    length(humann_cols),
    length(common),
    length(fvc_ids),
    length(fps_ids),
    sum(table(fvc_meta$ssn, fvc_meta$sample_type) == 1) * 0 + sum(apply(table(fvc_meta$ssn, fvc_meta$sample_type) > 0, 1, sum) == 2),
    sum(!meta$sampleID %in% colnames(path_mat) & meta$sampleID == "MAD059C")
  )
)
write.csv(flow, file.path(shared, "sample_flow.csv"), row.names = FALSE)

# Taxonomic OTU table (local kit only — do not fall back to /DATA01/MAD)
otu_paths <- c(file.path(data_dir, "otu_table.csv"))
otu_path <- otu_paths[file.exists(otu_paths)][1]
has_otu <- !is.na(otu_path)
tax_alpha <- NULL
if (has_otu) {
  message("Found OTU table: ", otu_path)
  otu <- read.csv(otu_path, row.names = 1, check.names = FALSE)
  otu <- otu[, intersect(colnames(otu), rownames(meta)), drop = FALSE]
  tax_alpha <- alpha_metrics(otu)
  tax_alpha <- tax_alpha %>% left_join(meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
  write.csv(tax_alpha, file.path(shared, "taxonomy_alpha_diversity.csv"), row.names = FALSE)
}

# --- Functional alpha diversity ---
func_alpha_fvc <- alpha_metrics(fvc_mat) %>%
  left_join(fvc_meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
func_alpha_fps <- alpha_metrics(fps_mat) %>%
  left_join(fps_meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
write.csv(func_alpha_fvc, file.path(p1_dir, "tables", "pathway_alpha_diversity_pain.csv"), row.names = FALSE)
write.csv(func_alpha_fps, file.path(p2_dir, "tables", "pathway_alpha_diversity_fecal.csv"), row.names = FALSE)

alpha_tests_p1 <- bind_rows(
  lapply(c("richness", "shannon"), function(m) {
    g1 <- func_alpha_fvc %>% filter(sample_type == "Fecal") %>% pull(.data[[m]])
    g2 <- func_alpha_fvc %>% filter(sample_type == "Cecal") %>% pull(.data[[m]])
    wt <- wilcox.test(g1, g2, paired = FALSE)
    paired_ssn <- func_alpha_fvc %>%
      select(ssn, sample_type, val = .data[[m]]) %>%
      pivot_wider(names_from = sample_type, values_from = val) %>%
      filter(!is.na(Fecal) & !is.na(Cecal))
    wtp <- if (nrow(paired_ssn) >= 3) {
      wilcox.test(paired_ssn$Fecal, paired_ssn$Cecal, paired = TRUE)
    } else {
      list(p.value = NA_real_)
    }
    tibble(metric = m, unpaired_p = wt$p.value, paired_p = wtp$p.value,
           median_fecal = median(g1), median_cecal = median(g2))
  })
)
write.csv(alpha_tests_p1, file.path(p1_dir, "tables", "pathway_alpha_tests.csv"), row.names = FALSE)

alpha_tests_p2 <- bind_rows(
  lapply(c("richness", "shannon"), function(m) {
    g1 <- func_alpha_fps %>% filter(pop == "no_pain") %>% pull(.data[[m]])
    g2 <- func_alpha_fps %>% filter(pop == "pain") %>% pull(.data[[m]])
    wt <- wilcox.test(g1, g2)
    tibble(metric = m, p_value = wt$p.value,
           median_no_pain = median(g1), median_pain = median(g2))
  })
)
write.csv(alpha_tests_p2, file.path(p2_dir, "tables", "pathway_alpha_pain_vs_control.csv"), row.names = FALSE)

# Alpha plots
p_alpha_p1 <- ggplot(func_alpha_fvc, aes(sample_type, shannon, fill = sample_type)) +
  geom_violin(trim = FALSE, alpha = 0.7) + geom_boxplot(width = 0.12, outlier.shape = NA) +
  labs(title = "Pathway Shannon diversity — pain cohort", y = "Shannon", x = NULL) +
  theme_bw() + theme(legend.position = "none")
save_plot(p_alpha_p1, file.path(p1_dir, "figures", "pathway_shannon_by_site"))

# --- PERMANOVA (pathway) extended ---
perm_vars_p1 <- c("sample_type", "calprotectin", "calprotectin_high", "phq_total", "gad_total", "bmi", "age")
perm_p1 <- bind_rows(lapply(perm_vars_p1, function(v) permanova_one(fvc_mat, fvc_meta, v, "pain_pathways")))
perm_vars_p2 <- c("calprotectin", "calprotectin_high", "age", "sex", "bmi", "phq_total",
                  "phq_difficulty", "gad_total", "gad_difficulty", "household_shared", "pop")
perm_p2 <- bind_rows(lapply(perm_vars_p2, function(v) permanova_one(fps_mat, fps_meta, v, "fecal_pathways")))
write.csv(bind_rows(perm_p1, perm_p2), file.path(shared, "pathway_permanova_extended.csv"), row.names = FALSE)

# Sensitivity: no recent antibiotics (fecal)
fps_abx0 <- fps_meta[fps_meta$antibiotics_recent == "no", , drop = FALSE]
perm_abx <- bind_rows(
  lapply(c("calprotectin", "bmi", "household_shared", "pop"), function(v) {
    permanova_one(fps_mat[, rownames(fps_abx0), drop = FALSE], fps_abx0, v, "fecal_no_antibiotics")
  })
)
write.csv(perm_abx, file.path(p2_dir, "tables", "permanova_antibiotics_excluded.csv"), row.names = FALSE)

# Pain-only fecal subset
fps_pain <- fps_meta[fps_meta$pop == "pain", , drop = FALSE]
perm_pain_only <- bind_rows(
  lapply(c("calprotectin", "bmi", "phq_total", "gad_total"), function(v) {
    permanova_one(fps_mat[, rownames(fps_pain), drop = FALSE], fps_pain, v, "fecal_pain_only")
  })
)
write.csv(perm_pain_only, file.path(p2_dir, "tables", "permanova_pain_only_fecal.csv"), row.names = FALSE)

# --- PCoA plots (pathways) ---
make_pcoa <- function(mat, meta, color_var, title, out_base) {
  common <- intersect(colnames(mat), rownames(meta))
  d <- vegdist(t(mat[, common, drop = FALSE]), method = "bray")
  pco <- cmdscale(d, k = 2, eig = TRUE)
  pct <- round(100 * pco$eig / sum(pco$eig), 1)
  df <- tibble(PC1 = pco$points[, 1], PC2 = pco$points[, 2], sampleID = common) %>%
    left_join(meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
  p <- ggplot(df, aes(PC1, PC2, color = .data[[color_var]])) +
    geom_point(size = 2.5, alpha = 0.85) +
    stat_ellipse(linewidth = 0.4) +
    labs(
      title = title,
      x = paste0("PC1 (", pct[1], "%)"),
      y = paste0("PC2 (", pct[2], "%)"),
      color = color_var
    ) + theme_bw()
  save_plot(p, out_base, width = 8, height = 6)
  df
}
make_pcoa(fvc_mat, fvc_meta, "sample_type", "Pathway PCoA — pain (fecal vs cecal)",
          file.path(p1_dir, "figures", "pcoa_pathways_by_site"))
make_pcoa(fps_mat, fps_meta, "pop", "Pathway PCoA — fecal (pain vs no pain)",
          file.path(p2_dir, "figures", "pcoa_pathways_by_pop"))
make_pcoa(fps_mat, fps_meta, "household_shared", "Pathway PCoA — household sharing",
          file.path(p2_dir, "figures", "pcoa_pathways_by_household"))

dir.create(file.path(shared, "tables"), showWarnings = FALSE, recursive = TRUE)

# --- Strict MaAsLin tables (q %%< 0.05) ---
maaslin_dirs <- list(
  pain_site = file.path(func_dir, "maaslin2/pain_fecal_vs_cecal/all_results.tsv"),
  pain_paired = file.path(func_dir, "maaslin2/pain_fecal_vs_cecal_paired_ssn/all_results.tsv"),
  fecal_meta = file.path(func_dir, "maaslin2/fecal_metadata_associations/all_results.tsv"),
  tax_pain = file.path(analysis_dir, "maaslin2/pain_fecal_vs_cecal/all_results.tsv"),
  tax_fecal = file.path(analysis_dir, "maaslin2/fecal_metadata_associations/all_results.tsv")
)
for (nm in names(maaslin_dirs)) {
  df <- read_maaslin_all(maaslin_dirs[[nm]]) %>% annotate_maaslin_direction()
  if (!nrow(df)) next
  write.csv(strict_maaslin(df, 0.05), file.path(shared, "tables", paste0("strict_q05_", nm, ".csv")), row.names = FALSE)
  write.csv(strict_maaslin(df, 0.25), file.path(shared, "tables", paste0("exploratory_q25_", nm, ".csv")), row.names = FALSE)
}

# --- Effect sizes and direction (fecal vs cecal) ---
cons_path <- read.csv(file.path(func_dir, "consensus_pain_gi_site_pathways.csv"),
                      stringsAsFactors = FALSE)
cons_tax <- read.csv(file.path(analysis_dir, "consensus/consensus_pain_gi_site.csv"),
                     stringsAsFactors = FALSE)
paired_pw <- read.csv(file.path(func_dir, "paired_wilcox_pain_gi_pathways.csv"),
                      stringsAsFactors = FALSE)
paired_tax <- read.csv(file.path(analysis_dir, "paired_wilcox_pain_gi.csv"),
                       stringsAsFactors = FALSE)
if (!"direction_label" %in% names(paired_pw)) paired_pw <- annotate_paired_wilcox_direction(paired_pw)
if (!"direction_label" %in% names(paired_tax)) paired_tax <- annotate_paired_wilcox_direction(paired_tax)
if (!"direction_label" %in% names(cons_path)) {
  maas_all_pw <- read_maaslin_all(maaslin_dirs$pain_site) %>% annotate_maaslin_direction()
  cons_path <- add_consensus_direction(cons_path, maas_all_pw, paired_pw)
}
if (!"direction_label" %in% names(cons_tax)) {
  maas_all_tax <- read_maaslin_all(maaslin_dirs$tax_pain) %>% annotate_maaslin_direction()
  cons_tax <- add_consensus_direction(cons_tax, maas_all_tax, paired_tax)
}
maas_pw <- read_maaslin_all(maaslin_dirs$pain_site) %>%
  filter(metadata == "sample_type", value == "Cecal") %>%
  annotate_maaslin_direction()
maas_tax <- read_maaslin_all(maaslin_dirs$tax_pain) %>%
  filter(metadata == "sample_type", value == "Cecal") %>%
  annotate_maaslin_direction()

# Combined consensus list with direction for manuscript / supplementary tables
write.csv(
  bind_rows(
    cons_tax %>% transmute(layer = "taxon", feature, n_methods, methods, cohort,
                           enriched_in, direction_label, direction_agreement = NA_character_),
    cons_path %>% transmute(layer = "pathway", feature, n_methods, methods, cohort,
                            enriched_in, direction_label, direction_agreement)
  ),
  file.path(shared, "tables", "pain_site_consensus_features.csv"),
  row.names = FALSE
)

top_pw <- head(cons_path$feature, 12)
effect_pw <- lapply(top_pw, function(f) {
  fmat <- match_maaslin_feature(f, rownames(fvc_mat))
  if (is.na(fmat)) return(NULL)
  mrow <- maas_pw %>% filter(feature == f | make.names(feature) == make.names(f)) %>% slice(1)
  prow <- paired_pw %>% filter(feature == f | make.names(feature) == make.names(f)) %>% slice(1)
  tibble(
    feature = fmat,
    maaslin_feature = f,
    maaslin_coef = if (nrow(mrow)) mrow$coef else NA_real_,
    maaslin_q = if (nrow(mrow)) mrow$qval else NA_real_,
    enriched_in = if (nrow(mrow)) mrow$enriched_in else NA_character_,
    direction_label = if (nrow(mrow)) mrow$direction_label else NA_character_,
    median_diff_cecal_minus_fecal = if (nrow(prow)) prow$median_diff else NA_real_,
    wilcox_enriched_in = if (nrow(prow)) prow$enriched_in %||% site_enriched_label(prow$median_diff) else NA_character_,
    wilcox_direction = if (nrow(prow)) prow$direction_label %||% site_direction_label(prow$median_diff) else NA_character_,
    wilcox_q = if (nrow(prow)) prow$qval else NA_real_,
    mean_rel_abund_fecal = mean(fvc_mat[fmat, fvc_meta$sample_type == "Fecal"]),
    mean_rel_abund_cecal = mean(fvc_mat[fmat, fvc_meta$sample_type == "Cecal"]),
    prevalence_fecal = mean(fvc_mat[fmat, fvc_meta$sample_type == "Fecal"] > 0),
    prevalence_cecal = mean(fvc_mat[fmat, fvc_meta$sample_type == "Cecal"] > 0)
  )
}) %>% bind_rows()
write.csv(effect_pw, file.path(p1_dir, "tables", "top_pathway_effect_sizes.csv"), row.names = FALSE)

cons_tax_strong <- read.csv(file.path(analysis_dir, "consensus/consensus_pain_gi_site_strong.csv"),
                            stringsAsFactors = FALSE)
effect_tax <- maas_tax %>%
  filter(feature %in% cons_tax_strong$feature) %>%
  transmute(
    feature, taxon_label = taxon_label(feature),
    maaslin_coef = coef, maaslin_q = qval,
    enriched_in, direction_label
  ) %>%
  left_join(
    paired_tax %>%
      transmute(
        feature,
        median_diff_cecal_minus_fecal = median_diff,
        wilcox_enriched_in = enriched_in %||% site_enriched_label(median_diff),
        wilcox_direction = direction_label %||% site_direction_label(median_diff),
        wilcox_q = qval
      ),
    by = "feature"
  )
write.csv(effect_tax, file.path(p1_dir, "tables", "top_taxon_effect_sizes.csv"), row.names = FALSE)

# --- Forest plots (coef > 0 => higher in cecum; fecal is MaAsLin reference) ---
forest_plot <- function(df, title, out_base, n = 20) {
  if (!nrow(df)) return(invisible(NULL))
  d <- df %>% arrange(qval) %>% slice_head(n = n) %>%
    mutate(
      dir_tag = ifelse(!is.na(enriched_in) & enriched_in %in% c("Cecal", "Fecal"),
                       paste0(" [", enriched_in, ">other]"), ""),
      label = paste0(substr(feature, 1, 52), dir_tag)
    )
  p <- ggplot(d, aes(x = coef, y = reorder(label, coef), color = qval < 0.05)) +
    geom_point(size = 2.5) + geom_vline(xintercept = 0, linetype = 2) +
    scale_color_manual(
      values = c("TRUE" = "#B2182B", "FALSE" = "grey50"),
      labels = c("TRUE" = "q < 0.05", "FALSE" = "q >= 0.05")
    ) +
    labs(
      title = title,
      subtitle = "Positive coef = higher in cecum; negative = higher in fecal (stool); Fecal is reference",
      x = "MaAsLin2 coefficient (Cecal vs Fecal)",
      y = NULL, color = NULL
    ) +
    theme_bw()
  save_plot(p, out_base, width = 11, height = max(5, n * 0.28))
}

forest_plot(maas_pw, "Pathways — fecal vs cecal (pain)", file.path(p1_dir, "figures", "forest_pathways_site"))
forest_plot(maas_tax, "Taxa — fecal vs cecal (pain)", file.path(p1_dir, "figures", "forest_taxa_site"))

maas_fecal <- read_maaslin_all(maaslin_dirs$tax_fecal) %>%
  filter(metadata %in% c("calprotectin", "gad_difficulty", "bmi", "household_shared")) %>%
  annotate_maaslin_direction()
forest_plot(maas_fecal %>% filter(metadata == "calprotectin"),
            "Taxa — calprotectin (fecal)", file.path(p2_dir, "figures", "forest_taxa_calprotectin"))
forest_plot(maas_fecal %>% filter(metadata == "gad_difficulty"),
            "Taxa — GAD difficulty (fecal)", file.path(p2_dir, "figures", "forest_taxa_gad_difficulty"))

# --- Boxplots top pathways ---
plot_pathway_boxes <- function(features, mat, meta, group_var, out_base) {
  long <- lapply(features, function(f) {
    fmat <- match_maaslin_feature(f, rownames(mat))
    if (is.na(fmat)) return(NULL)
    tibble(sampleID = colnames(mat), abundance = as.numeric(mat[fmat, ]), feature = fmat)
  }) %>% bind_rows() %>%
    left_join(meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
  long$feature <- substr(long$feature, 1, 50)
  p <- ggplot(long, aes(.data[[group_var]], abundance, fill = .data[[group_var]])) +
    geom_boxplot(outlier.size = 0.6) + facet_wrap(~feature, scales = "free_y", ncol = 3) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Top pathway abundances", x = NULL, y = "Relative abundance")
  save_plot(p, out_base, width = 11, height = 9)
}
plot_pathway_boxes(top_pw[1:min(9, length(top_pw))], fvc_mat, fvc_meta, "sample_type",
                 file.path(p1_dir, "figures", "boxplots_top_pathways_by_site"))

# Exploratory pathway-metadata scatter (fecal; Spearman screen)
sp <- read.csv(file.path(func_dir, "fecal_pathway_spearman_significant.csv"))
scatter_pathways <- function(features, xvar, title, out_base) {
  if (!length(features)) return(invisible(NULL))
  long <- lapply(features, function(f) {
    frow <- match_maaslin_feature(f, rownames(fps_mat))
    if (is.na(frow)) return(NULL)
    tibble(sampleID = colnames(fps_mat), abundance = as.numeric(fps_mat[frow, ]),
           feature = substr(frow, 1, 45))
  }) %>% bind_rows() %>%
    left_join(fps_meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row"))
  p <- ggplot(long, aes(.data[[xvar]], abundance)) + geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "#2166AC") +
    facet_wrap(~feature, scales = "free_y") +
    labs(title = title, x = xvar, y = "Relative abundance") + theme_bw()
  save_plot(p, out_base, width = 10, height = 4)
}
top_cal_pw <- sp %>% filter(metadata == "calprotectin") %>% arrange(qval) %>% slice_head(n = 3) %>% pull(feature)
scatter_pathways(top_cal_pw, "calprotectin", "Pathways vs fecal calprotectin (exploratory)",
                 file.path(p2_dir, "figures", "scatter_pathways_calprotectin"))
top_bmi_pw <- sp %>% filter(metadata == "bmi") %>% arrange(qval) %>% slice_head(n = 3) %>% pull(feature)
scatter_pathways(top_bmi_pw, "bmi", "Pathways vs BMI (exploratory Spearman)",
                 file.path(p2_dir, "figures", "scatter_pathways_bmi"))

# Taxon calprotectin boxplots if OTU available
if (has_otu) {
  tax_cal <- read.csv(file.path(analysis_dir, "fecal_maaslin_significant_by_metadata.csv")) %>%
    filter(metadata == "calprotectin") %>% arrange(qval) %>% slice_head(n = 3)
  if (nrow(tax_cal)) {
    long <- lapply(tax_cal$feature, function(f) {
      f_otu <- match_maaslin_feature(f, rownames(otu))
      if (is.na(f_otu)) return(NULL)
      tibble(sampleID = colnames(otu), abundance = as.numeric(otu[f_otu, ]), taxon = taxon_label(f_otu))
    }) %>% bind_rows()
    if (nrow(long) && "sampleID" %in% names(long)) {
      long <- long %>%
        left_join(fps_meta %>% meta_with_id(), by = c("sampleID" = "sample_id_row")) %>%
        dplyr::filter(!is.na(calprotectin_high))
      if (nrow(long)) {
        long$calprotectin_high <- factor(
          long$calprotectin_high,
          levels = c(0, 1),
          labels = c("normal", "high")
        )
        p <- ggplot(long, aes(calprotectin_high, abundance, fill = calprotectin_high)) +
          geom_boxplot() + facet_wrap(~taxon, scales = "free_y") +
          labs(title = "Taxa vs calprotectin elevation", x = "Calprotectin high", y = "Rel. abundance") +
          theme_bw()
        save_plot(p, file.path(p2_dir, "figures", "boxplots_taxa_calprotectin"), width = 9, height = 4)
      }
    }
  }
}

# --- Random forest (pathways) ---
run_rf <- function(mat, meta, response, label) {
  if (!has_rf) return(list(oob = NA_real_, importance = tibble(), label = label))
  common <- intersect(colnames(mat), rownames(meta))
  x <- t(mat[, common, drop = FALSE])
  colnames(x) <- make.names(colnames(x), unique = TRUE)
  y <- factor(meta[common, response])
  prev <- colSums(x > 0)
  keep <- names(sort(prev, decreasing = TRUE))[seq_len(min(80, ncol(x)))]
  df <- data.frame(y = y, x[, keep, drop = FALSE], check.names = FALSE)
  set.seed(42)
  rf <- randomForest(y ~ ., data = df, ntree = 1500, importance = TRUE)
  imp <- as.data.frame(importance(rf)) %>% tibble::rownames_to_column("feature")
  list(oob = rf$err.rate[nrow(rf$err.rate), "OOB"], importance = imp, label = label)
}
rf_p1 <- run_rf(fvc_mat, fvc_meta, "sample_type", "pathways_site")
rf_p2 <- run_rf(fps_mat, fps_meta, "pop", "pathways_pop")
rf_summary <- tibble(
  analysis = c(rf_p1$label, rf_p2$label),
  OOB_error = c(rf_p1$oob, rf_p2$oob),
  n_samples = c(ncol(fvc_mat), ncol(fps_mat))
)
write.csv(rf_summary, file.path(shared, "random_forest_pathway_summary.csv"), row.names = FALSE)
write.csv(rf_p1$importance, file.path(p1_dir, "tables", "rf_pathway_importance_site.csv"), row.names = FALSE)
rf_tax_dir <- file.path(shared, "random_forest_taxonomy")
dir.create(rf_tax_dir, showWarnings = FALSE)
for (rf_csv in list.files(file.path(analysis_dir, "random_forest"), pattern = "\\.csv$", full.names = TRUE)) {
  file.copy(rf_csv, file.path(rf_tax_dir, basename(rf_csv)), overwrite = TRUE)
}

# Taxonomic RF OOB from re-run if otu present
if (has_otu) {
  rf_tax_site <- run_rf(otu[, fvc_ids[fvc_ids %in% colnames(otu)], drop = FALSE], fvc_meta, "sample_type", "taxa_site")
  rf_tax_pop <- run_rf(otu[, fps_ids[fps_ids %in% colnames(otu)], drop = FALSE], fps_meta, "pop", "taxa_pop")
  write.csv(tibble(analysis = c("taxa_site", "taxa_pop"), OOB_error = c(rf_tax_site$oob, rf_tax_pop$oob)),
            file.path(shared, "random_forest_taxonomy_summary.csv"), row.names = FALSE)
}

# --- Spearman summaries (fecal) ---
if (file.exists(file.path(analysis_dir, "fecal_spearman_significant.csv"))) {
  file.copy(file.path(analysis_dir, "fecal_spearman_significant.csv"),
            file.path(p2_dir, "tables", "taxonomy_spearman_q25.csv"), overwrite = TRUE)
}

# --- Collect stats for Results text ---
perm_combined <- file.path(analysis_dir, "permanova_results.csv")
if (!file.exists(perm_combined)) {
  parts <- c(
    file.path(analysis_dir, "permanova_pain_fecal_vs_cecal.csv"),
    file.path(analysis_dir, "permanova_fecal.csv")
  )
  parts <- parts[file.exists(parts)]
  if (!length(parts)) stop("Missing taxonomic PERMANOVA tables in ", analysis_dir)
  write.csv(dplyr::bind_rows(lapply(parts, read.csv, stringsAsFactors = FALSE)),
            perm_combined, row.names = FALSE)
}
perm_tax <- read.csv(perm_combined)
perm_fecal_tax <- read.csv(file.path(analysis_dir, "fecal_permanova_by_metadata.csv"))
perm_tax_p <- perm_tax[[grep("^Pr", names(perm_tax), value = TRUE)[1]]]
n_cons_4 <- sum(read.csv(file.path(analysis_dir, "consensus/method_overlap_pain.csv"))$n == 4)

stats <- list(
  n_pain_paired = sum(apply(table(fvc_meta$ssn, fvc_meta$sample_type) > 0, 1, sum) == 2),
  n_fvc = length(fvc_ids),
  n_fps = length(fps_ids),
  n_fps_pain = sum(fps_meta$pop == "pain"),
  n_fps_control = sum(fps_meta$pop == "no_pain"),
  perm_tax_site_R2 = perm_tax$R2[perm_tax$cohort == "pain_fvc"][1],
  perm_tax_site_p = perm_tax_p[perm_tax$cohort == "pain_fvc"][1],
  perm_path_site_R2 = perm_p1$R2[perm_p1$metadata == "sample_type"],
  perm_path_site_p = perm_p1$p_value[perm_p1$metadata == "sample_type"],
  perm_fecal_pop_p = perm_tax_p[perm_tax$cohort == "fecal"][2],
  perm_fecal_hh_p = perm_fecal_tax$p_value[perm_fecal_tax$metadata == "household_shared"],
  perm_fecal_hh_R2 = perm_fecal_tax$R2[perm_fecal_tax$metadata == "household_shared"],
  perm_fecal_bmi_p = perm_fecal_tax$p_value[perm_fecal_tax$metadata == "bmi"],
  n_consensus_tax_4method = n_cons_4,
  n_consensus_path_2maaslin = nrow(cons_path),
  alpha_paired_shannon = alpha_tests_p1$paired_p[alpha_tests_p1$metric == "shannon"],
  rf_path_oob_site = rf_p1$oob,
  has_otu = has_otu
)

saveRDS(stats, file.path(shared, "results_stats.rds"))

# --- Generate Results markdown ---
fmt_p <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))

otu_note <- if (has_otu) {
  "Taxonomic alpha diversity statistics are in `shared/taxonomy_alpha_diversity.csv`."
} else {
  "**Note:** Place `otu_table.csv` in the replication kit directory to generate taxonomic alpha diversity and calprotectin taxon boxplots."
}

paper1_lines <- c(
  "# Results - Paper 1: Fecal versus cecal microbiome in adolescents with abdominal pain",
  "",
  "## Participants and profiling",
  "",
  paste0(
    "We analyzed **", stats$n_pain_paired, " adolescents with abdominal pain** with paired fecal and cecal specimens (**",
    stats$n_fvc, " samples**; **", stats$n_pain_paired, " subjects** with both sites). Taxonomic profiles were generated with SingleM. Functional profiles used HUMAnN3 MetaCyc pathways collapsed to **",
    nrow(path_mat), " community pathways** (", ncol(fvc_mat), " pain samples with HUMAnN; MAD059C lacked a pathway profile)."
  ),
  "",
  "## Community-level composition",
  "",
  paste0(
    "**Taxonomic beta diversity** (Aitchison/CLR) differed by sampling site (R2 = ", fmt_p(stats$perm_tax_site_R2, 3),
    ", P = ", fmt_p(stats$perm_tax_site_p, 4), "). **Pathway beta diversity** (Bray-Curtis) showed a concordant effect (R2 = ",
    fmt_p(stats$perm_path_site_R2, 3), ", P = ", fmt_p(stats$perm_path_site_p, 4), "). GAD-7 total was associated with pathway communities (P = ",
    fmt_p(perm_p1$p_value[perm_p1$metadata == "gad_total"], 3), ")."
  ),
  "",
  "See `figures/pcoa_pathways_by_site.pdf`.",
  "",
  "## Alpha diversity",
  "",
  paste0("Pathway Shannon paired Wilcoxon P = ", fmt_p(stats$alpha_paired_shannon, 4), ". ", otu_note),
  "",
  "## Differential taxa (fecal vs cecal)",
  "",
  paste0(
    "**", stats$n_consensus_tax_4method, " taxa** were supported by all four methods (MaAsLin2, paired MaAsLin2, paired Wilcoxon, LEfSe-style; q < 0.25). ",
    "Examples: *Phocaeicola vulgatus*, *Alistipes* spp., *Faecalibacterium* spp., *Oscillibacter* spp., *Ruminococcus gnavus*, *Agathobacter rectalis*."
  ),
  "",
  "## Differential pathways",
  "",
  paste0(
    "**", stats$n_consensus_path_2maaslin, " pathways** met MaAsLin2 consensus (unpaired + paired). ",
    "Random forest pathway site classification OOB error = ", fmt_p(stats$rf_path_oob_site, 3), "."
  ),
  "",
  "Strict q < 0.05 tables: `shared/tables/strict_q05_*`. Effect sizes: `tables/top_pathway_effect_sizes.csv`, `tables/top_taxon_effect_sizes.csv`.",
  "",
  "## Limitations",
  "",
  "One cecal sample (MAD059C) lacked HUMAnN output. Taxonomic alpha/PCoA boxplots require `otu_table.csv` in the project root. Pathway and taxonomy PERMANOVA R2 values are not directly comparable (Bray-Curtis vs Aitchison/CLR)."
)

paper2_lines <- c(
  "# Results - Paper 2: Fecal microbiome and clinical metadata",
  "",
  "## Participants",
  "",
  paste0(
    "Fecal cohort: **", stats$n_fps, " participants** (", stats$n_fps_pain, " pain, ",
    stats$n_fps_control, " controls)."
  ),
  "",
  "## Abdominal pain",
  "",
  paste0(
    "No significant pain vs no_pain compositional difference (taxonomy PERMANOVA P = ",
    fmt_p(stats$perm_fecal_pop_p, 3), "; pathway P = ", fmt_p(perm_p2$p_value[perm_p2$metadata == "pop"], 3), ")."
  ),
  "",
  "## Metadata and community structure",
  "",
  paste0(
    "Household sharing was significant in taxonomy PERMANOVA (R2 = ", fmt_p(stats$perm_fecal_hh_R2, 3),
    ", P = ", fmt_p(stats$perm_fecal_hh_p, 3), ") but not at the pathway level (P = ",
    fmt_p(perm_p2$p_value[perm_p2$metadata == "household_shared"], 3), "). BMI trend P = ",
    fmt_p(stats$perm_fecal_bmi_p, 3), "; pathway BMI PERMANOVA P = ",
    fmt_p(perm_p2$p_value[perm_p2$metadata == "bmi"], 3), ". Calprotectin PERMANOVA P = ",
    fmt_p(perm_fecal_tax$p_value[perm_fecal_tax$metadata == "calprotectin"], 3), "."
  ),
  "",
  paste0(
    "After excluding recent antibiotics (n = ", nrow(fps_abx0), "), household sharing P = ",
    fmt_p(perm_abx$p_value[perm_abx$metadata == "household_shared"], 3), " and BMI P = ",
    fmt_p(perm_abx$p_value[perm_abx$metadata == "bmi"], 3), " at the pathway level."
  ),
  "",
  "Sensitivity: `tables/permanova_antibiotics_excluded.csv`, `tables/permanova_pain_only_fecal.csv`.",
  "",
  "## Taxa and clinical metadata (MaAsLin2, q < 0.25)",
  "",
  "Calprotectin: *Pseudoruminococcus massiliensis* and others. GAD difficulty (impairment item): *Frisingicoccus caecimuris*, *Faecousia intestinalis*.",
  "",
  "## Functional pathways",
  "",
  "No multivariable pathway associations at q < 0.25 after full sample mapping; exploratory BMI Spearman hits in `fecal_pathway_spearman_significant.csv`.",
  "",
  paste0("Pathway pain-classification OOB error = ", fmt_p(rf_p2$oob, 3), " (chance 0.50)."),
  "",
  "Exploratory BMI pathway associations: `figures/scatter_pathways_bmi.pdf`. Taxonomic RF importances: `shared/random_forest_taxonomy/`.",
  "",
  "## Limitations",
  "",
  paste0(
    "Multivariable pathway MaAsLin2 had no q < 0.25 hits after sample-ID harmonization; univariate Spearman/BMI findings are exploratory. ",
    "Household effects were stronger for taxonomy than pathways (pathway PERMANOVA P = ",
    fmt_p(perm_p2$p_value[perm_p2$metadata == "household_shared"], 3), ")."
  )
)

writeLines(paper1_lines, file.path(p1_dir, "Results.md"))
writeLines(paper2_lines, file.path(p2_dir, "Results.md"))

writeLines(c(
  "# Methods harmonization (taxonomy vs function)",
  "",
  "| Layer | Distance / model | Transform |",
  "|-------|------------------|-----------|",
  "| Taxonomy (beta) | Aitchison | CLR (microViz) |",
  "| Pathways (beta) | Bray–Curtis | Collapsed MetaCyc rel. abund. |",
  "| MaAsLin2 | LM, NONE/NONE | Raw rel. abund. |",
  "| Alpha (pathways) | Shannon, richness | On collapsed pathways |",
  "",
  "Do not directly compare PERMANOVA R² between taxonomy and pathways."
), file.path(shared, "methods_harmonization.md"))

# Session info
writeLines(capture.output(sessionInfo()), file.path(shared, "session_info.txt"))

# Pandoc to docx
for (pair in list(c(file.path(p1_dir, "Results.md"), file.path(p1_dir, "Results.docx")),
                  c(file.path(p2_dir, "Results.md"), file.path(p2_dir, "Results.docx")))) {
  if (nzchar(Sys.which("pandoc"))) {
    system2("pandoc", c(pair[1], "-o", pair[2], "--from=markdown", "--to=docx"), stdout = FALSE)
  }
}

readme <- c(
  "# MAD manuscript outputs",
  "",
  "## Structure",
  "- `paper1_fecal_vs_cecal/` - pain cohort site comparison (figures, tables, Results.md/.docx)",
  "- `paper2_fecal_metadata/` - fecal cohort metadata (figures, tables, Results.md/.docx)",
  "- `shared/` - sample flow, strict/exploratory MaAsLin tables, pathway PERMANOVA, RF summaries",
  "",
  "## Regenerate",
  "```bash",
  "Rscript build_manuscript_package.R",
  "```",
  "",
  "## Optional inputs",
  "- `otu_table.csv` in the replication kit enables taxonomic alpha diversity and calprotectin taxon boxplots.",
  "- `combined_genefamilies_relab.tsv` not required for current package.",
  "",
  "## Known gaps",
  "- MAD059C missing from HUMAnN combined profiles",
  "- Gene-family and antibiotics-adjusted MaAsLin2 are partial (PERMANOVA sensitivity only)"
)
writeLines(readme, file.path(out_root, "README.md"))

message("Done. Outputs in ", out_root)
