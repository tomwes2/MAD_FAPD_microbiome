#!/usr/bin/env Rscript
# Rebuild Paper 1 figures/tables that appear in tables_and_figures.docx
# but were assembled outside the stripped replication kit.
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
  library(ggpubr)
  library(vegan)
})
source(file.path(data_dir, "R", "microbiome_helpers.R"))
source(file.path(data_dir, "R", "manuscript_utils.R"))
source(file.path(data_dir, "R", "direction_helpers.R"))
if (file.exists(file.path(data_dir, "R", "rf_importance_helpers.R"))) {
  source(file.path(data_dir, "R", "rf_importance_helpers.R"))
}

analysis_dir <- file.path(output_root, "analysis_output")
func_dir <- file.path(output_root, "functional_analysis_output")
p1_fig <- file.path(output_root, "manuscript_outputs", "paper1_fecal_vs_cecal", "figures")
p1_tab <- file.path(output_root, "manuscript_outputs", "paper1_fecal_vs_cecal", "tables")
net_dir <- file.path(func_dir, "pathway_network")
dir.create(p1_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(p1_tab, recursive = TRUE, showWarnings = FALSE)

save_both <- function(plot, base, width, height) {
  ggplot2::ggsave(paste0(base, ".pdf"), plot, width = width, height = height)
  ggplot2::ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 300)
}

normf <- function(x) gsub("[._ ]+", ".", tolower(trimws(x)))
sp_label <- function(f) sub("^OTU[0-9]+[._ ]+", "", gsub("\\.", " ", f))

gi_palette <- c(Fecal = "gold", Cecal = "darkred")

# --- Figure 1: alpha + beta diversity panel ---
message("Figure 1: diversity panel")
mad_objs <- load_mad_phyloseq(data_dir)
fvc_ps <- mad_objs$fvc_ps
fvc_amp <- ampvis2::amp_load(fvc_ps)
fvc_alpha <- ampvis2::amp_alphadiv(fvc_amp, measure = c("uniqueotus", "shannon"))

p_richness <- ggviolin(
  fvc_alpha, x = "sample_type", y = "uniqueOTUs",
  fill = "sample_type", palette = gi_palette, add = "boxplot"
) + stat_compare_means(method = "wilcox.test") +
  ggtitle("Richness") + theme(legend.position = "none")
p_shannon <- ggviolin(
  fvc_alpha, x = "sample_type", y = "Shannon",
  fill = "sample_type", palette = gi_palette, add = "boxplot"
) + stat_compare_means(method = "wilcox.test") +
  ggtitle("Shannon diversity") + theme(legend.position = "none")

perm_csv <- file.path(analysis_dir, "permanova_pain_fecal_vs_cecal.csv")
if (file.exists(perm_csv)) {
  perm_tab <- read.csv(perm_csv, stringsAsFactors = FALSE)
  perm_p <- perm_tab[[grep("^Pr", names(perm_tab), value = TRUE)[1]]][1]
  pca_sub <- sprintf(
    "PERMANOVA: R2 = %.3f, p %s",
    perm_tab$R2[1],
    format.pval(perm_p, digits = 2, eps = 0.001)
  )
} else {
  pca_sub <- "CLR PCA (species-level)"
}
p_pca <- fvc_ps %>%
  microViz::tax_agg("species") %>%
  microViz::tax_transform("clr") %>%
  microViz::ord_calc() %>%
  microViz::ord_plot(plot_taxa = 1:5, colour = "sample_type", size = 1.5) +
  ggplot2::stat_ellipse(ggplot2::aes(colour = sample_type), linewidth = 0.3) +
  ggplot2::scale_color_manual(values = gi_palette, name = "GI site") +
  ggplot2::labs(title = "Beta diversity (CLR PCA)", subtitle = pca_sub) +
  ggplot2::theme(legend.position = "bottom")

fig1 <- ggarrange(
  p_richness, p_shannon, p_pca,
  ncol = 3, nrow = 1, labels = c("A", "B", "C"),
  font.label = list(size = 14, face = "bold"),
  widths = c(1, 1, 1.2)
)
save_both(fig1, file.path(p1_fig, "figure1_diversity_panel"), width = 14, height = 5)
save_both(fig1, file.path(analysis_dir, "diversity_panel_pain_gi"), width = 14, height = 5)

# --- Table 1: MaAsLin2 + paired Wilcoxon concordant taxa ---
message("Table 1: concordant taxa")
maas <- read.delim(file.path(analysis_dir, "maaslin2/pain_fecal_vs_cecal/all_results.tsv"), check.names = FALSE) %>%
  annotate_maaslin_direction() %>%
  filter(metadata == "sample_type", value == "Cecal")
wilcox <- read.csv(file.path(analysis_dir, "paired_wilcox_pain_gi.csv"), stringsAsFactors = FALSE)
if (!"enriched_in" %in% names(wilcox)) wilcox <- annotate_paired_wilcox_direction(wilcox)
maas$feature_norm <- normf(maas$feature)
wilcox$feature_norm <- normf(wilcox$feature)
tab1 <- inner_join(
  maas %>% transmute(feature_norm, maaslin_feature = feature, maaslin_coef = coef, maaslin_q = qval, maas_site = enriched_in),
  wilcox %>% transmute(feature_norm, wilcox_feature = feature, median_diff = median_diff, wilcox_q = qval, wilcox_site = enriched_in),
  by = "feature_norm"
) %>%
  filter(!is.na(maas_site), maas_site == wilcox_site, maas_site %in% c("Cecal", "Fecal")) %>%
  mutate(species = vapply(maaslin_feature, sp_label, character(1))) %>%
  arrange(maas_site, desc(abs(maaslin_coef))) %>%
  select(species, enriched_in = maas_site, maaslin_coef, median_diff, maaslin_q, wilcox_q)
write.csv(tab1, file.path(p1_tab, "table1_concordant_taxa.csv"), row.names = FALSE)

# --- Figure 2: strong consensus |MaAsLin2 coef| ---
message("Figure 2: strong consensus taxa")
strong <- read.csv(file.path(analysis_dir, "consensus/consensus_pain_gi_site_strong.csv"), stringsAsFactors = FALSE)
plot_df <- strong %>%
  mutate(feature_norm = normf(feature)) %>%
  left_join(
    maas %>% group_by(feature_norm) %>% summarise(coef = coef[1], enriched_in = enriched_in[1], .groups = "drop"),
    by = "feature_norm"
  ) %>%
  mutate(
    species = vapply(feature, sp_label, character(1)),
    effect = abs(coef)
  ) %>%
  filter(!is.na(effect))
p2 <- ggplot(plot_df, aes(x = reorder(species, effect), y = effect, fill = enriched_in)) +
  geom_col(width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(Cecal = "#2166AC", Fecal = "#B2182B"), name = "Enriched in") +
  labs(
    title = "Strong consensus taxa (>=3 differential-abundance methods)",
    subtitle = "Bar height = |MaAsLin2 coefficient| (Cecal vs Fecal reference); paired Wilcoxon direction concordant",
    x = NULL, y = "|MaAsLin2 coefficient|"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
save_both(p2, file.path(p1_fig, "figure2_strong_consensus_direction"), width = 9, height = 6)

# --- Copy network figures if present (Figure X / Y) ---
copy_if <- function(src, dest) {
  if (file.exists(src)) file.copy(src, dest, overwrite = TRUE)
}
copy_if(file.path(net_dir, "community_by_enrichment.png"), file.path(p1_fig, "figureX_louvain_modules.png"))
copy_if(file.path(net_dir, "community_by_enrichment.pdf"), file.path(p1_fig, "figureX_louvain_modules.pdf"))
copy_if(file.path(net_dir, "network_full.png"), file.path(p1_fig, "figureY_pathway_network.png"))
copy_if(file.path(net_dir, "network_full.pdf"), file.path(p1_fig, "figureY_pathway_network.pdf"))
copy_if(file.path(net_dir, "community_labels.csv"), file.path(p1_tab, "tableY_louvain_modules.csv"))
copy_if(file.path(net_dir, "edge_type_counts.csv"), file.path(p1_tab, "tableX_network_edge_counts.csv"))

# --- Table 2 from RF permutation / drop-column if those CSVs exist ---
rf_tbl <- file.path(analysis_dir, "random_forest", "rf_de_consensus_table_strong.csv")
if (file.exists(rf_tbl)) {
  file.copy(rf_tbl, file.path(p1_tab, "table2_rf_strong_consensus.csv"), overwrite = TRUE)
} else {
  perm_f <- file.path(analysis_dir, "random_forest", "rf_permutation_importance.csv")
  drop_f <- file.path(analysis_dir, "random_forest", "rf_drop_column_importance.csv")
  dep_f <- file.path(analysis_dir, "random_forest", "rf_feature_dependence.csv")
  strong_f <- file.path(analysis_dir, "consensus", "consensus_pain_gi_site_strong.csv")
  if (file.exists(perm_f) && file.exists(strong_f)) {
    strong <- read.csv(strong_f, stringsAsFactors = FALSE)
    perm <- read.csv(perm_f, stringsAsFactors = FALSE)
    strong$feature_norm <- normf(strong$feature)
    perm$feature_norm <- normf(perm$taxon)
    tab2 <- merge(strong, perm[, c("feature_norm", "perm_importance")], by = "feature_norm", all.x = TRUE)
    if (file.exists(drop_f)) {
      drop <- read.csv(drop_f, stringsAsFactors = FALSE)
      drop$feature_norm <- normf(drop$taxon)
      tab2 <- merge(tab2, drop[, c("feature_norm", "drop_col_importance")], by = "feature_norm", all.x = TRUE)
    }
    if (file.exists(dep_f)) {
      dep <- read.csv(dep_f, stringsAsFactors = FALSE)
      dep$feature_norm <- normf(dep$taxon)
      tab2 <- merge(tab2, dep[, c("feature_norm", "dependence_r2_oob")], by = "feature_norm", all.x = TRUE)
    }
    tab2$species <- vapply(tab2$feature, sp_label, character(1))
    tab2$in_rf_model <- !is.na(tab2$perm_importance)
    write.csv(tab2, file.path(p1_tab, "table2_rf_strong_consensus.csv"), row.names = FALSE)
  }
}

message("Paper 1 figures written to ", p1_fig)
