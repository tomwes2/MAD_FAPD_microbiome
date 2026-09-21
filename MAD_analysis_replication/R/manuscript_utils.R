# Shared helpers for MAD manuscript package (tables, plots, HUMAnN loading).
# For differential *direction* (cecum vs fecal), see R/direction_helpers.R.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Map HUMAnN file column names (e.g. A1_clean_Abundance) to MAD sample IDs via sample_ids.csv.

load_sample_id_map <- function(csv_path) {
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

pathway_id <- function(feature) sub("\\|.*", "", feature)

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

# Recode metadata for stats: factor level order for sample_type is Fecal then Cecal (Wilcoxon sign).
prepare_meta <- function(meta) {
  meta$obese <- factor(ifelse(meta$bmi > 30, "obese", "non_obese"), levels = c("non_obese", "obese"))
  meta$calprotectin_binary <- factor(
    ifelse(meta$calprotectin_high %in% c(1, "1", TRUE), "high", "normal"),
    levels = c("normal", "high")
  )
  meta$pop <- factor(meta$pop, levels = c("no_pain", "pain"))
  meta$sample_type <- factor(meta$sample_type, levels = c("Fecal", "Cecal"))
  meta$sex <- factor(meta$sex)
  meta$phq_difficulty <- factor(meta$phq_difficulty, levels = 0:3)
  meta$gad_difficulty <- factor(meta$gad_difficulty, levels = 0:3)
  meta$household_id <- ifelse(
    meta$household == 0, paste0("unique_", meta$ssn), paste0("hh_", meta$household)
  )
  meta$household_shared <- factor(
    ifelse(meta$household == 0, "no_shared", "shared"), levels = c("no_shared", "shared")
  )
  meta$antibiotics_recent <- factor(ifelse(meta$antibiotics %in% c(1, "1", TRUE), "yes", "no"))
  meta
}

shannon <- function(x) {
  x <- x[x > 0]
  if (!length(x)) return(0)
  p <- x / sum(x)
  -sum(p * log(p))
}

alpha_metrics <- function(mat) {
  tibble(
    sampleID = colnames(mat),
    richness = colSums(mat > 0),
    shannon = vapply(seq_len(ncol(mat)), function(i) shannon(mat[, i]), numeric(1))
  )
}

permanova_one <- function(mat, meta, var, label) {
  common <- intersect(colnames(mat), rownames(meta))
  d <- vegdist(t(mat[, common, drop = FALSE]), method = "bray")
  f <- as.formula(paste("d ~", var))
  fit <- adonis2(f, data = meta[common, , drop = FALSE], permutations = 9999)
  tibble(cohort = label, metadata = var, R2 = fit$R2[1], F = fit$F[1], p_value = fit$`Pr(>F)`[1])
}

taxon_label <- function(x) {
  gsub("_", " ", sub("^OTU[0-9]+_", "", x))
}

# MaAsLin2 replaces special characters in feature names; match back to HUMAnN rownames.
match_maaslin_feature <- function(feature, rownames_vec) {
  if (feature %in% rownames_vec) return(feature)
  idx <- match(make.names(feature), make.names(rownames_vec))
  if (!is.na(idx)) return(rownames_vec[idx])
  NA_character_
}

read_maaslin_all <- function(path) {
  if (!file.exists(path)) return(tibble())
  read.delim(path, check.names = FALSE) %>% as_tibble()
}

strict_maaslin <- function(df, q = 0.05) {
  if (!nrow(df)) return(df)
  df %>% filter(qval < q)
}

save_plot <- function(p, path_base, width = 8, height = 6) {
  ggplot2::ggsave(paste0(path_base, ".pdf"), p, width = width, height = height)
  ggplot2::ggsave(paste0(path_base, ".png"), p, width = width, height = height, dpi = 300)
}
