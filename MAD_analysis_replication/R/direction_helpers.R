# Direction helpers for differential abundance results
#
# HOW TO READ GI-SITE COMPARISONS (fecal vs cecal)
# ------------------------------------------------
# MaAsLin2 models use reference = "sample_type,Fecal" (stool is the baseline).
# Result rows with metadata == "sample_type" and value == "Cecal" report the
# Cecal-vs-Fecal contrast:
#   coef > 0  -> pathway/taxon is HIGHER in CECUM than in FECAL (stool)
#   coef < 0  -> pathway/taxon is HIGHER in FECAL (stool) than in CECUM
#
# Paired Wilcoxon uses the same factor order (Fecal, then Cecal) and stores:
#   median_diff = median(cecal abundance - fecal abundance) within subject
#   median_diff > 0 -> higher in cecum; < 0 -> higher in fecal
#
# LEfSe-style Kruskal screens only test "difference by group" — they do not
# assign direction unless paired with MaAsLin or Wilcoxon.

#' Null-coalescing operator (used like x %||% default).
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Normalize feature names for matching across tools (MaAsLin dots vs HUMAnN colons/dashes).
normalize_feature <- function(x) {
  x <- tolower(trimws(x))
  gsub("[^a-z0-9]+", ".", x)
}

#' Short label for which GI site has higher abundance.
site_enriched_label <- function(coef_or_diff, tol = 1e-12) {
  ifelse(
    is.na(coef_or_diff), NA_character_,
    ifelse(coef_or_diff > tol, "Cecal", ifelse(coef_or_diff < -tol, "Fecal", "Tie/near-zero"))
  )
}

#' Plain-English direction for fecal vs cecal comparisons.
site_direction_label <- function(coef_or_diff, tol = 1e-12) {
  enriched <- site_enriched_label(coef_or_diff, tol)
  ifelse(
    enriched == "Cecal", "Higher in cecum than fecal (stool)",
    ifelse(enriched == "Fecal", "Higher in fecal (stool) than cecum", "No clear direction")
  )
}

#' Add direction columns to MaAsLin2 result tables (all_results or significant_results).
#'
#' @param df data.frame with columns metadata, value, coef
#' @return same data.frame with enriched_in, direction_label, coef_interpretation
annotate_maaslin_direction <- function(df) {
  if (!nrow(df)) return(df)
  required <- c("metadata", "value", "coef")
  if (!all(required %in% names(df))) {
    warning("annotate_maaslin_direction: missing columns; returning df unchanged")
    return(df)
  }

  df$enriched_in <- NA_character_
  df$direction_label <- NA_character_
  df$coef_interpretation <- NA_character_

  # GI site: only Cecal contrast rows are interpretable (Fecal is reference).
  idx_site <- df$metadata == "sample_type" & df$value == "Cecal"
  df$enriched_in[idx_site] <- site_enriched_label(df$coef[idx_site])
  df$direction_label[idx_site] <- site_direction_label(df$coef[idx_site])
  df$coef_interpretation[idx_site] <- paste0(
    "MaAsLin2 coef (Cecal vs Fecal reference): ",
    signif(df$coef[idx_site], 3)
  )

  # Pain vs control (when reference is no_pain).
  idx_pain <- df$metadata == "pop" & df$value == "pain"
  df$enriched_in[idx_pain] <- ifelse(df$coef[idx_pain] > 0, "pain", "no_pain")
  df$direction_label[idx_pain] <- ifelse(
    df$coef[idx_pain] > 0,
    "Higher in abdominal-pain group than controls",
    "Higher in controls than abdominal-pain group"
  )
  df$coef_interpretation[idx_pain] <- "MaAsLin2 coef (pain vs no_pain reference)"

  # Binary / categorical metadata (value is a level name).
  idx_cat <- df$metadata %in% c(
    "calprotectin_high", "sex", "phq_difficulty", "gad_difficulty", "household_shared", "obese"
  ) & !idx_site & !idx_pain
  if (any(idx_cat)) {
    df$enriched_in[idx_cat] <- paste0(df$metadata[idx_cat], "=", df$value[idx_cat])
    df$direction_label[idx_cat] <- ifelse(
      df$coef[idx_cat] > 0,
      paste0("Higher abundance when ", df$metadata[idx_cat], " = ", df$value[idx_cat]),
      paste0("Lower abundance when ", df$metadata[idx_cat], " = ", df$value[idx_cat])
    )
    df$coef_interpretation[idx_cat] <- paste0("MaAsLin2 coef for level ", df$value[idx_cat])
  }

  # Continuous metadata (value column repeats variable name, e.g. value == "bmi").
  cont_vars <- c("calprotectin", "age", "bmi", "phq_total", "gad_total")
  idx_cont <- df$metadata %in% cont_vars & (df$value == df$metadata | is.na(df$value))
  if (any(idx_cont)) {
    df$enriched_in[idx_cont] <- ifelse(df$coef[idx_cont] > 0, "increases_with_metadata", "decreases_with_metadata")
    df$direction_label[idx_cont] <- ifelse(
      df$coef[idx_cont] > 0,
      paste0("Higher abundance with increasing ", df$metadata[idx_cont]),
      paste0("Lower abundance with increasing ", df$metadata[idx_cont])
    )
    df$coef_interpretation[idx_cont] <- paste0(
      "MaAsLin2 coef per unit ", df$metadata[idx_cont], " (model-standardized if enabled)"
    )
  }

  df
}

#' Direction from paired Wilcoxon (median cecal - median fecal within subject).
annotate_paired_wilcox_direction <- function(df) {
  if (!nrow(df) || !"median_diff" %in% names(df)) return(df)
  df$enriched_in <- site_enriched_label(df$median_diff)
  df$direction_label <- site_direction_label(df$median_diff)
  df$coef_interpretation <- paste0(
    "Paired median difference (cecal - fecal): ", signif(df$median_diff, 3)
  )
  df
}

#' Direction from Spearman rho (exploratory univariate).
annotate_spearman_direction <- function(df) {
  if (!nrow(df) || !"rho" %in% names(df)) return(df)
  df$direction_label <- ifelse(
    df$rho > 0,
    paste0("Increases with ", df$metadata),
    paste0("Decreases with ", df$metadata)
  )
  df$enriched_in <- ifelse(df$rho > 0, "positive_association", "negative_association")
  df
}

#' Merge MaAsLin + Wilcoxon direction onto a consensus feature table.
add_consensus_direction <- function(consensus_df, maaslin_df, wilcox_df = NULL) {
  if (!nrow(consensus_df)) return(consensus_df)

  if (!"feature_norm" %in% names(consensus_df)) {
    consensus_df$feature_norm <- normalize_feature(consensus_df$feature)
  }

  maas_site <- maaslin_df %>%
    dplyr::filter(metadata == "sample_type", value == "Cecal") %>%
    dplyr::mutate(feature_norm = normalize_feature(feature)) %>%
    dplyr::select(feature_norm, maaslin_coef = coef, maaslin_q = qval,
                  maaslin_enriched_in = enriched_in, maaslin_direction = direction_label)

  out <- consensus_df %>%
    dplyr::left_join(maas_site, by = "feature_norm")

  if (!is.null(wilcox_df) && nrow(wilcox_df)) {
    wx <- wilcox_df %>%
      dplyr::mutate(feature_norm = normalize_feature(feature)) %>%
      dplyr::select(feature_norm, wilcox_median_diff = median_diff,
                    wilcox_enriched_in = enriched_in, wilcox_direction = direction_label)
    out <- out %>% dplyr::left_join(wx, by = "feature_norm")
  }

  # Consensus call: prefer MaAsLin; fall back to Wilcoxon if MaAsLin missing.
  out$enriched_in <- out$maaslin_enriched_in
  if ("wilcox_enriched_in" %in% names(out)) {
    miss <- is.na(out$enriched_in)
    out$enriched_in[miss] <- out$wilcox_enriched_in[miss]
  }
  out$direction_label <- out$maaslin_direction
  if ("wilcox_direction" %in% names(out)) {
    miss <- is.na(out$direction_label)
    out$direction_label[miss] <- out$wilcox_direction[miss]
  }
  out$direction_agreement <- NA_character_
  if ("wilcox_enriched_in" %in% names(out)) {
    both <- !is.na(out$maaslin_enriched_in) & !is.na(out$wilcox_enriched_in)
    agree <- both & out$maaslin_enriched_in == out$wilcox_enriched_in
    disagree <- both & out$maaslin_enriched_in != out$wilcox_enriched_in
    out$direction_agreement[agree] <- "MaAsLin and Wilcoxon agree"
    out$direction_agreement[disagree] <- "MaAsLin and Wilcoxon disagree — inspect abundances"
  }
  out
}

#' Write annotated MaAsLin TSV back to disk (optional convenience).
write_maaslin_with_direction <- function(path_in, path_out = path_in) {
  if (!file.exists(path_in)) return(invisible(FALSE))
  df <- read.delim(path_in, check.names = FALSE)
  df <- annotate_maaslin_direction(df)
  write.table(df, path_out, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(TRUE)
}
