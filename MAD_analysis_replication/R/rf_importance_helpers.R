# Random forest feature importance — permutation, drop-column, collinearity.
# Methods follow https://explained.ai/rf-importance/index.html

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Save a ggplot to PDF and high-resolution PNG.
save_gg_pdf_png <- function(plot, base_path, width = 8, height = 7, dpi = 300) {
  base <- sub("\\.(pdf|png)$", "", base_path, ignore.case = TRUE)
  pdf_path <- paste0(base, ".pdf")
  png_path <- paste0(base, ".png")
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi)
  invisible(list(pdf = pdf_path, png = png_path))
}

#' Save a ComplexHeatmap (or other grid draw) to PDF and PNG.
save_heatmap_pdf_png <- function(ht, base_path, width = 10, height = 8, dpi = 300) {
  base <- sub("\\.(pdf|png)$", "", base_path, ignore.case = TRUE)
  pdf_path <- paste0(base, ".pdf")
  png_path <- paste0(base, ".png")
  grDevices::pdf(pdf_path, width = width, height = height)
  ComplexHeatmap::draw(ht, merge_legends = TRUE)
  grDevices::dev.off()
  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
  ComplexHeatmap::draw(ht, merge_legends = TRUE)
  grDevices::dev.off()
  invisible(list(pdf = pdf_path, png = png_path))
}

#' Build RF feature matrix from phyloseq (top taxa by prevalence).
prepare_rf_matrix <- function(ps, response_var, n_top_taxa = 100, add_random_feature = TRUE,
                              random_feature_name = "random_feature") {
  mat <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) mat <- t(mat)
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  y <- meta[[response_var]]
  if (length(unique(na.omit(y))) < 2) {
    return(NULL)
  }
  prev <- rowMeans(mat > 0)
  keep_taxa <- names(sort(prev, decreasing = TRUE))[seq_len(min(n_top_taxa, length(prev)))]
  x <- t(mat[keep_taxa, , drop = FALSE])
  safe_names <- make.names(colnames(x), unique = TRUE)
  name_map <- stats::setNames(colnames(x), safe_names)
  colnames(x) <- safe_names
  df <- data.frame(y = factor(y), x, check.names = FALSE)
  if (add_random_feature) {
    set.seed(42)
    df[[random_feature_name]] <- stats::runif(nrow(df))
    name_map[random_feature_name] <- random_feature_name
  }
  df <- df[stats::complete.cases(df), ]
  list(
    df = df,
    name_map = name_map,
    taxa_safe = setdiff(names(name_map), random_feature_name),
    random_feature = random_feature_name
  )
}

#' Train random forest classifier; returns model + importance table.
fit_rf_classifier <- function(rf_data, label, response_var, ntree = 2000, seed = 42) {
  if (is.null(rf_data)) return(NULL)
  df <- rf_data$df
  set.seed(seed)
  rf <- randomForest::randomForest(y ~ ., data = df, ntree = ntree, importance = TRUE)
  imp_gini <- as.data.frame(randomForest::importance(rf, type = 2))
  imp_perm <- as.data.frame(randomForest::importance(rf, type = 1))
  imp <- imp_gini
  imp$MeanDecreaseAccuracy <- imp_perm$MeanDecreaseAccuracy
  imp$taxon <- unname(rf_data$name_map[rownames(imp)])
  imp$cohort <- label
  imp$response <- response_var
  imp$is_random_feature <- imp$taxon == rf_data$random_feature
  oob_err <- rf$err.rate[nrow(rf$err.rate), "OOB"]
  list(
    model = rf,
    importance = imp,
    oob_error = oob_err,
    oob_accuracy = 1 - oob_err,
    rf_data = rf_data
  )
}

#' Replication-based SE for permutation importance (stability across forest fits).
rf_permutation_importance_se <- function(rf_data, n_repeats = 10, ntree = 500, seed = 42) {
  if (is.null(rf_data)) return(NULL)
  df <- rf_data$df
  feats <- rownames(rf_data$name_map)
  mats <- vector("list", n_repeats)
  for (i in seq_len(n_repeats)) {
    set.seed(seed + i)
    rf <- randomForest::randomForest(y ~ ., data = df, ntree = ntree, importance = TRUE)
    imp <- as.data.frame(randomForest::importance(rf, type = 1))
    mats[[i]] <- imp$MeanDecreaseAccuracy
    names(mats[[i]]) <- rownames(imp)
  }
  mat <- do.call(cbind, mats)
  tibble::tibble(
    feature_safe = rownames(mat),
    taxon = unname(rf_data$name_map[rownames(mat)]),
    perm_importance = rowMeans(mat),
    perm_importance_sd = apply(mat, 1, stats::sd),
    is_random_feature = taxon == rf_data$random_feature
  ) %>%
    dplyr::arrange(dplyr::desc(perm_importance))
}

#' Drop-column importance: retrain without each feature; positive = feature helps OOB accuracy.
rf_drop_column_importance <- function(rf_data, features = NULL, ntree = 500, seed = 42,
                                      n_cores = 1L) {
  if (is.null(rf_data)) return(NULL)
  df <- rf_data$df
  y_col <- "y"
  if (is.null(features)) {
    features <- setdiff(names(df), y_col)
  } else {
    features <- intersect(features, setdiff(names(df), y_col))
  }
  form_full <- stats::as.formula(paste(y_col, "~ ."))
  set.seed(seed)
  rf_full <- randomForest::randomForest(form_full, data = df, ntree = ntree, importance = FALSE)
  baseline_oob <- rf_full$err.rate[nrow(rf_full$err.rate), "OOB"]

  drop_one <- function(feat) {
    cols <- c(y_col, setdiff(names(df), c(y_col, feat)))
    sub_df <- df[, cols, drop = FALSE]
    set.seed(seed)
    rf_drop <- randomForest::randomForest(
      stats::as.formula(paste(y_col, "~ .")),
      data = sub_df, ntree = ntree, importance = FALSE
    )
    oob_drop <- rf_drop$err.rate[nrow(rf_drop$err.rate), "OOB"]
    tibble::tibble(
      feature_safe = feat,
      taxon = unname(rf_data$name_map[feat]),
      oob_error_full = baseline_oob,
      oob_error_dropped = oob_drop,
      drop_col_importance = oob_drop - baseline_oob,
      is_random_feature = unname(rf_data$name_map[feat]) == rf_data$random_feature
    )
  }

  if (n_cores > 1L && .Platform$OS.type != "windows") {
    res <- parallel::mclapply(features, drop_one, mc.cores = min(n_cores, length(features)))
  } else {
    res <- lapply(features, drop_one)
  }
  dplyr::bind_rows(res) %>% dplyr::arrange(dplyr::desc(drop_col_importance))
}

#' Spearman correlation among RF features (samples x features matrix).
rf_spearman_correlation <- function(rf_data, features = NULL) {
  if (is.null(rf_data)) return(NULL)
  df <- rf_data$df
  x_cols <- setdiff(names(df), "y")
  if (!is.null(features)) x_cols <- intersect(features, x_cols)
  x_mat <- as.matrix(df[, x_cols, drop = FALSE])
  cor_mat <- stats::cor(x_mat, method = "spearman", use = "pairwise.complete.obs")
  cor_df <- as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
  names(cor_df) <- c("feature1", "feature2", "spearman_rho")
  taxon_map <- rf_data$name_map
  cor_df$taxon1 <- unname(taxon_map[cor_df$feature1])
  cor_df$taxon2 <- unname(taxon_map[cor_df$feature2])
  list(matrix = cor_mat, long = cor_df)
}

#' Feature dependence: OOB R² predicting each feature from all others (collinearity diagnostic).
rf_feature_dependence <- function(rf_data, features = NULL, ntree = 500, seed = 42) {
  if (is.null(rf_data)) return(NULL)
  df <- rf_data$df
  predictors <- setdiff(names(df), "y")
  if (!is.null(features)) predictors <- intersect(features, predictors)
  dep_one <- function(feat) {
    others <- setdiff(names(df), c("y", feat))
    sub_df <- df[, c(feat, others), drop = FALSE]
    names(sub_df)[1] <- "target"
    set.seed(seed)
    rf <- randomForest::randomForest(target ~ ., data = sub_df, ntree = ntree, importance = FALSE)
    oob_mse <- mean((df[[feat]] - rf$predicted)^2, na.rm = TRUE)
    var_y <- stats::var(df[[feat]], na.rm = TRUE)
    r2_oob <- if (var_y > 0) 1 - oob_mse / var_y else NA_real_
    tibble::tibble(
      feature_safe = feat,
      taxon = unname(rf_data$name_map[feat]),
      dependence_r2_oob = r2_oob,
      is_random_feature = unname(rf_data$name_map[feat]) == rf_data$random_feature
    )
  }
  dplyr::bind_rows(lapply(predictors, dep_one)) %>%
    dplyr::arrange(dplyr::desc(dependence_r2_oob))
}

#' Hierarchical clusters of correlated features (|rho| >= threshold).
rf_correlation_clusters <- function(cor_mat, threshold = 0.7) {
  if (is.null(cor_mat) || nrow(cor_mat) < 2) return(list(clusters = list(), membership = character()))
  d <- stats::as.dist(1 - abs(cor_mat))
  hc <- stats::hclust(d, method = "average")
  clusters <- stats::cutree(hc, h = 1 - threshold)
  cluster_list <- split(names(clusters), clusters)
  list(clusters = cluster_list, membership = clusters, hclust = hc)
}

#' Grouped permutation importance: permute correlated clusters together.
rf_grouped_permutation_importance <- function(rf_data, clusters, ntree = 500, seed = 42) {
  if (is.null(rf_data) || length(clusters) == 0) return(NULL)
  df <- rf_data$df
  set.seed(seed)
  rf_full <- randomForest::randomForest(y ~ ., data = df, ntree = ntree, importance = FALSE)
  baseline_oob <- rf_full$err.rate[nrow(rf_full$err.rate), "OOB"]

  perm_group <- function(group_feats) {
    df_perm <- df
    for (f in group_feats) {
      df_perm[[f]] <- sample(df_perm[[f]])
    }
    set.seed(seed)
    rf_perm <- randomForest::randomForest(y ~ ., data = df_perm, ntree = ntree, importance = FALSE)
    oob_perm <- rf_perm$err.rate[nrow(rf_perm$err.rate), "OOB"]
    tibble::tibble(
      group_id = paste(group_feats, collapse = " + "),
      n_features = length(group_feats),
      taxa = paste(unname(rf_data$name_map[group_feats]), collapse = " + "),
      grouped_perm_importance = oob_perm - baseline_oob
    )
  }
  dplyr::bind_rows(lapply(clusters, perm_group)) %>%
    dplyr::arrange(dplyr::desc(grouped_perm_importance))
}

plot_rf_importance_perm <- function(rf_list, top_n = 20, metric = "MeanDecreaseAccuracy") {
  imp <- rf_list$importance %>%
    dplyr::filter(!is_random_feature) %>%
    dplyr::arrange(dplyr::desc(.data[[metric]])) %>%
    dplyr::slice_head(n = top_n)
  ylab <- if (metric == "MeanDecreaseAccuracy") {
    "Permutation importance (mean decrease accuracy)"
  } else {
    "Mean decrease Gini"
  }
  ggplot2::ggplot(imp, ggplot2::aes(x = reorder(taxon, .data[[metric]]), y = .data[[metric]])) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("RF importance - ", rf_list$importance$cohort[1]),
      subtitle = sprintf(
        "OOB accuracy: %.1f%% | metric: %s",
        100 * rf_list$oob_accuracy,
        if (metric == "MeanDecreaseAccuracy") "permutation" else "Gini"
      ),
      x = NULL, y = ylab
    ) +
    ggplot2::theme_bw()
}

plot_rf_importance_comparison <- function(imp_combined, top_n = 20) {
  imp_combined %>%
    dplyr::filter(!is_random_feature) %>%
    dplyr::slice_head(n = top_n) %>%
    tidyr::pivot_longer(
      c(MeanDecreaseGini, MeanDecreaseAccuracy, drop_col_importance),
      names_to = "metric", values_to = "value"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        MeanDecreaseGini = "Gini",
        MeanDecreaseAccuracy = "Permutation",
        drop_col_importance = "Drop-column"
      ),
      metric = factor(metric, levels = c("Gini", "Permutation", "Drop-column"))
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(taxon, value), y = value, fill = metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "RF importance comparison (top taxa)",
      subtitle = "Permutation and drop-column preferred over Gini (Explained.ai)",
      x = NULL, y = "Importance", fill = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

plot_rf_collinearity_heatmap <- function(cor_mat, taxon_labels, title = "Spearman correlation") {
  if (is.null(cor_mat) || nrow(cor_mat) < 2) return(NULL)
  labels <- taxon_labels[rownames(cor_mat)]
  labels[is.na(labels)] <- rownames(cor_mat)[is.na(labels)]
  ord <- stats::hclust(stats::as.dist(1 - abs(cor_mat)), method = "average")$order
  cor_ord <- cor_mat[ord, ord, drop = FALSE]
  lab_ord <- labels[ord]
  ComplexHeatmap::Heatmap(
    cor_ord,
    name = "Spearman",
    col = viridisLite::viridis(100),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_side = "top",
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    row_labels = lab_ord,
    column_labels = lab_ord,
    column_title = title,
    rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
}

normalize_rf_feature <- function(x) {
  gsub("[._ ]+", ".", tolower(trimws(x)))
}

species_display_name <- function(feature) {
  x <- gsub("\\.", " ", as.character(feature))
  sub("^OTU[0-9]+[._ ]+", "", x)
}

#' Integrate RF importance metrics with differential-abundance consensus taxa.
#'
#' Builds a publication-oriented table: DE method support, three RF metrics,
#' robustness score, and collinearity / module annotations.
build_rf_de_consensus_table <- function(
    consensus_df,
    rf_importance_combined,
    cor_long = NULL,
    cor_clusters = NULL,
    grouped_perm = NULL,
    maaslin_hits = NULL,
    wilcox_hits = NULL,
    rf_top_n = 30,
    cor_threshold = 0.7,
    dependence_threshold = 0.35,
    min_de_methods = 2) {
  if (!nrow(consensus_df)) return(tibble::tibble())

  de <- consensus_df %>%
    dplyr::filter(n_methods >= min_de_methods) %>%
    dplyr::mutate(
      feature_norm = normalize_rf_feature(feature),
      species = vapply(feature, species_display_name, character(1))
    )

  rf <- rf_importance_combined %>%
    dplyr::filter(!is_random_feature) %>%
    dplyr::mutate(feature_norm = normalize_rf_feature(taxon)) %>%
    dplyr::group_by(feature_norm) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  rf_ranked <- rf %>%
    dplyr::mutate(
      rank_permutation = rank(-MeanDecreaseAccuracy, ties.method = "min"),
      rank_gini = rank(-MeanDecreaseGini, ties.method = "min"),
      rank_drop_col = dplyr::if_else(
        is.na(drop_col_importance),
        NA_real_,
        rank(-drop_col_importance, ties.method = "min")
      ),
      perm_top = rank_permutation <= rf_top_n,
      gini_top = rank_gini <= rf_top_n,
      drop_col_tested = !is.na(drop_col_importance),
      drop_col_support = dplyr::if_else(drop_col_tested, drop_col_importance > 0, NA)
    )

  maaslin_dir_tbl <- if (!is.null(maaslin_hits) && nrow(maaslin_hits)) {
    maaslin_hits %>%
      dplyr::mutate(feature_norm = normalize_rf_feature(feature)) %>%
      dplyr::group_by(feature_norm) %>%
      dplyr::summarise(
        maaslin_coef = coef[which.min(qval)][1],
        maaslin_qval = min(qval, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        enriched_in = dplyr::if_else(
          maaslin_coef > 0, "Cecal",
          dplyr::if_else(maaslin_coef < 0, "Fecal", "Tie")
        )
      )
  } else {
    tibble::tibble(
      feature_norm = character(), maaslin_coef = numeric(),
      maaslin_qval = numeric(), enriched_in = character()
    )
  }

  wilcox_dir_tbl <- if (!is.null(wilcox_hits) && nrow(wilcox_hits)) {
    wilcox_hits %>%
      dplyr::mutate(feature_norm = normalize_rf_feature(feature)) %>%
      dplyr::group_by(feature_norm) %>%
      dplyr::summarise(
        wilcox_median_diff = median_diff[1],
        wilcox_qval = min(qval, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble::tibble(
      feature_norm = character(), wilcox_median_diff = numeric(), wilcox_qval = numeric()
    )
  }

  partner_lookup <- function(feature_norm, restrict_norms = NULL) {
    if (is.null(cor_long) || !nrow(cor_long)) return(character())
    partners <- cor_long %>%
      dplyr::mutate(
        f1 = normalize_rf_feature(taxon1),
        f2 = normalize_rf_feature(taxon2)
      ) %>%
      dplyr::filter(
        abs(spearman_rho) >= cor_threshold,
        f1 == feature_norm | f2 == feature_norm
      ) %>%
      dplyr::mutate(
        partner_norm = dplyr::if_else(f1 == feature_norm, f2, f1),
        partner_label = dplyr::if_else(f1 == feature_norm, taxon2, taxon1)
      )
    partners <- partners %>% dplyr::filter(partner_norm != feature_norm)
    if (!is.null(restrict_norms)) {
      partners <- partners %>% dplyr::filter(partner_norm %in% restrict_norms)
    }
    partners %>%
      dplyr::distinct(partner_label) %>%
      dplyr::pull(partner_label) %>%
      sort()
  }

  cluster_id_for <- function(feature_norm) {
    if (is.null(cor_clusters) || !length(cor_clusters$membership)) return(NA_integer_)
    idx <- which(normalize_rf_feature(names(cor_clusters$membership)) == feature_norm)
    if (!length(idx)) return(NA_integer_)
    unname(cor_clusters$membership[idx[1]])
  }

  module_for <- function(taxon_label) {
    if (is.null(grouped_perm) || !nrow(grouped_perm)) return(NA_character_)
    hit <- grouped_perm %>%
      dplyr::filter(grepl(taxon_label, taxa, fixed = TRUE), n_features > 1) %>%
      dplyr::slice(1)
    if (!nrow(hit)) return(NA_character_)
    hit$taxa[1]
  }

  de_norms <- de$feature_norm
  tbl <- de %>%
    dplyr::left_join(
      rf_ranked %>%
        dplyr::select(
          feature_norm, taxon,
          MeanDecreaseAccuracy, MeanDecreaseGini, perm_importance, perm_importance_sd,
          drop_col_importance, dependence_r2_oob,
          rank_permutation, rank_gini, rank_drop_col,
          perm_top, gini_top, drop_col_tested, drop_col_support
        ),
      by = "feature_norm"
    ) %>%
    dplyr::left_join(maaslin_dir_tbl, by = "feature_norm") %>%
    dplyr::left_join(wilcox_dir_tbl, by = "feature_norm") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      in_rf_model = !is.na(MeanDecreaseAccuracy),
      rf_robustness_score = sum(c(
        perm_top, gini_top,
        if (is.na(drop_col_support)) FALSE else drop_col_support
      ), na.rm = TRUE),
      rf_robustness_label = dplyr::case_when(
        !in_rf_model ~ "Not in RF model",
        rf_robustness_score >= 3 ~ "High (all metrics)",
        rf_robustness_score == 2 ~ "Moderate (2/3 metrics)",
        perm_top ~ "Moderate (permutation only)",
        TRUE ~ "Low"
      ),
      cor_cluster_id = cluster_id_for(feature_norm),
      collinear_partners = {
        p <- partner_lookup(feature_norm)
        if (!length(p)) NA_character_ else paste(vapply(p, species_display_name, character(1)), collapse = "; ")
      },
      collinear_de_partners = {
        p <- partner_lookup(feature_norm, restrict_norms = de_norms)
        if (!length(p)) NA_character_ else paste(vapply(p, species_display_name, character(1)), collapse = "; ")
      },
      collinearity_flag = !is.na(collinear_partners) & collinear_partners != "",
      high_dependence = !is.na(dependence_r2_oob) & dependence_r2_oob >= dependence_threshold,
      correlation_module = if (collinearity_flag) module_for(species) else NA_character_
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      dplyr::desc(n_methods),
      dplyr::desc(rf_robustness_score),
      dplyr::desc(MeanDecreaseAccuracy)
    ) %>%
    dplyr::select(
      species, feature, n_methods, methods, enriched_in, maaslin_coef, maaslin_qval,
      wilcox_median_diff, wilcox_qval,
      MeanDecreaseAccuracy, rank_permutation, perm_importance, perm_importance_sd,
      MeanDecreaseGini, rank_gini,
      drop_col_importance, rank_drop_col, drop_col_tested,
      rf_robustness_score, rf_robustness_label,
      dependence_r2_oob, high_dependence,
      cor_cluster_id, collinearity_flag, collinear_partners, collinear_de_partners,
      correlation_module, in_rf_model
    )

  tbl
}

#' Dot-panel figure: RF metric support for DE-highlighted taxa.
plot_rf_de_support_figure <- function(tbl, min_de_methods = 3, max_taxa = 25,
                                      title = "RF support for differential-abundance taxa") {
  if (!nrow(tbl)) return(NULL)
  show <- tbl %>%
    dplyr::filter(n_methods >= min_de_methods) %>%
    dplyr::slice_head(n = max_taxa)

  if (!nrow(show)) {
    show <- tbl %>% dplyr::slice_head(n = max_taxa)
  }

  met_long <- show %>%
    dplyr::mutate(
      perm_norm = MeanDecreaseAccuracy / max(MeanDecreaseAccuracy, na.rm = TRUE),
      gini_norm = MeanDecreaseGini / max(MeanDecreaseGini, na.rm = TRUE),
      drop_norm = dplyr::if_else(
        is.na(drop_col_importance), NA_real_,
        drop_col_importance / max(drop_col_importance, na.rm = TRUE)
      ),
      collinearity = dplyr::if_else(collinearity_flag, "Collinear module", "Singleton")
    ) %>%
    tidyr::pivot_longer(
      c(perm_norm, gini_norm, drop_norm),
      names_to = "metric",
      values_to = "value_norm"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(
        metric,
        perm_norm = "Permutation",
        gini_norm = "Gini",
        drop_norm = "Drop-column"
      ),
      metric = factor(metric, levels = c("Permutation", "Gini", "Drop-column")),
      species = factor(species, levels = rev(unique(show$species)))
    )

  ggplot2::ggplot(
    met_long,
    ggplot2::aes(x = metric, y = species, size = value_norm, colour = collinearity)
  ) +
    ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::scale_size_area(max_size = 8, na.value = 0) +
    ggplot2::scale_colour_manual(values = c("Collinear module" = "#D55E00", "Singleton" = "#0072B2")) +
    ggplot2::labs(
      title = title,
      subtitle = "Point size = scaled metric within panel; orange = collinear with other RF/DE taxa",
      x = NULL, y = NULL, size = "Relative importance", colour = NULL
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")
}

#' Write RF x DE integration table and optional figure.
export_rf_de_consensus_table <- function(
    tbl,
    out_dir,
    prefix = "rf_de_consensus",
    min_de_methods_strong = 3) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(tbl, file.path(out_dir, paste0(prefix, "_table.csv")), row.names = FALSE)
  strong <- tbl %>% dplyr::filter(n_methods >= min_de_methods_strong)
  write.csv(
    strong,
    file.path(out_dir, paste0(prefix, "_table_strong.csv")),
    row.names = FALSE
  )
  p <- plot_rf_de_support_figure(tbl, min_de_methods = min_de_methods_strong)
  if (!is.null(p)) {
    save_gg_pdf_png(
      p, file.path(out_dir, paste0(prefix, "_support_figure")),
      width = 10, height = max(6, 0.35 * nrow(strong) + 2)
    )
  }
  invisible(list(table = tbl, strong = strong, plot = p))
}

plot_rf_dependence_heatmap <- function(dep_mat, taxon_labels, title = "Feature dependence (OOB R²)") {
  if (is.null(dep_mat) || nrow(dep_mat) < 2) return(NULL)
  labels <- taxon_labels[rownames(dep_mat)]
  labels[is.na(labels)] <- rownames(dep_mat)[is.na(labels)]
  ComplexHeatmap::Heatmap(
    dep_mat,
    name = "OOB R²",
    col = viridisLite::inferno(100),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    row_labels = labels,
    column_labels = labels,
    column_title = title
  )
}

#' Full robust RF importance workflow.
run_rf_deep_analysis <- function(ps, response_var, label, out_dir,
                               n_top_taxa = 100,
                               drop_col_top_n = 30,
                               cor_top_n = 40,
                               ntree_full = 2000,
                               ntree_fast = 500,
                               perm_repeats = 10,
                               cor_threshold = 0.7,
                               n_cores = 1L) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  rf_data <- prepare_rf_matrix(ps, response_var, n_top_taxa = n_top_taxa, add_random_feature = TRUE)
  if (is.null(rf_data)) {
    message("Skipping RF deep analysis for ", label)
    return(NULL)
  }

  rf_fit <- fit_rf_classifier(rf_data, label, response_var, ntree = ntree_full, seed = 42)
  perm_se <- rf_permutation_importance_se(rf_data, n_repeats = perm_repeats, ntree = ntree_fast, seed = 42)

  top_for_drop <- perm_se %>%
    dplyr::filter(!is_random_feature) %>%
    dplyr::slice_head(n = drop_col_top_n) %>%
    dplyr::pull(feature_safe)
  drop_col <- rf_drop_column_importance(
    rf_data, features = top_for_drop, ntree = ntree_fast, seed = 42, n_cores = n_cores
  )

  top_for_cor <- perm_se %>%
    dplyr::filter(!is_random_feature) %>%
    dplyr::slice_head(n = cor_top_n) %>%
    dplyr::pull(feature_safe)
  cor_res <- rf_spearman_correlation(rf_data, features = top_for_cor)
  dep <- rf_feature_dependence(rf_data, features = top_for_cor, ntree = ntree_fast, seed = 42)

  clust <- rf_correlation_clusters(cor_res$matrix, threshold = cor_threshold)
  grouped_perm <- if (length(clust$clusters) > 0) {
    rf_grouped_permutation_importance(rf_data, clust$clusters, ntree = ntree_fast, seed = 42)
  } else {
    NULL
  }

  imp_combined <- rf_fit$importance %>%
    dplyr::select(taxon, MeanDecreaseGini, MeanDecreaseAccuracy, is_random_feature) %>%
    dplyr::left_join(
      perm_se %>% dplyr::select(taxon, perm_importance, perm_importance_sd),
      by = "taxon"
    ) %>%
    dplyr::left_join(
      drop_col %>% dplyr::select(taxon, drop_col_importance),
      by = "taxon"
    ) %>%
    dplyr::left_join(
      dep %>% dplyr::select(taxon, dependence_r2_oob),
      by = "taxon"
    ) %>%
    dplyr::arrange(dplyr::desc(MeanDecreaseAccuracy))

  write.csv(rf_fit$importance, file.path(out_dir, "rf_importance_pain_gi_site.csv"), row.names = FALSE)
  write.csv(perm_se, file.path(out_dir, "rf_permutation_importance.csv"), row.names = FALSE)
  write.csv(drop_col, file.path(out_dir, "rf_drop_column_importance.csv"), row.names = FALSE)
  write.csv(cor_res$long, file.path(out_dir, "rf_collinearity_spearman.csv"), row.names = FALSE)
  write.csv(dep, file.path(out_dir, "rf_feature_dependence.csv"), row.names = FALSE)
  write.csv(imp_combined, file.path(out_dir, "rf_importance_combined.csv"), row.names = FALSE)
  if (!is.null(grouped_perm)) {
    write.csv(grouped_perm, file.path(out_dir, "rf_grouped_permutation_importance.csv"), row.names = FALSE)
  }

  summary_row <- tibble::tibble(
    cohort = label,
    response = response_var,
    n_samples = nrow(rf_data$df),
    n_features = length(rf_data$taxa_safe),
    oob_error = rf_fit$oob_error,
    oob_accuracy = rf_fit$oob_accuracy,
    random_feature_perm_importance = perm_se$perm_importance[perm_se$is_random_feature][1],
    random_feature_rank = which(perm_se$taxon == rf_data$random_feature),
    n_perm_repeats = perm_repeats,
    drop_col_features = length(top_for_drop),
    cor_threshold = cor_threshold,
    n_cor_clusters = length(clust$clusters)
  )
  write.csv(summary_row, file.path(out_dir, "rf_model_summary.csv"), row.names = FALSE)

  p_perm <- plot_rf_importance_perm(rf_fit, top_n = 20, metric = "MeanDecreaseAccuracy")
  save_gg_pdf_png(p_perm, file.path(out_dir, "rf_importance_permutation"), width = 9, height = 8)
  save_gg_pdf_png(p_perm, file.path(out_dir, "rf_importance_pain_gi_site"), width = 9, height = 8)

  p_gini <- plot_rf_importance_perm(rf_fit, top_n = 20, metric = "MeanDecreaseGini")
  save_gg_pdf_png(p_gini, file.path(out_dir, "rf_importance_gini"), width = 9, height = 8)

  top_compare <- imp_combined %>%
    dplyr::filter(!is_random_feature, !is.na(drop_col_importance)) %>%
    dplyr::slice_max(order_by = MeanDecreaseAccuracy, n = min(20, drop_col_top_n), with_ties = FALSE)
  if (nrow(top_compare) > 0) {
    p_cmp <- plot_rf_importance_comparison(top_compare, top_n = nrow(top_compare))
    save_gg_pdf_png(p_cmp, file.path(out_dir, "rf_importance_comparison"), width = 11, height = 8)
  }

  if (!is.null(cor_res$matrix) && nrow(cor_res$matrix) >= 2) {
    ht_cor <- plot_rf_collinearity_heatmap(
      cor_res$matrix, rf_data$name_map,
      title = sprintf("Spearman correlation (top %d RF taxa)", length(top_for_cor))
    )
    save_heatmap_pdf_png(
      ht_cor, file.path(out_dir, "rf_collinearity_heatmap"),
      width = 12, height = 10
    )
  }

  list(
    rf_fit = rf_fit,
    perm_se = perm_se,
    drop_col = drop_col,
    cor = cor_res,
    dependence = dep,
    clusters = clust,
    grouped_perm = grouped_perm,
    importance_combined = imp_combined,
    summary = summary_row
  )
}
