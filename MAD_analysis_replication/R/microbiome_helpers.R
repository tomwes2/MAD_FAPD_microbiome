# Shared helpers for MAD taxonomic analysis Rmd pipelines.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Load OTU table, taxonomy, metadata; return phyloseq objects.
load_mad_phyloseq <- function(data_dir = ".") {
  otu_path <- file.path(data_dir, "otu_table.csv")
  tax_path <- file.path(data_dir, "tax_table.csv")
  meta_path <- file.path(data_dir, "subject_data.csv")

  mad_otu_tab <- read.csv(otu_path, row.names = 1, check.names = FALSE)
  mad_tax_tab <- as(read.csv(tax_path, row.names = 1), "matrix")
  mad_samdat_tab <- read.csv(meta_path, row.names = 1, check.names = FALSE)

  common <- intersect(colnames(mad_otu_tab), rownames(mad_samdat_tab))
  mad_otu_tab <- mad_otu_tab[, common, drop = FALSE]
  mad_samdat_tab <- mad_samdat_tab[common, , drop = FALSE]

  mad_samdat_tab$obese <- factor(
    ifelse(mad_samdat_tab$bmi > 30, "obese", "non_obese"),
    levels = c("non_obese", "obese")
  )
  mad_samdat_tab$calprotectin_binary <- factor(
    ifelse(mad_samdat_tab$calprotectin_high %in% c(1, "1", TRUE), "high", "normal"),
    levels = c("normal", "high")
  )
  mad_samdat_tab$pop <- factor(mad_samdat_tab$pop, levels = c("no_pain", "pain"))
  mad_samdat_tab$sample_type <- factor(mad_samdat_tab$sample_type, levels = c("Fecal", "Cecal"))
  mad_samdat_tab$sex <- factor(mad_samdat_tab$sex)
  mad_samdat_tab$phq_depression_severity <- factor(
    mad_samdat_tab$phq_depression_binary,
    levels = c(0, 1, 2),
    labels = c("PHQ minimal (<5)", "PHQ mild (5-9)", "PHQ moderate+ (>=10)")
  )
  mad_samdat_tab$gad_anxiety_severity <- factor(
    mad_samdat_tab$gad_anxiety_binary,
    levels = c(0, 1, 2),
    labels = c("GAD minimal (<5)", "GAD mild (5-9)", "GAD moderate+ (>=10)")
  )
  mad_samdat_tab$phq_difficulty <- factor(
    mad_samdat_tab$phq_difficulty,
    levels = 0:3,
    labels = c(
      "PHQ: Not difficult", "PHQ: Somewhat difficult",
      "PHQ: Very difficult", "PHQ: Extremely difficult"
    )
  )
  mad_samdat_tab$gad_difficulty <- factor(
    mad_samdat_tab$gad_difficulty,
    levels = 0:3,
    labels = c(
      "GAD: Not difficult", "GAD: Somewhat difficult",
      "GAD: Very difficult", "GAD: Extremely difficult"
    )
  )
  mad_samdat_tab$household_id <- ifelse(
    mad_samdat_tab$household == 0,
    paste0("unique_", mad_samdat_tab$ssn),
    paste0("hh_", mad_samdat_tab$household)
  )
  mad_samdat_tab$household_shared <- factor(
    ifelse(mad_samdat_tab$household == 0, "no_shared", "shared"),
    levels = c("no_shared", "shared")
  )

  meta_priority <- intersect(
    c(
      "pop", "calprotectin", "calprotectin_high", "calprotectin_binary",
      "age", "sex", "bmi", "obese",
      "phq_total", "phq_difficulty", "gad_total", "gad_difficulty",
      "household_shared"
    ),
    colnames(mad_samdat_tab)
  )

  mad_ps <- phyloseq(
    otu_table(as.matrix(mad_otu_tab), taxa_are_rows = TRUE),
    tax_table(mad_tax_tab),
    sample_data(mad_samdat_tab)
  )
  fvc_ps <- subset_samples(mad_ps, pop == "pain")
  fps_ps <- subset_samples(mad_ps, sample_type == "Fecal")

  list(
    mad_ps = mad_ps,
    fvc_ps = fvc_ps,
    fps_ps = fps_ps,
    meta_priority = meta_priority
  )
}

prepare_maaslin_metadata <- function(ps, fixed_effects, extra_cols = c("ssn", "household_id")) {
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  rownames(meta) <- sample_names(ps)
  cols <- intersect(c(fixed_effects, extra_cols), names(meta))
  meta <- meta[, cols, drop = FALSE]
  meta <- as.data.frame(lapply(meta, function(x) {
    if (is.data.frame(x)) as.vector(x[[1]]) else as.vector(x)
  }), stringsAsFactors = FALSE)
  cat_vars <- c(
    "sample_type", "pop", "sex", "obese", "household_shared",
    "phq_difficulty", "gad_difficulty",
    "phq_depression_severity", "gad_anxiety_severity",
    "calprotectin_high", "ssn"
  )
  for (nm in intersect(cat_vars, names(meta))) meta[[nm]] <- factor(meta[[nm]])
  keep <- vapply(meta, function(x) length(unique(na.omit(x))) > 1, logical(1))
  meta <- meta[, keep, drop = FALSE]
  rownames(meta) <- sample_names(ps)
  meta
}

run_maaslin <- function(ps, fixed_effects, output_subdir, maaslin_dir,
                        reference = NULL, random_effects = NULL) {
  out <- file.path(maaslin_dir, output_subdir)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dat <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) dat <- t(dat)
  if (ncol(dat) == nsamples(ps) && nrow(dat) == ntaxa(ps)) dat <- t(dat)
  meta <- prepare_maaslin_metadata(ps, fixed_effects)
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
  do.call(Maaslin2::Maaslin2, args)
}

plot_maaslin_volcano <- function(fit, title, out_file, meta_priority) {
  res <- fit$results %>%
    dplyr::filter(metadata %in% meta_priority | metadata %in% c("sample_type", "pop")) %>%
    dplyr::mutate(neglog10q = -log10(pmax(qval, 1e-300)))
  if (!nrow(res)) return(invisible(NULL))
  p <- ggplot2::ggplot(res, ggplot2::aes(x = coef, y = neglog10q, colour = qval < 0.25)) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::facet_wrap(~ metadata, scales = "free_x") +
    ggplot2::scale_colour_manual(values = c("grey60", "firebrick")) +
    ggplot2::labs(title = title, x = "Coefficient", y = "-log10(q-value)") +
    ggplot2::theme_bw()
  ggplot2::ggsave(out_file, p, width = 12, height = max(4, 2 * length(unique(res$metadata))))
}

paired_wilcox_screen <- function(ps, group_var = "sample_type", p_thresh = 0.25, min_prev = 0.1) {
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  mat <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) mat <- t(mat)
  paired_ssn <- meta %>%
    dplyr::group_by(ssn, .data[[group_var]]) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(ssn) %>%
    dplyr::summarise(levels = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(levels == 2) %>%
    dplyr::pull(ssn)
  meta <- meta[meta$ssn %in% paired_ssn, , drop = FALSE]
  mat <- mat[, rownames(meta), drop = FALSE]
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  lv <- levels(factor(meta[[group_var]]))
  purrr::map_dfr(taxa, function(taxon) {
    d <- meta %>%
      dplyr::select(ssn, group = dplyr::all_of(group_var)) %>%
      dplyr::mutate(abund = as.numeric(mat[taxon, rownames(meta)]))
    g1 <- d %>% dplyr::filter(group == lv[1]) %>% dplyr::arrange(ssn)
    g2 <- d %>% dplyr::filter(group == lv[2]) %>% dplyr::arrange(ssn)
    if (nrow(g1) != nrow(g2)) return(NULL)
    tt <- wilcox.test(g1$abund, g2$abund, paired = TRUE)
    tibble::tibble(
      feature = taxon,
      comparison = paste(lv, collapse = "_vs_"),
      pval = tt$p.value,
      median_diff = median(g2$abund - g1$abund, na.rm = TRUE)
    )
  }) %>%
    dplyr::mutate(qval = p.adjust(pval, "BH")) %>%
    dplyr::filter(qval < p_thresh) %>%
    dplyr::arrange(qval)
}

lefse_screen <- function(ps, group_var, p_thresh = 0.25, min_prev = 0.1) {
  mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) mat <- t(mat)
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  grp <- factor(meta[[group_var]])
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  purrr::map_dfr(taxa, function(taxon) {
    x <- as.numeric(mat[taxon, rownames(meta)])
    kw <- kruskal.test(x ~ grp)
    tibble::tibble(feature = taxon, kruskal_p = kw$p.value)
  }) %>%
    dplyr::mutate(comparison = group_var, qval = p.adjust(kruskal_p, "BH")) %>%
    dplyr::filter(qval < p_thresh) %>%
    dplyr::arrange(qval)
}

spearman_screen <- function(ps, meta_var, p_thresh = 0.25, min_prev = 0.1) {
  mat <- as.matrix(otu_table(ps))
  if (!taxa_are_rows(ps)) mat <- t(mat)
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  rownames(meta) <- sample_names(ps)
  smp <- intersect(colnames(mat), rownames(meta))
  x <- as.numeric(meta[smp, meta_var, drop = TRUE])
  prev <- rowMeans(mat > 0)
  taxa <- names(prev[prev >= min_prev])
  purrr::map_dfr(taxa, function(taxon) {
    y <- as.numeric(mat[taxon, smp])
    ct <- cor.test(y, x, method = "spearman", exact = FALSE)
    tibble::tibble(feature = taxon, rho = unname(ct$estimate), pval = ct$p.value)
  }) %>%
    dplyr::mutate(metadata = meta_var, qval = p.adjust(pval, "BH")) %>%
    dplyr::filter(qval < p_thresh) %>%
    dplyr::arrange(qval)
}

normalize_feature <- function(x) {
  gsub("[._ ]+", ".", tolower(trimws(x)))
}

build_consensus <- function(cohort_label, comparison, method_hits, rf_imp = NULL,
                            min_methods = 2, rf_top_n = 30) {
  base <- dplyr::bind_rows(lapply(names(method_hits), function(m) {
    method_hits[[m]] %>%
      dplyr::transmute(feature, feature_norm = normalize_feature(feature), method = m, comparison)
  }))
  if (!nrow(base)) return(tibble::tibble())
  display <- base %>%
    dplyr::group_by(feature_norm) %>%
    dplyr::summarise(feature = feature[which.max(nchar(feature))][1], .groups = "drop")
  tally <- base %>%
    dplyr::group_by(feature_norm, comparison) %>%
    dplyr::summarise(
      n_methods = dplyr::n_distinct(method),
      methods = paste(sort(unique(method)), collapse = "+"),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_methods >= min_methods) %>%
    dplyr::left_join(display, by = "feature_norm") %>%
    dplyr::select(feature, feature_norm, comparison, n_methods, methods) %>%
    dplyr::mutate(cohort = cohort_label) %>%
    dplyr::arrange(dplyr::desc(n_methods), feature)
  if (!is.null(rf_imp) && nrow(rf_imp) > 0) {
    rf_top <- rf_imp %>%
      dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>%
      dplyr::slice_head(n = rf_top_n) %>%
      dplyr::mutate(feature_norm = normalize_feature(taxon))
    tally <- tally %>% dplyr::mutate(rf_top_support = feature_norm %in% rf_top$feature_norm)
  } else {
    tally <- tally %>% dplyr::mutate(rf_top_support = NA)
  }
  tally
}

read_maaslin_hits <- function(maaslin_dir, out_subdir, metadata_terms, q_thresh = 0.25) {
  f <- file.path(maaslin_dir, out_subdir, "significant_results.tsv")
  if (!file.exists(f)) {
    alt <- file.path(maaslin_dir, out_subdir, "all_results.tsv")
    if (!file.exists(alt)) return(tibble::tibble())
    res <- read.delim(alt, check.names = FALSE)
  } else {
    res <- read.delim(f, check.names = FALSE)
  }
  res %>%
    tibble::as_tibble() %>%
    dplyr::filter(metadata %in% metadata_terms, qval < q_thresh) %>%
    dplyr::transmute(feature = feature, metadata, qval, coef, comparison = metadata)
}

method_overlap <- function(hits_list) {
  if (length(hits_list) == 0) {
    return(tibble::tibble(feature = character(), feature_norm = character(),
                          methods = character(), n = integer()))
  }
  base <- dplyr::bind_rows(lapply(names(hits_list), function(m) {
    hits_list[[m]] %>%
      dplyr::transmute(feature, feature_norm = normalize_feature(feature), method = m)
  }))
  if (!nrow(base)) {
    return(tibble::tibble(feature = character(), feature_norm = character(),
                          methods = character(), n = integer()))
  }
  base %>%
    dplyr::group_by(feature_norm) %>%
    dplyr::summarise(
      feature = feature[which.max(nchar(feature))][1],
      methods = paste(sort(unique(method)), collapse = " + "),
      n = dplyr::n_distinct(method),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n))
}

run_permanova_microviz <- function(ps, vars, label) {
  perm <- ps %>%
    microViz::tax_transform("identity", rank = "species") %>%
    microViz::dist_calc("aitchison") %>%
    microViz::dist_permanova(variables = vars, seed = 1234, n_perms = 9999, n_processes = 1)
  list(label = label, table = microViz::perm_get(perm) %>% as.data.frame())
}

run_rf_classifier <- function(ps, response_var, label, n_top_taxa = 100) {
  mat <- as.data.frame(otu_table(ps))
  if (!taxa_are_rows(ps)) mat <- t(mat)
  meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
  y <- meta[[response_var]]
  if (length(unique(na.omit(y))) < 2) {
    message("Skipping RF for ", label, ": <2 classes in ", response_var)
    return(NULL)
  }
  prev <- rowMeans(mat > 0)
  keep_taxa <- names(sort(prev, decreasing = TRUE))[seq_len(min(n_top_taxa, length(prev)))]
  x <- t(mat[keep_taxa, , drop = FALSE])
  safe_names <- make.names(colnames(x), unique = TRUE)
  name_map <- stats::setNames(colnames(x), safe_names)
  colnames(x) <- safe_names
  df <- data.frame(y = factor(y), x, check.names = FALSE)
  df <- df[stats::complete.cases(df), ]
  set.seed(42)
  rf <- randomForest::randomForest(y ~ ., data = df, ntree = 2000, importance = TRUE)
  imp <- as.data.frame(importance(rf))
  imp$taxon <- unname(name_map[rownames(imp)])
  imp$cohort <- label
  imp$response <- response_var
  list(model = rf, importance = imp, oob = rf$err.rate[nrow(rf$err.rate), "OOB"])
}

plot_rf_importance <- function(rf_list, top_n = 20, out_file) {
  if (is.null(rf_list)) return(invisible(NULL))
  imp <- rf_list$importance %>%
    dplyr::arrange(dplyr::desc(MeanDecreaseGini)) %>%
    dplyr::slice_head(n = top_n)
  p <- ggplot2::ggplot(imp, ggplot2::aes(x = reorder(taxon, MeanDecreaseGini), y = MeanDecreaseGini)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("RF importance — ", rf_list$importance$cohort[1]),
      subtitle = paste0("OOB error: ", round(rf_list$oob, 3)),
      x = NULL, y = "Mean decrease Gini"
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(out_file, plot = p, width = 8, height = 7)
  p
}
