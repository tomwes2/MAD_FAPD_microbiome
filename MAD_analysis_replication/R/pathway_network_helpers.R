# Pathway network helpers (no MetaCyc/BioCyc API required).
# HUMAnN already reports MetaCyc pathway names; we link pathways by:
#   - name hierarchy / token overlap (knowledge)
#   - shared functional theme
#   - shared UniRef gene families (from precomputed map or genefamilies scan)
#   - abundance correlation (co-abundance)

`%||%` <- function(x, y) if (is.null(x)) y else x

#' HUMAnN collapsed row name -> MetaCyc-style label (for UniRef map matching).
humann_to_metacyc_label <- function(x) {
  vapply(x, function(s) {
    parts <- strsplit(s, "\\.\\.", fixed = TRUE)[[1]]
    if (length(parts) < 2L) return(s)
    id <- gsub("\\.", "-", parts[1])
    name <- gsub("\\.", " ", parts[2])
    paste0(id, ": ", name)
  }, character(1), USE.NAMES = FALSE)
}

#' Map consensus/MaAsLin feature IDs to rownames in load_pathway_matrix output.
match_pathway_rownames <- function(features, rownames_vec) {
  vapply(features, function(f) {
    if (f %in% rownames_vec) return(f)
    lab <- humann_to_metacyc_label(f)
    if (lab %in% rownames_vec) return(lab)
    idx <- match(make.names(f), make.names(rownames_vec))
    if (!is.na(idx)) return(rownames_vec[idx])
    NA_character_
  }, character(1), USE.NAMES = FALSE)
}

pathway_display_label <- function(x) {
  lab <- humann_to_metacyc_label(x)
  sub("^[^:]+: ", "", lab)
}

pathway_tokens <- function(x, min_len = 4L) {
  lab <- tolower(pathway_display_label(x))
  lab <- gsub("superpathway of ", " ", lab, fixed = TRUE)
  tok <- unique(unlist(strsplit(gsub("[^a-z0-9]+", " ", lab), "\\s+")))
  tok[nchar(tok) >= min_len & !tok %in% c("pathway", "biosynthesis", "degradation", "metabolism", "anaerobic", "aerobic")]
}

token_jaccard <- function(a, b) {
  ta <- pathway_tokens(a)
  tb <- pathway_tokens(b)
  if (!length(ta) || !length(tb)) return(0)
  length(intersect(ta, tb)) / length(union(ta, tb))
}

#' Edges when one pathway name is nested in another (superpathway / sub-pathway).
build_superpathway_edges <- function(features, min_chars = 12L) {
  feats <- unique(features)
  labels <- pathway_display_label(feats)
  edges <- list()
  for (i in seq_along(feats)) {
    for (j in seq_along(feats)) {
      if (i >= j) next
      li <- labels[i]
      lj <- labels[j]
      if (nchar(li) < min_chars && nchar(lj) < min_chars) next
      if (grepl(li, lj, fixed = TRUE) || grepl(lj, li, fixed = TRUE)) {
        edges[[length(edges) + 1L]] <- data.frame(
          from = feats[i], to = feats[j],
          edge_type = "superpathway_name",
          weight = 1,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(edges)) return(data.frame(from = character(), to = character(), edge_type = character(), weight = numeric()))
  do.call(rbind, edges)
}

#' Theme edges only when pathways share a theme AND substantive name/token overlap.
build_theme_edges <- function(features, classify_fn, min_jaccard = 0.12) {
  feats <- unique(features)
  themes <- vapply(feats, classify_fn, character(1))
  edges <- list()
  for (i in seq_len(length(feats) - 1L)) {
    for (j in (i + 1L):length(feats)) {
      if (themes[i] != themes[j]) next
      jac <- token_jaccard(feats[i], feats[j])
      if (jac < min_jaccard) next
      edges[[length(edges) + 1L]] <- data.frame(
        from = feats[i], to = feats[j],
        edge_type = paste0("shared_theme:", themes[i]),
        weight = jac,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(edges)) return(data.frame(from = character(), to = character(), edge_type = character(), weight = numeric()))
  do.call(rbind, edges)
}

#' Token overlap edges (pathway "product/substrate" similarity proxy without API).
build_token_edges <- function(features, min_jaccard = 0.2) {
  feats <- unique(features)
  n <- length(feats)
  edges <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      jac <- token_jaccard(feats[i], feats[j])
      if (jac >= min_jaccard) {
        edges[[length(edges) + 1L]] <- data.frame(
          from = feats[i], to = feats[j],
          edge_type = "shared_tokens",
          weight = jac,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(edges)) return(data.frame(from = character(), to = character(), edge_type = character(), weight = numeric()))
  do.call(rbind, edges)
}

load_uniref_map <- function(csv_path, features = NULL) {
  if (!file.exists(csv_path)) return(NULL)
  m <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!all(c("pathway", "uniref") %in% names(m))) return(NULL)
  if (!is.null(features)) {
    labels <- humann_to_metacyc_label(features)
    m <- m[m$pathway %in% labels, , drop = FALSE]
  }
  split(m$uniref, m$pathway)
}

build_uniref_jaccard_edges <- function(features, uniref_by_pathway_label,
                                      min_jaccard = 0.02, min_shared = 5L) {
  label <- humann_to_metacyc_label(features)
  sets <- lapply(label, function(l) unique(uniref_by_pathway_label[[l]] %||% character(0)))
  names(sets) <- features
  feats <- features[vapply(sets, length, integer(1)) > 0L]
  edges <- list()
  for (i in seq_along(feats)) {
    for (j in seq_along(feats)) {
      if (j <= i) next
      a <- feats[i]
      b <- feats[j]
      sa <- sets[[a]]
      sb <- sets[[b]]
      inter <- length(intersect(sa, sb))
      if (inter < min_shared) next
      jac <- inter / length(union(sa, sb))
      if (jac >= min_jaccard) {
        edges[[length(edges) + 1L]] <- data.frame(
          from = a, to = b,
          edge_type = "shared_uniref",
          weight = jac,
          shared_unirefs = inter,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(edges)) {
    return(data.frame(from = character(), to = character(), edge_type = character(),
                      weight = numeric(), shared_unirefs = integer()))
  }
  do.call(rbind, edges)
}

clr_transform <- function(mat) {
  mat <- mat + 1e-6
  log_mat <- log(mat)
  log_mat[!is.finite(log_mat)] <- NA
  cm <- colMeans(log_mat, na.rm = TRUE)
  cm[!is.finite(cm)] <- 0
  out <- sweep(log_mat, 2, cm, "-")
  out[!is.finite(out)] <- 0
  out
}

build_coabundance_edges <- function(mat, features, method = "spearman",
                                  min_abs_r = 0.45, max_p = 0.05,
                                  edge_type = "coabundance") {
  mapped <- match_pathway_rownames(features, rownames(mat))
  ok <- !is.na(mapped)
  if (sum(ok) < 3L) {
    return(data.frame(from = character(), to = character(), edge_type = character(),
                      weight = numeric(), cor = numeric(), p_value = numeric()))
  }
  id_map <- setNames(features[ok], mapped[ok])
  feats <- mapped[ok]
  if (length(feats) < 3L) {
    return(data.frame(from = character(), to = character(), edge_type = character(),
                      weight = numeric(), p_value = numeric()))
  }
  sub <- mat[feats, , drop = FALSE]
  sub <- clr_transform(sub)
  cm <- cor(t(sub), method = method)
  n <- nrow(cm)
  edges <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      r <- cm[i, j]
      if (is.na(r) || abs(r) < min_abs_r) next
      p <- tryCatch(
        cor.test(sub[i, ], sub[j, ], method = method, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
      if (is.na(p) || p > max_p) next
      edges[[length(edges) + 1L]] <- data.frame(
        from = id_map[rownames(cm)[i]],
        to = id_map[rownames(cm)[j]],
        edge_type = edge_type,
        weight = abs(r),
        cor = r,
        p_value = p,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(edges)) {
    return(data.frame(from = character(), to = character(), edge_type = character(),
                      weight = numeric(), cor = numeric(), p_value = numeric()))
  }
  do.call(rbind, edges)
}

build_paired_delta_edges <- function(mat, pair_map, features,
                                     method = "spearman", min_abs_r = 0.5, max_p = 0.1) {
  mapped <- match_pathway_rownames(features, rownames(mat))
  ok <- !is.na(mapped)
  id_map <- setNames(features[ok], mapped[ok])
  feats <- mapped[ok]
  if (length(feats) < 3L) {
    return(data.frame(from = character(), to = character(), edge_type = character(),
                      weight = numeric(), cor = numeric(), p_value = numeric()))
  }
  deltas <- matrix(NA_real_, nrow = length(feats), ncol = nrow(pair_map),
                   dimnames = list(feats, pair_map$ssn))
  for (k in seq_len(nrow(pair_map))) {
    f <- pair_map$fecal_id[k]
    c <- pair_map$cecal_id[k]
    if (!f %in% colnames(mat) || !c %in% colnames(mat)) next
    deltas[, k] <- mat[feats, c] - mat[feats, f]
  }
  rownames(deltas) <- id_map[rownames(deltas)]
  build_coabundance_edges(deltas, rownames(deltas), method = method,
                          min_abs_r = min_abs_r, max_p = max_p,
                          edge_type = "paired_delta")
}

prune_edges_top_k <- function(edges, k = 6L) {
  if (!nrow(edges) || k <= 0L) return(edges)
  edges <- edges[order(-edges$weight), , drop = FALSE]
  keep <- rep(FALSE, nrow(edges))
  count <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(nrow(edges))) {
    a <- edges$from[i]
    b <- edges$to[i]
    ca <- get0(a, envir = count, ifnotfound = 0L)
    cb <- get0(b, envir = count, ifnotfound = 0L)
    if (ca < k && cb < k) {
      keep[i] <- TRUE
      assign(a, ca + 1L, envir = count)
      assign(b, cb + 1L, envir = count)
    }
  }
  edges[keep, , drop = FALSE]
}

dedupe_edges <- function(edges) {
  if (!nrow(edges)) return(edges)
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  edges$key <- paste(pmin(edges$from, edges$to), pmax(edges$from, edges$to), edges$edge_type, sep = "|")
  edges <- edges[!duplicated(edges$key), , drop = FALSE]
  edges$key <- NULL
  edges
}

combine_edge_layers <- function(...) {
  dfs <- lapply(list(...), function(x) {
    if (is.null(x) || !nrow(x)) return(NULL)
    x
  })
  dfs <- dfs[!vapply(dfs, is.null, logical(1))]
  if (!length(dfs)) {
    return(data.frame(from = character(), to = character(), edge_type = character(), weight = numeric()))
  }
  all_cols <- unique(unlist(lapply(dfs, names)))
  dfs <- lapply(dfs, function(x) {
    for (cn in setdiff(all_cols, names(x))) x[[cn]] <- NA
    x[, all_cols, drop = FALSE]
  })
  dedupe_edges(do.call(rbind, dfs))
}

annotate_nodes <- function(features, node_meta) {
  nm <- node_meta[match(features, node_meta$feature), , drop = FALSE]
  nm$feature <- features
  nm$label <- pathway_display_label(features)
  nm$metacyc_label <- humann_to_metacyc_label(features)
  nm
}

run_communities <- function(edges, nodes, knowledge_only = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    message("Package 'igraph' not installed — using hclust fallback for communities.")
    nodes$community <- fallback_hclust_communities(edges, nodes$feature)
    return(list(graph = NULL, nodes = nodes))
  }
  if (knowledge_only) {
    e <- edges[!grepl("^coabundance$|^paired_delta$", edges$edge_type), , drop = FALSE]
  } else {
    e <- edges
  }
  if (!nrow(e)) {
    nodes$community <- 1L
    return(list(graph = NULL, nodes = nodes))
  }
  g <- igraph::graph_from_data_frame(e[, c("from", "to")], directed = FALSE, vertices = nodes$feature)
  igraph::E(g)$weight <- e$weight
  igraph::E(g)$type <- e$edge_type
  comm <- igraph::cluster_louvain(g, weights = igraph::E(g)$weight)
  nodes$community <- igraph::membership(comm)[nodes$feature]
  list(graph = g, nodes = nodes, membership = comm)
}

fallback_hclust_communities <- function(edges, features, k = 8L) {
  n <- length(features)
  if (n < 2L) return(rep(1L, n))
  adj <- matrix(0, n, n, dimnames = list(features, features))
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      adj[edges$from[i], edges$to[i]] <- edges$weight[i]
      adj[edges$to[i], edges$from[i]] <- edges$weight[i]
    }
  }
  d <- as.dist(1 - adj / max(adj))
  hc <- hclust(d, method = "average")
  k <- min(k, n)
  as.integer(cutree(hc, k = k))
}

#' Default human-readable labels (edit community_label_overrides.csv to change).
default_community_label_overrides <- function() {
  data.frame(
    community = c(1L, 2L, 3L, 4L, 5L, 7L, 8L),
    community_label = c(
      "Arginine, sialic acid & respiration",
      "Fatty acid beta-oxidation, TCA & glyoxylate bypass",
      "Fucose/rhamnose, hexitol, glyoxylate, menaquinol/heme superpathways",
      "Heme and fatty-acid biosynthesis",
      "Fecal-enriched: fermentation, reductive TCA, thiamine/cofactor salvage",
      "Purine degradation & salvage",
      "Pyrimidine/purine nucleotide de novo biosynthesis"
    ),
    stringsAsFactors = FALSE
  )
}

#' Auto-label a module from pathway display names (keyword scoring).
auto_community_label <- function(pathway_labels) {
  text <- tolower(paste(pathway_labels, collapse = " "))
  rules <- list(
    "Mucin sugars, glyoxylate & quinones" = "fucose|rhamnose|hexitol|glyoxyl|menaquinol",
    "Heme & lipid biosynthesis" = "heme|uroporphyrin|palmitate|fatty acid.*biosynthesis",
    "Fatty acid oxidation & TCA" = "beta.oxidation|tca|glyoxylate bypass|glycolysis",
    "Fecal fermentation & anaerobic energy" = "ferment|reductive tca|thiamine|stickland|butanoate",
    "Purine metabolism" = "purine|guanosine|adenosine.*degrad|nucleotide salvage",
    "Nucleotide de novo biosynthesis" = "de novo biosynthesis|deoxyribonucleotide",
    "Amino acid & arginine" = "arginine|ornithine|citrulline|lysine",
    "Mucin-related aminosugars" = "neuramin|acetylglucosamine|chondroitin|sial"
  )
  score_one <- function(pat) {
    m <- gregexpr(pat, text, perl = TRUE)[[1]]
    if (length(m) == 1L && m[1] == -1L) return(0L)
    length(m)
  }
  scores <- vapply(rules, score_one, integer(1))
  if (max(scores) > 0) return(names(which.max(scores))[1])
  "Mixed metabolism"
}

load_community_label_overrides <- function(path, use_defaults = TRUE) {
  if (!is.null(path) && file.exists(path)) {
    o <- read.csv(path, stringsAsFactors = FALSE)
    if (all(c("community", "community_label") %in% names(o))) return(o)
  }
  if (use_defaults) default_community_label_overrides()
  else data.frame(community = integer(), community_label = character(), stringsAsFactors = FALSE)
}

#' Assign community_label to nodes; write labels table.
assign_community_labels <- function(nodes, override_path, use_defaults = TRUE) {
  overrides <- load_community_label_overrides(override_path, use_defaults = use_defaults)

  if (!"community" %in% names(nodes)) stop("nodes missing 'community' column")
  comm_ids <- sort(unique(as.integer(nodes$community)))
  comm_ids <- comm_ids[!is.na(comm_ids)]
  rows <- lapply(comm_ids, function(cid) {
    sub <- nodes[nodes$community == cid, , drop = FALSE]
    n_cecal <- sum(sub$enriched_in == "Cecal", na.rm = TRUE)
    n_fecal <- sum(sub$enriched_in == "Fecal", na.rm = TRUE)
    dom <- if (n_fecal > n_cecal) "Fecal" else "Cecal"
    auto <- auto_community_label(sub$label)
    ov <- overrides$community_label[match(cid, overrides$community)]
    lab <- if (!is.na(ov) && nzchar(ov)) ov else auto
    data.frame(
      community = cid,
      n_pathways = nrow(sub),
      n_cecal = n_cecal,
      n_fecal = n_fecal,
      dominant_site = dom,
      auto_label = auto,
      community_label = lab,
      community_plot_label = paste0("C", cid, ": ", lab),
      stringsAsFactors = FALSE
    )
  })
  labels_df <- do.call(rbind, rows)
  idx <- match(as.integer(nodes$community), labels_df$community)
  nodes$community_label <- labels_df$community_label[idx]
  nodes$community_plot_label <- labels_df$community_plot_label[idx]
  list(nodes = nodes, labels = labels_df)
}

edge_style <- function(edge_types) {
  # Knowledge = dark; co-abundance / paired delta = grey (statistical layers).
  data.frame(
    color = ifelse(
      grepl("^coabundance$", edge_types), "#b3b3b3",
      ifelse(grepl("^paired_delta$", edge_types), "#888888", "#1a1a1a")
    ),
    lty = ifelse(grepl("^coabundance$|^paired_delta$", edge_types), 2L, 1L),
    stringsAsFactors = FALSE
  )
}

community_hull_palette <- function(n) {
  cols <- c(
    "#4393c3", "#d6604d", "#4daf4a", "#984ea3",
    "#ff7f00", "#a65628", "#f781bf", "#999999",
    "#66c2a5", "#fc8d62"
  )
  if (n <= length(cols)) return(cols[seq_len(n)])
  grDevices::colorRampPalette(cols)(n)
}

community_groups_on_graph <- function(g, nodes) {
  comm_ids <- sort(unique(as.integer(nodes$community)))
  comm_ids <- comm_ids[!is.na(comm_ids)]
  lapply(comm_ids, function(cid) {
    which(igraph::V(g)$name %in% nodes$feature[nodes$community == cid])
  })
}

wrap_label <- function(x, width = 26L) {
  words <- strsplit(trimws(x), "\\s+")[[1]]
  if (!length(words)) return(x)
  lines <- character()
  line <- words[1]
  if (length(words) > 1L) {
    for (w in words[-1]) {
      trial <- paste(line, w)
      if (nchar(trial) > width) {
        lines <- c(lines, line)
        line <- w
      } else {
        line <- trial
      }
    }
  }
  c(lines, line)
}

shrink_layout <- function(lay, fraction = 0.72) {
  ctr <- colMeans(lay)
  sweep(sweep(lay, 2, ctr, "-") * fraction, 2, ctr, "+")
}

layout_plot_limits <- function(lay, pad_frac = 0.38) {
  xr <- range(lay[, 1])
  yr <- range(lay[, 2])
  pad_x <- pad_frac * diff(xr)
  pad_y <- pad_frac * diff(yr)
  list(
    xlim = c(xr[1] - pad_x, xr[2] + pad_x),
    ylim = c(yr[1] - pad_y, yr[2] + pad_y)
  )
}

label_community_centroids <- function(lay, g, nodes, labs, hull_cols, xlim, ylim) {
  graphics::par(xpd = NA)
  for (i in seq_len(nrow(labs))) {
    cid <- labs$community[i]
    vidx <- which(igraph::V(g)$name %in% nodes$feature[nodes$community == cid])
    if (length(vidx) < 1L) next
    cx <- mean(lay[vidx, 1])
    cy <- mean(lay[vidx, 2])
    short <- sub("^C[0-9]+: ", "", labs$community_plot_label[i])
    lines <- wrap_label(short, width = 24L)
    if (length(lines) > 3L) {
      lines <- c(lines[1:2], paste0(substr(lines[3], 1, 20), "..."))
    }
    header <- paste0("C", cid)
    cex <- 0.68
    all_lines <- c(header, lines)
    w <- max(strwidth(all_lines, cex = cex, font = c(2, rep(1, length(lines)))))
    h <- sum(strheight(all_lines, cex = cex, font = c(2, rep(1, length(lines))))) * 1.2
    # Skip labels that would sit mostly outside the plot frame
    if (cx - w * 0.6 < xlim[1] || cx + w * 0.6 > xlim[2] ||
        cy - h * 0.6 < ylim[1] || cy + h * 0.6 > ylim[2]) {
      lines <- character()
      all_lines <- header
      w <- strwidth(header, cex = cex, font = 2)
      h <- strheight(header, cex = cex, font = 2) * 1.3
    }
    rect(
      cx - w * 0.58, cy - h * 0.52, cx + w * 0.58, cy + h * 0.52,
      col = grDevices::adjustcolor("white", alpha.f = 0.92),
      border = hull_cols[i], lwd = 2
    )
    y_top <- cy + h * 0.28
    graphics::text(cx, y_top, header, cex = cex, font = 2, col = hull_cols[i])
    if (length(lines)) {
      line_h <- strheight("Ag", cex = cex * 0.88)
      for (j in seq_along(lines)) {
        graphics::text(cx, y_top - j * line_h * 1.15, lines[j], cex = cex * 0.88, font = 1)
      }
    }
  }
  graphics::par(xpd = FALSE)
}

plot_network <- function(g, nodes, out_pdf, width = 18, height = 14,
                         title = "Pathway network",
                         label_communities = TRUE, draw_hulls = TRUE) {
  if (is.null(g) || !requireNamespace("igraph", quietly = TRUE)) return(invisible(NULL))
  pal <- c(Cecal = "#2166ac", Fecal = "#b2182b", Unknown = "#999999")
  vcol <- pal[nodes$enriched_in[match(igraph::V(g)$name, nodes$feature)]]
  vcol[is.na(vcol)] <- pal["Unknown"]
  sizes <- pmax(5, pmin(14, -log10(pmax(nodes$maaslin_q[match(igraph::V(g)$name, nodes$feature)], 1e-10)) * 0.9))

  est <- edge_style(igraph::E(g)$type)

  if (grepl("\\.png$", out_pdf, ignore.case = TRUE)) {
    grDevices::png(out_pdf, width = width, height = height, units = "in", res = 300)
  } else {
    grDevices::pdf(out_pdf, width = width, height = height)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(42)
  lay <- shrink_layout(igraph::layout_with_fr(g), fraction = 0.72)
  lim <- layout_plot_limits(lay, pad_frac = 0.38)

  groups <- NULL
  mark_col <- NULL
  mark_border <- NULL
  labs <- NULL
  if (label_communities && "community" %in% names(nodes)) {
    groups <- community_groups_on_graph(g, nodes)
    n_grp <- length(groups)
    hull_cols <- community_hull_palette(n_grp)
    mark_col <- grDevices::adjustcolor(hull_cols, alpha.f = 0.18)
    mark_border <- grDevices::adjustcolor(hull_cols, alpha.f = 0.75)
    labs <- unique(nodes[, c("community", "community_plot_label", "community_label")])
    labs <- labs[order(labs$community), , drop = FALSE]
  }

  # Extra margins: bottom for module legend, top for title.
  graphics::par(mar = c(11, 3, 4, 3), mgp = c(2, 0.5, 0), xpd = FALSE)
  plot_args <- list(
    g, layout = lay,
    rescale = FALSE,
    xlim = lim$xlim,
    ylim = lim$ylim,
    vertex.color = vcol,
    vertex.size = sizes,
    vertex.label = NA,
    edge.color = est$color,
    edge.lty = est$lty,
    edge.width = 0.4 + 1.8 * igraph::E(g)$weight,
    asp = 0,
    main = title
  )
  if (draw_hulls && !is.null(groups) && length(groups)) {
    plot_args$mark.groups <- groups
    plot_args$mark.col <- mark_col
    plot_args$mark.border <- mark_border
    plot_args$mark.expand <- 8
    plot_args$mark.lwd <- 2
  }
  do.call(igraph::plot.igraph, plot_args)

  if (label_communities && !is.null(labs) && nrow(labs)) {
    label_community_centroids(lay, g, nodes, labs, hull_cols, lim$xlim, lim$ylim)
  }

  if (!is.null(labs) && nrow(labs)) {
    comm_legend <- paste0("C", labs$community, ": ", labs$community_label)
    graphics::legend(
      "bottom", bty = "n", cex = 0.56, ncol = 2, title = "Louvain modules",
      legend = comm_legend,
      fill = grDevices::adjustcolor(hull_cols, alpha.f = 0.35),
      border = mark_border, x.intersp = 0.35, y.intersp = 0.9,
      inset = c(0, -0.02)
    )
  }

  graphics::legend(
    "topleft", bty = "n", cex = 0.72,
    legend = c(
      "Cecal-enriched", "Fecal-enriched",
      "Knowledge (black)", "Co-abundance (grey, dashed)",
      "Paired site shift (dark grey, dashed)"
    ),
    col = c(pal["Cecal"], pal["Fecal"], "#1a1a1a", "#b3b3b3", "#888888"),
    lty = c(NA, NA, 1, 2, 2), pch = c(19, 19, NA, NA, NA),
    lwd = c(NA, NA, 2, 2, 2)
  )
  invisible(NULL)
}

plot_community_bar <- function(nodes) {
  if (!requireNamespace("ggplot2", quietly = TRUE) || all(is.na(nodes$community))) return(NULL)
  df <- as.data.frame(table(
    if ("community_plot_label" %in% names(nodes)) nodes$community_plot_label else nodes$community,
    nodes$enriched_in
  ))
  names(df) <- c("community", "enriched_in", "n")
  df <- df[df$n > 0, , drop = FALSE]
  df$community <- factor(df$community, levels = unique(df$community))
  ggplot2::ggplot(df, ggplot2::aes(community, n, fill = enriched_in)) +
    ggplot2::geom_col(position = "stack") +
    ggplot2::scale_fill_manual(values = c(Cecal = "#2166ac", Fecal = "#b2182b")) +
    ggplot2::labs(x = NULL, y = "Pathways", fill = "Enriched in") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, size = 7.5)
    )
}
