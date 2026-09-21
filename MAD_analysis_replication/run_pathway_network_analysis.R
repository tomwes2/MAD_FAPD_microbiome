#!/usr/bin/env Rscript
# Pathway network analysis: fecal vs cecal (pain cohort).
# No MetaCyc/BioCyc API — uses HUMAnN pathway names, themes, optional UniRef map,
# and sample-level co-abundance on the pain subset.
#
# Prerequisite: run_functional_pathway_analysis.R (consensus + pathway matrix).
# Optional: functional_analysis_output/gprofiler_pain_site/pathway_to_uniref_map.csv
#   (built earlier from combined_genefamilies.tsv; speeds UniRef edges).
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
  library(ggplot2)
})
source(file.path(data_dir, "R", "manuscript_utils.R"))
if (file.exists(file.path(data_dir, "R", "pathway_functional_categories.R"))) {
  source(file.path(data_dir, "R", "pathway_functional_categories.R"))
}
source(file.path(data_dir, "R", "pathway_network_helpers.R"))

out_dir <- file.path(output_root, "functional_analysis_output", "pathway_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

consensus_path <- file.path(output_root, "functional_analysis_output", "consensus_pain_gi_site_pathways.csv")
if (!file.exists(consensus_path)) {
  stop("Missing ", consensus_path, " — run run_functional_pathway_analysis.R first.")
}

nodes_raw <- read.csv(consensus_path, stringsAsFactors = FALSE)
# Manuscript network uses MaAsLin unpaired + paired consensus (n_methods >= 3 → 93 pathways)
if ("n_methods" %in% names(nodes_raw)) {
  before <- nrow(nodes_raw)
  nodes_raw <- nodes_raw[nodes_raw$n_methods >= 3, , drop = FALSE]
  message("Filtered consensus n_methods >= 3: ", nrow(nodes_raw), " / ", before, " pathways")
}
features <- unique(nodes_raw$feature)
message("Network nodes: ", length(features), " consensus pathways")

# --- Load abundances (pain cohort) ---
id_lookup <- load_sample_id_map(file.path(data_dir, "sample_ids.csv"))
path_tsv <- file.path(data_dir, "combined_pathabundance_relab.tsv")
if (!file.exists(path_tsv)) stop("Missing ", path_tsv)

message("Loading HUMAnN pathway matrix...")
path_mat <- load_pathway_matrix(path_tsv, id_lookup)

meta <- read.csv(file.path(data_dir, "subject_data.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$sampleID
meta <- prepare_meta(meta)

pain_ids <- rownames(meta)[meta$pop == "pain" & rownames(meta) %in% colnames(path_mat)]
path_pain <- path_mat[, pain_ids, drop = FALSE]

# Paired fecal + cecal
pain_long <- meta[meta$pop == "pain", c("ssn", "sampleID", "sample_type"), drop = FALSE]
fecal <- pain_long[pain_long$sample_type == "Fecal", c("ssn", "sampleID")]
cecal <- pain_long[pain_long$sample_type == "Cecal", c("ssn", "sampleID")]
pair_map <- merge(fecal, cecal, by = "ssn", suffixes = c("_fecal", "_cecal"))
names(pair_map) <- c("ssn", "fecal_id", "cecal_id")

pair_map <- pair_map[
  pair_map$fecal_id %in% colnames(path_pain) & pair_map$cecal_id %in% colnames(path_pain),
  ,
  drop = FALSE
]
message("Pain samples with HUMAnN: ", ncol(path_pain), "; paired subjects: ", nrow(pair_map))

# --- Knowledge edges (no API) ---
message("Building knowledge edges...")
e_super <- build_superpathway_edges(features)
e_theme <- build_theme_edges(features, classify_pathway_category)
e_token <- build_token_edges(features, min_jaccard = 0.22)
knowledge_edges <- combine_edge_layers(e_super, e_theme, e_token)

uniref_map_path <- file.path(output_root, "functional_analysis_output", "gprofiler_pain_site", "pathway_to_uniref_map.csv")
e_uniref <- NULL
if (file.exists(uniref_map_path)) {
  message("UniRef edges from ", uniref_map_path)
  ur_list <- load_uniref_map(uniref_map_path, features = features)
  e_uniref <- build_uniref_jaccard_edges(features, ur_list, min_jaccard = 0.03, min_shared = 15L)
  knowledge_edges <- combine_edge_layers(knowledge_edges, e_uniref)
} else {
  message("No pathway_to_uniref_map.csv — skipping UniRef layer. ",
          "Re-run run_gprofiler_pain_pathways.R once (or we can build a small map).")
}

# --- Co-abundance edges ---
message("Building co-abundance edges (CLR + Spearman)...")
e_co <- build_coabundance_edges(
  path_pain, features,
  min_abs_r = 0.55, max_p = 0.05, edge_type = "coabundance"
)
e_co <- prune_edges_top_k(e_co, k = 5L)
e_delta <- build_paired_delta_edges(
  path_pain, pair_map, features,
  min_abs_r = 0.55, max_p = 0.1
)
e_delta <- prune_edges_top_k(e_delta, k = 5L)

all_edges <- combine_edge_layers(knowledge_edges, e_co, e_delta)

# Edges with both biological and statistical support
if (nrow(knowledge_edges) && nrow(combine_edge_layers(e_co, e_delta))) {
  kn <- paste(pmin(knowledge_edges$from, knowledge_edges$to),
              pmax(knowledge_edges$from, knowledge_edges$to), sep = "|")
  co <- combine_edge_layers(e_co, e_delta)
  cn <- paste(pmin(co$from, co$to), pmax(co$from, co$to), sep = "|")
  supported <- all_edges[paste(pmin(all_edges$from, all_edges$to),
                               pmax(all_edges$from, all_edges$to), sep = "|") %in% intersect(kn, cn), , drop = FALSE]
  write.csv(supported, file.path(out_dir, "edges_supported.csv"), row.names = FALSE)
}

# --- Nodes + communities ---
nodes <- annotate_nodes(features, nodes_raw)
nodes$functional_theme <- vapply(features, classify_pathway_category, character(1))
res_k <- run_communities(knowledge_edges, nodes, knowledge_only = TRUE)
nodes_k <- res_k$nodes

res_all <- run_communities(all_edges, nodes_k, knowledge_only = FALSE)
nodes_out <- res_all$nodes
nodes_out$community_knowledge <- nodes_k$community[match(nodes_out$feature, nodes_k$feature)]

kit_override <- file.path(data_dir, "functional_analysis_output", "pathway_network", "community_label_overrides.csv")
override_csv <- file.path(out_dir, "community_label_overrides.csv")
if (!file.exists(override_csv) && file.exists(kit_override)) {
  file.copy(kit_override, override_csv, overwrite = FALSE)
}
if (!file.exists(override_csv)) {
  write.csv(default_community_label_overrides(), override_csv, row.names = FALSE)
  message("Wrote default labels to ", override_csv, " (edit and re-run to customize)")
}
lab <- assign_community_labels(nodes_out, override_path = override_csv)
nodes_out <- lab$nodes
write.csv(lab$labels, file.path(out_dir, "community_labels.csv"), row.names = FALSE)

# --- Write tables (Cytoscape-ready) ---
write.csv(nodes_out, file.path(out_dir, "nodes.csv"), row.names = FALSE)
write.csv(all_edges, file.path(out_dir, "edges_all.csv"), row.names = FALSE)
write.csv(knowledge_edges, file.path(out_dir, "edges_knowledge.csv"), row.names = FALSE)
if (nrow(e_co)) write.csv(e_co, file.path(out_dir, "edges_coabundance.csv"), row.names = FALSE)
if (nrow(e_delta)) write.csv(e_delta, file.path(out_dir, "edges_paired_delta.csv"), row.names = FALSE)

mod_tab <- nodes_out %>%
  group_by(community, community_label, enriched_in) %>%
  summarise(n = dplyr::n(), pathways = paste(label, collapse = "; "), .groups = "drop") %>%
  arrange(community, desc(n))
write.csv(mod_tab, file.path(out_dir, "community_summary.csv"), row.names = FALSE)

# Edge counts
edge_counts <- as.data.frame(table(all_edges$edge_type))
names(edge_counts) <- c("edge_type", "n")
write.csv(edge_counts, file.path(out_dir, "edge_type_counts.csv"), row.names = FALSE)

# --- Plots ---
if (!is.null(res_all$graph)) {
  plot_network(
    res_all$graph, nodes_out, file.path(out_dir, "network_full.pdf"),
    title = "Pathway network: knowledge (black) + co-abundance (grey)",
    label_communities = TRUE, draw_hulls = TRUE
  )
  plot_network(
    res_all$graph, nodes_out, file.path(out_dir, "network_full.png"),
    title = "Pathway network: knowledge (black) + co-abundance (grey)",
    label_communities = TRUE, draw_hulls = TRUE
  )
}
if (!is.null(res_k$graph)) {
  plot_network(
    res_k$graph, nodes_k, file.path(out_dir, "network_knowledge_only.pdf"),
    title = "Knowledge-only pathway network", label_communities = FALSE
  )
}

p_bar <- plot_community_bar(nodes_out)
if (!is.null(p_bar)) {
  ggsave(file.path(out_dir, "community_by_enrichment.pdf"), p_bar, width = 9, height = 5)
  ggsave(file.path(out_dir, "community_by_enrichment.png"), p_bar, width = 9, height = 5, dpi = 300)
}

# --- README ---
writeLines(c(
  "# Pathway network (pain: fecal vs cecal)",
  "",
  "Built without MetaCyc/BioCyc API. HUMAnN already used MetaCyc to define pathways.",
  "",
  "## Edge types",
  "- **superpathway_name**: one pathway name nested in another",
  "- **shared_theme:***: same broad category (amino acid, carbohydrate, etc.)",
  "- **shared_tokens**: word overlap on pathway titles (substrate/product proxy)",
  "- **shared_uniref**: overlap of HUMAnN gene families (if pathway_to_uniref_map.csv exists)",
  "- **coabundance**: Spearman |r| on CLR abundances across pain samples",
  "- **paired_delta**: Spearman on (cecal - fecal) per subject",
  "",
  "## Files",
  "- `nodes.csv` — import to Cytoscape; color by `enriched_in`, size by `maaslin_q`",
  "- `edges_all.csv` — full network",
  "- `edges_knowledge.csv` — biochemistry/theme layer only",
  "- `community_summary.csv` — Louvain modules",
  "",
  "## Re-run",
  "```bash",
  "Rscript run_pathway_network_analysis.R",
  "```",
  "",
  paste0("Nodes: ", nrow(nodes_out), "; edges (all): ", nrow(all_edges),
         "; knowledge: ", nrow(knowledge_edges), "; paired: ", nrow(pair_map))
), file.path(out_dir, "README.md"))

message("\nDone. Outputs in ", out_dir)
message("Edges — knowledge: ", nrow(knowledge_edges),
        " | coabundance: ", nrow(e_co),
        " | paired delta: ", nrow(e_delta),
        " | total: ", nrow(all_edges))
