# Paper 1 replication kit

Inputs and R scripts for the **paired fecal vs cecal** analysis (adolescents with abdominal pain).

## Inputs

| File | Role |
|------|------|
| `subject_data.csv` | Sample metadata (85 rows; study IDs, site, covariates) |
| `sample_ids.csv` | HUMAnN column name → MAD `sampleID` |
| `otu_table.csv` | SingleM relative-abundance profiles |
| `tax_table.csv` | Taxonomy for OTU rows |
| `combined_pathabundance_relab.tsv` | HUMAnN MetaCyc pathway relative abundances |

## Scripts (run via `../run_local_replication.sh` or `Rscript run_replication.R`)

| Script | Role |
|--------|------|
| `microbiome_pain_fecal_vs_cecal.Rmd` | Taxonomic site comparison |
| `run_functional_pathway_analysis.R` | Pathway site comparison |
| `annotate_direction_outputs.R` | Direction labels on saved tables |
| `run_rf_importance_analysis.R` | Random forest importance |
| `run_rf_subject_grouped_cv.R` | Leave-one-subject-out CV + permutation |
| `run_pathway_network_analysis.R` | Pathway network + Louvain modules |
| `R/` | Shared helpers (`paths.R`, phyloseq loaders, network, RF) |

## Outputs (written under `MAD_OUTPUT_DIR`)

- `analysis_output/`
- `functional_analysis_output/` (includes `pathway_network/`)
- `reports/` (knitted HTML)

## R packages

```r
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggpubr", "vegan", "randomForest",
  "igraph", "rmarkdown", "knitr"
))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "Maaslin2"))
install.packages(c("microViz", "ampvis2", "pheatmap"))
```
