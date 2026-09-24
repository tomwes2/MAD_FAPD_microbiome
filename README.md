# Regional microbiome differences in adolescents with functional abdominal pain

Reproduce the **paired stool vs cecal aspirate** shotgun-metagenomic analyses that support the figures, tables, and supplementary statistics in the integrated manuscript draft.

This repository contains **analysis inputs and R code only**. It does not include manuscript Word/markdown sources or manuscript packaging scripts.

---

## Scientific scope

Within adolescents undergoing colonoscopy for abdominal pain (*n* = 25 paired participants), compare:

1. Species composition (SingleM relative abundances)  
2. MetaCyc pathway abundances (HUMAnN 3)  
3. Site discriminability (random forest + leave-one-subject-out CV)  
4. Pathway co-occurrence modules (integrated network + Louvain)

This is a **within-subject site comparison** (cecum vs stool), not pain vs healthy.

---

## Layout

```text
.
├── README.md
├── run_local_replication.sh          # entrypoint
├── scripts/
│   └── run_bbduk_pe.sh               # paired-end bbduk trim (QC used for this study)
├── MAD_analysis_replication/         # inputs + R pipeline
│   ├── subject_data.csv              # sample metadata
│   ├── sample_ids.csv
│   ├── otu_table.csv / tax_table.csv
│   ├── combined_pathabundance_relab.tsv
│   ├── microbiome_pain_fecal_vs_cecal.Rmd
│   ├── run_*.R
│   └── R/                            # helpers
└── supplementary_data/read_depth/    # library-depth tables (Suppl. depth stats)
```

Generated on your machine (gitignored): `analysis_output/`, `functional_analysis_output/`, `reports/`.

---

## Requirements

R with (among others): `phyloseq`, `microViz`, `ampvis2`, `vegan`, `Maaslin2`, `randomForest`, `igraph`, `ggplot2`, `ggpubr`, `dplyr`, `rmarkdown`.

See `MAD_analysis_replication/README.md` for install notes.

---

## Run

```bash
./run_local_replication.sh
```

Optional:

```bash
export MAD_DATA_DIR=/path/to/MAD_analysis_replication
export MAD_OUTPUT_DIR=/path/to/write/outputs
./run_local_replication.sh
```

---

## Key analysis outputs

| Analysis product | Primary output |
|---|---|
| Alpha/beta diversity | `analysis_output/` (diversity panels / PERMANOVA) |
| Strong-consensus taxa | `analysis_output/consensus/` |
| Pathway site models | `functional_analysis_output/` |
| Louvain modules / pathway network | `functional_analysis_output/pathway_network/` |
| RF importance + LOSO CV | `analysis_output/random_forest/` |
| Suppl. sequencing depth | `supplementary_data/read_depth/*.tsv` |

Manuscript figure assembly and Results packaging are kept outside this repository.

---

## Pipeline order

1. `microbiome_pain_fecal_vs_cecal.Rmd` — diversity, PERMANOVA, MaAsLin2, Wilcoxon, consensus  
2. `run_functional_pathway_analysis.R` — pathway site models  
3. `annotate_direction_outputs.R`  
4. `run_rf_importance_analysis.R` + `run_rf_subject_grouped_cv.R`  
5. `run_pathway_network_analysis.R` — network + Louvain  

---

## Notes

- SingleM tables are relative abundances; libraries are not rarefied.  
- Paired fecal and cecal sequencing depths do not differ meaningfully (see `supplementary_data/read_depth/`).  
- Louvain partitioning can vary slightly across runs; manuscript figures used the locked C1–C9 labeling from the primary analysis.  
- Column `ssn` in `subject_data.csv` is the study subject code (e.g. `MAD001`), used as the paired-subject ID.
- Upstream QC for this study was `bbduk` trimming only (`scripts/run_bbduk_pe.sh`). Assembly/mapping pipelines were not used for the site comparison.
