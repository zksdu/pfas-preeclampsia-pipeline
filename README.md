# PFAS × Preeclampsia — Analysis Code Package

**GitHub repository:** https://github.com/zksdu/pfas-preeclampsia-pipeline
**Zenodo archived version:** https://doi.org/10.5281/zenodo.22761408

Code accompanying the manuscript:

> **Unraveling the association between per- and polyfluoroalkyl substances (PFAS) and preeclampsia: a multidimensional strategy integrating network toxicology, machine learning, Mendelian randomization, and molecular docking**

An eight-step closed-loop computational pipeline: computational toxicology → placental transcriptomics → PPI/enrichment → three-algorithm machine learning → two-sample Mendelian randomization with colocalization → immune deconvolution → molecular docking → molecular dynamics.

## Author / Citation

**Bing Song** — The Third Affiliated Hospital of Guangzhou Medical University (bingsong2012683034@gzhmu.edu.cn)

If you use this code, please cite the corresponding manuscript and this archived version (Zenodo DOI: 10.5281/zenodo.22761409).

## Contents

```
scripts/                       All analysis code (see run order below)
DATA_SOURCES.md                Public data sources, accessions and URLs
LICENSE                        MIT license for original code
                               (CIBERSORT* files retain their original academic license)
```

## Environment

- **R** 4.6.1: limma, sva, data.table, clusterProfiler, org.Hs.eg.db, randomForest, glmnet, e1071, pROC, VennDiagram, future, furrr, purrr, dplyr
- **Python** 3.11: openmm 8.6.1 (conda-forge), openmmforcefields 0.15.1, openff-toolkit 0.18.0, openff-interchange 0.5.1, openff-units 0.3.2, pint 0.25.3, rdkit, mdtraj 1.11.1, networkx, lxml, xmltodict, pydantic, numpy, pandas
- **AutoDock Vina** 1.2.7; **Open Babel** 3.1.0 (Python bindings)
- Hardware used: NVIDIA P4 / RTX 4090 GPU (MD), 12-core CPU workstation

> Scripts hard-code a project base directory (`E:/Zcode/0911`); adjust `setwd()`/`BASE` to your local layout before running.

## Run order

| Step | Script(s) | Output |
|---|---|---|
| 1. PFAS target set | `ctd_targets.R` | CTD curated targets (5,063 genes) |
| 2. DEGs | `parse_series_matrix.R`, `deg_pipeline.R` | 65 DEGs (GSE75010); external cohort tables |
| 3. Intersection & enrichment | `enrich_ml.R`, `kegg_retry.R` | 24 intersection targets; GO (98) / KEGG (3) |
| 4. Machine learning | `enrich_ml.R` | 18 core genes; diagnostic model, AUC 0.915 / 0.830 |
| 5. MR + colocalization | `mr_full_analysis.R` (v2), `mr_fullwin_analysis.R` (v3), `mr_forest.R` | ENG OR 0.895, P = 0.037; coloc summaries |
| 6. Immune deconvolution | `immune_infiltration.R`, `immune_figures.R` | CIBERSORT fractions, PE-vs-CT tests, correlations |
| 7. Docking | `receptor_prep_v6.sh`, `receptor_convert.py`, `gen_docking_cmds.py`, `parse_docking.sh`, `docking_heatmap.R` | 36 Vina dockings; heatmap |
| 7b. Pose figures | `pose_figure.py` | 3-D pose renderings |
| 8. Molecular dynamics | `md_build_run.py`, `md_analyze.py` | 3 × 10 ns trajectories + metrics |
| 8b. MD verification | `md_verify_local.py`, `md_verify2.py`, `md_diag_eng.py` | Independent rechecks |
| Figures | `volcano.R`, `fig_md.R`, `validate_figures.R`, `validate_gse73374.R` | Publication figures |

Random seeds: 42 (ML), 123 (CIBERSORT permutations), 123 (MD velocities).

## Notes

- `CIBERSORT*.R` files are the standard CIBERSORT implementation (Newman et al., *Nat Methods* 2015) as distributed by the community wrapper (Moonerss/CIBERSORT), with two small patches documented in-file (base-R quantile normalization; sequential parallel stub). They retain the original academic-use license.
- MD used OpenFF Sage 2.2.1 (SMIRNOFF) with MMFF94-derived ligand charges, hydrogen-mass repartitioning (3 Da) and a 4 fs timestep; see manuscript Methods 2.8.
- MD trajectory files (~13 GB) are not included in this package; available from the corresponding author on reasonable request.
