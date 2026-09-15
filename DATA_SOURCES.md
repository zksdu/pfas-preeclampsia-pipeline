# Data Sources

All input data are from public repositories. No new human or animal data were generated.

## Gene expression (GEO)

| Accession | Tissue | Samples | Role |
|---|---|---|---|
| GSE75010 | Placenta | 77 normotensive / 80 preeclampsia | Main cohort (DEGs, ML training) |
| GSE73374 | Placenta | 17 / 19 | External validation |
| GSE48424 | Peripheral blood | 18 / 18 | Tissue-specificity control |

https://www.ncbi.nlm.nih.gov/geo/ (downloaded as `*_series_matrix.txt.gz`)

Platform annotations: GPL6244.annot.gz, GPL6480.annot.gz (NCBI GEO); GPL16686 via SOFT family.

## GWAS summary statistics

- **GCST90269903** — Tyrmi JS, et al. JAMA Cardiology 2023 (FinnGen R6 + Estonian Biobank + FINNPEC + InterPregGen meta-analysis; GRCh37; 16,743 cases / 296,824 total).
  https://www.ebi.ac.uk/gwas/studies/GCST90269903 — bulk FTP: `https://ftp.ebi.ac.uk/pub/databases/GWAS/summary_statistics/GCST90269xxx/GCST90269903/`

## eQTL (exposure instruments)

- eQTLGen consortium phase-1 cis-eQTL (Võsa et al., Nat Genet 2021), blood, n = 31,684:
  - significant cis-eQTL: `https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz`
  - full cis-eQTL: `https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz`
  - SNP AF: `2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz`

## Chemistry–gene interactions

- CTD (Comparative Toxicogenomics Database), curated human interactions for PFOA (CID 9554) and PFOS (CID 74483): `https://ctdbase.org/downloads/` (`CTD_chem_gene_ixns.tsv.gz`)

## Protein structures

- AlphaFold DB v6 monomer models, per-UniProt: `https://alphafold.ebi.ac.uk/files/AF-{UNIPROT}-F1-model_v6.pdb`

## Ligands

- PubChem SDF 3-D: PFOS CID 74483 (`https://pubchem.ncbi.nlm.nih.gov/compound/74483`), PFOA CID 9554. Deprotonated forms used throughout (pH 7.4 dominant species).

## Immune signature matrix

- LM22 (Newman et al., Nat Methods 2015), 547 genes × 22 cell types, community mirror: `https://github.com/Moonerss/CIBERSORT` (data via jsDelivr CDN mirror).
