# kegg_retry.R — KEGG 富集重试（rest.kegg.jp 已恢复）
options(stringsAsFactors = FALSE, timeout = 600)
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(data.table)})
setwd("E:/Zcode/0911")

inter <- readLines("results/targets/intersection_targets.txt")
cat("intersection genes:", length(inter), "\n")

sym2eg <- bitr(inter, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("mapped to Entrez:", nrow(sym2eg), "/", length(inter), "\n")

ekk <- enrichKEGG(gene = sym2eg$ENTREZID, organism = "hsa",
                  pAdjustMethod = "BH", pvalueCutoff = 1)   # 放宽拿全表，标注显著性
if (!is.null(ekk) && nrow(as.data.frame(ekk))) ekk <- setReadable(ekk, org.Hs.eg.db, "ENTREZID")
kk_df <- as.data.frame(ekk)
if (nrow(kk_df)) kk_df <- kk_df[order(kk_df$pvalue), ]
fwrite(kk_df, "results/enrich/kegg_intersection.csv")
cat("KEGG terms (all):", nrow(kk_df),
    " significant(p.adjust<0.05):", sum(kk_df$p.adjust < 0.05), "\n")
if (nrow(kk_df)) print(head(kk_df[, c("ID", "Description", "GeneRatio", "pvalue", "p.adjust")], 15))
cat("KEGG DONE\n")
