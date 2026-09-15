# immune_infiltration.R — Step6a: CIBERSORT(LM22) on GSE75010 胎盘主队列
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table); library(e1071); library(future); library(furrr); library(purrr); library(dplyr)})
setwd("E:/Zcode/0911")
enableParallel <- function(maxSize = 500, nThreads = NULL) {   # 打桩：顺序执行即可
  options(future.globals.maxSize = maxSize * 1024^2)
  plan(sequential)
}
source("scripts/CIBERSORT_doPerm.R")
source("scripts/CIBERSORT_CoreAlg.R")
source("scripts/CIBERSORT.R")

m <- readRDS("results/deg/merged_placenta_expr.rds")
expr <- as.matrix(m$expr)
storage.mode(expr) <- "double"
cat("mixture:", nrow(expr), "genes x", ncol(expr), "samples; max =", round(max(expr), 2), "\n")

lm22 <- as.matrix(read.delim("data/cibersort/LM22.txt", header = TRUE, sep = "\t",
                             row.names = 1, check.names = FALSE))
storage.mode(lm22) <- "double"
common <- intersect(rownames(lm22), rownames(expr))
cat("LM22 genes matched in mixture:", length(common), "/", nrow(lm22), "\n")
stopifnot(length(common) >= 400)

## ---- CIBERSORT（带缓存：fractions CSV 存在则跳过重算）----
if (!file.exists("results/immune/cibersort_fractions.csv")) {
  set.seed(123)
  res <- cibersort(sig_matrix = lm22, mixture_file = expr, perm = 100, QN = TRUE)
  frac <- res[, 1:22, drop = FALSE]
  pcol <- grep("^P[.-]value$", colnames(res)); ccol <- grep("^Correlation$", colnames(res))
  pval <- as.numeric(res[, pcol]); corr_qc <- as.numeric(res[, ccol])
  fwrite(cbind(sample = rownames(frac), group = as.character(
    factor(m$group, levels = c("CT", "PE"))), frac,
    P.value = pval, Correlation = corr_qc),
    "results/immune/cibersort_fractions.csv")
}
ftab <- fread("results/immune/cibersort_fractions.csv")
grp <- factor(ftab$group, levels = c("CT", "PE"))
frac <- as.matrix(ftab[, 3:24, with = FALSE]); storage.mode(frac) <- "double"
rownames(frac) <- ftab$sample
expr <- expr[, rownames(frac), drop = FALSE]
pval <- ftab$P.value
cat("CIBERSORT done:", ncol(frac), "cell types x", nrow(frac), "samples\n")
cat("samples with significant deconvolution (P<0.05):", sum(pval < 0.05, na.rm = TRUE), "/", length(pval), "\n")

## ---- PE vs CT（Wilcoxon）----
wt <- rbindlist(lapply(colnames(frac), function(cc) {
  w <- suppressWarnings(wilcox.test(frac[grp == "PE", cc], frac[grp == "CT", cc]))
  data.table(cell = cc, median_CT = median(frac[grp == "CT", cc]),
             median_PE = median(frac[grp == "PE", cc]),
             log2FC = log2((median(frac[grp == "PE", cc]) + 1e-6) /
                             (median(frac[grp == "CT", cc]) + 1e-6)),
             P = w$p.value)
}))
wt[, P.adj := p.adjust(P, "BH")]
wt <- wt[order(P)]
fwrite(wt, "results/immune/immune_pe_vs_ct.csv")
print(wt)

## ---- 核心基因 x 细胞比例（Spearman，长格式防错位）----
core <- readLines("results/ml/core_genes.txt")
core_in <- intersect(core, rownames(expr))
corr_dt <- rbindlist(lapply(core_in, function(g) {
  rbindlist(lapply(colnames(frac), function(cc) {
    ct <- suppressWarnings(cor.test(expr[g, ], frac[, cc], method = "spearman"))
    data.table(gene = g, cell = cc, rho = unname(ct$estimate), P = ct$p.value)
  }))
}))
corr_dt[, P.adj := p.adjust(P, "BH"), by = "gene"]
fwrite(corr_dt, "results/immune/immune_core_corr.csv")

focus <- c("ENG", "TREM1", "FLT1")
cat("\n--- Focus genes: cells with |rho|>0.2 & P<0.05 ---\n")
for (g in intersect(focus, core_in)) {
  h <- corr_dt[gene == g & abs(rho) > 0.2 & P < 0.05][order(-abs(rho))]
  if (nrow(h)) {
    cat(g, ":\n")
    cat(do.call(paste0, lapply(seq_len(nrow(h)), function(i)
      sprintf("  %-35s rho=%+.3f P=%.2g\n", h$cell[i], h$rho[i], h$P[i]))))
  }
}
cat("IMMUNE DONE\n")
