# immune_figures.R — Step6a 出图：免疫差异箱线图 + 核心基因-免疫细胞相关热图
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table)})
setwd("E:/Zcode/0911")

ftab <- fread("results/immune/cibersort_fractions.csv")
grp <- factor(ftab$group, levels = c("CT", "PE"))
frac <- as.matrix(ftab[, 3:24, with = FALSE]); storage.mode(frac) <- "double"
wt <- fread("results/immune/immune_pe_vs_ct.csv")
corr_dt <- fread("results/immune/immune_core_corr.csv")

## ---- Fig A: 显著差异细胞类型箱线图 ----
sig <- wt[P.adj < 0.05 & !is.na(P)][order(P.adj)]
sig <- sig[median_CT + median_PE > 0]   # 排除全零细胞
cells <- head(sig$cell, 9)
nc <- length(cells); ncol_p <- 3; nrow_p <- ceiling(nc / ncol_p)

png("results/figures/fig_immune_box.png", width = 2600, height = 2100, res = 200)
par(mfrow = c(nrow_p, ncol_p), mar = c(4.5, 4.5, 3.5, 1))
for (cc in cells) {
  v <- frac[, cc]
  boxplot(v ~ grp, names = c("Normal", "PE"), col = c("#4daf4a55", "#e41a1c55"),
          main = sprintf("%s\n(adj.P=%.2g)", cc, wt[cell == cc]$P.adj),
          ylab = "CIBERSOT fraction", cex.main = 1.0, outline = TRUE)
  beeswarm_jit <- jitter(rep(c(1, 2), table(grp)), amount = 0.06)
  points(beeswarm_jit, v, pch = 16, cex = 0.7,
         col = ifelse(grp == "PE", "#e41a1c", "#4daf4a"))
}
dev.off()

## ---- Fig B: 核心基因 x 免疫细胞 rho 热图 ----
core_order <- readLines("results/ml/core_genes.txt")
core_order <- core_order[nzchar(core_order)]
m <- dcast(corr_dt[gene %in% core_order], gene ~ cell, value.var = "rho")
setkey(m, gene); m <- m[core_order]
z <- as.matrix(m[, -1]); rownames(z) <- m$gene
colnames(z) <- gsub(" cells", "", colnames(z))

brk <- seq(-0.6, 0.6, length.out = 61)
cols <- c(colorRampPalette(c("#053061", "#2166ac", "#67a9cf"))(20),
          colorRampPalette(c("#f7f7f7"))(20),
          colorRampPalette(c("#ef8a62", "#b2182b", "#67001f"))(20))
png("results/figures/fig_immune_corr.png", width = 2600, height = 2000, res = 200)
layout(matrix(c(1, 2), nrow = 1), widths = c(4, 0.6))
par(mar = c(13, 13, 4, 0.5))
image(x = seq_len(nrow(z)), y = seq_len(ncol(z)), z = z, axes = FALSE,
      xlab = "", ylab = "", col = cols, breaks = brk,
      main = "Core genes x immune/stromal cells (Spearman rho)")
axis(1, at = seq_len(nrow(z)), labels = rownames(z), las = 2, cex.axis = 0.8, tick = FALSE)
axis(2, at = seq_len(ncol(z)), labels = colnames(z), las = 1, cex.axis = 0.85, tick = FALSE)
abline(h = seq(0.5, ncol(z), 1), v = seq(0.5, nrow(z), 1), col = "grey90", lwd = 0.4)
box()
par(mar = c(13, 1.5, 4, 3))
image(x = 1, y = seq(-0.6, 0.6, length.out = 60), z = matrix(seq(-0.6, 0.6, length.out = 60), nrow = 1),
      col = cols, axes = FALSE, xlab = "", ylab = "")
axis(2, at = seq(-0.6, 0.6, 0.2), cex.axis = 0.8, las = 1)
box()
dev.off()
cat("immune figures saved; significant cells:", nrow(sig), "\n")
