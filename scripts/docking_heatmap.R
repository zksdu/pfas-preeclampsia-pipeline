# docking_heatmap.R — 18 基因 x 2 配体 Vina 打分热图（越负结合越强）
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({library(data.table)})
setwd("E:/Zcode/0911")

d <- fread("results/docking/docking_summary.csv")
d[, gene := factor(gene, levels = unique(d[order(affinity_kcal_mol)]$gene))]
m <- dcast(d, gene ~ ligand, value.var = "affinity_kcal_mol")

png("results/figures/fig_docking_heatmap.png", width = 1700, height = 2200, res = 220)
op <- par(mar = c(5.5, 12, 4, 3))
z <- as.matrix(m[, .(PFOA, PFOS)])
brk <- seq(-10, -4, length.out = 61)
cols <- colorRampPalette(c("#08306b", "#2171b5", "#6baed6", "#f7fbff"))(60)
image(x = 1:2, y = seq_len(nrow(m)), z = t(z), axes = FALSE,
      xlab = "", ylab = "", col = cols, zlim = c(-10, -4), main = "Vina docking affinity (kcal/mol)")
axis(1, at = 1:2, labels = c("PFOA", "PFOS"), cex.axis = 1.3)
axis(2, at = seq_len(nrow(m)), labels = m$gene, las = 1, cex.axis = 0.9, tick = FALSE)
for (i in seq_len(nrow(m))) for (j in 1:2) {
  text(j, i, sprintf("%.2f", z[i, j]), cex = 0.85,
       col = if (z[i, j] < -8) "white" else "black")
}
box()
abline(h = seq(0.5, nrow(m) + 0.5, 1), col = "grey88", lwd = 0.6)
dev.off()
cat("heatmap saved;", nrow(m), "genes\n")
