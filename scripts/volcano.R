# volcano.R — GSE75010 主队列 DEG 火山图（标注 24 交集靶点）
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
setwd("E:/Zcode/0911")

deg <- fread("results/deg/deg_merged_placenta.csv")
inter <- readLines("results/targets/intersection_targets.txt")
inter <- inter[nzchar(inter)]

deg[, sig := ifelse(adj.P.Val < 0.05 & logFC >= 0.58, "up",
             ifelse(adj.P.Val < 0.05 & logFC <= -0.58, "down", "ns"))]
deg[, is_target := gene %in% inter]
deg[, padj_log := -log10(pmax(adj.P.Val, 1e-300))]

cols <- c(up = "#d6604d", down = "#4393c3", ns = "grey80")
png("results/figures/fig_volcano.png", width = 2000, height = 1700, res = 200)
par(mar = c(5, 5, 4, 2))
plot(deg$logFC, deg$padj_log, pch = 16, cex = 0.55,
     col = cols[deg$sig], xlab = "log2 fold change (PE vs normal)", ylab = "-log10 (adjusted P)",
     main = sprintf("Placental DEGs, GSE75010 (n up = %d, n down = %d)",
                    sum(deg$sig == "up"), sum(deg$sig == "down")))
abline(v = c(-0.58, 0.58), lty = 2, col = "grey55")
abline(h = -log10(0.05), lty = 2, col = "grey55")
lab <- deg[is_target == TRUE]
lab <- lab[order(-padj_log)][1:12]
lab <- lab[order(lab$logFC)]
points(lab$logFC, lab$padj_log, pch = 21, cex = 0.9, col = "black", bg = "#f4a582")
placed <- data.frame(x = numeric(0), y = numeric(0))
for (i in seq_len(nrow(lab))) {
  tx <- lab$logFC[i]; ty <- lab$padj_log[i]
  pos <- 4; dy <- 0
  for (j in seq_len(nrow(placed))) {
    if (abs(placed$x[j] - tx) < 0.14 & abs(placed$y[j] - ty) < 0.9) {
      dy <- dy - 0.55
    }
  }
  text(tx + 0.03, ty + dy, labels = lab$gene[i], pos = pos, cex = 0.75, offset = 0.35)
  placed <- rbind(placed, data.frame(x = tx, y = ty + dy))
}
legend("topright", c("up (PE)", "down (PE)", "ns", "PFAS-PE intersection target"),
       col = c(cols["up"], cols["down"], cols["ns"], "black"),
       pch = c(16, 16, 16, 21), pt.bg = c(NA, NA, NA, "#f4a582"), bty = "n", cex = 0.9)
dev.off()
cat("volcano saved; sig up:", sum(deg$sig == "up"), " down:", sum(deg$sig == "down"),
    " targets labelled:", nrow(lab), "\n")
