# flowchart.R — Fig.1 技术路线图（base graphics 版）
options(stringsAsFactors = FALSE)
setwd("E:/Zcode/0911")

png("results/figures/fig_overview.png", width = 2200, height = 2500, res = 200)
par(mar = c(1, 1, 3, 1))
plot(NA, NA, xlim = c(0, 10), ylim = c(0, 12), axes = FALSE, xlab = "", ylab = "",
     main = "Study design overview")

box_col <- c("#e8f0fe", "#e8f0fe", "#e8f0fe", "#e8f0fe", "#e8f0fe", "#e8f0fe", "#e8f0fe", "#e8f0fe")
tar_col <- "#fde9d9"

steps <- list(
  list(y = 11.2, main = "1  Computational toxicology & target profiling",
       sub = "PFOA / PFOS  |  CTD curated human targets (merged n = 5,063)"),
  list(y = 9.8, main = "2  Placental DEGs & intersection",
       sub = "GSE75010 limma (65 DEGs) intersect PFAS targets -> 24 genes"),
  list(y = 8.4, main = "3  PPI network & functional enrichment",
       sub = "STRING/MCODE  |  GO BP 98 terms, KEGG 3 terms"),
  list(y = 7.0, main = "4  Machine learning feature selection",
       sub = "RF + LASSO + SVM-RFE consensus -> 18 core genes; AUC 0.915 / 0.830 external"),
  list(y = 5.6, main = "5  Causal inference (MR + colocalization)",
       sub = "eQTLGen cis-eQTL -> PE GWAS (n = 296,824); ENG OR 0.895, P = 0.037"),
  list(y = 4.2, main = "6  Immune microenvironment",
       sub = "CIBERSORT/LM22: neutrophils/M2 down, CD8 T up; gene-immune correlations"),
  list(y = 2.8, main = "7  Molecular docking & MD",
       sub = "Vina 18 genes x 2 ligands (top -9.53); OpenMM Sage 10 ns"),
  list(y = 1.4, main = "8  Experimental validation (planned)",
       sub = "HTR-8/SVneo + HUVEC: PFOA/PFOS exposure, siRNA rescue")
)

for (s in steps) {
  y <- s$y
  colr <- if (grepl("^1", s$main)) tar_col else box_col[1]
  rect(1.2, y - 0.45, 8.8, y + 0.45, col = colr, border = "#3c78d8", lwd = 2)
  text(1.5, y + 0.13, s$main, adj = c(0, 0.5), cex = 1.05, font = 2)
  text(1.5, y - 0.22, s$sub, adj = c(0, 0.5), cex = 0.85, col = "grey25")
}
for (i in 1:(length(steps) - 1)) {
  y0 <- steps[[i]]$y - 0.45; y1 <- steps[[i + 1]]$y + 0.45
  arrows(5, y0, 5, y1, length = 0.12, lwd = 2, col = "#3c78d8")
}
dev.off()
cat("flowchart saved\n")
