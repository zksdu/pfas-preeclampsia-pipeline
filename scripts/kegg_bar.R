# kegg_bar.R — KEGG 富集 bar 图（3 条 nominal 显著通路，补充图）
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
setwd("E:/Zcode/0911")

kk <- fread("results/enrich/kegg_intersection.csv")
kk <- kk[order(pvalue)]
kk[, desc_short := ifelse(nchar(Description) > 42, paste0(substr(Description, 1, 40), "..."), Description)]
kk[, desc_short := factor(desc_short, levels = rev(desc_short))]

png("results/figures/fig_kegg_bar.png", width = 2000, height = 900, res = 200)
par(mar = c(5, 26, 4, 4))
bp <- barplot(-log10(pvalue) ~ desc_short, data = kk, horiz = TRUE,
              col = "#92c5de", border = "#2166ac", xlab = "-log10 (P)",
              xlim = c(0, max(-log10(kk$pvalue)) * 1.25),
              main = "KEGG pathways of 24 intersection targets (nominal P < 0.05)")
text(-log10(kk$pvalue) + 0.15, bp, labels = sprintf("P = %.4f", kk$pvalue), pos = 4, cex = 0.85)
dev.off()
cat("kegg bar saved;", nrow(kk), "terms\n")
