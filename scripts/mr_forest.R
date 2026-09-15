# mr_forest.R — 核心基因 cis-MR 森林图（优先金标准，退回 v2）
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
setwd("E:/Zcode/0911")

f <- "results/mr/mr_gold_summary.csv"
src <- "native Beta/SE"
if (!file.exists(f)) { f <- "results/mr/mr_full_summary.csv"; src <- "Z+MAF reconstructed" }
d <- fread(f)
d <- d[!is.na(OR)][order(OR)]
d[, ci_lo := exp(beta - 1.96 * se)]
d[, ci_hi := exp(beta + 1.96 * se)]
d[, lab := sprintf("%s (OR=%.3f, P=%.3g, nIV=%d)", gene, OR, P_MR, n_IV)]

png("results/figures/fig_mr_forest.png", width = 2400, height = 1600, res = 220)
par(mar = c(5, 14, 4, 2))
xlim <- range(c(d$ci_lo, d$ci_hi, 1))
plot(NA, NA, xlim = xlim, ylim = c(0.5, nrow(d) + 0.5),
     xlab = "OR for preeclampsia per 1-SD higher expression (95% CI)",
     yaxt = "n", ylab = "", log = "x", main = sprintf("Cis-MR of core genes on PE (%s)", src))
abline(v = 1, lty = 2, col = "grey55")
segments(d$ci_lo, seq_len(nrow(d)), d$ci_hi, seq_len(nrow(d)), lwd = 2)
points(d$OR, seq_len(nrow(d)), pch = 16,
       col = ifelse(d$P_MR < 0.05, "firebrick", "grey25"), cex = 1.2)
axis(2, at = seq_len(nrow(d)), labels = d$lab, las = 1, cex.axis = 0.75, tick = FALSE)
dev.off()
cat("forest plot saved:", nrow(d), "genes; source =", src, "\n")
