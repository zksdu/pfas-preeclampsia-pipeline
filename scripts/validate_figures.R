# validate_figures.R — 外部验证 + 论文图表生成
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
  library(clusterProfiler); library(org.Hs.eg.db)
})
setwd("E:/Zcode/0911")
dir.create("results/figures", showWarnings = FALSE)

core <- readLines("results/ml/core_genes.txt")
tr <- readRDS("results/deg/merged_placenta_expr.rds")     # 训练：GSE75010
va <- readRDS("results/deg/GSE48424_expr.rds")            # 外部：GSE48424 血液
deg <- fread("results/deg/deg_merged_placenta.csv")
inter <- readLines("results/targets/intersection_targets.txt")

## ---------- 1. 训练集拟合（与 enrich_ml.R 相同设定） ----------
Xtr <- t(tr$expr[core, , drop = FALSE])
ytr <- as.integer(tr$group == "PE")
df_tr <- data.frame(Xtr, PE = ytr)
fit <- glm(PE ~ ., data = df_tr, family = binomial)
p_tr <- as.numeric(predict(fit, type = "response"))
roc_tr <- pROC::roc(ytr, p_tr, quiet = TRUE)

## ---------- 2. GSE48424 外部验证 ----------
common <- intersect(core, rownames(va$expr))
cat("core genes present in GSE48424:", length(common), "/", length(core), "\n")
Xva <- t(va$expr[common, , drop = FALSE])
yva <- as.integer(va$group == "PE")
df_va <- as.data.frame(Xva)
# 用训练集系数构建打分（缺失基因按 0 权重处理）
cf <- coef(fit)
score <- rep(cf["(Intercept)"], nrow(Xva))
for (g in common) if (paste0("`", g, "`") %in% names(cf)) {
  score <- score + cf[paste0("`", g, "`")] * Xva[, g]
} else if (g %in% names(cf)) {
  score <- score + cf[[g]] * Xva[, g]
}
roc_va <- pROC::roc(yva, score, quiet = TRUE)
cat("train AUC:", round(auc(roc_tr), 3), "| external (blood) AUC:", round(auc(roc_va), 3), "\n")
cat("external p-value (Wilcoxon AUC vs 0.5):",
    wilcox.test(score ~ factor(yva))$p.value, "\n")

## ---------- 3. 图 1: 训练 ROC ----------
dfp <- data.frame(FPR = 1 - roc_tr$specificities, TPR = roc_tr$sensitivities,
                  set = sprintf("Training (GSE75010) AUC=%.3f", auc(roc_tr)))
dfp2 <- data.frame(FPR = 1 - roc_va$specificities, TPR = roc_va$sensitivities,
                   set = sprintf("External (GSE48424) AUC=%.3f", auc(roc_va)))
roc_dat <- rbind(dfp, dfp2)
p1 <- ggplot(roc_dat, aes(FPR, TPR, color = set)) +
  geom_line(linewidth = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  labs(x = "False positive rate", y = "True positive rate", color = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = c(0.62, 0.22))
ggsave("results/figures/fig_roc.png", p1, width = 5.5, height = 4.5, dpi = 300)

## ---------- 4. 图 2: 核心基因表达热图（训练集） ----------
em <- tr$expr[core, , drop = FALSE]
z <- t(scale(t(em)))   # 基因方向标准化
ord <- order(tr$group)
z <- z[, ord]
ann <- data.frame(sample = colnames(z), Group = tr$group[ord])
library(pheatmap)
png("results/figures/fig_heatmap_core.png", width = 2200, height = 1600, res = 260)
pheatmap::pheatmap(z, cluster_cols = TRUE, cluster_rows = TRUE,
                   annotation_col = data.frame(row.names = colnames(z), Group = ann$Group),
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   border_color = NA, fontsize_row = 7, fontsize_col = 4,
                   main = "PFAS-PE core genes (GSE75010)")
dev.off()

## ---------- 5. 图 3: GO 富集条图 ----------
ego <- fread("results/enrich/ego_intersection.csv")
if (nrow(ego)) {
  ego2 <- ego[order(p.adjust)][1:min(15, nrow(ego))]
  ego2$Description <- factor(ego2$Description, levels = rev(ego2$Description))
  p3 <- ggplot(ego2, aes(x = Count, y = Description, fill = -log10(p.adjust))) +
    geom_col(width = 0.7) +
    scale_fill_gradient(low = "#67A9CF", high = "#B2182B", name = expression(-log[10](adj.P))) +
    labs(x = "Gene count", y = NULL, title = "GO BP enrichment of PFAS-PE intersection targets") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(size = 11))
  ggsave("results/figures/fig_go_bar.png", p3, width = 9, height = 5.5, dpi = 300)
  cat("GO figure done:", nrow(ego2), "terms\n")
}

## ---------- 6. KEGG 重试 ----------
ekk <- tryCatch({
  k <- enrichKEGG(inter, organism = "hsa", pvalueCutoff = 0.1)
  if (!is.null(k) && nrow(as.data.frame(k))) setReadable(k, org.Hs.eg.db, keyType = "ENTREZID") else k
}, error = function(e) { cat("KEGG still unreachable:", conditionMessage(e), "\n"); NULL })
if (!is.null(ekk)) fwrite(as.data.frame(ekk), "results/enrich/kegg_intersection.csv")
cat("VALIDATION_FIGURES DONE\n")
