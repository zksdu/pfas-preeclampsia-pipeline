# enrich_ml.R — Step 3/4: 交集靶点富集分析 + 三算法机器学习核心基因筛选
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(glmnet)
  library(randomForest); library(e1071); library(pROC); library(data.table)
})
setwd("E:/Zcode/0911")
dir.create("results/ml", showWarnings = FALSE); dir.create("results/enrich", showWarnings = FALSE)

## ---------- 输入 ----------
merged <- readRDS("results/deg/merged_placenta_expr.rds")   # expr, batch, group
targets <- readLines("results/targets/pfas_ctd_targets.txt")
deg <- fread("results/deg/deg_merged_placenta.csv")
sig <- fread("results/deg/deg_sig_merged_placenta.csv")

## ---------- Step 3a: 交集靶点 ----------
inter <- intersect(targets, sig$gene)
cat("PFAS targets:", length(targets), " sig DEGs:", nrow(sig),
    " intersection:", length(inter), "\n")
writeLines(inter, "results/targets/intersection_targets.txt")

## ---------- Step 3b: GO/KEGG 富集（交集靶点） ----------
ego <- enrichGO(inter, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
fwrite(as.data.frame(ego), "results/enrich/ego_intersection.csv")
ekk <- enrichKEGG(inter, organism = "hsa", keyType = "kegg",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05)
if (!is.null(ekk) && nrow(as.data.frame(ekk))) {
  library(org.Hs.eg.db)
  ekk <- setReadable(ekk, org.Hs.eg.db, keyType = "ENTREZID")
}
fwrite(as.data.frame(ekk), "results/enrich/kegg_intersection.csv")
cat("GO terms:", nrow(as.data.frame(ego)), " KEGG terms:", nrow(as.data.frame(ekk)), "\n")

## ---------- Step 4: 机器学习 ----------
expr <- merged$expr; group <- factor(merged$group, levels = c("CT", "PE"))
feats <- inter
if (length(feats) < 5) stop("intersection too small for ML")
X <- t(expr[feats, , drop = FALSE])   # 样本 x 基因
y <- group
set.seed(42)

# ① 随机森林
rf <- randomForest(x = X, y = y, ntree = 1000, importance = TRUE)
imp <- as.data.frame(importance(rf)[, "MeanDecreaseGini", drop = FALSE])
imp$gene <- rownames(imp)
imp <- imp[order(-imp$MeanDecreaseGini), ]
rf_top <- head(imp$gene, 20)

# ② LASSO（cv.glmnet, lambda.1se）
Xm <- model.matrix(~ . - 1, data = as.data.frame(X))
cv <- cv.glmnet(Xm, y, family = "binomial", alpha = 1, nfolds = 10)
lasso_coef <- coef(cv, s = "lambda.1se")
lasso_top <- rownames(lasso_coef)[which(as.numeric(lasso_coef) != 0)]
lasso_top <- setdiff(lasso_top, "(Intercept)")
if (length(lasso_top) < 3) {
  lasso_coef <- coef(cv, s = "lambda.min")
  lasso_top <- setdiff(rownames(lasso_coef)[which(as.numeric(lasso_coef) != 0)], "(Intercept)")
}
cat("RF top20:", length(rf_top), " LASSO:", length(lasso_top), "\n")

# ③ SVM-RFE（线性核，简单实现）
svm_rfe <- function(X, y) {
  feats <- colnames(X); ret <- c()
  while (length(feats) > 1) {
    m <- svm(x = X[, feats, drop = FALSE], y = y, kernel = "linear", cost = 1)
    w <- sqrt(sum(m$coefs^2)) * abs(t(m$coefs) %*% m$SV)
    rk <- feats[order(as.numeric(abs(w)))]  # 最小权重先剔除
    ret <- c(ret, rk[1])
    feats <- setdiff(feats, rk[1])
  }
  ret  # 剔除顺序；反向即重要性排序
}
rfe_rank <- svm_rfe(as.data.frame(X), y)
svm_top <- rev(rfe_rank)[seq_len(min(20, length(rfe_rank)))]  # 重要性前20

venn_sets <- list(RF = rf_top, LASSO = lasso_top, SVM_RFE = svm_top)
core <- Reduce(intersect, venn_sets)
if (length(core) < 3) {  # 交集不足时降级为 ≥2 算法共享
  tab <- table(unlist(venn_sets))
  core <- names(tab)[tab >= 2]
}
cat("core genes:", paste(core, collapse = ", "), "\n")
writeLines(core, "results/ml/core_genes.txt")
fwrite(list(RF = paste(rf_top, collapse = ";"), LASSO = paste(lasso_top, collapse = ";"),
            SVM = paste(svm_top, collapse = ";"), CORE = paste(core, collapse = ";")),
       "results/ml/feature_sets.csv")

## ---------- 诊断模型 ----------
Xc <- X[, core, drop = FALSE]
df <- data.frame(Xc, PE = as.integer(y == "PE"))
fit_full <- glm(PE ~ ., data = df, family = binomial)
p_hat <- as.numeric(predict(fit_full, type = "response"))
roc <- pROC::roc(df$PE, p_hat, quiet = TRUE)
cat("model AUC:", round(auc(roc), 3), "\n")
capture.output(print(roc), file = "results/ml/roc_summary.txt")
fwrite(data.frame(sample = rownames(X), group = as.character(y), prob = p_hat),
       "results/ml/prediction.csv")
cat("ENRICH_ML DONE\n")
