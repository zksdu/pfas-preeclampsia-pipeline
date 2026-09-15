# validate_gse73374.R — 独立胎盘队列（GSE73374）外部验证
suppressPackageStartupMessages({library(data.table); library(pROC)})
setwd("E:/Zcode/0911")
core <- readLines("results/ml/core_genes.txt")
tr <- readRDS("results/deg/merged_placenta_expr.rds")

# 重建 GSE73374 表达（与 deg_pipeline 6b 相同流程）
src <- readLines("scripts/deg_pipeline.R")
s1 <- grep("parse_sm <- function", src); e1 <- grep("parse_gpl <- function", src)
eval(parse(text = paste(src[s1:(e1 - 2)], collapse = "\n")))
s1c <- grep("collapse_to_gene <- function", src); e1c <- grep("assign_group <- function", src)
eval(parse(text = paste(src[s1c:(e1c - 2)], collapse = "\n")))
p7 <- parse_sm("data/geo/GSE73374_series_matrix.txt.gz")
p7$pheno$group <- NULL
# 组别
v <- tolower(p7$pheno$Sample_characteristics_ch1)
is_nonpe <- grepl("non-pe|non pe|nonpe", v)
pe <- grepl("preeclam|pre-eclam", v) | (grepl("\\bpe\\b", v) & !is_nonpe)
ct <- (grepl("normotensive|control|healthy", v) | is_nonpe) & !pe
g7 <- rep(NA, nrow(p7$pheno)); g7[pe] <- "PE"; g7[ct] <- "CT"
# 注释
lines <- readLines(gzfile("data/GPL16686_head.soft.gz"), warn = FALSE); lines <- sub("\r$", "", lines)
pb <- grep("^!platform_table_begin$", lines); pe2 <- grep("^!platform_table_end$", lines)
ptab_lines <- sub("\t+$", "", lines[(pb + 1):(pe2 - 1)])
ptab <- read.delim(textConnection(ptab_lines), sep = "\t", header = TRUE,
                   check.names = FALSE, quote = "", fill = TRUE, colClasses = "character")
suppressPackageStartupMessages(library(org.Hs.eg.db))
acc <- ptab$GB_ACC; acc <- acc[nzchar(acc) & acc != "--unknown" & !is.na(acc)]
amap <- select(org.Hs.eg.db, keys = unique(acc), keytype = "ACCNUM", columns = c("ENTREZID", "SYMBOL"))
amap <- amap[!is.na(amap$SYMBOL) & nzchar(amap$SYMBOL), ]
ptab2 <- data.table(probe = ptab$ID, acc = ptab$GB_ACC)
ann <- merge(ptab2, amap, by.x = "acc", by.y = "ACCNUM", allow.cartesian = TRUE)
ann <- ann[!duplicated(ann$probe)]
e7 <- collapse_to_gene(p7$expr, data.frame(probe = ann$probe, symbol = ann$SYMBOL))
g7 <- g7[match(colnames(e7), p7$pheno$Sample_geo_accession)]
k <- !is.na(g7); e7 <- e7[, k]; g7 <- factor(g7[k], levels = c("CT", "PE"))
if (max(e7, na.rm = TRUE) > 100) e7 <- log2(e7 + 1)
cat("GSE73374:", nrow(e7), "genes x", ncol(e7), "samples |"); print(table(g7))

# 拟合训练模型并外推
common <- intersect(core, rownames(e7))
cat("core genes present:", length(common), "/", length(core), "\n")
Xtr <- t(tr$expr[core, , drop = FALSE]); ytr <- as.integer(tr$group == "PE")
fit <- glm(PE ~ ., data = data.frame(Xtr, PE = ytr), family = binomial)
cf <- coef(fit)
Xva <- e7[common, , drop = FALSE]
score <- rep(ifelse("(Intercept)" %in% names(cf), cf["(Intercept)"], 0), ncol(Xva))
for (g in common) {
  nm <- g; nm_bt <- paste0("`", g, "`")
  w <- if (nm %in% names(cf)) cf[[nm]] else if (nm_bt %in% names(cf)) cf[[nm_bt]] else 0
  score <- score + w * as.numeric(Xva[g, ])
}
roc <- pROC::roc(as.integer(g7 == "PE"), score, quiet = TRUE)
cat("EXTERNAL placenta (GSE73374) AUC:", round(auc(roc), 3),
    " p =", wilcox.test(score ~ g7)$p.value, "\n")
fwrite(data.frame(sample = colnames(Xva), group = as.character(g7), score = score),
       "results/ml/external_gse73374_prediction.csv")
png("results/figures/fig_roc_external_placenta.png", width = 1500, height = 1300, res = 260)
plot(roc, col = "#B2182B", main = sprintf("External validation: GSE73374 (placenta)\nAUC = %.3f", auc(roc)))
dev.off()
cat("DONE\n")
