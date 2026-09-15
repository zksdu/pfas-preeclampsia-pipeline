# deg_pipeline.R — Step 2: 三数据集解析、注释映射、合并去批次、limma 差异分析
# 用法: Rscript deg_pipeline.R
options(stringsAsFactors = FALSE, timeout = 600)
suppressPackageStartupMessages({library(limma); library(sva); library(data.table)})
setwd("E:/Zcode/0911")
dir.create("results", showWarnings = FALSE)
dir.create("results/deg", showWarnings = FALSE)

## ---------- 1. series matrix 解析（兼容单行式/两行式元数据） ----------
parse_sm <- function(gz) {
  lines <- readLines(gzfile(gz), warn = FALSE)
  lines <- sub("\r$", "", lines)   # 兼容 CRLF
  beg <- grep("^!series_matrix_table_begin$", lines)
  end <- grep("^!series_matrix_table_end$", lines)
  tab <- lines[(beg[1] + 1):(end[1] - 1)]
  expr <- read.delim(textConnection(tab), sep = "\t", header = TRUE,
                     check.names = FALSE, quote = "")
  rn <- as.character(expr[[1]])
  expr <- as.matrix(expr[, -1])
  colnames(expr) <- gsub('"', "", colnames(expr))
  rownames(expr) <- gsub('"', "", rn)

  n <- NA
  get_vals <- function(key) {
    # 两行式：key 单独一行，下一行为值
    idx <- which(lines == key)
    if (length(idx)) {
      v <- strsplit(gsub('"', "", lines[idx[1] + 1]), "\t", fixed = TRUE)[[1]]
      return(v[nzchar(v)])
    }
    # 单行式：key\tval\tval...
    idx2 <- grep(paste0("^", key, "\\t"), lines)
    if (length(idx2)) {
      f <- strsplit(gsub('"', "", lines[idx2[1]]), "\t", fixed = TRUE)[[1]]
      return(f[-1])
    }
    NULL
  }
  acc <- get_vals("!Sample_geo_accession")
  n <- length(acc)
  ch1_lines <- grep("^!Sample_characteristics_ch1\\t", lines)
  ch1_two <- which(lines == "!Sample_characteristics_ch1")
  ch1 <- c(
    if (length(ch1_two)) list(strsplit(gsub('"', "", lines[ch1_two[1] + 1]), "\t", fixed = TRUE)[[1]]),
    lapply(ch1_lines, function(i) {
      f <- strsplit(gsub('"', "", lines[i]), "\t", fixed = TRUE)[[1]]
      f[-1]
    })
  )
  ch1 <- lapply(ch1, function(x) if (!is.null(x) && length(x) == n) x else rep(NA, n))
  ch1_comb <- if (length(ch1)) apply(do.call(rbind, ch1), 2, paste, collapse = " ;; ") else rep(NA, n)

  gv <- function(key) { v <- get_vals(key); if (is.null(v) || !length(v)) rep(NA, n) else v }
  pheno <- data.frame(
    Sample_geo_accession = acc,
    Sample_title = gv("!Sample_title"),
    Sample_source_name = gv("!Sample_source_name"),
    Sample_characteristics_ch1 = ch1_comb,
    Sample_platform_id = gv("!Sample_platform_id"),
    stringsAsFactors = FALSE
  )
  list(expr = expr, pheno = pheno, platform = unique(pheno$Sample_platform_id))
}

## ---------- 2. 平台注释解析（兼容 GEO annot.gz 与 family.soft 平台表） ----------
parse_gpl <- function(gz) {
  lines <- readLines(gzfile(gz), warn = FALSE)
  lines <- sub("\r$", "", lines)
  # --- 格式 A: GEO annot.gz（列表头行以 ^ 开头且含 TAB） ---
  hdr_i <- grep("^\\^.*\\t", lines)
  if (length(hdr_i)) {
    hdr <- strsplit(sub("^\\^", "", lines[hdr_i[1]]), "\t", fixed = TRUE)[[1]]
    data_lines <- lines[(hdr_i[1] + 1):length(lines)]
    data_lines <- data_lines[nzchar(data_lines)]
    ann <- read.delim(textConnection(paste(data_lines, collapse = "\n")),
                      sep = "\t", header = FALSE, col.names = hdr,
                      check.names = FALSE, quote = "", fill = TRUE,
                      colClasses = "character")
    sym_col <- grep("gene symbol|gene symbol|symbol", hdr, ignore.case = TRUE, value = TRUE)[1]
    id_col <- hdr[1]
    out <- data.frame(probe = ann[[id_col]],
                      symbol = toupper(ann[[sym_col]]),
                      stringsAsFactors = FALSE)
  } else {
    # --- 格式 B: family.soft 平台表 ---
    beg <- grep("^!platform_table_begin$", lines)
    end <- grep("^!platform_table_end$", lines)
    stopifnot(length(beg) >= 1)
    if (!length(end)) end <- length(lines) else end <- end[1] - 1
    tab <- lines[(beg[1] + 1):end]
    ann <- read.delim(textConnection(tab), sep = "\t", header = TRUE,
                      check.names = FALSE, quote = "", fill = TRUE,
                      colClasses = "character")
    hdr <- names(ann)
    cat("GPL columns:", paste(hdr, collapse = " | "), "\n")
    ids <- as.character(ann[[1]])
    sym <- rep(NA_character_, nrow(ann))
    exact <- grep("^gene[ _]symbol$|^symbol$", hdr, ignore.case = TRUE, value = TRUE)
    if (length(exact)) {
      # annot.gz 风格：Gene symbol 列，多映射以 " /// " 分隔取第一个
      v <- ann[[exact[1]]]
      sym <- sub("\\s*///.*$", "", v)
    } else if (any(grepl("gene_assignment", hdr, ignore.case = TRUE))) {
      # Affy family.soft 风格："SYMBOL // name // location"
      ga <- grep("gene_assignment", hdr, ignore.case = TRUE, value = TRUE)[1]
      first <- sub("\\s*//.*$", "", ann[[ga]])
      sym <- ifelse(grepl("^[A-Z][A-Z0-9-]{0,11}$", first) & !grepl(" ", first), first, NA)
    } else {
      st <- grep("symbol|gene name", hdr, ignore.case = TRUE, value = TRUE)
      if (length(st)) sym <- ann[[st[1]]]
    }
    out <- data.frame(probe = ids, symbol = toupper(sym), stringsAsFactors = FALSE)
  }
  out$symbol[!nzchar(out$symbol) | out$symbol == "---" | is.na(out$symbol)] <- NA
  out <- out[!is.na(out$symbol) & !(out$symbol %in% c("NA", "")), ]
  out[!duplicated(out$probe), ]
}

## ---------- 3. 探针 → 基因折叠（按平均表达取最高探针） ----------
collapse_to_gene <- function(expr, ann) {
  m <- merge(data.table(probe = rownames(expr)), ann, by = "probe")
  m <- m[!is.na(symbol)]
  e <- expr[m$probe, , drop = FALSE]
  rownames(e) <- m$symbol
  # 多探针取平均表达最高的探针
  ord <- order(rowMeans(e), decreasing = TRUE)
  e <- e[ord, , drop = FALSE]
  e[!duplicated(rownames(e)), , drop = FALSE]
}

## ---------- 4. 分组判定 ----------
assign_group <- function(pheno) {
  # 在所有 characteristics 列中找疾病信息
  cols <- grep("characteristics|source_name|title", names(pheno), ignore.case = TRUE, value = TRUE)
  g <- rep(NA, nrow(pheno))
  for (cc in cols) {
    v <- tolower(pheno[[cc]])
    is_nonpe <- grepl("non-pe|non pe|nonpe", v)
    pe <- grepl("preeclam|pre-eclam", v) | (grepl("\\bpe\\b", v) & !is_nonpe)
    ct <- (grepl("normotensive|control|healthy", v) | is_nonpe) & !pe
    g[pe] <- "PE"; g[ct] <- "CT"
    if (sum(!is.na(g)) == nrow(pheno)) break
  }
  g
}

## ---------- 5. 逐数据集处理 ----------
gse_list <- list(
  GSE75010 = "data/geo/GSE75010_series_matrix.txt.gz",
  GSE73374 = "data/geo/GSE73374_series_matrix.txt.gz",
  GSE48424 = "data/geo/GSE48424_series_matrix.txt.gz"
)
gpl_list <- c(GSE75010 = "data/GPL6244.annot.gz",
              GSE73374 = "data/GPL16686_head.soft.gz",
              GSE48424 = "data/GPL6480.annot.gz")

parsed <- lapply(names(gse_list), function(g) {
  cat("==== parsing", g, "====\n")
  p <- parse_sm(gse_list[[g]])
  p$pheno$group <- assign_group(p$pheno)
  print(table(p$pheno$group, useNA = "ifany"))
  print(p$pheno[1:min(3, nrow(p$pheno)), grep("source_name|characteristics|group", names(p$pheno), ignore.case = TRUE), drop = FALSE])
  p
})
names(parsed) <- names(gse_list)

mapped <- lapply(names(parsed), function(g) {
  cat("==== annotating", g, "====\n")
  ann <- parse_gpl(gpl_list[[g]])
  cat("annotated probes:", nrow(ann), "\n")
  collapse_to_gene(parsed[[g]]$expr, ann)
})
names(mapped) <- names(parsed)

## ---------- 6. 主队列 GSE75010 limma 差异分析 ----------
e1 <- mapped$GSE75010
grp1 <- parsed$GSE75010$pheno$group[match(colnames(e1), parsed$GSE75010$pheno$Sample_geo_accession)]
keep <- !is.na(grp1)
e1 <- e1[, keep]; grp1 <- factor(grp1[keep], levels = c("CT", "PE"))
cat("GSE75010 primary cohort:", nrow(e1), "genes x", ncol(e1), "samples | groups:\n")
print(table(grp1))
if (max(e1, na.rm = TRUE) > 100) e1 <- log2(e1 + 1)

design <- model.matrix(~ 0 + grp1)
colnames(design) <- levels(grp1)
fit <- lmFit(e1, design)
cont <- makeContrasts(PE - CT, levels = design)
fit2 <- eBayes(contrasts.fit(fit, cont))
deg <- topTable(fit2, number = Inf, adjust = "BH")
fwrite(data.table(gene = rownames(deg), deg), "results/deg/deg_merged_placenta.csv")
sig <- subset(deg, adj.P.Val < 0.05 & abs(logFC) >= 0.58)
cat("significant DEGs (adj.P<0.05, |logFC|>=0.58):", nrow(sig),
    "| up:", sum(sig$logFC > 0), "down:", sum(sig$logFC < 0), "\n")
fwrite(data.table(gene = rownames(sig), sig), "results/deg/deg_sig_merged_placenta.csv")
saveRDS(list(expr = e1, group = grp1), "results/deg/merged_placenta_expr.rds")

## ---------- 6b. GSE73374 尽力映射（GB_ACC -> Entrez -> SYMBOL） ----------
if (requireNamespace("hugene21stprobeset.db", quietly = TRUE)) {
  lines <- readLines(gzfile("data/GPL16686_head.soft.gz"), warn = FALSE)
  lines <- sub("\r$", "", lines)
  pb <- grep("^!platform_table_begin$", lines); pe <- grep("^!platform_table_end$", lines)
  ptab_lines <- sub("\t+$", "", lines[(pb + 1):(pe - 1)])  # 去尾随制表符防空列名
  ptab <- read.delim(textConnection(ptab_lines), sep = "\t", header = TRUE,
                     check.names = FALSE, quote = "", fill = TRUE, colClasses = "character")
  acc <- ptab$GB_ACC
  acc <- acc[nzchar(acc) & acc != "--unknown" & !is.na(acc)]
  cat("GSE73374 GB_ACC available:", length(acc), "\n")
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  amap <- select(org.Hs.eg.db, keys = unique(acc), keytype = "ACCNUM",
                 columns = c("ENTREZID", "SYMBOL"))
  amap <- amap[!is.na(amap$SYMBOL) & nzchar(amap$SYMBOL), ]
  ptab2 <- data.table(probe = ptab$ID, acc = ptab$GB_ACC)
  ann73374 <- merge(ptab2, amap, by.x = "acc", by.y = "ACCNUM", allow.cartesian = TRUE)
  ann73374 <- ann73374[!duplicated(ann73374$probe)]
  ann73374 <- data.frame(probe = ann73374$probe, symbol = ann73374$SYMBOL)
  cat("GSE73374 best-effort annotated probes:", nrow(ann73374), "\n")
  if (nrow(ann73374) > 1000) {
    e2 <- collapse_to_gene(parsed$GSE73374$expr, ann73374)
    g2 <- parsed$GSE73374$pheno$group[match(colnames(e2), parsed$GSE73374$pheno$Sample_geo_accession)]
    k2 <- !is.na(g2)
    e2 <- e2[, k2]; g2 <- factor(g2[k2], levels = c("CT", "PE"))
    if (max(e2, na.rm = TRUE) > 100) e2 <- log2(e2 + 1)
    cat("GSE73374 cohort:", nrow(e2), "genes x", ncol(e2), "samples\n"); print(table(g2))
    design2 <- model.matrix(~ 0 + g2); colnames(design2) <- levels(g2)
    f2 <- eBayes(contrasts.fit(lmFit(e2, design2), makeContrasts(PE - CT, levels = design2)))
    deg2 <- topTable(f2, number = Inf, adjust = "BH")
    fwrite(data.table(gene = rownames(deg2), deg2), "results/deg/deg_gse73374.csv")
    sig2 <- subset(deg2, adj.P.Val < 0.05 & abs(logFC) >= 0.58)
    cat("GSE73374 sig DEGs:", nrow(sig2), "\n")
    fwrite(data.table(gene = rownames(sig2), sig2), "results/deg/deg_sig_gse73374.csv")
  }
}

## ---------- 7. GSE48424（血液）外部验证数据 ----------
g4 <- parsed$GSE48424
cat("GSE48424 tissue hints:\n")
print(table(grepl("blood|plasma|serum", tolower(paste(unlist(g4$pheno), collapse = " ")))))
e4 <- mapped$GSE48424
g4g <- g4$pheno$group[match(colnames(e4), g4$pheno$Sample_geo_accession)]
k4 <- !is.na(g4g)
e4 <- e4[, k4]; g4g <- factor(g4g[k4], levels = c("CT", "PE"))
if (max(e4, na.rm = TRUE) > 100) e4 <- log2(e4 + 1)
cat("GSE48424 cohort:", nrow(e4), "genes x", ncol(e4), "samples\n"); print(table(g4g))
saveRDS(list(expr = e4, group = g4g), "results/deg/GSE48424_expr.rds")
cat("PIPELINE DONE\n")
