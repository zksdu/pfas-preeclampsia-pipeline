# mr_analysis.R — Step 5: 双样本 MR（核心基因 cis-eQTL → PE GWAS）+ 共定位
# eQTLGen 显著 cis-eQTL 全表（含 SNPChr/SNPPos）× GWAS Catalog GCST90269903（chr:pos，无 rsID）
options(stringsAsFactors = FALSE, timeout = 1800)
suppressPackageStartupMessages({library(data.table); library(coloc)})
setwd("E:/Zcode/0911")
dir.create("results/mr", showWarnings = FALSE)

core <- readLines("results/ml/core_genes.txt")
cat("core genes:", paste(core, collapse = ", "), "\n")

## ---------- 1. 读入 PE GWAS（8列，chr/pos 无 rsID） ----------
gwas <- fread(cmd = "gzip -cd data/GCST90269903.tsv.gz", sep = "\t", showProgress = TRUE)
setnames(gwas, c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                 "beta", "standard_error", "effect_allele_frequency", "p_value"))
gwas[, chrpos := paste(chromosome, base_pair_location, sep = ":")]
cat("GWAS variants:", nrow(gwas), " genome-wide sig (p<5e-8):",
    sum(gwas$p_value < 5e-8), "\n")
gwas_key <- gwas[p_value < 1]  # 保留全部，索引加速：按 chrpos 建 key
setkey(gwas_key, chrpos)

## ---------- 2. 读入 eQTLGen 显著 cis-eQTL（流式按基因过滤） ----------
eqtl_file <- "data/eqtlgen_cis_significant.txt.gz"
hdr <- fread(cmd = paste("gzip -cd", shQuote(eqtl_file), "| head -2"), sep = "\t", header = TRUE)
cat("eQTLGen columns:", paste(names(hdr), collapse = ", "), "\n")
# 期望列: Pvalue SNP SNPChr SNPPos AssessedAllele OtherAllele Zscore Gene GeneSymbol ...
eqtl <- fread(cmd = paste("gzip -cd", shQuote(eqtl_file)),
              select = c("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                         "OtherAllele", "Zscore", "Gene", "GeneSymbol"),
              sep = "\t", showProgress = TRUE, quote = "")
eqtl[, chrpos := paste(SNPChr, SNPPos, sep = ":")]
cat("eQTL rows:", nrow(eqtl), " genes:", uniqueN(eqtl$GeneSymbol), "\n")

## ---------- 3. 逐基因 MR（cis 窗口内显著 eQTL 作为 IV，按 chr:pos 匹配 GWAS） ----------
mr_one_gene <- function(sym, fdr_th = 1e-5, gap = 1e5, max_iv = 20) {
  sub <- eqtl[GeneSymbol == sym & Pvalue < fdr_th]
  if (!nrow(sub)) return(list(ivs = 0))
  setorder(sub, Pvalue)
  if (nrow(sub) > max_iv) sub <- sub[1:max_iv]          # 取最强 20 个再做距离剪枝
  kept <- sub[SNPChr == sub$SNPChr[1]]
  ord <- kept[order(SNPPos)]
  gap_ok <- c(TRUE, diff(ord$SNPPos) > gap)
  pruned <- ord[gap_ok]
  n_iv <- nrow(pruned)
  m <- pruned[chrpos %in% gwas_key$chrpos]
  m <- cbind(m, gwas_key[match(m$chrpos, gwas_key$chrpos), .(g_beta = beta, g_se = standard_error,
                                                              g_ea = effect_allele, g_oa = other_allele,
                                                              g_p = p_value)])
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  flip <- (toupper(m$AssessedAllele) == toupper(m$g_oa)) & (toupper(m$OtherAllele) == toupper(m$g_ea))
  m[, beta_e := Zscore]
  m[flip, beta_e := -beta_e]
  m <- m[toupper(g_ea) %in% names(comp) & toupper(g_oa) %in% names(comp)]
  if (!nrow(m)) return(list(ivs = 0))
  b_exp <- m$beta_e; b_out <- m$g_beta; se_out <- m$g_se
  if (nrow(m) == 1) {                                   # 单 IV：Wald 比值
    ivw_b <- b_out / b_exp
    # 单 IV 的 SE 需 beta_exp 的 se，Z 无法精确传递；以 GWAS se/|Z| 近似
    ivw_se <- abs(se_out / b_exp)
    method <- "Wald ratio (Z-approx)"
  } else {
    ivw_b <- sum(b_exp * b_out / se_out^2) / sum(b_exp^2 / se_out^2)
    ivw_se <- sqrt(1 / sum(b_exp^2 / se_out^2))
    method <- "IVW (Z-approx)"
  }
  z <- ivw_b / ivw_se
  p <- 2 * pnorm(-abs(z))
  q <- NA; q_df <- NA; q_p <- NA
  if (nrow(m) > 1) {
    q <- sum(((b_out - b_exp * ivw_b) / se_out)^2); q_df <- nrow(m) - 1
    q_p <- pchisq(q, q_df, lower.tail = FALSE)
  }
  list(ivs = nrow(m), beta = ivw_b, se = ivw_se, p = p, OR = exp(ivw_b),
       Q = q, Q_p = q_p, method = method, detail = m)
}
res <- lapply(core, mr_one_gene)
names(res) <- core
out <- rbindlist(lapply(core, function(g) {
  r <- res[[g]]
  data.table(gene = g, n_IV = r$ivs,
             OR = if (!is.null(r$OR)) r$OR else NA,
             beta = if (!is.null(r$beta)) r$beta else NA,
             se = if (!is.null(r$se)) r$se else NA,
             P_MR = if (!is.null(r$p)) r$p else NA,
             Q_p = if (!is.null(r$Q_p)) r$Q_p else NA)
}), fill = TRUE)
fwrite(out, "results/mr/mr_summary.csv")
print(out)

## ---------- 4. 共定位（有 ≥50 个重叠 cis SNP 的基因） ----------
locus_one <- function(sym) {
  sub <- eqtl[GeneSymbol == sym]
  if (!nrow(sub)) return(NULL)
  chr <- sub$SNPChr[1]
  peak <- sub[which.min(Pvalue)]$SNPPos
  win <- sub[SNPChr == chr & SNPPos > peak - 5e5 & SNPPos < peak + 5e5]
  gw <- gwas[chromosome == chr & base_pair_location %in% win$SNPPos]
  if (nrow(gw) < 50) return(NULL)
  d <- merge(win[, .(chrpos, p_eqtl = Pvalue)], gw[, .(chrpos, p_gwas = p_value)], by = "chrpos")
  if (nrow(d) < 50) return(NULL)
  cl <- suppressWarnings(coloc.abf(
    dataset1 = list(snp = d$chrpos, pvalues = d$p_eqtl, type = "quant", N = 31684),
    dataset2 = list(snp = d$chrpos, pvalues = d$p_gwas, type = "cc",
                    N = 296824, s = 16743 / 296824)))  # Tyrmi 2023: PE 16743 例 / 共 296824
  list(gene = sym, PP.H0 = cl$summary["PP.H0.abf"], PP.H1 = cl$summary["PP.H1.abf"],
       PP.H3 = cl$summary["PP.H3.abf"], PP.H4 = cl$summary["PP.H4.abf"], n = nrow(d))
}
loc <- lapply(core, locus_one)
loc_res <- rbindlist(compact(loc), fill = TRUE)
if (nrow(loc_res)) { fwrite(loc_res, "results/mr/coloc_summary.csv"); print(loc_res) }
cat("MR SCRIPT DONE\n")
