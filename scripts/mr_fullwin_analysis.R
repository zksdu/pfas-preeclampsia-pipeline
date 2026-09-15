# mr_fullwin_analysis.R — v3 加强版：full 版 eQTL 全行（Z+MAF 重建）做定量 MR + coloc
# 依据：eQTLGen 官方不发布原生 Beta/SE（仅 Z），Z+MAF 重建为该数据标准做法；
# full 版含亚阈值变异 → IV 候选更多、coloc 窗口覆盖更全（v2 仅显著行）
options(stringsAsFactors = FALSE, timeout = 3600)
suppressPackageStartupMessages({library(data.table); library(coloc)})
setwd("E:/Zcode/0911")

core <- readLines("results/ml/core_genes.txt")

## ---------- 1. PE GWAS ----------
gwas <- fread(cmd = "gzip -cd data/GCST90269903.tsv.gz", sep = "\t", showProgress = TRUE)
setnames(gwas, c("chromosome", "base_pair_location", "effect_allele", "other_allele",
                 "beta", "standard_error", "effect_allele_frequency", "p_value"))
gwas[, chrpos := paste(chromosome, base_pair_location, sep = ":")]
cat("GWAS variants:", nrow(gwas), "\n")

## ---------- 2. full 版核心基因 eQTL ----------
eqtl <- fread(cmd = "gzip -cd data/mr_core_eqtl_beta.gz", sep = "\t", quote = "", showProgress = TRUE)
cat("core-gene eQTL rows:", nrow(eqtl), " genes:", uniqueN(eqtl$GeneSymbol), "\n")
stopifnot(all(c("SNP", "Zscore", "NrSamples", "GeneSymbol") %in% names(eqtl)))
eqtl[, chrpos := paste(SNPChr, SNPPos, sep = ":")]

af <- fread(cmd = "gzip -cd data/eqtlgen_snp_af.txt.gz", sep = "\t", quote = "",
            select = c("SNP", "AlleleB_all"), showProgress = TRUE)
setnames(af, "AlleleB_all", "AF_B")
eqtl <- merge(eqtl, af[, .(SNP, AF_B)], by = "SNP")
cat("after AF merge:", nrow(eqtl), "\n")
setorder(eqtl, GeneSymbol, SNP, Pvalue)
eqtl <- eqtl[!duplicated(eqtl[, .(GeneSymbol, SNP)])]
eqtl <- eqtl[!duplicated(eqtl[, .(GeneSymbol, chrpos)])]

## ---------- 3. Beta/SE 重建（eQTLGen 标准化表达，sdY=1） ----------
eqtl[, maf := pmin(AF_B, 1 - AF_B)]
eqtl[, se_rec := 1 / sqrt(NrSamples * 2 * maf * (1 - maf))]
eqtl[, beta_rec := Zscore * se_rec]

## ---------- 4. 逐基因定量 MR（p<1e-5，100kb 剪枝，max 30 IV） ----------
mr_one_gene <- function(sym, fdr_th = 1e-5, gap = 1e5, max_iv = 30) {
  sub <- eqtl[GeneSymbol == sym & Pvalue < fdr_th & is.finite(beta_rec) & is.finite(se_rec)]
  if (!nrow(sub)) return(list(ivs = 0))
  setorder(sub, Pvalue)
  if (nrow(sub) > max_iv) sub <- sub[1:max_iv]
  ord <- sub[order(SNPPos)]
  pruned <- ord[c(TRUE, diff(ord$SNPPos) > gap)]
  m <- pruned[chrpos %in% gwas$chrpos]
  m <- cbind(m, gwas[match(m$chrpos, gwas$chrpos), .(g_beta = beta, g_se = standard_error,
                                                     g_ea = effect_allele, g_oa = other_allele,
                                                     g_p = p_value)])
  comp <- c(A = "T", T = "A", C = "G", G = "C")
  ea <- toupper(m$AssessedAllele); oa <- toupper(m$OtherAllele)
  flip <- (ea == toupper(m$g_oa)) & (oa == toupper(m$g_ea))
  same <- (ea == toupper(m$g_ea)) & (oa == toupper(m$g_oa))
  m <- m[(same | flip) & ea %in% names(comp) & oa %in% names(comp)]
  if (!nrow(m)) return(list(ivs = 0))
  m[same, beta_x := beta_rec]; m[same, se_x := se_rec]
  m[flip & !same, beta_x := -beta_rec]; m[flip & !same, se_x := se_rec]
  b_x <- m$beta_x; se_x <- m$se_x; b_y <- m$g_beta; se_y <- m$g_se
  if (nrow(m) == 1) {
    b <- b_y / b_x; se <- sqrt((se_y / b_x)^2 + (se_x * b_y / b_x^2)^2)
    method <- "Wald ratio"
  } else {
    b <- sum(b_x * b_y / se_y^2) / sum(b_x^2 / se_y^2)
    se <- sqrt(1 / sum(b_x^2 / se_y^2))
    method <- "IVW"
  }
  p <- 2 * pnorm(-abs(b / se))
  q_p <- NA
  if (nrow(m) > 1) {
    qq <- sum(((b_y - b_x * b) / se_y)^2); q_p <- pchisq(qq, nrow(m) - 1, lower.tail = FALSE)
  }
  list(ivs = nrow(m), beta = b, se = se, p = p, OR = exp(b), Q_p = q_p, method = method)
}
out <- rbindlist(lapply(core, function(g) {
  r <- mr_one_gene(g)
  data.table(gene = g, n_IV = r$ivs,
             method = if (!is.null(r$method)) r$method else NA,
             OR = if (!is.null(r$OR)) r$OR else NA,
             beta = if (!is.null(r$beta)) r$beta else NA,
             se = if (!is.null(r$se)) r$se else NA,
             P_MR = if (!is.null(r$p)) r$p else NA,
             Q_p = if (!is.null(r$Q_p) && !is.na(r$Q_p)) r$Q_p else NA)
}), fill = TRUE)
out[, P_FDR := p.adjust(P_MR, method = "BH")]
fwrite(out, "results/mr/mr_v3_fullwin_summary.csv")
print(out)

## ---------- 5. 共定位（full 行，窗口更密） ----------
locus_one <- function(sym) {
  sub <- eqtl[GeneSymbol == sym & !is.na(beta_rec)]
  if (!nrow(sub)) return(NULL)
  chr <- sub$SNPChr[1]
  peak <- sub[which.min(Pvalue)]$SNPPos
  win <- sub[SNPChr == chr & SNPPos > peak - 5e5 & SNPPos < peak + 5e5]
  gw <- gwas[chromosome == chr & base_pair_location %in% win$SNPPos]
  if (nrow(gw) < 50) return(NULL)
  d <- merge(win[, .(chrpos, p_eqtl = Pvalue, b_eqtl = beta_rec, se_eqtl = se_rec, n_eq = NrSamples)],
             gw[, .(chrpos, p_gwas = p_value, b_gwas = beta, se_gwas = standard_error)], by = "chrpos")
  d <- d[!duplicated(d$chrpos)]
  if (nrow(d) < 50) return(NULL)
  cl <- suppressWarnings(coloc.abf(
    dataset1 = list(snp = d$chrpos, beta = d$b_eqtl, varbeta = d$se_eqtl^2,
                    N = median(d$n_eq), type = "quant", sdY = 1),
    dataset2 = list(snp = d$chrpos, beta = d$b_gwas, varbeta = d$se_gwas^2,
                    N = 296824, s = 16743 / 296824, type = "cc")))
  data.table(gene = sym, n = nrow(d),
             PP.H0 = cl$summary["PP.H0.abf"], PP.H1 = cl$summary["PP.H1.abf"],
             PP.H3 = cl$summary["PP.H3.abf"], PP.H4 = cl$summary["PP.H4.abf"])
}
loc <- lapply(core, locus_one)
loc_res <- rbindlist(Filter(Negate(is.null), loc), fill = TRUE)
if (nrow(loc_res)) { fwrite(loc_res, "results/mr/coloc_v3_fullwin_summary.csv"); print(loc_res) }

## ---------- 6. 与 v2 对照 ----------
if (file.exists("results/mr/mr_full_summary.csv")) {
  v2 <- fread("results/mr/mr_full_summary.csv")
  cmp <- merge(out[, .(gene, OR_v3 = OR, P_v3 = P_MR)],
               v2[, .(gene, OR_v2 = OR, P_v2 = P_MR)], by = "gene", all = TRUE)
  fwrite(cmp, "results/mr/mr_v3_vs_v2.csv")
  print(cmp)
}
cat("MR V3 FULLWIN DONE\n")
