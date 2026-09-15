# ctd_targets.R — Step 1a: 从 CTD 批量互作数据提取 PFOA/PFOS 的人源 curated 基因
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
setwd("E:/Zcode/0911")
dir.create("results/targets", showWarnings = FALSE)

f <- "data/CTD_chem_gene_ixns.tsv.gz"
stopifnot(file.exists(f))
# 表头在第 28 行（前 27 行为 '#' 元信息）
ctd <- fread(cmd = paste("gzip -cd", shQuote(f)), skip = 27, sep = "\t", quote = "",
             showProgress = FALSE, header = TRUE)
setnames(ctd, c("ChemicalName", "ChemicalID", "CasRN", "GeneSymbol", "GeneID",
                "GeneForms", "Organism", "OrganismID", "Interaction",
                "InteractionActions", "PubMedIDs"))
cat("total rows:", nrow(ctd), "\n")

hsa <- ctd[Organism == "Homo sapiens"]
cat("human rows:", nrow(hsa), "\n")
pfoa <- hsa[CasRN == "335-67-1"]
pfos <- hsa[CasRN == "1763-23-1"]
cat("PFOA rows:", nrow(pfoa), " PFOS rows:", nrow(pfos), "\n")

g_pfoa <- sort(unique(pfoa$GeneSymbol))
g_pfos <- sort(unique(pfos$GeneSymbol))
cat("PFOA genes:", length(g_pfoa), " PFOS genes:", length(g_pfos), "\n")
genes <- sort(unique(c(g_pfoa, g_pfos)))
cat("union PFAS genes:", length(genes), "\n")

writeLines(genes, "results/targets/pfas_ctd_targets.txt")
writeLines(g_pfoa, "results/targets/pfoa_genes.txt")
writeLines(g_pfos, "results/targets/pfos_genes.txt")
# 保留互作类型信息（供论文 Table S2）
ivn <- rbind(pfoa[, .(chem = "PFOA", GeneSymbol, InteractionActions)],
             pfos[, .(chem = "PFOS", GeneSymbol, InteractionActions)])
fwrite(unique(ivn), "results/targets/pfas_gene_interactions.csv")
print(head(genes, 40))
cat("CTD DONE\n")
