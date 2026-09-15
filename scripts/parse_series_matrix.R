# parse_series_matrix.R — 手动解析 GEO series matrix（不依赖 GEOquery）
# 用法: Rscript parse_series_matrix.R <file.txt.gz> <输出前缀>
args <- commandArgs(trailingOnly = TRUE)
gz <- args[1]; prefix <- args[2]

lines <- readLines(gzfile(gz), warn = FALSE)

meta <- lines[startsWith(lines, "!Series_") | startsWith(lines, "!Sample_")]
writeLines(meta, paste0(prefix, "_meta.txt"))

get_field <- function(name) {
  l <- lines[lines == paste0(name, "_table_begin") + 0]
}
# 提取表达矩阵
beg <- grep("^!series_matrix_table_begin$", lines)
end <- grep("^!series_matrix_table_end$", lines)
if (length(beg) == 0 || length(end) == 0) stop("no expression table found")
tab <- lines[(beg[1] + 1):(end[1] - 1)]
con <- textConnection(tab)
expr <- read.delim(con, sep = "\t", header = TRUE, check.names = FALSE,
                   quote = "", stringsAsFactors = FALSE)
close(con)
rownames(expr) <- expr$ID_REF
expr$ID_REF <- NULL
expr <- as.matrix(expr)
cat("expr dim:", nrow(expr), "x", ncol(expr), "\n")

# 提取样本注释
skeys <- c("!Sample_geo_accession", "!Sample_title", "!Sample_source_name",
           "!Sample_characteristics_ch1", "!Sample_platform_id", "!Sample_data_processing")
sdat <- lapply(skeys, function(k) {
  idx <- which(lines == k)
  if (!length(idx)) return(NULL)
  v <- lines[idx + 1]
  strsplit(gsub('"', "", v), "\t", fixed = TRUE)[[1]]
})
names(sdat) <- skeys
n <- length(sdat[["!Sample_geo_accession"]])
sdat <- lapply(sdat, function(x) if (is.null(x)) rep(NA, n) else x)
pheno <- as.data.frame(sdat, stringsAsFactors = FALSE)
names(pheno) <- sub("^!", "", names(pheno))
write.csv(pheno, paste0(prefix, "_pheno.csv"), row.names = FALSE)
cat("pheno dim:", nrow(pheno), "x", ncol(pheno), "\n")
print(head(pheno, 3))

saveRDS(list(expr = expr, pheno = pheno), paste0(prefix, ".rds"))
cat("saved", paste0(prefix, ".rds"), "\n")
