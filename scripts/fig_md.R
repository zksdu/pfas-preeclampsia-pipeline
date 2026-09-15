# fig_md.R — Fig.8: 三体系 MD 组图（RMSD蛋白 / Rg / RMSD配体 / 接触+氢键）
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
setwd("E:/Zcode/0911")

genes <- c("ENG", "FLT1", "TMEM45A")
cols3 <- c("#1b9e77", "#d95f02", "#7570b3")
names(cols3) <- genes

read_metrics <- function(g) {
  m <- fread(sprintf("results/md_gpu/%s_PFOS_md_metrics.csv", g))
  m[, gene := g]
  m[, t := time_ns]
  m
}
read_contacts <- function(g) {
  c <- fread(sprintf("results/md_gpu/%s_PFOS_md_contacts.csv", g))
  c[, gene := g]
  c[, t := time_ns]
  c
}
M <- rbindlist(lapply(genes, read_metrics))
C <- rbindlist(lapply(genes, read_contacts))

png("results/figures/fig_md.png", width = 3000, height = 2400, res = 210)
par(mfrow = c(2, 2), mar = c(4.5, 4.8, 3.5, 1.5))

rng <- function(v) c(0, max(v) * 1.08)

## (a) 蛋白 RMSD
plot(NA, NA, xlim = c(0, 10), ylim = rng(M$rmsd_protein_A),
     xlab = "time (ns)", ylab = "Protein RMSD (Å)", main = "(a) Protein RMSD")
for (g in genes) {
  d <- M[gene == g]
  lines(d$t, d$rmsd_protein_A, col = cols3[[g]], lwd = 1.6)
}
legend("bottomright", genes, col = cols3, lwd = 2, bty = "n", cex = 0.9)

## (b) Rg
plot(NA, NA, xlim = c(0, 10),
     ylim = range(M$rg_protein_A) + c(-1, 1) * max(M$rg_protein_A) * 0.06,
     xlab = "time (ns)", ylab = "Protein Rg (Å)", main = "(b) Radius of gyration")
for (g in genes) {
  d <- M[gene == g]
  lines(d$t, d$rg_protein_A, col = cols3[[g]], lwd = 1.6)
}
legend("bottomright", genes, col = cols3, lwd = 2, bty = "n", cex = 0.9)

## (c) 配体 RMSD
plot(NA, NA, xlim = c(0, 10), ylim = rng(M$rmsd_ligand_A),
     xlab = "time (ns)", ylab = "PFOS RMSD from docked pose (Å)", main = "(c) Ligand RMSD")
for (g in genes) {
  d <- M[gene == g]
  lines(d$t, d$rmsd_ligand_A, col = cols3[[g]], lwd = 1.6)
}
legend("bottomright", genes, col = cols3, lwd = 2, bty = "n", cex = 0.9)

## (d) 配体-蛋白接触与氢键（TMEM45A 有间歇接触；FLT1/ENG 为 0 接触、离子介导氢键）
ymax <- max(C$heavy_atom_contacts, 4)
plot(NA, NA, xlim = c(0, 10), ylim = c(-0.6, ymax * 1.1),
     xlab = "time (ns)", ylab = "Heavy-atom contacts (<4 Å) / H-bonds",
     main = "(d) PFOS-protein contacts")
for (g in genes) {
  d <- C[gene == g]
  lines(d$t, d$heavy_atom_contacts, col = cols3[[g]], lwd = 1.2)
}
points(C$heavy_atom_contacts == 0, type = "n")
mtext("H-bond means: ENG 0.05, FLT1 0.34, TMEM45A 0.22 per frame", side = 3,
      line = -2.2, cex = 0.75, adj = 0.02, col = "grey30")
legend("topright", c("contacts (<4 Å)", "zero-contact baseline"),
       col = c("grey30", "grey80"), lwd = c(1.5, 1), lty = c(1, 2), bty = "n", cex = 0.85)
abline(h = 0, col = "grey80", lty = 2)

dev.off()
cat("fig_md.png saved\n")
