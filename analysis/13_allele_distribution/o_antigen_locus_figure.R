# The O-antigen locus: one position, alternative clusters, one genealogy
#
# (a) The interval between gnd and hisF in each serotype. Every genome carries
#     a cluster here and the flanks are present in all 260, so what a pangenome
#     records as presence/absence of wzm, wzt, glf and wbbM is really which
#     cluster occupies the locus. Fill colour is the gene symbol, so a symbol
#     that recurs across serotypes keeps its colour even where the PIRATE
#     family differs; a heavy border marks exactly that case, a functional
#     equivalent scored as a different gene.
#
# (b,c) The galactan cluster alone, all genes on the same 182-185 carriers.
#     Synonymous diversity per site is near constant along the locus and every
#     gene is called high, while H ranges from 0.13 to 0.89. The genealogy is a
#     property of the locus; H is a property of each protein's constraint.
#
# Reads o_antigen_locus_table.csv (o_antigen_locus_table.R) and
# olocus_genes.csv.

library(data.table)

sero_min <- 8          # serotypes with at least this many genomes get a row
ref_sero <- "O2v2"     # the cluster the heavy borders are judged against
focal    <- c("g004051", "g000322")   # wzm and wbbM, the two focal AGs
foc_col  <- "#D62828"

tab  <- fread(paste0(outdir_dat, "/o_antigen_locus_table.csv"))
loci <- fread(paste0(outdir_dat, "/olocus_genes.csv"))

# 1. a representative gene order per serotype -------------------------------

main_anchor <- loci[anchor != "-", .N, by = anchor][order(-N)][1, anchor]
main_msu    <- loci[anchor == main_anchor, .N, by = msu][order(-N)][1, msu]
cat("locus syntenic block:", main_anchor, main_msu, "\n")

ser <- tab[n_genomes >= sero_min]
ser <- unique(loci[, .(gene_family, anchor, msu)])[ser, on = .(gene_family)]
cat("families drawn:", uniqueN(ser$gene_family),
    "  not placed in the block:", uniqueN(ser[anchor != main_anchor, gene_family]), "\n")
setorder(ser, -n_genomes, O_locus, pos)
sero_lev <- unique(ser$O_locus)

# the family each symbol takes in the reference cluster, to spot equivalents
ref_fam <- ser[O_locus == ref_sero, setNames(gene_family, tidy_name)]
ser[, swapped := !is.na(ref_fam[tidy_name]) & ref_fam[tidy_name] != gene_family]

# colours by symbol: recurring symbols get one, one-offs stay grey
sym_n <- ser[, uniqueN(O_locus), by = tidy_name][V1 >= 2][order(-V1), tidy_name]
pal   <- c("#B22222", "#3FA9F5", "#E8A33D", "#7FB069", "#7F5BA8", "#F26C2A",
           "#4C9A8F", "#C45CA0", "#8C6D3F", "#5B7FBF", "#9BBF5B", "#BF5B8C")
sym_col <- setNames(rep(pal, length.out = length(sym_n)), sym_n)
ser[, fill := fifelse(tidy_name %chin% sym_n, sym_col[tidy_name], "grey88")]

# 2. the galactan cluster, in locus order -----------------------------------

gal <- ser[O_locus == ref_sero & !is.na(H)]
setorder(gal, pos)
gal <- gal[tidy_name %chin% c("wzm", "wzt", "glf", "wbbN", "wbbO", "wbbM")]

# 3. figure ------------------------------------------------------------------

png(paste0(outdir_fig, "/o_antigen_locus.png"),
    width = 19, height = 12.6, units = "cm", res = 600,
    pointsize = 10, type = "cairo")

par(mar = c(1.2, 7.2, 4.4, 1.2))
nmax <- ser[, .N, by = O_locus][, max(N)]
ny   <- length(sero_lev)
plot(NA, xlim = c(-0.9, nmax + 0.9), ylim = c(ny + 0.6, -0.1),
     axes = FALSE, xlab = "", ylab = "")

for (i in seq_along(sero_lev)) {
  d <- ser[O_locus == sero_lev[i]]
  setorder(d, pos)
  # the backbone, with the two universal flanks
  segments(-0.75, i, nrow(d) + 0.75, i, col = "grey70", lwd = 1)
  text(-0.78, i, "gnd", adj = 1, cex = 0.6, col = "grey35", font = 3, xpd = NA)
  text(nrow(d) + 0.78, i, "hisF", adj = 0, cex = 0.6, col = "grey35",
       font = 3, xpd = NA)
  rect(seq_len(nrow(d)) - 0.44, i - 0.3, seq_len(nrow(d)) + 0.44, i + 0.3,
       col = d$fill, border = fifelse(d$swapped, "black", "white"),
       lwd = fifelse(d$swapped, 1.8, 0.6))
  text(seq_len(nrow(d)), i, substr(d$tidy_name, 1, 13), cex = 0.5, srt = 90,
       col = fifelse(d$fill == "grey88", "grey30", "white"))
  # the two focal families, wherever they occur
  fo <- which(d$gene_family %chin% focal)
  if (length(fo)) {
    rect(fo - 0.5, i - 0.37, fo + 0.5, i + 0.37, border = foc_col,
         lwd = 2.2, xpd = NA)
  }
  # a family PIRATE did not place in the locus anchor
  off <- which(d$anchor != main_anchor)
  if (length(off))
    points(off, rep(i + 0.44, length(off)), pch = 1, cex = 0.5, col = "grey45",
           xpd = NA)
  mtext(sprintf("%s  (n = %d)", sero_lev[i], d$n_genomes[1]), side = 2,
        at = i, las = 1, line = 0.4, cex = 0.68)
}
# the syntenic block the whole interval sits in
ytop <- 0.08
segments(0.5, ytop, nmax + 0.5, ytop, col = "grey35", lwd = 1.1, xpd = NA)
segments(c(0.5, nmax + 0.5), ytop, c(0.5, nmax + 0.5), ytop + 0.1,
         col = "grey35", lwd = 1.1, xpd = NA)
text(mean(c(0.5, nmax + 0.5)), ytop - 0.08,
     sprintf("syntenic block %s  (%s)", main_anchor, main_msu),
     cex = 0.68, col = "grey35", xpd = NA)
text(-0.78, ytop - 0.08, "anchor", adj = 1, cex = 0.6, col = "grey55", xpd = NA)
text(nmax + 0.78, ytop - 0.08, "anchor", adj = 0, cex = 0.6, col = "grey55",
     xpd = NA)

mtext("One locus, alternative clusters", side = 3, line = 3.0, adj = 0,
      cex = 0.95, font = 2, at = -0.9)
mtext("Every genome carries exactly one of these, so the pangenome records allelic variation as gene content",
      side = 3, line = 2.1, adj = 0, cex = 0.66, at = -0.9, col = "grey40")
mtext("Black border: same gene symbol, different PIRATE family.   Open circle: family not placed in the block",
      side = 3, line = 1.1, adj = 0, cex = 0.64, at = -0.9, col = "grey40")
mtext(expression("Red border: the focal AGs " * italic("wzm") *
                 " (g004051) and " * italic("wbbM") * " (g000322)"),
      side = 3, line = 0.2, adj = 0, cex = 0.64, at = -0.9, col = "#D62828")

dev.off()
cat("wrote", paste0(outdir_fig, "/o_antigen_locus.png"), "\n")
