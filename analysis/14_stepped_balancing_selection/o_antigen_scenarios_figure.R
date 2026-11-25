# The balancing-selection scenarios, drawn with the O-antigen data
#
# The three scenarios usually drawn for an accessory gene (genic selection
# alone; genic plus diversifying selection on sequence; genic plus long-term
# maintenance of alleles) all take absence to mean absence. The O antigen does
# not work that way. Every one of the 260 genomes carries a cluster in the
# gnd-hisF interval, and the 73 genomes scored as lacking wzm are exactly the
# genomes carrying a non-O1/O2 cluster there. Presence/absence of any galactan
# gene is therefore a readout of which locus occupies the position, not of that
# gene being gained and lost, which is scenario 0 below. Scenarios 1 to 3 then
# apply within the maintained locus.
#
# Panel (a) is the schematic. Panels (b) to (g) are the evidence: the protein
# allele spectrum of each galactan gene, partitioned by O subtype. The six genes
# share one genealogy (they are in one operon, in the same genomes, with
# s/site 0.227 to 0.310 and 3.1 to 4.3 times the segregating sites the
# coalescent null expects), so any difference between their panels is a
# difference in selection on protein sequence, not in genealogical depth.
#
# What the panels show: wzm collapses to a single allele in all four subtypes,
# while glf, wbbM, wbbN and wbbO carry subtype-specific alleles at the same
# genealogical depth. None of them is diversifying: piN/piS is 0.014 to 0.145
# across the operon, so scenario 2 is empty here.

# Parameters ---------------------------------------------------------------

min_sero  <- 5L     # O subtypes with fewer carriers are dropped from (b) to (g)
min_al    <- 3L     # alleles rarer than this are pooled into "other"
lab_cex   <- 0.62

# the operon, in locus order, with the tidy names used in the manuscript
focal_dt <- data.table(
  gene_family = c("g004051", "g003783", "g000322", "g001147", "g000817", "g000954"),
  short       = c("wzm", "wzt", "wbbM", "glf", "wbbN", "wbbO"))

# blues for O1, oranges for O2, as in o_antigen_variant_lollipops.R
sero_col <- c("O1v1" = "#2C6E8F", "O1v2" = "#8FC1D4",
              "O2v1" = "#C4631E", "O2v2" = "#F3BE8A")

# Data ---------------------------------------------------------------------

meta <- fread(paste0(outdir_dat, "/olocus.csv"), select = c("geno_id", "O_locus"))

# H, s/site and obs/exp come from the locus table, piN/piS from folder 14, so
# the figure quotes exactly the numbers the rest of the manuscript does
loc <- unique(fread(paste0(outdir_dat, "/o_antigen_locus_table.csv"),
                    select = c("gene_family", "k", "s_per_site", "obs_exp", "H", "n_al")),
              by = "gene_family")
pin <- fread(paste0(outdir_dat, "/ag_piN_piS_dt.csv"),
             select = c("gene_family", "pi_n", "pi_s", "piN_piS"))

focal_dt <- loc[focal_dt, on = .(gene_family)]
focal_dt <- pin[focal_dt, on = .(gene_family)]
setorder(focal_dt, -H)      # wzm first, so the contrast runs left to right

gc_tab <- Biostrings::GENETIC_CODE

# per-genome peptide, from the in-frame alignments written by folder 06. A codon
# holding a gap or an ambiguity code translates to "-" rather than to a residue
pep_of <- function(fam) {
  dna <- readDNAStringSet(paste0(outdir_dat, "/ag_fus_alns/", fam, ".fasta"))
  pep <- vapply(as.character(dna), function(x) {
    cd <- substring(x, seq(1, nchar(x) - 2, 3), seq(3, nchar(x), 3))
    a  <- gc_tab[cd]; a[is.na(a)] <- "-"
    paste(a, collapse = "")
  }, character(1))
  data.table(geno_id = sub("_[0-9]+$", "", names(dna)), pep = pep)
}

h_nei <- function(x) { n <- length(x); p <- prop.table(table(x)); (n/(n-1)) * (1 - sum(p^2)) }

# Nei's H partitioned within and between the O subtypes. Gst = 1 - Hw/Ht is the
# share of the allelic diversity that sits between subtypes rather than within
gst_of <- function(pep, sero) {
  w  <- data.table(sero, pep)[, .(h = h_nei(pep), n = .N), by = sero][n > 1]
  Ht <- h_nei(pep)
  c(Ht = Ht, Hw = sum(w$h * w$n) / sum(w$n), Gst = 1 - (sum(w$h * w$n) / sum(w$n)) / Ht)
}

dat <- vector("list", nrow(focal_dt))
names(dat) <- focal_dt$short

for (i in seq_len(nrow(focal_dt))) {
  d <- meta[pep_of(focal_dt$gene_family[i]), on = .(geno_id)]
  d <- d[O_locus %in% names(sero_col)]
  # alleles numbered by overall frequency, so allele 1 is the commonest
  d[, allele := as.integer(factor(pep, levels = names(sort(table(pep), decreasing = TRUE))))]
  g <- gst_of(d$pep, d$O_locus)
  focal_dt[i, c("Ht", "Hw", "Gst") := as.list(g)]
  dat[[i]] <- d
}

print(focal_dt[, .(short, k, n_al, H, s_per_site, obs_exp, pi_n, pi_s, piN_piS, Gst)])
fwrite(focal_dt, paste0(outdir_dat, "/o_antigen_scenario_stats.csv"))

# Figure --------------------------------------------------------------------

png(paste0(outdir_fig, "/o_antigen_scenarios.png"),
    width = 30, height = 24, units = "cm", res = 600,
    pointsize = 12, type = "cairo")

layout(matrix(c(1, 1, 1,
                2, 3, 4,
                5, 6, 7), nrow = 3, byrow = TRUE),
       heights = c(1.1, 1, 1))

# (a) the schematic ---------------------------------------------------------
# Four columns, each a handful of genomes drawn as a line with the two flanking
# core genes (gnd, hisF) as black blocks and the cluster between them. Colour is
# the allele, red ticks are non-synonymous differences, as in the original.

par(mar = c(0.4, 2.0, 2.2, 0.6))
plot(c(0, 100), c(0, 100), type = "n", axes = FALSE, xlab = "", ylab = "")

scen <- list(
  list(x = 2,  ttl = "0  Alternative loci at one position",
       sub = "several clusters compete for one\nchromosomal position; p/a of any\nconstituent gene reports which\ncluster is there",
       sig = "excess S, any H",
       obs = "the galactan locus itself\n185 genomes vs 73 carrying\nanother cluster"),
  list(x = 27, ttl = "1  Genic selection only",
       sub = "the cluster is maintained, its\nsequence is not diversified",
       sig = expression("excess " * italic(S) * ",  low " * italic(H) * ",  " * pi[N]/pi[S] * " << 1"),
       obs = "wzm"),
  list(x = 52, ttl = "2  Genic + diversifying selection",
       sub = "selection also makes many protein\nalleles (eg NFDS on an epitope)",
       sig = expression("excess " * italic(S) * ",  high " * italic(H) * ",  " * pi[N]/pi[S] * " > 1"),
       obs = "none"),
  list(x = 77, ttl = "3  Genic + maintenance of alleles",
       sub = "alleles are held long term on\ndistinct backgrounds (here the\nfour O subtypes, not habitats)",
       sig = expression("excess " * italic(S) * ",  high " * italic(H) * ",  " * pi[N]/pi[S] * " << 1"),
       obs = "wbbM, wbbN, wbbO, glf\nwzt, at intermediate H"))

# one genome: flanking core blocks, the backbone, and the cluster if carried
genome <- function(x0, y, w = 20, fill = NA, nsyn = 0, seed = 1, cw = 8.5) {
  segments(x0, y, x0 + w, y, lwd = 1.1, col = "grey25")
  rect(x0, y - 1.6, x0 + 2.2, y + 1.6, col = "grey15", border = NA)             # gnd
  rect(x0 + w - 2.2, y - 1.6, x0 + w, y + 1.6, col = "grey15", border = NA)     # hisF
  if (!is.na(fill)) {
    rect(x0 + 6.5, y - 1.9, x0 + 6.5 + cw, y + 1.9, col = fill, border = "grey40", lwd = 0.5)
    if (nsyn > 0) {
      set.seed(seed)
      tk <- sort(runif(nsyn, x0 + 7.2, x0 + 5.8 + cw))
      segments(tk, y - 1.9, tk, y + 1.9, col = "#D4342C", lwd = 1.1)
    }
  }
}

ys <- c(88, 80, 72, 64)

# 0: a different cluster in each genome, all in the one position
# different clusters hold different genes, so they differ in length too
for (j in 1:4) genome(scen[[1]]$x, ys[j], fill = c("#2C6E8F", "#C4631E", "#6BA368", "#9B7BB8")[j],
                      nsyn = 0, seed = j, cw = c(8.5, 10.5, 6.5, 9.5)[j])
# 1: one allele, carried by three of four
for (j in 1:4) genome(scen[[2]]$x, ys[j], fill = if (j == 2) NA else "#9CC3D5", nsyn = 0)
# 2: many alleles, many non-synonymous differences
for (j in 1:4) genome(scen[[3]]$x, ys[j], fill = if (j == 2) NA else "#9CC3D5",
                      nsyn = if (j == 2) 0 else 7, seed = 10 + j)
# 3: distinct alleles held on distinct backgrounds, few non-synonymous
for (j in 1:4) genome(scen[[4]]$x, ys[j],
                      fill = if (j == 2) NA else c("#2C6E8F", NA, "#C4631E", "#F3BE8A")[j],
                      nsyn = if (j == 2) 0 else 1, seed = 20 + j)

for (s in scen) {
  text(s$x, 98, s$ttl, adj = 0, font = 2, cex = 0.78, xpd = NA)
  text(s$x, 56, s$sub, adj = c(0, 1), cex = 0.62, col = "grey25")
  text(s$x, 36, s$sig, adj = c(0, 1), cex = 0.66)
  text(s$x, 26, "observed:", adj = c(0, 1), cex = 0.62, col = "grey45")
  text(s$x, 21, s$obs, adj = c(0, 1), cex = 0.66, font = 3)
}
segments(c(25.5, 50.5, 75.5), 2, c(25.5, 50.5, 75.5), 95, col = "grey85", lwd = 0.8)
mtext("a", side = 3, line = 0.6, adj = 0, font = 2, cex = 1.1,
      at = grconvertX(0.005, "nfc", "user"))

# (b) to (g) the allele spectra ---------------------------------------------
# One bar per O subtype, stacked by protein allele. If the alleles are
# indifferent to subtype the four bars carry the same colour; if alleles are
# held on subtype backgrounds the bars differ. Bar width tracks the number of
# genomes, so the subtypes are not given equal visual weight.

for (i in seq_len(nrow(focal_dt))) {

  d  <- dat[[i]]
  lv <- names(sero_col)[names(sero_col) %in% d$O_locus]
  n_s <- d[, .N, by = O_locus][N >= min_sero]
  lv <- lv[lv %in% n_s$O_locus]

  # alleles common enough to be worth a colour of their own; the rest are "other"
  big <- d[, .N, by = allele][N >= min_al][order(allele)]$allele
  d[, grp := fifelse(allele %in% big, as.character(allele), "other")]
  pal <- c(setNames(hcl.colors(max(length(big), 2), "Batlow")[seq_along(big)],
                    as.character(big)), other = "grey82")

  tb <- as.matrix(table(factor(d$grp, levels = names(pal)),
                        factor(d$O_locus, levels = lv)))
  tb <- prop.table(tb, 2)

  par(mar = c(3.2, 4.2, 2.6, 0.8))
  n_by <- as.vector(table(factor(d$O_locus, levels = lv)))
  w    <- sqrt(n_by) / max(sqrt(n_by))
  bp <- barplot(tb, col = pal, border = "white", space = 0.35, width = w,
                las = 1, ylim = c(0, 1.02), ylab = "", xlab = "", axisnames = FALSE)
  mtext("Proportion of carriers", side = 2, line = 2.8, cex = 0.62)
  mtext(lv, side = 1, at = bp, line = 0.45, cex = lab_cex)
  mtext(sprintf("n = %d", n_by), side = 1, at = bp, line = 1.45, cex = 0.5, col = "grey45")

  mtext(bquote(italic(.(focal_dt$short[i])) ~ "  " * .(focal_dt$n_al[i]) ~ "alleles"),
        side = 3, line = 1.1, adj = 0, cex = 0.72, font = 2)
  mtext(bquote(italic(H) == .(sprintf("%.2f", focal_dt$H[i])) ~ "  " *
               G[ST] == .(sprintf("%.2f", focal_dt$Gst[i])) ~ "  " *
               pi[N]/pi[S] == .(sprintf("%.3f", focal_dt$piN_piS[i])) ~ "  " *
               "o/e" == .(sprintf("%.1f", focal_dt$obs_exp[i]))),
        side = 3, line = 0.15, adj = 0, cex = 0.56, col = "grey30")
  mtext(letters[i + 1], side = 3, line = 1.1, adj = 0, font = 2, cex = 1.0,
        at = grconvertX(0.01, "nfc", "user"))
}

dev.off()
