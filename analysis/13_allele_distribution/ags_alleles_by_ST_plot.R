# Peptide alleles of one accessory gene, per genome
#
# 10_outlier_function/ags_alleles_by_ST_plot.R draws allele frequencies by ST
# from the cached allele x ST table. That table is aggregated, so it cannot be
# cut by anything other than ST. Here the alleles of a single focal family are
# rebuilt per genome from its fusion-aware alignment, which keeps the genome
# ID and so lets the same alleles be cut by ST and by isolation source in
# ags_alleles_by_env.R.
#
# Alleles are named the way the extraction script names them, by order of
# first appearance, so a_1 here is a_1 in allele_counts_by_ST.csv. The check
# below confirms that.
#
# Sourced for its objects: al (genome, ST, source, allele), al_cols (one
# colour per allele, ranked by overall frequency) and plot_allele_grp().
# With make_fig TRUE it also draws the two ST panels on their own.
#
# This is description. SpARK sampling is heavily uneven across both STs and
# sources, so the bars say how the alleles that were sampled are arranged,
# and nothing about how common any ST, source or allele is in the population.

library(data.table)
library(Biostrings)

fam      <- if (exists("fam")) fam else "g000990"
gene_lab <- if (exists("gene_lab")) gene_lab else "levR"
min_n    <- 2      # STs with fewer carriers than this are pooled, as in
                   # 10_outlier_function/ags_alleles_by_ST_plot.R
seedy    <- 21     # allele palette seed, matching allelic_diversity.R
n_top    <- 3      # major allele frequencies named in the header
make_fig <- if (exists("make_fig")) make_fig else TRUE

metafile <- "./input_data/spark_metadata.csv"
alndir   <- paste0(outdir_dat, "/ag_fus_alns")

# 1. peptide alleles per genome ---------------------------------------------

aln <- readDNAStringSet(paste0(alndir, "/", fam, ".fasta"))

# the alignment carries gaps; translate each CDS ungapped, as the extraction
# script does when it pulls the locus straight out of the genome
pep <- suppressWarnings(
  Biostrings::translate(DNAStringSet(gsub("-", "", as.character(aln))),
                        if.fuzzy.codon = "solve"))

# header is <geno_id>_<locus number>, so drop the last field
al <- data.table(geno_id = sub("_[^_]+$", "", names(aln)),
                 pep     = as.character(pep))
al[, allele := paste0("a_", match(pep, unique(pep)))]
al[, pep := NULL]

meta <- fread(metafile, select = c("id", "ST", "group_summary"))
setnames(meta, c("id", "group_summary"), c("geno_id", "env"))
meta[, env := trimws(tolower(env))]
al <- meta[al, on = .(geno_id)]

# alleles ranked by overall frequency, so colour tracks rank in every panel
al_rank <- al[, .N, by = allele][order(-N), allele]
al[, allele := factor(allele, levels = al_rank)]

# same palette and seed as allelic_diversity.R, so the panels read the same way
pal <- c(hcl.colors(40, "Batlow"), hcl.colors(40, "Mako"))
set.seed(seedy)
al_cols <- sample(pal, length(al_rank), replace = length(al_rank) > length(pal))
names(al_cols) <- al_rank

# 2. Nei's gene diversity, and the between-ST share of it -------------------

gene_div <- function(n) {
  n <- n[n > 0]; N <- sum(n)
  if (N < 2) return(NA_real_)
  (N^2 - sum(n^2)) / (N * (N - 1))
}

H_tot <- gene_div(al[, .N, by = allele]$N)
al[, n_st := .N, by = ST]
wst <- al[n_st >= 3, .(n = .N, H = gene_div(as.integer(table(allele)))), by = ST]
H_w <- wst[, sum(n * H) / sum(n)]

cat(sprintf("\n%s (%s): %d carriers, %d alleles, %d STs, %d sources\n",
            gene_lab, fam, nrow(al), uniqueN(al$allele),
            uniqueN(al$ST), uniqueN(al$env)))
cat(sprintf("overall H %.3f   mean within-ST H (n >= 3) %.3f   1 - Hw/Ht = %.3f\n",
            H_tot, H_w, 1 - H_w / H_tot))

# cross-check the per-genome alleles against the cached extraction
chk <- fread(paste0(outdir_dat, "/allele_counts_by_ST.csv"))[gene_family == fam]
if (nrow(chk)) {
  stc    <- setdiff(names(chk), c("gene_family", "allele_name"))
  cached <- melt(chk, id.vars = "allele_name", measure.vars = stc,
                 variable.name = "ST", value.name = "n")[n > 0]
  mine   <- al[, .N, by = .(allele_name = as.character(allele),
                            ST = factor(ST, levels = stc))]
  m <- merge(cached, mine, by = c("allele_name", "ST"), all = TRUE)
  setnafill(m, fill = 0L, cols = c("n", "N"))
  cat(sprintf("cross-check vs allele_counts_by_ST.csv: %s\n",
              if (all(m$n == m$N)) "identical"
              else paste(m[n != N, .N], "cells differ")))
}

# 3. the panel ---------------------------------------------------------------

# all but the ith piece replaced by a phantom of itself, so that coloured
# fragments of one plotmath string stay in register
# (as in 10_outlier_function/ags_alleles_by_ST_plot.R)
mask_parts <- function(parts, i) {
  out <- NULL
  for (j in seq_along(parts)) {
    pj  <- if (j == i) parts[[j]] else bquote(phantom(.(parts[[j]])))
    out <- if (is.null(out)) pj else bquote(.(out) * .(pj))
  }
  out
}

# One stacked bar per group. prop = TRUE gives the allele-frequency panel of
# allelic_diversity.R: bars scaled to 1, width proportional to sqrt(carriers)
# so the eye is not drawn to a group sampled once, and the n above each bar.
# prop = FALSE gives the same bars at their true height and equal width, which
# is where the sampling imbalance shows.
# lab_cols colours the group labels; pass src_cols for the source panels.
# fill: NULL stacks by allele in al_cols. Pass a named colour vector (such as
# src_cols) to stack by something else, in which case the rows of M are that
# something and the columns are the alleles.
plot_allele_grp <- function(M, prop = TRUE, ylab = NULL, hdr = NULL,
                            hdr_stats = TRUE, fill = NULL, lab_drop = 0,
                            lab_cols = NULL, lab_cex = 0.8, cex_hdr = 0.65,
                            lab_srt = 45) {

  if (is.null(fill)) M <- M[order(match(rownames(M), al_rank)), , drop = FALSE]
  np <- colSums(M); nb <- ncol(M); k <- sum(M)
  cols <- if (is.null(fill)) al_cols[rownames(M)] else fill[rownames(M)]

  if (prop) {
    wid <- sqrt(np) / max(sqrt(np))
    bp <- barplot(sweep(M, 2, np, "/"), col = cols, border = "white",
                  lwd = 0.4, width = wid, space = 0.35, ylim = c(0, 1),
                  yaxt = "n", names.arg = rep("", nb))
    axis(2, at = c(0, .5, 1), las = 1)
    text(bp, 1.02, np, cex = 0.65, col = "grey35", xpd = NA)
  } else {
    bp <- barplot(M, col = cols, border = "white", lwd = 0.4, space = 0.35,
                  ylim = c(0, max(np) * 1.08), las = 1,
                  names.arg = rep("", nb))
  }
  mtext(ylab %||% if (prop) "allele frequency" else "genomes", 2,
        line = 2.8, cex = 0.8)

  # Group labels take their group's colour, with a chip above each. Thorpe
  # et al.'s palette runs to cream, yellow and pale grey, which are unreadable
  # as text on white, so the text is darkened where the fill is light while
  # the chip keeps the published colour.
  darken <- function(x, floor_lum = 0.55) {
    r <- col2rgb(x) / 255
    lum <- colSums(r * c(0.2126, 0.7152, 0.0722))
    f <- pmin(1, floor_lum / pmax(lum, 1e-6))
    rgb(t(r * rep(f, each = 3)))
  }
  u  <- par("usr"); dy <- diff(u[3:4])
  y0 <- u[3] - dy * 0.055
  tc <- "black"
  if (!is.null(lab_cols)) {
    lc <- lab_cols[colnames(M)]; lc[is.na(lc)] <- "white"
    points(bp, rep(y0, length(bp)), pch = 22, bg = lc, col = "grey40",
           cex = 1.1, lwd = 0.5, xpd = NA)
    y0 <- y0 - dy * 0.045
    tc <- darken(lc)
  }
  # with sqrt widths a group of one is too thin to sit a label under, so
  # lab_drop suppresses the labels that would collide
  lab <- colnames(M); lab[np < lab_drop] <- ""
  text(bp, y0, lab, srt = lab_srt, xpd = NA, cex = lab_cex, col = tc,
       font = if (is.null(lab_cols)) 1 else 2,
       adj = if (lab_srt == 0) c(0.5, 1) else 1)

  # Header: name, then pangenome frequency, allele count and H, then the
  # frequency of the n_top major alleles in the colour of the allele each
  # refers to.
  if (!is.null(hdr)) {
    mtext(hdr, 3, line = 3.1, adj = 0, cex = cex_hdr + 0.03, font = 2)
  }
  if (!is.null(hdr) && hdr_stats) {
    nt  <- min(n_top, nrow(M))
    p_n <- sprintf("%1.2f", rowSums(M)[seq_len(nt)] / k)
    pre <- bquote("pan. freq =" ~ .(k) ~ "  |  " ~ .(nrow(M)) ~ "alleles  |  " ~
                    italic(H) ~ "=" ~ .(sprintf("%.2f", gene_div(rowSums(M)))) ~
                    "  |  ")
    parts <- c(list(pre), lapply(seq_len(nt), function(i)
      bquote(italic(p)[.(i)] ~ "=" ~ .(p_n[i]) *
               .(if (i < nt) ",  " else ""))))
    hc <- c("grey35", unname(cols[seq_len(nt)]))
    for (i in seq_along(parts))
      mtext(mask_parts(parts, i), 3, line = 1.4, adj = 0,
            cex = cex_hdr, col = hc[i])
  }
  invisible(bp)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# panel letter in user coordinates, so it survives any margin setting
panel_tag <- function(tag, dx = 0.07, dy = 0.28) {
  u <- par("usr")
  text(u[1] - diff(u[1:2]) * dx, u[4] + diff(u[3:4]) * dy,
       tag, font = 2, cex = 1.15, xpd = NA)
}

# allele x ST, rare STs pooled as the cached-table version does
st_keep <- al[, .N, by = ST][N >= min_n][order(-N), ST]
al[, st_g := factor(fifelse(ST %chin% st_keep, ST, "other"),
                    levels = c(st_keep, "other"))]
M_st <- as.matrix(dcast(al[, .N, by = .(allele, st_g)], allele ~ st_g,
                        value.var = "N", fill = 0, drop = FALSE),
                  rownames = "allele")
M_st <- M_st[rowSums(M_st) > 0, colSums(M_st) > 0, drop = FALSE]

# 4. the ST figure on its own ------------------------------------------------

if (make_fig) {
  png(paste0(outdir_fig, "/alleles_by_ST_", fam, ".png"),
      width = 18, height = 14, units = "cm", res = 600,
      pointsize = 11, type = "cairo")
  par(mfrow = c(2, 1), mar = c(4.6, 4.8, 4.8, 1.4), bty = "l")

  plot_allele_grp(M_st, prop = FALSE, hdr = paste0(gene_lab, ", by ST"))
  panel_tag("a")
  plot_allele_grp(M_st, prop = TRUE)
  panel_tag("b")

  dev.off()
  cat("wrote", paste0(outdir_fig, "/alleles_by_ST_", fam, ".png"), "\n")
}
