
# Allele frequencies by ST, for any accessory gene
#
# Sourced for its functions: it fixes no gene and draws nothing on source.
#   gene_div()       Nei's gene diversity from a vector of allele counts
#   collapse_st()    allele x ST counts for one gene, rare STs pooled
#   plot_allele_st() stacked barplot of allele frequency within each ST
#
# Counts are read from allele_counts_by_ST.csv (ags_alleles_by_ST.R), which
# translates every set 1 AG and tabulates peptide alleles against ST.

if (!exists("allele_dt")) {
  allele_dt <- fread(paste0(outdir_dat, "/allele_counts_by_ST.csv"))
}

# Nei's gene diversity ----------------------------------------------------
# H = (n / (n - 1)) * (1 - sum(p_i^2)), in count form (Nei 1987, eq. 8.4).
# The probability that two copies drawn without replacement differ.

gene_div <- function(n) {
  n <- n[n > 0]
  N <- sum(n)
  if (N < 2) return(NA_real_)
  (N^2 - sum(n^2)) / (N * (N - 1))
}

# Allele x ST counts ------------------------------------------------------
# min_n = 2 pools only the STs holding a single carrier of the gene.

collapse_st <- function(g_fam, min_n = 2, merge_lv = FALSE) {

  tab <- allele_dt[gene_family == g_fam]
  if (!nrow(tab)) stop("no allele counts for ", g_fam)

  # named explicitly, so that any column a caller has added to allele_dt is
  # not melted as though it were an ST
  st_v <- grep("^ST", names(tab), value = TRUE)

  long <- melt(tab, id.vars = "allele_name", measure.vars = st_v,
               variable.name = "ST", value.name = "n")

  long[, ST := as.character(ST)]
  if (merge_lv) long[, ST := sub("-[0-9]+LV$", "", ST)]   # merge single-locus variants

  long <- long[, .(n = sum(n)), by = .(allele_name, ST)][n > 0]

  tot  <- long[, .(N = sum(n)), by = ST][order(-N)]
  keep <- tot[N >= min_n, ST]

  long[, ST_grp := fifelse(ST %chin% keep, ST, "other")]

  out <- long[, .(n = sum(n)), by = .(allele_name, ST_grp)]

  lev <- c(tot[ST %chin% keep][order(-N), ST],
           if (any(long$ST_grp == "other")) "other")

  out[, ST_grp := factor(ST_grp, levels = lev)]

  setattr(out, "n_other_st", tot[!ST %chin% keep, .N])
  out[]
}

# Plot --------------------------------------------------------------------
# Bar width is proportional to sqrt(carriers), so the eye is not drawn to
# the allele composition of an ST sampled once or twice.

# all but the ith piece replaced by a phantom of itself, so that coloured
# fragments of one plotmath string stay in register
mask_parts <- function(parts, i) {
  out <- NULL
  for (j in seq_along(parts)) {
    pj  <- if (j == i) parts[[j]] else bquote(phantom(.(parts[[j]])))
    out <- if (is.null(out)) pj else bquote(.(out) * .(pj))
  }
  out
}

plot_allele_st <- function(g_fam, gname = g_fam, min_n = 2, seedy = 21,
                           n_top = 3, cex_hdr = 0.65) {

  long <- collapse_st(g_fam, min_n = min_n)

  M <- as.matrix(dcast(long, allele_name ~ ST_grp, value.var = "n", fill = 0),
                 rownames = "allele_name")

  # alleles ordered by overall frequency, so that colour tracks rank
  M <- M[order(rowSums(M), decreasing = TRUE), , drop = FALSE]

  n_all    <- rowSums(M)
  n_per_st <- colSums(M)
  nb       <- length(n_per_st)
  k        <- sum(M)

  Mp <- sweep(M, 2, n_per_st, "/")

  pal <- c(hcl.colors(40, "Batlow"),
           hcl.colors(40, "Mako"))

  set.seed(seedy)
  cols <- sample(pal, nrow(M), replace = nrow(M) > length(pal))

  wid <- sqrt(n_per_st) / max(sqrt(n_per_st))

  if (nb < 6) {
    bp <- barplot(Mp, col = cols, border = "white", lwd = 0.4,
                  width = wid, space = 0.35,
                  xlim = c(0, (sum(wid) + nb * 0.35 * mean(wid)) * 2),
                  ylim = c(0, 1), yaxt = "n", names.arg = rep("", nb))
  } else {
    bp <- barplot(Mp, col = cols, border = "white", width = wid,
                  space = 0.35, ylim = c(0, 1), yaxt = "n",
                  names.arg = rep("", nb))
  }

  axis(2, at = c(0, .5, 1), las = 1)
  mtext("allele frequency", 2, line = 2.8, cex = 0.8)
  text(bp, -0.06, colnames(Mp), srt = 45, adj = 1, xpd = NA, cex = 0.8)
  text(bp, 1.02, n_per_st, cex = 0.65, col = "grey35", xpd = NA)

  # Header: name, then pangenome frequency, allele count and H, then the
  # frequency of the n_top major alleles in the colour of the allele each
  # refers to.
  mtext(gname, 3, line = 3.1, adj = 0, cex = cex_hdr + 0.03, font = 2)

  n_top <- min(n_top, nrow(M))
  p_n   <- sprintf("%1.2f", n_all[seq_len(n_top)] / k)

  pre <- bquote("pan. freq =" ~ .(k) ~ "  |  " ~ .(nrow(M)) ~ "alleles  |  " ~
                  italic(H) ~ "=" ~ .(sprintf("%.2f", gene_div(n_all))) ~ "  |  ")

  parts <- c(list(pre),
             lapply(seq_len(n_top), function(i)
               bquote(italic(p)[.(i)] ~ "=" ~ .(p_n[i]) *
                        .(if (i < n_top) ",  " else ""))))

  cols_hdr <- c("grey35", cols[seq_len(n_top)])

  for (i in seq_along(parts)) {
    mtext(mask_parts(parts, i), 3, line = 1.4, adj = 0,
          cex = cex_hdr, col = cols_hdr[i])
  }

  invisible(list(bp = bp, M = M, cols = cols, H = gene_div(n_all)))
}
