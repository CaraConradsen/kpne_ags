# Stepped balancing selection: do wzm and wbbM protein variants track O type?
#
# wzm (g004051) and wbbM (g000322) sit in the same syntenic block, MSU_2,
# between the anchors A308:A309, and are carried by very nearly the same set of
# genomes (185 and 182 of 260). They are therefore present or absent together,
# yet their protein allelic diversity could hardly differ more: H = 0.13 over 8
# alleles for wzm against H = 0.89 over 50 for wbbM. Both are strongly purifying
# at the sequence level (piN/piS = 0.018 and 0.111), so the contrast is not one
# of diversifying against conservative selection.
#
# The remaining explanation is that wbbM, the galactan I galactosyltransferase,
# carries the residues that distinguish the O-antigen subtypes, while wzm, the
# ABC transporter permease that exports whatever galactan is built, does not.
# This script tests that directly: for every variable amino acid position in
# each protein it asks whether the minor residue is distributed across the O
# types at random, and draws the answer as a lolliplot in which each lollipop is
# a pie of the O types that carry the minor residue.
#
# If variants track O type the pies are close to single-coloured. If they do not
# they reproduce the overall composition of the carriers, shown in the legend.

# Parameters ---------------------------------------------------------------

min_minor  <- 1L           # a residue is drawn if this many genomes carry it
min_report <- 2L           # ... but only labelled, and tabulated, from this many
q_cut     <- 0.05          # BH threshold for the per-site tests
B_sim     <- 2e4           # replicates for the simulated Fisher p-values
anchor_keep <- "A308:A309" # the syntenic interval the two genes share
msu_keep    <- "MSU_2"

focal_dt <- data.table(
  gene_family = c("g004051", "g000322"),
  short       = c("wzm", "wbbM"),
  descr       = c("O-antigen ABC transporter permease",
                  "galactan I galactosyltransferase"))

# the four galactan subtypes plus the single ambiguous call. O1 and O2 share the
# galactan I cluster and differ by wbbY/wbbZ elsewhere on the chromosome, so the
# v1/v2 split is the one that belongs to this locus: blues for O1, oranges for O2
sero_col <- c("O1v1" = "#2C6E8F", "O1v2" = "#8FC1D4",
              "O2v1" = "#C4631E", "O2v2" = "#F3BE8A",
              "O1/O2v1" = "grey55")

lab_cex <- 0.55

# lollipopPlot is the current name for trackViewer::lolliplot; older installs
# export only the latter
lollipop <- if (exists("lollipopPlot")) lollipopPlot else lolliplot

# Data ---------------------------------------------------------------------

meta <- fread(paste0(outdir_dat, "/olocus.csv"),
              select = c("geno_id", "ST", "O_locus"))

# confirm both families really do sit in the one block before going further
msu <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
             select = c("gene_family", "geno_id", "anchor", "msu"))
msu <- msu[gene_family %in% focal_dt$gene_family]
print(msu[, .N, by = .(gene_family, anchor, msu)])
stopifnot(all(msu$anchor == anchor_keep), all(msu$msu == msu_keep))

# the peptide matrix, one row per genome. The alignments are the nucleotide
# alignments written by 06_piS_estimates_seg_sites/obs_ag_vals.R, so they are in
# frame and gapped by codon; a codon holding a gap translates to "-" rather than
# to an amino acid, which keeps indels visible without inventing residues
gc_tab <- Biostrings::GENETIC_CODE

pep_mat <- function(fam) {
  dna <- readDNAStringSet(paste0(outdir_dat, "/ag_fus_alns/", fam, ".fasta"))
  m <- do.call(rbind, lapply(as.character(dna), function(x) {
    cd <- substring(x, seq(1, nchar(x) - 2, 3), seq(3, nchar(x), 3))
    a  <- gc_tab[cd]
    a[is.na(a)] <- "-"      # gap, or any codon holding an ambiguity code
    a
  }))
  rownames(m) <- sub("_[0-9]+$", "", names(dna))     # fus locus tag -> geno_id
  if (all(m[, ncol(m)] %in% c("*", "-"))) m <- m[, -ncol(m), drop = FALSE]
  m
}

# Per-site variation -------------------------------------------------------
# One row per non-majority residue, so a site holding three residues gives two
# rows. Gaps and stops are not residues and are excluded from both the counts
# and the test, rather than being treated as a state of their own.

site_variants <- function(m, sero) {
  lv <- names(sero_col)[names(sero_col) %in% sero]
  rbindlist(lapply(seq_len(ncol(m)), function(j) {
    z  <- m[, j]
    ok <- !z %in% c("-", "*")
    tb <- sort(table(z[ok]), decreasing = TRUE)
    if (length(tb) < 2) return(NULL)
    rbindlist(lapply(names(tb)[-1], function(r) {
      car <- ok & z == r
      ct  <- table(factor(sero[car], levels = lv))
      # 2 x n_sero: carrying the minor residue against not, by O type
      tj  <- table(sero[ok], z[ok] == r)
      data.table(site = j, major = names(tb)[1], minor = r,
                 variant = paste0(names(tb)[1], j, r),
                 n_minor = sum(car),
                 P = fisher.test(tj, simulate.p.value = TRUE, B = B_sim)$p.value,
                 as.data.table(as.list(ct)))
    }))
  }))
}

# Allelic diversity, within and between O types ----------------------------
# Nei's unbiased H over whole peptides, partitioned the usual way: Gst = 1 - Hw/Ht
# says how much of the allelic diversity sits between the O types rather than
# within them. This is the same statistic as the per-site tests, read over the
# whole protein instead of one position at a time.

h_nei <- function(x) {
  n <- length(x)
  p <- prop.table(table(x))
  (n / (n - 1)) * (1 - sum(p^2))
}

part_h <- function(pep, sero) {
  w <- data.table(sero = sero, pep = pep)[, .(h = h_nei(pep), n = .N), by = sero][n > 1]
  Ht <- h_nei(pep)
  Hw <- sum(w$h * w$n) / sum(w$n)
  list(Ht = Ht, Hw = Hw, Gst = 1 - Hw / Ht, by_sero = w)
}

# Run ----------------------------------------------------------------------

res <- vector("list", nrow(focal_dt))
names(res) <- focal_dt$short

for (i in seq_len(nrow(focal_dt))) {

  fam  <- focal_dt$gene_family[i]
  m    <- pep_mat(fam)
  sero <- meta$O_locus[match(rownames(m), meta$geno_id)]
  pep  <- apply(m, 1, paste, collapse = "")

  sv <- site_variants(m, sero)
  sv[, q := p.adjust(P, "BH")]
  setorder(sv, -n_minor)

  ph <- part_h(pep, sero)
  al <- table(sero, factor(as.integer(factor(pep,
         levels = names(sort(table(pep), decreasing = TRUE))))))

  cat(sprintf("\n== %s (%s), %s: n = %d genomes, L = %d aa, %d alleles\n",
              focal_dt$short[i], fam, focal_dt$descr[i],
              nrow(m), ncol(m), uniqueN(pep)))
  cat(sprintf("   variable positions %d, of which %d have a minor residue in >= %d genomes\n",
              uniqueN(sv$site), uniqueN(sv[n_minor >= min_report]$site), min_report))
  cat(sprintf("   Ht = %.3f, Hw = %.3f, Gst = %.3f\n", ph$Ht, ph$Hw, ph$Gst))
  cat(sprintf("   variants whose frequency differs by O type (BH q < %.2f): %d of %d\n",
              q_cut, sv[n_minor >= min_report & q < q_cut, .N],
              sv[n_minor >= min_report, .N]))
  print(ph$by_sero[order(-n)])
  print(sv[n_minor >= min_report])

  res[[i]] <- list(m = m, sero = sero, pep = pep, sv = sv, ph = ph, al = al)

  fwrite(sv, paste0(outdir_dat, "/", focal_dt$short[i], "_site_variants_by_O.csv"))
  fwrite(as.data.table(al), paste0(outdir_dat, "/", focal_dt$short[i], "_alleles_by_O.csv"))
}

# Figure -------------------------------------------------------------------
# One panel per O type per gene, one figure per gene. The variant set is the
# same in every panel of a figure and the x axis is the same protein, so the
# panels read as a column: a position that is diagnostic of an O type is a full
# pie in that panel and an empty one in the others.
#
# Each lollipop is a non-majority residue; the pie is the frequency of that
# residue within the panel's O type, filled in the O type's colour against
# grey for the genomes of that O type carrying the majority residue. A black
# border marks a variant whose frequency differs by O type over all carriers
# (BH q < q_cut), tested once on the whole set rather than per panel.
#
# trackViewer draws in grid, so the panels go in a grid layout rather than
# through layout(), and newpage = FALSE keeps them all on the one page.

min_sero <- 5L     # an O type gets a panel only with this many carriers

lolli_gr <- function(sv, ot, n_ot, show_lab) {

  setorder(sv, site)
  n_in <- sv[[ot]]                       # carriers of the minor residue in this O type

  gr <- GRanges("prot", IRanges(sv$site, width = 1, names = sv$variant))
  # for type = "pie" the slices are the value* columns: here, has the minor
  # residue against has the majority one, within this O type only
  gr$value1 <- n_in
  gr$value2 <- n_ot - n_in
  gr$color  <- rep(list(c(unname(sero_col[ot]), "grey90")), length(gr))
  gr$border <- ifelse(sv$q < q_cut, "black", "grey75")
  # pie area, not radius, tracks how many carriers hold the residue overall, so
  # the same variant is the same size in every panel
  gr$cex <- 0.55 + 0.95 * sqrt(sv$n_minor) / sqrt(max(sv$n_minor))
  gr$label.parameter.rot <- 45
  # names only on the top panel, and only for the shared variants, or the
  # labels are repeated four times over and bury the figure
  gr$label.parameter.label <- if (show_lab)
    ifelse(sv$n_minor >= min_report, sv$variant, NA) else NA
  gr
}

feat_gr <- function(L, nm, fill) {
  f <- GRanges("prot", IRanges(1, width = L, names = nm))
  f$fill   <- fill
  f$height <- 0.035
  f
}

for (i in seq_len(nrow(focal_dt))) {

  r  <- res[[i]]
  L  <- ncol(r$m)
  sv <- r$sv[n_minor >= min_minor]

  # the O types with enough carriers to be worth a panel, in sero_col order
  cmp <- table(factor(r$sero, levels = names(sero_col)))
  ots <- names(cmp)[cmp >= min_sero]
  cat(sprintf("%s: panels for %s; dropped %s (fewer than %d carriers)\n",
              focal_dt$short[i], paste(ots, collapse = ", "),
              paste(names(cmp)[cmp > 0 & cmp < min_sero], collapse = ", "), min_sero))

  png(paste0(outdir_fig, "/o_antigen_variants_", focal_dt$short[i], ".png"),
      width = 30, height = 4.2 * length(ots) + 2.5, units = "cm", res = 600,
      pointsize = 12, type = "cairo")

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(length(ots), 1)))

  for (j in seq_along(ots)) {

    gr <- lolli_gr(sv, ots[j], cmp[[ots[j]]], show_lab = j == 1)
    ft <- feat_gr(L, paste0(focal_dt$short[i], " (", L, " aa)"),
                  if (i == 1) "#BFD9E6" else "#E8D2B8")

    # diversity of this O type on its own, for the panel header
    in_ot <- r$sero == ots[j]
    h_ot  <- h_nei(r$pep[in_ot])

    pushViewport(viewport(layout.pos.row = j, layout.pos.col = 1))
    lollipop(gr, ft, ranges = GRanges("prot", IRanges(1, L)),
             type = "pie", newpage = FALSE, cex = lab_cex + 0.35,
             ylab = FALSE, yaxis = FALSE, xaxis = j == length(ots))
    grid.text(letters[j], x = unit(0.5, "lines"),
              y = unit(1, "npc") - unit(0.8, "lines"),
              gp = gpar(fontface = 2, cex = 1.1))
    grid.text(sprintf("%s   n = %d, %d alleles, H = %.2f",
                      ots[j], sum(in_ot), uniqueN(r$pep[in_ot]), h_ot),
              x = unit(2, "lines"), y = unit(1, "npc") - unit(0.8, "lines"),
              just = "left", gp = gpar(cex = 0.8))
    popViewport()
  }

  popViewport()
  dev.off()
}
