# Allelic diversity of outlier accessory genes (manuscript figure)
#
# Synonymous diversity measures the depth of a gene's genealogy; it says
# nothing about whether that depth is expressed as protein variation. A gene
# can carry an excess of synonymous variation while remaining effectively
# monomorphic in peptide, or while segregating many peptides at intermediate
# frequency. Separating the two is what distinguishes alleles maintained by
# selection from a locus that has merely inherited the divergence of the
# lineages it is confined to.
#
# Allelic diversity is Nei's gene diversity, the probability that two
# randomly drawn copies of a gene carry different peptides,
#     H = (n / (n - 1)) * (1 - sum(p_i^2))
# with p_i the frequency of peptide allele i among the n carriers
# (Nei 1987, eq. 8.4). Peptides are the right unit for a functional or
# antigenic hypothesis and avoid counting synonymous variation twice.
#
# Allele counts per ST are read from allele_counts_by_ST.csv
# (ags_alleles_by_ST.R); the per-gene panels are drawn by
# ags_alleles_by_ST_plot.R.
#
# (a)   Nei's H against synonymous diversity, every outlier AG
# (b-f) allele frequencies by ST for abiEi, chrA, wzm, wbbM and tdh

# Parameters --------------------------------------------------------------

q_cut  <- 0.05
min_n  <- 2      # STs carrying fewer copies than this are pooled as "other"
seedy  <- 21     # allele palette seed

# Read the cached allele table fresh, so that a subset left in the session
# by an earlier run cannot be picked up silently.
allele_dt <- fread(paste0(outdir_dat, "/allele_counts_by_ST.csv"))

# gene_div(), collapse_st() and plot_allele_st()
source("./analysis/10_outlier_function/ags_alleles_by_ST_plot.R")

# Focal AGs: robust outliers (see AG outlier candidates doc, 2026-10-14).
# Labels use corrected gene names, not PIRATE consensus names.
focal_dt <- data.table(
  gene_family = c("g004051", "g003783", "g001147", "g000322",
                  "g003390", "g003914",
                  "g003299", "g001254",
                  "g002556", "g001058",
                  "g002655", "g003239",
                  "g002664", "g002490",
                  "g003377", "g004320",
                  "g002126_1", "g015591", "g006452", "g003974"),
  short = c("wzm", "wzt", "glf", "wbbM", "vgrG", "tssA", "abiEi", "abiEii",
            "asr", "chrA", "acoA", "acoB", "phyH", "phyH2", "sorF", "sorC",
            "glpT", "rplR", "ttdA", "tdh"),
  label = c("wzm, O-antigen export permease",
            "wzt, O-antigen export ATP-binding protein",
            "glf, UDP-galactopyranose mutase (O-antigen)",
            "wbbM, O-antigen glycosyltransferase (DUF4422)",
            "vgrG, T6SS tip protein",
            "tssA, T6SS assembly protein",
            "abiEi, Phage defence antitoxin",
            "abiEii, Phage defence toxin",
            "asr, Acid resistance protein",
            "chrA, Chromate transporter",
            "acoA, Acetoin dehydrogenase E1 alpha",
            "acoB, Acetoin dehydrogenase E1 beta",
            "phyH, Ectoine dioxygenase",
            "phyH, Phytanoyl-CoA dioxygenase",
            "sorF, Sorbose PTS EIIA",
            "sorC, Sorbitol operon regulator",
            "glpT, Glycerol-3-phosphate transporter (fosfomycin uptake)",
            "rplR, 50S ribosomal protein L18",
            "ttdA, Fumarate hydratase",
            "tdh, NAD(P)-dependent alcohol dehydrogenase"))

# genes shown as allele-frequency panels, in plotting order
panel_genes <- c(  "g004051", "g000322","g001058", "g003299", "g003974")

# Outlier calls -----------------------------------------------------------

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites", "syn_sites",
                               "expected_S", "p_bh"))

gene_pvals[, s_per_site := seg_sites / syn_sites]
gene_pvals[, call := fcase(p_bh < q_cut & seg_sites > expected_S, "high",
                           p_bh < q_cut & seg_sites < expected_S, "low",
                           default = "ns")]

# Allelic diversity -------------------------------------------------------
# Column sums of allele_counts_by_ST.csv reproduce the pangenome frequency
# of every set 1 AG exactly, so the alleles and the test are counted over
# the same carriers.

st_cols <- grep("^ST", names(allele_dt), value = TRUE)

allele_counts <- allele_dt[, .(gene_family, allele_name,
                               n_all = rowSums(.SD)), .SDcols = st_cols]

allele_stats <- allele_counts[n_all > 0,
                              .(k         = sum(n_all),
                                n_alleles = .N,
                                H         = gene_div(n_all)),
                              by = gene_family]

allele_stats <- gene_pvals[allele_stats, on = .(gene_family)]
allele_stats[focal_dt, on = .(gene_family), `:=`(short = i.short, label = i.label)]

stopifnot(allele_stats[, all(k == freq)],
          nrow(allele_stats) == uniqueN(gene_pvals$gene_family))

setorder(allele_stats, -H)
allele_stats[!is.na(short)]

allele_stats[, is_cloud := fcase(freq < (0.15*260), "y",
                                  default = "n")]

fwrite(allele_stats, paste0(outdir_dat, "/ag_allele_stats.csv"))

with(allele_stats[call == "ns" & is_cloud !="y"],
     suppressWarnings(cor.test(s_per_site, H, method = "spearman")))


# Panel helper ------------------------------------------------------------

# panel letter in user coordinates, so it survives any margin setting
panel_tag <- function(tag, dx = 0.11, dy = 0.12) {
  u <- par("usr")
  text(u[1] - diff(u[1:2]) * dx, u[4] + diff(u[3:4]) * dy,
       tag, font = 2, cex = 1.15, xpd = NA)
}

# Figure ------------------------------------------------------------------

png(paste0(outdir_fig, "/allelic_diversity.png"),
    width = 18, height = 28, units = "cm", res = 600,
    pointsize = 11, type = "cairo")

layout(matrix( c(1,2,3,rep(4:8, each = 3)), byrow = TRUE, ncol = 3),
       widths = c(1.7,1.5,0.5),
       heights = c(1.55, 1, 1, 1, 1, 1))

## (a) all outliers, highlighted genes labelled ---------------------------
# H is estimated from k carriers, so a gene held by two genomes carries far
# less information than one held by 200. Point area tracks pangenome
# frequency on a sqrt scale, as the bar widths in b-f do.


call_col <- c(high = "#3FA9F5", low = "#F26C2A")

all_ags <- allele_stats[call != "ns" & !is.na(H)]

cex_k <- function(k, lo = 0.30, hi = 1.45) {
  lo + (hi - lo) * (sqrt(k) - sqrt(min(all_ags$k))) /
    (sqrt(max(all_ags$k)) - sqrt(min(all_ags$k)))
}

for (cld in c("y", "n")) {
  if (cld == "y"){
    par(mar = c(4.2, 4, 2.6, 0), bty = "o")
  } else{
    par(mar = c(4.2, 0, 2.6, 0), bty = "o")
  }
  
  d  <- allele_stats[call != "ns" & is_cloud == cld]
  hl <- d[!is.na(short)]
  
  with(d,
       plot(s_per_site, H, type = "n", las = 1,
            xlim = c(0, max(all_ags$s_per_site) * 1.05),
            ylim = c(0, 1.04),
            xlab = "", yaxt = "n",
            ylab = ""))
  
  # largest first, so that common genes do not hide behind rare ones
  setorder(d, -k)
  
  with(d,
       points(s_per_site, H, pch = 16, cex = cex_k(k),
              col = adjustcolor(call_col[call], 0.30)))
  
  if (cld == "y") panel_tag("a", dx = 0.15, dy = 0.07)
  
  if (nrow(hl) > 0){
    with(hl,
         points(s_per_site, H, pch = 21, cex = cex_k(k) + 0.35,
                bg = adjustcolor(call_col[call], 0.95), 
                col = "grey15", lwd = 0.7))
    # labels of crowded points alternate above and below
    xs <- hl$s_per_site / diff(range(d$s_per_site))
    ys <- hl$H / diff(range(d$H))
    dm <- as.matrix(dist(cbind(xs, ys))); diag(dm) <- Inf
    crowded <- apply(dm, 1, min) < 0.10
    lab_pos <- rep(4, nrow(hl))
    lab_pos[crowded] <- rep(c(3, 1), length.out = sum(crowded))
    lab_pos[lab_pos == 1 & hl$H < 0.08] <- 3    # keep labels off the axis
    text(hl$s_per_site, hl$H, hl$short, pos = lab_pos, offset = 0.45,
         cex = 0.6, font = 3, xpd = NA)
    axis(side = 1, at = 0-0.1, labels = "Synonymous segregating sites per site",
         tick = FALSE, line = 1)
    
  } else{
    axis(side = 2, at = 0.5, labels = expression("Protein allelic diversity (" * italic(H) * ")"),
         tick = FALSE, line = 2)
    axis(side = 2, at = seq(0, 1.04, 0.2), labels = seq(0, 1.04, 0.2), las = 2)
    
    
  }
  
  if (cld == "y") {
    legend("top", inset = c(0, -0.12), xpd = TRUE,
           legend = expression(bold("Cloud AGs, ") * bolditalic(k) < bold("15%")),
           bty = "n", col = "grey40")
  }else{
    legend("top", inset = c(0, -0.12), xpd = TRUE,
           legend = expression(bold("Shell & Soft-core,") ~ bolditalic(k) >= bold("15%")),
           bty = "n", col = "grey40")
  }
  
  rho <- suppressWarnings(cor.test(d$s_per_site, d$H, method = "spearman"))
  
  p_lab <- if (rho$p.value < 0.001) "< 0.001" else sprintf("= %.3f", rho$p.value)
  
  legend("topright", xpd = NA, bty = "n", cex = 0.66,
    legend = as.expression(bquote(rho ~ "=" ~ .(sprintf("%.2f", rho$estimate)) *
                                    "," ~ italic(P) ~ .(p_lab))))
  
}

# Legends sit together on the right, where the cloud thins out: the call
# above, then the frequency scale, then the rank correlation.
par(mar = c(0, 0, 0, 0))
plot.new()

k_key <- c(10, 50, 150, 250)
k_key <- k_key[k_key <= max(d$k)]

# anchored off par("usr") rather than fixed user coordinates, so the block
# stays in the right margin whatever the x variable is
legend("topleft", inset = c(0.05, 0.15), xpd = TRUE,
       bty = "n", cex = 0.66, pch = c(16, 16, 21),
       pt.cex = c(0.9, 0.9, 1.2), y.intersp = 1.25,
       col = c(adjustcolor(call_col, 0.55), "grey15"),
       pt.bg = c(NA, NA, "white"),
       legend = c(sprintf("Excess diversity (n = %d)", sum(d$call == "high")),
                  sprintf("Deficit diversity (n = %d)", sum(d$call == "low")),
                  "Highlighted genes"))

legend("left",inset = c(0.1, 0.12), xpd = TRUE,
       bty = "n", cex = 0.66, pch = 16, col = "grey70",
       pt.cex = cex_k(k_key), y.intersp = 1.45,
       title = "Pangenome frequency", title.adj = 0, title.cex = 0.66,
       legend = k_key)






## (b-f) allele frequencies by ST ----------------------------------------

for (i in seq_along(panel_genes)) {

  par(mar = c(3.9, 4.8, 4.8, 1.4))

  plot_allele_st(panel_genes[i],
                 gname = focal_dt[gene_family == panel_genes[i], label],
                 min_n = min_n, seedy = seedy)

  panel_tag(letters[i + 1], dx = 0.075, dy = 0.85)
}

dev.off()
