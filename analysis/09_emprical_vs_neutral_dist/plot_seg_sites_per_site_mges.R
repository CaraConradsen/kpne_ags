
# Two-panel version of plot_seg_sites_per_site.R
# Panel a: all syntenic AGs (as before)
# Panel b: same, with MGE-associated AGs removed
# BH correction and Spearman's rho are recomputed within each panel, so panel b
# is a standalone test on the reduced set rather than a subset of panel a's calls.

alpha_fdr = 0.05

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

# Identify MGEs -----------------------------------------------------------
# Same COG source as plot_outliers_by_COG.R: egng_cog, renamed COG_funct_cat.

anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                     select = c("gene_family", "egng_cog",
                                "consensus_gene_name", "consensus_product")))

setnames(anno, old = "egng_cog", new = "COG_funct_cat")

# collapse to one row per family; flag any family carrying conflicting annotation
dup_fams <- anno[, .N, by = gene_family][N > 1, gene_family]
if (length(dup_fams))
  message("Families with >1 annotation row, first row taken: ", length(dup_fams))
anno <- unique(anno, by = "gene_family")

# COG X = Mobilome (prophages, transposons). COG strings can hold several
# letters, so search for X rather than matching exactly.
anno[, is_mge := grepl("X", COG_funct_cat, fixed = TRUE)]

# eggNOG assigns X conservatively, so also screen the product text.
# Set to FALSE to use the COG call alone.
use_product_screen <- TRUE

mge_regex <- paste0("transposase|transposon|integrase|excisionase|resolvase|",
                    "recombinase|insertion sequence|IS[0-9]|phage|prophage|",
                    "conjugal|conjugative|relaxase|mobilization|mobile element|",
                    "plasmid|replicase|type IV secretion")

# phage-defence genes are not themselves MGEs; keep them
mge_keep <- "anti-phage|phage resistance|phage defen|abortive infection|restriction"

if (use_product_screen) {
  anno[grepl(mge_regex, consensus_product, ignore.case = TRUE), is_mge := TRUE]
  anno[grepl(mge_keep, consensus_product, ignore.case = TRUE), is_mge := FALSE]
}

gene_pvals <- merge(gene_pvals, anno, by = "gene_family", all.x = TRUE)

# unannotated families are retained in panel b; note this in the legend
gene_pvals[is.na(is_mge), is_mge := FALSE]

cat("AGs tested:", nrow(gene_pvals),
    "| flagged MGE:", gene_pvals[is_mge == TRUE, .N],
    "| retained in panel b:", gene_pvals[is_mge == FALSE, .N], "\n")

# Recompute calls within a set --------------------------------------------
# p_two is unchanged; only the BH correction depends on the set being tested.

prep_calls <- function(dt) {
  dt <- copy(dt)
  dt[, p_bh := p.adjust(p_two, method = "BH")]
  dt[, s_per_site  := seg_sites / syn_sites]
  dt[, exp_per_site := expected_S / syn_sites]
  dt[, call := "ns"]
  dt[p_bh < alpha_fdr & seg_sites > expected_S, call := "high"]
  dt[p_bh < alpha_fdr & seg_sites < expected_S, call := "low"]
  dt[, call := factor(call, levels = c("ns", "high", "low"))]
  dt[]
}

dt_a <- prep_calls(gene_pvals)
dt_b <- prep_calls(gene_pvals[is_mge == FALSE])

cols <- c(ns = rgb(0.1, 0.1, 0.1, alpha = 0.5), high = "#3FA9F5", low = "#F26C2A")

# Genes of interest -------------------------------------------------------

trouble_children = data.table(gene_family = c("g000145", "g002611", "g004058",#Outer membrane
                                              "g000322", "g001147", "g004051", "g003783", "g000817",#5 co-inherited, O-antigen is textbook NFDS
                                              "g002097",#interbacterial competition
                                              "g000382","g000849_1", # Efflux and envelope stress, adaptive variation in efflux capacity
                                              "g003299"
                                              ),
                              consensus_gene_name = c("OmpA", "OmpK36","ompR",
                                                      "wbbM1", "rfbD", "rfbA", "rfbB", "glycosyltransferase",
                                                      "impA",
                                                      "acrA","baeS",
                                                      "AbiEi"
                                                      ),
                              dx = c(-20, -20, -15,
                                     15, -20, -20,  15, -20,
                                     15,
                                     -15, -12,
                                     -15),
                              dy = c( 0.000,  0.005, -0.010,
                                      -0.005,  0.010,  0.010,  0.015, -0.015,
                                      0.010,
                                      -0.005,  0.020,
                                      0.005))

# Panel drawing -----------------------------------------------------------
# Panels share both axes; only panel a carries tick labels and the axis title.
# The legend is drawn once in the outer right margin, so no panel draws its own.

draw_panel <- function(dt, xlim, ylim, panel_lab,
                       show_y = TRUE) {

  pt_col <- adjustcolor(cols[as.character(dt$call)], alpha.f = 0.6)

  plot(dt$freq, dt$s_per_site, type = "n",
       xlab = "",
       ylab = "", xaxt = "n", yaxt = "n",
       xlim = xlim, ylim = ylim)

  axis(side = 1, at = pretty(xlim), labels = pretty(xlim))

  if (show_y) {
    axis(side = 2, at = pretty(ylim), labels = pretty(ylim), las = 2)
    axis(side = 2, at = mean(ylim),
         labels = "",
         line = 2.5, tick = FALSE)
  } else {
    axis(side = 2, at = pretty(ylim), labels = FALSE)
  }

  points(dt$freq, dt$s_per_site, pch = 16, col = pt_col)

  # outlier count for this panel
  n_out = nrow(dt[call != "ns", ])
  perc_AGs_out = (n_out / nrow(dt)) * 100

  mtext(bquote(.(n_out)*" ("*.(sprintf("%2.2f", perc_AGs_out))*"%) of "*
                 .(nrow(dt))*" AGs"),
        side = 3, line = -1, cex = 0.7)

  spearmans <- with(dt, suppressWarnings(cor.test(freq, s_per_site,
                                                  method = "spearman")))

  rho_val <- sprintf("%.2f", spearmans$estimate)

  p_part <- if (spearmans$p.value < 0.001) {
    quote(italic(P) < 0.001)
  } else {
    bquote(italic(P) == .(sprintf("%.3f", spearmans$p.value)))
  }

  mtext(bquote(rho == .(rho_val) * ";" ~ .(p_part)),
        side = 3, adj = 0.5, line = -2.0, cex = 0.7)

  # ring and label the genes of interest
  tc <- merge(trouble_children,
              dt[, .(gene_family, s_per_site, freq)],
              by = "gene_family", all.x = TRUE)

  if (tc[is.na(s_per_site), .N] > 0)
    message("Panel ", panel_lab, " not plotted: ",
            paste(tc[is.na(s_per_site), consensus_gene_name], collapse = ", "))

  tc <- tc[!is.na(s_per_site)]

  points(tc$freq, tc$s_per_site, pch = 21, cex = 1.3,
         col = "black", bg = NA, lwd = 1.2)

  segments(tc$freq, tc$s_per_site,
           tc$freq + tc$dx, tc$s_per_site + tc$dy,
           col = "grey40")

  text(tc$freq + tc$dx, tc$s_per_site + tc$dy,
       labels = tc$consensus_gene_name,
       pos = ifelse(tc$dx < 0, 2, 4), offset = 0.2,
       cex = 0.6, font = 3, xpd = NA)

  mtext(panel_lab, side = 3, line = 0.3, adj = -0.02, font = 2, cex = 1)
}

# Draw --------------------------------------------------------------------
# shared limits taken from the full set so the two panels are directly comparable

xlim_shared <- range(dt_a$freq)
ylim_shared <- range(dt_a$s_per_site)

png(paste0(outdir_fig, "/outliers_seg_sites_per_site_mges.png"),
    width = 23.4, height = 12, units = "cm", res = 600,
    pointsize = 12, type = "cairo")

# the right-hand oma reserves the strip the legend sits in
par(mfrow = c(1, 2), family = "sans", bty = "o", tcl = -0.3,
    mgp = c(2.2, 0.6, 0), oma = c(1, 1, 0, 6.5))

# outer axis titles centre on the two panels, not on the whole device,
# so discount the oma strip given over to the legend
din <- par("din"); omi <- par("omi")
x_mid <- (omi[2] + (din[1] - omi[2] - omi[4]) / 2) / din[1]
y_mid <- (omi[1] + (din[2] - omi[1] - omi[3]) / 2) / din[2]

par(mar = c(2, 3, 2, 0.5))
draw_panel(dt_a, xlim_shared, ylim_shared, panel_lab = "a",
           show_y = TRUE)
mtext("Segregating sites per site", side = 2, outer = TRUE, adj = y_mid)

par(mar = c(2, 1.5, 2, 2))
draw_panel(dt_b, xlim_shared, ylim_shared, panel_lab = "b",
           show_y = FALSE)
mtext(expression("Pangenome frequency ("*italic(k)*")"),
      side = 1, outer = TRUE, adj = x_mid)

# Legend, outer right -----------------------------------------------------
# overlay a plot region spanning the whole device with the outer margins
# zeroed, so the legend can be dropped into the strip reserved by oma above

op <- par(no.readonly = TRUE)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

legend("right", inset = 0.005,
       title = expression(bold("BH FDR < 0.05")), title.adj = 0,
       legend = c("High outlier", "Low outlier", "Non-significant"),
       col = c(rgb(0.2, 0.6, 1, alpha = 1), "#F26C2A", "grey60"),
       pch = c(16, 16, 16),
       bty = "n", cex = 0.7, y.intersp = 1.2, xpd = NA)

par(op)
par(mfrow = c(1, 1))
dev.off()

# Outlier candidates with MGEs removed ------------------------------------

out_b <- merge(dt_b[call != "ns",
                    .(gene_family, freq, gene_length, syn_sites,
                      expected_S, seg_sites, p_two, p_bh, call)],
               anno[, .(gene_family, consensus_gene_name, consensus_product,
                        COG_funct_cat)],
               by = "gene_family", all.x = TRUE)

out_b[, obs_exp_ratio := fifelse(expected_S > 0, seg_sites / expected_S, NA_real_)]

setorder(out_b, call, -obs_exp_ratio)

fwrite(out_b, paste0(outdir_dat, "/Ls_outlier_candidates_no_mge.csv"))

# how the outlier set changes once MGEs are dropped
cat("\nPanel a: ", dt_a[call == "high", .N], " high, ",
    dt_a[call == "low", .N], " low\n",
    "Panel b: ", dt_b[call == "high", .N], " high, ",
    dt_b[call == "low", .N], " low\n", sep = "")
