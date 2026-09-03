
alpha_fdr = 0.05

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

gene_pvals[, s_per_site := seg_sites / gene_length]
gene_pvals[, exp_per_site := expected_S / gene_length]

## three-way call: high, low, or not significant
gene_pvals[, call := "ns"]
gene_pvals[p_bh < alpha_fdr & seg_sites > expected_S, call := "high"]
gene_pvals[p_bh < alpha_fdr & seg_sites < expected_S, call := "low"]
gene_pvals[, call := factor(call, levels = c("ns", "high", "low"))]

cols <- c(ns = rgb(0.1,0.1,0.1, alpha = 0.5), high = "#3FA9F5", low = "#F26C2A")
pt_col <- adjustcolor(cols[as.character(gene_pvals$call)], alpha.f = 0.6)


## expectation curve, one value per k (identical across m, so take any)
exp_curve <- unique(gene_pvals[, .(freq, exp_per_site)])[, .(y = mean(exp_per_site)), by = freq]
setorder(exp_curve, freq)

# ## point size carries gene length, the variable no longer on either axis
# cex_pt <- 0.4 + 1.0 * (rank(gene_pvals$gene_length) / nrow(gene_pvals))

# ---------------------------------------------------------------------------

# pdf("segsites_per_site.pdf", wigene_pvalsh = 3.3, height = 3.0, pointsize = 8)

par(mar = c(5, 5, 0.5, 0.5), bty = "n")

plot(gene_pvals$freq, gene_pvals$s_per_site, type = "n",
     xlab = "pangenome frequency, k",
     ylab = "", xaxt = "n", las = 2)

axis(side = 2, at = max(gene_pvals$s_per_site)/2,
     labels = "segregating sites per site", 
     line = 2.5, tick = FALSE)

axis(side = 1, at = pretty(gene_pvals$freq),
     labels = pretty(gene_pvals$freq))

points(gene_pvals$freq, gene_pvals$s_per_site,
       pch = 16, col = pt_col) #  cex = cex_pt,


legend("topright", bty = "n", pch = 16, pt.cex = 0.8, cex = 0.9,
       col = c("black", cols[["high"]], cols[["low"]], cols[["ns"]]),
       lty = c(1, NA, NA, NA), lwd = c(1.5, NA, NA, NA),
       legend = c("neutral expectation",
                  sprintf("excess (BH < %.2f)", alpha_fdr),
                  sprintf("deficit (BH < %.2f)", alpha_fdr),
                  "not significant"))

dev.off()

# ---------------------------------------------------------------------------
# Diagnostic worth looking at before trusting the panel above: observed
# against expected on the raw count scale. If the cloud sits systematically
# off the 1:1 line, the null is miscalibrated for the accessory genome
# rather than individual genes being under selection.

# pdf("obs_vs_exp_diagnostic.pdf", wigene_pvalsh = 3.3, height = 3.3, pointsize = 8)
par(mar = c(3.2, 3.4, 0.6, 0.6), mgp = c(2.0, 0.6, 0), bty = "n",
    family = "sans", las = 1, tcl = -0.3)

lim <- range(c(gene_pvals$expected_S, gene_pvals$seg_sites), finite = TRUE)
plot(gene_pvals$expected_S, gene_pvals$seg_sites, pch = 16, cex = 0.5, col = pt_col,
     xlim = lim, ylim = lim, xlab = "expected S", ylab = "observed S")
abline(0, 1, lty = 2, col = "grey40")

# dev.off()