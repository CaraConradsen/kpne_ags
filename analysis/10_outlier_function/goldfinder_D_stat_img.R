DEstimates <- fread(paste0(outdir_dat, "/DEstimates.csv"))

# Screens are supplied by the caller. They used to be hard-coded to 2 and 3,
# which are the screens holding the cluster heatmaps, so these panels drew
# straight over them and screens 4 and 5 were left empty.
if (!exists("scr_sim")) scr_sim <- 4L
if (!exists("scr_sub")) scr_sub <- 5L

# The direction column was renamed when the outlier set was rebuilt; accept
# either so an older DEstimates.csv still plots.
if (!"w_theta_sel" %in% names(DEstimates) && "pi_sel" %in% names(DEstimates))
  setnames(DEstimates, "pi_sel", "w_theta_sel")

d_cols <- c(balancing = rgb(0.1, 0.48, 0.8, alpha = 0.6),
            directional = rgb(0.95, 0.42, 0.16, alpha = 0.6))

# A gene found by both methods carries "simultaneous & subsequent", so the
# panels are selected by whether the label CONTAINS the method. Testing for
# equality would drop those genes from both panels.
in_sim <- DEstimates[grepl("simultaneous", score)]
in_sub <- DEstimates[grepl("subsequent",   score)]

# genes with no neutrality test have no direction to colour by
in_sim <- in_sim[!is.na(w_theta_sel)]
in_sub <- in_sub[!is.na(w_theta_sel)]

# --- Simultaneous (left) ---
screen(scr_sim)
par(mar = c(4, 4.5, 3, 0.1))

with(in_sim,
     plot(d, -log10(p_brow),
          pch = 16,
          ylim = c(0, 3),
          yaxt = "n",
          ylab = expression("-log"[10] ~ italic("P") * "(Brownian)"),
          xlab = "", main = "",
          col = d_cols[w_theta_sel]))
axis(side = 2, at = seq(0, 3, 0.5),
     labels = sprintf("%1.1f", seq(0, 3, 0.5)), las = 2)
with(in_sim,
     thigmophobe.labels(d, -log10(p_brow),
                        labels = clust_id, cex = 0.5))
mtext("Simultaneous", side = 3, adj = 0.5, line = 0.5)
mtext("d", side = 3, font = 2, cex = 1.1, line = 1.5, adj = -0.05)
abline(v = 0, lty = 2, col = "grey40")
abline(h = -log10(0.05), lty = 2, col = "grey40")

close.screen(2)

# --- Subsequent (right) ---
screen(scr_sub)
par(mar = c(4, 2, 3, 2.6))

with(in_sub,
     plot(d, -log10(p_brow),
          pch = 16,
          ylim = c(0, 3),
          yaxt = "n",
          ylab = "", xlab = "", main = "",
          col = d_cols[w_theta_sel]))
with(in_sub,
     thigmophobe.labels(d, -log10(p_brow),
                        labels = clust_id, cex = 0.5))
mtext("Subsequent", side = 3, adj = 0.5, line = 0.5)
mtext("e", side = 3, font = 2, cex = 1.1, line = 1.5, adj = -0.05)
abline(v = 0, lty = 2, col = "grey40")
abline(h = -log10(0.05), lty = 2, col = "grey40")

mtext("Phylogenetic D statistic", side = 1, 
      adj = -1.5, line = 2.5)

# --- Outer axis labels ---
# Need to return to full-device coordinates to draw outer mtext
close.screen(scr_sub)





