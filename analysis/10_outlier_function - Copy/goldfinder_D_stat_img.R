DEstimates <- fread(paste0(outdir_dat, "/DEstimates.csv"))

# --- Screen 1: Simultaneous (left) ---
screen(2)
par(mar = c(4, 4.5, 3, 0.1))

with(DEstimates[score == "simultaneous"],
     plot(d, -log10(p_brow),
          pch = 16,
          ylim = c(0, 3),
          yaxt = "n",
          ylab = expression("-log"[10] ~ italic("P") * "(Brownian)"),
          xlab = "", main = "",
          col = ifelse(w_theta_sel == "balancing",
                       rgb(0.1, 0.48, 0.8, alpha = 0.6),
                       rgb(0.95, 0.42, 0.16, alpha = 0.6))))
axis(side = 2, at = seq(0, 3, 0.5),
     labels = sprintf("%1.1f", seq(0, 3, 0.5)), las = 2)
with(DEstimates[score == "simultaneous"],
     thigmophobe.labels(d, -log10(p_brow),
                        labels = clust_id, cex = 0.5))
mtext("Simultaneous", side = 3, adj = 0.5, line = 0.5)
mtext("f", side = 3, font = 2, cex = 1.1, line = 1.5, adj = -0.05)
abline(v = 0, lty = 2, col = "grey40")
abline(h = -log10(0.05), lty = 2, col = "grey40")

close.screen(2)

# --- Screen 2: Subsequent (right) ---
screen(3)
par(mar = c(4, 2, 3, 2.6))

with(DEstimates[score == "subsequent"],
     plot(d, -log10(p_brow),
          pch = 16,
          ylim = c(0, 3),
          yaxt = "n",
          ylab = "", xlab = "", main = "",
          col = ifelse(w_theta_sel == "balancing",
                       rgb(0.1, 0.48, 0.8, alpha = 0.6),
                       rgb(0.95, 0.42, 0.16, alpha = 0.6))))
with(DEstimates[score == "subsequent"],
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
close.screen(3)





