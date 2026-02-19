
ci <- fread(paste0(outdir_dat,"consistencyindex.csv"))
gl <- fread(paste0(outdir_dat,"gain_loss_results.csv"))

gain_boot <- fread(paste0(outdir_dat, "/gain_loss_boot_results.csv"))



phylo_info = merge(ci, gl, by = "gene_family")

with(phylo_info,
     plot(gains, consistencyindex,
          bty = "L", pch= 19, 
          ylab = "Consistency index",
          xlab = "Number of gains",
          yaxt = "n",
          col = rgb(0.1,0.1,0.8, alpha = 0.9)))

axis(side = 2, at = seq(0,1, 0.2), 
     labels = seq(0,1, 0.2),
     las = 2)


box_plots <- phylo_info[,.(gene_family, consistencyindex, gains)]

box_plots[, ci_bin:= sprintf("%1.2f",round(consistencyindex / 0.05) * 0.05)]


box_plots[, ci_bin:= factor(ci_bin, levels = sprintf("%1.2f", seq(0,1,0.05)))]

with(box_plots,
     boxplot(gains ~ ci_bin, 
             yaxt = "n",
             xlab = "Consistency index bin", 
             ylab = "Number of gains",
             cex= 0.75, bty = "o",
             col = "white",
             border = "white",
             pch = 16, boxwex = 0.5))
abline(h=1, col ="firebrick3")     
axis(side = 2, at = seq(0,50, 10), 
     labels = seq(0,50, 10),
     las = 2)
with(box_plots,
     boxplot(gains ~ ci_bin, add = TRUE,
             xlab = "", ylab = "",
             yaxt = "n", xaxt = "n",
             col = "dodgerblue",
             border = "dodgerblue4",
             pch = 16, boxwex = 0.5))
text(16, 1,"Number of gains = 1", cex = 0.75,
     col = "firebrick3", pos = 3)

