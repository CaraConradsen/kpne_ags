# Plot Ci/parsimony against gene frequency (Gambin's plot) --------------------------------------

phylo_info <- fread(paste0(outdir_dat, "/phylo_info.csv"))

png(paste0(outdir_fig,"/gambin_par_gains.png"),
    width = 16, height = 9, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

#syntenic

gene_set <- c(
  "0" = "Unplaced",
  "1" = "Syntenic",
  "2" = "Non-syntenic"
)

par(mfrow = c(1,3),
    mar = c(3, 0,4,0),
    oma = c(2,4,0,0))

for (i in c(2,0,1)) {
  with(phylo_info[pis_out!=1 & syn_jun==i],
       plot(number_genomes, par_gains,
            xlim = c(0,260),
            pch = 19, #bty = "L",
            cex = 0.75, yaxt ="n",
            col = rgb(0.28, 0.52, 0.89, alpha = 0.5),
            xlab = "", 
            ylab = ""))
  
  if(i==2){
    axis(side=2, at = seq(0,40, 10),
         las = 2,
         labels = seq(0,40, 10))
  }
  
  axis(side = 3, at = 125, tick = FALSE,
       labels = bquote(
         .(gene_set[as.character(i)])~" (" * italic(n) ~ "=" ~ .(nrow(phylo_info[pis_out!=1 & syn_jun == i])) * ")"
       ))
  
  
}

mtext("Inferred parsimony gains", side = 2, line =2.75, outer = TRUE)
mtext("Pangenome frequency", side = 1, line = 0.75, outer = TRUE)

dev.off()
