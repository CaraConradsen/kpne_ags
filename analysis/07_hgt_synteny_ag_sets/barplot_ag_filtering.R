# pi outliers
piS_3IQR_outliers <- fread(paste0(outdir_dat, "/piS_3IQR_outliers.csv"))

# id set 2 outliers
set_1_ags <- fread(paste0(outdir_dat, "/phylo_info.csv"))[syn_jun == 1, gene_family]

pre_filter_ags <- fread(paste0(outdir_dat, "/pre_filter_ags.csv"))



# Barplot to show before / after filtering --------------------------------

pre_filter_ags$set = "All"
pre_filter_ags[gene_family %chin% piS_3IQR_outliers$gene_family, set:= "pi_out"]
pre_filter_ags[set!= "pi_out" & gene_family %chin% set_1_ags, set:= "set_1"]

pre_filter_ags[, bin := cut(number_genomes,
                       breaks = seq(min(number_genomes), max(number_genomes) + 3, by = 3),
                       right = FALSE)]

pre_filter_ags <- pre_filter_ags[,.(n = .N), by = c("set", "bin")]

pre_filter_ags[, n := log10(n)]

pre_filter_ags <- dcast(pre_filter_ags,  set ~ bin,
                   value.var = "n", fill=0)

pre_filter_ags_plot <- as.matrix(pre_filter_ags[,-1])
rownames(pre_filter_ags_plot) <- pre_filter_ags$set

# reverse row order
pre_filter_ags_plot <- pre_filter_ags_plot[nrow(pre_filter_ags_plot):1, ]


# Get the stacked barplot
par(mar=c(4,3,1,0.5))

bp <- barplot(
  pre_filter_ags_plot,
  col = c("#243040", "#5C7086", "#C3CDE8"), # "#8FA0B3"
  # for og
  # col = c(rgb(0.05, 0.05, 0.15),      # Very dark blue/navy (almost black)
  #         rgb(0.05, 0.4, 0.6),      # Stronger cyan-blue
  #         rgb(0.60, 0.70, 0.80),      # Lighter blue-gray
  #         rgb(0.80, 0.80, 0.80)),      # Light gray
  xaxt ="n",
  yaxt ="n",
  ylim = c(0,12),
  border = "white",
  space = 0.04,
  font.axis = 2,
  ylab = "",
  xlab = ""
)

axis(side = 1, at = seq(0, max(bp) + 0.5, length.out = 7),
labels = round(seq(0,257, length.out = 7), digits = 0))


axis(side = 1, at = length(bp)/2, labels = expression("Pangenome frequency ("*italic(k)*")"), 
     tick = FALSE, line = 1.5)

axis(side = 2, at = seq(0,12, 2),
     labels = seq(0,12, 2),las =1)

axis(side = 2, at = 6,
     tick = FALSE, line = 1,
     labels = expression("log"[10]~"Number of genes"))


legend("topright",
       legend = c("All AGs",
                  expression("After"~pi[S]~"filter"),
                  "After synteny filter"),
       fill = rev(c("#243040", "#5C7086", "#C3CDE8")), 
       border = rev(c("#243040", "#5C7086", "#C3CDE8")), 
       # border = c(rgb(0.80, 0.80, 0.80),
       #            rgb(0.60, 0.70, 0.80),
       #            rgb(0.05, 0.4, 0.6),
       #            rgb(0.05, 0.05, 0.15)),
       # fill = c(rgb(0.80, 0.80, 0.80),
       #          rgb(0.60, 0.70, 0.80),
       #          rgb(0.05, 0.4, 0.6),
       #          rgb(0.05, 0.05, 0.15)),
       bty = "n")

