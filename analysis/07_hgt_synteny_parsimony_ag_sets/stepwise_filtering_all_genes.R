
id_ags <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("gene_family","number_genomes", "average_dose", 
                                  "ag_type", "acrs_msu", "acrs_jun")))

dt <- copy(id_ags)

# pi outliers
piS_3IQR_outliers <- fread(paste0(outdir_dat, "/piS_3IQR_outliers.csv"))

# id set 2 outliers
set_1_ags <- fread(paste0(outdir_dat, "/phylo_info.csv"))[syn_jun == 1, gene_family]
set_2_ags <- fread(paste0(outdir_dat, "/phylo_info.csv"))[phy_par_keep == 1, gene_family]


# Step 1 — Total
dt1 <- dt
n1 <- nrow(dt1)

# Step 2 — Remove core
dt2 <- dt1[ag_type != "core"]
n2 <- nrow(dt2)

# Step 3 — Remove paralogs
dt3 <- dt2[average_dose <= 1]
n3 <- nrow(dt3)

# Step 4 — Remove singletons
dt4 <- dt3[number_genomes != 1]
n4 <- nrow(dt4)

# Step 5 — Remove piS outliers
dt5 <- dt4[!gene_family %chin% piS_3IQR_outliers$gene_family]
n5 <- nrow(dt5)

# Step 6 — Remove non syntenic
dt6 <- dt5[gene_family %in% set_1_ags]
n6 <- nrow(dt6)
fwrite(dt6[,.(gene_family)], paste0(outdir_dat, "/set_1_names.csv"))


# Step 7 — Remove genes with multiple parsimonious gains
dt7 <- dt6[gene_family %in% set_2_ags]
n7 <- nrow(dt7)
fwrite(dt7[,.(gene_family)], paste0(outdir_dat, "/set_2_names.csv"))

step_summary <- data.table(
  step = c(
    "Total gene families",
    "After removing core",
    "After removing paralogs",
    "After removing singletons",
    "After removing piS outliers",
    "After removing non-syntenic",
    "After removing mutliple\n parsimonious gainst"
  ),
  n_genes = c(n1, n1-n2, n2-n3, n3-n4, n4-n5, n5-n6, n6-n7),
  cum_ag = c(NA, n2, n3, n4, n5,n6, n7)
)

step_summary

# save data 
fwrite(step_summary, paste0(outdir_dat, "/step_summary.csv"))
       


# Barplot to show before / after filtering --------------------------------

ag_counts <- copy(dt4)

ag_counts$set = "All"
ag_counts[gene_family %chin% piS_3IQR_outliers$gene_family, set:= "pi_out"]
ag_counts[set!= "pi_out" & gene_family %chin% set_1_ags, set:= "set_1"]
ag_counts[gene_family %chin% dt7$gene_family, set:= "set_2"]

ag_counts[, bin := cut(number_genomes,
                       breaks = seq(min(number_genomes), max(number_genomes) + 4, by = 4),
                       right = FALSE)]

ag_counts <- ag_counts[,.(n = .N), by = c("set", "bin")]

ag_counts[, n := log10(n)]

ag_counts <- dcast(ag_counts,  set ~ bin,
                   value.var = "n", fill=0)

ag_counts_plot <- as.matrix(ag_counts[,-1])
rownames(ag_counts_plot) <- ag_counts$set

# reverse row order
ag_counts_plot <- ag_counts_plot[nrow(ag_counts_plot):1, ]


png(paste0(outdir_fig,"/AG_barplot.png"),
    width = 13, height = 15, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

# Get the stacked barplot
par(mar=c(4,3,0.5,0.5))

bp <- barplot(
  ag_counts_plot,
  col = c(rgb(0.05, 0.05, 0.15),      # Very dark blue/navy (almost black)
          rgb(0.05, 0.4, 0.6),      # Stronger cyan-blue
          rgb(0.60, 0.70, 0.80),      # Lighter blue-gray
          rgb(0.80, 0.80, 0.80)),      # Light gray
  xaxt ="n",
  yaxt ="n",
  ylim = c(0,12),
  border = "white",
  space = 0.04,
  font.axis = 2,
  ylab = "",
  xlab = "Number of genomes\n(binned in intervals of 2)"
)

axis(side = 1, at = seq(0, max(bp) + 0.5, length.out = 7),
       labels = round(seq(0,257, length.out = 7), digits = 0))

axis(side = 2, at = seq(0,12, 2),
     labels = seq(0,12, 2),las =1)

axis(side = 2, at = 6,
     tick = FALSE, line = 1,
     labels = expression("log"[10]~"Number of genes"))

legend("topright",
       legend = c("All", expression(pi[S]~"outliers"), "Set 1","Set 2"),
       border = c(rgb(0.80, 0.80, 0.80),
                  rgb(0.60, 0.70, 0.80),
                  rgb(0.05, 0.4, 0.6),
                  rgb(0.05, 0.05, 0.15)),
       fill = c(rgb(0.80, 0.80, 0.80),
                rgb(0.60, 0.70, 0.80),
                rgb(0.05, 0.4, 0.6),
                rgb(0.05, 0.05, 0.15)),
       bty = "n")

dev.off()
