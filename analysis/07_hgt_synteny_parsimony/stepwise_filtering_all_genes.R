
id_ags <- unique(fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                select = c("gene_family","number_genomes",
                           "ag_type", "acrs_msu", "acrs_jun")))


dt <- copy(id_ags)

# id set 2 outliers
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

dt[, set := fcase(gene_family %chin% set_1_ags$gene_family,1,
                   default = 0)]

set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

dt[gene_family %chin% set_2_ags$gene_family, 
   set := 2]

# Step 1 — Total
dt1 <- dt
n1 <- nrow(dt1)

# Step 2 — Remove core
dt2 <- dt1[ag_type != "core"]
n2 <- nrow(dt2)

# Step 3 — Remove paralogs
dt3 <- dt2[!grepl("_", gene_family)]
n3 <- nrow(dt3)

# Step 4 — Remove singletons
dt4 <- dt3[number_genomes != 1]
n4 <- nrow(dt4)

# Step 5 — Remove non syntenic
dt5 <- dt4[acrs_msu!=1][acrs_jun!=1]
n5 <- nrow(dt5)

# Step 7 — Remove genes with mutltiple parsimonious gains
dt6 <- dt5[set==2]
n6 <- nrow(dt6)

step_summary <- data.table(
  step = c(
    "Total gene families",
    "After removing core",
    "After removing paralogs",
    "After removing singletons",
    "After removing non-syntenic",
    "After removing mutliple\n parsimonious gainst"
  ),
  n_genes = c(n1, n1-n2, n2-n3, n3-n4, n4-n5, n5-n6),
  cum_ag = c(NA, n2, n3, n4, n5,n6)
)

step_summary

# save data 
fwrite(step_summary, paste0(outdir_dat, "/step_summary.csv"))
       


# Barplot to show before / after filtering --------------------------------

ag_counts <- copy(dt2)

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
  col = c(rgb(0.002,0.265,0.367),
          rgb(0.514,0.592,0.663),
          rgb(0.775,0.775,0.775)),
  xaxt ="n",
  yaxt ="n",
  border = "white",
  space = 0.04,
  font.axis = 2,
  ylab = "",
  xlab = "Number of genomes\n(binned in intervals of 2)"
)

axis(side = 1, at = seq(0, max(bp) + 0.5, length.out = 7),
       labels = round(seq(0,257, length.out = 7), digits = 0))

axis(side = 2, at = seq(0,9, 1),
     labels = seq(0,9, 1),las =1)

axis(side = 2, at = 4.5,
     tick = FALSE, line = 1,
     labels = expression("log"[10]~"Number of genes"))

legend("topright",
       legend = c("All", "Set 1","Set 2"),
       border = c(rgb(0.775,0.775,0.775),
                rgb(0.514,0.592,0.663),
                rgb(0.002,0.265,0.367)),
       fill = c(rgb(0.775,0.775,0.775),
                rgb(0.514,0.592,0.663),
                rgb(0.002,0.265,0.367)),
       bty = "n")

dev.off()
