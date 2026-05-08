
piS_dist <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

IQR_threshold = median(piS_dist$pi_s) + (IQR(piS_dist$pi_s) * 3)


ag_gene_dat <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                            select = c("gene_family","number_genomes", "average_dose", 
                                       "ag_type", "acrs_msu", "acrs_jun")))

ag_gene_dat <- ag_gene_dat[ag_type != "core"][average_dose <= 1][number_genomes != 1]

# non syntenic AGs

non_syn <- ag_gene_dat[acrs_msu==1 | acrs_jun ==1, gene_family]

# pi outliers
pis_out <- fread(paste0(outdir_dat, "/piS_3IQR_outliers.csv"))$gene_family

# genes without multiple parsimonious gains
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))$gene_family

multi_gains = ag_gene_dat[!gene_family %chin% set_2_ags, gene_family]

# calcualte euler
fit <- euler(
  list(
    "non-syntenic" = non_syn,
    "piS outliers"  = pis_out,
    "multiple gains"   = multi_gains
  ),
  shape = "ellipse"  # allows ellipses for a better fit when circles can't be perfect
)



# parsimony information ---------------------------------------------------

phylo_info <- fread(, paste0(outdir_dat, "/phylo_info.csv"))

# remove paralogs
phylo_info <- phylo_info[pis_out==0]

box_plots <- phylo_info[syn_jun == 1,.(gene_family, consistencyindex, par_gains)]

box_plots[, ci_bin:=sprintf("%1.2f",round(consistencyindex / 0.05) * 0.05)]

box_plots[, ci_bin:= factor(ci_bin, levels = sprintf("%1.2f", seq(0,1,0.05)))]

# Split gains by ci_bin (drops empty levels)
gain_list <- split(box_plots$par_gains, box_plots$ci_bin, drop = TRUE)

# Remove empty groups
gain_list_clean <- gain_list[sapply(gain_list, length) > 0]

positions <- match(names(gain_list_clean),
                   levels(box_plots$ci_bin))

# plot -----------------------------------------------------------------


png(paste0(outdir_fig,"/filtering_plot.png"),
    width = 16, height = 18.5, units = "cm", res = 300,
    pointsize = 11, type = "cairo")



split.screen(rbind(c(0.2, 0.7, 0.7, 1),
                   c(0, 0.55, 0.53, 0.7),
                   c(0, 0.55, 0.36, 0.53),
                   c(0, 0.55, 0.19, 0.36),
                   c(0, 0.55, 0, 0.19),
                   c(0.5, 1, 0.3, 0.7),
                   c(0.6, 1, 0, 0.3))
)

screen(1)

par(mar = c(4,4,0.5,0.5),
    oma = c(0,1,1,0))

hist(piS_dist$pi_s, 
     breaks = 500,
     xlim = c(0, 0.5),
     ylim = c(0, 3000),
     yaxt = "n",
     ylab = "",
     xlab = "",
     col = "black", 
     main = "")

segments(IQR_threshold,0,
         lwd = 1.1, lty = 3,
         IQR_threshold, 2500,
         col = "blue")

axis(side = 2, at = seq(0,3000, 500),
     labels = seq(0,3000, 500), las = 2)

axis(side = 2, at = 1500,
     line = 2.5, tick = FALSE,
     labels = "Frequency")

text(IQR_threshold,1500,
     paste0("median \u00B1 3*IQR = ",
            sprintf("%1.2f", IQR_threshold)),
     pos = 4, col = "#4444FF", font=2,
     cex = 0.75)

# u <- par("usr")
# xpad <- 0.1 * (u[2] - u[1])
# ypad <- 0.025 * (u[4] - u[3])
# 
# segments(0,0, ((u[1]+u[2])/3) - xpad, (u[3]+u[4])/3 + ypad,
#          col = "grey70", lty = 2)
# 
# segments(0.2,0, u[2]-xpad/8, (u[3]+u[4])/3 + ypad,
#          col = "grey70", lty = 2)

mtext(expression(pi[S]), side =1, line = 2)

mtext("a", side = 3, outer = TRUE,
      adj = 0, font = 2)

# v <- c(
#   grconvertX(u[1:2], "user", "ndc"),
#   grconvertY(u[3:4], "user", "ndc")
# )
# v <- c( (v[1]+v[2])/3 + 0.025, v[2], 0.86, 1 )
# 
# par(fig=v, new=TRUE, mar=c(0,0,0,0.5) )
# 
# #inset
# hist(piS_dist$pi_s, 
#      breaks = 200,
#      xlim = c(0, 0.2),
#      yaxt = "n",
#      ylab = "",
#      xlab = expression(pi[S]),
#      col = "black", 
#      main = "")
# 
# segments(IQR_threshold,0,
#          lwd = 1.1, lty = 3,
#          IQR_threshold, 2500,
#          col = "blue")
# 
# 
# axis(side = 2, at = seq(0,3000, 1000), 
#      labels = seq(0,3000, 1000), las = 2)
# 
# box()

close.screen(1)

# the synteny schematic ---------------------------------------------------

source("./analysis/ms_figures/synteny_schematic.R")


# Parsimony plot ----------------------------------------------------------

screen(6)

# parsimonious gains vs. CI
par(mar = c(4,4,0.25,0.25),
    bty = "n")

plot(NULL, bty = "n",
     xlim = c(1, 21),
     ylim = c(0, 45),
     yaxt = "n", 
     xaxt = "n",
     xlab = "",
     ylab = "")

axis(at  = 20, labels = "Number of gains",
     side = 2, tick = FALSE,
     line = 1.5)

axis(at  = 10.5, labels = "Consistency index", 
     side = 1, tick = FALSE, 
     line = 1.5)

axis(side = 2, at = seq(0,40,10),
     labels = seq(0,40,10), las = 2)

axis(side = 1, at = seq(1,21,4),
     labels = seq(0,1,0.2))

abline(h = 1, lty = 2)

# Add violins
for (i in seq_along(gain_list)) {
  vioplot(gain_list[[i]],
          drawRect = FALSE,
          col = rgb(0.8, 0.54, 0.18, alpha = 0.75),
          border = rgb(0.8, 0.54, 0.18, alpha = 1),
          at = positions[i], 
          bty = "n",   
          boxwex = 1,
          add = TRUE)
  
  boxplot(gain_list[[i]], 
          at = positions[i],
          yaxt = "n",
          outline = FALSE,
          boxwex = 0.5,
          staplewex = 0,
          outwex = 1.5,
          col = "gray20",
          border = "gray20",
          lty = 1,
          pch = 16,
          add = TRUE)
  
  points(positions[i],
         median(gain_list[[i]]),
         pch = 21,
         cex = 0.5,
         bg = "white")
}

# text(16.2, 1, "Number of gains = 1", 
#      col = "grey40", pos = 3, cex = 0.5)

mtext("c", side = 3, 
      font = 2,
      adj = -0.28)

close.screen(6)


# Venn diagram ------------------------------------------------------------


# Prepare a palette of 3 colors with R colorbrewer:
myCol <- brewer.pal(3, "Pastel2")

screen(7)

# Map base-R screen coordinates to grid viewport
vp <- baseViewports()
pushViewport(vp$figure)

pushViewport(viewport(
  x      = 0.5,
  y      = 0.5,
  width  = 0.85,   # adjust for left/right margin
  height = 0.85,   # adjust for top/bottom margin
  just   = c("centre", "centre")
))

venn_cols <- c("grey30","#E63946", "grey70")

venn <- plot(fit,
     fills = list(
       fill  = venn_cols,
       alpha = 0.5
     ),
     
     edges = list(
       lwd = 0.5,
       col = venn_cols,
       alpha = 0.55
     ),
     
     labels = list(
       labels = c(
         "non-syntenic",
         "\u03C0s outliers",
         "parsimonious gains >1"
       ),
       fontfamily = "sans",
       cex        = 0.6
     ),
     
     quantities = list(
       fontface   = "bold",
       fontfamily = "sans",
       cex        = 0.6
     ),
     
     margin = 1
)

grid.draw(venn)

grid.text("d", x = unit(0.02, "npc"), y = unit(0.95, "npc"),
          just = "left", gp =  gpar(fontface = "bold"))


close.screen(7)

close.screen(all.screens=TRUE)

dev.off()

