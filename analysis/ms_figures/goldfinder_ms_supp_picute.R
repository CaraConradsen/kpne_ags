# Supplementary figure: Goldfinder non-independence ------------------------
#   a  co-occurrence network, by scoring method
#   b  cluster 1 presence/absence heatmap
#   c  cluster 147 presence/absence heatmap
#   d  phylogenetic D, simultaneous
#   e  phylogenetic D, subsequent
#
# The D panels used to be drawn into screens 2 and 3, which already hold the
# heatmaps, so they overwrote them and screens 4 and 5 stayed blank. The
# screens are now passed explicitly.

heat_1   <- "./output/figures/clusters_heatmaps/cluster_1.png"
heat_2   <- "./output/figures/clusters_heatmaps/cluster_147.png"

stopifnot(file.exists(heat_1), file.exists(heat_2))

# one raster panel, scaled to fit without distorting the aspect ratio
draw_panel_png <- function(path, tag) {
  img <- readPNG(path)
  par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  d      <- dim(img)
  ar_img <- d[1] / d[2]
  pin    <- par("pin")
  ar_reg <- pin[2] / pin[1]
  if (ar_img > ar_reg) { h <- 1; w <- ar_reg / ar_img } else { w <- 1; h <- ar_img / ar_reg }
  rasterImage(img, (1 - w)/2, (1 - h)/2, (1 + w)/2, (1 + h)/2)
  mtext(tag, side = 3, adj = 0.1, font = 2)
}

# plot -----------------------------------------------------------------


png(paste0(outdir_fig,"/goldfinder_results.png"),
    width = 17, height = 22,  units = "cm", res = 300,
    pointsize = 11, type = "cairo")


split.screen(rbind(c(0,1,0.66,1),
                   c(0,0.5,0.33,0.6),
                   c(0.5,1,0.33,0.6),
                   c(0,0.5,0,0.3),
                   c(0.5,1,0,0.3)))

# the synteny schematic ---------------------------------------------------

screen(1)

source("./analysis/10_outlier_function/goldfinder_network_img.R")

close.screen(1)


screen(2)
draw_panel_png(heat_1, "b")
close.screen(2)

screen(3)
draw_panel_png(heat_2, "c")
close.screen(3)


# D stat ------------------------------------------------------------------
# into the empty lower screens, not over the heatmaps

scr_sim <- 4L
scr_sub <- 5L

source("./analysis/10_outlier_function/goldfinder_D_stat_img.R")



close.screen(all.screens=TRUE)

dev.off()
