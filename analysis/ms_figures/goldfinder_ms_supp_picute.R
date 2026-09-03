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

img <- readPNG("./output/figures/clusters_heatmaps/cluster_1.png")

par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
d      <- dim(img)
ar_img <- d[1] / d[2]
pin    <- par("pin")
ar_reg <- pin[2] / pin[1]
if (ar_img > ar_reg) { h <- 1; w <- ar_reg / ar_img } else { w <- 1; h <- ar_img / ar_reg }
rasterImage(img, (1 - w)/2, (1 - h)/2, (1 + w)/2, (1 + h)/2)

mtext("b.", side = 3, adj = 0.1, font = 2)

close.screen(2)

screen(3)

img <- readPNG("./output/figures/clusters_heatmaps/cluster_147.png")

par(mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i")
plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
d      <- dim(img)
ar_img <- d[1] / d[2]
pin    <- par("pin")
ar_reg <- pin[2] / pin[1]
if (ar_img > ar_reg) { h <- 1; w <- ar_reg / ar_img } else { w <- 1; h <- ar_img / ar_reg }
rasterImage(img, (1 - w)/2, (1 - h)/2, (1 + w)/2, (1 + h)/2)

mtext("c.", side = 3, adj = 0.1, font = 2)

close.screen(3)


# D stat ------------------------------------------------------------------

source("./analysis/10_outlier_function/goldfinder_D_stat_img.R")



close.screen(all.screens=TRUE)

dev.off()
