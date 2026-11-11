png(paste0(outdir_fig,"/dist_seg_sites_vs_neutral.png"),
    width = 266, height = 126, units = "mm", res = 600,
    pointsize = 12, type = "cairo", family = "sans")

par(mar = c(4, 5, 2, 1),   # lines, smaller because pointsize is 7
    tcl = -0.25, xpd = TRUE,                    # shorter ticks
    cex.axis = 0.9, cex.lab = 1)


# dev.off()


split.screen(rbind(c(0, 0.35, 0, 1),
                   c(0.35, 1, 0, 1))
)


screen(1)

source("./analysis/09_emprical_vs_neutral_dist/plot_alpha_effect.R")
mtext("a", side = 3, adj = -0.4, font = 2)

close.screen(1)

screen(2)

source("./analysis/09_emprical_vs_neutral_dist/plot_seg_sites_per_site.R")
mtext("b", side = 3, adj = -0.15, font = 2)

close.screen(2)


close.screen(all.screens=TRUE)

dev.off()
