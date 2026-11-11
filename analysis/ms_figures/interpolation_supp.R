png(paste0(outdir_fig,"/interpolation_alpha_&_theta.png"),
    width = 180, height = 80, units = "mm", res = 600,
    pointsize = 12, type = "cairo", family = "sans")

par(mfrow = c(1, 3),
    mar = c(4, 3.5, 1.5, 1.5),   # lines, smaller because pointsize is 7
    mgp = c(1.9, 0.55, 0),          # title, labels, line
    tcl = -0.25, xpd = TRUE,                    # shorter ticks
    cex.axis = 0.9, cex.lab = 1,
    lwd = 0.75)

source("./analysis/08_msprime_sim_alpha_theta/plot_alpha_interprolation.R")
mtext("a", side = 3, adj = -0.27, line = 0.2, font = 2)

source("./analysis/08_msprime_sim_alpha_theta/plot_theta_interpolation.R")

dev.off()