# output taken from determine_alpha_and_theta.py

interpo_data <- fread(paste0(outdir_dat, "/alpha_theta_curve.csv"))

R_obs  = 0.600
a_hat  = 3.395
pi_obs = 0.012
pib_hat = 0.7258                      # E[pi_branch | alpha_hat]
th_hat  = 2 * pi_obs / pib_hat        # 0.03307


# major ticks in alpha, evenly spread on the log scale
major <- c(0, 1, 3, 10, 30, 100, 300)
minor <- c(2, 5, 20, 50, 200, 500)

with(interpo_data,
     plot(log10_1p_alpha, R,type = "l", col = "transparent",
          bty = "L", lwd = 1.5, yaxt = "n",
          yaxs = "i", xaxt = "n",
          xlab = expression("Growth rate, " * alpha),
          ylab = expression(pi[S] / theta[W])))

axis(2, at = seq(0.2, 1.0, 0.2), las = 2)

axis(1, at = log10(1 + minor), labels = FALSE, tcl = -0.18)
axis(1, at = log10(1 + major), labels = major)

with(interpo_data,
     polygon(c(log10_1p_alpha, rev(log10_1p_alpha)),
             c(ci_lo, rev(ci_hi)),
             col = grey(0.85), border = NA))

with(interpo_data,
     lines(log10_1p_alpha, R, lwd = 1.25))

# observed R vs alpha
segments(0, R_obs, log10(1 + a_hat), R_obs, lty = 2, lwd = 1)
segments(log10(1 + a_hat), 0.10, log10(1 + a_hat), R_obs, lty = 2, lwd = 1)
points(log10(1 + a_hat), R_obs, pch = 21, bg = "white", cex = 1.1, lwd = 1.2)

text(log10(1 + a_hat), 0.2, bquote(hat(alpha) == .(round(a_hat, 2))),
     pos = 4, cex = 0.9)

text(log10(1 + a_hat), R_obs +0.02, bquote("Core" ~ pi[S] / theta[W] == .(sprintf("%.2f", R_obs))),
     pos = 4,cex = 0.9, offset = 0.2)

