# output taken from determine_alpha.py

interpo_data <- fread(paste0(outdir_dat, "/alpha_theta_curve.csv"))
alpha_dt <- fread(paste0(outdir_dat, "/alpha_hat.csv"))

R_obs  = alpha_dt$R_obs
a_hat  = alpha_dt$alpha_hat
a_lo   = alpha_dt$ci_lo
a_hi   = alpha_dt$ci_hi

# major ticks in alpha, evenly spread on the log scale
major <- c(0, 1, 3, 10, 30, 100, 300)
minor <- c(2, 5, 20, 50, 200, 500)

with(interpo_data,
     plot(log10_1p_alpha, R,type = "l", col = "transparent",
          bty = "L", lwd = 1.5, yaxt = "n",
          yaxs = "i", xaxt = "n",
          xlab = expression("\nGrowth rate, " * alpha),
          ylab = ""))

axis(side = 2, at = 0.5, labels = expression(pi[S] / theta[W]), line = 2,
     tick = FALSE)

axis(2, at = seq(0.2, 1.0, 0.2), las = 2)

axis(1, at = log10(1 + minor), labels = FALSE, tcl = -0.18)
axis(1, at = log10(1 + major), labels = major)

# bootstrap 95% CI on alpha_hat, drawn first so the curve sits on top
if (is.finite(a_lo) && is.finite(a_hi)) {
  rect(log10(1 + a_lo), par("usr")[3], log10(1 + a_hi), R_obs,
       col = grey(0.7, alpha = 0.5), border = NA)
}

# 95% CI on the simulated means
with(interpo_data,
     polygon(c(log10_1p_alpha, rev(log10_1p_alpha)),
             c(ci_lo, rev(ci_hi)),
             col = grey(0.85), border = NA))

with(interpo_data,
     lines(log10_1p_alpha, R, lwd = 1.25))

# observed R vs alpha
segments(0, R_obs, log10(1 + a_hat), R_obs, lty = 2, lwd = 1)
segments(log10(1 + a_hat), par("usr")[3], log10(1 + a_hat), R_obs, lty = 2, lwd = 1)
points(log10(1 + a_hat), R_obs, pch = 21, bg = "white", cex = 1.1, lwd = 1.2)

text(log10(1 + a_hat), 0.2,
     bquote(hat(alpha) == .(sprintf("%.2f", a_hat)) ~
              "[" * .(sprintf("%.2f", a_lo)) * "," ~
              .(sprintf("%.2f", a_hi)) * "]"),
     pos = 4, cex = 0.9)

text(log10(1 + a_hat), R_obs +0.02, bquote("Core" ~ pi[S] / theta[W] == .(sprintf("%.2f", R_obs))),
     pos = 4,cex = 0.9, offset = 0.2)

