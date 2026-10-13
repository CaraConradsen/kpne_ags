# output taken from determine_theta.py

curve_dt <- fread(paste0(outdir_dat, "/theta_pi_curve.csv"))
theta_dt <- fread(paste0(outdir_dat, "/theta_hat.csv"))

th_pi <- theta_dt[method == "quadratic_pi"]       # main estimate
th_tw <- theta_dt[method == "quadratic_thetaw"]   # consistency check

pi_obs     <- th_pi$target
thetaw_obs <- th_tw$target
R_obs      <- pi_obs / thetaw_obs
a_hat      <- th_pi$alpha

# quadratic fits through the origin, weights 1/se^2, theta = 0 excluded
# (same fit as fit_quad0 in determine_theta.py)
fit_pi <- lm(mean_pi ~ 0 + theta + I(theta^2),
             data = curve_dt[theta > 0], weights = 1 / se_pi^2)
fit_tw <- lm(mean_thetaw ~ 0 + theta + I(theta^2),
             data = curve_dt[theta > 0], weights = 1 / se_thetaw^2)

theta_grid <- seq(0, max(curve_dt$theta), length.out = 200)


# panel function ------------------------------------------------------------
plot_theta_panel <- function(x, y, lo, hi, fit, obs, est, ylab, obs_lab, panel) {

  ylim <- c(0, max(hi) * 1.05)

  plot(x, y, type = "n", bty = "L", las = 1,
       xaxs = "i", yaxs = "i",
       xlim = c(0, max(x) * 1.02), ylim = ylim,
       xlab = expression("Input " * theta * " (per site)"),
       ylab = ylab)

  mtext(panel, side = 3, adj = 0, line = 0.2, font = 2)

  # 95% CI on the simulated means
  polygon(c(x, rev(x)), c(lo, rev(hi)), col = grey(0.85), border = NA)

  # bootstrap 95% CI on theta_hat
  rect(est$ci_lo, 0, est$ci_hi, obs, col = grey(0.7, alpha = 0.5), border = NA)

  # quadratic fit and grid means
  lines(theta_grid, predict(fit, newdata = data.frame(theta = theta_grid)),
        lwd = 1.25)
  points(x, y, pch = 16, cex = 0.6)

  # observed value vs theta_hat
  segments(0, obs, est$theta_hat, obs, lty = 2, lwd = 1)
  segments(est$theta_hat, 0, est$theta_hat, obs, lty = 2, lwd = 1)
  points(est$theta_hat, obs, pch = 21, bg = "white", cex = 1.1, lwd = 1.2)

  text(est$theta_hat, ylim[2] * 0.08,
       bquote(hat(theta) == .(sprintf("%.4f", est$theta_hat)) ~
                "[" * .(sprintf("%.4f", est$ci_lo)) * "," ~
                .(sprintf("%.4f", est$ci_hi)) * "]"),
       pos = 4, cex = 0.9)

  text(est$theta_hat, obs, obs_lab, pos = 4, cex = 0.9, offset = 0.5)
}


# figure 1: interpolation of theta_hat --------------------------------------
png(paste0(outdir_fig, "/interpolation_theta.png"),
    width = 30, height = 14, units = "cm", res = 300,
    pointsize = 12, type = "cairo")
par(mfrow = c(1, 2), mar = c(4, 4.5, 1.5, 0.5))

with(curve_dt,
     plot_theta_panel(theta, mean_pi, ci_lo_pi, ci_hi_pi, fit_pi,
                      obs = pi_obs, est = th_pi,
                      ylab = expression("Mean simulated " * pi[S]),
                      obs_lab = bquote("Core" ~ pi[S] == .(sprintf("%.4f", pi_obs))),
                      panel = "a"))

with(curve_dt,
     plot_theta_panel(theta, mean_thetaw,
                      mean_thetaw - 1.96 * se_thetaw,
                      mean_thetaw + 1.96 * se_thetaw, fit_tw,
                      obs = thetaw_obs, est = th_tw,
                      ylab = expression("Mean simulated " * theta[W]),
                      obs_lab = bquote("Core" ~ theta[W] == .(sprintf("%.4f", thetaw_obs))),
                      panel = "b"))

dev.off()


# figure 2: diagnostics -----------------------------------------------------
diag_dt <- curve_dt[theta > 0]
diag_dt[, exp_pi := (theta / 2) * mean_pi_branch]
diag_dt[, ratio_lo := (mean_pi - 1.96 * se_pi) / exp_pi]
diag_dt[, ratio_hi := (mean_pi + 1.96 * se_pi) / exp_pi]

png(paste0(outdir_fig, "/interpolation_theta_diagnostics.png"),
    width = 30, height = 14, units = "cm", res = 300,
    pointsize = 12, type = "cairo")
par(mfrow = c(1, 2), mar = c(4, 4.5, 1.5, 0.5))

# (a) simulated pi against the analytic expectation (theta/2) x pi_branch
with(diag_dt,
     plot(theta, pi_over_expected, type = "n", bty = "L", las = 1,
          xlim = c(0, max(theta) * 1.02), xaxs = "i",
          ylim = range(c(ratio_lo, ratio_hi, 1)) + c(-0.01, 0.01),
          xlab = expression("Input " * theta * " (per site)"),
          ylab = expression("Simulated " * pi[S] / "expected " * pi[S])))
mtext("a", side = 3, adj = 0, line = 0.2, font = 2)
abline(h = 1, lty = 2)
with(diag_dt, segments(theta, ratio_lo, theta, ratio_hi, col = grey(0.5)))
with(diag_dt, points(theta, pi_over_expected, pch = 16, cex = 0.7))

# (b) R across theta; should be flat at R_obs if alpha is calibrated
R_range <- range(c(diag_dt$R_ratio_of_means, diag_dt$R_mean_of_ratios, R_obs))

with(diag_dt,
     plot(theta, R_ratio_of_means, type = "n", bty = "L", las = 1,
          xlim = c(0, max(theta) * 1.02), xaxs = "i",
          ylim = R_range + c(-0.02, 0.02),
          xlab = expression("Input " * theta * " (per site)"),
          ylab = expression(pi[S] / theta[W])))
mtext("b", side = 3, adj = 0, line = 0.2, font = 2)
abline(h = R_obs, lty = 2)
with(diag_dt, lines(theta, R_ratio_of_means, lwd = 1.25))
with(diag_dt, points(theta, R_ratio_of_means, pch = 16, cex = 0.7))
with(diag_dt, lines(theta, R_mean_of_ratios, lwd = 1.25, col = grey(0.55)))
with(diag_dt, points(theta, R_mean_of_ratios, pch = 1, cex = 0.7, col = grey(0.55)))

legend("bottomright", bty = "n", cex = 0.85,
       legend = c("Ratio of means", "Mean of ratios",
                  bquote("Core" ~ pi[S] / theta[W] == .(sprintf("%.2f", R_obs)))),
       lty = c(1, 1, 2), pch = c(16, 1, NA),
       col = c("black", grey(0.55), "black"))

text(0, R_range[2] + 0.015, bquote(alpha == .(round(a_hat, 2))), pos = 4, cex = 0.9)

dev.off()
