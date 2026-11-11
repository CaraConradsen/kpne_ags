
# Does the deficit rate fall with pangenome frequency? ---------------------
#
# The two-band contrast (shell vs soft-core) rests on the 95% boundary, and
# a reviewer can always ask why that line and not another. This script
# replaces it with a continuous model: among AGs where a deficit is
# detectable at all, does the probability of being a deficit outlier
# decline with k?
#
# Power is the confound throughout. E[S] = (theta/2) * t(k) * L_syn, so
# both frequency and gene length buy detectability, and the deficit rate
# would rise with k on power grounds alone. Two things handle it. First,
# the analysis is restricted to AGs where the most extreme possible
# observation (S = 0) could reach the FDR threshold; below that the test is
# unattainable, not merely weak. Second, log(syn_sites) enters the model, so
# the coefficient on log(k) is the frequency effect at fixed gene length.
#
# If the coefficient is negative and holds across a range of detectability
# thresholds, the two-band result is a corollary rather than an artefact of
# where the boundary was drawn.

q_cut     <- 0.05   # FDR level used for the outlier calls
es_grid   <- seq(4, 20, by = 1)   # detectability thresholds for the sweep
n_bin     <- 8                    # k bins drawn in panel a

# Outlier calls -----------------------------------------------------------

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites",
                               "syn_sites", "expected_S", "p_two", "p_bh"))

gene_pvals[, call := fcase(p_bh <= q_cut & seg_sites > expected_S, "high",
                           p_bh <= q_cut & seg_sites < expected_S, "low",
                           default = "ns")]

gene_pvals[, `:=`(deficit = as.integer(call == "low"),
                  excess  = as.integer(call == "high"))]

# Detectability -----------------------------------------------------------
# The most extreme deficit observable is S = 0, so the smallest two-sided p
# a gene can ever return is 2 * P(S = 0). Taking the null as Poisson with
# mean E[S] gives 2 * exp(-E[S]). The real null is a Poisson mixture over
# genealogies and so is overdispersed, which makes P(S = 0) larger than
# this and the bound optimistic; the threshold sweep below is what guards
# against that rather than a correction here.
#
# The cut-off is the realised BH critical value, the largest raw p that was
# rejected, rather than q_cut itself, which would be too generous.

bh_cut <- gene_pvals[p_bh <= q_cut, max(p_two)]

gene_pvals[, p_min   := pmin(1, 2 * exp(-expected_S))]
gene_pvals[, testable := p_min <= bh_cut]

cat(sprintf("\nBH critical value in use: %.4g\n", bh_cut))
cat(sprintf("equivalent to E[S] >= %.1f\n", -log(bh_cut / 2)))
cat(sprintf("deficit-testable AGs: %d of %d (%.1f%%)\n",
            gene_pvals[testable == TRUE, .N], nrow(gene_pvals),
            100 * gene_pvals[, mean(testable)]))

print(gene_pvals[, .(AGs = .N, deficit = sum(deficit), excess = sum(excess)),
                 by = testable])

# Any deficit outlier that falls outside the testable set would mean the
# bound is wrong, so check rather than assume.
stopifnot(gene_pvals[deficit == 1 & testable == FALSE, .N] == 0)

# Logistic model ----------------------------------------------------------
# Fitted on the testable AGs only. log(syn_sites) carries the length
# component of power, so the log(freq) term is the frequency effect with
# gene length held fixed.

det <- gene_pvals[testable == TRUE]

m_def <- glm(deficit ~ log(freq) + log(syn_sites),
             data = det, family = binomial())

# the same model for excess-diversity outliers, as a control: if the
# frequency effect is specific to deficits it should not appear here
m_exc <- glm(excess ~ log(freq) + log(syn_sites),
             data = det, family = binomial())

tidy_glm <- function(m, label) {
  s  <- summary(m)$coefficients
  ci <- suppressMessages(confint(m))
  data.table(model = label,
             term  = rownames(s),
             est   = s[, "Estimate"],
             se    = s[, "Std. Error"],
             OR    = exp(s[, "Estimate"]),
             lo    = exp(ci[, 1]),
             hi    = exp(ci[, 2]),
             p     = s[, "Pr(>|z|)"])
}

fit_dt <- rbind(tidy_glm(m_def, "deficit"), tidy_glm(m_exc, "excess"))

cat("\nlogistic models, deficit-testable AGs (n =", nrow(det), ")\n")
print(fit_dt[, .(model, term, OR = round(OR, 3), lo = round(lo, 3),
                 hi = round(hi, 3), p = signif(p, 3))])

fwrite(fit_dt, paste0(outdir_dat, "/ag_detectability_models.csv"))

# A one-unit change on the log scale is hard to read, so express the
# frequency effect as the odds ratio per doubling of k.
ci_def <- fit_dt[model == "deficit" & term == "log(freq)"]
cat(sprintf("\ndeficit odds per doubling of k: %.2f (%.2f-%.2f)\n",
            ci_def$OR^log(2), ci_def$lo^log(2), ci_def$hi^log(2)))

# Threshold sweep ---------------------------------------------------------
# The whole point of the exercise: if the coefficient only appears at one
# detectability cut-off then the objection has simply moved from the 95%
# boundary to this threshold instead.

sweep <- rbindlist(lapply(es_grid, function(e) {
  d <- gene_pvals[expected_S >= e]
  if (d[, sum(deficit)] < 5L || d[, uniqueN(deficit)] < 2L) return(NULL)
  m  <- glm(deficit ~ log(freq) + log(syn_sites), data = d, family = binomial())
  ci <- suppressMessages(confint(m))["log(freq)", ]
  data.table(min_ES = e, n = nrow(d), n_deficit = d[, sum(deficit)],
             est = coef(m)[["log(freq)"]], lo = ci[1], hi = ci[2],
             p = summary(m)$coefficients["log(freq)", "Pr(>|z|)"])
}))

cat("\nsensitivity to the detectability threshold:\n")
print(sweep[, .(min_ES, n, n_deficit,
                OR_per_doubling = round(2^est, 2),
                lo = round(2^lo, 2), hi = round(2^hi, 2),
                p = signif(p, 2))])

fwrite(sweep, paste0(outdir_dat, "/ag_detectability_sweep.csv"))

# Binned rates ------------------------------------------------------------
# Wilson intervals rather than normal-approximation ones, since several
# bins hold no deficits at all and the normal interval is useless there.

wilson <- function(x, n, conf = 0.95) {
  z <- qnorm(1 - (1 - conf) / 2)
  p <- x / n
  d <- 1 + z^2 / n
  c(centre = p,
    lo = ((p + z^2 / (2 * n)) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    hi = ((p + z^2 / (2 * n)) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}

det[, k_bin := cut(freq, unique(quantile(freq, seq(0, 1, length.out = n_bin + 1))),
                   include.lowest = TRUE)]

binned <- det[, {
  w <- wilson(sum(deficit), .N)
  .(n = .N, n_def = sum(deficit), k_med = as.numeric(median(freq)),
    rate = w[["centre"]], lo = w[["lo"]], hi = w[["hi"]])
}, by = k_bin][order(k_med)]

cat("\ndeficit rate by k, deficit-testable AGs only:\n")
print(binned[, .(k_bin, n, n_def, k_med, rate = round(100 * rate, 2))])

# The model fits a monotone curve whether or not the data are monotone, so
# check the bins directly. A strong negative rank correlation here means the
# decline is in the data rather than imposed by the link function.
rho_bin <- suppressWarnings(cor.test(binned$k_med, binned$rate,
                                     method = "spearman"))
cat(sprintf("binned rate vs k: Spearman rho = %.2f, P = %.3f (%d bins)\n",
            rho_bin$estimate, rho_bin$p.value, nrow(binned)))

fwrite(binned, paste0(outdir_dat, "/ag_detectability_binned.csv"))

# The two-band contrast this replaces, kept for comparison ----------------

det[, ag_type := fifelse(freq / 260 > 0.95, "soft", "shell")]

band_tab <- table(factor(det$ag_type, levels = c("shell", "soft")),
                  factor(det$deficit, levels = c(0, 1)))
band_f <- fisher.test(band_tab)

cat("\ntwo-band contrast among deficit-testable AGs:\n")
print(band_tab)
cat(sprintf("OR = %.2f (%.2f-%.2f), P = %.4g\n\n",
            band_f$estimate, band_f$conf.int[1], band_f$conf.int[2],
            band_f$p.value))

# Figure ------------------------------------------------------------------

pal <- c(deficit = "#B22222", excess = "#3FB0E8", ref = "#00285B")

png(paste0(outdir_fig, "/ag_detectability.png"),
    width = 20, height = 9, units = "cm", res = 300,
    pointsize = 11, type = "cairo")

par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.2, 1), bty = "L")

# a. observed rate against the fitted curve
plot(binned$k_med, 100 * binned$rate, type = "n",
     xlab = "Pangenome frequency, k", ylab = "Deficit outliers (%)",
     xlim = range(det$freq), ylim = c(0, 100 * max(binned$hi)))

arrows(binned$k_med, 100 * binned$lo, binned$k_med, 100 * binned$hi,
       code = 3, angle = 90, length = 0.03, col = "grey55", lwd = 1.1)
points(binned$k_med, 100 * binned$rate, pch = 21, bg = pal[["deficit"]],
       col = "white", cex = 1.3, lwd = 0.8)

# fitted curve at the median gene length, so the line is comparable across k
nd <- data.table(freq = seq(min(det$freq), max(det$freq), length.out = 200),
                 syn_sites = det[, median(syn_sites)])
nd[, fit := predict(m_def, newdata = nd, type = "response")]
lines(nd$freq, 100 * nd$fit, col = pal[["ref"]], lwd = 2)

legend("topright", bty = "n", cex = 0.85,
       legend = c("observed (Wilson 95% CI)", "logistic fit, median L"),
       pch = c(21, NA), lty = c(NA, 1), lwd = c(NA, 2), pt.cex = 1.3,
       pt.bg = c(pal[["deficit"]], NA), col = c("white", pal[["ref"]]))

mtext("a", side = 3, line = 0.8, adj = 0, font = 2, cex = 1.1)

# b. is the coefficient an artefact of where the threshold was drawn?
plot(sweep$min_ES, 2^sweep$est, type = "n",
     xlab = expression("Detectability threshold, minimum E["*italic(S)*"]"),
     ylab = "Deficit odds per doubling of k",
     ylim = range(2^sweep$lo, 2^sweep$hi, 1))

abline(h = 1, lty = 2, col = "grey40")
polygon(c(sweep$min_ES, rev(sweep$min_ES)),
        c(2^sweep$lo, rev(2^sweep$hi)),
        col = adjustcolor(pal[["ref"]], 0.15), border = NA)
lines(sweep$min_ES, 2^sweep$est, col = pal[["ref"]], lwd = 2)
points(sweep$min_ES, 2^sweep$est, pch = 16, col = pal[["ref"]], cex = 0.7)

mtext("b", side = 3, line = 0.8, adj = 0, font = 2, cex = 1.1)

dev.off()

cat("wrote", paste0(outdir_fig, "/ag_detectability.png"), "\n")
