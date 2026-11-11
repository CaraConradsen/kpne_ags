
# Multivariate structure in accessory gene diversity ----------------------
#
# Two questions, kept separate because they need different machinery:
#
#   1. Does protein allelic diversity (H) predict departure from the neutral
#      synonymous expectation, once power and non-independence are accounted
#      for? This is a regression question, not an ordination one.
#
#   2. Is there structure among gene-level diversity statistics that a
#      bivariate correlation misses? This is the ordination question, and it
#      needs more than two variables to be worth asking.
#
# On the response. s_per_site is a ratio of counts and H is bounded on [0, 1]
# with 41% of AGs at exactly 0 and 11% at exactly 1, so neither is a sensible
# Gaussian response. Instead the count itself is modelled, negative binomial
# with offset(log(expected_S)). The offset means every coefficient is a log
# fold-departure from the neutral expectation, and the frequency-dependence of
# power is absorbed exactly rather than approximated, because E[S] already
# contains t(k) and L_syn. H then enters as a predictor, where its boundary
# mass is harmless.
#
# On the random effect. A random intercept treats levels as exchangeable
# draws; k is continuous, ordered, and its effect on E[S] is already known and
# in the offset. Fitting (1 | k) therefore asks the model to absorb residual
# frequency structure non-parametrically across ~250 levels, many holding a
# single gene. It is fitted below anyway and compared, because the comparison
# is more informative than the assertion. The grouping that does carry real
# non-independence is the syntenic junction: genes at one junction are
# physically linked, co-inherited, and share a carrier set.

q_cut      <- 0.05
set.seed(42)

n_pc_show  <- 4          # principal components to tabulate
min_junc   <- 2          # junctions with fewer genes are pooled into "other"
rebuild    <- FALSE      # TRUE to re-derive the feature table from annotation

sky  <- "#3FA9F5"        # excess diversity, as used throughout the ms
fire <- "#F26C2A"        # reduced diversity
grey <- "grey70"

feat_file <- paste0(outdir_dat, "/ag_feature_table.csv")

# 1. Feature table --------------------------------------------------------
# Built once and cached, because the annotation file is large and nothing
# here changes it.

if (rebuild || !file.exists(feat_file)) {

  gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                      select = c("gene_family", "freq", "seg_sites", "syn_sites",
                                 "expected_S", "pi_s", "gene_length",
                                 "p_two", "p_bh"))

  allele_st <- fread(paste0(outdir_dat, "/ag_allele_stats.csv"),
                     select = c("gene_family", "H", "n_alleles", "s_per_site"))

  anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                select = c("gene_family", "ag_type", "msu", "anchor", "ST",
                           "threshold", "egng_cog", "ice", "is_fam",
                           "phg_state", "integron"))

  anno <- anno[gene_family %chin% gene_pvals$gene_family]

  # one row per gene: the junction is unique within set 1 by construction,
  # but ST varies across carriers and the MGE columns are per locus
  anno_g <- anno[, .(
    ag_type   = ag_type[1],
    junction  = paste(msu[1], anchor[1], sep = "_"),
    pir_thr   = as.numeric(threshold[1]),
    cog       = egng_cog[1],
    n_st      = uniqueN(ST),
    mge_any   = as.integer(any(nzchar(trimws(c(ice, is_fam, phg_state, integron)))))
  ), by = gene_family]

  dt <- merge(gene_pvals, allele_st, by = "gene_family")
  dt <- merge(dt, anno_g,            by = "gene_family")

  fwrite(dt, feat_file)

} else {
  dt <- fread(feat_file)
}

# Derived quantities ------------------------------------------------------
# log2 fold-departure uses a half-count correction because 30 AGs observed
# zero segregating sites.

dt[, `:=`(
  k            = freq,
  log2_obs_exp = log2((seg_sites + 0.5) / (expected_S + 0.5)),
  call         = fcase(p_bh <= q_cut & seg_sites > expected_S, "high",
                       p_bh <= q_cut & seg_sites < expected_S, "low",
                       default = "ns"),
  junc_f       = NA_character_)]

dt[, junc_f := fifelse(.N >= min_junc, junction, "other"), by = junction]
dt[, `:=`(junc_f = factor(junc_f),
          k_fac  = factor(freq),
          ag_type = factor(ag_type, levels = c("cloud", "shell", "soft")))]

# centred and scaled predictors, so the intercept is interpretable and the
# two models below are on a common footing
dt[, `:=`(H_z     = as.numeric(scale(H)),
          logk_z  = as.numeric(scale(log(k))),
          logL_z  = as.numeric(scale(log(syn_sites))),
          thr_z   = as.numeric(scale(pir_thr)),
          nst_z   = as.numeric(scale(log(n_st))))]

cat(sprintf("\nfeature table: %d AGs, %d junctions (%d retained as levels)\n",
            nrow(dt), uniqueN(dt$junction), uniqueN(dt$junc_f)))
print(dt[, .N, by = call][order(-N)])

# 2. Is k a usable random effect? -----------------------------------------
# Fitted rather than argued. M0 has no random effect, M1 is the proposal,
# M2 groups by syntenic junction, M3 has both.

form_fix <- seg_sites ~ H_z + logk_z + ag_type + offset(log(expected_S + 0.5))

m0 <- glmmTMB(form_fix, family = nbinom2, data = dt)
m1 <- update(m0, . ~ . + (1 | k_fac))
m2 <- update(m0, . ~ . + (1 | junc_f))
m3 <- update(m0, . ~ . + (1 | junc_f) + (1 | k_fac))

mods <- list(`no random effect` = m0, `(1|k)` = m1,
             `(1|junction)` = m2, `(1|junction)+(1|k)` = m3)

vc_sd <- function(m, grp) {
  v <- glmmTMB::VarCorr(m)$cond
  if (is.null(v[[grp]])) NA_real_ else sqrt(as.numeric(v[[grp]]))
}

comp <- data.table(
  model   = names(mods),
  AIC     = vapply(mods, AIC, numeric(1)),
  logLik  = vapply(mods, function(m) as.numeric(logLik(m)), numeric(1)),
  df      = vapply(mods, function(m) attr(logLik(m), "df"), numeric(1)),
  sd_k    = vapply(mods, vc_sd, numeric(1), grp = "k_fac"),
  sd_junc = vapply(mods, vc_sd, numeric(1), grp = "junc_f"),
  conv    = vapply(mods, function(m) m$fit$convergence == 0, logical(1)))

comp[, dAIC := AIC - min(AIC)]

cat("\nrandom-effect comparison (negative binomial, offset = log E[S]):\n")
print(comp[, .(model, df, AIC = round(AIC, 1), dAIC = round(dAIC, 1),
               sd_k = round(sd_k, 3), sd_junc = round(sd_junc, 3), conv)])

fwrite(comp, paste0(outdir_dat, "/mv_random_effect_comparison.csv"))

# number of k levels holding a single gene, which is why (1|k) is poorly
# identified whatever the AIC says
cat(sprintf("\nk levels: %d, of which %d contain a single AG\n",
            uniqueN(dt$k_fac), dt[, .N, by = k_fac][N == 1, .N]))

# Zero inflation. 44% of AGs have no segregating sites at all, most of them
# because E[S] is small rather than because anything suppressed variation, so
# whether the negative binomial can absorb them on its own is an empirical
# question.

best_nb <- mods[[comp[which.min(AIC), model]]]
best_zi <- update(best_nb, ziformula = ~ logk_z)

zi_cmp <- data.table(model = c("negative binomial", "zero-inflated nb"),
                     AIC   = c(AIC(best_nb), AIC(best_zi)))
zi_cmp[, dAIC := AIC - min(AIC)]

cat("\nzero inflation:\n")
print(zi_cmp[, .(model, AIC = round(AIC, 1), dAIC = round(dAIC, 1))])

best <- if (AIC(best_zi) < AIC(best_nb) - 2) best_zi else best_nb
cat(sprintf("carrying forward: %s\n",
            if (identical(best, best_zi)) "zero-inflated nb" else "negative binomial"))

# 3. The fitted model -----------------------------------------------------

cf <- summary(best)$coefficients$cond
fit_dt <- data.table(term = rownames(cf),
                     est  = cf[, "Estimate"],
                     se   = cf[, "Std. Error"],
                     z    = cf[, "z value"],
                     p    = cf[, "Pr(>|z|)"])
fit_dt[, `:=`(fold = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))]

cat(sprintf("\nselected model: %s\n", comp[which.min(AIC), model]))
print(fit_dt[, .(term, fold = round(fold, 3), lo = round(lo, 3),
                 hi = round(hi, 3), p = signif(p, 3))])

fwrite(fit_dt, paste0(outdir_dat, "/mv_glmm_coefficients.csv"))

# 4. Diagnostics ----------------------------------------------------------
# Simulation-based quantile residuals, computed directly rather than through
# DHARMa so the script carries no extra dependency. Under a correct model
# these are uniform on [0, 1]; the randomised correction handles the
# discreteness of a count response.

sim_quantile_resid <- function(m, n_sim = 250) {
  y   <- model.response(model.frame(m))
  sim <- simulate(m, nsim = n_sim)
  sim <- as.matrix(as.data.frame(sim))
  lo  <- rowMeans(sim <  y)
  eq  <- rowMeans(sim == y)
  pmin(pmax(lo + runif(length(y)) * eq, 1e-6), 1 - 1e-6)
}

qres  <- sim_quantile_resid(best)
fitv  <- fitted(best)
pres  <- residuals(best, type = "pearson")
disp  <- sum(pres^2) / df.residual(best)

ks <- suppressWarnings(ks.test(qres, "punif"))

cat(sprintf("\noverdispersion (Pearson chi-sq / df) = %.3f\n", disp))
# At n = 4,071 the KS test rejects on deviations too small to matter, so the
# statistic D is the quantity to read, not the p-value.
cat(sprintf("quantile residuals vs uniform: KS D = %.3f (P = %.3g, see note)\n",
            ks$statistic, ks$p.value))
cat(sprintf("observed zeros %d, mean simulated zeros %.0f\n",
            sum(model.response(model.frame(best)) == 0),
            mean(colSums(as.matrix(as.data.frame(simulate(best, nsim = 100))) == 0))))

png(paste0(outdir_fig, "/mv_glmm_diagnostics.png"),
    width = 24, height = 16, units = "cm", res = 300, pointsize = 11,
    type = "cairo")

par(mfrow = c(2, 3), mar = c(4.2, 4.2, 2.6, 1), bty = "L")

# a. are the quantile residuals uniform?
qqplot(qunif(ppoints(length(qres))), qres, pch = 16, cex = 0.3,
       col = adjustcolor("grey25", 0.4),
       xlab = "Expected quantile", ylab = "Observed quantile")
abline(0, 1, col = fire, lwd = 2)
mtext("a  Quantile residuals", side = 3, line = 0.8, adj = 0, font = 2)

# b. residual against fitted: should be flat
plot(log10(fitv), qres, pch = 16, cex = 0.3, col = adjustcolor("grey25", 0.4),
     xlab = expression(log[10]~"fitted"), ylab = "Quantile residual")
abline(h = 0.5, lty = 2, col = "grey40")
lines(lowess(log10(fitv), qres, f = 0.3), col = fire, lwd = 2)
mtext("b  Residual vs fitted", side = 3, line = 0.8, adj = 0, font = 2)

# c. residual against k: curvature here means the offset has not absorbed
# the frequency dependence
plot(log10(dt$k), qres, pch = 16, cex = 0.3, col = adjustcolor("grey25", 0.4),
     xlab = expression(log[10]~italic(k)), ylab = "Quantile residual")
abline(h = 0.5, lty = 2, col = "grey40")
lines(lowess(log10(dt$k), qres, f = 0.3), col = fire, lwd = 2)
mtext("c  Residual vs frequency", side = 3, line = 0.8, adj = 0, font = 2)

# d. residual against H, the predictor of interest
plot(dt$H, qres, pch = 16, cex = 0.3, col = adjustcolor("grey25", 0.4),
     xlab = "Protein allelic diversity (H)", ylab = "Quantile residual")
abline(h = 0.5, lty = 2, col = "grey40")
lines(lowess(dt$H, qres, f = 0.3), col = fire, lwd = 2)
mtext("d  Residual vs H", side = 3, line = 0.8, adj = 0, font = 2)

# e. are the random intercepts normal?
re <- ranef(best)$cond
re <- re[[1]][, 1]
qqnorm(re, pch = 16, cex = 0.5, col = "grey25", main = "")
qqline(re, col = fire, lwd = 2)
mtext("e  Random intercepts", side = 3, line = 0.8, adj = 0, font = 2)

# f. observed against fitted counts. A half-count offset keeps the 1,792
# AGs with no segregating sites on the plot; dropping them would hide the
# part of the fit that matters most.
yobs <- model.response(model.frame(best))
plot(fitv + 0.5, yobs + 0.5, pch = 16, cex = 0.3,
     col = adjustcolor("grey25", 0.4), log = "xy",
     xlab = "Fitted + 0.5", ylab = "Observed + 0.5")
abline(0, 1, col = fire, lwd = 2)
mtext(sprintf("f  Observed vs fitted (%d zeros)", sum(yobs == 0)),
      side = 3, line = 0.8, adj = 0, font = 2)

dev.off()
cat("wrote", paste0(outdir_fig, "/mv_glmm_diagnostics.png"), "\n")

# 5. Ordination -----------------------------------------------------------
# A PCA of H against s_per_site alone would be a rotated scatter plot: with
# two variables the first component is just the major axis of their
# correlation and carries no information the correlation does not. The
# ordination is therefore run over the wider set of gene-level statistics,
# and the bivariate relationship is reported as a correlation.

cat(sprintf("\nH vs s_per_site: Spearman rho = %.3f (n = %d)\n",
            cor(dt$H, dt$s_per_site, method = "spearman"), nrow(dt)))

pca_vars <- c("H", "n_alleles", "s_per_site", "pi_s", "log2_obs_exp",
              "k", "syn_sites", "pir_thr", "n_st")

X <- dt[, ..pca_vars]

# counts and frequencies are right-skewed over orders of magnitude
for (v in c("n_alleles", "k", "syn_sites", "n_st")) X[[v]] <- log10(X[[v]])

ok <- complete.cases(X) & is.finite(rowSums(as.matrix(X)))
X  <- X[ok]
d  <- dt[ok]

cat(sprintf("PCA on %d AGs x %d variables (%d dropped for missing or non-finite)\n",
            nrow(X), ncol(X), sum(!ok)))

pca <- prcomp(X, center = TRUE, scale. = TRUE)
ev  <- pca$sdev^2
pv  <- 100 * ev / sum(ev)

cat("\nvariance explained:\n")
print(data.table(PC = paste0("PC", seq_along(pv)),
                 eigenvalue = round(ev, 3),
                 percent = round(pv, 1),
                 cumulative = round(cumsum(pv), 1))[1:n_pc_show])

load_dt <- as.data.table(pca$rotation[, 1:n_pc_show], keep.rownames = "variable")
cat("\nloadings:\n")
print(load_dt[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)])

fwrite(load_dt, paste0(outdir_dat, "/mv_pca_loadings.csv"))
fwrite(cbind(d[, .(gene_family, call, ag_type)],
             as.data.table(pca$x[, 1:n_pc_show])),
       paste0(outdir_dat, "/mv_pca_scores.csv"))

# Broken-stick: components worth interpreting are those whose eigenvalue
# exceeds what a random partition of the variance would give.
p  <- length(ev)
bs <- rev(cumsum(1 / p:1)) * (sum(ev) / p)
cat(sprintf("\ncomponents exceeding broken-stick expectation: %d\n", sum(ev > bs)))

png(paste0(outdir_fig, "/mv_pca.png"),
    width = 26, height = 9, units = "cm", res = 300, pointsize = 11,
    type = "cairo")

par(mfrow = c(1, 3), mar = c(4.2, 4.2, 2.6, 1), bty = "L")

# a. scree against the broken-stick null
bp <- barplot(pv, names.arg = seq_along(pv), col = "grey80", border = "white",
              xlab = "Component", ylab = "Variance explained (%)",
              ylim = c(0, max(pv, 100 * bs / sum(ev)) * 1.1))
lines(bp, 100 * bs / sum(ev), type = "b", pch = 16, col = fire, lwd = 2, cex = 0.7)
legend("topright", bty = "n", cex = 0.85, lty = c(NA, 1), pch = c(15, 16),
       col = c("grey80", fire), legend = c("observed", "broken stick"))
mtext("a  Scree", side = 3, line = 0.8, adj = 0, font = 2)

# b. scores, coloured by call
pal <- c(ns = adjustcolor(grey, 0.35), high = adjustcolor(sky, 0.55),
         low = adjustcolor(fire, 0.9))
ord <- order(factor(d$call, levels = c("ns", "high", "low")))
plot(pca$x[ord, 1], pca$x[ord, 2], pch = 16,
     cex = ifelse(d$call[ord] == "ns", 0.3, 0.5),
     col = pal[d$call[ord]],
     xlab = sprintf("PC1 (%.1f%%)", pv[1]), ylab = sprintf("PC2 (%.1f%%)", pv[2]))
abline(h = 0, v = 0, lty = 3, col = "grey70")
legend("topright", bty = "n", cex = 0.85, pch = 16,
       col = c(grey, sky, fire), legend = c("not significant", "excess", "reduced"))
mtext("b  Scores", side = 3, line = 0.8, adj = 0, font = 2)

# c. loadings
plot(pca$rotation[, 1], pca$rotation[, 2], type = "n",
     xlab = sprintf("PC1 (%.1f%%)", pv[1]), ylab = sprintf("PC2 (%.1f%%)", pv[2]),
     xlim = range(pca$rotation[, 1]) * 1.45,
     ylim = range(pca$rotation[, 2]) * 1.45)
abline(h = 0, v = 0, lty = 3, col = "grey70")
arrows(0, 0, pca$rotation[, 1], pca$rotation[, 2], length = 0.06,
       col = "grey40", lwd = 1.3)
text(pca$rotation[, 1] * 1.18, pca$rotation[, 2] * 1.18,
     labels = rownames(pca$rotation), cex = 0.72, col = "grey15")
mtext("c  Loadings", side = 3, line = 0.8, adj = 0, font = 2)

dev.off()
cat("wrote", paste0(outdir_fig, "/mv_pca.png"), "\n")

# Do the components separate the calls? -----------------------------------
# PCA is unsupervised, so this is a description of the ordination rather
# than a test of it.

cat("\nPC scores by call (median):\n")
print(cbind(d[, .(gene_family, call)],
            as.data.table(pca$x[, 1:n_pc_show]))[
              , lapply(.SD, function(x) round(median(x), 2)),
              by = call, .SDcols = paste0("PC", 1:n_pc_show)][order(call)])

for (i in 1:2) {
  kw <- kruskal.test(pca$x[, i] ~ factor(d$call))
  cat(sprintf("PC%d by call: Kruskal-Wallis chi-sq = %.1f, df = %d, P = %.3g\n",
              i, kw$statistic, kw$parameter, kw$p.value))
}
