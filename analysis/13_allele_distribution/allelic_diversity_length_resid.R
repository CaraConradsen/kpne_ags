# Panel allelic_diversity with H residualised on gene length

span       <- 0.75     # loess span for the length fit
fit_within <- FALSE    # TRUE: residualise separately in cloud and shell
q_cut      <- 0.05

call_col <- c(high = "#3FA9F5", low = "#F26C2A")

# Focal AGs, as in allelic_diversity.R. Labels use corrected gene names.
focal_dt <- data.table(
  gene_family = c("g004051", "g003783", "g001147", "g000322",
                  "g003390", "g003914",
                  "g003299", "g001254",
                  "g002556", "g001058",
                  "g002655", "g003239",
                  "g002664", "g002490",
                  "g003377", "g004320",
                  "g002126_1", "g015591", "g006452", "g003974"),
  short = c("wzm", "wzt", "glf", "wbbM", "vgrG", "tssA", "abiEi", "abiEii",
            "asr", "chrA", "acoA", "acoB", "phyH", "phyH2", "sorF", "sorC",
            "glpT", "rplR", "ttdA", "tdh"))

# Data ---------------------------------------------------------------------
# gene_length lives in gene_pvals_fdr.csv, H and the call in ag_allele_stats.csv

alle <- fread(paste0(outdir_dat, "/ag_allele_stats.csv"))
leng <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
              select = c("gene_family", "gene_length"))

d <- leng[alle, on = .(gene_family)][!is.na(H) & !is.na(gene_length) & gene_length > 0]
# ag_allele_stats.csv already carries short and label columns; drop them so
# the highlight set is exactly focal_dt and not whatever that file holds
d[, c("short", "label") := NULL]
d[focal_dt, on = .(gene_family), short := i.short]

# Residualise H on log gene length ------------------------------------------

if (fit_within) {
  d[, rH := residuals(loess(H ~ log(gene_length), span = span)), by = is_cloud]
} else {
  d[, rH := residuals(loess(H ~ log(gene_length), .SD, span = span))]
}

cat(sprintf("H vs log(length): rho = %.3f over %d AGs\n",
            cor(d$H, log(d$gene_length), method = "spearman"), nrow(d)))
cat(sprintf("length accounts for %.1f%% of the variance in H\n",
            100 * (1 - var(d$rH) / var(d$H))))

fwrite(d[, .(gene_family, short, k, gene_length, s_per_site, H, rH, call, is_cloud)],
       paste0(outdir_dat, "/allelic_diversity_length_resid.csv"))

# Figure --------------------------------------------------------------------

all_ags <- d[call != "ns"]

# point area tracks pangenome frequency on a sqrt scale, as in the original
cex_k <- function(k, lo = 0.30, hi = 1.45) {
  lo + (hi - lo) * (sqrt(k) - sqrt(min(all_ags$k))) /
    (sqrt(max(all_ags$k)) - sqrt(min(all_ags$k)))
}

panel_tag <- function(tag, dx = 0.11, dy = 0.12) {
  u <- par("usr")
  text(u[1] - diff(u[1:2]) * dx, u[4] + diff(u[3:4]) * dy,
       tag, font = 2, cex = 1.15, xpd = NA)
}

# Marker colours for the three named genes. Deliberately not the blue and
# orange of the call palette above, which mean excess and deficit.
mk_col <- c(wzm = "#7F5BA8", wbbM = "#2E8B57", chrA = "#C45CA0")

png(paste0(outdir_fig, "/allelic_diversity_length_resid.png"),
    width = 18, height = 16, units = "cm", res = 600,
    pointsize = 11, type = "cairo")

# (a) the fit the residuals come from, (b) the residuals it produces,
# (c, d) those residuals against synonymous diversity; legend down the side
layout(matrix(c(1, 2, 5,
                3, 4, 5), nrow = 2, byrow = TRUE),
       widths = c(1.7, 1.5, 0.5), heights = c(1, 1.15))

## (a) H against gene length, with the loess that defines the residual ------

par(mar = c(4.2, 5.2, 3.0, 1.2), bty = "o")

plot(d$gene_length / 3, d$H, log = "x", pch = 16, cex = 0.3, las = 1,
     col = adjustcolor("grey30", 0.28),
     xlab = "", ylab = "", ylim = c(0, 1.02))
mtext("Protein length (aa)", side = 1, line = 2.6, cex = 0.72)
mtext(expression("Protein allelic diversity (" * italic(H) * ")"),
      side = 2, line = 3.0, cex = 0.72)

gg  <- exp(seq(log(min(d$gene_length)), log(max(d$gene_length)), length.out = 200))
fit <- loess(H ~ log(gene_length), d, span = span)
lines(gg / 3, predict(fit, newdata = data.frame(gene_length = gg)),
      col = "#F26C2A", lwd = 2.5)

mark <- d[gene_family %chin% c("g004051", "g000322", "g001058")]
points(mark$gene_length / 3, mark$H, pch = 21, cex = 1.1, lwd = 0.8,
       bg = mk_col[mark$short], col = "grey15")
text(mark$gene_length / 3, mark$H, mark$short,
     pos = thigmophobe(log(mark$gene_length), mark$H), offset = 0.45,
     cex = 0.58, font = 3, col = mk_col[mark$short])

mtext("a", side = 3, line = 1.4, adj = 0, font = 2, cex = 1.15,
      at = grconvertX(0.005, "nfc", "user"))
mtext(sprintf("loess fit: rho = %.2f, %.1f%% of the variance in H",
              cor(d$H, log(d$gene_length), method = "spearman"),
              100 * (1 - var(d$rH) / var(d$H))),
      side = 3, line = 0.5, adj = 0, cex = 0.66, col = "#F26C2A",
      at = grconvertX(0.06, "nfc", "user"))

## (b) the distribution of residuals -----------------------------------------

par(mar = c(4.2, 3.6, 3.0, 1.2), bty = "o")

hh <- hist(d$rH, breaks = 45, plot = FALSE)
plot(hh, col = "grey85", border = NA, main = "", las = 1,
     xlab = "", ylab = "", ylim = c(0, max(hh$counts) * 1.26))
abline(v = 0, lty = 3, col = "grey45")
mtext(expression("Residual " * italic(H) * " (gene length removed)"),
      side = 1, line = 2.6, cex = 0.72)
mtext("AGs", side = 2, line = 2.6, cex = 0.72)

# tick heights alternate: wbbM and chrA sit 0.04 apart
setorder(mark, rH)
ht <- c(0.72, 1.02)[1 + (seq_len(nrow(mark)) %% 2)]
for (i in seq_len(nrow(mark))) {
  segments(mark$rH[i], 0, mark$rH[i], max(hh$counts) * ht[i],
           col = mk_col[mark$short[i]], lwd = 2)
  # a label sitting directly above a near neighbour's tick is pushed aside
  near <- any(abs(mark$rH[-i] - mark$rH[i]) < 0.08)
  text(mark$rH[i], max(hh$counts) * ht[i], mark$short[i],
       pos = if (!near) 3 else if (ht[i] < 1) 2 else 4,
       cex = 0.58, font = 3, col = mk_col[mark$short[i]], offset = 0.25)
}

mtext("b", side = 3, line = 1.4, adj = 0, font = 2, cex = 1.15,
      at = grconvertX(0.005, "nfc", "user"))
mtext("the two modes are the AGs pinned at H = 0 and H = 1",
      side = 3, line = 0.5, adj = 0, cex = 0.66, col = "grey45",
      at = grconvertX(0.08, "nfc", "user"))

## (c, d) residual H against synonymous diversity ---------------------------

ylim <- range(all_ags$rH) * 1.06
xlim <- c(0, max(all_ags$s_per_site) * 1.05)

for (cld in c("y", "n")) {

  par(mar = if (cld == "y") c(4.2, 5.2, 2.6, 0) else c(4.2, 0, 2.6, 0), bty = "o")

  dd <- d[call != "ns" & is_cloud == cld]
  hl <- dd[!is.na(short)]

  plot(dd$s_per_site, dd$rH, type = "n", las = 1,
       xlim = xlim, ylim = ylim, xlab = "", ylab = "", yaxt = "n")
  abline(h = 0, lty = 3, col = "grey60")

  setorder(dd, -k)      # largest first, so common genes are not hidden
  points(dd$s_per_site, dd$rH, pch = 16, cex = cex_k(dd$k),
         col = adjustcolor(call_col[dd$call], 0.30))

  if (cld == "y") {
    panel_tag("c", dx = 0.15, dy = 0.07)
    axis(2, at = pretty(ylim), las = 2)
    axis(2, at = mean(ylim), tick = FALSE, line = 2.4,
         labels = expression("Residual " * italic(H) * " (gene length removed)"))
  }
  if (cld == "y")
    mtext("Synonymous segregating sites per site", side = 1, line = 2.6,
          cex = 0.72, at = par("usr")[2])

  if (nrow(hl)) {
    points(hl$s_per_site, hl$rH, pch = 21, cex = cex_k(hl$k) + 0.35,
           bg = adjustcolor(call_col[hl$call], 0.95), col = "grey15", lwd = 0.7)
    # the focal genes crowd badly in the upper middle, so each label is
    # pushed away from its nearest neighbours rather than alternated
    text(hl$s_per_site, hl$rH, hl$short,
         pos = thigmophobe(hl$s_per_site, hl$rH), offset = 0.5,
         cex = 0.6, font = 3, xpd = NA)
  }

  legend("top", inset = c(0, -0.12), xpd = TRUE, bty = "n", col = "grey40",
         legend = if (cld == "y")
           expression(bold("Cloud AGs, ") * bolditalic(k) < bold("15%"))
         else
           expression(bold("Shell & Soft-core,") ~ bolditalic(k) >= bold("15%")))

  if (cld == "n") panel_tag("d", dx = 0.03, dy = 0.07)

  rho <- suppressWarnings(cor.test(dd$s_per_site, dd$rH, method = "spearman"))
  p_lab <- if (rho$p.value < 0.001) "< 0.001" else sprintf("= %.3f", rho$p.value)
  legend("topright", xpd = NA, bty = "n", cex = 0.66,
         legend = as.expression(bquote(rho ~ "=" ~ .(sprintf("%.2f", rho$estimate)) *
                                         "," ~ italic(P) ~ .(p_lab))))
}

# legends, as in the original
par(mar = c(0, 0, 0, 0)); plot.new()

legend("topleft", inset = c(0.05, 0.15), xpd = TRUE, bty = "n", cex = 0.66,
       pch = c(16, 16, 21), pt.cex = c(0.9, 0.9, 1.2), y.intersp = 1.25,
       col = c(adjustcolor(call_col, 0.55), "grey15"),
       pt.bg = c(NA, NA, "white"),
       legend = c(sprintf("Excess diversity (n = %d)",
                          all_ags[call == "high" & is_cloud == "n", .N]),
                  sprintf("Deficit diversity (n = %d)",
                          all_ags[call == "low" & is_cloud == "n", .N]),
                  "Highlighted genes"))

k_key <- c(10, 50, 150, 250)
k_key <- k_key[k_key <= max(all_ags$k)]
legend("left", inset = c(0.1, 0.12), xpd = TRUE, bty = "n", cex = 0.66,
       pch = 16, col = "grey70", pt.cex = cex_k(k_key), y.intersp = 1.45,
       title = "Pangenome frequency", title.adj = 0, title.cex = 0.66,
       legend = k_key)

dev.off()

