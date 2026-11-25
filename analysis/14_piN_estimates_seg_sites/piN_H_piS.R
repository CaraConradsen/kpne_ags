# Protein variation against the synonymous backdrop

q_cut   <- 0.05
# H counts protein haplotypes and so rises with gene length; piN does not,
# being per-site. use_resid replaces H with its residual from a loess on log
# gene length, which puts the two protein measures on a comparable footing.
# piN_H_piS_resid.R sets this and sources this file.
use_resid <- if (exists("use_resid")) use_resid else FALSE
span      <- 0.75
iso     <- c(1, 0.5, 0.1)     # piN/piS lines drawn in (b)
lab_cex <- 0.6

call_col <- c(high = "#3FA9F5", low = "#F26C2A", ns = "grey70")

# Focal AGs, as in ms_figures/allelic_diversity.R
focal_dt <- data.table(
  gene_family = c("g004051", "g003783", "g001147", "g000322",
                  "g003299", "g001254",
                  "g001058",
                  "g003974", "g000990"),
  short = c("wzm", "wzt", "glf", "wbbM", "abiEi", "abiEii",
            "chrA", "tdh", "levR"))

# Data ---------------------------------------------------------------------

alle <- fread(paste0(outdir_dat, "/ag_allele_stats.csv"),
              select = c("gene_family","k","H","n_alleles","call","is_cloud"))
pin  <- fread(paste0(outdir_dat, "/ag_piN_piS_dt.csv"),
              select = c("gene_family","pi_n","pi_s","gene_length","nonsyn_sites",
                         "syn_sites","seg_n_sites","seg_s_sites","piN_piS"))

d <- pin[alle, on = .(gene_family), nomatch = 0]
d[focal_dt, on = .(gene_family), short := i.short]

# the cloud is left out, as in the H panels: k < 15% gives too few carriers
# for either statistic to be worth plotting
d <- d[is_cloud == "n" & !is.na(H) & !is.na(pi_n) & !is.na(pi_s)]

cat(sprintf("shell & soft-core AGs with H and piN: %d  (excess %d, deficit %d, ns %d)\n",
            nrow(d), d[call=="high",.N], d[call=="low",.N], d[call=="ns",.N]))

# the plotted protein-count measure, raw or length-residualised. The loess is
# fitted over every AG in the panel set, then each row is drawn from it, so the
# residual means the same thing in all three rows.
if (use_resid) {
  d[, Hp := residuals(loess(H ~ log(gene_length), .SD, span = span))]
  H_lab  <- expression("Residual " * italic(H) * " (gene length removed)")
  suffix <- "_resid"
  cat(sprintf("length accounts for %.1f%% of the variance in H\n",
              100 * (1 - var(d$Hp) / var(d$H))))
} else {
  d[, Hp := H]
  H_lab  <- expression("Protein allelic diversity (" * italic(H) * ")")
  suffix <- ""
}

rho <- function(x, y) suppressWarnings(cor.test(x, y, method = "spearman"))

fwrite(d, paste0(outdir_dat, "/ag_H_piN_piS_dt.csv"))

# Figure --------------------------------------------------------------------
# Three rows: the non-significant AGs as the reference, then the two kinds of
# departure separately, because they behave nothing alike. Axes are free
# between rows, since the deficit set spans a tenth of the diversity the other
# two do and would otherwise sit in one corner. Ranges are therefore NOT
# comparable row to row: read each row's shape, not its position.

# point area tracks pangenome frequency on a sqrt scale, as in allelic_diversity.R
cex_k <- function(k, lo = 0.30, hi = 1.45) {
  lo + (hi - lo) * (sqrt(k) - sqrt(min(d$k))) / (sqrt(max(d$k)) - sqrt(min(d$k)))
}

panel_tag <- function(tag, dx = 0.145, dy = 0.10) {
  u <- par("usr")
  text(u[1] - diff(u[1:2]) * dx, u[4] + diff(u[3:4]) * dy,
       tag, font = 2, cex = 1.1, xpd = NA)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# one panel for one subset, on axis ranges fixed by the whole dataset
draw_panel <- function(dd, xv, yv, xlab, ylab, xlim, ylim, row_lab = NULL) {
  
  x <- dd[[xv]]; y <- dd[[yv]]
  plot(x, y, type = "n", las = 1, xlab = "", ylab = "", xlim = xlim, ylim = ylim)
  mtext(xlab, side = 1, line = 2.7, cex = 0.68)
  mtext(ylab, side = 2, line = 3.0, cex = 0.68)
  if (!is.null(row_lab)) mtext(row_lab, side = 2, line = 5.0, cex = 0.78, font = 2)
  
  o <- order(dd$k)                       # largest last, so common genes show
  points(x[o], y[o], pch = 16, cex = cex_k(dd$k[o]),
         col = adjustcolor(call_col[dd$call[o]], 0.35))
  
  hl <- which(!is.na(dd$short))
  if (length(hl)) {
    points(x[hl], y[hl], pch = 21, cex = cex_k(dd$k[hl]) + 0.35,
           bg = adjustcolor(call_col[dd$call[hl]], 0.95), col = "grey15", lwd = 0.7)
    text(x[hl], y[hl], dd$short[hl], pos = thigmophobe(x[hl], y[hl]),
         offset = 0.45, cex = lab_cex, font = 3, xpd = NA)
  }
  
  legend("topright", bty = "n", cex = 0.9, inset = c(0.0, -0.1),xpd = TRUE,
         legend = as.expression(bquote(rho ~ "=" ~ .(sprintf("%.2f", rho(x, y)$estimate)) *
                                         "," ~ italic(n) ~ "=" ~ .(nrow(dd))*
                                         "," ~ italic(P) ~ "=" ~ .(sprintf("%.2f", rho(x, y)$p.value)))))
}

# constant piN/piS is a line through the origin, so no log axis is needed and
# the AGs with zero diversity stay on the plot
add_iso <- function() {
  u <- par("usr")
  for (r in iso) {
    abline(a = 0, b = r, col = "grey45", lty = if (r == 1) 2 else 3, lwd = 1.1)
    xr <- min(u[2] * 0.67, (u[4] / r) * 0.83)
    text(xr, r * xr, bquote(pi[N] / pi[S] == .(r)), adj = c(0.5, -0.45),
         cex = 0.8, col = "grey45")
  }
}

rows <- list(list(sub = d[call == "ns"],   lab = "Not significant",   tags = c("a","b","c")),
             list(sub = d[call == "high"], lab = "High synonymous diversity",  tags = c("d","e","f")),
             list(sub = d[call == "low"],  lab = "Low synonymous diversity", tags = c("g","h","i")))


# plot --------------------------------------------------------------------

png(paste0(outdir_fig, "/piN_H_piS", suffix, ".png"),
    width = 28, height = 20, units = "cm", res = 600,
    pointsize = 12, type = "cairo")

layout(matrix(c(1,  2,  3, 10,
                4,  5,  6, 10,
                7,  8,  9, 10), nrow = 3, byrow = TRUE),
       widths = c(1.4, 1.4, 1.4, 0.6), heights = c(1, 1, 1))

for (i in seq_along(rows)) {
  
  rw  <- rows[[i]]
  top <- i == 1
  sb  <- rw$sub
  
  # each row is scaled to its own data
  xs <- c(0, max(sb$pi_s) * 1.04)
  yn <- c(0, max(sb$pi_n) * 1.06)
  xc <- c(0, max(sb$pi_n) * 1.04)
  yh <- if (use_resid) range(d$Hp) * 1.04 else c(0, 1.02)
  
  ## how many variants
  par(mar = c(4.4, 6.4, if (top) 2.2 else 1.2, 1.0), bty = "o")
  draw_panel(sb, "pi_s", "Hp",
             expression("Synonymous diversity (" * pi[S] * ")"),
             H_lab, xs, yh, row_lab = rw$lab)
  if (use_resid) abline(h = 0, lty = 3, col = "grey60")
  panel_tag(rw$tags[1])

  ## how different they are
  par(mar = c(4.4, 5.2, if (top) 2.2 else 1.2, 1.0))
  draw_panel(sb, "pi_s", "pi_n",
             expression("Synonymous diversity (" * pi[S] * ")"),
             bquote(atop("Non-synonymous diversity ( "~ pi[N] ~ ")",
                         "")),
             # expression("Non-synonymous diversity (" * pi[N] * ")"),
             xs, yn)
  add_iso()
  panel_tag(rw$tags[2])

  ## the two protein measures against each other
  par(mar = c(4.4, 5.2, if (top) 2.2 else 1.2, 1.0))
  draw_panel(sb, "pi_n", "Hp",
             expression("Non-synonymous diversity (" * pi[N] * ")"),
             H_lab, xc, yh)
  if (use_resid) abline(h = 0, lty = 3, col = "grey60")
  panel_tag(rw$tags[3])
}

## legends
par(mar = c(0, 0, 0, 0)); plot.new()

legend("topleft", inset = c(0.02, 0.12), xpd = TRUE, bty = "n", cex = 0.62,
       pch = c(16, 16, 16, 21), pt.cex = c(0.9, 0.9, 0.9, 1.1), y.intersp = 1.3,
       col = c(adjustcolor(call_col[c("high","low","ns")], 0.6), "grey15"),
       pt.bg = c(NA, NA, NA, "white"),
       legend = c(sprintf("Excess (n = %d)",  d[call=="high", .N]),
                  sprintf("Deficit (n = %d)", d[call=="low",  .N]),
                  sprintf("ns (n = %d)",      d[call=="ns",   .N]),
                  "Highlighted"))

k_key <- c(50, 150, 250)
k_key <- k_key[k_key <= max(d$k)]
legend("left", inset = c(0.06, 0.0), xpd = TRUE, bty = "n", cex = 0.62,
       pch = 16, col = "grey70", pt.cex = cex_k(k_key), y.intersp = 1.5,
       title = "Pangenome\nfrequency", title.adj = 0, title.cex = 0.62,
       legend = k_key)

dev.off()

