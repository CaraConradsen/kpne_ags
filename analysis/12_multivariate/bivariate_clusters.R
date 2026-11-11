
# Do H and synonymous departure fall into groups once k is removed? ------
#
# Low-frequency AGs are excluded, because H is a function of the integer
# allele counts among k carriers and at small k its support is a few discrete
# points. Where to cut is an empirical question rather than the 15% cloud
# boundary. The proportion of AGs pinned at H = 0 or H = 1 falls steeply and
# then levels off:
#
#   k >= 5    24% at a boundary      k >= 20    10%
#   k >= 10   14%                    k >= 39     6%
#   k >= 13   12%
#
# and the cluster structure is unchanged across that range (mean silhouette
# 0.39 to 0.40 at G = 4 for every cut from 13 to 39; it only "improves" at
# k >= 5 because the pile of H = 0 genes forms a tight artificial group).
#
# k >= 13 is therefore used: 12% boundary mass, 1,259 AGs rather than 786,
# it is the same 5% threshold Goldfinder already imposes elsewhere in the
# manuscript, and it admits both CRISPR-Cas loci.
#
# Both axes are residualised on log(k) with a loess, because H rises with
# frequency for the same reason power does, and a cluster that simply
# recovered the frequency distribution would say nothing about selection.
# The departure is already conditioned on k through E[S]; residualising it
# as well removes any curvature the offset leaves behind.

q_cut    <- 0.05
k_min    <- 13          # 5% of 260; see note below
# which AGs go into the clustering:
#   "all"  every AG above k_min
#   "bh"   BH-significant departures only
#   "raw"  uncorrected p <= q_cut only
sig_set  <- if (exists("sig_set")) sig_set else "all"
suffix   <- c(all = "", bh = "_bh", raw = "_raw")[sig_set]
span     <- 0.75        # loess span for the k control
g_try    <- 2:8
# the silhouette scan prefers 4 over all AGs but 3 within the outliers, where
# the non-significant middle that cluster 2 occupies has been removed
g_show   <- if (exists("g_show")) g_show else if (sig_set == "all") 4 else 3
set.seed(42)

sky  <- "#3FA9F5"; fire <- "#F26C2A"; grey <- "grey70"

# Genes discussed in the manuscript, grouped by the system they belong to so
# that a whole operon can be seen to travel together (or not).
focal <- data.table(
  gene_family = c("g003179","g005162","g006892","g005811","g005291","g006512","g002579",
                  "g004051","g000322","g003783",
                  "g003299","g001254",
                  "g003390","g000115","g003914","g004223","g005072",
                  "g003974","g005024","g005695","g004284","g008275",
                  "g001557","g001847","g005305",
                  "g001058","g003679","g006237","g007471"),
  gene        = c("cas3(B)","cas1e(B)","cas2e(B)","cas5e(B)","cse2(A)","cas1e(A)","cas3(A)",
                  "wzm","wbbM","rfbB",
                  "abiEi","abiEii",
                  # vgrG-2 (g007173) and vasK-2 (g004984) are dropped: the
                  # annotation is wrong, they belong to other families
                  "vgrG-1","vasK-1","tssA","tssF","tssG",
                  "tdh","ydjK","ydjI","ydjL","glpR",
                  "mrkC","fimD","fimA",
                  "chrA","ccmA","intA","mchF"),
  system      = c(rep("CRISPR-Cas", 7), rep("O-antigen", 3), rep("Defence / TA", 2), rep("T6SS", 5),
                  rep("ydj operon", 5), rep("Fimbrial", 3), rep("Other", 4)))

# Scenario read off each zone's position. Zones are the clusters, except that
# cluster 1 (the synonymous-deficit cloud) is split by allele count: the sweep
# signature is shared, the number of surviving protein alleles is what
# separates a hard sweep from a soft one. Edit freely: these are
# interpretations of where a zone sits, not outputs of the clustering.
if (!exists("scenario")) {
  scenario <- if (sig_set == "all") {
    c("Directional sweep",
      "No synonymous excess\ndespite protein\nvariation",
      "Selection on presence\n(protein conserved)",
      "Allelic NFDS")
  } else {
    c("Directional sweep",
      "Selection on presence\n(protein conserved)",
      "Allelic NFDS")
  }
}
lab_cex  <- 0.58         # scenario label size
lab_pad  <- 0.05         # how far past the hull each label sits, as a
                         # fraction of the axis range
# Manual adjustment per cluster, in axis fractions, applied after the
# automatic placement. Cluster 2 sits near the middle of the plane, so the
# ray that pushes the others outward has no clear direction for it and throws
# it into the right margin; it is dropped into the empty ground below instead.
lab_nudge <- if (sig_set == "all") {
  rbind(c(0, 0), c(-0.03, -0.24), c(0, 0), c(0, 0))
} else {
  rbind(c(0, 0), c(0, 0), c(0, 0))
}

sys_col <- c(`CRISPR-Cas`   = "#E05C9E",
             `O-antigen`    = "#00285B", `Defence / TA` = "#B22222",
             `T6SS`         = "#1B7F5C", `ydj operon`   = "#F26C2A",
             `Fimbrial`     = "#7F5BA8", `Other`        = "#C9A227")

# Data --------------------------------------------------------------------

d <- fread(paste0(outdir_dat, "/ag_feature_table.csv"))

d[, k := freq]
d[, pp  := if (sig_set == "raw") p_two else p_bh]
d[, grp := fcase(pp <= q_cut & seg_sites > expected_S, "high",
                 pp <= q_cut & seg_sites < expected_S, "low",
                 default = "ns")]

d <- d[k >= k_min]

# half-count correction, because 30 AGs observed no segregating sites
d[, dep := log2((seg_sites + 0.5) / (expected_S + 0.5))]

# The k control is fitted on every AG above the cut, then the subset is taken.
# Fitting it inside a significant-only subset would remove a trend estimated
# from genes already chosen for departing, and the subset's k distribution is
# shaped by detectability rather than by biology.
d[, `:=`(rH   = residuals(loess(H   ~ log(k), span = span)),
         rdep = residuals(loess(dep ~ log(k), span = span)))]

d_all <- copy(d)         # panel a always shows the k control as it was fitted

if (sig_set != "all") {
  cat(sprintf("\nrestricted to %s outliers: %d of %d AGs\n",
              if (sig_set == "raw") "uncorrected p <= 0.05" else "BH-significant",
              d[grp != "ns", .N], nrow(d)))
  d <- d[grp != "ns"]
}

cat(sprintf("\n%d non-cloud AGs (k >= %d)\n", nrow(d), k_min))
cat(sprintf("correlation with log(k)  before: H %.3f, departure %.3f\n",
            cor(d$H, log(d$k)), cor(d$dep, log(d$k))))
cat(sprintf("                          after: H %.3f, departure %.3f\n",
            cor(d$rH, log(d$k)), cor(d$rdep, log(d$k))))
cat(sprintf("residual H vs residual departure: r = %.3f\n",
            cor(d$rH, d$rdep)))

X <- scale(as.matrix(d[, .(rH, rdep)]))

# How many groups? --------------------------------------------------------
# Silhouette rather than BIC alone: a Gaussian mixture will keep adding
# components to approximate a continuous density, so BIC can favour a
# partition that does not separate anything.

sw <- rbindlist(lapply(g_try, function(G) {
  km <- kmeans(X, G, nstart = 50)
  data.table(G = G, silhouette = mean(silhouette(km$cluster, dist(X))[, 3]),
             wss = km$tot.withinss)
}))

mc <- Mclust(X, G = g_try, verbose = FALSE)

cat("\ncluster number:\n")
print(sw[, .(G, silhouette = round(silhouette, 3), wss = round(wss))])
cat(sprintf("mclust prefers %s with G = %d; its silhouette is %.3f\n",
            mc$modelName, mc$G,
            mean(silhouette(mc$classification, dist(X))[, 3])))
cat("(silhouette >0.5 strong, 0.25-0.5 weak to moderate, <0.25 none)\n")

fwrite(sw, paste0(outdir_dat, "/bivar_cluster_selection", suffix, ".csv"))

# Groups ------------------------------------------------------------------

km <- kmeans(X, g_show, nstart = 200)

# order clusters by residual departure so the labelling is stable across runs
ord <- order(tapply(d$rdep, km$cluster, median))
d[, cl := factor(match(km$cluster, ord))]

summ <- d[, .(n = .N,
              rH      = round(median(rH), 2),
              rdep    = round(median(rdep), 2),
              H       = round(median(H), 2),
              alleles = as.numeric(median(n_alleles)),
              k       = as.numeric(median(k)),
              high    = sum(grp == "high"),
              low     = sum(grp == "low"),
              ns      = sum(grp == "ns")), by = cl][order(cl)]

cat(sprintf("\nk-means, G = %d (silhouette %.3f):\n", g_show,
            mean(silhouette(as.integer(d$cl), dist(X))[, 3])))
print(summ)

fwrite(cbind(d[, .(gene_family, k, H, n_alleles, dep, rH, rdep, grp)],
             cluster = d$cl),
       paste0(outdir_dat, "/bivar_cluster_assignments", suffix, ".csv"))

fx <- focal[d, on = .(gene_family), nomatch = 0]

cat(sprintf("\nfocal genes: %d of %d retained (cloud AGs are excluded)\n",
            nrow(fx), nrow(focal)))
drop <- focal[!gene_family %chin% d$gene_family, gene]
if (length(drop)) cat("  excluded as cloud:", paste(drop, collapse = ", "), "\n")

cat("\nfocal genes by cluster:\n")
print(fx[, .(gene, system, cluster = cl, k, H = round(H, 2),
             alleles = n_alleles, rdep = round(rdep, 2), call = grp)][order(cluster, system)])

cat("\nfocal genes per cluster, by system:\n")
print(dcast(fx, system ~ cl, fun.aggregate = length, value.var = "gene"))

# Cluster 1 holds both sweep shapes. A hard sweep leaves one protein allele;
# alleles replacing in turn leave several with one dominant, so allele count
# separates them where H cannot.
cat("\ncluster 1, split by allele count:\n")
print(d[cl == 1, .(n = .N, alleles = as.numeric(median(n_alleles)),
                   H = round(median(H), 2), rH = round(median(rH), 2),
                   rdep = round(median(rdep), 2), k = as.numeric(median(k))),
        by = .(shape = fifelse(n_alleles <= 2, "1-2 alleles (hard)",
                               "3+ alleles (soft)"))])

# The presence/absence NFDS scenario also predicts intermediate prevalence.
cat("\nprevalence per cluster:\n")
print(d[, .(n = .N, median_prevalence = round(median(k) / 260, 2),
            pct_between_0.3_0.7 = round(100 * mean(k / 260 >= 0.3 & k / 260 <= 0.7))),
        by = cl][order(cl)])

# Figure ------------------------------------------------------------------

png(paste0(outdir_fig, "/bivar_clusters", suffix, ".png"),
    width = 32, height = 11, units = "cm", res = 300, pointsize = 11,
    type = "cairo")

# panel c carries two dozen gene labels, so it gets the extra width
layout(matrix(1:3, nrow = 1), widths = c(1, 1, 1.75))
par(mar = c(4.3, 4.3, 2.8, 1), bty = "L")

# a. why k has to be removed
plot(d_all$k, d_all$H, log = "x", pch = 16, cex = 0.35,
     col = adjustcolor("grey30", 0.35),
     xlab = expression("Pangenome frequency, "*italic(k)),
     ylab = "Protein allelic diversity (H)")
kk <- exp(seq(log(min(d_all$k)), log(max(d_all$k)), length.out = 100))
lines(kk, predict(loess(H ~ log(k), data = d_all, span = span),
                  newdata = data.frame(k = kk)), col = fire, lwd = 2.5)
legend("bottomright", bty = "n", cex = 0.85, lty = 1, lwd = 2.5, col = fire,
       legend = sprintf("loess (r = %.2f)", cor(d_all$H, log(d_all$k))))
mtext("a  H rises with frequency", side = 3, line = 0.8, adj = 0, font = 2)

# b. how many groups
plot(sw$G, sw$silhouette, type = "b", pch = 16, col = "#00285B", lwd = 2,
     xlab = "Number of clusters", ylab = "Mean silhouette width",
     ylim = c(0, max(sw$silhouette) * 1.25))
abline(h = 0.25, lty = 2, col = "grey50")
text(max(g_try), 0.25, "no structure", pos = 3, cex = 0.72, col = "grey40")
points(g_show, sw[G == g_show, silhouette], pch = 21, bg = fire, cex = 1.5)
mtext("b  Cluster number", side = 3, line = 0.8, adj = 0, font = 2)

# c. the residual plane. Clusters are drawn as outlines so that the points
# themselves are free to carry the call, which is the thing a reader wants to
# check a cluster against.
pal <- c("#B22222", "#E8A33D", "#7FB069", "#3FA9F5", "#7F5BA8", "#888888")
# drop the orange so that a scenario keeps the same colour whether or not the
# non-significant middle is in the plot
if (sig_set != "all") pal <- pal[-2]
call_col <- c(high = sky, low = fire, ns = "grey72")

# a convex hull over the inner hull_q of each cluster: the full hull is
# stretched into a sliver by one or two stray genes
hull_q <- 0.92
draw_hull <- function(x, y, col) {
  dd <- sqrt(((x - median(x)) / sd(x))^2 + ((y - median(y)) / sd(y))^2)
  keep <- dd <= quantile(dd, hull_q)
  x <- x[keep]; y <- y[keep]
  h <- chull(x, y)
  polygon(x[h], y[h], border = col, col = adjustcolor(col, 0.06),
          lwd = 1.6, xpd = NA)
}

# the limits are padded so the scenario labels have somewhere to sit outside
# the hulls without being clipped
plot(d$rH, d$rdep, type = "n",
     xlim = extendrange(d$rH, f = 0.10),
     ylim = extendrange(d$rdep, f = 0.08),
     xlab = "Residual H  (k removed)",
     ylab = "Residual departure, log2 obs/exp  (k removed)")
abline(h = 0, v = 0, lty = 3, col = "grey65")

for (i in seq_len(g_show)) d[cl == i, draw_hull(rH, rdep, pal[i])]

points(d$rH, d$rdep, pch = 16, cex = 0.45,
       col = adjustcolor(call_col[d$grp], 0.7))

points(fx$rH, fx$rdep, pch = 21, bg = sys_col[fx$system], col = "white",
       cex = 1.25, lwd = 1)

# thigmophobe picks a label side pointing away from the nearest neighbour,
# which is the cheapest way to keep two dozen labels legible
lab_pos <- thigmophobe(fx$rH, fx$rdep)
text(fx$rH, fx$rdep, labels = fx$gene, pos = lab_pos, offset = 0.35,
     cex = 0.6, font = 3, col = "grey10", xpd = NA)

# Scenario labels sit outside their own hull. Each is placed on the ray
# running from the centre of the plot through the cluster centroid, at the
# hull vertex furthest along that ray plus a pad, so the label clears the
# cloud it names instead of sitting on top of the points.
usr <- par("usr")
dx <- usr[2] - usr[1]; dy <- usr[4] - usr[3]
gx <- mean(usr[1:2]);  gy <- mean(usr[3:4])

cen <- rbindlist(lapply(seq_len(g_show), function(i) {
  x <- d[cl == i, rH]; y <- d[cl == i, rdep]
  ux <- (median(x) - gx) / dx; uy <- (median(y) - gy) / dy
  r  <- sqrt(ux^2 + uy^2); ux <- ux / r; uy <- uy / r
  p  <- ((x - gx) / dx) * ux + ((y - gy) / dy) * uy   # distance along the ray
  j  <- which.max(p)
  data.table(cl = i, n = length(x),
             x = x[j] + (lab_pad * ux + lab_nudge[i, 1]) * dx,
             y = y[j] + (lab_pad * uy + lab_nudge[i, 2]) * dy)
}))

for (i in seq_len(nrow(cen))) {
  txt <- paste0(scenario[i], "\nn = ", cen$n[i])
  # a cluster sitting near the middle gets an unstable ray and can be thrown
  # off the panel, so keep the whole box inside the plotting region
  hw <- max(strwidth(strsplit(txt, "\n")[[1]], cex = lab_cex)) / 2 * 1.35
  hh <- strheight(txt, cex = lab_cex) / 2 * 1.2
  boxed.labels(min(max(cen$x[i], usr[1] + hw), usr[2] - hw),
               min(max(cen$y[i], usr[3] + hh), usr[4] - hh),
               txt, cex = lab_cex, font = 2,
               col = pal[i], bg = "white", border = NA,
               xpad = 1.1, ypad = 1.2)
}

call_lab <- c(high = "excess (high outlier)", low = "deficit (low outlier)",
              ns = "not significant")
keep_lab <- names(call_lab) %chin% unique(d$grp)
legend("bottomright", bty = "n", cex = 0.68, pch = 16, pt.cex = 1.1,
       col = call_col[keep_lab], legend = call_lab[keep_lab])
legend("topleft", bty = "n", cex = 0.68, pch = 21, pt.cex = 1.1,
       pt.bg = sys_col, col = "white", legend = names(sys_col))
mtext(sprintf("c  Clusters (silhouette %.2f)",
              mean(silhouette(as.integer(d$cl), dist(X))[, 3])),
      side = 3, line = 0.8, adj = 0, font = 2)

dev.off()
cat("\nwrote", paste0(outdir_fig, "/bivar_clusters", suffix, ".png"), "\n")
