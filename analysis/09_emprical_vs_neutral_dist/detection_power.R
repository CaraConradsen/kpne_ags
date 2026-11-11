# Poisson:   P(S = s)  =  exp(-lambda) * lambda^s / s!
#   at s = 0:  lambda^0 = 1,  0! = 1
#   so         P(S = 0)  =  exp(-lambda)
#   
# min p  =  2 × P(S = 0)  =  2 × exp(-E[S])
# 
# 2 × exp(-E[S])  <=  bh
# exp(-E[S])  <=  bh / 2
# -E[S]   <=  log(bh / 2)
# E[S]   >=  -log(bh / 2)
# 
# with bh = 0.0146:   E[S] >= 4.92

theta_dt <- round(fread(paste0(outdir_dat, "/theta_hat.csv"))[method=="analytic_pi", theta_hat], 
                  digits = 3)

alpha_fdr = 0.05

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))
# gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_m_fdr.csv"))

ag_set_size = length(fread("./output/data/set_1_names.csv")$gene_family)

gene_pvals[, s_per_site := seg_sites / syn_sites]
gene_pvals[, exp_per_site := expected_S / syn_sites]

## three-way call: high, low, or not significant
gene_pvals[, call := "ns"]
gene_pvals[p_bh < alpha_fdr & seg_sites > expected_S, call := "high"]
gene_pvals[p_bh < alpha_fdr & seg_sites < expected_S, call := "low"]
gene_pvals[, call := factor(call, levels = c("ns", "high", "low"))]

cols <- c(ns = rgb(0.1,0.1,0.1, alpha = 0.5), high = "#3FA9F5", low = "#F26C2A")
pt_col <- adjustcolor(cols[as.character(gene_pvals$call)], alpha.f = 0.6)

# The most extreme deficit observable is S = 0, so the floor on the two-sided
# p is 2 * P(S = 0). Under a Poisson null that is 2 * exp(-E[S]) <= bh

bh_threshold   <- gene_pvals[p_bh <= 0.05, max(p_two)]
wall <- -log(bh_threshold / 2) # calculated from poisson E[S]    

sky <- "#3FB0E8"; fire <- "#B22222"


# start plot --------------------------------------------------------------
png(paste0(outdir_fig, "/deficit_detectability_range.png"),
    width = 30, height = 14.25, units = "cm", res = 300, pointsize = 11,
    type = "cairo")

par(mfrow = c(1, 2), mar = c(4, 4, 3, 0.5), bty = "n")

##  how many AGs sit either side of the wall ------------------------
h <- hist(log10(gene_pvals$expected_S), breaks = 60, plot = FALSE)

plot(h, col = ifelse(h$mids < log10(wall), "grey80", sky), border = "white",
     xaxt = "n", main = "",
     xlab = expression("Expected segregating sites, E["*italic(S)*"]"),
     ylab = "Number of accessory genes")

axis(1, at = log10(c(0.01, 0.1, 1, 10, 100)), labels = c(0.01, 0.1, 1, 10, 100))
abline(v = log10(wall), lty = 2,lwd = 1.5, col = fire)
# 
# text(log10(0.1), par("usr")[4] * 0.95,
#      sprintf("%.0f%% of AGs\nno deficit detectable", 100 * mean(gene_pvals$expected_S < wall)),
#      adj = 0, col = "grey30")

text(log10(8), par("usr")[4] * 0.45,
     bquote(bold("E[S] =" ~ .(round(wall, 1)))),
     adj = 0, col = "grey30")

legend("topleft", bty="n",
       col = "grey30", cex = 1.1,
       "P(S = 0)")

mtext("a", side = 3, line = 0, adj = -0.15, font = 2, cex = 1.1)


# detectability limit -----------------------------------------------------
det_lim <- copy(gene_pvals)

det_lim[, testable := expected_S >= wall]
det_lim[, deficit  := p_bh <= alpha_fdr & seg_sites < expected_S]

# t(k) recovered from the gene table
det_lim[, t_k := 2 * expected_S / (theta_hat * syn_sites)]
tk <- det_lim[, .(t_k = as.numeric(median(t_k))), by = freq][order(freq)]
tk[, L_min := 2 * wall / (theta_hat * t_k)]   # boundary: min length at each k

L_med <- det_lim[, median(syn_sites)]
k_med <- tk[L_min <= L_med][1]$freq

cat(sprintf("BH critical value %.4g, E[S] = %.2f\n", bh_threshold, wall))
cat(sprintf("median gene (%d syn sites) testable from k = %d\n", L_med, k_med))

pal <- c(no = "grey80", yes = "grey10", def = "#F26C2A")



# a. the boundary in (k, L_syn) space
with(det_lim,
     plot(freq, syn_sites, type = "n", log = "y",
          yaxt = "n",
          xlab = expression("Pangenome frequency,"~italic(k)),
          ylab = expression("Synonymous sites,"~italic(L[s])~"\n")))

# mtext(side = 3, adj = 0,cex = 0.95, sprintf("BH critical value %.4g, E[S] = %.2f, \nmedian gene (%d synonymous sites) testable from k  = %d\n", 
#                                  bh_threshold, wall, L_med, k_med))

axis(side = 2,
     at     = c(50, 100, 200, 500, 1000, 2000),
     labels = c(50, 100, 200, 500, 1000, 2000),
     las    = 1)

with(det_lim[testable == FALSE, ],
     points(freq, syn_sites,
            pch = 16, cex = 0.6, col = pal[["no"]]))

with(det_lim[testable == TRUE],
     points(freq,  syn_sites,
            pch = 16, cex = 0.6, 
            col = adjustcolor(pal[["yes"]], 0.6)))

with(det_lim[deficit == TRUE],
     points(freq, syn_sites,
            pch = 21, cex = 1, 
            bg = pal[["def"]]))

lines(tk$freq, tk$L_min, col = "black", lwd = 2)

segments(0,L_med, k_med, L_med, lwd =2,
         lty = 3, col = "grey40")

segments(k_med, 1, k_med, L_med, lwd =2,
         lty = 3, col = "grey40")

legend("topright", bty = "n", cex = 0.95,
       pt.cex = c(0.6,0.6,1, NA), pch = c(16, 16, 21, NA),
       lty = c(NA, NA, NA, 1), lwd = c(NA, NA, NA, 2),
       col = c(pal[["no"]], pal[["yes"]], "black", "black"),
       pt.bg = c(NA, NA, pal[["def"]], NA),
       legend = c("Low outlier undetectable", "Low outlier detectable",
                  "Low outlier", "detection boundary"))
mtext("b", side = 3, line = 0, adj = -0.15, font = 2, cex = 1.1)

dev.off()
