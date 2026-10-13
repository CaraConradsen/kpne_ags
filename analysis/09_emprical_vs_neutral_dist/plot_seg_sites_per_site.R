
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


## expectation curve, one value per k (identical across m, so take any)
exp_curve <- unique(gene_pvals[, .(freq, exp_per_site)])[, .(y = mean(exp_per_site)), by = freq]
setorder(exp_curve, freq)

# ## point size carries gene length, the variable no longer on either axis
# cex_pt <- 0.4 + 1.0 * (rank(gene_pvals$gene_length) / nrow(gene_pvals))

# ---------------------------------------------------------------------------
png(paste0(outdir_fig,"/outliers_seg_sites_per_site.png"),
    width = 17, height = 13, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

par(mar = c(5, 5, 2, 1), bty = "o")

plot(gene_pvals$freq, gene_pvals$s_per_site, type = "n",
     xlab = expression("Pangenome frequency ("*italic(k)*")"),
     ylab = "", xaxt = "n", las = 2, ylim = c(0,0.6))

axis(side = 2, at = max(gene_pvals$s_per_site)/2,
     labels = "segregating sites per site", 
     line = 2.5, tick = FALSE)

axis(side = 1, at = pretty(gene_pvals$freq),
     labels = pretty(gene_pvals$freq))

points(gene_pvals$freq, gene_pvals$s_per_site,
       pch = 16, col = pt_col) #  cex = cex_pt,

# # temp
# acq_res <- fread(paste0(outdir_dat, "/ag_acquisition_ceiling.csv"))
# single_acq_ags <- acq_res[status == "single", gene_family]
# with(gene_pvals[!gene_family %chin% single_acq_ags],
#      points(freq, s_per_site,
#             pch = 16, col = "red"))

# Ed's method
# high_flanks <- fread(paste0(outdir_dat, "/family_2d.csv"))[hgt_class %chin% c("delta_only", "strong"), gene_family]
# with(gene_pvals[gene_family %chin% high_flanks],
#      points(freq, s_per_site,
#             pch = 16, col = "red"))

perc_AGs_out = (nrow(gene_pvals[call!="ns",])/nrow(gene_pvals))*100

mtext(bquote(atop(.((nrow(gene_pvals[call!="ns",])))*" ("*.(sprintf("%2.2f",perc_AGs_out))*"%) AGs",
                  "")), 
      side = 3, line = -1, cex = 0.75)

spearmans <- with(gene_pvals, suppressWarnings(cor.test(freq, s_per_site, method = "spearman")))

rho_val <- sprintf("%.2f", spearmans$estimate)

p_part <- if (spearmans$p.value < 0.001) {
  quote(italic(P) < 0.001)
} else {
  bquote(italic(P) == .(sprintf("%.3f", spearmans$p.value)))
}

mtext(bquote(rho == .(rho_val) * ";" ~ .(p_part)),
      side = 3, adj = 0.5, line = -1.9, cex = 0.8)

legend("topright", #horiz = TRUE,
       legend = c("High outlier (BH FDR < 0.05)","Low outlier (BH FDR < 0.05)",
                  "Non-significant"),
       col = c(rgb(0.2, 0.6, 1, alpha = 1),"#F26C2A", "grey60"),
       pch = c(16, 16,16),
       bty = "n", cex = 0.8)

# sanity checks (new)
trouble_children = data.table(gene_family = c("g000145", "g002611", "g004058",#Outer membrane 
                                              "g000322", "g001147", "g004051", "g003783", "g000817",#5 co-inherited, O-antigen is textbook NFDS
                                              "g002097",#interbacterial competition
                                              "g000382","g000849_1", # Efflux and envelope stress, adaptive variation in efflux capacity 
                                              "g003299"
                                              ), 
                              consensus_gene_name = c("OmpA", "OmpK36","ompR", 
                                                      "wbbM1", "rfbD", "rfbA", "rfbB", "glycosyltransferase",
                                                      "impA",
                                                      "acrA","baeS",
                                                      "AbiEi"
                                                      ),
                              dx = c(-20, -20, -15,
                                     15, -20, -20,  15, -20,
                                     15,
                                     -15, -12,
                                     -15),
                              dy = c( 0.000,  0.005, -0.010,
                                      -0.005,  0.010,  0.010,  0.015, -0.015,
                                      0.010,
                                      -0.005,  0.020,
                                      0.005))

tc <- merge(trouble_children,
            gene_pvals[, .(gene_family, s_per_site, freq)],
            by = "gene_family", all.x = TRUE)

# flag any genes not found in gene_pvals
if (tc[is.na(s_per_site), .N] > 0)
  message("Not plotted: ", paste(tc[is.na(s_per_site), consensus_gene_name], collapse = ", "))

tc <- tc[!is.na(s_per_site)]

# ring the genes of interest
points(tc$freq, tc$s_per_site, pch = 21, cex = 1.6,
       col = "black", bg = NA, lwd = 1.5)

# leader lines run all the way to the label anchor
segments(tc$freq, tc$s_per_site,
         tc$freq + tc$dx, tc$s_per_site + tc$dy,
         col = "grey40")

# pos = 2 puts text left of the anchor, pos = 4 puts it right,
# so the line always meets the near edge of the label
text(tc$freq + tc$dx, tc$s_per_site + tc$dy,
     labels = tc$consensus_gene_name,
     pos = ifelse(tc$dx < 0, 2, 4), offset = 0.2,
     cex = 0.8, font = 3, xpd = NA)
# dev.off()

# ---------------------------------------------------------------------------
# Diagnostic worth looking at before trusting the panel above: observed
# against expected on the raw count scale. If the cloud sits systematically
# off the 1:1 line, the null is miscalibrated for the accessory genome
# rather than individual genes being under selection.

# pdf("obs_vs_exp_diagnostic.pdf", wigene_pvalsh = 3.3, height = 3.3, pointsize = 8)
par(mfrow = c(1,2), mar = c(3.2, 3.4, 0.6, 0.6), mgp = c(2.0, 0.6, 0), bty = "n",
    family = "sans", las = 1, tcl = -0.3)

lim <- range(c(gene_pvals$expected_S, gene_pvals$seg_sites), finite = TRUE)
plot(gene_pvals$expected_S, gene_pvals$seg_sites, pch = 16, cex = 0.5, col = pt_col,
     xlim = lim, ylim = lim, xlab = "expected S", ylab = "observed S")
abline(0, 1, lty = 2, col = "grey40")
with(gene_pvals[gene_family %chin% high_flanks],
     points(expected_S, seg_sites, cex = 0.5,
            pch = 16, col = "red"))

lim <- range(c(gene_pvals$exp_per_site, gene_pvals$s_per_site), finite = TRUE)
plot(gene_pvals$exp_per_site, gene_pvals$s_per_site, pch = 16, cex = 0.5, col = pt_col,
     xlim = lim, ylim = lim, xlab = "expected S per site", ylab = "observed S per site")
abline(0, 1, lty = 2, col = "grey40")
par(mfrow = c(1,1))
# dev.off()


anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                     select = c("gene_family", "consensus_gene_name","consensus_product",
                                "alleles_at_maximum_threshold")))

gene_pvals[call=="low", .(gene_family, freq,gene_length, syn_sites, expected_S, seg_sites)][order(expected_S)]-> temp_low; temp_low 
names_low <- anno[gene_family %chin% gene_pvals[call=="low", gene_family]]
names_low[order(match(gene_family, temp_low$gene_family))]

gene_pvals[call=="high" & freq >50, .(gene_family, freq,gene_length, syn_sites, expected_S, seg_sites)][order(expected_S)]-> temp_high; temp_high  
names_high <- anno[gene_family %chin% gene_pvals[call=="high" & freq >50, gene_family]]
names_high[order(match(gene_family, temp_high$gene_family))]


# output these candidates:
low_cands <- merge(temp_low, names_low, by = "gene_family")
low_cands[, outlier_set := "low"]

high_cands <- merge(temp_high, names_high, by = "gene_family")
high_cands[, outlier_set := "high"]

out <- rbind(low_cands, high_cands)

out[, seg_per_allele := fifelse(alleles_at_maximum_threshold > 0,
                                seg_sites / alleles_at_maximum_threshold,
                                NA_real_)]

out[, obs_exp_ratio := fifelse(expected_S > 0, seg_sites / expected_S, NA_real_)]

fwrite(out,
       paste0(outdir_dat, "/Ls_outlier_candidates.csv"))



# every family whose name or product mentions ompR or envZ
check <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                      select = c("gene_family", "freq", "max_dose", "ag_type",
                                 "consensus_gene_name", "consensus_product")))

check[grepl("ompR|envZ", consensus_gene_name, ignore.case = TRUE) |
        grepl("ompR|envZ|osmolarity|EnvZ", consensus_product, ignore.case = TRUE)][
          order(-freq)]



