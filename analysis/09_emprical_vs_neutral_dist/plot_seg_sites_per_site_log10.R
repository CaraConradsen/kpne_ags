
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
# png(paste0(outdir_fig,"/outliers_seg_sites_per_site.png"),
#     width = 17, height = 13, units = "cm", res = 300,
#     pointsize = 10, type = "cairo")

par(bty = "o")

plot(gene_pvals$freq, log10(gene_pvals$s_per_site), type = "n", lwd = 1.1,
     xlab = expression("Pangenome frequency ("*italic(k)*")"),
     ylab = "Log10 Segregating sites per site", las=2, xaxt = "n")

# axis(side = 2, at = max(gene_pvals$s_per_site)/2,
#      labels = "Segregating sites per site", 
#      line = 2.5, tick = FALSE, cex = 1)

axis(side = 1, at = seq(0,250,50),
     labels = seq(0,250,50))

points(gene_pvals$freq, log10(gene_pvals$s_per_site),
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
      side = 3, line = -1, cex = 0.9)

spearmans <- with(gene_pvals, suppressWarnings(cor.test(freq, log10(s_per_site), method = "spearman")))

rho_val <- sprintf("%.2f", spearmans$estimate)

p_part <- if (spearmans$p.value < 0.001) {
  quote(italic(P) < 0.001)
} else {
  bquote(italic(P) == .(sprintf("%.3f", spearmans$p.value)))
}

mtext(bquote(rho == .(rho_val) * ";" ~ .(p_part)),
      side = 3, adj = 0.5, line = -1.9, cex = 0.9)

legend("bottomright", #horiz = TRUE,
       title = expression("BH FDR,"~italic(q) < 0.05),
       legend = c("High outlier","Low outlier",
                  "Non-significant"),
       col = c(rgb(0.2, 0.6, 1, alpha = 1),"#F26C2A", "grey60"),
       pch = c(16, 16,16),
       bty = "n", cex = 0.9)

legend("bottom", expression(hat(italic(alpha)) == 2.89*","~hat(italic(theta)) == 0.031),
       bty = "n", cex = 0.9)

# sanity checks (new)
trouble_children = data.table(gene_family = c("g004051", "g003783", "g001147", "g000322", # O1/O2 O-antigen locus (one locus, A308:A309)
                                              "g003390", "g003914",                       # T6SS (one cluster, A164:A165)
                                              "g003299", "g001254",                       # AbiE toxin-antitoxin operon (A22:A23)
                                              "g002556",                                  # acid resistance
                                              "g001058",                                  # chromate efflux
                                              "g002655", "g003239",                       # acetoin catabolism operon (A169:A170)
                                              "g002126_1",                                # glycerol-3-phosphate uptake (fosfomycin entry)
                                              "g006452"),                                 # fumarate hydratase
                              consensus_gene_name = c("wzm", "wzt", "glf", "wbbM",
                                                      "vgrG", "tssA",
                                                      "abiEi", "abiEii",
                                                      "asr",
                                                      "chrA",
                                                      "acoA", "acoB",
                                                      "glpT",
                                                      "fumarate hydratase"),
                              dx = c(12, 12, -6, -8,
                                     -10, -8,
                                     -9, -8,
                                     6, 
                                     -6,
                                     -9,
                                     9, -6,
                                     4), 
                                # 
                                # 
                                # 16,  16, -16, -18,
                                #      -10, -10,
                                #      -13, -9,
                                #       10,
                                #      -10,
                                #      -13,  13,
                                #      -10,
                                #       8),
                              dy = c( 0.060, -0.005,  0.040,  0.000,
                                      0.030, -0.030,
                                      0.000,  -0.015,
                                     -0.040,
                                      0.050,
                                     -0.050, -0.020,
                                      0.040,
                                      0.005))

tc <- merge(trouble_children,
            gene_pvals[, .(gene_family, log10(s_per_site), freq)],
            by = "gene_family", all.x = TRUE)


tc <- tc[!is.na(V2)]

# ring the genes of interest
points(tc$freq, tc$V2, pch = 21, cex = 1.6,
       col = "black", bg = NA, lwd = 1.5)

# leader lines run all the way to the label anchor
segments(tc$freq, tc$V2,
         tc$freq + tc$dx, tc$V2 + tc$dy,
         col = "grey40")

# pos = 2 puts text left of the anchor, pos = 4 puts it right,
# so the line always meets the near edge of the label
text(tc$freq + tc$dx, tc$V2 + tc$dy,
     labels = tc$consensus_gene_name,
     pos = ifelse(tc$dx < 0, 2, 4), offset = 0.2,
     cex = 0.8, font = 3, xpd = NA)
# dev.off()

# 
# gene_pvals[call == "low" & freq < 100] -> ad_track_candidates
# 
# anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
#               select = c("locus_tag","gene_family", "product","gene_length", "gene",
#                          "freq", "ag_type","ST",
#                          "consensus_gene_name","consensus_product",
#                          "alleles_at_maximum_threshold"))
# 
# anno <- anno[gene_family %chin% ad_track_candidates$gene_family]
# anno <- anno[, n_ST_allele := .N, by = c("gene_family", "ST")]
# 
# anno <- unique(anno[,.(gene_family, freq, consensus_gene_name, consensus_product,ST, n_ST_allele, alleles_at_maximum_threshold)])
# 
# anno <- anno[, n_ST := uniqueN(ST), gene_family]
# 
# anno <- dcast(anno, gene_family + freq + consensus_gene_name + consensus_product + alleles_at_maximum_threshold + n_ST ~ ST,  value.var = "n_ST_allele" )
# 
# fwrite(anno, "C:/Users/carac/Desktop/low_ags.csv")
