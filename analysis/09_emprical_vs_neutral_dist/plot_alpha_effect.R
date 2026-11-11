
# Plot the effect of demography against observed AGs

# Get FDRs ----------------------------------------------------------------
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_sites),
                               freq = mean(freq)), gene_family]

# get piS theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]



# get pvals against different expected distributions--------------------
list_prod_lng <- list.files(outdir_dat, pattern = "pi_prod_lng_alpha_[0-9.]+\\.csv",
                            recursive = TRUE, full.names = TRUE)

gene_pvals_alpha <- lapply(list_prod_lng, function(prod_lng){
  
  # get alpha
  alpha <- as.numeric(sub(".*_alpha_([0-9.]+)\\.csv$", "\\1", prod_lng))
  
  theta_prod_lng <- fread(prod_lng)
  
  # calculate p
  # add expected S to gene_table
  expected_by_freq <- theta_prod_lng[, .(expected_S = sum(seg_sites * p)), by = freq]
  
  ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = "freq"]
  
  # two-sided p via the "minimum tail doubled" 
  gene_pvals_alpha <- ag_age_S_dt[theta_prod_lng, on = .(freq, seg_s_sites = seg_sites), 
                            nomatch = NULL]
  
  # for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
  gene_pvals_alpha[, p_lower := cdf]                    # P(X <= observed)
  gene_pvals_alpha[, p_upper := 1 - cdf + p]          # P(X >= observed)
  gene_pvals_alpha[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]
  
  
  # Benjamini-Hochberg
  gene_pvals_alpha[, p_bh := p.adjust(p_two, method = "BH")]
  
  gene_pvals_alpha[, alpha := alpha]
  
  gene_pvals_alpha[,.(gene_family, freq, seg_s_sites, expected_S, p_two, p_bh, alpha)]
  
})

gene_pvals_alpha <- rbindlist(gene_pvals_alpha)

# # merge with constant assumptions
# constant_gene_pvals_alpha = fread(paste0(outdir_dat, "/gene_pvals_alpha_fdr.csv"),
#                             select = c("gene_family", "freq", "seg_s_sites", 
#                                        "expected_S", "p_two", "p_bh"))[,alpha := 0.0]
# 
# gene_pvals_alpha <- rbind(gene_pvals_alpha, constant_gene_pvals_alpha)


# get_cis -----------------------------------------------------------------
# get 95% quantiles for the neutral distribution from msprime
ci_summary <- fread(paste0(outdir_dat, "/ci_summary_expanding_pop.csv"))

# # get 95% quantiles for the neutral distribution from rcoal.R
# ci_rcoal <- fread(paste0(outdir_dat, "/ci_summary.csv"))[,alpha := 0.0][,s_set := NULL]
# 
# ci_summary <- rbind(ci_summary,ci_rcoal)

# Examine outliers --------------------------------------------------------

# Join neutral quantiles and define outliers
gene_pvals_alpha[ci_summary, on = c("alpha","freq"),
        pi_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

gene_pvals_alpha[ci_summary, on = c("alpha","freq"),
        pi_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]



# Plot --------------------------------------------------------------------

# png(paste0(outdir_fig,"/demography_2.0.png"),
#     width = 25.7, height = 9.4, units = "cm", res = 300,
#     pointsize = 12, type = "cairo")

#Panel 1
# CI_area -----------------------------------------------------------------
# ci_summary quantiles are counts at the simulated gene length
m_target <- 435
max_y <- 40 / m_target

#theta
with(gene_pvals_alpha[alpha==0],
     plot(freq, seg_s_sites, pch = 16,#
          col = "white",
          xlim = c(0, 260),  
          ylim = c(0, max_y),
          cex = 0.45,
          # yaxs = "i",
          bty = "n",
          yaxt = "n",
          ylab = "Synonymous segregating sites per site\n",
          xlab = expression("Pangenome frequency ("*italic(k)*")"))
)
axis(side = 2, at = seq(0, max_y, length.out = 7),
     labels = sprintf("%.3f", seq(0, max_y, length.out = 7)), las = 2)
alpha_vals <- sort(unique(ci_summary$alpha))
line_cols <- c(
  "#000000",
  "#228833",   # green
  "#AA3377",   # purple
  "#EE6677",   # red-pink
  "#CCBB44"    # yellow-olive
)
names(line_cols) <- as.character(sort(unique(ci_summary$alpha)))

sm <- function(x, y, span = 0.2) {
  predict(loess(y ~ x, span = span), newdata = data.frame(x = x))
}

for (i in alpha_vals) {
  a <- as.character(i)
  d <- ci_summary[alpha == i][order(freq)]
  
  up <- sm(d$freq, d$upper_95 / m_target)
  lo <- sm(d$freq, d$lower_95 / m_target)
  
  # keep the original alpha = 0 as a filled band
  if (i == 0) {
    lines(d$freq, up, col = line_cols[a], lwd = 1.65)
    lines(d$freq, lo, col = line_cols[a], lwd = 1.65)
  }else{
  lines(d$freq, up, col = line_cols[a], lwd = 1.6, lty = "44")
  lines(d$freq, lo, col = line_cols[a], lwd = 1.6, lty = "44")
  }
}

alpha_vals[alpha_vals==2.4] <- 2.8 # Quick fix

legend("topleft", bty = "n",
       legend = sapply(alpha_vals, function(x) bquote(alpha == .(x))),
       col = line_cols, lty = c(1,2,2,2), lwd = 1.8)

# dev.off()





