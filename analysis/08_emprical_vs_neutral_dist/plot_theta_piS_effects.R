
# Plot the effect of varying piS against observed AGs

# Get FDRs ----------------------------------------------------------------
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))
setnames(ag_age_S_dt, "number_genomes", "freq")

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_s_sites),
                               freq = mean(freq)), gene_family]

# get piS theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]



# get pvals against different expected distributions--------------------
list_prod_lng <- list.files(paste0(outdir_dat,"/theta_pis_out"), pattern = "^prob_pi_lng_0\\.[0-9]{2}\\.csv$",
                            recursive = TRUE, full.names = TRUE)

gene_pvals <- lapply(list_prod_lng, function(prod_lng){
  
  # get piS
  piS <- as.numeric(gsub("prob_pi_lng_", "", gsub(".csv","", basename(prod_lng))))
  
  theta_prod_lng <- fread(prod_lng)
  
  # calculate p
  # add expected S to gene_table
  expected_by_freq <- theta_prod_lng[, .(expected_S = sum(seg_sites * p)), by = freq]
  
  ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = "freq"]
  
  # two-sided p via the "minimum tail doubled" 
  gene_pvals <- ag_age_S_dt[theta_prod_lng, on = .(freq, seg_s_sites = seg_sites), 
                            nomatch = NULL]
  
  # for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
  gene_pvals[, p_lower := cdf]                    # P(X <= observed)
  gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
  gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]
  
  
  # Benjamini-Hochberg
  gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]
  
  gene_pvals[, theta := piS]
  
  gene_pvals[,.(gene_family, freq, seg_s_sites, expected_S, p_two, p_bh, theta)]
  
})

gene_pvals <- rbindlist(gene_pvals)

# get_cis -----------------------------------------------------------------
# get 95% quantiles for the neutral distribution

ci_summary_files <- list.files(paste0(outdir_dat,"/theta_pis_out"), pattern = "^ci_summary_0\\.[0-9]{2}\\.csv$",
                            recursive = TRUE, full.names = TRUE)

ci_summary <- lapply(ci_summary_files, function(ci_sum){
  
  temp_ci_sum <- fread(ci_sum)
  
  temp_ci_sum

})

ci_summary <- rbindlist(ci_summary)

# Examine outliers --------------------------------------------------------

# Join neutral quantiles and define outliers
gene_pvals[ci_summary, on = c("theta","freq"),
        pi_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

gene_pvals[ci_summary, on = c("theta","freq"),
        pi_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]



# Plot --------------------------------------------------------------------

png(paste0(outdir_fig,"/varying_piS.png"),
    width = 17, height = 22, units = "cm", res = 300,
    pointsize = 10, type = "cairo")

par(mfrow=c(5,2), oma=c(5,5,0.5,0.5))
max_y = 30

sm <- function(x, y, span = 0.2) {
  predict(loess(y ~ x, span = span), newdata = data.frame(x = x))
}

thetas <- sort(unique(ci_summary$theta))
for (i in thetas) {
  
  par(mar = c(0,0,0,0))
  
  plot(NULL, pch = 16,#
       col = "grey40",
       xlim = c(0, 260),  
       ylim = c(0, max_y),
       cex = 0.45,
       yaxt = "n",
       xaxt = "n",
       ylab = "",
       xlab = "")
  
  piS_vals <- ci_summary[theta==i][order(freq)]

  up <- sm(piS_vals$freq, piS_vals$upper_95)
  lo <- sm(piS_vals$freq, piS_vals$lower_95)
  
  polygon(c(piS_vals$freq, rev(piS_vals$freq)), c(up, rev(lo)),
              col = "grey80", border = NA)
    
  with(gene_pvals[theta==i], 
       points(freq, seg_s_sites, pch = 16,
              ylim = c(0, 1),
              cex = 0.55,
              yaxt = "n", bty="L",
              ylab = "",
              xlab = "",
              col = ifelse(pi_outlier ==1 & pi_sel  == "balancing" & p_bh <= 0.05,
                           rgb(0.2, 0.6, 1, alpha = 0.5),
                           ifelse(pi_outlier ==1 & pi_sel  != "balancing" & p_bh <= 0.05,
                                  rgb(0.949, 0.424, 0.165, alpha = 0.5), 
                                  rgb(0.1,0.1,0.1, alpha = 0.5))
              )))
  
  bal_AGs = nrow(gene_pvals[theta==i][pi_sel != "directional" & p_bh <= 0.05])
  dir_AGs = nrow(gene_pvals[theta==i][pi_sel == "directional" & p_bh <= 0.05])
  
  mtext(bquote("  "~pi[S]~"="~ .(i)), side = 3, adj = 0, line = -2, cex = 1)
  
  mtext(bquote(.(bal_AGs)~"(bal)"~":"~.(dir_AGs)~"(dir)"), side = 3, line = -2, cex = 0.85)

  
  if (which(thetas == i) %% 2 == 1){
    axis(side = 2, at = seq(0,max_y, length.out = 7), xpd = TRUE,
         labels = sprintf("%3.0f", seq(0,max_y, length.out = 7)), las = 2)
  }
  
  if (which(thetas == i) > length(thetas)-2) {
    axis(side = 1, at = seq(0,250, length.out = 6), xpd = TRUE,
         labels = seq(0,250, length.out = 6))
  }
  
}

mtext("Number of synonymous segregating sites",
      side = 2, outer = TRUE, line = 3.5)
mtext(expression("Pangenome frequency ("*italic(k)*")"),
      side = 1, outer = TRUE, line = 3.5)

dev.off()







