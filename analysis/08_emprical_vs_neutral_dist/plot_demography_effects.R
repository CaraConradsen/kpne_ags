
# Plot the effect of demography against observed AGs

# Get FDRs ----------------------------------------------------------------
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))
setnames(ag_age_S_dt, "number_genomes", "freq")

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_s_sites),
                               freq = mean(freq)), gene_family]

# get piS theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]



# get pvals against different expected distributions--------------------
list_prod_lng <- list.files(outdir_dat, pattern = "pi_prod_lng_alpha_[0-9.]+\\.csv",
                            recursive = TRUE, full.names = TRUE)

gene_pvals <- lapply(list_prod_lng, function(prod_lng){
  
  # get alpha
  alpha <- as.numeric(sub(".*_alpha_([0-9.]+)\\.csv$", "\\1", prod_lng))
  
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
  
  gene_pvals[, alpha := alpha]
  
  gene_pvals[,.(gene_family, freq, seg_s_sites, expected_S, p_two, p_bh, alpha)]
  
})

gene_pvals <- rbindlist(gene_pvals)

# # merge with constant assumptions
# constant_gene_pvals = fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
#                             select = c("gene_family", "freq", "seg_s_sites", 
#                                        "expected_S", "p_two", "p_bh"))[,alpha := 0.0]
# 
# gene_pvals <- rbind(gene_pvals, constant_gene_pvals)


# get_cis -----------------------------------------------------------------
# get 95% quantiles for the neutral distribution from msprime
ci_summary <- fread(paste0(outdir_dat, "/ci_summary_expanding_pop.csv"))

# # get 95% quantiles for the neutral distribution from rcoal.R
# ci_rcoal <- fread(paste0(outdir_dat, "/ci_summary.csv"))[,alpha := 0.0][,s_set := NULL]
# 
# ci_summary <- rbind(ci_summary,ci_rcoal)

# Examine outliers --------------------------------------------------------




# Join neutral quantiles and define outliers
gene_pvals[ci_summary, on = c("alpha","freq"),
        pi_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

gene_pvals[ci_summary, on = c("alpha","freq"),
        pi_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]



# Plot --------------------------------------------------------------------

png(paste0(outdir_fig,"/demography.png"),
    width = 17, height = 10, units = "cm", res = 300,
    pointsize = 10, type = "cairo")

layout(matrix(c(1, 2), nrow = 1), widths = c(3, 2))
par(oma = c(0, 0, 1, 0), mar = c(4, 5, 1, 1))

# Panel 1 
# CI_area -----------------------------------------------------------------

max_y = 40

#theta
# with(gene_pvals[alpha==0],
plot(freq, seg_s_sites, pch = 16,#
     col = "grey40",
     xlim = c(0, 260),  
     ylim = c(0, max_y),
     cex = 0.45,
     yaxt = "n",
     ylab = "Number of synonymous segregating sites",
     xlab = expression("Pangenome frequency ("*italic(k)*")"))
# )
axis(side = 2, at = seq(0,max_y, length.out = 7),
     labels = sprintf("%3.0f", seq(0,max_y, length.out = 7)), las = 2)
alpha_vals <- sort(unique(ci_summary$alpha))
line_cols <- c(
  rgb(0.1,0.1,0.1, alpha = 0.5),   # grey        alpha = 0 (neutral reference)
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
  
  up <- sm(d$freq, d$upper_95)
  lo <- sm(d$freq, d$lower_95)
  
  # keep the original alpha = 0 as a filled band
  if (i == 0) {
    polygon(c(d$freq, rev(d$freq)), c(up, rev(lo)),
            col = line_cols[a], border = NA)
  }else{
  lines(d$freq, up, col = line_cols[a], lwd = 1.6, lty = 2)
  lines(d$freq, lo, col = line_cols[a], lwd = 1.6, lty = 2)
  }
}


legend("topleft", bty = "n",
       legend = sapply(alpha_vals, function(x) bquote(alpha == .(x))),
       col = line_cols, lty = c(NA,2,2,2), 
       pch = c(15, NA,NA,NA), lwd = 1.8)

mtext("a", side = 3, outer = TRUE,
      adj = 0, font = 2)  



# Panel 2 


# Robustness plot ---------------------------------------------------------

alpha_robustness <- gene_pvals[p_bh <=0.05, .(n=.N), by=c("alpha", "pi_sel")][!is.na(pi_sel)]
alpha_robustness <- dcast(alpha_robustness, pi_sel~alpha)
alpha_robustness <- alpha_robustness[order(-pi_sel)]
alpha_bp = as.matrix(alpha_robustness[,-1])
rownames(alpha_bp) = alpha_robustness$pi_sel

barplot(
  alpha_bp,
  col = c("#F26C2A", 
          rgb(0.2, 0.6, 1, alpha = 1)), 
  yaxt ="n",
  ylim = c(0,600),
  border = "white",
  space = 0.5,
  # font.axis = 2,
  ylab = "Number of outliers",
  xlab = expression(alpha~"value")
)
axis(side = 2, at = seq(0,600,100),
     labels = seq(0,600,100), las = 2)

legend("topright", 
       legend = c("High outlier (BH FDR < 0.05)","Low outlier (BH FDR < 0.05)"),
       col = c(rgb(0.2, 0.6, 1, alpha = 1),"#F26C2A"),
       pch = c(15,15),cex = 0.75,
       bty = "n")

mtext("b", side = 3, outer = TRUE,
      adj = 0.6, font = 2) 

dev.off()

layout(1) #reset






