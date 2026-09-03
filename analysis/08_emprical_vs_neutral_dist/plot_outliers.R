
# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
ci_summary <- fread(paste0(outdir_dat, "/ci_summary.csv"))

# fdr corrected p-values
# Pool all p-values, one global BH correction
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

# subset data
dt_set1 <- gene_pvals[gene_family %chin% set_1_ags$gene_family]


# Join neutral quantiles and define outliers
dt_set1[ci_summary[s_set=="pi"], on = "freq",
          pi_outlier := fifelse(
            seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
          )
]

dt_set1[ci_summary[s_set=="pi"], on = "freq",
        pi_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]

plot_dat_list <- list(dt_set1)

fwrite(cbind(dt_set1, set=1),
       na = NA,
       paste0(outdir_dat, "/set_outliers.csv"))

# Plot graph --------------------------------------------------------------
# nice light blue: 0.400, 0.714, 0.922
# 
png(paste0(outdir_fig,"/theta_piS_outliers_seg_sites.png"),
    width = 25.7, height = 14.2875, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

par(oma=c(0.5,2,2,0.5)) # Set outer margins

# mat <- matrix(c(
#   rep(1,5),rep(2,6),
#   rep(3,5),
#   rep(4,4),
#   rep(c(5,6,7,8), each = 3)), nrow = 2, byrow = TRUE)
# layout(mat, heights = c(0.55, 0.45))
# #layout.show(n=8)

# SEVEN LAYOUT
mat <- matrix(c(rep(1,5),
                rep(2,5),
                rep(3:7, each = 2)), nrow = 2, byrow = TRUE)
layout(mat, heights = c(0.55, 0.45))
#layout.show(n=7)

max_y = 35

# barplot of retained AGs
source("./analysis/07_hgt_synteny_ag_sets/barplot_ag_filtering.R")
mtext("a", side = 3, outer = TRUE,
      adj = 0, font = 2)



# The results! ------------------------------------------------------------

sm <- function(x, y, span = 0.2) {
  predict(loess(y ~ x, span = span), newdata = data.frame(x = x))
}

# theta neutral distribution (constant population)
for (i in 1) {#seq_along(plot_dat_list)
  
  if(i==1){
    par(mar=c(4,6,1,0))
  }else{
    par(mar=c(4,1,1,0.5))
  }
  
  #theta
  plot(NULL, pch = 16,
       xlim = c(0, 260),  
       ylim = c(0, max_y),
       cex = 0.45,
       yaxt = "n",
       ylab = "",
       xlab = "")
  
  line_vals_w_theta <- ci_summary[s_set== "pi"][order(freq)]
  
  up <- sm(line_vals_w_theta$freq, line_vals_w_theta$upper_95)
  lo <- sm(line_vals_w_theta$freq, line_vals_w_theta$lower_95)
  
  polygon(
    c(line_vals_w_theta$freq, rev(line_vals_w_theta$freq)),
    c(up, rev(lo)),
    col = "grey80",
    border = "grey80"
  )

  with(plot_dat_list[[i]], 
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
  if(i==1){
    axis(side = 2, at = seq(0,max_y, length.out = 7),
         labels = sprintf("%3.0f", seq(0,max_y, length.out = 7)), las = 2)
    axis(side = 2, at = max_y/2, tick = FALSE,line = 1.5,
         labels = "Number of synonymous segregating sites")
    mtext("b", side = 3, outer = TRUE,
          adj = 0.5, font = 2)
  }else{
    axis(side = 1, at = -20, tick = FALSE,line = 1.5,
         labels = expression("Pangenome frequency ("*italic(k)*")"), xpd = TRUE)
    # mtext("c", side = 3, outer = TRUE,
    #       adj = 0.69, font = 2)
  }
  
  perc_AGs_out = nrow(plot_dat_list[[i]][pi_outlier==1 & p_bh <= 0.05])/nrow(plot_dat_list[[i]])*100
  
  # mtext(bquote(atop("Set "~.(i),
  #                   .(nrow(plot_dat_list[[i]][theta_outlier==1& p_bh <= 0.05]))*" ("*.(sprintf("%2.2f",perc_AGs_out))*"%) AGs")), 
  #       side = 3, line = -1.75, cex = 0.75)
  
  mtext(bquote(atop(.(nrow(plot_dat_list[[i]][pi_outlier==1& p_bh <= 0.05]))*" ("*.(sprintf("%2.2f",perc_AGs_out))*"%) AGs",
                    "")), 
        side = 3, line = -1.75, cex = 0.75)
  
  spearmans <- with(plot_dat_list[[i]], suppressWarnings(cor.test(freq, seg_s_sites, method = "spearman")))
  
  rho_val <- sprintf("%.2f", spearmans$estimate)
  
  p_part <- if (spearmans$p.value < 0.001) {
    quote(italic(P) < 0.001)
  } else {
    bquote(italic(P) == .(sprintf("%.3f", spearmans$p.value)))
  }
  
  mtext(bquote(rho == .(rho_val) * ";" ~ .(p_part)),
        side = 3, adj = 0.9, line = -2, cex = 0.8)
  
}

legend(-5,40, #horiz = TRUE,
       legend = c("","","High outlier (BH FDR < 0.05)","Low outlier (BH FDR < 0.05)", "Non-significant", "Neutral 95% CI"),
       col = c("white","white",rgb(0.2, 0.6, 1, alpha = 1),"#F26C2A", "black", "grey80"),
       pch = c(16, 16, 16, 16,16,15),
       bty = "n")

# add in slices
# source("./analysis/08_emprical_vs_neutral_dist/neut_w_theta_probs_sliced.R")
source("./analysis/08_emprical_vs_neutral_dist/neut_piS_theta_probs_sliced.R")

mtext("c", side = 3, outer = TRUE,
      adj = 0, font = 2, line = -25)


dev.off()

# 
# dat <- merge(unique(anno[,.(gene_family, number_genomes, consensus_gene_name, consensus_product, alleles_at_maximum_threshold)]),
#              plot_dat_list[[i]][, .(gene_family,freq, seg_s_sites)],
#              all.x = TRUE, by = "gene_family")
# dat[4,3] = "AbiEii/AbiGii"
# dat[7,3] = "AbiEi"
# 
# with(dat,
#      points(freq, seg_s_sites, pch = 1, cex = 1.5, col= "red", lwd = 2))
# with(dat,
#      text(freq, seg_s_sites, consensus_gene_name, cex = 0.75, pos = 4, col = "red"))
# 
# 

