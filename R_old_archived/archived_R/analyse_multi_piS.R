

# Does theta cause outliers to reach neutral? -----------------------------
# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
ci_theta_list <- list.files(paste0(outdir_dat,"/theta_pis_out"),
                            pattern ="ci_summary", 
                            full.names = TRUE, recursive = TRUE)


ci_summary_list <- lapply(ci_theta_list, function(lst){
  ci_dat <- fread(lst)
  ci_dat$theta_val <- as.numeric(gsub("ci_summary_", "", 
                                      gsub(".csv", "", basename(lst))))
  
  ci_dat
})

ci_summary <- rbindlist(ci_summary_list)


# fdr corrected p-values
# Pool all p-values, one global BH correction
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_mutli_pis_fdr.csv"))

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# subset data
dt_set1 <- gene_pvals[gene_family %chin% set_1_ags$gene_family]
dt_set2 <- gene_pvals[gene_family %chin% set_2_ags$gene_family]


# Join neutral quantiles and define outliers
dt_set1[ci_summary[s_set=="piS_theta"], on = c("theta_val", "freq"),
        theta_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

dt_set1[ci_summary[s_set=="piS_theta"], on = c("theta_val", "freq"),
        piS_theta_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]

dt_set2[ci_summary[s_set=="piS_theta"], on = c("theta_val", "freq"),
        theta_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

dt_set2[ci_summary[s_set=="piS_theta"], on = c("theta_val", "freq"),
        piS_theta_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]



# find central theta ------------------------------------------------------
theta_vals = sort(unique(dt_set1$theta_val))

dot_spread <- lapply(theta_vals, function(theta){
  fit <- with(dt_set1[theta_val == theta],
              lm(expected_S~freq))
  
  resid <- dt_set1[theta_val == theta, seg_s_sites] - predict(fit)      # fit = your lm object
  n_above <- sum(resid > 0)
  n_below <- sum(resid < 0)
  n_on    <- sum(resid == 0)         # exact ties, usually rare
  
  c(theta, n_above, n_on, n_below, (n_above/n_below))
  
})

dot_spread <- do.call(rbind, dot_spread)
dot_spread




# Plot data ---------------------------------------------------------------

dt_set1$set = 1

plot_theta <- rbind(dt_set1)

theta_vals = sort(unique(plot_theta$theta_val))

png(paste0(outdir_fig,"/pis_outliers_seg_sites_multi_set1.png"),
    width = 210, height = 297, units = "mm", res = 300,
    pointsize = 12, type = "cairo")

par(oma=c(4,2,2,0)) # Set outer margins

mat <- matrix(1:10, nrow = 5, byrow = TRUE)
layout(mat)
#layout.show(n=14)

# Wattersons' theta neutral distribution
for (t in theta_vals) {
  
    if(round(t, 2) %in% round(seq(0.01, 0.1, 0.02), 2)){
      par(mar=c(0,4,0,0))
    }else{
      par(mar=c(0,1,0,0.5))
    }
    
    #theta
    plot(NULL, pch = 16,
         xlim = c(0, 260),  
         ylim = c(0, 25),
         cex = 0.45,
         yaxt = "n",xpd = TRUE,
         xaxt = ifelse(t %in% c(0.09,0.1), "s","n"),
         ylab = "",
         xlab = "")
    
    line_vals_piS_theta <- ci_summary[theta_val == t][order(freq)]
    
    polygon(
      c(line_vals_piS_theta$freq, rev(line_vals_piS_theta$freq)),
      c(line_vals_piS_theta$upper_95, rev(line_vals_piS_theta$lower_95)),
      col = "grey80",
      border = "grey80"
    )
    
    with(plot_theta[theta_val == t], 
         points(freq, seg_s_sites, pch = 16,
                ylim = c(0, 1),
                cex = 0.55,
                yaxt = "n", bty="L",
                ylab = "",
                xlab = "",
                col = ifelse(theta_outlier ==1 & piS_theta_sel == "balancing" & p_bh <= 0.05,
                             rgb(0.2, 0.6, 1, alpha = 0.5),
                             ifelse(theta_outlier ==1 & piS_theta_sel != "balancing" & p_bh <= 0.05,"#F26C2A", 
                                    rgb(0.1,0.1,0.1, alpha = 0.5))
                )))
    
    if(round(t, 2) %in% round(seq(0.01, 0.1, 0.02), 2)){
      axis(side = 2, at = seq(0,25, length.out = 7),
           labels = sprintf("%3.0f", seq(0,25, length.out = 7)), las = 2)
    }
    
    perc_AGs_out = nrow(plot_theta[theta_val == t][theta_outlier==1 & p_bh <= 0.05])/nrow(plot_theta[theta_val == t])*100
    
    mtext(bquote(.(nrow(plot_theta[theta_val == t][theta_outlier==1& p_bh <= 0.05]))*" ("*.(sprintf("%2.2f",perc_AGs_out))*"%) AGs" *
                   ","~pi[S]~"="~.(t)), 
          side = 3, line = -1.75, cex = 0.75)
  }

legend("topleft", #horiz = TRUE,
       legend = c("","","High outlier (BH FDR < 0.05)","Low outlier (BH FDR < 0.05)", "Non-significant", "Neutral 95% CI"),
       col = c("white","white",rgb(0.2, 0.6, 1, alpha = 1),"#F26C2A", "black", "grey80"),
       pch = c(16, 16, 16, 16,16,15),
       bty = "n")

mtext("Number of synonymous segregating sites",
      outer = TRUE, side = 2)

mtext(expression("Pangenome frequency ("*italic(k)*")"),
      outer = TRUE, side = 1, line = 3)

mtext(expression("Increasing"~pi[S]),
      outer = TRUE, side = 3, line = 0, cex = 1.1)

dev.off()

















