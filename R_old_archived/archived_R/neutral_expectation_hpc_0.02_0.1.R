# Import observed data ----------------------------------------------------

# Import core mutation rate -----------------------------------------------
for (i in seq(0.01, 0.02,0.01)) {
  
  core_mut_rate = i
  
  cat("\n\nStarting with theta W = ", core_mut_rate,"\n")
  
  # import m target
  
  m_target = 435
  
  
  # Scale trees to gene length (m) ------------------------------------------
  # only use 500K trees
  res_unique <- fread(paste0(outdir_dat, "/res_unique.csv"))
  
  # calculate the expected number of segregating sites ----------------------
  # we multiply the total branch length of each sub-tree 
  # that produces x descendants by half our estimate of 2Nu from the core genes. 
  
  w_theta = (core_mut_rate)/2
  
  res_unique[, theta_exp_segsites := subtree_br_len * w_theta]
  
  # put on a scale of segregating sites per gene using median AG length, 700
  
  res_unique[, theta_exp_segsites := theta_exp_segsites * m_target]
  
  # Then using this expectation we generate the probability of actually 
  # observing 0, 1, 2….etc segregating sites using the Poisson distribution.
  
  # Then we sum these probabilities across sub-trees to generate the final 
  # distribution of the number of segregating sites. We then compare our 
  # observations against this distribution.
  
  # Determine maximum number of segregating sites to consider
  max_seg_sites <- round(max(res_unique$theta_exp_segsites) + 3) 
  
  site_rng <- 0:max_seg_sites
  
  s_cols <- paste0("s_", site_rng)
  
  
  # Using θw-based expectation ----------------------------------------------
  prob_theta <- matrix(0, nrow(res_unique), length(site_rng))
  
  start.time <- Sys.time()
  
  for (j in seq_along(site_rng)) {
    prob_theta[, j] <- dpois(site_rng[j], res_unique$theta_exp_segsites)
  }
  
  # get duration
  end.time <- Sys.time()
  
  cat("Loop took: ", end.time - start.time) # Time difference of 11.89873 mins
  
  colnames(prob_theta) = paste0("s_", site_rng)
  
  prob_theta <- as.data.frame(prob_theta)
  
  setDT(prob_theta)
  
  prob_theta <- cbind(res_unique[, .(freq,tree_id,theta_exp_segsites)], prob_theta)
  
  prob_theta <- prob_theta[, lapply(.SD, sum, na.rm = TRUE), .SDcols = s_cols, by = freq]
  
  # -------------------------------------------------------------------------
  
  #theta
  
  setorderv(prob_theta, cols = "freq", order = 1L)
  
  # Row sums for normalisation
  prob_theta[, row_total := rowSums(.SD), .SDcols = s_cols]
  
  # Normalise each column by row total
  prob_theta[, (s_cols) := lapply(.SD, function(x) x / row_total), 
             .SDcols = s_cols]
  
  prob_theta[, row_total:= NULL]
  
  # save normalised prob_theta
  fwrite(prob_theta, paste0("./output/data/theta_w_out/prob_theta_norm_",core_mut_rate,".csv"))
  # prob_theta <- fread(paste0(outdir_dat, "/prob_theta_norm.csv"))
  
  theta_prod_lng <- melt(prob_theta,
                         id.vars = "freq",
                         variable.name = "seg_sites", 
                         value.name = "p")
  
  # convert seg_sites to numeric
  theta_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]
  
  
  # Calculate cdf -----------------------------------------------------------
  # compute the cumulative distribution and 95% CI boundaries
  
  # Compute CDF per frequency class
  setorder(theta_prod_lng, freq, seg_sites)
  
  theta_prod_lng[, cdf := cumsum(p), by = freq]
  
  fwrite(theta_prod_lng, paste0("./output/data/theta_w_out/theta_prod_lng_",core_mut_rate,".csv"))
  
  # Extract 2.5% and 97.5% quantile boundaries per k
  theta_ci_bounds <- theta_prod_lng[, .(
    lower_95 = seg_sites[which.max(cdf >= 0.025)],  # 2.5th percentile
    median_S = seg_sites[which.max(cdf >= 0.500)],  # 50th percentile
    upper_95 = seg_sites[which.max(cdf >= 0.975)],   # 97.5th percentile
    iqr      = seg_sites[which.max(cdf >= 0.750)] - seg_sites[which.max(cdf >= 0.250)]
  ), by = freq]
  
  theta_ci_bounds$s_set = "w_theta"
  
  fwrite(theta_ci_bounds, paste0("./output/data/theta_w_out/ci_summary_",core_mut_rate,".csv"))
  
  cat("Done with theta W = ", core_mut_rate)
  
}



# Does theta cause outliers to reach neutral? -----------------------------
# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
ci_theta_list <- list.files(paste0(outdir_dat, "/theta_w_out"),
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
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_mutli_theta_fdr.csv"))

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# subset data
dt_set1 <- gene_pvals[gene_family %chin% set_1_ags$gene_family]
dt_set2 <- gene_pvals[gene_family %chin% set_2_ags$gene_family]


# Join neutral quantiles and define outliers
dt_set1[ci_summary[s_set=="w_theta"], on = c("theta_val", "freq"),
        theta_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

dt_set1[ci_summary[s_set=="w_theta"], on = c("theta_val", "freq"),
        w_theta_sel := fcase(
          seg_s_sites  < lower_95, "directional", 
          seg_s_sites  > upper_95, "balancing",
          default = NA
        )
]

dt_set2[ci_summary[s_set=="w_theta"], on = c("theta_val", "freq"),
        theta_outlier := fifelse(
          seg_s_sites  < lower_95 | seg_s_sites  > upper_95, 1L, 0L
        )
]

dt_set2[ci_summary[s_set=="w_theta"], on = c("theta_val", "freq"),
        w_theta_sel := fcase(
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
dt_set2$set = 2


plot_theta <- rbind(dt_set1, dt_set2)

theta_vals = sort(unique(plot_theta[!theta_val %chin% c(0.08,0.09)]$theta_val))

png(paste0(outdir_fig,"/thetaw_outliers_seg_sites_multi.png"),
    width = 210, height = 297, units = "mm", res = 300,
    pointsize = 12, type = "cairo")

par(oma=c(4,2,2,0)) # Set outer margins

mat <- matrix(1:14, nrow = 7, byrow = TRUE)
layout(mat)
#layout.show(n=14)

# Wattersons' theta neutral distribution
for (t in theta_vals) {
  
  for (i in 1:2) {
    
  if(i==1){
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
       xaxt = ifelse(t == 0.07, "s","n"),
       ylab = "",
       xlab = "")
  
  line_vals_w_theta <- ci_summary[theta_val == t][order(freq)]
  
  polygon(
    c(line_vals_w_theta$freq, rev(line_vals_w_theta$freq)),
    c(line_vals_w_theta$upper_95, rev(line_vals_w_theta$lower_95)),
    col = "grey80",
    border = "grey80"
  )
  
  with(plot_theta[theta_val == t & set == i], 
       points(freq, seg_s_sites, pch = 16,
              ylim = c(0, 1),
              cex = 0.55,
              yaxt = "n", bty="L",
              ylab = "",
              xlab = "",
              col = ifelse(theta_outlier ==1 & w_theta_sel == "balancing" & p_bh <= 0.05,
                           rgb(0.2, 0.6, 1, alpha = 0.5),
                           ifelse(theta_outlier ==1 & w_theta_sel != "balancing" & p_bh <= 0.05,"#F26C2A", 
                                  rgb(0.1,0.1,0.1, alpha = 0.5))
              )))
  if(i==1){
    axis(side = 2, at = seq(0,25, length.out = 7),
         labels = sprintf("%3.0f", seq(0,25, length.out = 7)), las = 2)
    }

  perc_AGs_out = nrow(plot_theta[theta_val == t & set == i][theta_outlier==1 & p_bh <= 0.05])/nrow(plot_theta[theta_val == t & set == i])*100
  
  mtext(bquote(.(nrow(plot_theta[theta_val == t & set == i][theta_outlier==1& p_bh <= 0.05]))*" ("*.(sprintf("%2.2f",perc_AGs_out))*"%) AGs" *
                 ","~theta[W]~"="~.(t)), 
        side = 3, line = -1.75, cex = 0.75)
}
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

mtext(expression("Increasing"~theta[W]),
      outer = TRUE, side = 3, line = 0, cex = 1.1)

dev.off()

















