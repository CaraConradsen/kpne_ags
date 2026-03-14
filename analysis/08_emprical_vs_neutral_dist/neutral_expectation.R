# Import observed data ----------------------------------------------------

# from get_ag_pis.R
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

# remove outliers
ag_age_S_dt <- ag_age_S_dt[pi_out==0][, pi_out:=NULL]


# Scale trees -------------------------------------------------------------
# res <- fread(paste0(outdir_dat,"/neutral_sim_260_ntree_10000.csv"))

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# Calculate empirical mean for frequency class k
pi_emp_k_set_1 = ag_age_S_dt[gene_family %chin% set_1_ags$gene_family, .(mean_Sm = mean(Sm_std)), by = freq]
pi_emp_k_set_2 = ag_age_S_dt[gene_family %chin% set_2_ags$gene_family, .(mean_Sm = mean(Sm_std)), by = freq]

# Anchor genealogy estimates by empirical means ---------------------------
res_1 <- copy(res)
res_2 <- copy(res)

# Join pi_emp_k onto res by freq and create new column
res_1[pi_emp_k_set_1, on = "freq", pi_sim := t_hat * mean_Sm]
res_2[pi_emp_k_set_2, on = "freq", pi_sim := t_hat * mean_Sm]

# Outlier detections: summarise the neutral distribution per frequency --------

neutral_summary_1 <- res_1[pi_sim <= 1,.(
  pi_lower = quantile(pi_sim, 0.025, na.rm = TRUE),
  pi_median = median(pi_sim, na.rm = TRUE),
  pi_upper = quantile(pi_sim, 0.975, na.rm = TRUE),
  iqr = IQR(pi_sim, na.rm = TRUE, type = 7)
), by = freq]

neutral_summary_1$set = 1

neutral_summary_2 <- res_2[pi_sim <= 1,.(
  pi_lower = quantile(pi_sim, 0.025, na.rm = TRUE),
  pi_median = median(pi_sim, na.rm = TRUE),
  pi_upper = quantile(pi_sim, 0.975, na.rm = TRUE),
  iqr = IQR(pi_sim, na.rm = TRUE, type = 7)
), by = freq]

neutral_summary_2$set = 2

neutral_summary <- rbind(neutral_summary_1[!is.na(pi_median)],
                         neutral_summary_2[!is.na(pi_median)])

fwrite(neutral_summary, paste0(outdir_dat, "/neutral_summary.csv"))
