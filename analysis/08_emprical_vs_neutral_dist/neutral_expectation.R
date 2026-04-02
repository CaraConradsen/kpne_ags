# Import observed data ----------------------------------------------------

# from get_ag_pis.R
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

# remove outliers
ag_age_S_dt <- ag_age_S_dt[pi_out==0][, pi_out:=NULL]


# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# Import core mutation rate -----------------------------------------------

core_mut_rate <- fread(paste0(outdir_dat, "/core_mut_rate.csv"))


# Scale trees to gene length (m) ------------------------------------------
file_out_name = paste0(outdir_dat, "/neutral_sim_260_ntree_1000000.csv")
res <- fread(file_out_name)


# Reduce res size ---------------------------------------------------------
# only use 500K trees

# take first occurance
res_unique <- res[tree_id <= 500000]; rm(res)

fwrite(res_unique, paste0(outdir_fig, "/res_unique.csv"))
# res_unique <- fread(paste0(outdir_fig, "/res_unique.csv"))

# calculate the expected number of segregating sites ----------------------
# we multiply the total branch length of each sub-tree 
# that produces x descendants by half our estimate of 2Nu from the core genes. 

pi_theta = (core_mut_rate$mean_pi)/2
w_theta = (core_mut_rate$mean_theta_w)/2

res_unique[, pi_exp_segsites := subtree_br_len * pi_theta]
res_unique[, theta_exp_segsites := subtree_br_len * w_theta]

# put on a scale of segregating sites per gene using median AG length, 700

res_unique[, pi_exp_segsites := pi_exp_segsites * 700]
res_unique[, theta_exp_segsites := theta_exp_segsites * 700]

# Then using this expectation we generate the probability of actually 
# observing 0, 1, 2….etc segregating sites using the Poisson distribution.

# Then we sum these probabilities across sub-trees to generate the final 
# distribution of the number of segregating sites. We then compare our 
# observations against this distribution.

# Determine maximum number of segregating sites to consider
max_seg_sites <- round(max(res_unique$pi_exp_segsites , res_unique$theta_exp_segsites) + 1) 

site_rng <- 0:max_seg_sites

s_cols <- paste0("s_", site_rng)

# Look at pi expectation first ----------------------------------------------

prob_pi <- matrix(0, nrow(res_unique), length(site_rng))

start.time <- Sys.time()

for (j in seq_along(site_rng)) {
  prob_pi[, j] <- dpois(site_rng[j], res_unique$pi_exp_segsites)
}

# get duration
end.time <- Sys.time()
end.time - start.time # Time difference of 8.102087 mins

colnames(prob_pi) = paste0("s_", site_rng)

setDT(as.data.frame(prob_pi))

prob_pi <- cbind(res_unique[, .(freq,tree_id,pi_exp_segsites)], prob_pi)

prob_pi <- prob_pi[, lapply(.SD, sum, na.rm = TRUE), .SDcols = s_cols, by = freq]


fwrite(prob_pi, paste0(outdir_dat, "/prob_pi.csv"))
# prob_pi <- fread(paste0(outdir_dat, "/prob_pi.csv"))


# Using θw-based expectation ----------------------------------------------
prob_theta <- matrix(0, nrow(res_unique), length(site_rng))

start.time <- Sys.time()

for (j in seq_along(site_rng)) {
  prob_theta[, j] <- dpois(site_rng[j], res_unique$theta_exp_segsites)
}

# get duration
end.time <- Sys.time()
end.time - start.time # Time difference of 8.102087 mins

colnames(prob_theta) = paste0("s_", site_rng)

prob_theta <- as.data.frame(prob_theta)

setDT(prob_theta)

prob_theta <- cbind(res_unique[, .(freq,tree_id,theta_exp_segsites)], prob_theta)

prob_theta <- prob_theta[, lapply(.SD, sum, na.rm = TRUE), .SDcols = s_cols, by = freq]

fwrite(prob_theta, paste0(outdir_dat, "/prob_theta.csv"))
# prob_theta <- fread(paste0(outdir_dat, "/prob_theta.csv"))


# -------------------------------------------------------------------------


setorderv(prob_pi, cols = "freq", order = 1L)

# Row sums for normalisation
prob_pi[, row_total := rowSums(.SD), .SDcols = s_cols]

# Normalise each column by row total
prob_pi[, (s_cols) := lapply(.SD, function(x) x / row_total), 
        .SDcols = s_cols]

prob_pi[, row_total:= NULL]

pi_prod_lng <- melt(prob_pi,
                    id.vars = "freq",
                    variable.name = "seg_sites", 
                    value.name = "p")

# convert seg_sites to numeric
pi_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]


#theta

setorderv(prob_theta, cols = "freq", order = 1L)

# Row sums for normalisation
prob_theta[, row_total := rowSums(.SD), .SDcols = s_cols]

# Normalise each column by row total
prob_theta[, (s_cols) := lapply(.SD, function(x) x / row_total), 
            .SDcols = s_cols]

prob_theta[, row_total:= NULL]

theta_prod_lng <- melt(prob_theta,
                    id.vars = "freq",
                    variable.name = "seg_sites", 
                    value.name = "p")

# convert seg_sites to numeric
theta_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]


# Calculate cdf -----------------------------------------------------------
# compute the cumulative distribution and 95% CI boundaries

# Compute CDF per frequency class
setorder(pi_prod_lng, freq, seg_sites)
setorder(theta_prod_lng, freq, seg_sites)

pi_prod_lng[, cdf := cumsum(p), by = freq]
theta_prod_lng[, cdf := cumsum(p), by = freq]

# Extract 2.5% and 97.5% quantile boundaries per k
pi_ci_bounds <- pi_prod_lng[, .(
  lower_95 = seg_sites[which.max(cdf >= 0.025)],  # 2.5th percentile
  median_S = seg_sites[which.max(cdf >= 0.500)],  # 50th percentile
  upper_95 = seg_sites[which.max(cdf >= 0.975)],   # 97.5th percentile
  iqr      = seg_sites[which.max(cdf >= 0.750)] - seg_sites[which.max(cdf >= 0.250)]
), by = freq]

pi_ci_bounds$s_set = "pi"

theta_ci_bounds <- theta_prod_lng[, .(
  lower_95 = seg_sites[which.max(cdf >= 0.025)],  # 2.5th percentile
  median_S = seg_sites[which.max(cdf >= 0.500)],  # 50th percentile
  upper_95 = seg_sites[which.max(cdf >= 0.975)],   # 97.5th percentile
  iqr      = seg_sites[which.max(cdf >= 0.750)] - seg_sites[which.max(cdf >= 0.250)]
), by = freq]

theta_ci_bounds$s_set = "w_theta"

# Outlier detections: summarise the neutral distribution per frequency --------

ci_summary <- rbind(pi_ci_bounds, theta_ci_bounds)

fwrite(ci_summary, paste0(outdir_dat, "/ci_summary.csv"))
