# Import observed data ----------------------------------------------------

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

# Import core mutation rate -----------------------------------------------

core_mut_rate <- fread(paste0(outdir_dat, "/core_mut_rate.csv"))[3,]

# import m target

m_target <- as.numeric(readLines(paste0(outdir_dat, "/m_target.txt")))


# Scale trees to gene length (m) ------------------------------------------
# only use 500K trees
# file_out_name = paste0(outdir_dat, "/neutral_sim_260_ntree_1000000.csv")
# # take first occurance
# res <- fread(file_out_name)
# res_unique <- res[tree_id <= 500000]; rm(res)
# fwrite(res_unique, paste0(outdir_dat, "/res_unique.csv"))
# res_unique <- fread(paste0(outdir_dat, "/res_unique.csv"))

exp_0_tr_file = paste0(outdir_dat, "/expanding_sim_260_alpha0_ntree_500000.csv")
exp_0_tr <- fread(exp_0_tr_file, fill=TRUE)

exp_2.4_tr_file = paste0(outdir_dat, "/expanding_sim_260_alpha2.4_ntree_500000.csv")
exp_2.4_tr <- fread(exp_2.4_tr_file, fill=TRUE)

exp_5_tr_file = paste0(outdir_dat, "/expanding_sim_260_alpha5_ntree_500000.csv")
exp_5_tr <- fread(exp_5_tr_file, fill=TRUE)

exp_10_tr_file = paste0(outdir_dat, "/expanding_sim_260_alpha10_ntree_500000.csv")
exp_10_tr <- fread(exp_10_tr_file, fill=TRUE)

# calculate the expected number of segregating sites ----------------------
# we multiply the total branch length of each sub-tree 
# that produces x descendants by half our estimate of 2Nu from the core genes. 

pi_theta = (core_mut_rate$mean_pi_s)/2

exp_0_tr[, pi_exp_segsites := subtree_br_len * pi_theta]
exp_2.4_tr[, pi_exp_segsites := subtree_br_len * pi_theta]
exp_5_tr[, pi_exp_segsites := subtree_br_len * pi_theta]
exp_10_tr[, pi_exp_segsites := subtree_br_len * pi_theta]

# put on a scale of segregating sites per gene using median AG length, 700

exp_0_tr[, pi_exp_segsites := pi_exp_segsites * m_target]
exp_2.4_tr[, pi_exp_segsites := pi_exp_segsites * m_target]
exp_5_tr[, pi_exp_segsites := pi_exp_segsites * m_target]
exp_10_tr[, pi_exp_segsites := pi_exp_segsites * m_target]

# Then using this expectation we generate the probability of actually 
# observing 0, 1, 2….etc segregating sites using the Poisson distribution.

# Then we sum these probabilities across sub-trees to generate the final 
# distribution of the number of segregating sites. We then compare our 
# observations against this distribution.


# Function ----------------------------------------------------------------

# Want to be able to dynamically caluclate the neutral expectation 
# for summed subtrees simmulated with msprime under varaious demography schemes


get_exp_expectations <- function(exp_tr_dt){
  
  # get alpha
  nm <- deparse(substitute(exp_tr_dt))
  alpha <- suppressWarnings(as.numeric(sub("^exp_([0-9.]+)_tr$", "\\1", nm)))
  if (is.na(alpha)) stop("No alpha value found in name: ", nm)
  cat(paste(rep("-",50), collapse = ""), "\n","alpha = ", alpha, "\n")

  
  # Determine maximum number of segregating sites to consider
  lam_max  <- max(exp_tr_dt$pi_exp_segsites)
  site_rng <- 0:qpois(1 - 1e-10, lam_max)
  cat("site_rng now 0 to", max(site_rng),
      "| tail mass =", ppois(max(site_rng), lam_max, lower.tail = FALSE), "\n")
  
  n_s      <- length(site_rng)
  
  s_cols   <- paste0("s_", site_rng)
  
  n_k      <- max(exp_tr_dt$freq)
  
  
  ## === chunked accumulation=============================
  
  lam_v <- exp_tr_dt$pi_exp_segsites
  k_v   <- exp_tr_dt$freq
  
  accumulate_chunk <- function(idx) {
    p  <- vapply(site_rng, function(s) dpois(s, lam_v[idx]), numeric(length(idx)))
    rs <- rowsum(p, k_v[idx], reorder = TRUE)
    out <- matrix(0, n_k, n_s)
    out[as.integer(rownames(rs)), ] <- rs
    out
  }
  
  chunk    <- 2e5
  starts   <- seq(1, nrow(exp_tr_dt), by = chunk)
  idx_list <- lapply(starts, function(s) s:min(s + chunk - 1, nrow(exp_tr_dt)))
  
  cl <- makeCluster(5)
  on.exit(stopCluster(cl), add = TRUE)
  clusterExport(cl, c("lam_v", "k_v", "site_rng", "n_k", "n_s"),
                envir = environment())
  acc <- Reduce(`+`, parLapply(cl, idx_list, accumulate_chunk))

  ## === Downstream ===================================
  
  prob_pi <- data.table(freq = seq_len(n_k))
  prob_pi[, (s_cols) := as.data.table(acc)]
  prob_pi <- prob_pi[rowSums(acc) > 0]             # drop empty k classes
  setorderv(prob_pi, "freq", 1L)
  
  prob_pi[, row_total := rowSums(.SD), .SDcols = s_cols]
  prob_pi[, (s_cols) := lapply(.SD, function(x) x / row_total), .SDcols = s_cols]
  
  ## ACHTUNG!!: if site_rng does not span the bulk of the Poisson mass for the
  ## largest lambda, normalising by row_total silently rescales a truncated
  ## distribution and the quantiles below will be biased. Check first:
  ##   max(exp_tr_dt$pi_exp_segsites); range(site_rng)
  ## The tail mass beyond max(site_rng) should be negligible.
  cat(sprintf("max lambda = %.2f, max site_rng = %d, tail mass = %.3g\n",
              max(exp_tr_dt$pi_exp_segsites), max(site_rng),
              ppois(max(site_rng), max(exp_tr_dt$pi_exp_segsites),
                    lower.tail = FALSE)))
  
  prob_pi[, row_total := NULL]
  fwrite(prob_pi, file.path(outdir_dat, paste0("prob_pi_norm_alpha_",alpha,".csv")))
  
  pi_prod_lng <- melt(prob_pi, id.vars = "freq",
                      variable.name = "seg_sites", value.name = "p")
  pi_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]
  setorder(pi_prod_lng, freq, seg_sites)
  pi_prod_lng[, cdf := cumsum(p), by = freq]
  fwrite(pi_prod_lng, file.path(outdir_dat, paste0("pi_prod_lng_alpha_",alpha,".csv")))
  
  ## Quantiles. NOTE: which.max(cdf >= q) returns 1 when the condition is never
  ## met, silently giving the smallest seg_sites instead of flagging the problem.
  ## first_at() returns NA instead, so truncation shows up rather than hiding.
  first_at <- function(x, cdf, q) {
    i <- which(cdf >= q)
    if (length(i)) x[i[1L]] else NA_integer_
  }
  
  pi_ci_bounds <- pi_prod_lng[, .(
    lower_95 = first_at(seg_sites, cdf, 0.025),
    median_S = first_at(seg_sites, cdf, 0.500),
    upper_95 = first_at(seg_sites, cdf, 0.975),
    iqr      = first_at(seg_sites, cdf, 0.750) - first_at(seg_sites, cdf, 0.250)
  ), by = freq]
  pi_ci_bounds[, alpha := alpha]
  
  if (anyNA(pi_ci_bounds))
    message("Some quantiles unreachable -- site_rng is too narrow for large k.")
  
  # Outlier detections: summarise the neutral distribution per frequency --------
  
  # ci_summary <- rbind(pi_ci_bounds, theta_ci_bounds)
  
  return(pi_ci_bounds)
  
}


# implement neutral expectation function -----------------------------------

# 0
start.time <- Sys.time()

alpha_0_ci_bounds = get_exp_expectations(exp_0_tr)

Sys.time() - start.time
#Time difference of 10.99774 mins


# 2.4
start.time <- Sys.time()

alpha_2.4_ci_bounds = get_exp_expectations(exp_2.4_tr)

Sys.time() - start.time
#Time difference of 6.877202 mins

# 5
start.time <- Sys.time()

alpha_5_ci_bounds = get_exp_expectations(exp_5_tr)

Sys.time() - start.time
#Time difference of 6.377719 mins

# 10
start.time <- Sys.time()

alpha_10_ci_bounds = get_exp_expectations(exp_10_tr)

Sys.time() - start.time
#Time difference of 5.686743 mins


ci_summary <- rbind(alpha_0_ci_bounds, alpha_2.4_ci_bounds, alpha_5_ci_bounds, alpha_10_ci_bounds)

fwrite(ci_summary, paste0(outdir_dat, "/ci_summary_expanding_pop.csv"))

