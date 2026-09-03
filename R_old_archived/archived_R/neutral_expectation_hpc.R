library(data.table)
library(foreach)
library(doParallel)

setDTthreads(1) 

# CPUs from SLURM, not hardcoded
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
outdir_dat = "./ag_expected_sites_per_len"

# Import observed data ----------------------------------------------------

# Import core mutation rate -----------------------------------------------
# from alpha = 3.395; expected value of pi_branch = 0.7258, therefore 

theta_hat  = 0.033
m = x# needs to read in, not sure how

# Scale trees to gene length (m) ------------------------------------------
# only use 500K trees

file_out_name = paste0(outdir_dat, "/expanding_sim_260_alpha3.395_ntree_500000.csv")
res_unique <- fread(file_out_name, fill=TRUE)[!is.na(tree_id)]

# calculate the expected number of segregating sites ----------------------
# we multiply the total branch length of each sub-tree 
# that produces x descendants by half our estimate of 2Nu from the core genes. 

pi_theta = (theta_hat)/2

res_unique[, pi_exp_segsites := subtree_br_len * pi_theta]

# put on a scale of segregating sites per gene using median AG length, 700

res_unique[, pi_exp_segsites := pi_exp_segsites * m]

# Then using this expectation we generate the probability of actually 
# observing 0, 1, 2….etc segregating sites using the Poisson distribution.

# Then we sum these probabilities across sub-trees to generate the final 
# distribution of the number of segregating sites. We then compare our 
# observations against this distribution.

# Determine maximum number of segregating sites to consider
lam_max  <- max(res_unique$pi_exp_segsites)
site_rng <- 0:qpois(1 - 1e-10, lam_max)
cat("site_rng now 0 to", max(site_rng),
    "| tail mass =", ppois(max(site_rng), lam_max, lower.tail = FALSE), "\n")

n_s      <- length(site_rng)

s_cols   <- paste0("s_", site_rng)

n_k      <- max(res_unique$freq)


## === chunked accumulation=============================

start.time <- Sys.time()

accumulate_chunk <- function(idx) {
  lam <- res_unique$pi_exp_segsites[idx]
  k   <- res_unique$freq[idx]
  p   <- vapply(site_rng, function(s) dpois(s, lam), numeric(length(lam)))
  rs  <- rowsum(p, k, reorder = TRUE)
  out <- matrix(0, n_k, n_s)
  out[as.integer(rownames(rs)), ] <- rs
  out
}

chunk    <- 2e5
starts   <- seq(1, nrow(res_unique), by = chunk)
idx_list <- lapply(starts, function(s) s:min(s + chunk - 1, nrow(res_unique)))

cl <- makeCluster(6)
clusterExport(cl, c("res_unique", "site_rng", "n_k", "n_s", "accumulate_chunk"))
acc <- Reduce(`+`, parLapply(cl, idx_list, accumulate_chunk))
stopCluster(cl)

Sys.time() - start.time
#Time difference of 9.573947 mins

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
##   max(res_unique$pi_exp_segsites); range(site_rng)
## The tail mass beyond max(site_rng) should be negligible.
cat(sprintf("max lambda = %.2f, max site_rng = %d, tail mass = %.3g\n",
            max(res_unique$pi_exp_segsites), max(site_rng),
            ppois(max(site_rng), max(res_unique$pi_exp_segsites),
                  lower.tail = FALSE)))

prob_pi[, row_total := NULL]
fwrite(prob_pi, file.path(outdir_dat, "prob_pi_norm_",m,".csv"))

pi_prod_lng <- melt(prob_pi, id.vars = "freq",
                    variable.name = "seg_sites", value.name = "p")
pi_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]
setorder(pi_prod_lng, freq, seg_sites)
pi_prod_lng[, cdf := cumsum(p), by = freq]
fwrite(pi_prod_lng, file.path(outdir_dat, "pi_prod_lng_",m,".csv"))

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
pi_ci_bounds[, s_set := "pi"]

if (anyNA(pi_ci_bounds))
  warning("Some quantiles unreachable -- site_rng is too narrow for large k.")

# Outlier detections: summarise the neutral distribution per frequency --------

# ci_summary <- rbind(pi_ci_bounds, theta_ci_bounds)

fwrite(pi_ci_bounds, paste0(outdir_dat, "/ci_summary_",m,".csv"))
