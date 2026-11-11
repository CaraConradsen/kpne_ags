#!/usr/bin/env Rscript

# Sensitivity of outlier detection to theta (Figure S5)
#
# Same machinery as neutral_expectation_hpc.R, but holding gene length fixed
# and sweeping theta instead of sweeping gene length at the fitted theta.
#
# Three things differ from the earlier version of this script:
#   1. It reads the EXPANDING genealogies (alpha_hat = 2.89), not the
#      constant-size set. The earlier script scaled res_unique.csv, which
#      was simulated under alpha = 0, so its intervals were the wrong shape
#      and the outlier counts were not comparable with Figure 2.
#   2. Quantiles use first_at(), which returns NA when a quantile is not
#      reached, rather than which.max(cdf >= q), which silently returns the
#      first element and hides a truncated site range.
#   3. Probabilities accumulate in chunks with rowsum() rather than building
#      one nrow(res_unique) x n_sites matrix, which does not fit in memory
#      at 500,000 trees.
#
# One array task per theta value.

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

# Array bookkeeping -------------------------------------------------------

outdir_dat <- Sys.getenv("OUTDIR_DAT", unset = ".")
n_cores    <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))

# Gene length at which the sensitivity is shown, in synonymous sites. The
# median over the syntenic AGs is 400; the earlier script used 435, which
# was the median of an earlier AG set.
m_target   <- as.integer(Sys.getenv("M_TARGET", unset = "400"))

# theta grid. These are INPUT theta per site, not pi_S: under expansion the
# fitted input theta_hat = 0.031 reproduces a core pi_S of 0.0121 and a core
# theta_W of 0.0195 (Figure S4), so the panels should be labelled by input
# theta and the corresponding pi_S stated in the legend.
theta_vals <- as.numeric(strsplit(
  Sys.getenv("THETA_GRID", unset = "0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10"),
  ",")[[1]])

task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
if (is.na(task)) stop("SLURM_ARRAY_TASK_ID not set")
if (task > length(theta_vals)) { cat("nothing to do\n"); quit(save = "no") }

theta_hat <- theta_vals[task]

outdir <- file.path(outdir_dat, "theta_grid_out")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("task", task, "| theta =", theta_hat, "| m =", m_target, "\n")

# Import the expanding genealogies ----------------------------------------
# alpha_hat = 2.89, the value at which the simulated pi_S / theta_W ratio
# matches the observed core ratio of 0.62 (Figure S4a).

tree_file <- file.path(outdir_dat,
                       "expanding_sim_260_alpha2.89_ntree_500000.csv")

res_unique <- fread(tree_file, fill = TRUE,
                    select = c("tree_id", "freq", "subtree_br_len"))[!is.na(tree_id)]
res_unique[, tree_id := NULL]

# Expected segregating sites ----------------------------------------------
# Total branch length of each sub-tree producing k descendants, scaled by
# half the mutation rate, then by the number of synonymous sites.

pi_theta <- theta_hat / 2

res_unique[, pi_exp_segsites := subtree_br_len * pi_theta * m_target]

n_k <- max(res_unique$freq)

start.time <- Sys.time()

lam_max  <- max(res_unique$pi_exp_segsites)
site_rng <- 0:qpois(1 - 1e-10, lam_max)
n_s      <- length(site_rng)
s_cols   <- paste0("s_", site_rng)

cat("site_rng 0 to", max(site_rng),
    "| tail mass =", ppois(max(site_rng), lam_max, lower.tail = FALSE), "\n")

# Chunked accumulation ----------------------------------------------------

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

acc <- Reduce(`+`, mclapply(idx_list, accumulate_chunk, mc.cores = n_cores))

print(Sys.time() - start.time)

# Normalise ---------------------------------------------------------------

prob_pi <- data.table(freq = seq_len(n_k))
prob_pi[, (s_cols) := as.data.table(acc)]
prob_pi <- prob_pi[rowSums(acc) > 0]
setorderv(prob_pi, "freq", 1L)

prob_pi[, row_total := rowSums(.SD), .SDcols = s_cols]
prob_pi[, (s_cols) := lapply(.SD, function(x) x / row_total), .SDcols = s_cols]

## ACHTUNG!!: if site_rng does not span the bulk of the Poisson mass for the
## largest lambda, normalising by row_total silently rescales a truncated
## distribution and the quantiles below will be biased. This matters more
## here than in the main run, because the largest theta on the grid is three
## times the fitted value and lambda scales with it.
cat(sprintf("max lambda = %.2f, max site_rng = %d, tail mass = %.3g\n",
            lam_max, max(site_rng),
            ppois(max(site_rng), lam_max, lower.tail = FALSE)))

prob_pi[, row_total := NULL]
fwrite(prob_pi, file.path(outdir,
                          sprintf("prob_pi_norm_theta%s.csv", theta_hat)))

pi_prod_lng <- melt(prob_pi, id.vars = "freq",
                    variable.name = "seg_sites", value.name = "p")
pi_prod_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]
setorder(pi_prod_lng, freq, seg_sites)
pi_prod_lng[, cdf := cumsum(p), by = freq]
fwrite(pi_prod_lng, file.path(outdir,
                              sprintf("pi_prod_lng_theta%s.csv", theta_hat)))

# Quantiles ---------------------------------------------------------------
## NOTE: which.max(cdf >= q) returns 1 when the condition is never met,
## silently giving the smallest seg_sites instead of flagging the problem.
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
pi_ci_bounds[, theta := theta_hat]
pi_ci_bounds[, m     := m_target]

if (anyNA(pi_ci_bounds))
  warning(sprintf("theta = %s: some quantiles unreachable -- site_rng too narrow for large k.",
                  theta_hat))

fwrite(pi_ci_bounds, file.path(outdir,
                               sprintf("ci_summary_theta%s.csv", theta_hat)))

cat("\ntask", task, "done, theta =", theta_hat, "\n")
