
# the loglikelihood of growth rate alpha is the sum of all gene 
# probabilities, given the gene sits in k with growth rate alpha,
# how likley are we to see s segregating sites
# ℓ(α) = ∑[g] log P_hat (S=sg | kg, α)



# Get AG observed segregating sites per frequency -------------------------
ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

colnames(ag_age_S_dt)[10] = "freq"

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_s_sites),
                               freq = mean(freq)), gene_family]


# get piS theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]

k_obs = ag_age_S_dt$freq
s_obs = ag_age_S_dt$seg_s_sites

# Read in lookup tables ---------------------------------------------------


list_gene_prob_tables <- list.files(outdir_dat, pattern = "^prob_pi_norm_alpha_\\d+(\\.\\d+)?\\.csv$",
                                    full.names = TRUE, recursive = TRUE)

alphas <- gsub("^prob_pi_norm_alpha_|\\.csv$", "", basename(list_gene_prob_tables))

pmf_list <- setNames(
  lapply(list_gene_prob_tables, function(pmf_file) {
    pmf_dt <- fread(pmf_file)
    pmf <- as.matrix(pmf_dt[, -1])
    rownames(pmf) <- pmf_dt$freq
    pmf
  }),
  alphas
)

# pmf is a matrix; rows indexed by k, columns by s (starting at s = 0)
loglik <- function(pmf, k_obs, s_obs) {
  p <- pmf[cbind(k_obs - 1, s_obs + 1)]   # +1 because s starts at 0, at k starts at 2
  sum(log(p))
}

ll <- sapply(alphas, function(a) loglik(pmf_list[[a]], k_obs, s_obs))

loglike_dt <- data.table(loglike = ll, alpha = names(ll))
