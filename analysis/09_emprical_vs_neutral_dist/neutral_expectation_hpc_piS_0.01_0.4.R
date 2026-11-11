# Import observed data ----------------------------------------------------
library(data.table)
library(foreach)
library(doParallel)

setDTthreads(1) 

# CPUs from SLURM, not hardcoded
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

res_unique_master <- fread("./res_unique.csv")

# Import core mutation rate -----------------------------------------------
theta_vals <- c(seq(0.01, 0.04, 0.01), seq(0.1, 0.4, 0.05))

cl <- makeForkCluster(ncores)   # fork shares res_unique_master (copy-on-write)
registerDoParallel(cl)

foreach(core_mut_rate = theta_vals,
        .packages = "data.table",
        .errorhandling = "pass") %dopar% {
          
          cat("\n\nStarting with theta piS = ", core_mut_rate,"\n")
          
          res_unique <- copy(res_unique_master)
          
          # import m target
          
          m_target = 435
          
          
          # Scale trees to gene length (m) ------------------------------------------
          # only use 500K trees
          
          
          # calculate the expected number of segregating sites ----------------------
          # we multiply the total branch length of each sub-tree 
          # that produces x descendants by half our estimate of 2Nu from the core genes. 
          
          piS_theta = (core_mut_rate)/2
          
          res_unique[, pi_exp_segsites := subtree_br_len * piS_theta]
          
          # put on a scale of segregating sites per gene using median AG length
          
          res_unique[, pi_exp_segsites := pi_exp_segsites * m_target]
          
          # Then using this expectation we generate the probability of actually 
          # observing 0, 1, 2….etc segregating sites using the Poisson distribution.
          
          # Then we sum these probabilities across sub-trees to generate the final 
          # distribution of the number of segregating sites. We then compare our 
          # observations against this distribution.
          
          # Determine maximum number of segregating sites to consider
          max_seg_sites <- qpois(0.99999, max(res_unique$pi_exp_segsites))
          
          site_rng <- 0:max_seg_sites
          
          s_cols <- paste0("s_", site_rng)
          
          
          # Using piS-based expectation ----------------------------------------------
          prob_pi <- matrix(0, nrow(res_unique), length(site_rng))
          
          for (j in seq_along(site_rng)) {
            prob_pi[, j] <- dpois(site_rng[j], res_unique$pi_exp_segsites)
          }
          
          colnames(prob_pi) = paste0("s_", site_rng)
          
          prob_pi <- as.data.frame(prob_pi)
          
          setDT(prob_pi)
          
          prob_pi <- cbind(res_unique[, .(freq,tree_id,pi_exp_segsites)], prob_pi)
          
          prob_pi <- prob_pi[, lapply(.SD, sum, na.rm = TRUE), .SDcols = s_cols, by = freq]
          
          # -------------------------------------------------------------------------
          
          setorderv(prob_pi, cols = "freq", order = 1L)
          
          # Row sums for normalisation
          prob_pi[, row_total := rowSums(.SD), .SDcols = s_cols]
          
          # Normalise each column by row total
          prob_pi[, (s_cols) := lapply(.SD, function(x) x / row_total), 
                  .SDcols = s_cols]
          
          prob_pi[, row_total:= NULL]
          
          # save normalised prob_theta
          fwrite(prob_pi, paste0("./theta_pis_out/prob_theta_pis_norm_",core_mut_rate,".csv"))
          
          prob_pi_lng <- melt(prob_pi,
                              id.vars = "freq",
                              variable.name = "seg_sites", 
                              value.name = "p")
          
          # convert seg_sites to numeric
          prob_pi_lng[, seg_sites := as.integer(gsub("s_", "", seg_sites))]
          
          
          # Calculate cdf -----------------------------------------------------------
          # compute the cumulative distribution and 95% CI boundaries
          
          # Compute CDF per frequency class
          setorder(prob_pi_lng, freq, seg_sites)
          
          prob_pi_lng[, cdf := cumsum(p), by = freq]
          
          fwrite(prob_pi_lng, paste0("./theta_pis_out/prob_pi_lng_",core_mut_rate,".csv"))
          
          # Extract 2.5% and 97.5% quantile boundaries per k
          pis_ci_bounds <- prob_pi_lng[, .(
            lower_95 = seg_sites[which.max(cdf >= 0.025)],  # 2.5th percentile
            median_S = seg_sites[which.max(cdf >= 0.500)],  # 50th percentile
            upper_95 = seg_sites[which.max(cdf >= 0.975)],   # 97.5th percentile
            iqr      = seg_sites[which.max(cdf >= 0.750)] - seg_sites[which.max(cdf >= 0.250)]
          ), by = freq]
          
          pis_ci_bounds$s_set = "piS_theta"
          
          fwrite(pis_ci_bounds, paste0("./theta_pis_out/ci_summary_",core_mut_rate,".csv"))
          
          cat("\n\nDone with theta piS  =", core_mut_rate)
          
        }

stopCluster(cl)