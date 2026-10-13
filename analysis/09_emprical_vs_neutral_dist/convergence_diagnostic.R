# Using Welford's (1962)

# Stabilising the variance ------------------------------------------------
file_out_name = paste0(outdir_dat, "/neutral_sim_260_ntree_1000000.csv")
# res <- fread(file_out_name)


# For each frequency class k, compute variance of subtree_br_len
n_genomes <- 260
n_trees   <- 1000000
step_size <- 10000
k_reps    <- round(seq(2, n_genomes - 1, length.out = 10), digits = 0)
breaks    <- c(0, seq(step_size, n_trees, by = step_size))

# --- Key res by freq then tree_id for fast subsetting ---
setkey(res, freq, tree_id)

# --- Initialise one Welford accumulator per k ---
accum <- data.table(
  freq = k_reps,
  n    = 0L,
  mean = 0,
  M2   = 0
)
setkey(accum, freq)

# --- Pre-filter res to only k_reps rows once, outside the loop ---
res_k <- res[freq %in% k_reps]
setkey(res_k, tree_id, freq)

# --- Single sequential pass over batches ---
var_convergence <- vector("list", length(breaks) - 1)

for (b in seq_len(length(breaks) - 1)) {
  
  lo <- breaks[b] + 1
  hi <- breaks[b + 1]
  
  # Extract this batch - fast keyed subset, no copying of full res
  batch <- res_k[tree_id %between% c(lo, hi)]
  
  # Compute batch-level Welford statistics per k in one vectorised step
  batch_stats <- batch[, .(
    n_new    = .N,
    mean_new = mean(subtree_br_len),
    M2_new   = fifelse(.N > 1, var(subtree_br_len) * (.N - 1L), 0)
  ), by = freq]
  
  # Update accumulators for all k simultaneously using parallel Welford
  accum[batch_stats, on = "freq", `:=`(
    n    = {
      nA <- n; nB <- i.n_new
      nA + nB
    },
    mean = {
      nA <- n; nB <- i.n_new
      (nA * mean + nB * i.mean_new) / (nA + nB)
    },
    M2   = {
      nA <- n; nB <- i.n_new; delta <- i.mean_new - mean
      M2 + i.M2_new + delta^2 * nA * nB / (nA + nB)
    }
  )]
  
  # Snapshot variance for all k at this checkpoint
  var_convergence[[b]] <- accum[n > 1L, .(
    freq         = freq,
    var_t        = M2 / (n - 1L),
    n_obs        = n,
    n_trees_done = hi
  )]
  
  if (b %% 10 == 0) cat("Batch", b, "of", length(breaks) - 1, "done\n")
}

var_convergence <- rbindlist(var_convergence)

# Plot --------------------------------------------------------------------

png(paste0(outdir_fig,"/variance_of_subtree_variance.png"),
width = 15.9, height = 9.25, units = "cm", res = 300,
pointsize = 12, type = "cairo")

par(mar = c(6,5,1,1))

plot(NULL,
     xlim = c(0, max(var_convergence$n_trees_done)),
     ylim = c(0, max(var_convergence$var_t)),
     bty = "L",
     yaxt = "n",
     xaxt = "n",
     xlab  = "Number of simulated trees",
     ylab  = "")

labels <- c("0", "200K", "400K", "600K", "800K", "1M")
axis(1, at = seq(0, 1e6, 2e5), labels = labels)

axis(side = 2,
     at = seq(0, 0.08, 0.02),
     labels = sprintf("%1.2f", seq(0, 0.08, 0.02)),
     las = 2)
mtext("Variance of subtree branch length",
      side = 2, line = 3.5)

colrs = rainbow(length(k_reps))

for (i in seq_along(k_reps)) {
  k = k_reps[i]
  
  with(var_convergence[freq == k],
       lines(n_trees_done, var_t,
             col = colrs[i],
             lwd = 2))
  
}

legend("bottom", 
       legend = k_reps,
       col = colrs,
       lwd = 2,
       ncol = length(k_reps), 
       xpd = TRUE,    
       inset = c(0, -0.5),     # Moves it below the x-axis
       bty = "n",
       cex = 0.5)

dev.off()
