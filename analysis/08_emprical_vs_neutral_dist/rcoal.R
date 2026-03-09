
# 1. For a tree, compute freq & subtree length for each internal node
analyse_pi_desc_subtrees  <- function(tree, n_genomes) {
  
  # tree metrics
  Ntip <- length(tree$tip.label)
  Nnode <- tree$Nnode
  internal_nodes <- (Ntip + 1):(Ntip + Nnode)
  
  out <- lapply(internal_nodes, function(node){
    
    # get descendants nodes
    desc_tips <- phytools::getDescendants(tree, node)
    # return tips
    desc_tips <- desc_tips[desc_tips <= Ntip]
    
    freq <- length(desc_tips)
    
    # subtree: keep only those tips and sum branch lengths
    # use keep.tip; it preserves edge.length for the pruned tree
    subtree <- keep.tip(tree, tree$tip.label[desc_tips])
    subtree_br_len <- sum(subtree$edge.length)
    data.frame(node = node, freq = freq, subtree_br_len = subtree_br_len)
    
  })
  
  # remove the root node
  return(rbindlist(out)[freq != n_genomes])
  
}




# 2. Main simulation: many trees, collect results
simulate_freq_vs_diversity <- function(n_genomes, n_trees, seed = 42, num_core = 10) {
  set.seed(seed)
  
  # Set up parallel cluster
  cl <- parallel::makeCluster(num_core)
  doParallel::registerDoParallel(cl)
  
  # Export custom functions to workers
  parallel::clusterExport(cl, varlist = c("analyse_pi_desc_subtrees"), envir = environment())
  
  # Run in parallel
  results <- foreach(t = seq_len(n_trees), .combine = "rbind",
                     .packages = c("ape", "phytools", "data.table")) %dopar% {
                       tree <- rcoal(n_genomes)
                       tree_nodes_dt <- analyse_pi_desc_subtrees (tree, n_genomes)
                       tree_nodes_dt$tree_id <- t
                       # wrap inside a list so .combine=rbindlist can concatenate properly
                       tree_nodes_dt[, c("tree_id", "node", "freq", "subtree_br_len")]
                     }
  
  # Clean up
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  return(results)
}




# Global variables --------------------------------------------------------
tot_pangenome_size = 260
n_sim_trees = 10000

# Run simulations (adjust n_tips, n_trees to taste; beware runtime if n_trees large)

# start timer
start.time <- Sys.time()

res <- simulate_freq_vs_diversity(n_genomes = tot_pangenome_size, n_trees = n_sim_trees)

# get duration
end.time <- Sys.time()
end.time - start.time# 260 genomes, 10,000; Time difference of 5.397505 mins

# Normalise/scaled results by mean freq subtree_br_len
res[, t_hat := subtree_br_len/mean(subtree_br_len), freq]

fwrite(res, paste0(outdir_dat,"/neutral_sim_260_ntree_10000.csv"))

# Scale trees -------------------------------------------------------------
# res <- fread(paste0(outdir_dat,"/neutral_sim_260_ntree_10000.csv"))

# from get_ag_pis.R
ag_age_S_dt <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

# THERE ARE MISSING GENE FAMILY ALIGNMENTS??? >> increased dosage (n = 21), not included in features_alignment

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# Calculate empirical mean for frequency class k
pi_emp_k_set_1 = ag_age_S_dt[gene_family %chin% set_1_ags$gene_family, .(mean_Sm = mean(Sm)), by = freq]
pi_emp_k_set_2 = ag_age_S_dt[gene_family %chin% set_2_ags$gene_family, .(mean_Sm = mean(Sm)), by = freq]

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
