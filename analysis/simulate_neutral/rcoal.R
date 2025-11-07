
# 1. For a tree, compute freq & subtree length for each internal node
analyse_length_desc_subtrees <- function(tree) {
  
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
  
  return(rbindlist(out))
  
}


# 2. Main simulation: many trees, collect results
simulate_freq_vs_diversity <- function(n_genomes, n_trees, seed = 42, num_core = 2) {
  set.seed(seed)
  
  # Set up parallel cluster
  cl <- parallel::makeCluster(num_core)
  doParallel::registerDoParallel(cl)
  
  # Export custom functions to workers
  parallel::clusterExport(cl, varlist = c("analyse_length_desc_subtrees"), envir = environment())
  
  # Run in parallel
  results <- foreach(t = seq_len(n_trees), .combine = rbindlist,
                     .packages = c("ape", "phytools", "data.table")) %dopar% {
                       tree <- rcoal(n_genomes)
                       tree_nodes_dt <- analyse_length_desc_subtrees(tree)
                       tree_nodes_dt$tree_id <- t
                       # wrap inside a list so .combine=rbindlist can concatenate properly
                       list(tree_nodes_dt[, c("tree_id", "node", "freq", "subtree_br_len")])
                     }
  
  # Clean up
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  
  return(results)
}

# Example run (adjust n_tips, n_trees to taste; beware runtime if n_trees large)
res <- simulate_freq_vs_diversity(n_genomes = 93, n_trees = 500)

# Base-R plot: frequency vs subtree_br_len (scatter + LOWESS smooth)
# jitter freq slightly for visualization
plot(res$subtree_br_len, res$freq, #jitter(res$freq, amount = 0.2),
     ylab = "Accessory gene frequency",
     bty = "l", yaxt = "n",
     col = rgb(0.25,0.25,0.8, alpha = 0.5),
     xlab = "Subtree total branch length",
     main = "Accessory gene frequency vs subtree length",
     pch = 19, cex = 0.6)

axis(side = 2, at = pretty(c(0, max(res$freq))), 
     las = 2, labels = pretty(c(0, max(res$freq))))
