
# 1. For a tree, compute freq & subtree length for each internal node

analyse_pi_desc_subtrees <- function(tree, n_genomes, internal_nodes){

  children <- split(tree$edge[,2], tree$edge[,1])
  
  # preallocate
  freq <- numeric(length(internal_nodes))
  subtree_br_len <- numeric(length(internal_nodes))
  
  for(i in seq_along(internal_nodes)){
    node <- internal_nodes[i]
    # recursively collect all descendant tips
    desc <- node
    tips <- integer(0)
    stack <- node
    while(length(stack) > 0){
      n <- stack[1]; stack <- stack[-1]
      if(n <= n_genomes) tips <- c(tips, n) else stack <- c(stack, children[[as.character(n)]])
    }
    freq[i] <- length(tips)
    subtree_br_len[i] <- sum(tree$edge.length[tree$edge[,2] %in% tips])
  }
  
  data.table(internal_node = internal_nodes, freq = freq, subtree_br_len = subtree_br_len)
}

# 2. Main simulation: many trees, collect results
simulate_freq_vs_diversity_disk <- function(n_genomes, n_trees, out_file, chunk_size = 1000, seed = 42) {
  set.seed(seed)
  
  # define columns
  dt_empty <- data.table(
    internal_node = integer(),
    freq = integer(),
    subtree_br_len = numeric(),
    tree_id = integer()
  )
  
  # Open file connection for appending results
  fwrite(dt_empty, file = out_file)

  # get internal nodes
  n_node <- n_genomes - 1
  internal_nodes <- (n_genomes + 2):(n_genomes + n_node) # exclude the root node (AGs)
  
  for(t_start in seq(1, n_trees, by = chunk_size)) {
    t_end <- min(t_start + chunk_size - 1, n_trees)
    trees_chunk <- lapply(t_start:t_end, function(t) rcoal(n_genomes))
    
    chunk_res <- rbindlist(lapply(seq_along(trees_chunk), function(i) {
      tree <- trees_chunk[[i]]
      tree_nodes_dt <- analyse_pi_desc_subtrees(tree, n_genomes, internal_nodes)
      tree_nodes_dt[, tree_id := t_start + i - 1]
    }))
    
    # append to disk
    fwrite(chunk_res, file = out_file, append = TRUE)
    rm(chunk_res, trees_chunk)
    gc()
    
    # Print progress
    cat(sprintf("Trees %d - %d written to disk\n", t_start, t_end))
    
  }
  
  # read final file if needed
  fread(out_file)
}


# Global variables --------------------------------------------------------
tot_pangenome_size = 260
n_sim_trees = 1000000
file_out_name = paste0(outdir_dat, "/neutral_sim_260_ntree_1000000.csv")

# Run simulations (adjust n_tips, n_trees to taste; beware runtime if n_trees large)

# start timer
start.time <- Sys.time()

res <- simulate_freq_vs_diversity_disk(n_genomes = tot_pangenome_size, 
                                       n_trees = n_sim_trees,
                                       out_file = file_out_name, 
                                       chunk_size = 10000)
# get duration
end.time <- Sys.time()
end.time - start.time
# 260 genomes, 10,000; Time difference of 2.04299 mins
# 260 genomes, 1,000,000; Time difference of 3.405551 hours

# units are in 4neu -> need to halve for haploids

# Normalise/scaled results by mean freq subtree_br_len
res[, t_hat := subtree_br_len/mean(subtree_br_len), freq]

# export tree
fwrite(res, file_out_name)



