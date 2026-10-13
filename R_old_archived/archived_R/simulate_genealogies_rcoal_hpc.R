library(ape)
library(data.table)

# Subtree tip count (k) and summed subtree branch length, per internal node.
# Only want AGs, so the root is dropped, so k is always < n_genomes.
# Set include_tips = TRUE if you also want k = 1.

analyse_pi_desc_subtrees <- function(tree, n_genomes, include_tips = FALSE) {

  tr     <- reorder(tree, "postorder")   # children always precede their parent
  parent <- tr$edge[, 1]
  child  <- tr$edge[, 2]
  el     <- tr$edge.length   

  tot <- n_genomes + tr$Nnode
  k   <- integer(tot)    # descendant tip count
  bl  <- numeric(tot)    # summed branch length below the node
  k[seq_len(n_genomes)] <- 1L

  for (i in seq_along(child)) {
    p <- parent[i]; ch <- child[i]
    k[p]  <- k[p]  + k[ch]
    bl[p] <- bl[p] + bl[ch] + el[i]     # child's subtree + the edge to it
  }

  root  <- parent[length(parent)]        # last post-order edge belongs to the root
  nodes <- if (include_tips) seq_len(tot) else (n_genomes + 1):tot
  nodes <- setdiff(nodes, root)

  out <- data.table(internal_node = nodes,
                    freq          = k[nodes],
                    subtree_br_len = bl[nodes])
  out[freq < n_genomes]                  # explicit guard: keep k < n only
}

# 2. Main simulation: many trees, streamed to disk in chunks.

simulate_freq_vs_diversity_disk <- function(n_genomes, n_trees, out_file,
                                            chunk_size = 1000, seed = 42) {
  set.seed(seed)

  fwrite(data.table(internal_node = integer(),
                    freq          = integer(),
                    subtree_br_len = numeric(),
                    tree_id       = integer()),
         file = out_file)

  for (t_start in seq(1, n_trees, by = chunk_size)) {
    t_end <- min(t_start + chunk_size - 1, n_trees)

    chunk_res <- rbindlist(lapply(t_start:t_end, function(t) {
      dt <- analyse_pi_desc_subtrees(rcoal(n_genomes), n_genomes)
      dt[, tree_id := t]
    }))

    fwrite(chunk_res, file = out_file, append = TRUE)
    rm(chunk_res); gc()
    cat(sprintf("Trees %d - %d written to disk\n", t_start, t_end))
  }
  invisible(out_file)
}


# Globals -----------------------------------------------------------------
tot_pangenome_size <- 260
n_sim_trees        <- 500000
file_out_name      <- file.path(outdir_dat, "neutral_sim_260_ntree_500K.csv")

start.time <- Sys.time()
simulate_freq_vs_diversity_disk(n_genomes  = tot_pangenome_size,
                                n_trees    = n_sim_trees,
                                out_file   = file_out_name,
                                chunk_size = 10000)
Sys.time() - start.time

# Normalise by the mean subtree length within each k class.
res <- fread(file_out_name)           
res[, t_hat := subtree_br_len / mean(subtree_br_len), by = freq]
fwrite(res, file_out_name)
