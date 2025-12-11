
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
simulate_freq_vs_diversity <- function(n_genomes, n_trees, seed = 42, num_core = 2) {
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

# Example run (adjust n_tips, n_trees to taste; beware runtime if n_trees large)
res <- simulate_freq_vs_diversity(n_genomes = 93, n_trees = 50)




# Base-R plot: frequency vs subtree_br_len (scatter + LOWESS smooth)
# jitter freq slightly for visualization
plot(res$freq/93,res$subtree_br_len,  #jitter(res$freq, amount = 0.2),
     xlab = "Accessory gene frequency",
     bty = "l", yaxt = "n",
     col = rgb(0.25,0.25,0.8, alpha = 0.5),
     ylab = "Subtree total branch length (Age 4N generations)",
     main = "Accessory gene frequency vs subtree length",
     pch = 19, cex = 0.6)

axis(side = 2, at = pretty(c(0, max(res$subtree_br_len))),
     las = 2, labels = pretty(c(0, max(res$subtree_br_len))))

# axis(side = 2, at = pretty(c(0, max(res$freq))),
#      las = 2, labels = pretty(c(0, max(res$freq))))


# Scale trees -------------------------------------------------------------
# THERE ARE MISSING GENE FAMILY ALIGNMENTS???

pangraph_anno <- fread(paste0(outdir_dat,"/pangraph_anno.csv"))

# Get non-paralog genes
focal_gene_families = unique(pangraph_anno$gene_family)

focal_gene_families =  c(focal_gene_families[!grepl("_", focal_gene_families)],
                         focal_gene_families[grepl("_1", focal_gene_families)])

gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_485_lng_rds_out/feature_sequences/"

ag_seg_sites = lapply(focal_gene_families, function(i) {
  
  # --- Try reading the FASTA (skip if missing) ---
  temp_string <- tryCatch(
    readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
    error = function(e) return(NULL)
  )
  
  # If file missing → skip
  if (is.null(temp_string)) return(NULL)
  
  # --- Subset to loci present in pangraph annotation ---
  sampled_loci = pangraph_anno[gene_family == i, locus_tag]
  
  # Subset by sequence names; avoid errors if some names not found
  temp_string <- temp_string[names(temp_string) %in% sampled_loci]
  
  # Skip if 
  if (length(temp_string) == 0) return(NULL)
  if (length(temp_string) %in% c(1, 93)) return(NULL)
  
  # --- Compute gene length and frequency ---
  gene_len = width(temp_string)[1]
  freq = length(temp_string)
  
  # --- Convert to DNAbin ---
  dna_bin <- as.DNAbin(temp_string)
  
  # --- Compute segregating sites ---
  segpos <- seg.sites(dna_bin)
  S <- length(segpos)
  
  data.frame(
    gene_family = i,
    gene_len = gene_len,
    freq = freq,
    S = S
  )
})

ag_seg_sites = rbindlist(ag_seg_sites)

ag_seg_sites[, .(mean_S = mean(S)), by = freq]




