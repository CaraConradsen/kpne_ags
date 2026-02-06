
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
tot_pangenome_size = 214
n_sim_trees = 10000

# Example run (adjust n_tips, n_trees to taste; beware runtime if n_trees large)
res <- simulate_freq_vs_diversity(n_genomes = tot_pangenome_size, n_trees = n_sim_trees)

# Normalise/scaled results by mean freq subtree_br_len
res[, t_hat := subtree_br_len/mean(subtree_br_len), freq]

fwrite(res, paste0(outdir_dat,"/neutral_sim_214_ntree_10000.csv"))

# Scale trees -------------------------------------------------------------
# res <- fread(paste0(outdir_dat,"/neutral_sim_214_ntree_10000.csv"))
# THERE ARE MISSING GENE FAMILY ALIGNMENTS???

AGs_post_cat <- fread(paste0(outdir_dat, "/AGs_post_cat.csv"))

AGs_post_cat <- AGs_post_cat[multi_gain!=1]

# Get non-paralog genes
focal_gene_families = unique(AGs_post_cat$gene_family)

# focal_gene_families =  c(focal_gene_families[!grepl("_", focal_gene_families)],
#                          focal_gene_families[grepl("_1", focal_gene_families)])

gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_485_lng_rds_out/feature_sequences/"

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallel foreach loop
ag_seg_sites <- foreach(i = focal_gene_families, 
                        .combine = rbind, 
                        .packages = c("Biostrings","ape","pegas","data.table")) %dopar% {
  
  # Read the FASTA, skip if missing
  temp_string <- tryCatch(
    readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
    error = function(e) return(NULL)
  )
  
  # skip if file is missing
  if (is.null(temp_string)) return(NULL)
  
  #  Subset to loci present in pangraph annotation
  sampled_loci = AGs_post_cat[gene_family == i, locus_tag]
  
  # Subset by sequence names
  temp_string <- temp_string[names(temp_string) %in% sampled_loci]
  
  # Skip if no sequence, ORFan or core gene
  if (length(temp_string) == 0) return(NULL)
  if (length(temp_string) %in% c(1, tot_pangenome_size)) return(NULL)
  
  # Compute gene length and frequency 
  gene_len = width(temp_string)[1]
  freq = length(temp_string)
  
  # Convert to DNAbin
  dna_bin <- as.DNAbin(temp_string)
  
  # Calculate segregating sites 
  seg_site <- seg.sites(dna_bin)
  S <- length(seg_site)
  
  data.frame(
    gene_family = i,
    m = gene_len,
    freq = freq,
    S = S
  )
}

# Stop cluster when done
stopCluster(cl)

setDT(ag_seg_sites)

# Calculate segregating site per site
ag_seg_sites[, Sm := S / m]

# Calculate empirical mean for frequency class k
pi_emp_k = ag_seg_sites[, .(mean_Sm = mean(Sm)), by = freq]



# Anchor genealogy estimates by empirical means ---------------------------

# Join pi_emp_k onto res by freq and create new column
res[pi_emp_k, on = "freq", pi_sim := t_hat * mean_Sm]


# # Simulated frequency vs pi_sim
# plot(res$freq/tot_pangenome_size,res$pi_sim,
#      xlab = "Accessory gene frequency",
#      bty = "l", yaxt = "n",
#      col = rgb(0.25,0.25,0.8, alpha = 0.5),
#      ylab = "Simulated pi",
#      main = paste0("Neutral distribution,\nn = ",n_sim_trees,
#                    ", pangenome = ",tot_pangenome_size),
#      pch = 19, cex = 0.6)
# 
# axis(side = 2, at = pretty(c(0, max(res$pi_sim, na.rm = TRUE))),
#      las = 2, labels = pretty(c(0, max(res$pi_sim, na.rm = TRUE))))
# 

# Outlier detections: summarise the neutral distribution per frequency --------

neutral_summary <- res[, .(
  pi_lower = quantile(pi_sim, 0.025, na.rm = TRUE),
  pi_median = median(pi_sim, na.rm = TRUE),
  pi_upper = quantile(pi_sim, 0.975, na.rm = TRUE),
  iqr = IQR(pi_sim, na.rm = TRUE, type = 7)
), by = freq]

fwrite(neutral_summary, paste0(outdir_dat, "/neutral_summary.csv"))
