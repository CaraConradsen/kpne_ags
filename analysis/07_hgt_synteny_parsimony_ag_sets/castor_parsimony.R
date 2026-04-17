# tree: rooted phylo object
# gene_matrix: matrix of 0/1, rows = taxa, columns = genes
# taxa names must match tree$tip.label

core_tree <- read.tree("./input_data/bootstrapped_pirate_gubbins/boot_gub_pir_100.final_bootstrapped_tree.tre")

core_tree_resolved <- multi2di(core_tree)

# Get remaining accessory genes after synteny analysis --------------------
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

ag_homo_dat <- msu_regions_anchored[ag_type != "core"][number_genomes != 1]#[acrs_msu!=1][acrs_jun != 1]

# create presence/absence for gene_families -------------------------------
# Create binary ag presence_absence data.table 
# needs a dummy start and end, where rows are each gene family and columns are genomes

# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_homo_dat[,.(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

gene_matrix <- as.matrix(ag_presence_absence[,-1])

# for castor need 1 -> n states
gene_matrix <- gene_matrix + 1

row.names(gene_matrix) <- ag_presence_absence$geno_id



# Run analysis ------------------------------------------------------------

cl <- makeCluster(num_cores - 2)
registerDoParallel(cl)

genes <- colnames(gene_matrix)
n_genes <- length(genes)
ntip  <- Ntip(core_tree_resolved)
nnode <- Nnode(core_tree_resolved)

gain_loss_results <- foreach(
  g = genes,
  .combine = rbind,
  .packages = c("castor", "ape")
) %dopar% {
  
  trait <- gene_matrix[, g]
  
  trait <- trait[core_tree_resolved$tip.label]
  

  fit <- asr_max_parsimony(core_tree_resolved, 
                               trait, Nstates=2, 
                               transition_costs="all_equal")

  ntip  <- Ntip(core_tree_resolved)
  nnode <- Nnode(core_tree_resolved)
  
  internal_node_numbers <- (ntip + 1):(ntip + nnode)
  
  # Convert list to matrix: rows = nodes, columns = possible states
  node_states <- max.col(fit$ancestral_likelihoods) - 1  # get binary states
  
  node_states_named <- setNames(node_states, internal_node_numbers)
  
  states <- numeric(ntip + nnode)
  states[1:ntip] <- as.numeric(as.character(trait - 1))
  states[(ntip + 1):(ntip + nnode)] <- node_states
  
  edges <- core_tree_resolved$edge
  
  gains  <- sum(states[edges[,1]] == 0 &
                  states[edges[,2]] == 1,
                na.rm = TRUE)
  
  losses <- sum(states[edges[,1]] == 1 &
                  states[edges[,2]] == 0,
                na.rm = TRUE)
  
  data.frame(
    gene_family = g,
    parsimony_score = fit$total_cost,
    gains = gains,
    losses = losses,
    anc_state = node_states_named["261"][[1]]
  )
}


stopCluster(cl)

fwrite(gain_loss_results, paste0(outdir_dat, "/maximum_parsimony_gain_loss_results.csv"))

