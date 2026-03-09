# tree: rooted phylo object
# gene_matrix: matrix of 0/1, rows = taxa, columns = genes
# taxa names must match tree$tip.label

core_tree <- read.tree("./input_data/bootstrapped_gubbins/tmp8yl_0c9w/RAxML_bestTree.core_genome_aln.iteration_20")

core_tree_resolved <- multi2di(core_tree)

# Get remaining accessory genes after synteny analysis --------------------
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

ag_homo_dat <- msu_regions_anchored[acrs_msu!=1][acrs_jun != 1][ag_type != "core"][number_genomes != 1][!grepl("_",gene_family)]

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
  .packages = c("ape")
) %dopar% {
  
  trait <- factor(gene_matrix[, g], levels = c(0,1))
  
  trait <- trait[core_tree_resolved$tip.label]
  
  # Use a matrix instead of data.frame to avoid the 1-sequence problem
  phy_data <- phyDat(as.matrix(trait), type = "USER", levels = c(0,1))
  
  
  fit <- try(ancestral.pars(
    core_tree_resolved, 
    phy_data, type = "MPR"),
             silent = TRUE)# MPR = most parsimonious reconstruction

  
  # Catch failed fits
  if(inherits(fit, "try-error")) {
    return(data.frame(
      gene_family = g,
      gains = NA,
      losses = NA,
      anc_state = NA
    ))
  }
  
  ntip  <- Ntip(core_tree_resolved)
  nnode <- Nnode(core_tree_resolved)
  
  # Convert list to matrix: rows = nodes, columns = possible states
  node_states <- sapply(fit, function(x) x[1])  # just pick first state
  
  states <- numeric(ntip + nnode)
  states[1:ntip] <- as.numeric(as.character(trait))
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
    gains = gains,
    losses = losses,
    anc_state = node_states["261"][[1]]
  )
}


stopCluster(cl)

fwrite(gain_loss_results, paste0(outdir_dat, "/gain_loss_results.csv"))

