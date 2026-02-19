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

  fit <- try(ace(trait, core_tree_resolved,
                 type="discrete", model="ARD"),
             silent = TRUE)

  # Catch failed fits
  if(inherits(fit, "try-error") || any(is.na(fit$lik.anc))) {
    return(data.frame(
      gene_family = g,
      gains = NA,
      losses = NA,
      anc_state = NA
    ))
  }

  ntip  <- Ntip(core_tree_resolved)
  nnode <- Nnode(core_tree_resolved)

  node_states <- apply(fit$lik.anc, 1, function(x) {
    as.numeric(names(x)[which.max(x)])
  })

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



# Bootstrap estimates -----------------------------------------------------

# get 100 bootstrapped trees
boot_trees <- read.tree("./input_data/bootstrapped_gubbins/tmp8yl_0c9w/RAxML_bootstrap.core_genome_aln.iteration_20.bootstrapped_trees")

boot_trees <- boot_trees[1:100]

# tidy trees
boot_trees <- lapply(boot_trees, function(tr) {
  # assign length 1 to all branches
  tr$edge.length <- rep(1, nrow(tr$edge))
  
  # remove polytomies
  tr <- multi2di(tr)
  
  tr <- midpoint.root(tr)
  
})

class(boot_trees) <- "multiPhylo"



cl <- makeCluster(num_cores - 2)
registerDoParallel(cl)

genes <- colnames(gene_matrix)
n_genes <- length(genes)
ntip  <- Ntip(core_tree_resolved)
nnode <- Nnode(core_tree_resolved)

gain_loss_boot_results <- foreach(
  g = genes,
  .combine = rbind,
  .packages = c("ape", "data.table")
) %dopar% {
  
  trait <- factor(gene_matrix[, g], levels = c(0,1))
  
  trait <- trait[core_tree_resolved$tip.label]
  
  fits <- lapply(boot_trees, function(tree){
                try(ace(trait, tree,
                 type="discrete", model="ARD"),
             silent = TRUE)})
  
  # Catch failed fits
  check_fit <- lapply(fits, function(fit){
    if(inherits(fit, "try-error") || any(is.na(fit$lik.anc))) {
      c(FALSE)
    }else{c(TRUE)}
  })
  
  # Convert list of length-1 logicals to a logical vector
  keep_trees <- sapply(check_fit, `[[`, 1)
  
  fits <- fits[keep_trees]
  
  # Catch failed fits
  if(length(fits)==0) {
    return(data.frame(
      gains = NA,
      losses = NA,
      gene_family = g,
      rep = NA
    ))
  }
  
  
  all_states <- lapply(fits, function(fit){
    node_states <- apply(fit$lik.anc, 1, function(x) {
      as.numeric(names(x)[which.max(x)])
    })
    states <- numeric(ntip + nnode)
    states[1:ntip] <- as.numeric(as.character(trait))
    states[(ntip + 1):(ntip + nnode)] <- node_states
    
    states
  })

  all_gains_losses <- lapply(seq_along(boot_trees[keep_trees]), function(i){
    
    edges <- boot_trees[[i]]$edge
    
    state <- all_states[[i]]
    
    gains  <- sum(state[edges[,1]] == 0 &
                    state[edges[,2]] == 1,
                  na.rm = TRUE)
    
    losses <- sum(state[edges[,1]] == 1 &
                    state[edges[,2]] == 0,
                  na.rm = TRUE)
    
    data.frame(gains, losses)
    
  })

  all_gains_losses <- rbindlist(all_gains_losses)
  
  all_gains_losses$gene_family = g
  
  all_gains_losses$rep = 1:nrow(all_gains_losses)

  all_gains_losses
}


stopCluster(cl)

fwrite(gain_loss_boot_results, paste0(outdir_dat, "/gain_loss_boot_results.csv"))

