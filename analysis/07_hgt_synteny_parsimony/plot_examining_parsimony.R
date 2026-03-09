# Will need:
#   1. the recombination free core genome phylogeny
#   2. a binary presence / absence file to recode

# Check convert gubbins final tre file to newick is rooted -----------------------------
core_gub_tree <- read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre") 

is.rooted(core_gub_tree)

# Create binary ag presence_absence data.table------------------------------------------------------------
# needs a dummy start and end, where rows are each gene family and columns are genomes

ag_info <- fread(paste0(outdir_dat,"/pangraph_anno.csv"),
                 select = c("geno_id","gene_family", "ag_type"))

# id orfans
ag_info[, n := .N, gene_family]

# remove paralogs (retain genes with only _1 or no underscore)

keep_ags = unique(c(grep("_", ag_info$gene_family, invert = TRUE, value = TRUE),
                    grep("_1", ag_info$gene_family, value = TRUE)))

ag_info <- ag_info[gene_family %chin% keep_ags]


# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_info[n!=1 & ag_type !="core",.(geno_id, gene_family)],
  formula = gene_family ~ geno_id,
  fun.aggregate = length,
  value.var = "gene_family"
)

# make dummy regions
dummy_regions <- data.frame(
  start = seq(1, by = 2, length.out = nrow(ag_presence_absence)),
  end   = seq(2, by = 2, length.out = nrow(ag_presence_absence))
)

ag_presence_absence <- cbind(dummy_regions, ag_presence_absence)

setDT(ag_presence_absence)


# Maximum parsimony -------------------------------------------------------
htree = core_gub_tree

gene_tip <- t(ag_presence_absence[gene_family=="g000076"][1,-(1:3)])

gene_mat <- factor(gene_tip)

names(gene_mat) <- row.names(gene_tip)

# Convert binary presence/absence to phyDat
gene_phy <- phyDat(gene_mat, type="USER", levels=c(0,1))

# Perform parsimony ancestral reconstruction
anc <- ancestral.pars(core_gub_tree, gene_phy, type="MPR")  # MPR = most parsimonious reconstruction

# Extract most parsimonious state for each node/tip
states <- sapply(anc, function(x) which.min(x) - 2)*-1


# 1. Create a vector of node states in order of node numbers
# ape numbers tips 1:Ntip and internal nodes Ntip+1:Ntip+Nnode
all_nodes <- c(htree$tip.label, htree$node.label)
node_states <- states[all_nodes]

# 2. Initialize edge colors
edge_colors <- rep("grey40", nrow(htree$edge))

# 3. Loop through edges and assign colors based on outer and inner node states
for(i in 1:nrow(htree$edge)){
  parent_node <- htree$edge[i,1]
  child_node  <- htree$edge[i,2]
  
  parent_state <- node_states[parent_node]
  child_state  <- node_states[child_node]
  
  # Apply coloring rules
  if(parent_state == 1 & child_state == 1){
    edge_colors[i] <- "dodgerblue1"
  } else {
    edge_colors[i] <- "grey40"
  }
}

# 4. Plot the tree
plot(htree, edge.color = edge_colors, 
     cex = 0.5,
     align.tip.label = TRUE,
     edge.width = 3,
     show.node.label = FALSE,
     show.tip.label = TRUE)
tiplabels(pch=21,
          bg=ifelse(node_states[1:Ntip(htree)]==1,"dodgerblue1","grey40"),
          cex=0.5)


# 1. Convert tip states into a matrix (phangorn needs a matrix)
char_matrix <- as.matrix(states[htree$tip.label])
rownames(char_matrix) <- htree$tip.label

# 2. Use parsimony score function
pscore <- parsimony(htree, gene_phy, method = "fitch")

# 3. Calculate minimum number of steps (for binary character)
min_steps <- max(sum(gene_phy == 0), sum(gene_phy == 1))

# 4. Consistency index
CI <- min_steps / pscore
CI




