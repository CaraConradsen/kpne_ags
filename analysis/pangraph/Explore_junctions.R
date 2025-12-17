

msu_regions_anchored <- fread(paste0(outdir_tab, "/msu_regions_anchored.csv"))



# Examine the msu blocks ---------------------------------------------------

msu_loci <- msu_regions_anchored[grepl(":", anchor)]

# Are there any AGs across junctions?
msu_loci_grp = unique(msu_loci[,.(ag_type, gene_family, msu)])

across_msu_loci = msu_loci_grp[,.(n = .N), by = c("ag_type", "gene_family")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_msu = 0
msu_regions_anchored[gene_family %chin% across_msu_loci, acrs_msu := 1]

# Examine the syntenic junctions ---------------------------------------------------

junction_loci <- msu_regions_anchored[grepl(":", anchor)]

# Are there any AGs across junctions?
junction_loci_grp = unique(junction_loci[,.(ag_type, gene_family, msu, anchor)])

across_junction_loci = junction_loci_grp[,.(n = .N), by = c("ag_type", "gene_family", "msu")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_jun = 0
msu_regions_anchored[gene_family %chin% across_junction_loci, acrs_jun := 1]




# Gene Parsimony ---------------------------------------------------------------

# Generate a binary matrix giving gene presence/absence
PA <- dcast(
  data = msu_regions_anchored[,.(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

PA_mat <- as.matrix(PA[, -1])
rownames(PA_mat) <- PA$geno_id

# get tree
par_tree <- ape::read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre")


n_genes <- ncol(PA_mat)
res_ml <- data.frame(gene = colnames(PA_mat), n_ag = colSums(PA_mat), gains_ml = NA_integer_, stringsAsFactors = FALSE)

for(j in seq_len(n_genes)) {
  x <- PA_mat[, j]
  # skip trivial loci
  if(length(unique(x)) == 1) {
    res_ml$gains_ml[j] <- 0
    next
  }
  # fit discrete model (ARD allows different gain/loss rates)
  fit <- tryCatch(ace(x, par_tree, type = "discrete", model = "ARD",
                      ip = c(1,1), CI = TRUE), error = function(e) NULL)
  if(is.null(fit)) {
    res_ml$gains_ml[j] <- NA
    next
  }
  # node likelihoods: fit$lik.anc is matrix nodes x states
  probs <- fit$lik.anc   # rows = nodes, cols = states (state order follows levels of x -> 0 then 1)
  # get ML state per node
  ml_state_node <- apply(probs, 1, which.max) - 1  # convert to 0/1
  # build full vector of states for parent/child comparison: tips + nodes
  tip_states <- x
  node_states <- ml_state_node
  # name mapping: ape node numbering: tips = 1:Ntip, nodes = Ntip+1 : Ntip+Nnode
  Ntip <- length(par_tree$tip.label)
  full_states <- numeric(Ntip + par_tree$Nnode)
  full_states[1:Ntip] <- tip_states
  full_states[(Ntip+1):(Ntip+par_tree$Nnode)] <- node_states
  
  # count 0->1 transitions on edges
  edge_par <- par_tree$edge[,1]
  edge_child <- par_tree$edge[,2]
  parent_state <- full_states[edge_par]
  child_state  <- full_states[edge_child]
  gains <- sum(parent_state == 0 & child_state == 1, na.rm = TRUE)
  res_ml$gains_ml[j] <- gains
}

# Rank by inferred gains
res_ml <- res_ml[order(-res_ml$gains_ml, -res_ml$n_ag), ]
head(res_ml, 20)

# ace() or even fitMk() will give you the number of gains implicitly (via reconstructed states), 
# but they do not directly annotate which branches those gains occurred on. To locate where on 
# the tree a gene likely appeared (and thus identify candidate HGT events), you need stochastic 
# mapping. This is the standard approach for discrete traits like gene presence/absence.

sim_maps <- make.simmap(par_tree, x, model = "ARD", nsim = 100, pi = "estimated")


# Illustrate jumping gene -------------------------------------------------
#g000077
plot_77 = unique(msu_regions_anchored[gene_family == "g000077", .(anchor, msu)])

plot_77_dat  = rbind(
  msu_regions_anchored[msu == "MSU_4" & anchor %chin% c("A29:A30", "A29", "A30",
                                                        "A88:A89", "A88", "A89")],
  msu_regions_anchored[msu == "MSU_2" & anchor %chin% c("A101:A102", "A101", "A102")],
  msu_regions_anchored[msu == "MSU_0" & anchor %chin% c("A24:A25", "A24","A25",
                                                        "A5:A6", "A5", "A6",
                                                        "A20:A21", "A20","A21")]
)

plot_77_dat <- plot_77_dat[geno_id %chin% msu_regions_anchored[gene_family == "g000077", geno_id]]

plot_77_dat <- plot_77_dat[!ST %chin% c("ST258", "ST258-1LV")]



# synteny plot data
# Part 1. get tree
core_gub_tree <- read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre")

ST_labels <- unique(pangraph_anno[geno_id %chin% core_gub_tree$tip.label,
                             .(geno_id, ST)])

# re-label with ST
ST_labels <- ST_labels[match(core_gub_tree$tip.label, geno_id)]  # reorder

core_gub_tree$tip.label <- ST_labels$ST

ST_gene_tree = keep.tip(core_gub_tree, unique(plot_77_dat$ST))
ST_gene_tree$node.label <- NULL

# Convert ape phylo → Newick string
nwk <- write.tree(ST_gene_tree)

# Convert Newick string → ade4 phylog
tree_phylog <- ade4::newick2phylog(nwk)

# Part 2. get gene annotations
generate_dna_segs <- function(anno_file, synt_dat, geno_ord){
  # get anno min max
  subset_anno_minmax <- synt_dat[, .(
    start = min(start),
    end   = max(end)
  ), by = geno_id]
  
  # get annotations
  subset_anno <- subset_anno_minmax[
    anno_file,
    on = .(geno_id, start <= start, end >= end),
    nomatch = 0
  ]
  
  # set locus_tag as the name column
  colnames(subset_anno) = gsub("locus_tag", "name", colnames(subset_anno))
  
  # fix strand
  subset_anno[, strand:= fcase(strand == "+", 1,
                               default = -1)]
  
  # add plot values
  subset_anno$fill = "grey50"
  subset_anno$col = "grey50"
  subset_anno$gene_type = "headless_arrows"#"arrows"
  subset_anno$lty = 1
  subset_anno$lwd = 1
  subset_anno$pch = 8
  subset_anno$cex = 1
  
  # highlight core
  subset_anno[ag_type == "soft", 
              `:=`(fill = "white", col = "black", gene_type = "blocks")]
  subset_anno[ag_type == "shell", 
              `:=`(fill = "grey20", col = "black", gene_type = "blocks")]
  subset_anno[ag_type == "cloud", 
              `:=`(fill = "brown1", col = "darkred", gene_type = "blocks")]
  
  # order by genome order
  subset_anno[, ST := factor(ST, levels = geno_ord)]
  
  # order the data.table by that factor order
  setorder(subset_anno, ST)
  
  # order columns
  first_cols <- c("name","start", "end", "strand")
  
  # reassign column order
  setcolorder(subset_anno, c(intersect(first_cols, names(subset_anno)), setdiff(names(subset_anno), first_cols)))
  
  
  # Split by ST, drop ST column in each sub-table
  
  seq_list = split(subset_anno, by = "ST")
  
  # seq_list = lapply(seq_list, dna_seg)
  
  return(seq_list)
  
  
}

dna_seqs <- generate_dna_segs(anno_file = plot_77_dat, geno_ord = ST_gene_tree$tip.label,
                              synt_dat = plot_77_dat)

pdf(paste0(outdir_fig,
           "/jumping_g000078_highlighted.pdf"), 
    width = 40, height = 16, paper = "special")

layout(matrix(c(1:4), nrow=1), widths=c(1,2,3,3))  # left = tree, right = other plot
par(mar = c(4.5, 0.05,0.05,0))
plot(ST_gene_tree,
     edge.width = 1,                 
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.65,
     no.margin = FALSE, main="")


# then synteny plot
# get y-range
n_genomes = length(dna_seqs)

# get x-range
# x_bp_minmax <- rbindlist(
#   lapply(dna_seqs, function(g) g[, .(start = min(start), end = max(end))]),
#   idcol = "gene_set"
# )[, .(start = min(start), end = max(end))]

x_bp_minmax <- data.frame(start = 1043500,
                          end = 1415000)

par(mar = c(4.5,5,2,0.2))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = unlist(x_bp_minmax),
  ylim = c(0.5, n_genomes + 0.5),
  xlab = "Position (bp)",
  ylab = "", xaxs = "i",
  yaxt = "n", yaxs = "i",
  xaxt = "n"
)

axis(1)

height = 1

# Then add some genes
invisible(
  lapply(seq_len(length(dna_seqs)), function(i){
    genes_tmp = dna_seqs[[i]]
    
    ST_temp = as.character(unique(dna_seqs[[i]]$ST))
    
    gene_width <- mean(genes_tmp$end - genes_tmp$start)
    
    invisible(
      lapply(seq_len(nrow(genes_tmp)), function(j) {
        berryFunctions::roundedRect(genes_tmp$start[j],
                                    i - ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2),
                                    genes_tmp$end[j],
                                    i + ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2), 
                                    col = genes_tmp$fill[j],
                                    rounding = ifelse(genes_tmp$ag_type[j]=="core",
                                                      0, 0.05),
                                    border = genes_tmp$col[j])
      })
    )
    
    # highlight focal gene
    x_focal = genes_tmp[gene_family == "g000077", start]
    points(x_focal, i, pch = 15, cex = 8,
           col = "#00A9CF")
    
    axis(side = 2, at = i, las = 2, col.axis = "navyblue", cex.axis = 0.75,
         labels = ST_temp, tick = FALSE)
    
    
  }
  )
)

x_bp_minmax <- data.frame(start = 1750000,
                          end = 2060000)

par(mar = c(4.5,5,2,0.2))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = unlist(x_bp_minmax),
  ylim = c(0.5, n_genomes + 0.5),
  xlab = "Position (bp)",
  ylab = "", xaxs = "i",
  yaxt = "n", yaxs = "i",
  xaxt = "n"
)

axis(1)

height = 1

# Then add some genes
invisible(
  lapply(seq_len(length(dna_seqs)), function(i){
    genes_tmp = dna_seqs[[i]]
    
    ST_temp = as.character(unique(dna_seqs[[i]]$ST))
    
    gene_width <- mean(genes_tmp$end - genes_tmp$start)
    
    invisible(
      lapply(seq_len(nrow(genes_tmp)), function(j) {
        berryFunctions::roundedRect(genes_tmp$start[j],
                                    i - ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2),
                                    genes_tmp$end[j],
                                    i + ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2), 
                                    col = genes_tmp$fill[j],
                                    rounding = ifelse(genes_tmp$ag_type[j]=="core",
                                                      0, 0.05),
                                    border = genes_tmp$col[j])
      })
    )
    
    # highlight focal gene
    x_focal = genes_tmp[gene_family == "g000077", start]
    points(x_focal, i, pch = 15, cex = 8,
           col = "#00A9CF")
    
    axis(side = 2, at = i, las = 2, col.axis = "navyblue", cex.axis = 0.75,
         labels = ST_temp, tick = FALSE)
    
    
  }
  )
)

# Then add some genes
invisible(
  lapply(seq_len(length(dna_seqs)), function(i){
    genes_tmp = dna_seqs[[i]]
    
    ST_temp = as.character(unique(dna_seqs[[i]]$ST))
    
    gene_width <- mean(genes_tmp$end - genes_tmp$start)
    
    invisible(
      lapply(seq_len(nrow(genes_tmp)), function(j) {
        berryFunctions::roundedRect(genes_tmp$start[j],
                                    i - ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2),
                                    genes_tmp$end[j],
                                    i + ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2), 
                                    col = genes_tmp$fill[j],
                                    rounding = ifelse(genes_tmp$ag_type[j]=="core",
                                                      0, 0.05),
                                    border = genes_tmp$col[j])
      })
    )
    
    # highlight focal gene
    x_focal = genes_tmp[gene_family == "g000077", start]
    points(x_focal, i, pch = 15, cex = 8,
           col = "#00A9CF")
    
    axis(side = 2, at = i, las = 2, col.axis = "navyblue", cex.axis = 0.75,
         labels = ST_temp, tick = FALSE)
    
    
  }
  )
)

x_bp_minmax <- data.frame(start = 2250000,
                          end = 2800000)

par(mar = c(4.5,5,2,0.2))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = unlist(x_bp_minmax),
  ylim = c(0.5, n_genomes + 0.5),
  xlab = "Position (bp)",
  ylab = "", xaxs = "i",
  yaxt = "n", yaxs = "i",
  xaxt = "n"
)

axis(1)

height = 1

# Then add some genes
invisible(
  lapply(seq_len(length(dna_seqs)), function(i){
    genes_tmp = dna_seqs[[i]]
    
    ST_temp = as.character(unique(dna_seqs[[i]]$ST))
    
    gene_width <- mean(genes_tmp$end - genes_tmp$start)
    
    invisible(
      lapply(seq_len(nrow(genes_tmp)), function(j) {
        berryFunctions::roundedRect(genes_tmp$start[j],
                                    i - ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2),
                                    genes_tmp$end[j],
                                    i + ifelse(genes_tmp$ag_type[j]=="core",
                                               height/3,height/2), 
                                    col = genes_tmp$fill[j],
                                    rounding = ifelse(genes_tmp$ag_type[j]=="core",
                                                      0, 0.05),
                                    border = genes_tmp$col[j])
      })
    )
    
    # highlight focal gene
    x_focal = genes_tmp[gene_family == "g000077", start]
    points(x_focal, i, pch = 15, cex = 8,
           col = "#00A9CF")
    
    axis(side = 2, at = i, las = 2, col.axis = "navyblue", cex.axis = 0.75,
         labels = ST_temp, tick = FALSE)
    
    
  }
  )
)

dev.off()
