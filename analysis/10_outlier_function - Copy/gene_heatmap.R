# gene heatmaps for simultaneous genes

# input data --------------------------------------------------------------
clean_coin_ags <- fread(paste0(outdir_dat,"/clean_coin_ags.csv"),
                        select = c("gene_family","p_bh", "w_theta_sel",
                                   "n_genes_assoc_cluster", "score", "clust_id"))

# Pool all p-values, one global BH correction
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

sel_type <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                  select = c("gene_family", "set", "w_theta_sel"))[set==1][, set:=NULL]

gene_pvals <- merge(gene_pvals, sel_type,
                    by = "gene_family")

# set colour map up
gene_pvals[, heat_val:= fcase(w_theta_sel == "balancing" & p_bh <= 0.05, 1,
                              w_theta_sel == "balancing" & p_bh > 0.05, 0.5,
                              w_theta_sel == "directional" & p_bh <= 0.05,-1,
                              w_theta_sel == "directional" & p_bh > 0.05, -0.5,
                              default = 2)]

# load clusters
sim_clus_dt <- fread(paste0(outdir_dat, "/simul_clust_fams.csv"))
sub_clus_dt <- fread(paste0(outdir_dat, "/sub_clust_fams.csv"))

# prescence / absence

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))

# # Generate a binary matrix giving gene presence/absence
# pa <- dcast(
#   data = pan_anno[,.(geno_id, gene_family)],
#   formula = geno_id ~ gene_family,
#   fun.aggregate = length,
#   value.var = "gene_family"
# )
# 
# # checks
# # # Which rows have any value > 1
# # pa[rowSums(pa[, -1] > 1) > 0, 1]
# # 
# # # Which columns have any value > 1
# # names(pa[, -1])[colSums(pa[, -1] > 1) > 0]

# ascribe heatmap colours
gene_cols <- setdiff(colnames(pa), "geno_id")

# # convert pa
# pa[, (gene_cols) := lapply(.SD, as.double), .SDcols = gene_cols]
# 
# for (g in gene_cols) {
#   hv <- gene_pvals[gene_family == g, heat_val]
#   if (length(hv) == 1) {
#     pa[get(g) == 1, (g) := hv]
#   }
#   if (length(hv) == 0) {
#     pa[get(g) == 1, (g) := 2]
#   }
# }

# fwrite(pa, paste0(outdir_dat, "/pa_gf_gene_heatmap.csv"))

pa <- fread(paste0(outdir_dat, "/pa_gf_gene_heatmap.csv"))

# Prepare tree as a dendogram ---------------------------------------------

# Load tree
tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")
if (!is.rooted(tree)) {
  tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
}

# Log  compression
tree$edge.length <- log1p(tree$edge.length)

tree_ultra <- phytools::force.ultrametric(tree) 

hc <- as.hclust(tree_ultra)
# Convert to dendrogram
dend <- as.dendrogram(hc)


# Cluster heat map --------------------------------------------------------

plot_gene_heatmap <- function(clst_id, dend, gene_pvals, pa, tree){
  
  # set colour conditions
  breaks <- c(-1.5, -0.75, -0.25, 0.25, 0.75, 1.5, 2.5)
  cols   <- c("#F26C2A", "#F9A876", "white", "#4DB3FF", "#1A7ACC", "grey70")
  
  out_grp_cols <- setNames(cols, c(-1,-0.5,0,0.5,1,2))

  temp_clst_dt <- clean_coin_ags[clust_id == clst_id]
  
  # focal gene
  focal_out_ag <- temp_clst_dt[p_bh == min(p_bh), gene_family][1]
  
  # focal gene colour
  focal_ag_col <- gene_pvals[gene_family==focal_out_ag, heat_val]
  
  present <- pan_anno[gene_family == focal_out_ag, geno_id]
  
  dend %>% set("by_labels_branches_col", 
               value = present, 
               TF_values = c(out_grp_cols[[as.character(focal_ag_col)]],Inf)) -> dend_col
  
  # plot(dend_col, horiz = TRUE, leaflab = "none")
  
  # get simultaneous ags first
  sim_ags <- if (length( temp_clst_dt[score=="simultaneous", gene_family])>0){
    temp_clst_dt[score=="simultaneous"][order(p_bh), gene_family]
  } else{NULL}
  
  sim_clus_num <- if (length(sim_ags) > 0) {
    unique(sim_clus_dt[gene_family %chin% sim_ags, assoc_cluster])
  } else { NULL }
  
  othr_sim_ags <- if (!is.null(sim_clus_num) && length(sim_clus_num) > 0) {
    setdiff(sim_clus_dt[assoc_cluster %chin% as.character(sim_clus_num), gene_family], sim_ags)
  } else { NULL }
  
  sub_ags <- if (length(temp_clst_dt[score!="simultaneous", gene_family])>0){
    temp_clst_dt[score!="simultaneous"][order(p_bh), gene_family]
  } else{NULL}
  
  sub_clus_num <- if (length(sub_ags) > 0) {
    unique(sub_clus_dt[gene_family %chin% sub_ags, assoc_cluster])
  } else { NULL }
  
  othr_sub_ags <- if (!is.null(sub_clus_num) && length(sub_clus_num) > 0) {
    setdiff(sub_clus_dt[assoc_cluster %chin% as.character(sub_clus_num), gene_family], sub_ags)
  } else { NULL }
  
  # put everything together, priotitising sim
  temp_coin_genes <- unique(c(sim_ags, othr_sim_ags, 
                              setdiff(sub_ags, c(sim_ags, othr_sim_ags)),
                              setdiff(othr_sub_ags, c(sim_ags, othr_sim_ags))))
  
  
  temp_mat <- pa[, c("geno_id", temp_coin_genes), with = FALSE]
  
  temp_mat[, geno_id := factor(geno_id, levels = tree$tip.label)]
  setorder(temp_mat, geno_id)
  
  mat <- as.matrix(temp_mat[,-1])
  row.names(mat) <- temp_mat$geno_id
  
  heatmap(mat, Colv = NA, Rowv = dend_col,
          scale = "none",
          labRow = NA, 
          col    = cols,        
          cexCol = 1, 
          main = paste("Cluster", clst_id, ": type", paste(unique(temp_clst_dt[order(score), score]), collapse = "/")), 
          margins = c(4, 2),
          breaks = breaks)

  
}

# for (i in unique(clean_coin_ags$clust_id)) {
#   png(paste0(outdir_fig, "/clusters_heatmaps/cluster_",i,".png"),
#       width  = 13.33 * (4/7),   # 7.62 inches
#       height = 7.5  * (4/5),    # 6 inches
#       units  = "in",
#       res    = 150,              # 300 for print quality
#       bg     = "white")
#   
#   plot_gene_heatmap(i, dend, gene_pvals, pa, tree)
#   
#   dev.off()
# }
# 
