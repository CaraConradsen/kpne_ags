
all_associations <- unique(fread("./output/data/sorted_gf_clust.csv",
                                 select = c("gene_family", "score", "clust_id",  "p_bh", 
                                            "w_theta_sel")))

all_associations <- all_associations[p_bh <= 0.05]

all_associations <- all_associations[, .(score = unlist(strsplit(score, ":"))), 
                                     by = .(gene_family, clust_id, p_bh, w_theta_sel)]


# get pa and tree ---------------------------------------------------------------

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))

# Generate a binary matrix giving gene presence/absence
pa <- dcast(
  data = pan_anno[,.(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

# Load tree
tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")


# D stat ------------------------------------------------------------------


ag_presence_absence <- as.data.frame(pa)

com_dat_ag_ap <- comparative.data(tree, ag_presence_absence, geno_id)

# simultaneous
DEstimates_sim_out <- lapply(unique(all_associations[score== "simultaneous", gene_family]), function(g){
  d_info <- do.call(caper::phylo.d,
                    list(data = com_dat_ag_ap, binvar = as.name(g)))
  data.table(gene_family = d_info$binvar, d = d_info$DEstimate,
             p_rand = d_info$Pval1, p_brow = d_info$Pval0)
})

DEstimates_sim_out <- rbindlist(DEstimates_sim_out)

DEstimates_sim_out <- merge(DEstimates_sim_out,all_associations[score== "simultaneous",],
                            by = "gene_family")




# subsequent
DEstimates_sub_out <- lapply(unique(all_associations[score== "subsequent", gene_family]), function(g){
  d_info <- do.call(caper::phylo.d,
                    list(data = com_dat_ag_ap, binvar = as.name(g)))
  data.table(gene_family = d_info$binvar, d = d_info$DEstimate,
             p_rand = d_info$Pval1, p_brow = d_info$Pval0)
})

DEstimates_sub_out <- rbindlist(DEstimates_sub_out)

DEstimates_sub_out <- merge(DEstimates_sub_out,all_associations[score== "subsequent",],
                            by = "gene_family")

DEstimates_sub_out$score = "subsequent"

fwrite(rbind(DEstimates_sim_out, DEstimates_sub_out),
       paste0(outdir_dat, "/DEstimates.csv"))