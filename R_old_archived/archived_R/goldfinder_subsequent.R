# Analyse the output from Goldfinder


# Part 1 Simultaneous estimates of association ----------------------------

ass_clusters <- read.delim("./input_data/goldfinder/sub_output/association_clusters.txt",
                           header = FALSE)

header_idx <- which(startsWith(ass_clusters$V1, ">"))

# end of each cluster = one before the next header; last cluster runs to end of df
starts <- header_idx + 1
ends   <- c(header_idx[-1] - 1, nrow(ass_clusters))

# list of row indices per cluster
member_rows <- Map(seq, starts, ends)

a_clus_dt <- lapply(seq_along(member_rows), function(g) {
  assoc_cluster = gsub(">","", ass_clusters$V1[header_idx[g]])
  gene_family = ass_clusters$V1[member_rows[[g]]]
  
  data.table(cbind(assoc_cluster,gene_family))
})

a_clus_dt <- rbindlist(a_clus_dt)

a_clus_dt[,c("assoc_cluster", "n_genes_assoc_cluster"):= tstrsplit(assoc_cluster, ",", fill = TRUE)]
a_clus_dt[, c("assoc_cluster", "n_genes_assoc_cluster") := lapply(.SD, as.numeric), .SDcols = c("assoc_cluster", "n_genes_assoc_cluster")]
a_clus_dt[, gene_family:= gsub(",", "", gene_family)]
a_clus_dt[, gene_family:= tstrsplit(gene_family, "__", fill = TRUE, keep = 2)]

a_clus_dt$assc_type = "associate"


# Find gene outlier associations 
ag_outliers <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]

ag_outliers <- merge(ag_outliers, a_clus_dt,
                     all.x = TRUE, by = "gene_family")

ag_synteny <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

associated_outliers <- ag_outliers[!is.na(assoc_cluster)][order(assoc_cluster)]

associated_outliers[, c("expected_S", "seg_s_sites","p",
                        "cdf","p_lower","p_upper") := NULL]

setnames(associated_outliers, "freq", "number_genomes")

# find synteny
syn_prod_info <- lapply(unique(associated_outliers$assoc_cluster), function(clust) {
  
  # find associated gene families 
  temp_dat <- unique(ag_synteny[gene_family %in% a_clus_dt[assoc_cluster==clust, gene_family], 
                                .(gene_family, number_genomes, ag_type, msu, anchor, consensus_product)])
  
  # return info 
  cbind(assoc_cluster = clust, temp_dat)
})

syn_prod_info <- rbindlist(syn_prod_info)


associated_outliers <- merge(associated_outliers, syn_prod_info,
                             all = TRUE, 
                             by = c("assoc_cluster", "gene_family", "number_genomes"))

associated_outliers[, n_genes_assoc_cluster := max(n_genes_assoc_cluster, na.rm = TRUE), assoc_cluster]

write.xlsx(associated_outliers, paste0(outdir_xlsx, "/AG_outliers.xlsx"),
           sheetName="sub_association", append = TRUE)


# Add in Cogs -------------------------------------------------------------
# don't double count genes!!

pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("geno_id", "gene_family", 
                                         "fus_locus_tag","egng_cog", 
                                         "consensus_gene_name", "consensus_product")))

# rename eggnog_cog to COG_funct_cat
setnames(pangraph_anno, old = "egng_cog", new = "COG_funct_cat")

gene_pvals <- fread(paste0(outdir_dat, "/set_outliers.csv"))[theta_outlier==1]


# Cluster vs D ------------------------------------------------------------


# Get remaining accessory genes after synteny analysis --------------------
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

ag_out_dat <- msu_regions_anchored[gene_family %chin% ag_outliers$gene_family]

core_gub_tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")

# create presence/absence for gene_families -------------------------------
# Create binary ag presence_absence data.table 
# needs a dummy start and end, where rows are each gene family and columns are genomes

# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_out_dat[,.(geno_id, gene_family)],
  formula = geno_id~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

ag_presence_absence <- as.data.frame(ag_presence_absence)

com_dat_ag_ap <- comparative.data(core_gub_tree, ag_presence_absence, geno_id)

DEstimates_ag_out <- lapply(unique(ag_out_dat$gene_family), function(g){
  d_info <- do.call(caper::phylo.d,
                    list(data = com_dat_ag_ap, binvar = as.name(g)))
  data.table(gene_family = d_info$binvar, d = d_info$DEstimate,
             p_rand = d_info$Pval1, p_brow = d_info$Pval0)
})

DEstimates_ag_out <- rbindlist(DEstimates_ag_out)

DEstimates_ag_out_cl <- merge(DEstimates_ag_out,a_clus_dt,
                           by = "gene_family")


# plot --------------------------------------------------------------------

tab <- merge(DEstimates_ag_out_cl[,.(gene_family,assoc_cluster, p_brow)],
             unique(ag_synteny[,.(gene_family,number_genomes, consensus_product)]),
             all.x= TRUE, by ="gene_family")

tab <- merge(tab, ag_outliers[,.(gene_family, seg_s_sites, w_theta_sel, set)],
             by ="gene_family")

tab <- rbind(tab[set==2],
             tab[!gene_family %chin% tab[set==2, gene_family]])

tab <- unique(tab[, .(gene_family, set, assoc_cluster,number_genomes,   
               consensus_product, seg_s_sites, w_theta_sel, p_brow)])

tab <- tab[order(p_brow)]

colnames(tab) <- c("AG", "Set", "Clust.", "Pan. freq.", "Gene product",
                   "Syn. seg. sites", "Selection", "p_brow")


#plot

pdf(paste0(outdir_fig, "/d_vs_associations_sub_clusters.pdf"), 
    width = 12, height = 9, 
    pointsize = 12)

layout(matrix(c(1, 2), nrow = 1), widths = c(0.45, 0.55))
with(DEstimates_ag_out_cl,
     plot(d, -log10(p_brow),
     pch = 16, bty = "L",
     ylab = expression("-log"[10]~italic("P")*"(Brownian)"),
     xlab = "Phylogenetic D statistic",
     col = ifelse(-log10(p_brow) <= -log10(0.05),
                  rgb(0.1,0.1,0.1, alpha = 0.6),
                  rgb(0.1,0.8,0.1, alpha = 0.6))))

with(DEstimates_ag_out_cl,
     thigmophobe.labels(d, -log10(p_brow),
                        labels = assoc_cluster, 
                        cex = 0.5))

abline(v = 0, lty= 2, col = "grey40")
abline(h = -log10(0.05), lty= 2, col = "grey40")

arrows(x0 = -0.6, y0 = 1.25,
       x1 = -0.05, y1 = 1.25,    
       code = 1,             
       angle = 20,           
       length = 0.1,         
       col = "grey30",
       lwd = 2,
       xpd = NA)            

text(-0.35, 1.25, cex = 0.75, 
     col = "grey30", pos= 1,
     "stronger\nphylognetic signal", xpd = TRUE)

arrows(x0 = 0.3, y0 = 1.35,
       x1 = 0.55, y1 = 1.35,    
       code = 3,             
       angle = 20,           
       length = 0.1,         
       col = rgb(0.1,0.8,0.1, alpha = 1),
       lwd = 2,
       xpd = NA)            

text(0.45, 1.35, cex = 0.75, 
     col = "grey30", pos= 3,
     "HGT", xpd = TRUE)



plot.new()
sig_cl <- nrow(tab[p_brow<=0.05])
ns_cl <- nrow(tab)- sig_cl
num_cols = ncol(tab)-1

abg <- matrix(c(rep(rgb(0.1,0.8,0.1, alpha = 0.6), sig_cl*num_cols),
                rep(rgb(1,1,1, alpha = 0.6), ns_cl*num_cols)),
              nrow = nrow(tab), byrow = TRUE)

plot.window(xlim = c(0, 1), ylim = c(0, 1))
addtable2plot(-0.1,-0.15, tab[,-"p_brow"],
              bty = "o", display.rownames = FALSE,
              cex = 0.55, bg = abg, xpad = 0,ypad = 0.95,
              hlines = TRUE, vlines = FALSE,
              title = expression("Outlier genes ranked by -log"[10]~italic("P")*"(Brownian)"))

dev.off()
