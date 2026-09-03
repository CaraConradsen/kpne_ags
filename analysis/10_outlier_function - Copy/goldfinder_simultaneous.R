# Analyse the output from Goldfinder


# Part 1 Simultaneous estimates of association ----------------------------

ass_clusters <- read.delim("./input_data/goldfinder/simul_output/association_clusters.txt",
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

# assign whether dis
dis_pairs <- fread("./input_data/goldfinder/simul_output/simultaneous_dissociation_significant_pairs.csv")
dis_pairs[, Gene_1:= tstrsplit(Gene_1, "__", fill = TRUE, keep = 2)]
dis_pairs[, Gene_2:= tstrsplit(Gene_2, "__", fill = TRUE, keep = 2)]

all_clus_dt <- rbind(a_clus_dt, data.table(assoc_cluster = max(a_clus_dt$assoc_cluster)+1,
                                           gene_family = c(dis_pairs$Gene_1, dis_pairs$Gene_2),
                                           n_genes_assoc_cluster = 2,
                                           assc_type = "disassociate"))

# Find gene outlier associations 
ag_outliers <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]

ag_outliers <- merge(ag_outliers, all_clus_dt,
                     all.x = TRUE, by = "gene_family")

ag_synteny <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

associated_outliers <- ag_outliers[!is.na(assoc_cluster)][order(assoc_cluster)]

associated_outliers[, c("expected_S", "seg_s_sites","p",
                        "cdf","p_lower","p_upper") := NULL]

setnames(associated_outliers, "freq", "number_genomes")

# find synteny
syn_prod_info <- lapply(unique(associated_outliers$assoc_cluster), function(clust) {
  
  # find associated gene families 
  temp_dat <- unique(ag_synteny[gene_family %in% all_clus_dt[assoc_cluster==clust, gene_family], 
                                .(gene_family, number_genomes, ag_type, msu, anchor, consensus_product)])
  
  # return info 
  cbind(assoc_cluster = clust, temp_dat)
})

syn_prod_info <- rbindlist(syn_prod_info)


associated_outliers <- merge(associated_outliers, syn_prod_info,
                             all = TRUE, 
                             by = c("assoc_cluster", "gene_family", "number_genomes"))

associated_outliers[, n_genes_assoc_cluster := max(n_genes_assoc_cluster, na.rm = TRUE), assoc_cluster]

# write.xlsx(associated_outliers, paste0(outdir_xlsx, "/AG_outliers.xlsx"), sheetName="sim_association", append = TRUE)


# Add in Cogs -------------------------------------------------------------
# don't double count genes!!

pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("geno_id", "gene_family", 
                                         "fus_locus_tag","egng_cog", 
                                         "consensus_gene_name", "consensus_product")))

# rename eggnog_cog to COG_funct_cat
setnames(pangraph_anno, old = "egng_cog", new = "COG_funct_cat")

gene_pvals <- fread(paste0(outdir_dat, "/set_outliers.csv"))[theta_outlier==1]


# Calculate signed log10(p) -----------------------------------------------
# Join neutral quantiles and define outliers
# Poisson-Based Approaches
# Since segregating sites are count data, Poisson tests are statistically appropriate:
# Calculate the mean from neutral genes (λ)
# Use Poisson probability: P(X ≥ observed | λ)
# Convert to a score: -log10(p-value)

# signed magnitude: positive = excess, negative = deficit
gene_pvals[, direction := fifelse(seg_s_sites > expected_S, 1, -1)]
gene_pvals[, signed_logp := direction * -log10(p_two)]

# Plot by functional category -----------------------------------------
cog_set <- gene_pvals[,.(gene_family, set, signed_logp)]

# add cogs
cog_set[pangraph_anno[,.(gene_family, COG_funct_cat)], on = "gene_family",
        COG := COG_funct_cat
]

cog_set <- rbind(cog_set[, .(COG_letter = unlist(strsplit(COG, ""))), by = c("gene_family","set","signed_logp")],
                 cog_set[COG=="",.(gene_family, set, signed_logp, COG_letter= COG)])

# add cog annotation from COG_dt
cog_set <- merge(cog_set, COG_dt[,.(COG_letter,COG_function, desc)],
                 all.x = TRUE, by = "COG_letter")

cog_set <- merge(cog_set,
                 gene_pvals[, .(gene_family, set, seg_s_sites, p_two, p_bh, w_theta_sel)],
                 all.x = TRUE, by = c("gene_family", "set"))


cog_set <- merge(cog_set,
                 unique(pangraph_anno[,.(gene_family, consensus_gene_name, consensus_product)]),
                 all.x = TRUE, by = "gene_family")

# write.xlsx(cog_set, paste0(outdir_xlsx, "/AG_outliers.xlsx"), sheetName="cog_set")



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

DEstimates_ag_out_cl <- merge(DEstimates_ag_out,all_clus_dt,
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

pdf(paste0(outdir_fig, "/d_vs_associations_clusters.pdf"), 
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

arrows(x0 = -0.35, y0 = 1.25,
       x1 = -0.05, y1 = 1.25,    
       code = 1,             
       angle = 20,           
       length = 0.1,         
       col = "grey30",
       lwd = 2,
       xpd = NA)            

text(-0.2, 1.25, cex = 0.75, 
     col = "grey30", pos= 1,
     "stronger\nphylognetic signal", xpd = TRUE)

arrows(x0 = 0.2, y0 = 1.35,
       x1 = 0.7, y1 = 1.35,    
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
