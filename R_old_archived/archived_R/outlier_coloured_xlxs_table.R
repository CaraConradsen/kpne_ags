# create and export an xlxs worksheet for the supplemental

outlier_ags <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                     select = c("gene_family", "freq", "set", "seg_s_sites", 
                                "p_two", "p_bh", "theta_outlier", "w_theta_sel"))

# keep only outliers
outlier_ags <- outlier_ags[theta_outlier==1 & p_bh <= 0.05][,theta_outlier := NULL]

outlier_ags <- outlier_ags[, if (any(set == 2)) .SD[set == 2] else .SD, by = gene_family]

# add in gene information
pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("gene_family", "egng_cog", "ag_type","gene_length",
                                         "consensus_gene_name", "consensus_product")))

# take longest gene
pangraph_anno <- pangraph_anno[, .SD[which.max(gene_length)], by = gene_family]

outlier_ags <- merge(outlier_ags,
                     pangraph_anno,
                     all.x  = TRUE, 
                     by = "gene_family")

setcolorder(outlier_ags, c("w_theta_sel","gene_family","freq","gene_length",
                           "ag_type","set","seg_s_sites", "p_two","p_bh",        
                           "egng_cog", "consensus_gene_name","consensus_product"))

setorderv(outlier_ags, cols = c("w_theta_sel", "set", "freq"),
          order = c(-1, 1, 1))

outlier_ags[, p_two := sprintf("%1.4f", p_two)]
outlier_ags[, p_bh  := sprintf("%1.4f", p_bh)]
# fix small vals
outlier_ags[p_two == "0.0000", p_two := "<0.0001"]
outlier_ags[p_bh == "0.0000", p_bh := "<0.0001"]

setnames(outlier_ags,
         c("w_theta_sel","gene_family","freq","ag_type","set","seg_s_sites","gene_length",
           "p_two","p_bh","egng_cog","consensus_gene_name","consensus_product"),
         c("Selection type","Gene family","Gene frequency",
           "Pangenome category","Set","Downsampled Seg. Sites","Longest gene length",
           "P-value","FDR P","COG category","Gene name","Annotation"))


# COGS --------------------------------------------------------------------
cog_set <- fread(paste0(outdir_dat, "/bh_p_cog_set.csv"),
                 select = c("COG_letter", "gene_family","signed_logp", "COG_function"))[COG_letter!=""]

setnames(cog_set,c("COG_letter","gene_family","signed_logp", "COG_function"),
         c("COG", "Gene family","Signed log10(p)","COG function"))

cog_set <- merge(cog_set, outlier_ags[,-"COG category"],
                 all.x = TRUE, by = "Gene family")


# associations ------------------------------------------------------------
clean_coin_ags <- fread(paste0(outdir_dat,"/clean_coin_ags.csv"),
                        select = c("gene_family","p_bh", "w_theta_sel", "set","perc_in_set1",
                                   "n_genes_assoc_cluster", "score", "clust_id"))
# Id the set 2 AG
coin_ag_set2 <- clean_coin_ags[set == 2, gene_family]

clean_coin_ags[, set:=NULL]

# fix overlappers
clean_coin_ags[, score := paste(sort(unique(score)), collapse = ":"), 
               by = .(gene_family, clust_id)]

# Collapse to one row per gene_family/clust_id
clean_coin_ags <- unique(clean_coin_ags, by = c("gene_family", "clust_id"))

# Add other genes
all_gf_ags <- rbind(fread(paste0(outdir_dat, "/sub_clust_fams.csv")), 
                    fread(paste0(outdir_dat, "/simul_clust_fams.csv")))
coin_genes <- unique(clean_coin_ags$gene_family)

sorted_gf_clust <- lapply(coin_genes, function(i) {
  temp_focal <- clean_coin_ags[gene_family == i]
  
  temp <- all_gf_ags[gene_family == i]
  
  temp <- split(temp, seq_len(nrow(temp)))
  
  assoc_genes <- lapply(temp, function(clust_vals){
    all_gf_ags[assoc_cluster == clust_vals$assoc_cluster & score == clust_vals$score]
  })
  
  assoc_genes <- rbindlist(assoc_genes)[,.(gene_family, score)]
  
  assoc_genes$clust_id = temp_focal$clust_id
  
  # don't self-match
  assoc_genes <- assoc_genes[!gene_family %chin% temp_focal$gene_family]
  
  rbindlist(list(temp_focal[,-c("p_bh","w_theta_sel")], assoc_genes), fill = TRUE)
})

sorted_gf_clust <- rbindlist(sorted_gf_clust)

# add p-vals and sets
gene_pvals <- fread("./output/data/gene_pvals_fdr.csv", 
                    select = c("gene_family", "seg_s_sites","p_bh"))

sorted_gf_clust <- merge(sorted_gf_clust, gene_pvals,
                         all.x = TRUE, by = "gene_family")

outlier_ags_infor <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                     select = c("gene_family", "set", "w_theta_sel"))

outlier_ags_infor <- outlier_ags_infor[, .SD[which.max(set)], by = gene_family]

sorted_gf_clust <- merge(sorted_gf_clust, outlier_ags_infor,
                         all.x = TRUE, by = "gene_family")

# add gene and product
sorted_gf_clust <- merge(sorted_gf_clust, pangraph_anno[, .(gene_family, gene_length, consensus_gene_name, consensus_product)],
                         all.x = TRUE, by = "gene_family")

msu_anchor <- unique(fread("./output/data/msu_regions_anchored.csv", 
                    select = c("gene_family", "number_genomes", "msu", "anchor")))

sorted_gf_clust <- merge(sorted_gf_clust, msu_anchor,
                         all.x = TRUE, by = "gene_family")

# remove dups
sorted_gf_clust <- unique(sorted_gf_clust)

setorderv(sorted_gf_clust, cols = c("clust_id", "p_bh"),
          order = c(1, 1))

fwrite(sorted_gf_clust, paste0(outdir_dat,"/sorted_gf_clust.csv"))

setnames(sorted_gf_clust, 
         old = c("gene_family", "perc_in_set1", "n_genes_assoc_cluster", "score", 
                 "clust_id",  "seg_s_sites", "p_bh", "set", "w_theta_sel", 
                 "gene_length", "consensus_gene_name", "consensus_product", 
                 "number_genomes", "msu", "anchor"),
         new = c("Gene family", "% clust in Set 1", "Cluster size", "Scoring method", 
                 "Cluster", "Downsampled Seg. Sites", "FDR p-value", 
                 "Set", "Selection type", "Longest gene length", "Gene name", "Annotation", 
                 "Gene frequency",  "MSU", "Anchor"))


# Create a coloured xlxs workbook for analysis ----------------------------

wb <- openxlsx::createWorkbook()

# Define styles
orange_sig    <- createStyle(fgFill = "#F0A855")
blue_sig      <- createStyle(fgFill = "#6FA8D4")
orange_ns_sig <- createStyle(fgFill = "#F5C897")
blue_ns_sig   <- createStyle(fgFill = "#A8C8E8")
ns_in_set     <- createStyle(fgFill = "grey80")
bold_header   <- createStyle(textDecoration = "bold")


# Sheet 1: Data dictionary ------------------------------------------------
dict <- data.table(
  Sheet = c(
    rep("BH corrected selection outliers", 12),
    rep("COGs", 15),
    rep("Gene associations", 17)
  ),
  Column = c(
    "Selection type", "Gene family", "Gene frequency", "Longest gene length",
    "Pangenome category", "Set", "Downsampled Seg. Sites", "P-value", "FDR P",
    "COG category", "Gene name", "Annotation",
    "COG", "Gene family", "Signed log10(p)", "COG function",
    "Selection type", "Gene frequency", "Longest gene length",
    "Pangenome category", "Set", "Downsampled Seg. Sites", "P-value", "FDR P",
    "Gene name", "Annotation",
    "Gene family", "% clust in Set 1", "Cluster size", "Scoring method",
    "Cluster", "Downsampled Seg. Sites", "FDR p-value", "Set", "Selection type",
    "Longest gene length", "Gene name", "Annotation", "Gene frequency",
    "Sequence type", "MSU", "Anchor", "Gene frequency"
  ),
  Description = c(
    "Mode of selection detected: balancing or directional",
    "PIRATE gene family identifier",
    "Number of isolates carrying this gene family",
    "Length (bp) of the longest representative gene in the family",
    "Pangenome category: core, shell, or cloud",
    "Filter set gene belongs to: 1 (syntenic) or 2 (single-gain)",
    "Synonymous segregating sites after downsampling",
    "Two-tailed p-value from the selection test",
    "Benjamini-Hochberg corrected p-value",
    "COG functional category letter code",
    "Consensus gene name from PIRATE annotation",
    "Consensus product annotation from PIRATE",
    "COG functional category letter code",
    "PIRATE gene family identifier",
    "Signed log10 p-value; sign indicates balancing (+) or directional (-)",
    "Full name of COG functional category",
    "Mode of selection detected: balancing or directional",
    "Number of isolates carrying this gene family",
    "Length (bp) of the longest representative gene in the family",
    "Pangenome category: core, shell, or cloud",
    "Filter set gene belongs to: 1 (syntenic) or 2 (single-gain)",
    "Synonymous segregating sites after downsampling",
    "Two-tailed p-value from the selection test",
    "Benjamini-Hochberg corrected p-value",
    "Consensus gene name from PIRATE annotation",
    "Consensus product annotation from PIRATE",
    "PIRATE gene family identifier",
    "Percentage of genes in the co-occurrence cluster present in Set 1",
    "Number of gene families in the co-occurrence cluster",
    "Goldfinder scoring method: simultaneous, subsequent, or both",
    "Cluster identifier assigned after merging simultaneous and subsequent clusters",
    "Synonymous segregating sites after downsampling",
    "Benjamini-Hochberg corrected p-value from the selection test",
    "Filter set gene belongs to: 1 (syntenic) or 2 (single-gain)",
    "Mode of selection detected: balancing or directional",
    "Length (bp) of the longest representative gene in the family",
    "Consensus gene name from PIRATE annotation",
    "Consensus product annotation from PIRATE",
    "Number of isolates carrying this gene family",
    "MLST sequence type of the isolate",
    "Mobilisation/synteny unit region identifier",
    "Core-block flanking gene pair defining genomic anchor position",
    "Number of isolates carrying this gene family"
  ),
  Colour = c(
    rep("Orange (dark): directional outlier; Blue (dark): balancing outlier", 12),
    rep("Orange (dark): directional outlier; Blue (dark): balancing outlier", 15),
    rep("Orange (dark): directional, FDR <= 0.05; Orange (light): directional, FDR > 0.05; Blue (dark): balancing, FDR <= 0.05; Blue (light): balancing, FDR > 0.05; Grey: not a selection outlier", 17)
  )
)

addWorksheet(wb, "Data dictionary")
openxlsx::writeData(wb, "Data dictionary", dict)
addStyle(wb, "Data dictionary", bold_header,
         rows = 1, cols = 1:ncol(dict), gridExpand = TRUE)

# Sheet 2: BH corrected selection outliers --------------------------------
addWorksheet(wb, "BH corrected selection outliers")
openxlsx::writeData(wb, "BH corrected selection outliers", outlier_ags)

dir_rows <- which(outlier_ags[["Selection type"]] == "directional") + 1
bal_rows <- which(outlier_ags[["Selection type"]] != "directional") + 1

addStyle(wb, "BH corrected selection outliers", orange_sig, rows = dir_rows,
         cols = 1:ncol(outlier_ags), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "BH corrected selection outliers", blue_sig, rows = bal_rows,
         cols = 1:ncol(outlier_ags), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "BH corrected selection outliers", bold_header,
         rows = 1, cols = 1:ncol(outlier_ags), gridExpand = TRUE)

# Sheet 2: COGs -----------------------------------------------------------
addWorksheet(wb, "COGs")
openxlsx::writeData(wb, "COGs", cog_set)

cog_dir_rows <- which(cog_set[["Selection type"]] == "directional") + 1
cog_bal_rows <- which(cog_set[["Selection type"]] != "directional") + 1

addStyle(wb, "COGs", orange_sig, rows = cog_dir_rows,
         cols = 1:ncol(cog_set), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "COGs", blue_sig, rows = cog_bal_rows,
         cols = 1:ncol(cog_set), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "COGs", bold_header,
         rows = 1, cols = 1:ncol(cog_set), gridExpand = TRUE)

# Sheet 3: Gene associations ----------------------------------------------
addWorksheet(wb, "Gene associations")
openxlsx::writeData(wb, "Gene associations", sorted_gf_clust)

sel  <- sorted_gf_clust[["Selection type"]]
fdr  <- sorted_gf_clust[["FDR p-value"]]
set  <- sorted_gf_clust[["Set"]]

na_rows      <- which(is.na(sel) & !is.na(set)) + 1
sig_dir_rows <- which(!is.na(sel) & sel == "directional" & !is.na(fdr) & fdr <= 0.05) + 1
ns_dir_rows  <- which(!is.na(sel) & sel == "directional" & (is.na(fdr) | fdr >  0.05)) + 1
sig_bal_rows <- which(!is.na(sel) & sel != "directional" & !is.na(fdr) & fdr <= 0.05) + 1
ns_bal_rows  <- which(!is.na(sel) & sel != "directional" & (is.na(fdr) | fdr >  0.05)) + 1

addStyle(wb, "Gene associations", ns_in_set,     rows = na_rows,      cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Gene associations", orange_sig,    rows = sig_dir_rows, cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Gene associations", orange_ns_sig, rows = ns_dir_rows,  cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Gene associations", blue_sig,      rows = sig_bal_rows, cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Gene associations", blue_ns_sig,   rows = ns_bal_rows,  cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "Gene associations", bold_header,
         rows = 1, cols = 1:ncol(sorted_gf_clust), gridExpand = TRUE)


saveWorkbook(wb, paste0(outdir_xlsx, "/outliers_ag.xlsx"), overwrite = TRUE)
