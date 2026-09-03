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





# # Create a xlxs workbook for supp ----------------------------------------

wb <- openxlsx::createWorkbook()
addWorksheet(wb, "outliers")
openxlsx::writeData(wb, "outliers", outlier_ags)

# Bold headers
addStyle(wb, "outliers", createStyle(textDecoration = "bold"), rows = 1, cols = 1:ncol(outlier_ags), gridExpand = TRUE)

saveWorkbook(wb, paste0(outdir_xlsx, "/outliers_ag_supp.xlsx"), overwrite = TRUE)


