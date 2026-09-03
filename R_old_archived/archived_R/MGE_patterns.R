outlier_ags <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                     select = c("gene_family", "freq", "set", "seg_s_sites", 
                                "p_two", "p_bh", "theta_outlier", "w_theta_sel"))

# keep only outliers
outlier_ags <- outlier_ags[theta_outlier==1 & p_bh <= 0.05][,theta_outlier := NULL]

outlier_ags <- outlier_ags[, if (any(set == 2)) .SD[set == 2] else .SD, by = gene_family]

# add in gene information
pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("gene_family", "fus_locus_tag")))

outlier_ags <- merge(outlier_ags, pangraph_anno,
                     all.x = TRUE, by = c("gene_family"))

# get full ag info

pirate_full_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                                 select = c("fus_locus_tag","geno_id","gene_family","alleles_at_maximum_threshold",
                                            "phg_state", "most_common_phage", 
                                            "integron","ice","is_fam", "is_cluster","consensus_gene_name",
                                            "anchor","msu", "consensus_product", "egng_cog")))

anno <- unique(pirate_full_anno[,.(gene_family,consensus_gene_name, consensus_product, egng_cog)])

outlier_ags <- merge(outlier_ags, pirate_full_anno,
                     all.x = TRUE, by = c("fus_locus_tag"))

outlier_ags <- unique(outlier_ags[, -c("fus_locus_tag", "geno_id")])

fwrite(outlier_ags, "C:/Users/carac/Desktop/outlier_mobilome.csv")



outlier_ags[ice!=""]



outlier_ags <- merge(outlier_ags, anno,
                     all.x = TRUE, by = "gene_family")


fwrite(outlier_ags[, .(gene_family,consensus_gene_name,consensus_product, egng_cog )], "C:/Users/carac/Desktop/anno_outlier_mobilome.csv")


pirate_full_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                                 select = c("gene_family","anchor","msu")))

sites = unique(pirate_full_anno[,.(msu, anchor)])
sites = sites[grepl(":", anchor)]

sites = sites[,sites:= paste(msu, anchor, sep = "|"), .I]


hostspot <- fread("./output/data/msu_site_classification.csv")

hostspot[, phage_n := round(n_families * phage_frac, digits = 0)]

hostspot[, ag_n := n_families - (ice_n + phage_n)]

data <- as.matrix(hostspot[, .(ag_n, ice_n, phage_n)])
colnames(data) <- c("Non-MGEs", "ICE/IME", "Phages")
rownames(data) <- hostspot$site

col = c("#3A7CB5", 
        "#F0A860", 
        "#C0392B") 

# Get the stacked barplot
par(mar = c(9,4,5,2), xpd = TRUE)
barplot(t(data), 
        col= col, 
        border="white", 
        space=0.04, 
        font.axis=1,
        ylab = "Number of Accessory Genes",
        cex.axis = 0.85,
        las = 2,
        xlab="")

legend("topright",
        legend = colnames(data),              # the 3 stacked categories
        fill   = col,       # same colours, same order
        border = "white",
        bty    = "n",
        cex    = 0.85)


labs <- hostspot[,.(n_families, site_class, dominant_function)][order(-n_families),]
labs[, y:=.I]

with(labs[n_families %in% c(40,26,18,16,15,11,6,4) & !y %in% c(6, 18,20:21,23:24)],
     text(y-1, n_families+0.2,  labels = paste0(site_class," (",dominant_function, ")"),
          srt = 20, cex = 0.8, pos = 4))


