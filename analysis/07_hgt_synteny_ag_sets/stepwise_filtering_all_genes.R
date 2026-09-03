
id_ags <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("gene_family","number_genomes", "average_dose", 
                                  "ag_type", "acrs_msu", "acrs_jun")))

dt <- copy(id_ags)

# pi outliers
piS_3IQR_outliers <- fread(paste0(outdir_dat, "/piS_3IQR_outliers.csv"))

# syntenic ags
ag_names <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"), 
                   select = c("gene_family"))$gene_family 

syntenic_ags <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                             select = c("gene_family", "acrs_msu", "acrs_jun","msu"))

syntenic_ags <- syntenic_ags[gene_family %chin% ag_names][msu!=""]

syntenic_ags <- unique(syntenic_ags)

syntenic_ags <- syntenic_ags[acrs_msu!= 1][acrs_jun!=1]$gene_family

# Step 1 — Total
dt1 <- dt
n1 <- nrow(dt1)

# Step 2 — Remove core
dt2 <- dt1[ag_type != "core"]
n2 <- nrow(dt2)

# Step 3 — Remove paralogs
dt3 <- dt2[average_dose <= 1]
n3 <- nrow(dt3)

# Step 4 — Remove singletons
dt4 <- dt3[number_genomes != 1]
n4 <- nrow(dt4)
fwrite(dt4, paste0(outdir_dat, "/pre_filter_ags.csv"))

# Step 5 — Remove piS outliers
dt5 <- dt4[!gene_family %chin% piS_3IQR_outliers$gene_family]
n5 <- nrow(dt5)

# Step 6 — Remove non syntenic
dt6 <- dt5[gene_family %in% syntenic_ags]
n6 <- nrow(dt6)
fwrite(dt6[,.(gene_family)], paste0(outdir_dat, "/set_1_names.csv"))

step_summary <- data.table(
  step = c(
    "Total gene families",
    "After removing core",
    "After removing paralogs",
    "After removing singletons",
    "After removing piS outliers",
    "After removing non-syntenic"
  ),
  n_genes = c(n1, n1-n2, n2-n3, n3-n4, n4-n5, n5-n6),
  cum_ag = c(NA, n2, n3, n4, n5,n6)
)

step_summary

# save data 
fwrite(step_summary, paste0(outdir_dat, "/step_summary.csv"))
       

