
id_ags <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("gene_family","number_genomes", "average_dose", 
                                  "ag_type", "acrs_msu", "acrs_jun")))

dt <- copy(id_ags)

# pi outliers
piS_3IQR_outliers <- fread(paste0(outdir_dat, "/piS_3IQR_outliers.csv"))

# id set 2 outliers
set_1_ags <- fread(paste0(outdir_dat, "/phylo_info.csv"))[syn_jun == 1, gene_family]
set_2_ags <- fread(paste0(outdir_dat, "/phylo_info.csv"))[phy_par_keep == 1, gene_family]


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
dt6 <- dt5[gene_family %in% set_1_ags]
n6 <- nrow(dt6)
fwrite(dt6[,.(gene_family)], paste0(outdir_dat, "/set_1_names.csv"))


# Step 7 — Remove genes with multiple parsimonious gains
dt7 <- dt6[gene_family %in% set_2_ags]
n7 <- nrow(dt7)
fwrite(dt7[,.(gene_family)], paste0(outdir_dat, "/set_2_names.csv"))

step_summary <- data.table(
  step = c(
    "Total gene families",
    "After removing core",
    "After removing paralogs",
    "After removing singletons",
    "After removing piS outliers",
    "After removing non-syntenic",
    "After removing mutliple\n parsimonious gainst"
  ),
  n_genes = c(n1, n1-n2, n2-n3, n3-n4, n4-n5, n5-n6, n6-n7),
  cum_ag = c(NA, n2, n3, n4, n5,n6, n7)
)

step_summary

# save data 
fwrite(step_summary, paste0(outdir_dat, "/step_summary.csv"))
       

