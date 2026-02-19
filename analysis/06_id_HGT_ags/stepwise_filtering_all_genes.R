
id_ags <- unique(fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                select = c("gene_family","number_genomes",
                           "ag_type", "acrs_msu", "acrs_jun")))


dt <- copy(id_ags)

# Step 1 — Total
dt1 <- dt
n1 <- nrow(dt1)

# Step 2 — Remove core
dt2 <- dt1[ag_type != "core"]
n2 <- nrow(dt2)

# Step 3 — Remove paralogs
dt3 <- dt2[!grepl("_", gene_family)]
n3 <- nrow(dt3)

# Step 4 — Remove singletons
dt4 <- dt3[number_genomes != 1]
n4 <- nrow(dt4)

# Step 5 — Remove non syntenic
dt5 <- dt4[acrs_msu!=1][acrs_jun!=1]
n5 <- nrow(dt5)

step_summary <- data.table(
  step = c(
    "Total gene families",
    "After removing core",
    "After removing paralogs",
    "After removing singletons",
    "After removing non-syntenic"
  ),
  n_genes = c(n1, n1-n2, n2-n3, n3-n4, n4-n5),
  cum_ag = c(NA, n2, n3, n4, n5)
)


# save data 
fwrite(dt5[, .(gene_family)], paste0(outdir_dat, "/group_syn_gene_fams.csv"))
       
# save data 
consistencyindex <- fread(paste0(outdir_dat, "/consistencyindex.csv"))


