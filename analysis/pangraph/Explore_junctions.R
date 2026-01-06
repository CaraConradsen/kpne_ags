

msu_regions_anchored <- fread(paste0(outdir_tab, "/msu_regions_anchored.csv"))

pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"))

# Examine the msu blocks ---------------------------------------------------

msu_loci <- msu_regions_anchored[grepl(":", anchor)]

# Are there any AGs across junctions?
msu_loci_grp = unique(msu_loci[,.(ag_type, gene_family, msu)])

across_msu_loci = msu_loci_grp[,.(n = .N), by = c("ag_type", "gene_family")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_msu = 0
msu_regions_anchored[gene_family %chin% across_msu_loci, acrs_msu := 1]

# Examine the syntenic junctions ---------------------------------------------------

junction_loci <- msu_regions_anchored[grepl(":", anchor)]

# Are there any AGs across junctions?
junction_loci_grp = unique(junction_loci[,.(ag_type, gene_family, msu, anchor)])

across_junction_loci = junction_loci_grp[,.(n = .N), by = c("ag_type", "gene_family", "msu")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_jun = 0
msu_regions_anchored[gene_family %chin% across_junction_loci, acrs_jun := 1]

