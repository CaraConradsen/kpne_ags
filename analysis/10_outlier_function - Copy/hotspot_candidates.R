

pirate_anno_mges <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                     select = c("geno_id", "gene_family", "number_genomes", "consensus_gene_name", 
                                "consensus_product", "alleles_at_maximum_threshold", 
                                "phg_state", "most_common_phage", 
                                "integron", "ice", "is_fam", "is_cluster", "msu", "anchor"))


set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

pirate_anno_mges <- unique(pirate_anno_mges[gene_family %chin% set_1_ags$gene_family])


fwrite(pirate_anno_mges, "C:/Users/carac/Desktop/pirate_anno_mges.csv")
