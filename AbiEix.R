pangraph_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("gene_family", "number_genomes", "gene_length","egng_cog", 
                                         "consensus_gene_name", "consensus_product"))

unique(pangraph_anno[grepl("AbiE", consensus_product)]) 

AbiEix_ags <- c("g001254","g007857","g003299","g011206","g014364")

ag_pvals <- fread("./output/data/gene_pvals_fdr.csv")

# AbiEii
focal_ags <- ag_pvals[gene_family %chin% AbiEix_ags & freq >32, gene_family]

number_of_alleles <- fread("./input_data/PIRATE_260_hybrid_chr_hps0.6_out/PIRATE.gene_families.ordered.tsv")

number_of_alleles[gene_family %chin% focal_ags, .(gene_family, threshold, alleles_at_maximum_threshold)]

pangraph_anno[gene_family %chin% focal_ags][,.(n_alleles = .N), c("gene_family", "gene_length")]
