# Calcualte the age of AGs

PIRATE_dir = "./input_data/PIRATE_1695_out/"

# Get pirate pangenome 
kpne_pangome_dt = fread(paste0(PIRATE_dir, "PIRATE.gene_families.ordered.tsv"))

# list informative columns
cols_to_keep = c("gene_family", "threshold", "number_genomes")

cols_to_rmv = colnames(kpne_pangome_dt[,1:22])[!colnames(kpne_pangome_dt[,1:22]) %in% cols_to_keep]

kpne_pangome_dt[, (cols_to_rmv) := NULL]

kpne_pangome_dt[, pan_freq := round(number_genomes/1695, digits=3)]

kpne_pangome_dt_lng <- melt(kpne_pangome_dt,
                            id.vars = c("gene_family", "threshold", "number_genomes", "pan_freq"), 
                            variable.name = "geno_id", 
                            value.name = "loci_id")

# remove missing loci
kpne_pangome_dt_lng <- kpne_pangome_dt_lng[loci_id!="",]

# 8321877 unique loci
length(unique(kpne_pangome_dt_lng[pan_freq >= 0.95, gene_family]))# missing 40?
length(unique(kpne_pangome_dt_lng[pan_freq < 0.95, gene_family]))
length(unique(kpne_pangome_dt_lng[number_genomes == 1, gene_family]))

# NEED TO CORRECT FOR THE FISSION/FUSION AND DOASGE ISSUE


# Create AG set without ORfans --------------------------------------------
ag_lng <- kpne_pangome_dt_lng[pan_freq < 0.95 & number_genomes != 1, ]


# Map to loci names (to then extract .fnas)

kpne_pan_genes = fread(paste0(PIRATE_dir,"loci_list.tab"))

colnames(kpne_pan_genes) = c("loci_id", "gene_family", "threshold", "allele", "geno_id")

ag_lng_mapped <- merge(kpne_pan_genes, ag_lng, by=c("gene_family", "loci_id", "geno_id", "threshold"))

# 11,992 unique AG gene families
fwrite(ag_lng_mapped, paste0(outdir_dat,"/ag_lng_mapped.csv"))

