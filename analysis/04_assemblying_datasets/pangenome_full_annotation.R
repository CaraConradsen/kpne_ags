# This script integrates pangenome data
# analysis to refine the classification of core and accessory genes.
# 
# Functions:
# - Parse PIRATE output to import gene presence/absence data.
# - Merge GFF annotations and external metadata for downstream comparative genomics.
# 
# Outputs:
# - Updated core/accessory gene tables.
# - Integrated GFF annotation dataset for visualization or recombination analysis.


# Import gene annotations -------------------------------------------------
hybrid_gffs_files <- list.files("./input_data/PIRATE_260_hybrid_chr_out/modified_gffs/",
                                full.names = TRUE, recursive = TRUE, 
                                pattern = ".gff")

# Extract COGS 
extract_cogs <- function(db_ls){
  db_ls <- unlist(db_ls)
  cogs <- gsub("COG:","", db_ls[grepl("COG:", db_ls)])
  cogs <- cogs[order(nchar(cogs))]
  return(cogs)
}

extract_keggs <- function(db_ls){
  db_ls <- unlist(db_ls)
  kegg <- gsub("KEGG:","", db_ls[grepl("KEGG:", db_ls)])
  return(kegg)
}

cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

pirate_anno <- foreach(i = 1:length(hybrid_gffs_files), 
                      .packages = c("rtracklayer", "gUtils"),
                      .combine = "rbind") %dopar% {
                        
                        geno_id = gsub(".gff","", basename(hybrid_gffs_files[i]))
                        
                        # get pirate anno
                        anno <- rtracklayer::import.gff(hybrid_gffs_files[i])
                        
                        anno_dt <- gr2dt(anno)
                        
                        anno_dt$geno_id = geno_id
                        
                        # get cogs
                        res <- lapply(anno_dt$Dbxref, extract_cogs)
                        res_dt <- data.table::transpose(res)
                        anno_dt[, c("COG_funct_cat", "COG_ortho_grp") := res_dt]
                        anno_dt[, KEGG := extract_keggs(Dbxref), by = .I]
                        
                        anno_dt[type == "CDS", .(geno_id, seqnames, start,end, strand, 
                                                            locus_tag,gene, product,
                                                            COG_funct_cat, COG_ortho_grp, KEGG)]
                        

                      }

stopCluster(cl)


# add pirate & core syntenic groups -----------------------------------------
focal_genome_names = unique(pirate_anno$geno_id)
pirate <- fread("./input_data/PIRATE_260_hybrid_chr_out/PIRATE.gene_families.ordered.tsv",
                select = c("gene_family", "consensus_gene_name",
                           "consensus_product", "number_genomes", 
                           "cluster", "cluster_order", focal_genome_names))

pirate_lng <- melt(pirate, 
                   id.vars = c("gene_family", "number_genomes",
                               "consensus_gene_name",
                               "consensus_product",
                               "cluster", "cluster_order"),
                   variable.name = "geno_id", 
                   value.name = "fus_locus_tag"
                   )

pirate_lng <- pirate_lng[fus_locus_tag!=""]

pirate_lng[, locus_tag := fus_locus_tag]

#tidy fused pirate loci
pirate_lng[
  ,
  .(locus_tag = unlist(strsplit(locus_tag, ";"))),
  by = setdiff(names(pirate_lng), "locus_tag")
]


pirate_anno <- merge(pirate_anno, pirate_lng,
                     all.x = TRUE, by = c("geno_id","locus_tag"))


# Add ST ------------------------------------------------------------------
ST_info <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv", 
                 select = c("id", "ST"))

colnames(ST_info)[1] = "geno_id"

pirate_anno <- merge(pirate_anno, ST_info,
                     all.x = TRUE, by = "geno_id")

# fwrite(pirate_anno, paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))
# pirate_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))


# Add in MGE information --------------------------------------------------
#ISEScan 




