# This script integrates pangenome data, plasmid identification, and genome synteny
# analysis to refine the classification of core and accessory genes.
# 
# Functions:
# - Parse PIRATE output to import gene presence/absence data.
# - Identify and annotate plasmid-borne contigs from assembly metadata.
# - Recalculate gene categorization (core vs accessory) based on plasmid content and synteny.
# - Detect conserved syntenic blocks across genomes.
# - Merge GFF annotations and external metadata for downstream comparative genomics.
# 
# Outputs:
# - Updated core/accessory gene tables.
# - Annotated synteny map and plasmid summary files.
# - Integrated GFF annotation dataset for visualization or recombination analysis.


# Id plasmid contigs ------------------------------------------------------


lng_rd_fasta_files = list.files("./input_data/kpne_485_fasta/",
                                full.names = TRUE)

cl <- makePSOCKcluster(num_cores-6)
registerDoParallel(cl)

asmbly_dat <- foreach(i = 1:length(lng_rd_fasta_files), 
                      .packages = c("mlplasmids", "Biostrings"),
                      .combine = "rbind") %dopar% {
                        ls = lng_rd_fasta_files[i]
                        lr_fasta = Biostrings::readDNAStringSet(ls)
                        
                        # second line of evidence classifying plasmids
                        plas_prob <- plasmid_classification(path_input_file = ls, 
                                                            species = 'Klebsiella pneumoniae', 
                                                            full_output = TRUE)[, c("Contig_name", "Prob_Chromosome")]
                        
                        # tidy info
                        colnames(plas_prob) = c("seqnames", "prob_chr")
                        plas_prob$prob_chr = round(plas_prob$prob_chr, digits = 2)
                        
                        # get largest contig
                        if(length(lr_fasta)>1){
                          # assign chromosome or plasmid 
                          chromosomal_contig = names(lr_fasta[width(lr_fasta) == max(width(lr_fasta))])
                          plasmid_contig <- names(lr_fasta[width(lr_fasta) < max(width(lr_fasta))])
                          contig_dat = rbind(data.frame(seqnames = chromosomal_contig, asmbly_type = "chromosome"),
                                             data.frame(seqnames = plasmid_contig, asmbly_type = "plasmid"))
                        }else{
                          contig_dat = data.frame(seqnames = names(lr_fasta), asmbly_type = "chromosome")
                        }
                        
                        merge(contig_dat, plas_prob, by = "seqnames")
                      }

stopCluster(cl)

setDT(asmbly_dat)
# fwrite(asmbly_dat, paste0(outdir_dat, "/asmbly_dat.csv"))
# asmbly_dat <- fread(paste0(outdir_dat, "/asmbly_dat.csv"))

# Import gene annotations -------------------------------------------------
lng_rd_gffs_files <- list.files("./input_data/PIRATE_485_lng_rds_out/modified_gffs/",
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

pirate_anno <- foreach(i = 1:length(lng_rd_gffs_files), 
                      .packages = c("rtracklayer", "gUtils"),
                      .combine = "rbind") %dopar% {
                        
                        geno_id = gsub(".gff","", basename(lng_rd_gffs_files[i]))
                        
                        # get pirate anno
                        anno <- rtracklayer::import.gff(lng_rd_gffs_files[i])
                        
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

#add assembly info

pirate_anno <- merge(pirate_anno, asmbly_dat[,.(seqnames, asmbly_type)],
                     all.x = TRUE, by = "seqnames")



# add pirate & core syntenic groups -----------------------------------------
focal_genome_names = unique(pirate_anno$geno_id)
pirate <- fread("./input_data/PIRATE_485_lng_rds_out/PIRATE.gene_families.ordered.tsv",
                select = c("gene_family", "number_genomes", 
                           "cluster", "cluster_order", focal_genome_names))

pirate_lng <- melt(pirate, 
                   id.vars = c("gene_family", "number_genomes", "cluster", "cluster_order"),
                   variable.name = "geno_id", 
                   value.name = "locus_tag"
                   )

pirate_lng <- pirate_lng[locus_tag!=""]

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




# Check for dnaA, dnaN ... ------------------------------------------------
setorderv(pirate_anno, c("geno_id", "seqnames", "start"))

correct_start <- pirate_anno[gene=="dnaA" & strand == "+" & grepl("_1$", seqnames) & start == 70, geno_id]
pirate_anno[geno_id %chin% correct_start & grepl("_1$", seqnames), .(end = max(end)), by = geno_id][,quantile(end)]

# short contigs
short_contig = pirate_anno[grepl("_1$", seqnames), .(end = max(end)), by = geno_id][end < 4200000, geno_id]
pirate_anno[geno_id %chin% short_contig & gene=="dnaA" & grepl("_1$", seqnames)]#

# remainder?
pirate_anno[!geno_id %chin% c(short_contig, correct_start) & gene=="dnaA"] 

pirate_anno[!ST %chin% pirate_anno[geno_id %chin%  correct_start,ST], ST] %>% unique()

