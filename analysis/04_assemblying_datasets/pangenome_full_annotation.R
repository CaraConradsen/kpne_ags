# This script integrates pangenome data
# analysis to refine the classification of core and accessory genes.
# 
# Functions:
# - Parse PIRATE output to import gene presence/absence data.
# - take note of the fusion/fisson genes and use consensus
# - add Eggnog
# - Merge GFF annotations and external metadata for downstream comparative genomics.
# 
# Outputs:
# - Updated core/accessory gene tables.
# - Integrated GFF annotation dataset for visualization or recombination analysis.


# PART 1. Import gene annotations -------------------------------------------------
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
                        
                        anno_dt$gene_length = width(anno)
                        
                        # get cogs
                        res <- lapply(anno_dt$Dbxref, extract_cogs)
                        res_dt <- data.table::transpose(res)
                        anno_dt[, c("COG_funct_cat", "COG_ortho_grp") := res_dt]
                        anno_dt[, KEGG := extract_keggs(Dbxref), by = .I]
                        
                        anno_dt[type == "CDS", .(geno_id, seqnames, start,end, strand,gene_length, 
                                                            locus_tag,gene, product,
                                                            COG_funct_cat, COG_ortho_grp, KEGG)]
                        

                      }

stopCluster(cl)

# add pirate & core syntenic groups -----------------------------------------
focal_genome_names = unique(pirate_anno$geno_id)
pirate_lng <- fread("./input_data/PIRATE_260_hybrid_chr_out/PIRATE.gene_families.ordered.tsv",
                select = c("gene_family", "consensus_gene_name",
                           "consensus_product", "number_genomes", 
                           "cluster", "cluster_order", focal_genome_names))

pirate_lng <- melt(pirate_lng, 
                   id.vars = c("gene_family", "number_genomes",
                               "consensus_gene_name",
                               "consensus_product",
                               "cluster", "cluster_order"),
                   variable.name = "geno_id", 
                   value.name = "fus_locus_tag"
                   )

pirate_lng <- pirate_lng[fus_locus_tag!=""]

pirate_lng$locus_tag = pirate_lng$fus_locus_tag

# tidy fusion brackets
pirate_lng$locus_tag = gsub("[()]", "", pirate_lng$locus_tag)

#tidy fission pirate loci
pirate_lng <- pirate_lng[
  ,
  .(locus_tag = unlist(strsplit(locus_tag, ";"))),
  by = setdiff(names(pirate_lng), "locus_tag")
]

#tidy fusion pirate loci
pirate_lng <- pirate_lng[
  ,
  .(locus_tag = unlist(strsplit(locus_tag, ":"))),
  by = setdiff(names(pirate_lng), "locus_tag")
]

# add in annotation, removing all short genes not analysed in pirate
pirate_anno <- merge(pirate_lng, pirate_anno, 
                     all.x = TRUE, by = c("geno_id","locus_tag"))


# fix fission/fusion start
collapsed <- pirate_anno[, .(fus_locus_tag, start, end, strand, gene_length)]

collapsed <- collapsed[, .(
  fstart  = min(start),
  fend    = max(end),
  fstrand = strand[which.max(gene_length)]
), by = fus_locus_tag]


pirate_anno <- merge(pirate_anno, collapsed,
                     all.x = TRUE, by = c("fus_locus_tag"))


# Add ST ------------------------------------------------------------------
ST_info <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv", 
                 select = c("id", "ST"))

colnames(ST_info)[1] = "geno_id"

pirate_anno <- merge(pirate_anno, ST_info,
                     all.x = TRUE, by = "geno_id")


# Add eggnog --------------------------------------------------------------
eggnog24 <- fread("./input_data/cog-24.mapping.tab",
                  select = c("COG no.", "Func"))

colnames(eggnog24) = c("COGid","COG_category")

eggnog <- fread("./input_data/eggnog/klebsiella_eggnog.emapper.annotations", 
                skip = 4, fill = TRUE, 
                select = c("#query", "eggNOG_OGs", "score"))

# Pick best hit per protein (highest bitscore)
eggnog <- eggnog[order(-score), .SD[1], by = `#query`]

# extract first COG ID
eggnog[, COGid := sub(".*?(COG[0-9]+).*", "\\1", eggNOG_OGs)]

# replace non-COG rows with NA
eggnog[!grepl("COG[0-9]+", eggNOG_OGs), COGid := NA]

# map eggnog24 onto our focal loci
eggnog <- merge(eggnog[,.(`#query`, COGid)], eggnog24,
                all.x = TRUE, by = "COGid")

colnames(eggnog)[2:3] <- c("gene_family", "egng_cog")

pirate_anno <- merge(pirate_anno, eggnog[,.(gene_family, egng_cog )],
                     all.x = TRUE, by = "gene_family")


# PART 2. Add in MGE information --------------------------------------------------

# phastest

phastest_files <- list.files("./input_data/phastest_results/",
                             pattern = ".json", full.names = TRUE,
                             recursive = TRUE)

phastest_dt <- foreach(pf = phastest_files,
                    .packages = c("jsonlite"),
                    .combine = "rbind") %do% {
                      
                      geno_id = gsub(".json", "", basename(pf))
                      
                      pf_df <- fromJSON(pf)
                      
                      if (!is.data.frame(pf_df)) {
                        return(NULL)  
                      }
                      
                      cbind(pf_df,geno_id)
                    }

setDT(phastest_dt)

phastest_dt[,region := NULL]

colnames(phastest_dt)[3] = "phg_state"

phastest_dt <- phastest_dt[, .(
  geno_id,
  rstart = start,
  rend   = stop,
  phg_state,
  most_common_phage
)]

setkey(phastest_dt, geno_id, rstart, rend)

genes_dt <- pirate_anno[, .(
  gene_family,
  geno_id,
  locus_tag,
  gstart = start,
  gend   = end
)]

setkey(genes_dt, geno_id, gstart, gend)

gene_phage <- foverlaps(
  genes_dt,
  phastest_dt,
  by.x = c("geno_id", "gstart", "gend"),
  by.y = c("geno_id", "rstart", "rend"),
  type = "within",     # gene fully inside phage region
  nomatch = NA
)


pirate_anno <- merge(pirate_anno, 
                     gene_phage[!is.na(phg_state), .(locus_tag, phg_state, most_common_phage)],
                     all.x = TRUE, by = "locus_tag")

# intermediate 
# fwrite(pirate_anno, paste0(outdir_dat, "/pirate_anno.csv"))
# pirate_anno <- fread(paste0(outdir_dat, "/pirate_anno.csv"))


# Add integrons -----------------------------------------------------------

# determine if any integrons were found
integron_files <- list.files("./input_data/integron_finder_results/", 
                           recursive = TRUE, pattern = ".summary",
                           full.names = TRUE)

integrons_sum <- foreach(i = integron_files,
                     .packages = "data.table",
                     .combine = "rbind") %do%{
                       
                       geno_id = gsub(".summary", "", basename(i))
                       
                       cbind(geno_id,
                             fread(i,
                                   select = c("ID_replicon", "CALIN", "complete","In0")
                                   )
                             )
                     } 

# filter results
integrons_sum <- integrons_sum[CALIN!=0 | complete != 0 | In0 != 0]

# long results
integrons_sum <- melt(integrons_sum, id.vars=1:2,
                      variable.name = "integron")[value!=0]


integrons_sum <- foreach(i = integrons_sum$geno_id,
                         .packages = "data.table",
                         .combine = "rbind") %do%{
                           
                           int_dat <- fread(list.files(paste0("./input_data/integron_finder_results/",
                                                              i,"/"),
                                                       recursive = TRUE, pattern = ".integron",
                                                       full.names = TRUE),
                                            select = c("type", "pos_beg", "pos_end", "strand", "annotation"))
                           
                           colnames(int_dat)[c(1:3,5)] <- c("integron", "start", "end", "integron_anno")
                           int_dat$geno_id  = i
                           
                           int_dat
                         } 


# fix strand
integrons_sum[, strand:= fcase(strand == "-1", "-",
                               default = "+")]

integrons_dt <- integrons_sum[, .(
  geno_id,
  rstart = start,
  rend   = end,
  integron
)]

setkey(integrons_dt, geno_id, rstart, rend)
setkey(genes_dt, geno_id, gstart, gend)

gene_integrons <- foverlaps(
  genes_dt,
  integrons_dt,
  by.x = c("geno_id", "gstart", "gend"),
  by.y = c("geno_id", "rstart", "rend"),
  type = "any",     # gene fully inside phage region
  nomatch = NA
)

pirate_anno <- merge(pirate_anno, 
                     gene_integrons[!is.na(integron), .(locus_tag, integron)],
                     all.x = TRUE, by = "locus_tag")


# ICEfinder

ice_files <- list.files("./input_data/icefinder_results/", 
                        recursive = TRUE, full.names = TRUE,
                        pattern = "ICE_details.tsv")

ice_dt <- foreach(i = ice_files, 
                  .packages = c("data.table"),
                  .combine = "rbind") %do% {
                    fread(i, 
                          select = c("MGE", "Location"))
                  } 

# Split on last underscore
ice_dt[, c("geno_id", "ice") := tstrsplit(MGE, "_(?=[^_]+$)", perl = TRUE)]

# Split out the numbers inside parentheses
ice_dt[, c("start", "end") := tstrsplit(gsub(".*:\\((\\d+)\\.\\.(\\d+)\\)", "\\1_\\2", Location), "_")]

# Convert to integer
ice_dt[, c("start", "end") := .(as.integer(start), as.integer(end))]


ice_dt <- ice_dt[, .(
  geno_id,
  rstart = start,
  rend   = end,
  ice
)]

setkey(ice_dt, geno_id, rstart, rend)
setkey(genes_dt, geno_id, gstart, gend)

gene_ice <- foverlaps(
  genes_dt,
  ice_dt,
  by.x = c("geno_id", "gstart", "gend"),
  by.y = c("geno_id", "rstart", "rend"),
  type = "any",     # gene fully inside phage region
  nomatch = NA
)

pirate_anno <- merge(pirate_anno, 
                     gene_ice[!is.na(ice), .(locus_tag, ice)],
                     all.x = TRUE, by = "locus_tag")

# merge with pan_anno

fwrite(pirate_anno, paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))
# pirate_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))


# ISEScan  ----------------------------------------------------------------

# TnCentral ---------------------------------------------------------------



# Pangenome for analysis --------------------------------------------------
# need to account for fusion loci on different strands

pan_anno <- unique(pirate_anno[,.(gene_family,geno_id,fus_locus_tag, 
                                  number_genomes, consensus_product,
                                  fstart, fend, fstrand, ST)])

colnames(pan_anno)[6:8] <- c("start","end","strand")


# Assign Core vs accessory ------------------------------------------------

n_genomes = length(unique(pan_anno[, geno_id]))


pan_anno[, ag_type := fcase(number_genomes < 0.99 * n_genomes & number_genomes >=  0.95 * n_genomes, "soft",
                            number_genomes < 0.95 * n_genomes & number_genomes >=  0.15 * n_genomes, "shell",
                            number_genomes < 0.15 * n_genomes, "cloud",
                            default = "core")]
# get total length
fasta_dir <- list.files("./input_data/kpne_260_genohead_fasta/",
                         recursive = TRUE, full.names = TRUE,
                        pattern = ".fasta")

tot_length_dt <- foreach(i = fasta_dir,
                         .packages = c("Biostrings", "data.table"),
                         .combine = "rbind") %do% {
                           dna_str <- readDNAStringSet(i)
                           data.frame(geno_id = names(dna_str), tot_len = width(dna_str))
                         } 
setDT(tot_length_dt)

pan_anno <- merge(pan_anno, tot_length_dt,
                  all.x = TRUE, by =c("geno_id"))

fwrite(pan_anno, paste0(outdir_dat, "/pangenome_anno.csv"))






