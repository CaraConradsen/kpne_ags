# kpne genes
# in this file, all core + ag + and mge 
# regions will be put together in a single dt


# 1. Import gffs -------------------------------------------------------------
# get gff info
gff_files_list = list.files("./input_data/kpne_1695_gffs", 
                            pattern = ".gff", full.names = T)

# start timer
start.time <- Sys.time()

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Process assembly data
kpne_prokka_annotation <- foreach(i = 1:length(gff_files_list), 
                                .combine = rbind, 
                                .packages = c("data.table", "rtracklayer")) %dopar% { 
                                  # Try to import GFF file, skip if it fails
                                  gff_dt <- tryCatch({
                                    as.data.table(rtracklayer::import(gff_files_list[i], 
                                                                      colnames = c("source", "type", "ID", "locus_tag", "product", "Name")))
                                  }, error = function(e) {
                                    warning(paste("Skipping unreadable file:", gff_files_list[i]))
                                    return(NULL)  # Return NULL to skip this iteration
                                  })
                                  
                                  # Return processed data
                                  gff_dt
                                }




# Stop the parallel backend after the loop is done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 18.45449 mins

# get unique loci IDs
kpne_prokka_annotation$locus_id = paste0(gsub("SPARK_","K", kpne_prokka_annotation$seqnames),
                                         gsub("k__0","", kpne_prokka_annotation$locus_tag))

# There are 8835314 unique loci


#fwrite(kpne_prokka_annotation, paste0(outdir_dat, "/kpne_prokka_annotation.csv"))
# kpne_prokka_annotation <- fread(paste0(outdir_dat, "/kpne_prokka_annotation.csv"))


# 2. Get Gene family contig names from alns -------------------------------

# set directories
pirate_genes_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_1695_out/feature_sequences"

# get gene names
pirate_genes_faas <- list.files(path = pirate_genes_dir, pattern = ".aa.fasta", full.names = TRUE)

# start timer
start.time <- Sys.time()
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# foreach loop in parallel
pirate_pangenome_aln<- foreach(i = 1:length(pirate_genes_faas),
                               .combine = "rbind",
                               .packages = c("Biostrings")) %dopar% {
                                 gene_aln_temp <- readAAStringSet(pirate_genes_faas[i])
                                 gene_aln_temp <- lapply(names(gene_aln_temp),
                                                         function(x) strsplit(x, "\\|")[[1]])
                                 gene_aln_temp <- data.frame(gene_family_aln = gsub(".aa.fasta","", basename(pirate_genes_faas[i])),
                                                             ID_contig = sapply(gene_aln_temp, `[`, 1))
                                 # info to be outputted
                                 gene_aln_temp
                               }

#fwrite(pirate_pangenome_aln, paste0(outdir_dat, "/pirate_pangenome_aln.csv"))
# pirate_pangenome_aln <- fread(paste0(outdir_dat, "/pirate_pangenome_aln.csv"))



#3. Import integrons (IntegronFinder)  -------------------------------------------

integron_dir = list.files("./input_data/integron_finder_results/", full.names = TRUE,
                          recursive = TRUE, pattern = ".summary")

# start timer
start.time <- Sys.time()

cl <- makeCluster(4)#num_cores)
registerDoParallel(cl)

# Process assembly data
kpne_integrons <- foreach(i = 1:length(integron_dir), 
                                  .combine = rbind, 
                          .packages = c("data.table", "stringr", "S4Vectors")) %dopar% { 
                            
                            chunk <- readLines(integron_dir[i])
                            chunk <- chunk[3:length(chunk)]
                            
                            # Split the character vector by tab
                            split_data <- strsplit(chunk, "\t")
                            
                            # Convert to data.table
                            dt <- data.table(matrix(unlist(split_data), ncol = length(split_data[[1]]), byrow = TRUE))
                            
                            tmp_dt <- dt[2:nrow(dt),]
                            
                            # Assign column names
                            setnames(tmp_dt, c(unlist(dt[1,])))
                            
                            integron_contigs = tmp_dt[CALIN!=0 | complete!=0, ID_replicon]
                            
                            if(isEmpty(integron_contigs)){
                              return(NULL)  # Return NULL to skip this iteration
                            } else{
                              chunk <- readLines(gsub("summary","integrons",integron_dir[i]))
                              chunk <- chunk[3:length(chunk)]
                              
                              # Split the character vector by tab
                              split_data <- strsplit(chunk, "\t")
                              
                              # Convert to data.table
                              dt <- data.table(matrix(unlist(split_data), ncol = length(split_data[[1]]), byrow = TRUE))
                              
                              tmp_dt <- dt[2:nrow(dt),]
                              
                              # Assign column names
                              setnames(tmp_dt, c(unlist(dt[1,])))
                              
                              tmp_dt <- tmp_dt[ID_replicon %chin% integron_contigs, .(ID_replicon, pos_beg, pos_end, strand, type, type_elt, annotation)]
                              
                              tmp_dt
                            }
                          }

# Stop the parallel backend after the loop is done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 6.855873 secs

# tidy up integrons to merge with main data
kpne_integrons$type = gsub("complete", "integron", kpne_integrons$type)

kpne_integrons[type_elt=="protein", type_elt:= "CDS"]

kpne_integrons$product = paste(kpne_integrons$type, kpne_integrons$annotation, sep=" ")

kpne_integrons[, strand := fcase(strand == -1, "-",
                                 strand == 1, "+",
                                 default = NA)]

kpne_integrons[,c("type", "annotation") := NULL]

colnames(kpne_integrons)[c(1:3,5)] = c("seqnames", "start", "end", "type")

kpne_integrons$source = "IntegronFinder"

#fwrite(kpne_integrons, paste0(outdir_dat, "/kpne_integrons.csv"))
# kpne_integrons <- fread(paste0(outdir_dat, "/kpne_integrons.csv"))

