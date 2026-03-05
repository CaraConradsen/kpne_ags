# Get AG piS (age)

# input require data sets
# ags

pangraph_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                                  "number_genomes", "start","end",
                                  "strand","ST","ag_type"))


# estimates for set 1
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

pangraph_anno <- pangraph_anno[gene_family %chin% set_1_ags$gene_family]

# number of genomes in pangenome
tot_pangenome_size = length(unique(pangraph_anno$geno_id))

list_unique_ags = unique(pangraph_anno$gene_family)

# fission/fusion list -----------------------------------------------------

pangraph_anno[, n:=.N, by = c("gene_family", "geno_id")]

# Filter rows where n > 1
multi_loci <- pangraph_anno[n > 1]

# Group by gene_family and geno_id, collect locus_tags
f_loci <- multi_loci[, .(locus_tags = list(locus_tag)), by = .(gene_family, geno_id)]

# Now collapse by gene_family into a list of lists
f_loci <- f_loci[, .(loci_per_geno = list(locus_tags)), by = gene_family]

# Convert to named list
gene_family_list <- setNames(f_loci$loci_per_geno, f_loci$gene_family)


# Alignment loop ----------------------------------------------------------

# function to collapse fission/fusion loci
collapse_alignment <- function(dna) {
  
  mat <- as.matrix(dna)
  
  merged <- apply(mat, 2, function(col) {
    bases <- col[col != "-"]
    if (length(bases) == 0) "-" else bases[1]
  })
  
  DNAString(paste0(merged, collapse = ""))
}


gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_out/feature_sequences/"

# start timer
start.time <- Sys.time()

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallel foreach loop
ag_age_S_dt <- foreach(i = list_unique_ags, 
                        .combine = rbind, 
                        .packages = c("Biostrings","ape","pegas","data.table")) %dopar% {
                          
                          # Read the FASTA, skip if missing
                          temp_string <- tryCatch(
                            readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
                            error = function(e) return(NULL)
                          )
                          
                          # skip if file is missing
                          if (is.null(temp_string)) return(NULL)
                          
                          #  Subset to loci present in pangraph annotation
                          sampled_loci = pangraph_anno[gene_family == i, locus_tag]
                          
                          # Subset by sequence names
                          temp_string <- temp_string[names(temp_string) %in% sampled_loci]
                          
                          # fix fission / fussion
                          if(length(gene_family_list[[i]]) >= 1){
                            normies <- temp_string[!names(temp_string) %in% unlist(gene_family_list[[i]])]
                            
                            collapsed_loc <- lapply(gene_family_list[[i]], function(fus_loc){
                              collapse_alignment(temp_string[names(temp_string) %in% fus_loc])
                            })
                            
                            
                            collapsed_loc <- DNAStringSet(unlist(collapsed_loc))
                            
                            collapsed_loc_tag <- sapply(gene_family_list[[i]], `[`, 1)
                            
                            collapsed_loc <- DNAStringSet(setNames(collapsed_loc, collapsed_loc_tag))
                            
                            normies <- c(normies, DNAStringSet(collapsed_loc))
                            
                            temp_string <- copy(normies)
                          }
                          
                          
                          # Skip if no sequence, ORFan or core gene
                          if (length(temp_string) == 0) return(NULL)
                          if (length(temp_string) %in% c(1, tot_pangenome_size)) return(NULL)
                          
                          # save k strains
                          loc_freq = length(temp_string)
                          
                          # Convert to "alignment" class
                          alignment_obj <- list(
                            nb = length(temp_string),
                            nam = names(temp_string),
                            seq = as.character(temp_string),  # must be character vector
                            com = NULL
                          )
                          class(alignment_obj) <- "alignment"
                          
                          kaks_result <- try(seqinr::kaks(alignment_obj), silent = TRUE)
                          
                          if (inherits(kaks_result, "try-error") || is.null(kaks_result$ks)) {
                            mean_ks <- NA
                          } else {
                            mean_ks <- mean(kaks_result$ks, na.rm = TRUE)
                          }
                          
                          # Compute gene length
                          gene_len = width(temp_string)[1]

                          # Convert to DNAbin
                          dna_bin <- as.DNAbin(temp_string)
                          
                          # Calculate segregating sites 
                          seg_site <- seg.sites(dna_bin)
                          S <- length(seg_site)
                          
                          #return value
                          data.frame(
                            gene_family = i,
                            m = gene_len,
                            mean_ks = mean_ks,
                            freq = loc_freq,
                            S = S
                          )

                        }


# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 24.2675 secs

# Stop cluster when done
stopCluster(cl)

setDT(ag_age_S_dt)

# Calculate segregating site per site
# m is gene length, S is number of segregating sites
ag_age_S_dt[, Sm := S / m]

# # checks
# ag_age_S_dt <- merge(ag_age_S_dt, unique(pangraph_anno[,.(gene_family, number_genomes)]), all.x = TRUE, by = "gene_family")
# ag_age_S_dt[freq != number_genomes]# should be nothing

fwrite(ag_age_S_dt, paste0(outdir_dat, "/ag_age_S_dt.csv"))
# ag_age_S_dt <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))



