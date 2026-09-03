# input require data sets
# ags

fasta_dir <- paste0(outdir_dat, "/ag_fus_alns/")

ag_fasta_files = list.files(fasta_dir, pattern = ".fasta",
                            full.names = TRUE, recursive = TRUE)

cl <- makeCluster(num_cores)
registerDoParallel(cl)

count_missing_seq <- foreach(
  i = ag_fasta_files,
  .packages = c("Biostrings", "data.table"),
  .combine  = "rbind"
) %dopar% {
  
  temp_string <- tryCatch(
    readDNAStringSet(i),
    error = function(e) return(NULL)
  )
  
  if (is.null(temp_string)) return(NULL)
  
  gene_family <- gsub("\\.fasta$", "", basename(i))
  
  counts <- sapply(temp_string, function(fna) {
    sum(letterFrequency(fna, letters = c("N", "-"), as.prob = TRUE))
  })
  
  data.table(
    gene_family = gene_family,
    number_genomes = length(temp_string),
    n_isolates_gte_0.05 = length(counts[counts >= 0.05]),
    mean_missing = mean(counts)
  )
}

stopCluster(cl)

count_missing_seq[, gene_lvl_prevl := n_isolates_gte_0.05/number_genomes]

# fwrite(count_missing_seq, paste0(outdir_dat,"/count_missing_seq.csv"))


# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))[,gene_family]
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))[,gene_family]

count_missing_seq[, set := fcase(gene_family %in% set_1_ags, 1,
                                 default = 0)]
count_missing_seq[, set := fcase(gene_family %in% set_2_ags, 2,
                                 default = set)]

par(mforw = c(2,1))

with(count_missing_seq[set %in% c(1,2)],
     hist(mean_missing, breaks = 100, 
          main = "Set 1"))



