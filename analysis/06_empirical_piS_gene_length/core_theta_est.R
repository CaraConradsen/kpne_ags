# Get AG piS (age)

# input require data sets
# ags

core_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                   select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                              "number_genomes", "start","end", "average_dose",
                              "strand","ST","ag_type"))


# remove estimates of core, paralogs and singletons
core_anno <- core_anno[ag_type=="core"][average_dose <= 1]

# number of genomes in pangenome
tot_pangenome_size = length(unique(core_anno$geno_id))

list_unique_cores = unique(core_anno$gene_family)

# fission/fusion list -----------------------------------------------------

core_anno[, n:=.N, by = c("gene_family", "geno_id")]

# Filter rows where n > 1
multi_loci <- core_anno[n > 1]

# Group by gene_family and geno_id, collect locus_tags
f_loci <- multi_loci[, .(locus_tags = list(locus_tag)), by = .(gene_family, geno_id)]

# Now collapse by gene_family into a list of lists
f_loci <- f_loci[, .(loci_per_geno = list(locus_tags)), by = gene_family]

# Convert to named list
gene_family_list <- setNames(f_loci$loci_per_geno, f_loci$gene_family)

# Set hypergeometric sampling limits --------------------------------------
m_target = 372

set.seed(42)

downsample_segsites <- function(S, m, m_target) {
  if (m < m_target) {
    # Gene is below target length — no downsampling needed
    return(NA)
  } else if (m == m_target){
    return(S)
  }
  # Draw one realisation from the hypergeometric
  rhyper(nn = 1,
         m  = S,
         n  = m - S,
         k  = m_target)
}

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

# set up directories

gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_out/feature_sequences/"

file_dir = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data"
py_script = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/PolyFastA.py"

# Step 1: write all FASTAs first 
fasta_dir <- paste0(outdir_dat, "/polyfasta_input/")
dir.create(fasta_dir, showWarnings = FALSE)

for (i in list_unique_cores) {
  
  temp_string <- tryCatch(
    readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
    error = function(e) NULL
  )
  if (is.null(temp_string)) next
  
  sampled_loci <- core_anno[gene_family == i, locus_tag]
  temp_string  <- temp_string[names(temp_string) %in% sampled_loci]
  
  if (length(gene_family_list[[i]]) >= 1) {
    normies <- temp_string[!names(temp_string) %in% unlist(gene_family_list[[i]])]
    collapsed_loc <- lapply(gene_family_list[[i]], function(fus_loc) {
      collapse_alignment(temp_string[names(temp_string) %in% fus_loc])
    })
    collapsed_loc     <- DNAStringSet(unlist(collapsed_loc))
    collapsed_loc_tag <- sapply(gene_family_list[[i]], `[`, 1)
    collapsed_loc     <- DNAStringSet(setNames(collapsed_loc, collapsed_loc_tag))
    temp_string       <- c(normies, DNAStringSet(collapsed_loc))
  }
  
  if (length(temp_string) == 0) next
  
  writeXStringSet(temp_string,
                  paste0(fasta_dir, i, ".fasta"),
                  format = "fasta")
}

# Step 2: single PolyFastA call on the whole directory
fasta_dir_wsl <- gsub("C:", "/mnt/c", gsub("\\\\", "/", fasta_dir))
result_wsl    <- gsub("C:", "/mnt/c", gsub("\\\\", "/", 
                                           paste0(outdir_dat, "/polyfasta_core_results.tsv")))

cmd_polyfasta <- paste(
  "wsl",
  "/home/carac/anaconda3/bin/python",
  "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/PolyFastA.py",
  "--dir", fasta_dir_wsl,
  "--cds",
  "--out", result_wsl
)

system(cmd_polyfasta, wait = TRUE, 
       ignore.stdout = TRUE, 
       ignore.stderr = TRUE)

#Time difference of 51.88084 secs

# Step 3: read results
polyfasta_results <- fread(paste0(outdir_dat, "/polyfasta_core_results.tsv"),
                           select = c("file", "seg_sites_S"))


# the file column will be the safe_name.fasta — recover gene family name
polyfasta_results[, gene_family := gsub("\\.fasta$", "", basename(file))]

# Step 4: now do the rest of your loop without any system() calls

# start timer
start.time <- Sys.time()

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# foreach loop
core_2Nu_dt <- foreach(i = list_unique_cores, 
                       .combine = rbind, 
                       .packages = c("Biostrings","ape","pegas","data.table")) %do% {
                         
                         S_syn <- polyfasta_results[gene_family == i, seg_sites_S]
                         hyperg_S_syn <- if (length(S_syn) == 0) NA else S_syn
                         
                         # Read the fixed FASTA, skip if missing
                         temp_string <- tryCatch(
                           readDNAStringSet(paste0(fasta_dir, i, ".fasta")),
                           error = function(e) return(NULL)
                         )
                         
                         # skip if file is missing
                         if (is.null(temp_string)) return(NULL)
                         
                         # Convert to DNAbin
                         dna_bin <- as.DNAbin(temp_string)
                         
                         # Calculate segregating sites 
                         seg_site <- seg.sites(dna_bin)
                         S <- length(seg_site)
                         
                         # Calculate θw
                         # save k strains
                         n <- length(temp_string)
                         
                         L <- width(temp_string)[1]
                         
                         # If downsampling 
                         if(width(temp_string)[1] > m_target){
                           hyperg_S_syn = downsample_segsites(hyperg_S_syn, L, m_target)
                           L = m_target
                         }
                         
                         a_n <- sum(1 / (1:(n-1)))
                         
                         theta_w <- hyperg_S_syn / (L * a_n)
                         
                         # calculate pi here
                         
                         
                         #return value
                         data.frame(
                           gene_family = i,
                           m = width(temp_string)[1],
                           freq = n, 
                           S = S,
                           S_syn = S_syn,
                           hyperg_S_syn = hyperg_S_syn,
                           hyperg_m = L,
                           theta_w = theta_w
                           # pi = pi_value,
                         )
                         
                       }

# Stop cluster when done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 3.97464 mins

setDT(core_2Nu_dt)

# Save output
fwrite(core_2Nu_dt, paste0(outdir_dat, "/core_2Nu_dt.csv"))


#output 2Nu
# mean_pi <- mean(core_2Nu_dt$pi, na.rm = TRUE)
mean_theta_w <- mean(core_2Nu_dt$theta_w, na.rm = TRUE)
median_m <- median(core_2Nu_dt$m, na.rm = TRUE)

mut_rate = data.frame(cbind(ag_type = "core", mean_theta_w, median_m)) #mean_pi,

fwrite(mut_rate, paste0(outdir_dat, "/core_mut_rate.csv"))


# remove additional fasta files
unlink(fasta_dir, recursive = TRUE)





# Inspect values ----------------------------------------------------------

png(paste0(outdir_fig,"/2Nmu_estimators.png"),
    width = 13, height = 15, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

# Plot comparison
with(core_2Nu_dt,
     plot(pi, theta_w,
          pch = 16, col = rgb(0.1,0.2,0.4, alpha = 0.5),
          cex = 0.5, 
          xlab = "Nucleotide diversity (π)", 
          ylab = "Watterson's θw",
          main = "Comparison of 2Nμ estimates")
)
abline(0, 1, col = "red", lty = 2)
dev.off()



