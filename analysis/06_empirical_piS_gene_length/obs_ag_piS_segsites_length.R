# Get AG piS (age)

# input require data sets
# ags

pangraph_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                                  "number_genomes", "start","end", "average_dose",
                                  "strand","ST","ag_type"))


# remove estimates of core, paralogs and singletons
pangraph_anno <- pangraph_anno[ag_type!="core"][number_genomes > 1][average_dose <= 1]

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

# set up directories

gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_out/feature_sequences/"

file_dir = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data"
py_script = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/PolyFastA.py"



# Step 1: write all FASTAs first 
fasta_dir <- paste0(outdir_dat, "/polyfasta_input/")
dir.create(fasta_dir, showWarnings = FALSE)

for (i in list_unique_ags) {
  
  temp_string <- tryCatch(
    readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
    error = function(e) NULL
  )
  if (is.null(temp_string)) next
  
  sampled_loci <- pangraph_anno[gene_family == i, locus_tag]
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
  if (length(temp_string) %in% c(1, tot_pangenome_size)) next
  
  writeXStringSet(temp_string,
                  paste0(fasta_dir, i, ".fasta"),
                  format = "fasta")
}

# Step 2: single PolyFastA call on the whole directory
fasta_dir_wsl <- gsub("C:", "/mnt/c", gsub("\\\\", "/", fasta_dir))
result_wsl    <- gsub("C:", "/mnt/c", gsub("\\\\", "/", 
                                           paste0(outdir_dat, "/polyfasta_results.tsv")))

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
polyfasta_results <- fread(paste0(outdir_dat, "/polyfasta_results.tsv"),
                           select = c("file", "seg_sites_S", "pi_S"))


# the file column will be the safe_name.fasta — recover gene family name
polyfasta_results[, gene_family := gsub("\\.fasta$", "", basename(file))]

# Step 4: now do the rest of your loop without any system() calls

# start timer
start.time <- Sys.time()

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# foreach loop
ag_age_S_dt <- foreach(i = list_unique_ags, 
                        .combine = rbind, 
                        .packages = c("Biostrings","ape","pegas","data.table")) %do% {
                          
                          S_syn <- polyfasta_results[gene_family == i, seg_sites_S]
                          S_syn <- if (length(S_syn) == 0) NA else S_syn
                          
                          pi_S <- polyfasta_results[gene_family == i, pi_S]
                          pi_S <- if (length(pi_S) == 0) NA else pi_S
                          
                          # Read the fixed FASTA, skip if missing
                          temp_string <- tryCatch(
                            readDNAStringSet(paste0(fasta_dir, i, ".fasta")),
                            error = function(e) return(NULL)
                          )
                          
                          # skip if file is missing
                          if (is.null(temp_string)) return(NULL)
                          
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
                            S = S,
                            S_syn = S_syn,
                            pi_S = pi_S
                          )

                        }

# Stop cluster when done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 2.108996 mins

setDT(ag_age_S_dt)

# # checks
# ag_age_S_dt <- merge(ag_age_S_dt, unique(pangraph_anno[,.(gene_family, number_genomes)]), all.x = TRUE, by = "gene_family")
# ag_age_S_dt[freq != number_genomes]# should be nothing

# PolyFastA vs kaks piS
with(ag_age_S_dt,
     plot(mean_ks, pi_S,
          pch = 16, col = rgb(0.1,0.1, 0.4, alpha = 0.5),
          cex = 0.95, 
          xlab = expression("kaks'"~pi[S]),
          ylab = expression("PolyFastA's"~pi[S])))
abline(a=0, b=1, lty=2, col="red")

# Calcualte piS outliers --------------------------------------------------

with(ag_age_S_dt,
     hist(pi_S,
          breaks = 200,
          xlab = expression("PolyFastA's"~pi[S])))
abline(v = median(ag_age_S_dt$pi_S), lty=2, col="red")

IQR_threshold = median(ag_age_S_dt$pi_S) + (IQR(ag_age_S_dt$pi_S) * 3)

abline(v = IQR_threshold, lty=2, col="red")

ag_age_S_dt[, pi_out := fcase(pi_S > IQR_threshold, 1,
                              default = 0)]

fwrite(ag_age_S_dt[pi_S > IQR_threshold, .(gene_family)],
       paste0(outdir_dat, "/piS_3IQR_outliers.csv"))



# Downsample genes using a hypergeometric distribution ----------------------------------------

# Examine the distribution first
quantile(ag_age_S_dt[pi_out == 0, m], probs = c(0.1, 0.25, 0.5, 0.75, 0.9))
# 10%    25%    50%    75%    90% 
# 228.0  372.0  699.0 1125.0 1622.7 


# Define standardised target length 
m_target <- quantile(ag_age_S_dt[pi_out == 0, m], probs = c(0.25))


# For genes LONGER than m_target, downsample segregating sites
# by drawing m_target sites without replacement from m total sites,
# of which S are segregating — i.e. sample from hypergeometric

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

set.seed(42)

ag_age_S_dt$hyperg_S_syn <- mapply(
  downsample_segsites,
  S       = ag_age_S_dt$S_syn,
  m       = ag_age_S_dt$m,
  m_target = m_target
)

# Standardised length column, all genes now treated as m_target long
ag_age_S_dt$hyperg_m <- ifelse(ag_age_S_dt$m >= m_target, 
                            m_target, 
                            NA)


# Export observes seg sites, gene length and piS --------------------------

fwrite(ag_age_S_dt, paste0(outdir_dat, "/ag_age_S_dt.csv"))
# ag_age_S_dt <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

unlink(fasta_dir, recursive = TRUE)


# Examine gene length -----------------------------------------------------
IQR_gene_threshold = median(ag_age_S_dt[pi_out==0]$m) + (IQR(ag_age_S_dt[pi_out==0]$m) * 3)

ag_age_S_dt[, genel_out := fcase(m > IQR_gene_threshold, 1,
                                 default = 0)]

par(mfrow = c(2,1))

h1 <- hist(ag_age_S_dt$m,
           breaks = 50,
           main ="",
           col="dodgerblue",
           border = "dodgerblue4",
           xlab = "Gene length")

abline(v = IQR_gene_threshold, lwd =1.1, col = "red",
       xpd = FALSE)

u <- par("usr")

text(u[2],u[4], expression("With"~pi[S]~"outliers"),
     pos = 2,
     xpd = TRUE)

hist(ag_age_S_dt[pi_out==0]$m,
     breaks = 50,
     xlim = c(0, max(h1$breaks)),
     main ="",
     col="dodgerblue",
     border = "dodgerblue4",
     xlab = "Gene length")

u <- par("usr")

text(u[2],u[4], expression("Without"~pi[S]~"outliers"),
     pos = 2,
     xpd = TRUE)

abline(v = IQR_gene_threshold, lwd =1.1, col = "red",
       xpd = FALSE)

par(mfrow = c(1,1))

with(ag_age_S_dt[pi_out==0],
     plot(freq, m,
          pch = 16,
          bty = "L", 
          ylab = "Gene length",
          xlab = "Genome frequency",
          col = rgb(.1,.1,.1,alpha = 0.4)))

abline(h = median(ag_age_S_dt[pi_out==0, m]),
       lty = 2, 
       lwd =1.1, col = "red",
       xpd = FALSE)

text(max(ag_age_S_dt$freq)/2,
     median(ag_age_S_dt[pi_out==0, m]),
     "median",
     pos = 3, col = "red"
)

abline(h = IQR_gene_threshold,
       lty = 2, 
       lwd =1.1, col = "dodgerblue4",
       xpd = FALSE)

text(max(ag_age_S_dt$freq)/2,
     IQR_gene_threshold,
     "median + 3*IQR",
     pos = 3, col = "dodgerblue4"
)

