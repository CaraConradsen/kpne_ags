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

# step 1, write all FASTAs first 
fasta_dir <- paste0(outdir_dat, "/ag_fus_alns/")
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


# step 2, create codon walker ----------------------------------------

# codon classification 

fourfold_prefixes           <- c("CT","GT","TC","CC","AC","GC","CG","GG")
twofold_purine_prefixes     <- c("AA","AG","CA","GA","TT")              # A/G synonymous
twofold_pyrimidine_prefixes <- c("AA","AG","CA","GA","TT","TA","TG")    # C/T synonymous

# codon-walker ------------------------------------------------------------

codon_walker <- function(aln, mode = c("4fold","2fold")) {
  mode  <- match.arg(mode)
  L     <- ncol(aln)
  n_cod <- L / 3L
  out   <- vector("list", n_cod)
  
  for (i in seq_len(n_cod)) {
    sub  <- aln[, (3*i-2):(3*i)]
    keep <- rowSums(sub == "N" | sub == "-") == 0
    sub  <- sub[keep, , drop = FALSE]
    if (nrow(sub) < 2) next
    
    prefix <- paste0(sub[,1], sub[,2])
    if (length(unique(prefix)) > 1) next       # mixed 1st/2nd position -> skip
    p     <- prefix[1]
    third <- sub[, 3]
    
    if (mode == "4fold") {
      if (!(p %in% fourfold_prefixes)) next
    } else {                                   # 2fold
      bases  <- unique(third)
      is_pur <- all(bases %in% c("A","G")) && (p %in% twofold_purine_prefixes)
      is_pyr <- all(bases %in% c("C","T")) && (p %in% twofold_pyrimidine_prefixes)
      if (!is_pur && !is_pyr) next
    }
    
    counts   <- tabulate(match(third, c("A","C","G","T")), 4)
    out[[i]] <- c(n = nrow(sub), counts)
  }
  do.call(rbind, out)
}


walk_genes <- function(fasta_dir, unique_cores, 
                       mode = "4fold"){
  
  results <- vector("list", length(unique_cores))
  
  for(i in seq_along(unique_cores)){
    
    focal_gene = unique_cores[i]
    
    # Read the fixed FASTA, skip if missing
    temp_string <- tryCatch(
      readDNAStringSet(paste0(fasta_dir, focal_gene, ".fasta")),
      error = function(e) return(NULL)
    )
    
    # skip if file is missing
    if (is.null(temp_string)) return(NULL)
    
    gene_length  <- unique(width(temp_string))
    
    gene_mat <- as.matrix(temp_string)
    gene_mat[!gene_mat %in% c("A","C","G","T","-")] <- "N"
    
    results[[i]] <- cbind(codon_walker(gene_mat, mode = mode), gene_length, focal_gene)
  }
  do.call(rbind, results)
}



# Run ---------------------------------------------------------------------

start.time <- Sys.time()

sites_4 <- walk_genes(fasta_dir, list_unique_ags, mode = "4fold")
sites_2 <- walk_genes(fasta_dir, list_unique_ags, mode = "2fold")

sites_4_dt <- as.data.table(sites_4); setnames(sites_4_dt, c("n","A","C","G","T", "gene_length", "gene_family"))
cols = c("n","A","C","G","T", "gene_length")
sites_4_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_4_dt[, class := "4fold"]

sites_2_dt <- as.data.table(sites_2); setnames(sites_2_dt, c("n","A","C","G","T", "gene_length", "gene_family"))
cols = c("n","A","C","G","T", "gene_length")
sites_2_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_2_dt[, class := "2fold"]

sites_combined_dt <- rbind(sites_4_dt, sites_2_dt)

end.time <- Sys.time()
end.time - start.time # Time difference of 2.537345 mins

fwrite(sites_combined_dt, paste0(outdir_dat, "/ag_aln_sites.csv"))

# remove additional fasta files
unlink(fasta_dir, recursive = TRUE)

# Identify pi outliers ----------------------------------------------------
# per-gene processing 
per_ag <- sites_combined_dt[, {
  Ls <- .N
  n_j      <- n
  counts   <- cbind(A, C, G, T)
  
  p_mat    <- counts / n_j
  h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
  seg      <- rowSums(counts > 0) >= 2
  
  n_max    <- max(n_j)
  a_lookup <- c(NA_real_, cumsum(1 / seq_len(n_max - 1)))
  a_nj     <- a_lookup[n_j]
  
  .(Ls        = as.integer(Ls),
    pi_s      = as.numeric(mean(h)),
    thetaW_s  = as.numeric(mean(seg / a_nj)),
    n_median  = as.numeric(median(n_j)),
    seg_s_sites = as.integer(sum(seg)))
  
}, by = .(gene_family, class, gene_length)]


combined_pi_dt <- per_ag[, .(
  pi_s     = mean(pi_s),
  gene_length = mean(gene_length),
  Ls_total = sum(Ls),
  seg_total = sum(seg_s_sites)
), by = gene_family]

IQR_threshold = median(combined_pi_dt$pi_s) + 3*IQR(combined_pi_dt$pi_s)
quantile(combined_pi_dt$pi_s, probs = 0.9)

fwrite(combined_pi_dt[pi_s > IQR_threshold, .(gene_family)],
       paste0(outdir_dat, "/piS_3IQR_outliers.csv"))

fwrite(combined_pi_dt, paste0(outdir_dat, "/ag_age_S_dt.csv"))


# Inspect distributions ---------------------------------------------------

par(mfrow=c(1,2))

with(combined_pi_dt,
     hist(pi_s, 
          breaks = 500,
          yaxt = "n",
          ylab = "Frequency",
          xlab = expression(pi[S]),
          col = "black", 
          main = "Dist. of pis"))



segments(IQR_threshold,0,
         lwd = 1.1, lty = 3,
         IQR_threshold, 2000,
         col = "blue")
text(IQR_threshold, 1500, "median + 3*IQR: 0.10",
     col = "blue", pos = 4)

axis(side = 2, at = seq(0,3000, 500),
     labels = seq(0,3000, 500), las = 2)

red_ags = combined_pi_dt[pi_s <= IQR_threshold]
red_quantiles = quantile(red_ags$gene_length, probs = c(0.1, 0.25))

with(red_ags,
     hist(gene_length, 
          breaks = 500,
          yaxt = "n",
          ylab = "Frequency",
          xlab = "Gene length",
          col = "black", 
          main = "Dist. of gene length,\nafter removing pis outliers"))

axis(side = 2, at = seq(0,200, 50),
     labels = seq(0,200, 50), las = 2)

lapply(red_quantiles, function(x) {
  abline(v = x, col = "red", lty = 2)})

text(red_quantiles[[2]], 150, "Q2: 375bp",
     col = "red", pos = 4)


# Calculate AG segregating sites for target gene size ------------------------


ag_genes_hypergeometric <- function(sites_dt, m_target = 375, min_sites = 2) {
  
  # drop genes shorter than m_target
  sites_dt <- sites_dt[gene_length >= m_target]
  
  # per-gene processing 
  sites_dt[, {
    Ls <- .N
    
    Ls_keep <- rhyper(1, Ls, gene_length[1] - Ls, 375)
    
    if (Ls_keep < min_sites) {
      # return an empty row; filtered out below
      .(Ls = Ls, Ls_keep = Ls_keep, pi_s = NA_real_, 
        n_median = NA_real_, seg_s_sites = NA_integer_)
    } else {
      idx      <- sample.int(Ls, size = Ls_keep)
      n_j      <- n[idx]
      counts   <- cbind(A[idx], C[idx], G[idx], T[idx])
      
      p_mat    <- counts / n_j
      h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
      seg      <- rowSums(counts > 0) >= 2
      
      .(Ls        = as.integer(Ls),
        Ls_keep   = as.integer(Ls_keep),
        pi_s      = as.numeric(mean(h)),
        n_median  = as.numeric(median(n_j)),
        seg_s_sites = as.integer(sum(seg)))
    }
  }, by = .(gene_family, class, gene_length)]
}


ag_seg_sites_dt <- ag_genes_hypergeometric(sites_combined_dt)


# add back frequency
ag_seg_sites_dt <- merge(ag_seg_sites_dt,
                         unique(pangraph_anno[, .(gene_family, number_genomes)]),
                         all.x = TRUE, by = "gene_family")

# Save output
fwrite(ag_seg_sites_dt, paste0(outdir_dat, "/ag_seg_sites_dt.csv"))
