# Calculate theta from core genes

core_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                   select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                              "number_genomes", "start","end", "average_dose",
                              "strand","ST","ag_type"))


# remove get only core genes and remove paralogs
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

# Get gene alignments and control for fission/fusion -----------------------

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
fasta_dir <- paste0(outdir_dat, "/core_fus_alns/")
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

sites_4 <- walk_genes(fasta_dir, list_unique_cores, mode = "4fold")
sites_2 <- walk_genes(fasta_dir, list_unique_cores, mode = "2fold")

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
end.time - start.time # Time difference of 4.122814 mins

fwrite(sites_combined_dt, paste0(outdir_dat, "/core_aln_sites.csv"))


# per-gene binomial thinning + pi / thetaW summary -----------------------
# sites_dt: must contain columns n, A, C, G, T, gene_length, gene_family, class
# m_target: target gene length in bp (e.g. 375)
# min_sites: minimum retained sites for a gene to be kept (default 2)

summarise_genes_hypergeometric <- function(sites_dt, m_target = 375, min_sites = 2) {
  
  # drop genes shorter than m_target
  sites_dt <- sites_dt[gene_length >= m_target]
  
  # per-gene processing 
  sites_dt[, {
    Ls <- .N
    
    Ls_keep <- rhyper(1, Ls, gene_length[1] - Ls, 375)

    if (Ls_keep < min_sites) {
      # return an empty row; filtered out below
      .(Ls = Ls, Ls_keep = Ls_keep, pi_s = NA_real_, thetaW_s = NA_real_,
        n_median = NA_real_, seg_s_sites = NA_integer_)
    } else {
      idx      <- sample.int(Ls, size = Ls_keep)
      n_j      <- n[idx]
      counts   <- cbind(A[idx], C[idx], G[idx], T[idx])
      
      p_mat    <- counts / n_j
      h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
      seg      <- rowSums(counts > 0) >= 2
      
      n_max    <- max(n_j)
      a_lookup <- c(NA_real_, cumsum(1 / seq_len(n_max - 1)))
      a_nj     <- a_lookup[n_j]
      
      .(Ls        = as.integer(Ls),
        Ls_keep   = as.integer(Ls_keep),
        pi_s      = as.numeric(mean(h)),
        thetaW_s  = as.numeric(mean(seg / a_nj)),
        n_median  = as.numeric(median(n_j)),
        seg_s_sites = as.integer(sum(seg)))
    }
  }, by = .(gene_family, class, gene_length)]
}


core_theta_dt <- summarise_genes_hypergeometric(sites_combined_dt)


# Save output
fwrite(core_theta_dt, paste0(outdir_dat, "/core_theta_dt.csv"))

#output 2Nu

mean_theta <- core_theta_dt[, .(mean_pi_s = mean(pi_s),
                  mean_thetaW_s = mean(thetaW_s)), class]

mut_rate = rbind(mean_theta, 
      cbind(class = "unweighted_mean", 
            mean_theta[, lapply(.SD, mean), .SDcols = c("mean_pi_s", "mean_thetaW_s")]))

fwrite(mut_rate, paste0(outdir_dat, "/core_mut_rate.csv"))


# remove additional fasta files
unlink(fasta_dir, recursive = TRUE)

png(paste0(outdir_fig,"/2Nmu_estimators.png"),
    width = 13, height = 15, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

# Plot comparison
with(core_theta_dt,
     plot(pi_s,  thetaW_s,
          pch = 16, col = rgb(0.1,0.2,0.4, alpha = 0.5),
          cex = 0.5, 
          xlab = "Synonymous nucleotide diversity (πS)", 
          ylab = "Watterson's θS",
          main = "Comparison of 2Nμ estimates")
)
abline(0, 1, col = "black", lty = 2)
abline(with(core_theta_dt,
            lm(thetaW_s~pi_s)),
       col = "red", lty = 2)
dev.off()


