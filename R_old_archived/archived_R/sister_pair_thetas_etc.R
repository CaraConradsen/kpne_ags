# GET THOSE CHERRIES
# What to determine whether including higher frequency AGs includes more recombination
# than our core genes

# Get Gubbins Pangraph core tree -----------------------------------------------------------
# Load tree
tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")

ntip     <- Ntip(tree)
internal <- (ntip + 1):(ntip + tree$Nnode)   

desc      <- Descendants(tree, internal, type = "tips")  
size      <- lengths(desc)
clades    <- lapply(desc, function(i) tree$tip.label[i])

# group clades by their size
clades_by_size <- split(clades, size)
# 
# clades_by_size[["2"]]   # all cherries
# clades_by_size[["3"]]   # all triplets
# 


# Calculate theta from core genes -----------------------------------------

core_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                   select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                              "number_genomes", "start","end", "average_dose",
                              "strand","ST","ag_type"))


# remove get only core genes and remove paralogs
core_anno <- core_anno[ag_type=="core"][average_dose <= 1]

# number of genomes in pangenome
tot_pangenome_size = length(unique(core_anno$geno_id))

list_unique_cores = unique(core_anno$gene_family)


# import m target

m_target <- as.numeric(readLines(paste0(outdir_dat, "/m_target.txt")))


# Prefix tables -----------------------------------------------------------

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

# gene walker -- reads each gene ONCE, slices across all frequencies -------

walk_genes_cherries <- function(fasta_dir, unique_cores, clades_by_size,
                                mode = "4fold",
                                n_cores = parallel::detectCores() - 4){
  
  freqs <- names(clades_by_size)
  
  per_gene <- function(focal_gene){
    
    temp_string <- tryCatch(
      readDNAStringSet(paste0(fasta_dir, focal_gene, ".fasta")),
      error = function(e) return(NULL)
    )
    if (is.null(temp_string)) return(NULL)   # per-gene skip, not whole run
    
    gene_length <- unique(width(temp_string))
    
    # convert once
    gene_mat <- as.matrix(temp_string)
    gene_mat[!gene_mat %in% c("A","C","G","T","-")] <- "N"
    
    # strip trailing _<digits> suffix so headers match clade members exactly
    core_names <- sub("_[0-9]+$", "", rownames(gene_mat))
    
    # loop frequencies, then clades within each frequency
    gene_out <- list()
    for (k in freqs){
      focal_genomes_sets <- clades_by_size[[k]]
      
      codon_counts <- lapply(focal_genomes_sets, function(set){
        sub_mat <- gene_mat[core_names %in% set, , drop = FALSE]
        if (nrow(sub_mat) < 2) return(NULL)
        cbind(codon_walker(sub_mat, mode = mode), gene_length, focal_gene)
      })
      
      codon_counts <- Map(function(codon_counts, i){
        if (is.null(codon_counts)) return(NULL)
        cbind(codon_counts, rep = i)
      }, codon_counts, seq_along(codon_counts))
      
      codon_counts <- do.call(rbind, codon_counts)
      if (!is.null(codon_counts))
        gene_out[[k]] <- cbind(codon_counts, freq = k)
    }
    
    do.call(rbind, gene_out)
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(Biostrings))
    parallel::clusterExport(cl,
                            c("codon_walker", "fourfold_prefixes",
                              "twofold_purine_prefixes", "twofold_pyrimidine_prefixes",
                              "clades_by_size", "freqs", "mode", "fasta_dir"),
                            envir = environment())
    results <- parallel::parLapply(cl, unique_cores, per_gene)
  } else {
    results <- parallel::mclapply(unique_cores, per_gene, mc.cores = n_cores)
  }
  
  do.call(rbind, results)
}

# Run ---------------------------------------------------------------------

fasta_dir <- paste0(outdir_dat, "/core_fus_alns/")   # core fusion alignments must exist

start.time <- Sys.time()
sites_4 <- walk_genes_cherries(fasta_dir, list_unique_cores, clades_by_size,
                              mode = "4fold")
sites_2 <- walk_genes_cherries(fasta_dir, list_unique_cores, clades_by_size,
                               mode = "2fold")
end.time <- Sys.time()
end.time - start.time # Time difference of 41.7497 mins

sites_4_dt <- as.data.table(sites_4); setnames(sites_4_dt, c("n","A","C","G","T", "gene_length", "gene_family", "rep", "freq"))
cols = c("n","A","C","G","T", "gene_length")
sites_4_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_4_dt[, class := "4fold"]

sites_2_dt <- as.data.table(sites_2); setnames(sites_2_dt, c("n","A","C","G","T", "gene_length", "gene_family", "rep", "freq"))
cols = c("n","A","C","G","T", "gene_length")
sites_2_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_2_dt[, class := "2fold"]

sites_combined_dt <- rbind(sites_4_dt, sites_2_dt)



fwrite(sites_combined_dt, paste0(outdir_dat,"/sister_core_theta.csv"))
# sites_combined_dt <- fread(paste0(outdir_dat,"/sister_core_theta.csv"))


# Get frequency based core ------------------------------------------------


summarise_mutli_genes_hypergeometric <- function(sites_dt, m_target, min_sites = 2) {
  
  # drop genes shorter than m_target
  sites_dt <- sites_dt[gene_length >= m_target]
  
  # per-gene processing 
  sites_dt[, {
    Ls <- .N
    
    Ls_keep <- rhyper(1, Ls, gene_length[1] - Ls, m_target)
    
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
  }, by = .(gene_family, class,rep, freq, gene_length)]
}


core_theta_multi_dt <- summarise_mutli_genes_hypergeometric(sites_combined_dt, m_target)


# Save output
fwrite(core_theta_multi_dt, paste0(outdir_dat, "/core_theta_multi_dt.csv"))
# core_theta_multi_dt <- fread(paste0(outdir_dat, "/core_theta_multi_dt.csv"))

theta_plot = core_theta_multi_dt[, .(thetaW = mean(thetaW_s)), by = c("gene_family", "freq", "rep")]

temp_box <- with(theta_plot, 
                 boxplot(thetaW~freq, pch = 16, xlim=c(0,46),
                         cex =0.5, main = expression(theta[W]~"per core sub-clades"),
                         xlab = "Pangenome frequency", 
                         ylab = expression(theta[W]),
                         col = "dodgerblue"))#ylim = c(0,0.1),
abline(h=0.02, col = "red")
text(45, 0.02, expression(theta[W]~"= 0.02"),
     pos = 3, col = "red", xpd = TRUE)



