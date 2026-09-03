# Calculate theta from core genes

# codon classification ----------------------------------------------------

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


# the chunk maker ---------------------------------------------------------

walk_chunked <- function(fasta_path, mode = "4fold", mem_budget_mb = 600) {
  dna   <- readDNAStringSet(fasta_path)
  L_nt  <- unique(width(dna))
  n_seq <- length(dna)
  stopifnot(length(L_nt) == 1, L_nt %% 3 == 0)
  
  bytes_per_cell <- 16
  max_chunk_nt   <- (mem_budget_mb * 1024^2) %/% (n_seq * bytes_per_cell)
  chunk_codons   <- max(1L, max_chunk_nt %/% 3L)
  chunk_nt       <- min(chunk_codons * 3L, L_nt)
  chunk_starts   <- seq(1L, L_nt, by = chunk_nt)
  n_chunks       <- length(chunk_starts)
  message(sprintf("[%s] %d seq x %d nt -> %d chunks of up to %d codons",
                  mode, n_seq, L_nt, n_chunks, chunk_codons))
  
  results <- vector("list", n_chunks)
  for (ci in seq_along(chunk_starts)) {
    start <- chunk_starts[ci]
    end   <- min(start + chunk_nt - 1L, L_nt)
    chunk <- as.matrix(subseq(dna, start = start, end = end))
    chunk[!chunk %in% c("A","C","G","T","-")] <- "N"
    results[[ci]] <- codon_walker(chunk, mode = mode)
    rm(chunk); gc(verbose = FALSE)
    message(sprintf("  [%s] chunk %d / %d (positions %d-%d)",
                    mode, ci, n_chunks, start, end))
  }
  do.call(rbind, results)
}

# get summary stats for pi & thetaW, variance, SE, 95% CI -----------------

summarise_sites <- function(sites_dt, n_for_var, label) {
  n_j   <- sites_dt$n
  quantile_n_j = quantile(n_j)
  counts  <- as.matrix(sites_dt[, .(A, C, G, T)])
  p_mat <- counts / n_j
  h     <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
  seg   <- rowSums(counts > 0) >= 2
  
  n_max    <- max(n_j)
  a_lookup <- c(NA, cumsum(1 / seq_len(n_max - 1)))
  a_nj     <- a_lookup[n_j]
  
  pi_hat    <- mean(h)
  theta_hat <- mean(seg / a_nj)
  L         <- nrow(sites_dt)
  
  a1 <- sum(1 / seq_len(n_for_var - 1))
  a2 <- sum(1 / seq_len(n_for_var - 1)^2)
  var_theta <- theta_hat / (a1 * L) + theta_hat^2 * a2 / (a1^2 * L)
  var_pi    <- (n_for_var + 1) / (3 * L * (n_for_var - 1)) * theta_hat +
    2 * (n_for_var^2 + n_for_var + 3) /
    (9 * L * n_for_var * (n_for_var - 1)) * theta_hat^2
  se_pi    <- sqrt(var_pi)
  se_theta <- sqrt(var_theta)
  
  data.table(
    class      = label,
    L          = L,
    q0.25      = quantile_n_j[2],
    q0.25      = quantile_n_j[2],
    pi         = pi_hat,
    var_pi     = var_pi,
    se_pi      = se_pi,
    pi_lo95    = pi_hat - 1.96 * se_pi,
    pi_hi95    = pi_hat + 1.96 * se_pi,
    thetaW     = theta_hat,
    var_theta  = var_theta,
    se_theta   = se_theta,
    theta_lo95 = theta_hat - 1.96 * se_theta,
    theta_hi95 = theta_hat + 1.96 * se_theta
  )
}

# Run ---------------------------------------------------------------------
core_aln <- "./input_data/bootstrapped_pirate_gubbins/recomb_mask_core_alignment.fasta"

start.time <- Sys.time()

sites_4 <- walk_chunked(core_aln, mode = "4fold")
sites_2 <- walk_chunked(core_aln, mode = "2fold")

sites_4_dt <- as.data.table(sites_4); setnames(sites_4_dt, c("n","A","C","G","T"))
sites_4_dt[, class := "4fold"]

sites_2_dt <- as.data.table(sites_2); setnames(sites_2_dt, c("n","A","C","G","T"))
sites_2_dt[, class := "2fold"]

sites_combined_dt <- rbind(sites_4_dt, sites_2_dt)

end.time <- Sys.time()
end.time - start.time # Time difference of 4.701017 mins

# fwrite(sites_combined_dt, paste0(outdir_dat, "/core_aln_sites.csv"))
fwrite(sites_combined_dt, paste0(outdir_dat, "/no_recomb_core_aln_sites.csv"))

n_for_var <- 260   # max sample size; adjust if your alignment differs
summary_dt <- rbind(
  summarise_sites(sites_4_dt,        n_for_var, "4fold"),
  summarise_sites(sites_2_dt,        n_for_var, "2fold"),
  summarise_sites(sites_combined_dt, n_for_var, "combined")
)

print(summary_dt)

# fwrite(summary_dt, paste0(outdir_dat, "/core_aln_summary_dt.csv"))
fwrite(summary_dt, paste0(outdir_dat, "/no_recomb_core_aln_summary_dt.csv"))

