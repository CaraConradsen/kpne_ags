# Get AG piN (protein-level constraint)
#
# Companion to 06_piS_estimates_seg_sites/obs_ag_vals.R, and built the same
# way: walk every codon of every AG alignment, keep the sites whose class is
# unambiguous, and summarise base frequencies per gene.
#
# 06 keeps the third position of 4-fold and 2-fold degenerate codons, where a
# change is synonymous. Here the complement is kept: NONDEGENERATE (0-fold)
# positions, where every change alters the amino acid. piN is therefore the
# mirror of piS at 4-fold sites, and the two are on the same per-site scale,
# so piN/piS needs no site-counting correction of the Nei-Gojobori kind.
#
# input require data sets
# ags, and the FASTAs written by 06 (output/data/ag_fus_alns/)

pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                              select = c("geno_id", "gene_family", "locus_tag",
                                         "freq", "start","end", "max_dose",
                                         "strand","ST","ag_type")))


# remove estimates of core, paralogs and singletons
pangraph_anno <- pangraph_anno[ag_type!="core"][freq > 1][max_dose <= 1]

# number of genomes in pangenome
tot_pangenome_size = length(unique(pangraph_anno$geno_id))

list_unique_ags = unique(pangraph_anno$gene_family)

# one locus per genome per family after fission/fusion removal
stopifnot(pangraph_anno[, .N, by = .(gene_family, geno_id)][, all(N == 1)])

# set up directories

gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_hps0.6_out/feature_sequences/"

# step 1, the FASTAs ------------------------------------------------------
# 06 writes these; rebuilt here only if they are missing, so that the two
# scripts can never be reading different alignments

fasta_dir <- paste0(outdir_dat, "/ag_fus_alns/")

rebuild_fastas <- FALSE   # TRUE rewrites every alignment from the genomes

if (rebuild_fastas || !dir.exists(fasta_dir) ||
    length(list.files(fasta_dir, pattern = "\\.fasta$")) == 0) {

  dir.create(fasta_dir, showWarnings = FALSE)

  for (i in list_unique_ags) {

    temp_string <- tryCatch(
      readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
      error = function(e) NULL
    )
    if (is.null(temp_string)) next

    sampled_loci <- pangraph_anno[gene_family == i, locus_tag]
    temp_string  <- temp_string[names(temp_string) %in% sampled_loci]

    if (length(temp_string) < 2) next
    if (length(temp_string) %in% c(1, tot_pangenome_size)) next

    writeXStringSet(temp_string,
                    paste0(fasta_dir, i, ".fasta"),
                    format = "fasta")
  }
}


# step 2, codon degeneracy ------------------------------------------------
# Read off the standard genetic code rather than hard-coded, so the table
# cannot drift from Biostrings. fold[codon, position] is the number of bases
# at that position, the codon's own included, that leave the amino acid
# unchanged: 1 = nondegenerate, 2 = two-fold, 3 = three-fold (Ile), 4 = four-
# fold. Codons are indexed ci = (b1-1)*16 + (b2-1)*4 + b3 over A,C,G,T, so
# the walker can look them up with arithmetic instead of pasting strings.

base_idx <- c(A = 1L, C = 2L, G = 3L, T = 4L)

cod_tab <- local({
  gc    <- Biostrings::GENETIC_CODE
  b     <- c("A","C","G","T")
  fold  <- matrix(0L, 64L, 3L)
  stp   <- logical(64L)
  for (i1 in 1:4) for (i2 in 1:4) for (i3 in 1:4) {
    ci <- (i1-1L)*16L + (i2-1L)*4L + i3
    ch <- c(b[i1], b[i2], b[i3])
    cd <- paste0(ch, collapse = "")
    stp[ci] <- gc[[cd]] == "*"
    for (pos in 1:3) {
      alt <- vapply(b, function(x) { y <- ch; y[pos] <- x; paste0(y, collapse = "") }, "")
      fold[ci, pos] <- sum(gc[alt] == gc[[cd]])
    }
  }
  list(fold = fold, stop = stp)
})

# non-synonymous codon-walker ---------------------------------------------
# One pass per codon, then one test per position. A position is kept when
#   (i)  the other two positions are invariant across the retained sequences,
#        so the codon context, and with it the degeneracy, is well defined;
#   (ii) no observed codon at that site is a stop;
#   (iii) EVERY observed codon is nondegenerate at that position, so that any
#        change at the site is non-synonymous. Requiring all observed codons
#        rather than the majority one is the conservative choice: it drops
#        partially degenerate sites (first position of Leu and Arg) instead
#        of counting a synonymous change as non-synonymous.

codon_walker_n <- function(aln) {
  n_cod <- ncol(aln) / 3L
  out   <- vector("list", n_cod * 3L)
  k     <- 0L

  for (i in seq_len(n_cod)) {
    sub  <- aln[, (3*i-2):(3*i), drop = FALSE]
    keep <- rowSums(sub == "N" | sub == "-") == 0
    sub  <- sub[keep, , drop = FALSE]
    if (nrow(sub) < 2) next

    idx <- matrix(base_idx[sub], nrow = nrow(sub))
    if (anyNA(idx)) next

    for (pos in 1:3) {
      ctx <- switch(pos, c(2L,3L), c(1L,3L), c(1L,2L))
      if (!all(idx[, ctx[1]] == idx[1, ctx[1]])) next
      if (!all(idx[, ctx[2]] == idx[1, ctx[2]])) next

      u    <- unique(idx[, pos])
      tmpl <- idx[1, ]
      ci   <- integer(length(u))
      for (z in seq_along(u)) {
        t2 <- tmpl; t2[pos] <- u[z]
        ci[z] <- (t2[1]-1L)*16L + (t2[2]-1L)*4L + t2[3]
      }
      if (any(cod_tab$stop[ci]))            next
      if (any(cod_tab$fold[ci, pos] != 1L)) next

      k        <- k + 1L
      out[[k]] <- c(n = nrow(sub), tabulate(idx[, pos], 4L))
    }
  }
  do.call(rbind, out)
}


walk_genes_n <- function(fasta_dir, unique_cores){

  results <- vector("list", length(unique_cores))

  for(i in seq_along(unique_cores)){

    focal_gene = unique_cores[i]

    # Read the fixed FASTA, skip if missing
    temp_string <- tryCatch(
      readDNAStringSet(paste0(fasta_dir, focal_gene, ".fasta")),
      error = function(e) return(NULL)
    )

    # skip if file is missing
    if (is.null(temp_string)) next

    gene_length  <- unique(width(temp_string))
    if (length(gene_length) != 1) next
    if (gene_length %% 3 != 0) next

    gene_mat <- as.matrix(temp_string)
    gene_mat[!gene_mat %in% c("A","C","G","T","-")] <- "N"

    # a short or gap-heavy gene can leave no usable nondegenerate site, and
    # cbind() on a NULL would make a 2-column row that breaks the rbind below
    cw <- codon_walker_n(gene_mat)
    if (is.null(cw)) next

    results[[i]] <- cbind(cw, gene_length, focal_gene)
  }
  do.call(rbind, results)
}



# Run ---------------------------------------------------------------------

start.time <- Sys.time()

sites_0 <- walk_genes_n(fasta_dir, list_unique_ags)

sites_0_dt <- as.data.table(sites_0); setnames(sites_0_dt, c("n","A","C","G","T", "gene_length", "gene_family"))
cols = c("n","A","C","G","T", "gene_length")
sites_0_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_0_dt[, class := "0fold"]

end.time <- Sys.time()
end.time - start.time # Time difference of roughly 10 mins

fwrite(sites_0_dt, paste0(outdir_dat, "/ag_aln_sites_N.csv"))


# Per-gene piN ------------------------------------------------------------
# Same summary as 06, over nondegenerate sites instead of degenerate ones

summarise_genes_n <- function(sites_dt) {

  sites_dt[, {
    Ln <- .N
    n_j      <- as.integer(n)
    counts   <- cbind(A, C, G, T)

    p_mat    <- counts / n_j
    h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
    seg      <- rowSums(counts > 0) >= 2

    n_max    <- max(n_j)
    a_lookup <- c(NA_real_, cumsum(1 / seq_len(n_max - 1)))
    a_nj     <- a_lookup[n_j]

    .(Ln          = as.integer(Ln),
      pi_n        = as.numeric(mean(h)),
      thetaW_n    = as.numeric(mean(seg / a_nj)),
      n_median    = as.numeric(median(n_j)),
      seg_n_sites = as.integer(sum(seg)))
  }, by = .(gene_family, class, gene_length)]
}


ag_theta_n_dt <- summarise_genes_n(sites_0_dt)

# one row per gene already, since every site here is of one class
stopifnot(ag_theta_n_dt[, .N, by = gene_family][, all(N == 1)])

combined_pi_n_dt <- ag_theta_n_dt[, .(gene_family,
                                      pi_n,
                                      gene_length,
                                      nonsyn_sites = Ln,     # nondegenerate sites assayed
                                      seg_n_sites)]          # non-syn segregating sites

# add back freq
combined_pi_n_dt <- merge(combined_pi_n_dt,
                          unique(pangraph_anno[,.(gene_family,freq)]),
                          all.x = TRUE, by = "gene_family")

fwrite(combined_pi_n_dt, paste0(outdir_dat, "/ag_piN_dt.csv"))


# piN / piS ---------------------------------------------------------------
# 06 must have been run: ag_age_S_dt.csv holds pi_s, syn_sites and the
# synonymous segregating sites over the same alignments.

s_file <- paste0(outdir_dat, "/ag_age_S_dt.csv")

if (file.exists(s_file)) {

  pi_s_dt <- fread(s_file, select = c("gene_family","pi_s","syn_sites","seg_sites"))
  setnames(pi_s_dt, "seg_sites", "seg_s_sites")

  pi_dt <- merge(combined_pi_n_dt, pi_s_dt, all.x = TRUE, by = "gene_family")

  # undefined wherever the gene carries no synonymous variation, which is
  # most of the deficit set, so it is left NA rather than set to Inf
  pi_dt[, piN_piS := fifelse(pi_s > 0, pi_n / pi_s, NA_real_)]
  pi_dt[, seg_N_S := fifelse(seg_s_sites > 0, seg_n_sites / seg_s_sites, NA_real_)]

  fwrite(pi_dt, paste0(outdir_dat, "/ag_piN_piS_dt.csv"))

  cat(sprintf("AGs with piN: %d   with piN/piS: %d   median piN/piS: %.3f\n",
              nrow(pi_dt), pi_dt[!is.na(piN_piS), .N],
              pi_dt[, median(piN_piS, na.rm = TRUE)]))

} else {
  warning("ag_age_S_dt.csv not found: run 06 first for piN/piS")
  pi_dt <- combined_pi_n_dt
}


# Inspect distributions ---------------------------------------------------

par(mfrow=c(1,3))

with(pi_dt,
     hist(pi_n,
          breaks = 500,
          yaxt = "n",
          ylab = "Frequency",
          xlab = expression(pi[N]),
          col = "black",
          main = "Dist. of piN"))

axis(side = 2, at = seq(0,3000, 500),
     labels = seq(0,3000, 500), las = 2)

with(pi_dt,
     hist(nonsyn_sites,
          breaks = 500,
          yaxt = "n",
          ylab = "Frequency",
          xlab = "Nondegenerate sites assayed",
          col = "black",
          main = "Dist. of non-synonymous sites"))

axis(side = 2, at = seq(0,200, 50),
     labels = seq(0,200, 50), las = 2)

if ("piN_piS" %in% names(pi_dt)) {

  with(pi_dt[!is.na(piN_piS)],
       hist(log10(piN_piS + 1e-4),
            breaks = 200,
            yaxt = "n",
            ylab = "Frequency",
            xlab = expression(log[10](pi[N]/pi[S])),
            col = "black",
            main = "Dist. of piN/piS"))

  axis(side = 2, at = seq(0,1000, 200),
       labels = seq(0,1000, 200), las = 2)

  # neutrality, for reference: anything right of this is protein diversity
  # outrunning synonymous diversity
  abline(v = 0, col = "blue", lty = 3, lwd = 1.1)
  text(0, 100, "piN = piS", col = "blue", pos = 4)
}
