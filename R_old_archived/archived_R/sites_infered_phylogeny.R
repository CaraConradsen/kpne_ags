# --- Part 1: pull synonymous-site columns (bases, not counts) per gene ---

# reuse your prefix tables
fourfold_prefixes <- c("CT","GT","TC","CC","AC","GC","CG","GG")

extract_synonymous_cols <- function(aln_mat) {
  # aln_mat: strains x positions, values A/C/G/T/-/N
  L <- ncol(aln_mat)
  n_cod <- L %/% 3L
  keep_cols <- integer(0)
  
  for (i in seq_len(n_cod)) {
    cols <- (3*i-2):(3*i)
    sub  <- aln_mat[, cols, drop = FALSE]
    prefix <- paste0(sub[,1], sub[,2])
    # only fourfold codons where 1st/2nd positions are unambiguous and consistent
    good <- prefix %in% fourfold_prefixes
    if (!any(good)) next
    if (length(unique(prefix[good])) != 1) next   # mixed prefix -> skip
    keep_cols <- c(keep_cols, cols[3])            # 3rd position only
  }
  aln_mat[, keep_cols, drop = FALSE]
}

fasta_dir <- paste0(outdir_dat, "/core_fus_alns/")

# all strain names from the tree tips define the column set
all_strains <- tree$tip.label

per_gene_syn <- lapply(list_unique_cores, function(g) {
  f <- paste0(fasta_dir, g, ".fasta")
  if (!file.exists(f)) return(NULL)
  dna <- tryCatch(readDNAStringSet(f), error = function(e) NULL)
  if (is.null(dna)) return(NULL)
  
  m <- as.matrix(dna)
  m[!m %in% c("A","C","G","T","-")] <- "N"
  rownames(m) <- sub("_[0-9]+$", "", rownames(m))   # match tip labels
  
  syn <- extract_synonymous_cols(m)
  if (ncol(syn) == 0) return(NULL)
  
  # align rows to the full strain set; missing strains -> N
  out <- matrix("N", nrow = length(all_strains), ncol = ncol(syn),
                dimnames = list(all_strains, NULL))
  shared <- intersect(rownames(syn), all_strains)
  out[shared, ] <- syn[shared, ]
  out
})

per_gene_syn <- per_gene_syn[!vapply(per_gene_syn, is.null, logical(1))]

# concatenate across genes (strains aligned by row)
syn_concat <- do.call(cbind, per_gene_syn)
cat("synonymous alignment:", nrow(syn_concat), "strains x",
    ncol(syn_concat), "sites\n")

# --- Part 2: build the tree ---

syn_phydat <- phyDat(syn_concat, type = "DNA")

# distance + NJ for a fast look-see (what Adam asked for)
dm <- dist.ml(syn_phydat, model = "JC69")
nj_tree <- NJ(dm)

write.tree(nj_tree, paste0(outdir_dat, "/synonymous_NJ_tree.tre"))



plot(nj_tree, type = "fan", cex = 0.3, no.margin = TRUE)


# NJ can produce tiny negative branch lengths — set them to 0
nj_tree$edge.length[nj_tree$edge.length < 0] <- 0

# 260 tips: a fan (circular) layout is the most legible
plot(nj_tree, type = "fan", cex = 0.3, no.margin = TRUE, label.offset = 0.0002)
add.scale.bar(length = 0.001, cex = 0.7)
mtext("synonymous substitutions per site (dS)", side = 1, cex = 0.7)

