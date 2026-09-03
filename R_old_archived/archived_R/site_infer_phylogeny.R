fourfold_prefixes           <- c("CT","GT","TC","CC","AC","GC","CG","GG")
twofold_purine_prefixes     <- c("AA","AG","CA","GA","TT")
twofold_pyrimidine_prefixes <- c("AA","AG","CA","GA","TT","TA","TG")

extract_syn_cols <- function(aln_mat) {   # returns fold4 and fold2 column matrices
  L <- ncol(aln_mat); n_cod <- L %/% 3L
  c4 <- list(); c2 <- list()
  for (i in seq_len(n_cod)) {
    sub    <- aln_mat[, (3*i-2):(3*i), drop = FALSE]
    prefix <- paste0(sub[,1], sub[,2])
    pref_ok <- prefix[!grepl("[N-]", prefix)]
    if (length(unique(pref_ok)) != 1L) next     # ambiguous/mixed 1st-2nd pos -> skip
    p     <- pref_ok[1]
    third <- sub[, 3]                            # named by strain -> rownames survive cbind
    
    if (p %in% fourfold_prefixes) {
      c4[[length(c4)+1]] <- third               # every change here is synonymous
    } else if (p %in% c(twofold_purine_prefixes, twofold_pyrimidine_prefixes)) {
      obs <- unique(third[!third %in% c("N","-")])
      keep <-
        (all(obs %in% c("A","G")) && p %in% twofold_purine_prefixes) ||
        (all(obs %in% c("C","T")) && p %in% twofold_pyrimidine_prefixes)
      if (keep) c2[[length(c2)+1]] <- third
    }
  }
  list(fold4 = if (length(c4)) do.call(cbind, c4) else NULL,
       fold2 = if (length(c2)) do.call(cbind, c2) else NULL)
}

align_to_all <- function(mat, all_strains) {
  if (is.null(mat) || ncol(mat) == 0) return(NULL)
  out <- matrix("N", length(all_strains), ncol(mat),
                dimnames = list(all_strains, NULL))
  shared <- intersect(rownames(mat), all_strains)
  out[shared, ] <- mat[shared, ]
  out
}

per_gene <- lapply(list_unique_cores, function(g) {
  f <- paste0(fasta_dir, g, ".fasta")
  if (!file.exists(f)) return(NULL)
  dna <- tryCatch(readDNAStringSet(f), error = function(e) NULL)
  if (is.null(dna)) return(NULL)
  m <- as.matrix(dna)
  m[!m %in% c("A","C","G","T","-")] <- "N"
  rownames(m) <- sub("_[0-9]+$", "", rownames(m))
  sc <- extract_syn_cols(m)
  list(fold4 = align_to_all(sc$fold4, all_strains),
       fold2 = align_to_all(sc$fold2, all_strains))
})

bind_class <- function(pg, cls) {
  mats <- lapply(pg, function(x) if (is.null(x)) NULL else x[[cls]])
  mats <- mats[!vapply(mats, is.null, logical(1))]
  do.call(cbind, mats)
}
syn4 <- bind_class(per_gene, "fold4")
syn2 <- bind_class(per_gene, "fold2")
cat("4-fold sites:", ncol(syn4), " 2-fold sites:", ncol(syn2), "\n")

library(phangorn)
d4 <- as.matrix(dist.ml(phyDat(syn4, type = "DNA"), model = "JC69"))
d2 <- as.matrix(dist.ml(phyDat(syn2, type = "DNA"), model = "JC69"))

# both matrices are indexed by all_strains, so order matches
d_mean <- (d4 + d2) / 2
# if any pair lacks comparable sites in one class, fall back to the other:
d_mean[is.na(d_mean)] <- pmax(d4, d2, na.rm = TRUE)[is.na(d_mean)]

nj_tree <- NJ(as.dist(d_mean))
nj_tree$edge.length[nj_tree$edge.length < 0] <- 0

pairs <- list(
  c("SPARK_1138_C2", "SPARK_276_C1"),
  c("SPARK_776_C1",  "SPARK_1965_C1"),
  c("SPARK_1298_C1", "SPARK_262_C2"),
  c("SPARK_1867_C1", "SPARK_1618_C1"),
  c("SPARK_2810_C1", "SPARK_418_C1"),
  c("SPARK_355_C1",  "SPARK_195_C2"),
  c("SPARK_2013_C1", "SPARK_67_C1")
)

# default all tips black, then colour each pair
tipcols <- rep("black", Ntip(nj_tree))
names(tipcols) <- nj_tree$tip.label

pal <- rainbow(length(pairs))           # or a fixed vector of 7 colours
for (i in seq_along(pairs)) tipcols[pairs[[i]]] <- pal[i]

plot(nj_tree, type = "fan", cex = 0.3, no.margin = TRUE,
     label.offset = 0.0002, tip.color = tipcols)
add.scale.bar(length = 0.001, cex = 0.7)

# Worst offenders
offenders <- theta_plot[thetaW > 0.4]
offenders[, genome_pairs := clades_by_size[["2"]][rep]]

