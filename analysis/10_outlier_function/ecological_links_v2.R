# Ecological distribution (v2)
# The SPARK dataset records isolation origin (hospital vs non-hospital) and,
# for human isolates, disease vs carriage. treeWAS tests whether AG presence
# tracks either category beyond what the phylogeny explains, then we ask
# whether significant genes are outliers or sit in an outlier co-occurrence
# cluster (clean_coin_ags, from goldfinder.R).

library(data.table)
library(ape)
library(treeWAS)

if (!exists("outdir_dat")) outdir_dat <- "./output/data"

# Metadata ----------------------------------------------------------------

pan_anno_genos <- unique(fread(paste0(outdir_dat, "/pangenome_anno.csv"),
                               select = c("geno_id")))

spark_meta <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv",
                    select = c("id", "group_summary", "GROUP", "place_sampling",
                               "TYPE"))

spark_meta <- spark_meta[id %chin% pan_anno_genos$geno_id]
stopifnot(nrow(spark_meta) == 260)

spark_meta[, hospital_assoc := fcase(TYPE == "human hospital", "hospital",
                                     default = "non_hospital")]

# Hospital isolates by sampling site --------------------------------------
# confirms the 118/167 single-hospital figure in the text

hosp_sites <- spark_meta[hospital_assoc == "hospital",
                         .(n = .N), by = place_sampling][order(-n)]
hosp_sites[, prop := round(n / sum(n), 3)]
print(hosp_sites)
cat("Hospital isolates:", hosp_sites[, sum(n)],
    "| largest single site:", hosp_sites[1, n], "\n")

# Gene presence/absence matrix --------------------------------------------

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))[number_genomes < 260][number_genomes > 12]

gene_pa <- dcast(
  data = pan_anno[, .(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

pa_all <- as.matrix(gene_pa[, -1])
row.names(pa_all) <- gene_pa$geno_id
pa_all[pa_all > 1] <- 1L                # paralogues collapse to presence
storage.mode(pa_all) <- "integer"

tree <- ape::read.tree("./input_data/goldfinder/pang_gub_tree.newick")

# treeWAS wrapper ---------------------------------------------------------
# phen: 0/1 vector named by isolate. Tree is pruned to those isolates, snps
# rows matched to sub_tree$tip.label, invariant genes dropped, and
# n.snps.sim set from the filtered matrix.

run_treewas <- function(phen) {
  sub_tree <- keep.tip(tree, names(phen))
  
  snps <- pa_all[sub_tree$tip.label, , drop = FALSE]
  n_pres <- colSums(snps)
  snps <- snps[, n_pres > 0 & n_pres < nrow(snps), drop = FALSE]   # drop invariant
  
  phen <- phen[sub_tree$tip.label]
  stopifnot(identical(names(phen), row.names(snps)), !anyNA(phen))
  
  cat("isolates:", nrow(snps), "| variable genes:", ncol(snps),
      "| phen = 1:", sum(phen), "\n")
  
  set.seed(42)
  treeWAS(
    snps  = snps,
    phen  = phen,
    tree  = sub_tree,
    test  = c("terminal", "simultaneous", "subsequent"),
    n.snps.sim = 10 * ncol(snps),
    p.value = 0.05,
    p.value.correct = "fdr",
    snps.reconstruction = "parsimony",
    phen.reconstruction = "parsimony",
    plot.tree = FALSE
  )
}

# significant genes from each test; treeWAS returns a character message
# rather than a data.frame when a test finds nothing
sig_genes <- function(out) {
  tests <- c("terminal", "simultaneous", "subsequent")
  res <- rbindlist(lapply(tests, function(t) {
    s <- out[[t]]$sig.snps
    if (!is.data.frame(s) || nrow(s) == 0) return(NULL)
    data.table(gene_family = row.names(s), test = t)
  }))
  if (nrow(res) == 0) return(data.table(gene_family = character(), test = character()))
  res[, .(test = paste(sort(unique(test)), collapse = " & ")), by = gene_family]
}

# 1a. Hospital vs non-hospital (all 260) ----------------------------------

phen_hosp <- setNames(as.integer(spark_meta$hospital_assoc == "hospital"),
                      spark_meta$id)

out_hosp <- run_treewas(phen_hosp)
sig_hosp <- sig_genes(out_hosp)
print(sig_hosp)

# 1b. Disease vs carriage (171 human isolates) ----------------------------

meta_dc <- spark_meta[grepl("carriage|disease", group_summary)]
meta_dc[, disease := fifelse(grepl("carriage", group_summary), 0L, 1L)]
stopifnot(nrow(meta_dc) == 171)

phen_dc <- setNames(meta_dc$disease, meta_dc$id)

out_dc <- run_treewas(phen_dc)
sig_dc <- sig_genes(out_dc)
print(sig_dc)

# 2. Are significant genes outliers, or in an outlier cluster? ------------
# clean_coin_ags: one row per goldfinder gene, with the current outlier call
# (low / high / ns; NA if not a set 1 AG with a neutrality test)

clean_coin_ags <- fread(paste0(outdir_dat, "/clean_coin_ags.csv"),
                        select = c("gene_family", "p_bh", "call", "score", "clust_id"),
                        na.strings = c("", "NA"))   # blank call = no neutrality test

# guard against stray whitespace and call/score swaps
clean_coin_ags[, `:=`(call = trimws(call), score = trimws(score))]
stopifnot(clean_coin_ags[!is.na(call), all(call %chin% c("low", "high", "ns"))],
          !clean_coin_ags[, any(score %chin% c("low", "high", "ns"))])

# a cluster counts as an outlier group if any member is a low/high outlier
clust_out <- clean_coin_ags[, .(n_genes   = .N,
                                n_low     = sum(call == "low", na.rm = TRUE),
                                n_high    = sum(call == "high", na.rm = TRUE)),
                            by = clust_id]
clust_out[, outlier_clust := n_low + n_high > 0]

# outlier status of every tested AG, not only goldfinder genes
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites", "expected_S", "p_bh"))
gene_pvals[, call := fcase(p_bh <= 0.05 & seg_sites > expected_S, "high",
                           p_bh <= 0.05 & seg_sites < expected_S, "low",
                           default = "ns")]

check_sig <- function(sig) {
  if (nrow(sig) == 0) { cat("no significant genes\n"); return(invisible(NULL)) }
  x <- gene_pvals[, .(gene_family, freq, call)][sig, on = .(gene_family)]
  x <- clean_coin_ags[, .(gene_family, clust_id)][x, on = .(gene_family)]
  x <- clust_out[x, on = .(clust_id)]
  x[, .(gene_family, test, freq, call, clust_id, n_genes, n_low, n_high, outlier_clust)]
}

check_sig(sig_hosp)
check_sig(sig_dc)

# Raw (non-phylogenetic) per-gene Fisher's tests were run in v1 and found
# many hospital-associated genes; these are not reported because they do not
# control for clonal structure. See ecological_links.R if needed.
