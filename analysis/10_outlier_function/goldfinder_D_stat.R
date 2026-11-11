# Phylogenetic D for the genes Goldfinder placed in a cluster
#
# D is a property of a gene's presence/absence pattern, not of the scoring
# method that found it, so it is computed ONCE per gene family. The earlier
# version looped separately over the simultaneous and the subsequent gene
# lists and rbind()ed the results, which computed D twice for the 135 genes
# found by both and wrote them out as duplicate rows. It then overwrote the
# label on the second half with a flat score = "subsequent", so a gene found
# by both appeared once as simultaneous and once as subsequent.
#
# It also filtered to p_bh <= 0.05 before computing, so D only ever existed
# for genes that were outliers at the time, and went stale as soon as the
# outlier set changed. D is now computed for every cluster member, and the
# outlier call is joined on afterwards.

d_reps <- 1000L   # phylo.d permutations; the caper default

all_associations <- fread(paste0(outdir_dat, "/clean_coin_ags.csv"))

cat("computing D for", nrow(all_associations), "gene families\n")

# Presence/absence and tree -----------------------------------------------

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"),
                  select = c("geno_id", "gene_family"))

pa <- dcast(pan_anno, geno_id ~ gene_family,
            fun.aggregate = length, value.var = "gene_family")

# phylo.d needs a binary variable. fun.aggregate = length gives a count, so
# a gene present twice in a genome would enter as 2. Fission and fusion loci
# are collapsed upstream, which should make this impossible; check rather
# than assume.
stopifnot(max(as.matrix(pa[, -"geno_id"])) <= 1L)

tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")

ag_presence_absence <- as.data.frame(pa)

com_dat_ag_ap <- comparative.data(tree, ag_presence_absence, geno_id)

# genes dropped when the tree and the matrix were matched cannot be tested
testable <- intersect(all_associations$gene_family,
                      names(com_dat_ag_ap$data))

if (length(testable) < nrow(all_associations))
  warning(sprintf("%d of %d cluster members are absent from the matched data and are skipped",
                  nrow(all_associations) - length(testable), nrow(all_associations)))

# D statistic -------------------------------------------------------------

DEstimates <- rbindlist(lapply(testable, function(g) {
  d_info <- do.call(caper::phylo.d,
                    list(data = com_dat_ag_ap, binvar = as.name(g),
                         permut = d_reps))
  data.table(gene_family = d_info$binvar,
             d           = d_info$DEstimate,
             p_rand      = d_info$Pval1,
             p_brow      = d_info$Pval0)
}))

# Join the cluster and outlier information --------------------------------
# One row per gene family, so score keeps the collapsed
# "simultaneous & subsequent" label rather than being overwritten.

DEstimates <- all_associations[DEstimates, on = .(gene_family)]

stopifnot(DEstimates[, uniqueN(gene_family)] == nrow(DEstimates))

cat("\nD computed for", nrow(DEstimates), "genes\n")
print(DEstimates[, .(genes = .N), by = score][order(-genes)])
print(DEstimates[!is.na(call), .(genes = .N), by = .(score, call)][order(score, call)])

fwrite(DEstimates, paste0(outdir_dat, "/DEstimates.csv"))
