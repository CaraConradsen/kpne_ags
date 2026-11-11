# Goldfinder non-independence: map AG outlier calls onto the clusters
#
# Goldfinder itself is not re-run. This reads its cluster output and attaches
# the current outlier calls from gene_pvals_fdr.csv.
#
# Cluster numbering is continuous across the two scoring methods. The
# simultaneous association clusters keep their own ids, the simultaneous
# dissociation pairs form one further cluster, and the subsequent clusters
# are offset to continue the sequence, so no id is shared between methods.
#
# A gene family found by both methods is kept ONCE, labelled
# "simultaneous & subsequent", and carries the simultaneous cluster as its
# clust_id. Both ids are retained in clust_id_sim and clust_id_sub.

sep_score <- " & "   # how a gene found by both methods is labelled

# Parse Goldfinder cluster files ------------------------------------------
# Each cluster is a ">id,n_genes" header followed by one gene per line, as
# "<something>__<gene_family>".

parse_clusters <- function(f) {

  x <- read.delim(f, header = FALSE)

  h <- which(startsWith(x$V1, ">"))

  # end of each cluster = one before the next header; last runs to the end
  rows <- Map(seq, h + 1L, c(h[-1] - 1L, nrow(x)))

  dt <- rbindlist(lapply(seq_along(rows), function(g)
    data.table(assoc_cluster = gsub(">", "", x$V1[h[g]]),
               gene_family   = x$V1[rows[[g]]])))

  dt[, c("assoc_cluster", "n_genes_assoc_cluster") :=
       tstrsplit(assoc_cluster, ",", fill = TRUE)]
  dt[, c("assoc_cluster", "n_genes_assoc_cluster") :=
       lapply(.SD, as.integer), .SDcols = c("assoc_cluster", "n_genes_assoc_cluster")]

  dt[, gene_family := gsub(",", "", gene_family)]
  dt[, gene_family := tstrsplit(gene_family, "__", fill = TRUE, keep = 2)]

  dt[]
}

sim <- parse_clusters("./input_data/goldfinder/simul_output/association_clusters.txt")
sub <- parse_clusters("./input_data/goldfinder/sub_output/association_clusters.txt")

# dissociating pairs, simultaneous scoring only
dis_pairs <- fread("./input_data/goldfinder/simul_output/simultaneous_dissociation_significant_pairs.csv")
dis_pairs[, Gene_1 := tstrsplit(Gene_1, "__", fill = TRUE, keep = 2)]
dis_pairs[, Gene_2 := tstrsplit(Gene_2, "__", fill = TRUE, keep = 2)]

# Continuous cluster numbering --------------------------------------------

n_sim      <- max(sim$assoc_cluster)
dis_id     <- n_sim + 1L
sub_offset <- dis_id

sim[, `:=`(clust_id = assoc_cluster, method = "simultaneous",
           assc_type = "associate")]

dis <- data.table(
  assoc_cluster         = dis_id,
  gene_family           = unique(c(dis_pairs$Gene_1, dis_pairs$Gene_2)),
  n_genes_assoc_cluster = uniqueN(c(dis_pairs$Gene_1, dis_pairs$Gene_2)),
  clust_id              = dis_id,
  method                = "simultaneous",
  assc_type             = "disassociate")

sub[, `:=`(clust_id = assoc_cluster + sub_offset, method = "subsequent",
           assc_type = "associate")]

gf_long <- rbind(sim, dis, sub)

stopifnot(gf_long[, uniqueN(method), by = .(clust_id)][, all(V1 == 1)])

cat("clusters:", uniqueN(gf_long$clust_id), "numbered 1 -",
    max(gf_long$clust_id), "\n")
cat("  simultaneous:", gf_long[method == "simultaneous", uniqueN(clust_id)],
    " subsequent:",   gf_long[method == "subsequent",   uniqueN(clust_id)], "\n")

fwrite(gf_long, paste0(outdir_dat, "/goldfinder_clusters_long.csv"))

# simul_clust_fams.csv is still read by annotation_audit_supp_table.R as an
# annotation column, so it is kept, now written from the same parse rather
# than from a separate copy of this code.
fwrite(gf_long[method == "simultaneous",
               .(assoc_cluster, gene_family, n_genes_assoc_cluster,
                 assc_type, score = method)],
       paste0(outdir_dat, "/simul_clust_fams.csv"))

fwrite(gf_long[method == "subsequent",
               .(assoc_cluster, gene_family, n_genes_assoc_cluster,
                 assc_type, score = method)],
       paste0(outdir_dat, "/sub_clust_fams.csv"))

# One row per gene family -------------------------------------------------
# A gene found by both methods is collapsed rather than counted twice. It
# keeps the simultaneous cluster, since that is the stronger evidence, and
# the subsequent cluster is retained alongside it.

gf <- gf_long[, .(
  score        = paste(sort(unique(method)), collapse = sep_score),
  clust_id     = clust_id[method == "simultaneous"][1],
  clust_id_sim = clust_id[method == "simultaneous"][1],
  clust_id_sub = clust_id[method == "subsequent"][1],
  assc_type    = paste(sort(unique(assc_type)), collapse = sep_score)
), by = gene_family]

# genes found only by the subsequent method have no simultaneous cluster
gf[is.na(clust_id), clust_id := clust_id_sub]

stopifnot(!anyNA(gf$clust_id), gf[, uniqueN(gene_family)] == nrow(gf))

# cluster size, after collapsing
gf[, n_genes_assoc_cluster := .N, by = clust_id]

cat("\ngene families, counted once each:\n")
print(gf[, .(genes = .N), by = score][order(-genes)])

# Set 1 membership --------------------------------------------------------

set1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

gf[, n_in_set1   := sum(gene_family %chin% set1_ags), by = clust_id]
gf[, perc_in_set1 := round((n_in_set1 / n_genes_assoc_cluster) * 100)]

# Attach the current outlier calls ----------------------------------------

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites", "syn_sites",
                               "expected_S", "p_two", "p_bh"))

# Direction of the departure, and whether it passed the FDR threshold. This
# is what the network and D-statistic panels colour by. It must be derived
# from the CURRENT p-values: the file shipped before this carried a column
# pi_sel built from an earlier outlier set.
gene_pvals[, sel_type := fifelse(seg_sites > expected_S,
                                    "balancing", "directional")]
gene_pvals[, call := fcase(p_bh <= 0.05 & seg_sites > expected_S, "high",
                           p_bh <= 0.05 & seg_sites < expected_S, "low",
                           default = "ns")]

setnames(gene_pvals, "freq", "freq")

clean_coin_ags <- gene_pvals[gf, on = .(gene_family)]

# Goldfinder ran on the whole pangenome, so some cluster members are not
# set 1 AGs and have no neutrality test.
cat("\nGoldfinder genes:", nrow(clean_coin_ags),
    "| with a neutrality test:", clean_coin_ags[!is.na(p_bh), .N], "\n")
print(clean_coin_ags[!is.na(call), .(genes = .N), by = .(score, call)][order(score, call)])

fwrite(clean_coin_ags, paste0(outdir_dat, "/clean_coin_ags.csv"))

# Fisher's exact test: is coincidence associated with selection? ----------
# Run overall, then separately for each direction against the
# non-outliers. The two directions are opposite hypotheses, so pooling them
# asks only whether a cluster member is unusual in either direction, and a
# category enriched for both would cancel. Each direction is compared with
# the NON-OUTLIERS, never with the opposite direction.
#
# Over set 1 AGs above 5% frequency, as before.

pan_anno <- unique(fread(paste0(outdir_dat, "/pangenome_anno.csv"),
                         select = c("gene_family", "freq")))

pan_anno <- pan_anno[freq >= (260 * 0.05)]

pan_anno[, gf_hits := fifelse(gene_family %chin% gf$gene_family, "Y", "N")]

pan_anno <- gene_pvals[, .(gene_family, call)][pan_anno, on = .(gene_family)]

# set 1 AGs only, and only those with a neutrality test
pan_anno_set1 <- pan_anno[gene_family %chin% set1_ags & !is.na(call)]

fisher_coincidence <- function(dt, label) {

  tab <- table(factor(dt$gf_hits, levels = c("N", "Y")),
               factor(dt$sel_grp, levels = c("unsel", "sel")))

  res <- fisher.test(tab)

  data.table(comparison = label,
             n_sel      = sum(dt$sel_grp == "sel"),
             n_unsel    = sum(dt$sel_grp == "unsel"),
             sel_in_gf  = tab["Y", "sel"],
             unsel_in_gf = tab["Y", "unsel"],
             OR         = unname(res$estimate),
             CI_lower   = round(res$conf.int[1], 4),
             CI_upper   = round(res$conf.int[2], 4),
             pval       = round(res$p.value, 4))
}

# overall: any outlier against the non-outliers
overall <- copy(pan_anno_set1)[, sel_grp := fifelse(call != "ns", "sel", "unsel")]

# each direction against the non-outliers, the opposite direction dropped
high_v_ns <- pan_anno_set1[call != "low"][, sel_grp := fifelse(call == "high", "sel", "unsel")]
low_v_ns  <- pan_anno_set1[call != "high"][, sel_grp := fifelse(call == "low",  "sel", "unsel")]

fexact_coincident_genes <- rbind(
  fisher_coincidence(overall,   "any outlier vs non-outlier"),
  fisher_coincidence(high_v_ns, "high diversity vs non-outlier"),
  fisher_coincidence(low_v_ns,  "low diversity vs non-outlier"))

fexact_coincident_genes[, p_bh := round(p.adjust(pval, "BH"), 4)]

cat("\ncoincidence vs selection, set 1 AGs at >= 5% frequency:\n")
print(fexact_coincident_genes)

fwrite(fexact_coincident_genes, paste0(outdir_dat, "/fexact_coincident_genes.csv"))
