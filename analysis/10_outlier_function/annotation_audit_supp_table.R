# ########################################################## ####
# Annotation audit of outlier AGs -> Supplementary Table      ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####

# Script purpose ####
# Applies the annotation-quality criteria uniformly to every AG family tested,
# so that the families highlighted in Fig. X are selected by stated rules and
# not by inspection. Writes one row per tested family with each criterion as its
# own column, so a reader can see what passed, what failed, and why.

# Criteria (see Methods):
#  i.   >= 90% of per-locus annotations agree with each other and the consensus
#  ii.  product is specific, not a domain/family/hypothetical assignment
#  iii. no same-named or same-product family at >= 0.95 * n_genomes (paralogue of a core gene)
#  iv.  no sibling family with a near-disjoint genome set summing to near-universal (split family)
#  v.   < 10% of member loci carry a mobile element annotation

n_genomes   <- 260
hf_cut      <- 0.95 * n_genomes   # high-frequency sibling threshold, 247
split_cut   <- 0.95 * n_genomes   # union at which a set of siblings is called a split
jacc_cut    <- 0.05               # max genome-set overlap for siblings to count as a split
agree_cut   <- 0.90
mge_cut     <- 0.10
alpha_fdr   <- 0.05

# generic product patterns, criterion ii
generic_pat <- paste("domain-containing", "family protein", "hypothetical",
                     "DUF[0-9]+", "putative", "uncharacteri[sz]ed",
                     "-like protein", "^protein$",
                     sep = "|")


# Load ---------------------------------------------------------------------------------

anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
              select = c("geno_id", "locus_tag", "gene_family", "freq", "ag_type",
                         "consensus_gene_name", "consensus_product",
                         "gene", "product", "anchor", "msu",
                         "phg_state", "integron", "ice", "is_fam"))

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

# Goldfinder co-inheritance clusters, used as an annotation column only
simul_clust <- fread(paste0(outdir_dat, "/simul_clust_fams.csv"))


# Family-level summary -------------------------------------------------------------------
# one row per family, with the modal per-locus product and its share

fam <- anno[, .(n_loci        = .N,
                freq          = freq[1],
                ag_type       = ag_type[1],
                cons_name     = consensus_gene_name[1],
                cons_product  = consensus_product[1],
                anchor        = names(sort(table(anchor), decreasing = TRUE))[1],
                msu           = names(sort(table(msu),    decreasing = TRUE))[1]),
            by = gene_family]

# modal per-locus product and the proportion of loci carrying it
top_prod <- anno[!is.na(product),
                 .N, by = .(gene_family, product)
                 ][order(gene_family, -N)
                   ][, .(top_product = product[1],
                         top_n       = N[1],
                         n_distinct  = .N,
                         second      = fifelse(.N > 1, product[2], NA_character_),
                         second_n    = fifelse(.N > 1, N[2], NA_integer_)),
                     by = gene_family]

fam <- merge(fam, top_prod, by = "gene_family", all.x = TRUE)
fam[, top_frac := round(top_n / n_loci, 3)]

# modal per-locus gene symbol, kept for the audit trail only; Bakta symbols are
# transferred at protein-cluster level and are not used to assign identity
top_gene <- anno[!is.na(gene) & gene != "",
                 .N, by = .(gene_family, gene)
                 ][order(gene_family, -N)
                   ][, .(top_gene = gene[1]), by = gene_family]

fam <- merge(fam, top_gene, by = "gene_family", all.x = TRUE)


# Criterion i: per-locus annotation concordance -------------------------------------------

fam[, crit_i_concordant := top_frac >= agree_cut & top_product == cons_product]


# Criterion ii: product specificity --------------------------------------------------------

fam[, crit_ii_specific := !grepl(generic_pat, top_product,      ignore.case = TRUE) &
                          !grepl(generic_pat, cons_product,     ignore.case = TRUE) &
                          !is.na(cons_name) & cons_name != "" & cons_name != "unitary"]


# Criterion v: mobile element content ------------------------------------------------------
# computed here because criteria iii and iv need the genome sets built below

mge <- anno[, .(mge_frac   = round(mean(!is.na(phg_state) | !is.na(integron) |
                                        !is.na(ice) | !is.na(is_fam)), 3),
                phage_frac = round(mean(!is.na(phg_state)), 3)),
            by = gene_family]

fam <- merge(fam, mge, by = "gene_family", all.x = TRUE)
fam[, crit_v_not_mge := mge_frac < mge_cut]


# Criterion iii: high-frequency sibling -----------------------------------------------------
# a family whose gene name or product is also carried by a near-universal family is
# a paralogue of a core gene and cannot be labelled with that name

fam[, name_base := tstrsplit(cons_name, "_", fixed = TRUE, keep = 1L)]

hf <- fam[freq >= hf_cut, .(gene_family, name_base, cons_product, freq)]

# by shared gene-name stem
hf_by_name <- merge(fam[, .(gene_family, name_base)],
                    hf[!is.na(name_base) & name_base != "" & name_base != "unitary",
                       .(sib = gene_family, name_base, sib_freq = freq)],
                    by = "name_base", allow.cartesian = TRUE)
hf_by_name <- hf_by_name[gene_family != sib]

# by shared consensus product, excluding generic products
hf_by_prod <- merge(fam[!grepl(generic_pat, cons_product, ignore.case = TRUE),
                        .(gene_family, cons_product)],
                    hf[, .(sib = gene_family, cons_product, sib_freq = freq)],
                    by = "cons_product", allow.cartesian = TRUE)
hf_by_prod <- hf_by_prod[gene_family != sib]

hf_sibs <- rbindlist(list(hf_by_name[, .(gene_family, sib, sib_freq)],
                          hf_by_prod[, .(gene_family, sib, sib_freq)]),
                     use.names = TRUE)
hf_sibs <- unique(hf_sibs)

hf_sum <- hf_sibs[, .(hf_sib = paste0(sib, "(", sib_freq, ")", collapse = ";"),
                      n_hf_sib = .N),
                  by = gene_family]

fam <- merge(fam, hf_sum, by = "gene_family", all.x = TRUE)
fam[is.na(n_hf_sib), n_hf_sib := 0L]
fam[, crit_iii_no_hf_sib := n_hf_sib == 0L]


# Criterion iv: split family -----------------------------------------------------------------
# siblings sharing a name or product but occupying near-disjoint sets of genomes are
# fragments of one family; their union approaches the genome count

gsets <- anno[, .(genos = list(unique(geno_id))), by = gene_family]
setkey(gsets, gene_family)

# candidate sibling pairs: shared name stem or shared non-generic product, any frequency
cand_name <- merge(fam[!is.na(name_base) & name_base != "" & name_base != "unitary",
                       .(gene_family, name_base)],
                   fam[!is.na(name_base) & name_base != "" & name_base != "unitary",
                       .(sib = gene_family, name_base)],
                   by = "name_base", allow.cartesian = TRUE)[gene_family != sib]

cand_prod <- merge(fam[!grepl(generic_pat, cons_product, ignore.case = TRUE),
                       .(gene_family, cons_product)],
                   fam[!grepl(generic_pat, cons_product, ignore.case = TRUE),
                       .(sib = gene_family, cons_product)],
                   by = "cons_product", allow.cartesian = TRUE)[gene_family != sib]

cand <- unique(rbindlist(list(cand_name[, .(gene_family, sib)],
                              cand_prod[, .(gene_family, sib)]),
                         use.names = TRUE))

# overlap of genome sets for each candidate pair
cand[, `:=`(g_set = gsets[.SD$gene_family, genos],
            s_set = gsets[.SD$sib,         genos])]

cand[, `:=`(n_shared = mapply(function(a, b) length(intersect(a, b)), g_set, s_set),
            n_min    = mapply(function(a, b) min(length(a), length(b)),  g_set, s_set),
            n_sib    = lengths(s_set))]

# a sibling counts towards a split only if it is near-disjoint
cand[, disjoint := n_shared <= jacc_cut * n_min]

split_sum <- cand[disjoint == TRUE,
                  .(split_parts = paste0(sib, "(", n_sib, ")", collapse = ";"),
                    n_split_part = .N,
                    split_extra  = sum(n_sib)),
                  by = gene_family]

fam <- merge(fam, split_sum, by = "gene_family", all.x = TRUE)
fam[is.na(n_split_part), `:=`(n_split_part = 0L, split_extra = 0L)]
fam[, split_union := freq + split_extra]
fam[, crit_iv_not_split := !(split_union >= split_cut & split_union > freq)]


# Assemble the table ----------------------------------------------------------------------

fam[, eligible := crit_i_concordant & crit_ii_specific & crit_iii_no_hf_sib &
                  crit_iv_not_split & crit_v_not_mge]

# first failing criterion, so the table says why a family was excluded
fam[, fail_reason := fifelse(eligible, "",
                      fifelse(!crit_i_concordant,  "i: per-locus annotations disagree",
                       fifelse(!crit_ii_specific,  "ii: non-specific product or gene name",
                        fifelse(!crit_iii_no_hf_sib, "iii: high-frequency sibling present",
                         fifelse(!crit_iv_not_split, "iv: fragmented family",
                                                     "v: mobile element content")))))]

supp <- merge(gene_pvals[, .(gene_family, k = freq, seg_sites, syn_sites,
                             expected_S, pi_s, p_bh)],
              fam, by = "gene_family", all.x = TRUE)

supp[, s_over_e := round(seg_sites / expected_S, 2)]
supp[, direction := fifelse(p_bh >= alpha_fdr, "ns",
                     fifelse(seg_sites > expected_S, "high", "low"))]

# Goldfinder cluster as an annotation column, not a filter
supp <- merge(supp,
              unique(simul_clust[, .(gene_family,
                                     gf_cluster = assoc_cluster,
                                     gf_cluster_n = n_genes_assoc_cluster)]),
              by = "gene_family", all.x = TRUE)

setcolorder(supp, c("gene_family", "cons_name", "top_gene", "cons_product", "top_product",
                    "k", "ag_type", "anchor", "msu",
                    "syn_sites", "seg_sites", "expected_S", "s_over_e", "pi_s",
                    "p_bh", "direction",
                    "n_loci", "n_distinct", "top_n", "top_frac", "second", "second_n",
                    "crit_i_concordant", "crit_ii_specific", "crit_iii_no_hf_sib",
                    "crit_iv_not_split", "crit_v_not_mge",
                    "eligible", "fail_reason",
                    "hf_sib", "n_hf_sib", "split_parts", "n_split_part", "split_union",
                    "mge_frac", "phage_frac", "gf_cluster", "gf_cluster_n"))

supp[, c("name_base", "freq") := NULL]
setorderv(supp, c("direction", "p_bh", "k"), c(1, 1, -1))


# Counts for the Results text --------------------------------------------------------------

cat("families tested:              ", nrow(supp), "\n")
cat("significant (BH <", alpha_fdr, "):       ", supp[direction != "ns", .N], "\n")
cat("  high / low:                 ", supp[direction == "high", .N], "/",
                                      supp[direction == "low",  .N], "\n")
cat("significant and eligible:     ", supp[direction != "ns" & eligible == TRUE, .N], "\n")
cat("  high / low:                 ", supp[direction == "high" & eligible == TRUE, .N], "/",
                                      supp[direction == "low"  & eligible == TRUE, .N], "\n\n")

cat("exclusions among significant families:\n")
print(supp[direction != "ns" & eligible == FALSE, .N, by = fail_reason][order(-N)])


# Write -------------------------------------------------------------------------------------

fwrite(supp, paste0(outdir_dat, "/supp_table_annotation_audit.csv"))

wb <- createWorkbook()
addWorksheet(wb, "Annotation audit")
writeData(wb, "Annotation audit", supp, headerStyle = createStyle(textDecoration = "bold"))
freezePane(wb, "Annotation audit", firstActiveRow = 2, firstActiveCol = 2)
setColWidths(wb, "Annotation audit", cols = 1:ncol(supp), widths = "auto")
saveWorkbook(wb, paste0(outdir_dat, "/supp_table_annotation_audit.xlsx"), overwrite = TRUE)
