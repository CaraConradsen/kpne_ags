# O-antigen locus: gene content and diversity statistics, by serotype
#
# Every genome carries one O-antigen gene cluster in the interval between gnd
# (g000180) and hisF (g001818), both of which are present in all 260 genomes.
# The clusters occupying that interval differ by serotype and share little
# homology, so PIRATE scores them as different gene families. What looks in a
# pangenome like presence/absence of wzm, wzt, glf and wbbM is therefore the
# O1/O2 galactan cluster being present rather than one of the alternatives:
# the O9-family wzm/wzt pair, the O3 mannan genes, and so on.
#
# This script tabulates, for every serotype in the dataset, the gene families
# that make up its locus, with the neutrality-test statistics and the allele
# frequency spectrum of each. Gene names are the curated tidy_name, not the
# PIRATE consensus name, which disagrees with it for the galactan cluster
# (consensus calls wzm "rfbA", wzt "rfbB", glf "rfbD", wbbM "wbbM1") and which
# gives the name wzm to a second, unrelated family, g004721.
#
# Core genes shared by every serotype (ugd, hisIE, wcaG) are dropped, since
# they flank the variable cluster rather than belonging to it.

library(data.table)

min_prev  <- 0.5   # a family is part of a serotype's locus at this prevalence
min_n_sero <- 3    # serotypes with fewer genomes than this are reported but flagged

loci  <- fread(paste0(outdir_dat, "/olocus_genes.csv"))

# tidy_name holds the curated gene symbol for most families but the full
# product description for a few (g000264 is "Uridine diphosphate galacturonate
# 4-epimerase", not "wcaG"), so take the first thing that looks like a symbol
symb <- function(a, b, p) {
  ok <- function(x) !is.na(x) & nzchar(x) & nchar(x) <= 10
  fifelse(ok(a), a, fifelse(ok(b), b, substr(p, 1, 28)))
}
loci[, tidy_name := symb(tidy_name, consensus_gene_name, consensus_product)]
meta  <- fread(paste0(outdir_dat, "/olocus.csv"))        # geno_id, ST, O_locus, ...
pv    <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
               select = c("gene_family", "freq", "seg_sites", "syn_sites",
                          "expected_S", "p_bh"))
alst  <- fread(paste0(outdir_dat, "/ag_allele_stats.csv"),
               select = c("gene_family", "n_alleles", "H", "call"))
gf    <- fread(paste0(outdir_dat, "/goldfinder_clusters_long.csv"))
adt   <- fread(paste0(outdir_dat, "/allele_counts_by_ST.csv"))

# 1. allele frequency spectrum, p1 to p3 ------------------------------------
# Row sums of allele_counts_by_ST.csv give each allele's total across STs, so
# the spectrum is counted over exactly the carriers the test used.

st_cols <- grep("^ST", names(adt), value = TRUE)
spec <- adt[, .(gene_family, allele_name, n = rowSums(.SD)), .SDcols = st_cols][n > 0]
setorder(spec, gene_family, -n)
spec <- spec[, .(k_al = sum(n),
                 p1 = n[1] / sum(n),
                 p2 = if (.N >= 2) n[2] / sum(n) else NA_real_,
                 p3 = if (.N >= 3) n[3] / sum(n) else NA_real_), by = gene_family]

# 2. which families make up each serotype's locus ---------------------------

loci <- meta[, .(geno_id, O_locus)][loci, on = .(geno_id)]
n_sero <- loci[, .(n_genomes = uniqueN(geno_id)), by = O_locus]

fam <- loci[, .(n = uniqueN(geno_id),
                tidy_name = tidy_name[1],
                consensus_gene_name = consensus_gene_name[1],
                consensus_product = consensus_product[1],
                pos = round(mean(idx), 1)),      # mean rank from gnd, for ordering
            by = .(O_locus, gene_family)]
fam <- n_sero[fam, on = .(O_locus)]
fam[, prev := n / n_genomes]
fam <- fam[prev >= min_prev]

# genes present in every serotype are the flanks, not the cluster
n_sero_tot <- uniqueN(fam$O_locus)
fam[, n_sero_with := .N, by = gene_family]       # serotypes carrying this family
ubiq <- fam[n_sero_with == n_sero_tot, unique(gene_family)]
cat("dropped as shared by all", n_sero_tot, "serotypes:",
    paste(unique(fam[gene_family %chin% ubiq, tidy_name]), collapse = ", "), "\n")
fam <- fam[!gene_family %chin% ubiq]

# 3. join the statistics -----------------------------------------------------

tab <- pv[fam, on = .(gene_family)]
tab <- alst[tab, on = .(gene_family)]
tab <- spec[tab, on = .(gene_family)]
tab <- gf[method == "simultaneous", .(gene_family, assoc_cluster)][tab, on = .(gene_family),
                                                                   mult = "first"]

tab[, s_per_site := seg_sites / syn_sites]
tab[, obs_exp    := seg_sites / expected_S]

setnames(tab, c("freq", "seg_sites", "syn_sites", "expected_S"),
         c("k", "S", "syn_sites", "exp_S"))

setcolorder(tab, c("O_locus", "n_genomes", "pos", "gene_family", "tidy_name",
                   "consensus_gene_name", "n", "prev", "n_sero_with", "assoc_cluster",
                   "k", "S", "syn_sites", "s_per_site", "exp_S", "obs_exp",
                   "call", "p_bh", "n_alleles", "H", "p1", "p2", "p3"))
setorder(tab, -n_genomes, O_locus, pos)

fwrite(tab, paste0(outdir_dat, "/o_antigen_locus_table.csv"))
cat("wrote", paste0(outdir_dat, "/o_antigen_locus_table.csv"),
    "  rows:", nrow(tab), "  serotypes:", uniqueN(tab$O_locus), "\n")

# 4. printed summary ---------------------------------------------------------

# families carried by O-locus genomes but absent from the neutrality test:
# worth knowing, because they are part of the locus and were not assessed
untested <- unique(tab[is.na(k), .(gene_family, tidy_name,
                                   genomes = NA_integer_)])
if (nrow(untested)) {
  untested[, genomes := loci[.SD, on = .(gene_family), uniqueN(geno_id),
                             by = .EACHI]$V1]
  cat("\nin the locus but absent from gene_pvals_fdr.csv (not in set 1):\n")
  print(untested[order(-genomes)])
}

fmt <- function(d) d[, .(gene = tidy_name, fam = gene_family, k, S, syn = syn_sites,
                         `s/site` = round(s_per_site, 3), exp = round(exp_S, 1),
                         `o/e` = round(obs_exp, 1), call,
                         alleles = n_alleles, H = round(H, 2),
                         p1 = round(p1, 2), p2 = round(p2, 2), p3 = round(p3, 2))]

for (o in n_sero[order(-n_genomes), O_locus]) {
  d <- tab[O_locus == o]
  if (!nrow(d)) next
  cat(sprintf("\n=== %s   (%d genomes, %d cluster genes)%s\n",
              o, d$n_genomes[1], nrow(d),
              if (d$n_genomes[1] < min_n_sero) "  [few genomes]" else ""))
  print(fmt(d))
}

# the two focal genes and everything they associate with
cat("\n=== Goldfinder cluster of wzm (g004051) and wbbM (g000322)\n")
cl <- gf[gene_family %chin% c("g004051", "g000322"), unique(assoc_cluster)]
mem <- gf[assoc_cluster %chin% cl, unique(gene_family)]
d <- unique(loci[gene_family %chin% mem,
                 .(gene_family, tidy_name, consensus_product)])
d <- alst[pv[d, on = .(gene_family)], on = .(gene_family)]
d <- spec[d, on = .(gene_family)]
d[, `:=`(s_per_site = seg_sites / syn_sites, obs_exp = seg_sites / expected_S)]
setnames(d, c("freq", "seg_sites", "syn_sites", "expected_S"),
         c("k", "S", "syn_sites", "exp_S"))
print(d[order(-k), .(gene = tidy_name, fam = gene_family, k, S, syn = syn_sites,
                     `s/site` = round(s_per_site, 3), exp = round(exp_S, 1),
                     `o/e` = round(obs_exp, 1), call, alleles = n_alleles,
                     H = round(H, 2), p1 = round(p1, 2), p2 = round(p2, 2),
                     p3 = round(p3, 2))])
cat("\nGoldfinder cluster(s):", paste(cl, collapse = ", "),
    " members in the O locus interval:", nrow(d), "of", length(mem), "\n")
