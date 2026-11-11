# Ecological distribution
# The SPARK dataset contains information on isolated origin 
# (in the broad categories: hospital-associated, non-hospital-associated). 

pan_anno_genos <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                                  select = c("geno_id")))

spark_meta <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv", 
                    select = c("id", "group_summary", "group_summary_l2", 
                               "GROUP", "ENVIRONMENTAL_GROUP", "place_sampling", "ORIGIN",
                               "specific_group", "ASSOCIATED_GROUP",
                               "ASSOCIATED_SPECIES", "GENERAL_ISOLATION_SOURCE",
                               "TYPE", "IN_OUT","SAMPLE_TYPE", "SAMPLE_DESCRIPTION"))



spark_meta <- spark_meta[id %chin% pan_anno_genos$geno_id]

spark_meta[, hospital_assoc := fcase(TYPE == "human hospital", "hospital",
                                      default = "non_hospital")]

tab <- table(spark_meta$GROUP, spark_meta$hospital_assoc)   # 3 x 2

chi <- chisq.test(tab)$statistic
n <- sum(tab)
V <- sqrt(chi / (n * (min(dim(tab)) - 1)))
V

# A. Host domain (GROUP) — Human 171 / Animal 68 / Environmental 21. Clean three-way split, no missing values. Binarise to Human vs non-human (171/89) if you want a single 2×2.
# B. Hospital-associated vs non-hospital-associated (from TYPE) — "human hospital" 167 vs everything else 93. This is your headline example and the split is well-powered.
# C. Disease 95 vs carriage 76 (72 hospital + 4 community); both arms well-powered.

# Step 1 treeWAS ----------------------------------------------------------

#  gene presence/absence matrix
pan_anno<- fread(paste0(outdir_dat, "/pangenome_anno.csv"))[ag_type != "core"][number_genomes > 12]

# pan_anno <- pan_anno[geno_id %chin% spark_meta[grepl("carriage|disease", group_summary), id],]

# Generate a binary matrix giving gene presence/absence
gene_pa <- dcast(
  data = pan_anno[,.(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

dim(gene_pa)

# # checks
# # Which rows have any value > 1
# gene_pa[rowSums(gene_pa[, -1] > 1) > 0, 1]
# 
# # Which columns have any value > 1
# names(gene_pa[, -1])[colSums(gene_pa[, -1] > 1) > 0]
# 

meta_carriage <- spark_meta[grepl("carriage|disease", group_summary)]
meta_carriage[, carriage := fcase(grepl("carriage", group_summary), "carriage",
                                     default = "disease")]

snps <- as.matrix(gene_pa[,-1])          # ensure it's a matrix, not data.frame
row.names(snps) <- gene_pa$geno_id
storage.mode(snps) <- "integer"

# 2. Phenotype = hospital association, named by isolate
# phen <- setNames(spark_meta$hospital_assoc, spark_meta$id) 
phen <- setNames(meta_carriage$carriage, meta_carriage$id) 

# 3. Tree = your Gubbins recombination-corrected phylogeny
tree <- ape::read.tree("./input_data/goldfinder/pang_gub_tree.newick")

sub_tree <- keep.tip(tree, meta_carriage$id)

ids <- sub_tree$tip.label

# Check: isolates in snps not in the tree, and vice versa
setdiff(rownames(snps), ids)
setdiff(ids, rownames(snps))   # should be empty

# Subset and order snps to match the tree
snps_sub <- snps[ids, , drop = FALSE]

# Drop genes that are now invariant (present in all or none)
n <- colSums(snps_sub)
snps_sub <- snps_sub[, n > 0 & n < nrow(snps_sub)]

# Match phen to the same isolates, coded 0/1
phen_sub <- setNames(as.numeric(phen[ids] == "disease"), ids)
table(phen_sub)

set.seed(42)
out <- treeWAS(
  snps  = snps_sub,
  phen  = phen_sub,
  tree  = sub_tree,
  test  = c("terminal", "simultaneous", "subsequent"),
  n.snps.sim = 10 * ncol(snps_sub),
  p.value = 0.05,
  p.value.correct = "fdr",
  snps.reconstruction = "parsimony",
  phen.reconstruction = "parsimony",
  plot.tree = FALSE
)

assoc <- c(row.names(out[c("terminal")]$terminal$sig.snps),
           row.names(out[c("simultaneous")]$simultaneous$sig.snps))

set_1_ags <- fread("./output/data/set_1_names.csv")

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "p_bh"))

gene_pvals <- gene_pvals[freq > 12]

gene_pvals[, is_outlier := fcase(p_bh <= 0.05, "sel",
                                default = "non_sel")]

tab <- table(outlier = gene_pvals$is_outlier,
             category_assoc = gene_pvals$gene_family %in% assoc)
fisher.test(tab)


# Raw fisher's test -------------------------------------------------------
# separate the ID column from the gene columns

pa <- as.matrix(gene_pa[,-1])          # ensure it's a matrix, not data.frame
row.names(pa) <- gene_pa$geno_id
pa <- t(pa)
pa[1:5,1:5]

hosp <- spark_meta$hospital_assoc[match(colnames(pa), spark_meta$id)] == "hospital"
stopifnot(length(hosp) == ncol(pa), !anyNA(hosp))


res <- apply(pa, 1, function(g) {
  tab <- table(factor(g, c(0,1)), hosp)
  ft <- fisher.test(tab)
  
  c(
    p = ft$p.value,
    OR = unname(ft$estimate),
    prev_hosp = mean(g[hosp]),
    prev_nonhosp = mean(g[!hosp])
  )
})

res <- as.data.table(t(res), keep.rownames = "gene_family")
res[, bh_p := p.adjust(p, "BH")]

res[, hos_grp:= fcase(bh_p < 0.05 & OR > 1, "hosp",
                      bh_p < 0.05 & OR < 1, "non-hosp",
                      default = "NS")]

raw_ecol <- merge(gene_pvals, res,
             all.x = TRUE, by ="gene_family")

raw_ecol[hos_grp!="NS", .(n = .N), hos_grp]

# hos_grp     n
# 1: non-hosp   132
# 2:     hosp   256

raw_ecol_tab <- table(raw_ecol[hos_grp!="NS", .(is_outlier, hos_grp)])

fisher.test(raw_ecol_tab)

