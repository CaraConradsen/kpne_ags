

outlier_ags <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                     select = c("freq", "gene_family", "p_bh", 
                                "pi_outlier", "pi_sel", "set"))[set==1]

outlier_ags <- outlier_ags[p_bh < 0.05]

# AG presence / absence ---------------------------------------------------

# Get remaining accessory genes after synteny analysis --------------------
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

ag_dat <- msu_regions_anchored[ag_type != "core"][number_genomes != 1]#[!grepl("_",gene_family)][acrs_msu!=1][acrs_jun != 1]

set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))

ag_dat <- ag_dat[gene_family %chin% set_1_ags$gene_family]

# create presence/absence for gene_families -------------------------------
# Create binary ag presence_absence data.table 
# needs a dummy start and end, where rows are each gene family and columns are genomes

# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_dat[,.(geno_id, gene_family)],
  formula = geno_id ~ gene_family,
  fun.aggregate = length,
  value.var = "gene_family"
)

outlier_genos <- lapply(outlier_ags$gene_family, function(ag){
  data.table(gene_family = ag, 
             geno_id = ag_presence_absence[get(ag) == 1, geno_id])
})

outlier_genos <- rbindlist(outlier_genos)

outlier_genos_sum <- outlier_genos[,.(n =.N), geno_id]

# get total ags
gene_counts <- t(ag_presence_absence[,-1])
colnames(gene_counts) <- ag_presence_absence$geno_id

gene_counts <- colSums(gene_counts)

gene_counts_dt <- data.table(geno_id = names(gene_counts),
                             tot_ags = gene_counts)


outlier_genos_sum <- merge(outlier_genos_sum, gene_counts_dt,
                           by = "geno_id")

outlier_genos_sum[,perc := round((n/tot_ags)*100, digits = 2)]


# ST colours --------------------------------------------------------------
pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))

ST_groups <- unique(pan_anno[,.(geno_id, ST)])

# get colours
ST_cols = ST_groups[,.(n=.N), ST][n>1, .(ST)]
setorderv(ST_cols, "ST",1)

ST_cols[, cols := colorRampPalette(c("#420D55","blueviolet",
                                     "darkblue","#008080",
                                     "forestgreen","#799000",
                                     "yellow","#FF8000",
                                     "brown"))(nrow(ST_cols))]

ST_groups <- merge(ST_groups, ST_cols,
                   all.x = TRUE, by ="ST")
ST_groups[is.na(cols), cols := "black" ]

outlier_genos_sum <- merge(outlier_genos_sum, 
                           ST_groups, by ="geno_id")

# Genetic distance --------------------------------------------------------
core_tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")


D <- cophenetic(core_tree)          # or dist.dna on the alignment
bg_divergence <- rowMeans(D)        # per-genome background metric

bg_divergence <- data.table(geno_id = names(bg_divergence),
                            bg_divergence = bg_divergence)


outlier_genos_sum <- merge(outlier_genos_sum, bg_divergence,
                           by = "geno_id")

with(outlier_genos_sum,
     plot(perc, bg_divergence,
          bty = "L",
          xlim = c(0,20),
          pch = 16,
          col = cols,
          ylab = "Mean pairwise core distance",
          xlab = "Significant high-diversity AGs (% per genome)"))


# Burden vs accessory-genome size, coloured by source. 
# Directly tests your hospital-vs-water idea. x = number of tested AGs
# per genome (or total AG count), y = outlier %, colour = ecological 
# source from metadata. If outlier % is flat against accessory 
# size once you see the colours cluster, the percentage is 
# just tracking repertoire size/ecology.

with(outlier_genos_sum,
     plot(tot_ags, perc,
          bty = "L",
          ylim = c(0,20),
          pch = 16,
          col = cols,
          xlab = "Total AGs",
          ylab = "Significant high-diversity AGs (% per genome)"))



# Heatmaps ----------------------------------------------------------------


# Log  compression
core_tree$edge.length <- log1p(core_tree$edge.length)

tree_ultra <- phytools::force.ultrametric(core_tree) 

hc <- as.hclust(tree_ultra)
# Convert to dendrogram
dend <- as.dendrogram(hc)

# Cluster heat map 
temp_mat <- ag_presence_absence[, c("geno_id", outlier_ags$gene_family), with = FALSE]

temp_mat[, geno_id := factor(geno_id, levels = core_tree$tip.label)]
setorder(temp_mat, geno_id)

mat <- as.matrix(temp_mat[,-1])
row.names(mat) <- temp_mat$geno_id

# set colour conditions
breaks <- c(-0.5, 0.5, 1.5)
cols   <- c("white", "#4DB3FF")

heatmap(mat, Colv = TRUE, Rowv = dend,
        scale = "none",
        labRow = NA, 
        col    = cols,        
        cexCol = 1, 
        main = "", 
        margins = c(4, 2),
        breaks = breaks)

  

col_dist <- dist(t(mat), method = "binary")
col_dend <- as.dendrogram(hclust(col_dist))

heatmap(mat, Colv = col_dend, Rowv = dend,
          scale = "none",
          labRow = NA,
          col = cols,
          breaks = breaks,
          trace = "none",
          dendrogram = "column",
          key = FALSE)


# The heatmap partitions the 423 into the two categories 
# Adam's argument demands: phylogeny-congruent (most of them — 
# background, not attributable) and phylogeny-incongruent (the 
# minority — where attribution is even possible). This is the 
# presence-level version of the decoupling test; it doesn't 
# yet separate recurrent HGT from balancing selection among the 
# incongruent set, but it isolates which genes 
# are eligible for that question.
