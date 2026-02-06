###############################################################################
#  Klebsiella Genome Feature Integration
# 
# Description:
#   This script processes and integrates genomic features across Klebsiella
#   chromosomes, including:
#     - All chromosomal genes (core and accessory)
#     - Mobile genetic elements (MGEs) and functional categories
#     - Pangraph syntenic blocks (MSUs)
#
# Input:
#   - Genome annotation tables
#   - Core/accessory gene definitions
#   - Pangraph syntenic blocks
#   - MGE/recombination annotations
#
# Output:
#   - Integrated data tables linking genes, syntenic blocks, MGEs, and recombination
#   - Tables ready for visualization and comparative analysis
#
# Settings for core
# fcase(n < 0.99 * n_genomes & n >=  0.95 * n_genomes, "soft",
#       n < 0.95 * n_genomes & n >=  0.15 * n_genomes, "shell",
#       n < 0.15 * n_genomes, "cloud",
#       default = "core")]
###############################################################################
# Retain only correct oriented chromosomes --------------------------------
pan_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))[,c(1:11, 13,17)]


# subset anno by pangraph anno
pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"),
                       select = c("locus_tag", "ag_type"))# 214 genomes


pan_anno <- merge(pangraph_anno, pan_anno,
                  all.x = TRUE, by ="locus_tag")

# add syntenic identifiers
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                              select = c("locus_tag", "anchor", "msu"))# accounts for 99.68% loci 

pan_anno <- merge(pan_anno,msu_regions_anchored,
                  all.x = TRUE, by ="locus_tag")

# tidy binary cats
pan_anno[is.na(acrs_msu), acrs_msu := NA]
pan_anno[is.na(acrs_jun), acrs_jun := NA]


# Fix paralogs ------------------------------------------------------------
# check to see if paralogs are no longer paralogs after subsetting
paralogs <- unique(pan_anno[grepl("_",gene_family), .(gene_family)])

paralogs[, broad_fam := tstrsplit(gene_family, "_", fill = TRUE, keep = 1)]

paralogs[, n:=.N, broad_fam]

non_paralogs = paralogs[n<2, gene_family]

paralogs <- paralogs[n>1, gene_family]

# # update data
# pan_anno[gene_family %chin% non_paralogs,
#          gene_family := tstrsplit(gene_family, "_", fill = TRUE, keep = 1)]

pan_anno$paralog = 0

pan_anno[gene_family %chin% paralogs, paralog := 1]

pan_anno[, broad_fam := tstrsplit(gene_family, "_", fill = TRUE, keep = 1)]

# Add Eggnog Cogs ---------------------------------------------------------

eggnog_anno <- fread("./input_data/klebsiella_eggnog.emapper.annotations", fill=TRUE,skip = 4,
      select=c("#query","COG_category","Preferred_name"))

colnames(eggnog_anno)[1] ="gene_family"

eggnog_anno[, gene_family := tstrsplit(gene_family, ";", fill = TRUE, keep = 1)]

# Remove bottom lines
eggnog_anno <- eggnog_anno[grepl("g0", gene_family)]

# add to pan_anno
pan_anno <- merge(pan_anno, 
                  unique(eggnog_anno[,.(gene_family, COG_category)]),
                  all.x = TRUE, by = c("gene_family"))

pan_anno[is.na(COG_category), COG_category := "-"]

# Add MGEs ----------------------------------------------------------------

# Add in ISESCan IS



# Some counts -------------------------------------------------------------
unique(pan_anno[, broad_fam]) %>% length
unique(pan_anno[ag_type=="core", broad_fam]) %>% length
unique(pan_anno[ag_type!="core" & paralog == 1,broad_fam]) %>% length

# Some paralogs are still in the data, but family is core
AGs_post_cat = pan_anno[!broad_fam %chin% pan_anno[ag_type=="core", broad_fam]]
AGs_post_cat = AGs_post_cat[paralog == 0]

# remove singletons
AGs_post_cat[, n:=.N, gene_family]
AGs_post_cat <- AGs_post_cat[n!=1]

AGs_post_cat[,c("n", "broad_fam"):= NULL]


# Id across msu & junctions -----------------------------------------------

# id duplicates within and across msu/junctions ---------------------------

# Examine the msu blocks 

msu_loci <- AGs_post_cat[grepl(":", anchor)| (anchor=="" | is.na(anchor))]

# Are there any AGs across msu?
msu_loci_grp = unique(msu_loci[,.(ag_type, gene_family, msu)])

across_msu_loci = msu_loci_grp[,.(n = .N), by = c("ag_type", "gene_family")][n>1][,gene_family]


# code across junctions
AGs_post_cat$acrs_msu = 0
AGs_post_cat[gene_family %chin% across_msu_loci, acrs_msu := 1]

# Examine the syntenic junctions ---------------------------------------------------
# take into account null estimates

junction_loci <- AGs_post_cat[grepl(":", anchor)]

# Are there any AGs across junctions?
junction_loci_grp = unique(junction_loci[,.(ag_type, gene_family, msu, anchor)])

across_junction_loci = junction_loci_grp[,.(n = .N), by = c("ag_type", "gene_family", "msu")][n>1][,gene_family]


# code across junctions
AGs_post_cat$acrs_jun = 0
AGs_post_cat[gene_family %chin% across_junction_loci, acrs_jun := 1]

# remove non-syntenic blocks
AGs_post_cat <- AGs_post_cat[acrs_jun == 0 & acrs_msu == 0]

# Homoplasy within junction AGs ----------------------------------------
# Check convert gubbins final tre file to newick is rooted -----------------------------
# create dir
outdir_homoplasy <- "./output/homoplasyFinder"

if (file.exists(outdir_homoplasy)==FALSE){
  dir.create(outdir_homoplasy, recursive = TRUE)
}

core_gub_tree <- read.tree("./input_data/pangraph/gub_graph.node_labelled.final_tree.tre") 

is.rooted(core_gub_tree)

# save to homoplasy dir
write.tree(core_gub_tree, file = paste0(outdir_homoplasy, "/core_gub_tree.tre"))

# Create binary ag presence_absence data.table 
# needs a dummy start and end, where rows are each gene family and columns are genomes

# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = AGs_post_cat[,.(geno_id, gene_family)],
  formula = gene_family ~ geno_id,
  fun.aggregate = length,
  value.var = "gene_family"
)

# make dummy regions
dummy_regions <- data.frame(
  start = seq(1, by = 2, length.out = nrow(ag_presence_absence)),
  end   = seq(2, by = 2, length.out = nrow(ag_presence_absence))
)

ag_presence_absence <- cbind(dummy_regions, ag_presence_absence)

setDT(ag_presence_absence)

# save to homoplasy dir
fwrite(ag_presence_absence[,-3], paste0(outdir_homoplasy, "/presenceAbsence_INDELs.csv"))


# Run homoplasyFinder -----------------------------------------------------
old_dir <- getwd()
setwd(outdir_homoplasy)


# Find the csv and tree files
presenceAbsenceFile <- paste0(getwd(), "/presenceAbsence_INDELs.csv") 
treeFile <- paste0(getwd(), "/core_gub_tree.tre")


# Get the current working directory
workingDirectory <- paste0(getwd(), "/")

# Run the HomoplasyFinder jar tool
inconsistentPositions <- runHomoplasyFinderInJava(treeFile=treeFile, 
                                                  presenceAbsenceFile=presenceAbsenceFile, 
                                                  path=workingDirectory)
# clean file names
files <- list.files()
new_files <- gsub("_\\d{2}-\\d{2}-\\d{2}", "", files)
file.rename(from = files, to = new_files)

# Read in the output table
resultsFile <- paste0(workingDirectory, "consistencyIndexReport.txt")
homoplasy_res <- read.table(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)

setDT(homoplasy_res)

setwd(old_dir)

png(
  filename = paste0(outdir_fig, "/Homoplasy_214.png"),
  width = 800,  # total width in pixels
  height = 900, # total height in pixels
  res = 150      # resolution, 150–300 dpi
)

par(mar = c(4,4,0.1,0.1))
hist(homoplasy_res$MinimumNumberChangesOnTree, 
     breaks = 41, xlim = c(1,41), 
     main = "",
     yaxt = "n", xlab = "Minimum number of changes on Tree",
     col = "dodgerblue", border="dodgerblue4")
axis(side = 2, at = seq(0, 2000, length.out = 5),
     las =2, labels = seq(0, 2000, length.out = 5))

dev.off()

png(
  filename = paste0(outdir_fig, "/Homoplasy_CI_214.png"),
  width = 800,  # total width in pixels
  height = 900, # total height in pixels
  res = 150      # resolution, 150–300 dpi
)
par(mar = c(4,4,0.1,0.1))
hist(homoplasy_res$ConsistencyIndex, 
     breaks = 41, xlim = c(0,0.5), 
     main = "",
     yaxt = "n", xlab = "Consistency index",
     col = "firebrick", border="firebrick4")
axis(side = 2, at = seq(0, 1500, length.out = 4),
     las =2, labels = seq(0, 1500, length.out = 4))

dev.off()

homoplasy_res[, multi_gain:= fcase(ConsistencyIndex < 0.5 | MinimumNumberChangesOnTree > 3, 1,
                                   default = 0)]
colnames(homoplasy_res)[1:2] = c("start", "end")

homoplasy_res <- merge(ag_presence_absence[,1:3], homoplasy_res, 
                       all.x = TRUE, by = c("start", "end"))

AGs_post_cat <- merge(AGs_post_cat, homoplasy_res[, .(gene_family,ConsistencyIndex, MinimumNumberChangesOnTree, multi_gain)],
                      all.x = TRUE, by = "gene_family")

AGs_post_cat[is.na(multi_gain), multi_gain:= 0]

fwrite(AGs_post_cat, paste0(outdir_dat, "/AGs_post_cat.csv"))

# Ancestral reconstructions -----------------------------------------------

# Example: your binary presence/absence matrix
# PA: rows = genomes, columns = genes
# core_tree: recombination-corrected core genome tree

PA = t(ag_presence_absence[,4:ncol(ag_presence_absence)])
colnames(PA) = ag_presence_absence$gene_family


ag_gains_list <- foreach(ag_num = 1:ncol(PA),
                        .packages = c("ape", "phangorn")) %do% {
  # Extract gene column (tip states)
  gene_idx <- ag_num
  
  # gene_fam
  gene_family = colnames(PA)[ag_num]
  
  # Convert to phyDat for this gene
  gene_phy <- phyDat(t(t(PA[, gene_idx])), type="USER", levels=c(0,1))
  
  # ACCTRAN reconstruction
  gene_anc <- ancestral.pars(core_gub_tree, gene_phy, type="ACCTRAN")
  
  # The result is a list of tip matrices
  # Use first element to extract ancestral states
  node_names = c(row.names(PA), paste0("Node_", 1: Nnode(core_gub_tree)))
  anc_matrix <- t(sapply(node_names, function(node) gene_anc[[node]]))  # rows = nodes, columns = 0/1
  
  # Count 0→1 gains along edges
  edges <- core_gub_tree$edge
  gains <- 0
  for (i in 1:nrow(edges)) {
    parent <- edges[i,1]
    child  <- edges[i,2]
    
    parent_state <- which.max(anc_matrix[parent,]) - 1
    child_state  <- which.max(anc_matrix[child,]) - 1
    
    if (parent_state == 0 & child_state == 1) gains <- gains + 1
  }
  
  # return
  list(gene_family = gene_family, gains = gains)
}

ag_gains <- rbindlist(ag_gains_list)

png(
  filename = paste0(outdir_fig, "/Ancestral_parsimony_214.png"),
  width = 800,  # total width in pixels
  height = 900, # total height in pixels
  res = 150      # resolution, 150–300 dpi
)
par(mar = c(4,4,0.1,0.1))
hist(ag_gains$gains, breaks = 80, 
     col = "forestgreen", 
     yaxt = "n", main = "",
     xlab = "Number of gains (parent < child)")
axis(side = 2, at = seq(0, 800, 200),
     las =2, labels = seq(0, 800, 200))
text(80, 820, 'phangorn::ancestral.pars(type="ACCTRAN")', pos = 2)

dev.off()

AGs_post_cat <- merge(AGs_post_cat, ag_gains,
                      all.x = TRUE, by = "gene_family")

AGs_post_cat[, phy_multi:= fcase(multi_gain == 1 & gains > 1,1,
                                   default = 0)]



# COG plot for HGTs -------------------------------------------------------
pan_anno


