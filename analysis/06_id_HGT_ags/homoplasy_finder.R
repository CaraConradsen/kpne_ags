# Homoplasy within junction AGs ----------------------------------------
# create dir
outdir_homoplasy <- "./output/homoplasyFinder"

if (file.exists(outdir_homoplasy)==FALSE){
  dir.create(outdir_homoplasy, recursive = TRUE)
}

# get bootstrapped trees
boot_trees <- read.tree("./input_data/bootstrapped_gubbins/tmp8yl_0c9w/RAxML_bootstrap.core_genome_aln.iteration_20.bootstrapped_trees")

# Get remaining accessory genes after synteny analysis --------------------
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

ag_homo_dat <- msu_regions_anchored[acrs_msu!=1][acrs_jun != 1][ag_type != "core"][number_genomes != 1][!grepl("_",gene_family)]

# create presence/absence for gene_families -------------------------------
# Create binary ag presence_absence data.table 
# needs a dummy start and end, where rows are each gene family and columns are genomes

# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_homo_dat[,.(geno_id, gene_family)],
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

# save data 
fwrite(ag_presence_absence, paste0(outdir_dat, "/homoplasy_ag_presence_absence.csv"))


# Get bootstrapped homoplasy estiamtes ------------------------------------

# split each ag
ag_dat <- ag_presence_absence[, -3, with = FALSE][1:10,]# CHANGE THIS!!

ag_dat[, {
  
  ag_dir <- file.path(outdir_homoplasy, paste0("ag_", start))
  dir.create(ag_dir, recursive = TRUE, showWarnings = FALSE)
  
  fwrite(.SD,
         file.path(ag_dir, paste0("presenceAbsence_INDELs_", start, ".csv")))
  
  NULL
  
}, by = start]


# write trees
tree_rep_dir <- paste0(outdir_homoplasy,"/trees")

dir.create(tree_rep_dir, recursive = TRUE, showWarnings = FALSE)

for (i in 1:10) {#length(boot_trees) CHANGE THIS!!
  write.tree(boot_trees[[i]], file = paste0(tree_rep_dir, "/boot_tree_",i,".tre"))
}

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
