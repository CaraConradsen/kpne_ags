# Will need:
#   1. the recombination free core genome phylogeny
#   2. a binary presence / absence file to recode

# Check convert gubbins final tre file to newick is rooted -----------------------------
core_gub_tree <- read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre") 

is.rooted(core_gub_tree)

# Create binary ag presence_absence data.table------------------------------------------------------------
# needs a dummy start and end, where rows are each gene family and columns are genomes

ag_info <- fread(paste0(outdir_dat,"/pangraph_anno.csv"),
                 select = c("geno_id","gene_family", "ag_type"))

# id orfans
ag_info[, n := .N, gene_family]

# remove paralogs (retain genes with only _1 or no underscore)

keep_ags = unique(c(grep("_", ag_info$gene_family, invert = TRUE, value = TRUE),
                    grep("_1", ag_info$gene_family, value = TRUE)))

ag_info <- ag_info[gene_family %chin% keep_ags]


# Generate a binary matrix giving gene presence/absence
ag_presence_absence <- dcast(
  data = ag_info[n!=1 & ag_type !="core",.(geno_id, gene_family)],
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

# save dummy regions key
fwrite(ag_presence_absence[, .(start, end, gene_family)],
       file = paste0(outdir_dat, "/dummy_regions_gene_family_key.csv"))

ag_presence_absence[, gene_family := NULL]

# export csv for homoplasy
fwrite(ag_presence_absence,
       file = paste0(outdir_dat, "/ag_presence_absence.csv"))

# Run homoplasyfinder -----------------------------------------------------
# create directory
outdir_homoplasy <- "./output/homoplasy"

if (file.exists(outdir_homoplasy)==FALSE){
  dir.create(outdir_homoplasy, recursive = TRUE)
}


# Find the FASTA and tree files attached to package
presenceAbsenceFile <- list.files(outdir_dat, full.names = TRUE, pattern = "ag_presence_absence.csv")
treeFile <- "./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre"

# Set working directory
workingDirectory <- paste0(normalizePath(outdir_homoplasy, winslash = "/"), "/")

# Run the HomoplasyFinder jar tool
inconsistentPositions <- runHomoplasyFinderInJava(treeFile=treeFile, 
                                                  presenceAbsenceFile=presenceAbsenceFile, 
                                                  path=workingDirectory,
                                                  createFasta=FALSE)

# Read in the output table
resultsFile <- list.files(workingDirectory, full.names = TRUE, pattern =  ".txt")
CI_report <- fread(resultsFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)


# Read in the annotated tree
anno_tree <- read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre") 

anno_tree <- read.tree(list.files(workingDirectory, full.names = TRUE, pattern =  ".tree"))

# Plot the annotated tree
plotAnnotatedTree(tree, inconsistentPositions, fastaFile)



