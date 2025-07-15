# ########################################################## ####
# AG introductions in KLEBSIELLA                              ####
# Author:    Cara Conradsen                                  ####
# ########################################################## ####

# Script purpose ####
# To investigate accessory genes (AGs) in Klebsiella pneumonia

# Environment set-up -----------------------------------------------------------------
# Load in dependencies and create subdirectories

# import vcf R scripts from the analysis subdirectory
invisible(# removes verbose warnings
  AG_analysis <- list.files("./analysis/general", 
                             full.names = T,
                            pattern ="load_packages.R",
                             recursive = T)) 
invisible(sapply(
  AG_analysis, 
  source)
)


# Package dependencies
# this function attaches libraries, but also checks and installs missing packages
packages <- c("ape","Biostrings","BSgenome", "data.table","doParallel","dplyr","MASS", "jsonlite",
              "eulerr","foreach","GenomicRanges","ggplot2","ggtree","gUtils","ggnewscale",
              "ggiraph","htmlwidgets","igraph", "jsonlite","pegas", "phyloseq", "micropan","VennDiagram", 
              "phytools","RColorBrewer","rtracklayer", "seqinr","stats","stringr","treeio", "vegan")

# load in functions, uses BiocManager::install which handles both 
# CRAN and Bioconductor packages
load_packages(packages)

# set data.table threads
setDTthreads(threads = detectCores()-2)

# Set up parallel backend 
num_cores <- detectCores()-2

# create output directories
outdir_tab <- "./output/tables"
outdir_fig <- "./output/figures"
outdir_dat <- "./output/data"

if (file.exists(outdir_tab)==FALSE){
  dir.create(outdir_tab, recursive = TRUE)
}
if (file.exists(outdir_fig)==FALSE){
  dir.create(outdir_fig, recursive = TRUE)
}
if (file.exists(outdir_dat)==FALSE){
  dir.create(outdir_dat, recursive = TRUE)
}

