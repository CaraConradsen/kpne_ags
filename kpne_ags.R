# ########################################################## ####
# Selection on AG in KLEBSIELLA                              ####
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
                            pattern =".R",
                             recursive = T)) 
invisible(sapply(
  AG_analysis, 
  source)
)

# library(remotes)
# remotes::install_version("rlang", version = "1.1.2")
# remotes::install_version("dplyr", version = "1.1.3")
# Package dependencies
# this function attaches libraries, but also checks and installs missing packages
packages <- c("ape","Biostrings","BSgenome", "doParallel","data.table", "MASS", "yyjsonr",
              "eulerr","foreach","GenomicRanges","ggplot2","ggtree","gUtils","ggnewscale", "RCy3",
              "mlplasmids","msaR","shape", "ade4", "genoPlotR","MASS","phangorn",
              "ggiraph","htmlwidgets","igraph", "jsonlite","pegas", "phyloseq", "micropan","VennDiagram", 
              "phytools","RColorBrewer","rtracklayer", "seqinr","stats","stringr","treeio", "vegan",
              "vioplot", "homoplasyFinder")

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

