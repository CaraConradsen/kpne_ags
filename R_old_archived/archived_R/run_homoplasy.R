args <- commandArgs(trailingOnly = TRUE)

presence_file <- args[1]
tree_file     <- args[2]
out_dir       <- args[3]

options(java.parameters = "-Xmx55g")
library(homoplasyFinder)

runHomoplasyFinderInJava(
  treeFile = tree_file,
  presenceAbsenceFile = presence_file,
  path = out_dir,
  includeConsistentSitesInReport = TRUE,
  multithread = TRUE
)
