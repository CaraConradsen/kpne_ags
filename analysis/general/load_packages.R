load_packages <- function(packages){
  # takes a vector of package names
  
  # Install BiocManager if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  # Check for missing packages and install if needed
  missing_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if(length(missing_packages)) BiocManager::install(missing_packages)
  
  # Check and install gUtils from GitHub if not installed
  if (!requireNamespace("gUtils", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("mskilab/gUtils")
  }
  
  # Load the packages
  sapply(packages, require, character.only = TRUE)
}