# The piN / H / piS grid, with H residualised on gene length
#
# H counts protein haplotypes, so a longer gene reaches a higher H at the same
# genealogical depth; piN is per-site and carries no such dependence. Putting
# the two side by side therefore compares a length-sensitive measure with a
# length-free one. This replots piN_H_piS.R with H replaced by its residual
# from a loess on log gene length, fitted over all the AGs in the panel set.
#
# Everything else is unchanged, and the output is written alongside the raw
# version as piN_H_piS_resid.png.

use_resid <- TRUE

source("./analysis/14_piN_estimates_seg_sites/piN_H_piS.R")

use_resid <- FALSE   # leave the session as it was found
