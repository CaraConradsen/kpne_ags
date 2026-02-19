# Homoplasy within junction AGs ----------------------------------------
# create dir
outdir_homoplasy <- "./output/homoplasyFinder"

if (file.exists(outdir_homoplasy)==FALSE){
  dir.create(outdir_homoplasy, recursive = TRUE)
}

core_tree <- read.tree("./input_data/bootstrapped_gubbins/tmp8yl_0c9w/RAxML_bestTree.core_genome_aln.iteration_20")

core_tree_resolved <- multi2di(core_tree)

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


# Get bootstrapped homoplasy estimates ------------------------------------

# split each ag
ag_dat <- ag_presence_absence[, -3, with = FALSE]

fwrite(ag_dat,file.path(outdir_homoplasy, paste0("/presenceAbsence.csv")))


# # set into groups
# ag_dat[, group := (seq_len(.N) - 1) %/% 35 + 1]
# 
# ag_dat[, {
# 
#   ag_dir <- file.path(outdir_homoplasy, paste0("ag_", group))
#   dir.create(ag_dir, recursive = TRUE, showWarnings = FALSE)
# 
# 
#   NULL
# 
# }, by = group]

write.tree(core_tree_resolved, file = paste0(outdir_homoplasy, "/core_tree_resolved.tree"))


# Run homoplasyFinder -----------------------------------------------------

# Find the csv and tree files
presenceAbsenceFile <- list.files(outdir_homoplasy,
                                  recursive = TRUE, full.names = TRUE,
                                  pattern = ".csv")

treeFile <- list.files(outdir_homoplasy,
                       recursive = TRUE, full.names = TRUE,
                       pattern = ".tree")

system2("Rscript",
        args = c("./analysis/06_id_HGT_ags/run_homoplasy.R",
                 presenceAbsenceFile,
                 treeFile,
                 paste0(outdir_homoplasy,"/")))


# Import Consistency indices ---------------------------------------------
# get ci
list_homoplasies <- list.files(outdir_homoplasy,
                               recursive = TRUE, full.names = TRUE,
                               pattern = ".txt")

CI_est <- fread(list_homoplasies)

# Plot relationship -------------------------------------------------------

colnames(CI_est) <- tolower(colnames(CI_est))

CI_est <- merge(CI_est,
                ag_presence_absence[,.(start, end, gene_family)],
                all.x = TRUE, by = c("start", "end"))

# save data 
fwrite(CI_est[,-c(1,2)], paste0(outdir_dat, "/consistencyindex.csv"))


grping_dat <- unique(ag_homo_dat[,.(gene_family, anchor, number_genomes)])

grping_dat[, syn_jun := fcase(anchor == "", 0,
                              default = 1)]

grping_dat[, anchor := NULL]

plot_data <- merge(CI_est,grping_dat,
                   all.x=TRUE, by = c("gene_family") )


png(paste0(outdir_fig,"/ci_vs_pangenome.png"),
    width = 11, height = 10.5, units = "cm", res = 300,
    pointsize = 10, type = "cairo")


par(mar = c(4, 4,2,0.5))
with(plot_data,
     plot(number_genomes, consistencyindex ,
          pch = 19, bty = "L",
          cex = 0.5, yaxt ="n",
          col = ifelse(syn_jun==0,
                       rgb(0.1,0.1, 0.8, alpha = 0.5),
                       rgb(0.8,0.1, 0.1, alpha = 0.5)),
          xlab = "Number of genomes", 
          ylab = "Mean consistency index"))

axis(side=2, at = seq(0,1, 0.2),
     las = 2,
     labels = sprintf("%.1f", seq(0,1, 0.2)))


legend("top", pch = c(19,19), 
       col = c(rgb(0.8,0.1, 0.1, alpha = 0.5),
               rgb(0.1,0.1, 0.8, alpha = 0.5)
       ), cex = 0.75, bty="n",
       legend = c("Within syntenic block", "Unknown synteny"))

dev.off()
