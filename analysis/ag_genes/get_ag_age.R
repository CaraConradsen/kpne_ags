# Get AG Age


# input require data sets
# ags
ag_dt <- fread(paste0(outdir_dat,"/ag_lng_mapped.csv"))

# pirate aln
pirate_pangenome_aln = fread(paste0(outdir_dat,"/pirate_pangenome_aln.csv"))

colnames(pirate_pangenome_aln)[2] = "loci_id"

ag_dt <- merge(ag_dt, pirate_pangenome_aln, all.x = TRUE, by="loci_id")

list_unique_genomes = unique(ag_dt[!is.na(gene_family_aln), geno_id])

list_unique_ags = unique(ag_dt[!is.na(gene_family_aln), gene_family]) # 11,900 out of 11,992


# loop
# directory of algined genes
nt_aln_dir = "./input_data/PIRATE_1695_out/feature_sequences/"


# # start timer
# start.time <- Sys.time() 
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)

ag_age_dt <- foreach(i = list_unique_ags,
                     .packages = c("data.table", "seqinr", "Biostrings"),
                     .combine = rbind) %do% {
                       
                       
                       
                       n100_smpl_geno = sample(list_unique_genomes,1000)
                       
                       n100_smpl_loci_names = ag_dt[gene_family == i & geno_id %chin% n100_smpl_geno, loci_id]
                       
                       loc_freq = length(n100_smpl_loci_names)
                       
                       if(loc_freq >1){
                         focal_ag_aln <- Biostrings::readDNAStringSet(paste0(nt_aln_dir,i,".nucleotide.fasta"))
                         focal_ag_aln <- focal_ag_aln[names(focal_ag_aln) %in% n100_smpl_loci_names]
                         
                         if (any(grepl("-", as.character(focal_ag_aln)))==FALSE){
                         
                         # Convert to "alignment" class
                         alignment_obj <- list(
                           nb = length(focal_ag_aln),
                           nam = names(focal_ag_aln),
                           seq = as.character(focal_ag_aln),  # must be character vector
                           com = NULL
                         )
                         class(alignment_obj) <- "alignment"
                         
                         class(alignment_obj) <- "alignment"
                         
                         kaks_result <- try(seqinr::kaks(alignment_obj), silent = TRUE)
                         
                         if (inherits(kaks_result, "try-error") || is.null(kaks_result$ks)) {
                           max_ks <- NA
                         } else {
                           max_ks <- max(kaks_result$ks, na.rm = TRUE)
                         }
                         
                         }else{
                           max_ks = NA
                         }
                       } else{
                         max_ks = 0
                       }
                       
                       data.frame(gene_family = i, max_ks = max_ks, ag_count = loc_freq)
                     }

# # Stop the parallel backend after the loop is done
# stopCluster(cl)
# 
# # get duration
# end.time <- Sys.time()
# end.time - start.time# Time difference of 

fwrite(ag_age_dt, paste0(outdir_dat, "/ag_age_dt.csv"))
# ag_age_dt <-fread(paste0(outdir_dat, "/ag_age_dt.csv"))

# plot data
ag_age_dt = as.data.table(ag_age_dt)

# prop of non NAs
non_na = nrow(ag_age_dt[!is.na(max_ks)])

# set axes
pretty_x = c(0, 0.4)
pretty_y = pretty(range(ag_age_dt$ag_count/1000, na.rm = TRUE))

mat <- matrix(c(1,2,3,4), byrow = TRUE, ncol=2)


# start figure

png(filename = paste0(outdir_fig,"/Ag_age_vs_freq_2.png"),
    pointsize = 20,
    width = 9, height = 9, units = "in", res = 300)

# layout(mat, widths = c(2.25,0.75),
#        heights = c(0.75, 2.25))
# 
# par(mar = c(0.25,4,0.5,0.5))
# with(ag_age_dt[!is.na(max_ks),],
#      hist(max_ks,
#           xaxt = "n",
#           las=2,
#           border = "darkblue",
#           ylab = "count",
#           col = "royalblue1",
#           breaks = 80,
#           main="",
#           xlim = range(pretty_x)))
# 
# plot.new()

par(mar = c(4,4,0.25,0.5))
with(ag_age_dt[!is.na(max_ks),],
     plot( ag_count/1000,max_ks,
          pch=16, 
          las = 2,
          xaxt = "n",
          ylab = "AG frequency in 1000 genomes",
          xlab = "maximum dS",
          xlim = range(pretty_x),
          ylim = range(pretty_y),
          col = rgb(0.1,0.1,0.9, alpha = 0.5)))
axis(side = 1, at = seq(0,1,0.1), labels = seq(0,1,0.1))
# legend("topright", bty='n',
#        paste0(non_na, " non-NA AGs out of ", nrow(ag_age_dt)))
# 
# # barplot version
# par(mar = c(4,0.25,0.25,0.5))
# hist_ag_count <- with(ag_age_dt[!is.na(max_ks),],
#                       hist(ag_count/1000,
#                            breaks = 80, 
#                            plot = FALSE))
# 
# # Plot as a horizontal barplot
# barplot(hist_ag_count$counts,
#         horiz = TRUE,
#         yaxt = "n",
#         xlab = "count",
#         space = 0,
#         main= "",
#         col = "royalblue1",
#         names.arg = round(hist_ag_count$mids, 1))

dev.off()
