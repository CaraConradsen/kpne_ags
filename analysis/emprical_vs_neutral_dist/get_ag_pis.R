# Get AG piS (age)

# input require data sets
# ags

pangraph_anno <- fread(paste0(outdir_dat,"/pangraph_anno.csv"))

list_unique_ags = unique(pangraph_anno[ag_type!="core", gene_family]) 

# remove paralogs (retain genes with only _1 or no underscore)

list_unique_ags = c(grep("_", list_unique_ags, invert = TRUE, value = TRUE),
                    grep("_1", list_unique_ags, value = TRUE))# 14,330 ags

# numer of genomes in pangenome
tot_pangenome_size = length(unique(pangraph_anno$geno_id))

# loop
# directory of algined genes
gene_align_loc = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_485_lng_rds_out/feature_sequences/"


# start timer
start.time <- Sys.time()

# Set up parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallel foreach loop
ag_age_dt <- foreach(i = list_unique_ags, 
                        .combine = rbind, 
                        .packages = c("Biostrings","ape","pegas","data.table")) %dopar% {
                          
                          # Read the FASTA, skip if missing
                          temp_string <- tryCatch(
                            readDNAStringSet(paste0(gene_align_loc, i, ".nucleotide.fasta")),
                            error = function(e) return(NULL)
                          )
                          
                          # skip if file is missing
                          if (is.null(temp_string)) return(NULL)
                          
                          #  Subset to loci present in pangraph annotation
                          sampled_loci = pangraph_anno[gene_family == i, locus_tag]
                          
                          # Subset by sequence names
                          temp_string <- temp_string[names(temp_string) %in% sampled_loci]
                          
                          # Skip if no sequence, ORFan or core gene
                          if (length(temp_string) == 0) return(NULL)
                          if (length(temp_string) %in% c(1, tot_pangenome_size)) return(NULL)
                          
                          # save k strains
                          loc_freq = length(temp_string)
                          
                          # Convert to "alignment" class
                          alignment_obj <- list(
                            nb = length(temp_string),
                            nam = names(temp_string),
                            seq = as.character(temp_string),  # must be character vector
                            com = NULL
                          )
                          class(alignment_obj) <- "alignment"
                          
                          kaks_result <- try(seqinr::kaks(alignment_obj), silent = TRUE)
                          
                          if (inherits(kaks_result, "try-error") || is.null(kaks_result$ks)) {
                            max_ks <- NA
                          } else {
                            max_ks <- max(kaks_result$ks, na.rm = TRUE)
                          }
                          
                          #return value
                          data.frame(gene_family = i, max_ks = max_ks, freq = loc_freq)
                          

                        }


# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 49.81543 secs

# Stop cluster when done
stopCluster(cl)

setDT(ag_age_dt)

fwrite(ag_age_dt, paste0(outdir_dat, "/ag_age_dt.csv"))
# ag_age_dt <-fread(paste0(outdir_dat, "/ag_age_dt.csv"))



# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
neutral_summary <- fread(paste0(outdir_dat, "/neutral_summary.csv"))

# Join neutral quantiles and define outliers
ag_age_dt[neutral_summary, on = "freq",
          outlier := fifelse(
            max_ks < pi_lower | max_ks > pi_upper, 1L, 0L
          )
]



# Plot graph --------------------------------------------------------------

pdf(paste0(outdir_fig, "/emp_vs_neut_outliers.pdf"), width = 8.27, height = 5.845)

par(mar = c(4,4,1,1))
with(ag_age_dt, 
     plot(freq, max_ks, pch = 16,
          ylim = c(0, 2.5),
          cex = 0.85,
          yaxt = "n", bty="L",
          ylab = expression(pi[S]),
          xlab = "Frequency in the pangenome",
          col = ifelse(outlier==1, rgb(0.8,0.1,0.1, alpha = 0.75),
                       rgb(0.1,0.1,0.8, alpha = 0.75))))

axis(side = 2, at = seq(0,2.5,0.5),
     labels = sprintf("%1.2f", seq(0,2.5,0.5)), las = 2)

legend("top", horiz = TRUE,
       legend = c("Outliers", "Non-outlier"),
       col = c("red", "blue"),
       pch = 16,
       bty = "n")

dev.off()

# Deviation values ------------------------------------------------------
# Join neutral quantiles and define outliers
ag_age_dt[neutral_summary, on = "freq",
          D := (max_ks - pi_median)/iqr
]


# add cogs
cogs <- pangraph_anno[COG_funct_cat!="", .(
  COG_mode = names(which.max(table(COG_funct_cat)))
), by = gene_family]


ag_age_dt[cogs, on = "gene_family",
          COG_letter := COG_mode
]

ag_age_dt <- merge(ag_age_dt,COG_dt,
                   all.x = TRUE,
                   by = "COG_letter")


cog_bxplt_dt = ag_age_dt[!is.na(COG_num) & D <= 5 & D >= -5]


pdf(paste0(outdir_fig, "/cog_D.pdf"), width = 8.27, height = 11.69)

par(mar = c(7,14,1,1))                   
plot(NULL,
     bty = "L",yaxt = "n",
     ylab = "",
     xlab = "Robust deviation score",
     ylim = c(26, 1),
     xlim = c(-5, 5))

with(COG_dt,
     axis(side = 2, at = COG_num, cex.axis = 0.5,
          col = "grey30",
          labels = paste0(COG_letter,", ",COG_function), las = 2))

abline(v = 0, col = "grey30", lwd = 1)
abline(v = 1.5, col = "grey30", lty = 2)
abline(v = -1.5, col = "grey30", lty = 2)

for (i in 1:26) {
  if(nrow(cog_bxplt_dt[COG_num == i]) >= 10){
    with(cog_bxplt_dt[COG_num == i],
         vioplot(D, at = i, 
                 col = "salmon1",
                 border = NA,
                 horizontal = TRUE,
                 drawRect = FALSE,
                 wex = 0.75,
                 add = TRUE))
  }
}

for (i in 1:26) {
  with(cog_bxplt_dt[COG_num == i],
       points(D, COG_num, factor = 0.5, pch = 16, cex = 0.65,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  if(nrow(cog_bxplt_dt[COG_num == i]) >= 10){
    with(cog_bxplt_dt[COG_num == i],
         points(median(D), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}

# fewer segregating sites than expected → candidate positive selection

# more segregating sites than expected → candidate balancing selection

arrows(x0 = -5, y0 = 29.75,      # start of arrow
       x1 = -1.5, y1 = 29.75,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey20",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(-3.25, 30.25, cex = 0.95, col = "grey30",
     "positive selection candidate", xpd = TRUE)

arrows(x0 = 5, y0 = 29.75,      # start of arrow
       x1 = 1.5, y1 = 29.75,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey30",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(3.25, 30.25, cex = 0.95, col = "grey30",
     "balancing selection candidate", xpd = TRUE)

dev.off()


