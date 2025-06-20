# Fix up plots
temp_plamsid_dat <- fread(paste0(outdir_dat, "/kpneu_plasmids.csv"), select = c("contig_id", "origin"))
colnames(temp_plamsid_dat) = c("contig", "contig_origin")


# Import modified gffs to get loci and contig info-----------------------------------------
# get gff info
gff_files_list = list.files("./input_data/PIRATE_1695_out/modified_gffs", 
                            pattern = ".gff", full.names = T)

# start timer
start.time <- Sys.time()

cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Process assembly data
modified_annotation <- foreach(i = 1:length(gff_files_list), 
                               .combine = rbind, 
                               .packages = c("data.table", "rtracklayer")) %dopar% { 
                                 # Try to import GFF file, skip if it fails
                                 gff_dt <- 
                                   as.data.table(rtracklayer::import(gff_files_list[i], 
                                                                     colnames = c("type", "locus_tag", "prev_locus",
                                                                                  "product")))
                                 
                                 
                                 
                                 # get only focal CDS
                                 gff_dt <- gff_dt[grepl("CDS", type)]
                                 
                                 # Return processed data
                                 gff_dt
                               }




# Stop the parallel backend after the loop is done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 10.45449 mins


modified_annotation$type = NULL

colnames(modified_annotation)[c(1)] = c("contig")

fwrite(modified_annotation, paste0(outdir_dat, "/modified_annotation.csv"))
# modified_annotation <- fread(paste0(outdir_dat, "/modified_annotation.csv"))

# Import long pirate data and id loci -------------------------------------


kpne_pangome_dt_lng <- fread(paste0(outdir_dat, "/kpne_pangome_dt_lng.csv"))

colnames(kpne_pangome_dt_lng)[6] = "locus_tag"


# merge everything together -----------------------------------------------

pangenome_data  = merge(modified_annotation, kpne_pangome_dt_lng,
                        all.x = TRUE, by = "locus_tag")

pangenome_data  = merge(pangenome_data, temp_plamsid_dat,
                        all.x = TRUE, by = "contig")


pangenome_data_chr <- pangenome_data[contig_origin!="plasmid" & !is.na(gene_family)]

pangenome_data_chr[, new_n_geno := uniqueN(geno_id), by = gene_family]




# Assign groups -----------------------------------------------------------
pangenome_data_chr[, pan_freq := new_n_geno/1695]

pangenome_data_chr[, pan_grp := fcase(new_n_geno == 1695, "core",
                                      new_n_geno < 1695 & pan_freq >= 0.95, "soft-core",
                                      pan_freq < 0.95 & pan_freq >= 0.15, "shell",
                                      default = "cloud")]

# account for paralogs
pangenome_data_chr[, gene_group := tstrsplit(gene_family, "_", fill = TRUE,
                                             keep = 1L)]

temp_count = unique(pangenome_data_chr[, .(gene_group, gene_family, new_n_geno)])
temp_count[, para_rank := frank(-new_n_geno, ties.method = "min"), 
                   by = gene_group]

pangenome_data_chr = merge(pangenome_data_chr, temp_count[,.(gene_group, gene_family, para_rank)],
                           all.x = TRUE, by = c("gene_group", "gene_family"))

# Inspect paralogs - 
unique(pangenome_data_chr[para_rank > 1, .(product)]) %>% View()

fwrite(pangenome_data_chr, paste0(outdir_dat, "/pangenome_data_chr.csv"))
# pangenome_data_chr <- fread(paste0(outdir_dat, "/pangenome_data_chr.csv"))

# Select fasta to estimate diversity on -----------------------------------


# === Run on all files in directory ===
input_dir <- "./input_data/PIRATE_1695_out/feature_sequences"  # change to your actual folder
fasta_files <- list.files(input_dir, pattern = "nucleotide.fasta", full.names = TRUE)
fasta_files <- as.data.table(fasta_files)

# Subset to AGs
chunk_size <- 100

ag_names <- unique(pangenome_data_chr[pan_grp != "core" & new_n_geno > 1 & para_rank == 1, gene_family])# remove singletons
ag_locus_tag <- unique(pangenome_data_chr[pan_grp != "core" & new_n_geno > 1 & para_rank == 1, .(gene_family, locus_tag)])

chunks <- split(ag_names, ceiling(seq_along(ag_names) / chunk_size))

# Apply each chunk's pattern to subset dt
fasta_files_flt <- lapply(chunks, function(chunk) {
  ag_pattern <- paste(chunk, collapse = "|")
  fasta_files[grepl(ag_pattern, fasta_files, perl = TRUE)]
})

# Combine results into one data.table
fasta_files_flt <- rbindlist(fasta_files_flt, use.names = TRUE, fill = TRUE)

#detect and remove duplicates
fasta_files_flt[, n:=.N, by = fasta_files]

fasta_files_flt <- fasta_files_flt[n <2]

fasta_files_flt$n = NULL

fwrite(fasta_files_flt, paste0(outdir_dat, "/fasta_files_flt.csv"))
# fasta_files_flt <- fread(paste0(outdir_dat, "/fasta_files_flt.csv"))

fasta_files_flt <- fasta_files_flt[,fasta_files]

# Estimate nuc.div ----------------------------------------------------------------
library("silentsitesdiv")


# start timer
start.time <- Sys.time()

cl <- makeCluster(6)
registerDoParallel(cl)

# function
ag_pis_dt <- foreach(i = fasta_files_flt,
                     .packages = c("data.table", "silentsitesdiv", "Biostrings"),
                     .combine = rbind) %dopar% {
                       
                       gene_fam = gsub(".nucleotide.fasta", "",basename(i))
                       
                       focal_ag_aln <- Biostrings::readDNAStringSet(i)
                       
                       focal_ag_aln <- focal_ag_aln[names(focal_ag_aln) %in% ag_locus_tag[gene_family == gene_fam, locus_tag]]
                       
                       piS_dist <- silentsitesdiv::pairwise_piS(as.character(focal_ag_aln))
                       
                       max_piS <- max(piS_dist, na.rm = TRUE)
                       
                       avg_piS <- mean(piS_dist, na.rm = TRUE)
                       
                       print(paste(gene_fam))
                       
                       data.frame(gene_family = gene_fam, max_piS = max_piS, avg_piS = avg_piS)
                     }

ag_pis_dt <- as.data.table(ag_pis_dt)

# Stop the parallel backend after the loop is done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 53.17464 mins

fwrite(ag_pis_dt, paste0(outdir_dat, "/ag_pis_dt.csv"))
# ag_pis_dt <- fread(paste0(outdir_dat, "/ag_pis_dt.csv"))



# 
# 
# sixers <-unique(pangenome_data_chr[new_n_geno == 6, gene_group])
# vals <-  foreach(i = 1: length(sixers), 
#                  .combine = "rbind", 
#                  .packages = c("data.table")) %do% {
#   result <- file.exists(paste0("./input_data/PIRATE_1695_out/feature_sequences/",sixers[i],".nucleotide.fasta"))
#   data.table(data.frame(geno = sixers[i], exist = result[1]))
#                  }
# 
# nsixers <- vals[exist==TRUE, geno]
# 
# for (i in 1:length(nsixers)) {
#   dna_bin <- Biostrings::readDNAStringSet(paste0("./input_data/PIRATE_1695_out/feature_sequences/",nsixers[i],".nucleotide.fasta"))
#   dna_bin <- dna_bin[names(dna_bin) %in% ag_locus_tag[gene_family == nsixers[i], locus_tag]]
#   cat("\n\n", nsixers[i], "\n")
#   print(unique(dna_bin))
# }
# 
# explore synonymous estimates in  g017612 (no missing) and g017688

focal_ag_aln <- Biostrings::readDNAStringSet("./input_data/PIRATE_1695_out/feature_sequences/g017612.nucleotide.fasta")

focal_ag_aln <- focal_ag_aln[names(focal_ag_aln) %in% ag_locus_tag[gene_family == "g017612", locus_tag]]




# Run kaks ----------------------------------------------------------------

# start timer
start.time <- Sys.time()

cl <- makeCluster(12)
registerDoParallel(cl)

ag_age_dt <- foreach(i = fasta_files_flt,
                     .packages = c("data.table", "seqinr", "Biostrings"),
                     .combine = rbind) %dopar% {

                       gene_fam = gsub(".nucleotide.fasta", "",basename(i))

                       focal_ag_aln <- Biostrings::readDNAStringSet(i)

                       focal_ag_aln <- focal_ag_aln[names(focal_ag_aln) %in% ag_locus_tag[gene_family == gene_fam, locus_tag]]


                       # Convert to "alignment" class
                       alignment_obj <- list(
                         nb = length(focal_ag_aln),
                         nam = names(focal_ag_aln),
                         seq = as.character(focal_ag_aln),  # must be character vector
                         com = NULL
                       )
                       class(alignment_obj) <- "alignment"

                       kaks_result <- try(seqinr::kaks(alignment_obj, rmgap = TRUE), silent = TRUE)

                       max_ks <- max(kaks_result$ks, na.rm = TRUE)
                       
                       avg_ks <- mean(kaks_result$ks, na.rm = TRUE)

                       data.frame(gene_family = gene_fam, max_ks = max_ks, avg_ks = avg_ks)
                     }

ag_age_dt <- as.data.table(ag_age_dt)

# Stop the parallel backend after the loop is done
stopCluster(cl)

# get duration
end.time <- Sys.time()
end.time - start.time# Time difference of 22.69124 mins

fwrite(ag_age_dt, paste0(outdir_dat, "/ag_age_dt.csv"))
# ag_age_dt <- fread(paste0(outdir_dat, "/ag_age_dt.csv"))

# Add frequency -----------------------------------------------------------

freq_plot_dat <- merge(ag_age_dt, 
                       unique(pangenome_data_chr[, .(gene_family, pan_freq, pan_grp, threshold)]),
                       all.x = TRUE, by = "gene_family")

freq_plot_dat <- merge(freq_plot_dat, ag_pis_dt,
                       all.x = TRUE, by = "gene_family")

freq_plot_dat[, col := fcase(pan_grp == "soft-core", "#D55E00",
                             pan_grp == "shell", "#56B4E9", 
                             default = "#009E73"
                             )]

with(freq_plot_dat,
     plot(avg_ks, avg_piS, 
          pch = 16, col = col))
abline(a=0, b=1, lty = 2)

with(unique(freq_plot_dat[,.(pan_grp, col)]), 
     legend("top", bty = "n", 
            title = "AG type",
            legend = pan_grp, 
            fill = col)
)

with(freq_plot_dat,
     plot(pan_freq, avg_piS, 
          pch = 16, col = col))
abline(a=0, b=1, lty = 2)

fwrite(freq_plot_dat, paste0(outdir_dat, "/freq_plot_dat.csv"))


# Thresholds --------------------------------------------------------------
library(viridis)

# Example mapping
threshold <- c(98, 95, 90, 80, 70, 60, 50)

# Get colours
threshold_col <- data.frame(threshold_col = viridis(n = length(threshold))[rank(-threshold)],
                            threshold = threshold)

freq_plot_dat <- merge(freq_plot_dat, threshold_col, all.x = TRUE, by="threshold")

# plot --------------------------------------------------------------------
png(filename = paste0(outdir_fig, "/piS_vs_avgKs.png"),
    width = 3840, height = 4473, units = "px", pointsize = 11,
    res = 500, type = "cairo-png")

layout(matrix(c(1,1,2:5), nrow = 3, ncol = 2, byrow = TRUE),
       height = c(0.85,0.85, 1.5))

par(mar = c(4,12,4,20), xpd = TRUE)
with(freq_plot_dat, 
     plot(avg_ks, avg_piS, pch=16,bty="L",
          xlab = expression("Average " * K[S]),
          ylab = expression(pi[S]),
          col = threshold_col))
abline(a=0, b=1, lty = 2, xpd = FALSE)
text(1.95, 0.95, "1:1 line")

mtext(expression("A. My estimate of " * pi[S] * " vs. average "* K[S] * " (kaks)"), 
      line = 1, adj = -0.25, cex = 0.75)

with(threshold_col, 
     legend("topleft", bty = "n", inset = c(1.05, 0),
            title = "Percentage identity threshold\nof gene family",
            legend = threshold, cex=1.05,
            fill = threshold_col)
)

par(mar = c(4,4,4,4))
with(freq_plot_dat,
     plot(pan_freq, avg_piS, pch=16,bty="L",
          ylab = expression(pi[S]),
          xlab = 'Frequency in the pangenome',
          main = "",
          col = threshold_col))#rgb(0.15,0.1,0.7, alpha = 0.25)

mtext(expression("B. My estimate of " * pi[S] * " vs. pangenome frequency"), 
      line = 1, adj = -0.25, cex = 0.75)

with(freq_plot_dat,
     plot(pan_freq, avg_ks, pch=16,bty="L",
          ylab = expression("Average " * K[S]),
          xlab = 'Frequency in the pangenome',
          main = "",
          col = threshold_col))#rgb(0.15,0.1,0.7, alpha = 0.25)
mtext(expression("C. Kaks' average " * K[S] * " vs. pangenome frequency"), 
      line = 1, adj = -0.25, cex = 0.75)

par(mar = c(3,4,5,4))
# Density estimation
dens <- kde2d(y = freq_plot_dat$avg_piS, x = freq_plot_dat$pan_freq, n = 50)
max_den = round(max(dens$z))
trunc_val = 40
dens$z[dens$z > trunc_val] <- trunc_val # Cap values at 40
# plot
persp(dens, col = "lightblue",theta = 40, ticktype = "detailed", 
      ltheta = 120,shade = 0.5, xlab = "Pangenome frequency",
      zlim = c(0, trunc_val), 
      main = paste0("\n\n(data truncated at density = 40; max density = ",max_den,")"),
      zlab = "Density", ylab = "\u03C0S")

mtext(expression("D. Density distribution of " * pi[S]), 
      line = 2, adj = -0.15, cex = 0.75)


# Density estimation
dens <- kde2d(y = freq_plot_dat$avg_ks, x = freq_plot_dat$pan_freq, n = 50)
max_den = round(max(dens$z))
trunc_val = 40
dens$z[dens$z > trunc_val] <- trunc_val # Cap values at 40
# plot
persp(dens, col = "lightblue",theta = 40, ticktype = "detailed", 
      ltheta = 120,shade = 0.5, xlab = "Pangenome frequency",
      zlim = c(0, trunc_val), 
      main = paste0("\n\n(data truncated at density = 40; max density = ",max_den,")"),
      zlab = "Density", ylab = "Average KS")
mtext(expression("E. Density distribution of average " * K[S]), 
      line = 2, adj = -0.15, cex = 0.75)

dev.off()



# Looking at avgKs < 0.1 and threshold >95 --------------------------------
layout(matrix(c(1:2), nrow = 1, ncol = 2, byrow = TRUE))

low_vals = freq_plot_dat[threshold >= 95 & avg_ks <= 0.1]# & pan_grp != "soft-core"]

# Density estimation
dens <- kde2d(y = low_vals$avg_ks, x = low_vals$pan_freq, n = 50)
max_den = round(max(dens$z))
trunc_val = 100
dens$z[dens$z > trunc_val] <- trunc_val # Cap values at 40
# plot
persp(dens, col = "lightblue",theta = 40, ticktype = "detailed", 
      ltheta = 120,shade = 0.5, xlab = "Pangenome frequency",
      zlim = c(0, trunc_val),
      main = "",
      zlab = "Density", ylab = "\u03C0S")
mtext(paste0("A. all accessory genes \n(identity threshold \u2265 95%, \u03C0S \u2264 0.1, n = ", nrow(low_vals),")"),
      line = 2, adj = -0.15)

# no soft core
low_vals = freq_plot_dat[threshold >= 95 & avg_ks <= 0.1 & pan_grp != "soft-core"]

# Density estimation
dens <- kde2d(y = low_vals$avg_ks, x = low_vals$pan_freq, n = 50)
max_den = round(max(dens$z))
trunc_val = 100
dens$z[dens$z > trunc_val] <- trunc_val # Cap values at 40
# plot
persp(dens, col = "lightblue",theta = 40, ticktype = "detailed", 
      ltheta = 120,shade = 0.5, xlab = "Pangenome frequency",
      zlim = c(0, trunc_val),
      main = "",
      zlab = "Density", ylab = "\u03C0S")
mtext(paste0("B. no soft-core genes \n(identity threshold \u2265 95%, \u03C0S \u2264 0.1, n = ", nrow(low_vals),")"),
      line = 2, adj = -0.15)

# only soft core
low_vals = freq_plot_dat[threshold >= 95 & avg_ks <= 0.1 & pan_grp == "soft-core"]


