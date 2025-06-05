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

pangenome_data_chr[, pan_grp := fcase(pan_freq > 0.99, "core",
                                      pan_freq <= 0.99 & pan_freq >= 0.95, "soft-core",
                                      pan_freq < 0.95 & pan_freq >= 0.15, "shell",
                                      default = "cloud")]

# plot distributions ------------------------------------------------------

pie_dat = unique(pangenome_data_chr[,.(gene_family, pan_grp)])

pie_dat = pie_dat[,.(n = .N), by=pan_grp]

pie_dat[, lab := paste0(pan_grp, " (", n,")")]

# thresholds
pirate_thresholds = unique(pangenome_data_chr[,.(gene_family, threshold)])

pirate_thresholds <- pirate_thresholds[,.(n=.N), by = threshold]

setorderv(pirate_thresholds, "threshold", 1L)

# accessory genes
ag_freq = unique(pangenome_data_chr[pan_grp!="core",.(gene_family, pan_freq)])


# png("C:/Users/carac/Dropbox/Vos_Lab/working_pictures/pangenome_stats.png",
#     width = 12.25, height = 5.25, units = "in", res = 300, 
#     pointsize = 20,
#     bg = "transparent", type = "cairo")

# plots
mat <- matrix(c(1,2,3,4,5,3), byrow =TRUE, ncol = 3)

layout(mat, widths = c(0.75,0.75, 1),
       heights = c(1,0.55))
# layout.show(n=5)


par(mar = c(4,5,4,1))
bar_thres<- with(pirate_thresholds,
                 barplot(n~threshold,
                         border = NA,
                         yaxt = "n",
                         col = "grey40",
                         # main = "",
                         main = "Percentage identity of gene \nfamily at lowest threshold  \n",
                         xpd=TRUE,
                         ylab = "",
                         xlab = ""))
abline(h=0)
axis(side=2, at = seq(0,10000, 2000),
     labels = seq(0,10000, 2000),
     las = 2)
axis(side=2, at = 10000/2, line = 2.5,
     labels = "Count", tick = F)
axis(side=1, at = bar_thres[4], line = 1.5,
     labels = "Threshold", tick = F)


# Green-themed colours
green_palette <- c("grey70", "grey40",  "yellowgreen", "#1b5837")  # dark to light greens

par(mar = c(0.5,1,1,6))
order_idx <- c(3, 4, 1, 2)
# Create pie chart
with(pie_dat,
     pie(n[order_idx],
         labels = c(lab[order_idx][1:3], ""),
         col = green_palette,
         main = "\n\nPangenome",
         border = "#fafaf8",
         clockwise = TRUE
     )
)

# adjust label
angle_radians <- 2 * pi * 14258/2 / sum(pie_dat$n)

# Radius at which to place labels (1 is on the circumference)
r <- -1.2
x_pos <- r * cos(angle_radians)
y_pos <- r * sin(angle_radians)

# Add text at custom positions
text(x_pos, y_pos, xpd=TRUE,
     labels = pie_dat[pan_grp=="cloud", lab])


# ag frequency
par(mar = c(4,5,4,0.5))

ag_freq_hist <- hist(ag_freq[,pan_freq],
                     breaks = 20, plot = FALSE
)

# Assign colours based on left edge of each bin
bar_col <- cut(ag_freq_hist$breaks[-length(ag_freq_hist$breaks)],
               breaks = c(-Inf, 0.15, 0.95, Inf),
               labels = c("#1b5837", "yellowgreen",  "grey40"),
               right = FALSE)  # Include left edge, exclude right

# Plot histogram with coloured bars
plot(ag_freq_hist,
     main = "Frequency of acessory genes",
     ylab="", yaxt="n",
     xlab = "",
     col = as.character(bar_col)
)

legend("topright", bty="n", fill = c("#1b5837", "yellowgreen",  "grey40"),
       legend = c("cloud, <0.15", "shell, <0.95 ","soft-core, <0.99"))

pretty_axis <- pretty(range(ag_freq_hist$counts))

axis(side=2, at = pretty_axis,
     labels = pretty_axis,
     las = 2)
axis(side=2, at = 7000, line = 2.5,
     labels = "Number of genes", tick = F)
axis(side = 1, at = ag_freq_hist$breaks[11],
     labels = "Genome frequency", line =1.5, 
     tick = FALSE)

plot.new();plot.new()

# dev.off()




# Syntenic blocks for AGs --------------------------------------------------------

syntenic_blocks <- fread("./input_data/PIRATE_1695_out/pangenome.syntenic_blocks.tsv")
colnames(syntenic_blocks) = c("block_id", "n_syn_geno", "n_co_genes", 
                              "genes", "genes_up", "genes_down")

# Subset to AGs
chunk_size <- 1000

ag_names <- unique(pangenome_data_chr[pan_grp!="core", gene_family])

chunks <- split(ag_names, ceiling(seq_along(ag_names) / chunk_size))

# Apply each chunk's pattern to subset dt
ag_syntenic_blocks_ls <- lapply(chunks, function(chunk) {
  ag_pattern <- paste(chunk, collapse = "|")
  syntenic_blocks[grepl(ag_pattern, genes, perl = TRUE)]
})

# Combine results into one data.table
ag_syntenic_blocks <- rbindlist(ag_syntenic_blocks_ls, use.names = TRUE, fill = TRUE)

#detect and remove duplicates
ag_syntenic_blocks[, n:=.N, by = block_id]

ag_syntenic_blocks <- ag_syntenic_blocks[n <2]

ag_syntenic_blocks$n = NULL


# Multiple intros of AGs --------------------------------------------------
ag_gene_fam_dt <- unique(pangenome_data_chr[pan_grp!="core",
                                            .(gene_family, number_genomes,
                                              pan_freq)])

# count the ags
ag_syn_block_count <- foreach(i = ag_gene_fam_dt$gene_family, 
                              .combine = "rbind",
                              .packages = "data.table") %do% {
                                temp_syn <- ag_syntenic_blocks[grepl(i, genes)]
                                
                                count_temp_syn = nrow(temp_syn)
                                
                                temp_syn <- temp_syn[, lapply(.SD, function(x) paste(x, collapse = ",")), 
                                                     .SDcols = c("block_id", "n_syn_geno", "n_co_genes")]
                                
                                temp_syn$gene_family = i
                                
                                temp_syn$n_blocks = count_temp_syn
                                
                                temp_syn
                                
                              }



# Separate paralogs, single insertion and unlocalised ---------------------
ag_syn_block_count[, gene_fam := tstrsplit(gene_family, "_", fixed = TRUE, keep = 1L)]


ag_syn_block_count[, group := fcase(grepl("_", gene_family), "paralog",
                                    !grepl("_", gene_family) & n_blocks == 1, "single_syn_block",
                                    !grepl("_", gene_family) & n_blocks == 0, "unlocalised",
                                    default = NA)]


synten_hist <- unique(ag_syn_block_count[, .(gene_fam, group)])[, .(n=.N), by = group]

synten_hist[, perc := n / sum(n)]    
# Here we have the gene families groups as:
# "single-syntenic block", “Paralogous / Multi-copy (inc. possible HGT)”,
# and "unlocalised"
with(synten_hist,
     barplot(n~group
     ))



# NT plot -----------------------------------------------------------------

# === Run on all files in directory ===
input_dir <- "./input_data/PIRATE_1695_out/feature_sequences"  # change to your actual folder
fasta_files <- list.files(input_dir, pattern = "nucleotide.fasta", full.names = TRUE)
fasta_files <- as.data.table(fasta_files)

# Subset to AGs
chunk_size <- 100

ag_names <- unique(ag_syn_block_count[n_syn_geno !="1", gene_fam])# remove singletons

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

fasta_files_flt <- fasta_files_flt[,fasta_files]

ag_age_dt <- foreach(i = fasta_files_flt,
                     .packages = c("data.table", "seqinr", "Biostrings"),
                     .combine = rbind) %do% {
                       
                       

                       if(loc_freq >1){
                         
                         focal_ag_aln <- Biostrings::readDNAStringSet(i)
                         focal_ag_aln <- focal_ag_aln[!grepl("-", names(focal_ag_aln))]
                         
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

