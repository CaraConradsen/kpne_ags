# set set
set.seed(42)

mash_dist <- fread("./input_data/kpne_412_chr_fasta/mash_dist.tsv")
mash_dist[, V1:= gsub(".fasta", "", V1)]
mash_dist[, V2:= gsub(".fasta", "", V2)]
mash_dist <- mash_dist[V1!=V2]


# Use mash against all genomes --------------------------------------------
min(mash_dist$V3[mash_dist$V3!=0])# equals 2.38274e-05

mash_threshold = 2e-05

mash_dat <- dcast(mash_dist,
                  V1~V2, 
                  value.var = "V3",
                  drop = TRUE)

mash_mat <- as.matrix(mash_dat[,-1])

row.names(mash_mat) = mash_dat$V1


hc_all <- hclust(as.dist(mash_mat), method = "average")

# visulaise dendogram
par(mar = c(4,5,4,2))
plot(hc_all, cex = 0.5, labels = FALSE,
     xlab = "", ylab = "",
     yaxt = "n",xaxt = "n",
     main = "Cluster Dendrogram")

axis(side = 2, at = seq(0,0.012, length.out = 7), las = 2,
     labels = sprintf("%1.3f", seq(0,0.012, length.out = 7)))

axis(side = 2, at = 0.006, line = 3, 
     labels = "Cluster merge height (average Mash distance)",
     tick = FALSE)

axis(side = 1, at = 164, cex = 10,
     line = -1.75,
     labels = "Genomes",
     tick = FALSE)

# add ST points

hang_height = max(hc_all$height) * 0.1

n <- length(hc_all$labels)

# For each leaf, find the height of its first merge
terminal_height <- sapply(seq_len(n), function(i) {
  min(hc_all$height[hc_all$merge[,1] == -i | hc_all$merge[,2] == -i])
})

terminal_height = terminal_height - hang_height

# X positions from plotting order
x_pos <- match(seq_len(n), hc_all$order)

# create data.table
all_hc_plot_dat <- data.frame(x_pos = x_pos, y_pos = terminal_height,
                              geno_id = hc_all$labels)
setDT(all_hc_plot_dat)

ST_grps = fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv",
                select = c("id", "ST"))

colnames(ST_grps)[1] = "geno_id"

all_hc_plot_dat <- merge(all_hc_plot_dat, ST_grps, 
                         by = "geno_id")

all_hc_plot_dat[, ST_brd := tstrsplit(ST, "-", fill = TRUE, keep = 1)]

all_hc_plot_dat[, n_brd := .N , ST_brd]

plot_cols = sort(unique(all_hc_plot_dat[n_brd>1, ST]))

plot_cols = cbind(plot_cols, 
                  rainbow(length(plot_cols)))
colnames(plot_cols) = c("ST", "col")

all_hc_plot_dat <- merge(all_hc_plot_dat, plot_cols, 
                         all.x = TRUE, by = "ST")

all_hc_plot_dat[is.na(col), col := "black"]

# Plot dendrogram and overlay points
with(all_hc_plot_dat, 
     points(x_pos, y_pos, 
            pch = 16, cex = 0.8,
            col = col)
)

abline(h = mash_threshold, col = "red", lty = 2)
text(332, 0.00002,"< 2e-05", pos = 3, col = "red", xpd=TRUE)

# add legend
plot_info = unique(all_hc_plot_dat[col!="black", .(ST, col)])
with(plot_info[1:11,],
     legend(60, -0.002, legend = ST, col = col,
            x.intersp = 0.4,
            bty="n", horiz = TRUE,cex = 0.6,
            pch = 16, xpd = TRUE))
with(plot_info[12:22,],
     legend(60, -0.00225, legend = ST, col = col,
            x.intersp = 0.4,
            bty="n", horiz = TRUE,cex = 0.6,
            pch = 16, xpd = TRUE))
with(plot_info[23:31,],
     legend(60, -0.0025, legend = ST, col = col,
            x.intersp = 0.4,
            bty="n", horiz = TRUE,cex = 0.6,
            pch = 16, xpd = TRUE))
text(50, -0.0025,"ST (n>1)", xpd = TRUE)

# cut clones
clusters <- cutree(hc_all, h = mash_threshold)

# Get unique cluster IDs
unique_clusters <- unique(clusters)

# Sample one genome per cluster
representatives <- sapply(unique_clusters, function(clust_id) {
  genomes_in_cluster <- names(clusters)[clusters == clust_id]
  sample(genomes_in_cluster, 1)
})

write(representatives,
       paste0(outdir_dat, "/all_no_clones_mash2e-05.csv"))


# Output representative genomes -------------------------------------------

# Define source and destination directories
src_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_412_chr_fasta"
dest_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_260_chr_fasta"

for (i in representatives) {
  
  # Example file
  file = paste0(i,".fasta")
  
  # Construct full paths
  src_file  <- file.path(src_dir, file)
  dest_file <- file.path(dest_dir, file)

  # Move the file
  file.rename(src_file, dest_file)
  
}

convert_files = list.files("./input_data/kpne_260_chr_fasta/",
                           pattern = ".fasta", full.names = TRUE)

for (i in convert_files) {
  
  name = gsub(".fasta", ".fna", basename(i))
  
  cat("done ..", name)
  
  fna <- Biostrings::readDNAStringSet(i)
  
  Biostrings::writeXStringSet(fna, 
                              paste0("./input_data/kpne_260_chr_fna/", name))
}

