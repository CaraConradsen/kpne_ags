# compare to panX db

par(mfrow =c(1,2), oma = c(0,0,1,0))
panx_ec_ST = fread("C:/Users/carac/Dropbox/Eyre-Walker_work/data/ecoliST",
                   fill=TRUE)[,c(3)] 
colnames(panx_ec_ST) = "ST"
panx_ec_ST[,n:=.N, ST]

unique(panx_ec_ST)[,n] %>% hist(col = "dodgerblue2",
                                main = "ecoli",
                                ylab = "Number of STs",
                                xlab = "Number of genomes",
                                breaks = max(panx_ec_ST$n))

panx_sa_ST = fread("C:/Users/carac/Dropbox/Eyre-Walker_work/data/saureusST",
                   fill=TRUE)[,c(3)] 
colnames(panx_sa_ST) = "ST"
panx_sa_ST[,n:=.N, ST]

unique(panx_sa_ST)[,n] %>% hist(col = "firebrick2",
                                main = "saureus",
                                ylab = "Number of STs",
                                xlab = "Number of genomes",
                                breaks = max(panx_sa_ST$n))
mtext("ST sampling for PanX genomes", side = 3, outer = TRUE)

# set set
set.seed(42)

# Retain only correct oriented chromosomes --------------------------------
# 328 genomes
pan_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))
correct_start <- pan_anno[gene=="dnaA" & strand == "+" & 
                            grepl("_1$", seqnames) & start == 70, 
                          .(geno_id, ST)]

setorderv(correct_start, cols = "geno_id")

geno_list = correct_start[,geno_id]# 328 genomes

chr_contigs <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
                            select = c("seqnames", "asmbly_type")))

# REMOVE PLASMIDS
chr_contigs <- chr_contigs[asmbly_type!="plasmid"][,asmbly_type:=NULL]


# Export focal chromosomal fasta for mash & ANI------------------------------------------

for (i in geno_list) {
  full_fasta <- readDNAStringSet(lng_rd_fasta_files[grepl(i,lng_rd_fasta_files)])
  if(length(full_fasta)>1){
    contig = chr_contigs[grepl(i, seqnames), seqnames]
    
    #subset full fasta
    full_fasta = full_fasta[contig]
  }
  writeXStringSet(full_fasta, paste0("./input_data/kpne_328_chr_fasta/",i,".fasta"))
}


# create fastani ref list
# Windows paths from list.files
fasta_paths <- list.files(
  "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_328_chr_fasta", 
  pattern = "\\.fasta$", 
  full.names = TRUE
)

# Convert to WSL/Linux paths
fasta_paths <- gsub("^C:", "/mnt/c", fasta_paths)
fasta_paths <- gsub("\\\\", "/", fasta_paths)  # backslash -> forward slash

# Write REFERENCE_LIST
writeLines(fasta_paths, "REFERENCE_LIST")


# Split into STs; mash dist fastANI ---------------------------------------
ST_grps <- unique(pan_anno[,.(ST, geno_id)])

ST_grps <- ST_grps[geno_id %chin% geno_list]

ST_grps[, ST_brd:= tstrsplit(ST,"-", fill = TRUE, keep = 1)]

# plot STs
ST_grps[, n_brd:=.N, by =c("ST_brd")]
ST_grps[, n:=.N, by =c("ST_brd","ST")]

st_brd_grps = unique(ST_grps[,ST_brd])

mash_dist <- fread("./input_data/kpne_328_chr_fasta/mash_dist.tsv")
mash_dist[, V1:= gsub(".fasta", "", V1)]
mash_dist[, V2:= gsub(".fasta", "", V2)]
mash_dist <- mash_dist[V1!=V2]


# par(mfrow = c(2,2))
# for (t in c(5e-04, 1e-04, 5e-05, 2e-05)) {
  
mash_threshold = 2e-05

rep_genos <- foreach(i = st_brd_grps,
                     .packages = c("data.table", "stats"),
                     .combine = c) %do%  {
                       st_genos <- ST_grps[ST_brd==i, geno_id]
                       
                       if(length(st_genos)> 1){
                       
                       mash_dat <- dcast(mash_dist[V1 %chin% st_genos & V2 %chin% st_genos],
                                         V1~V2, value.var = "V3", drop = TRUE)
                       mash_mat <- as.matrix(mash_dat[,-1])
                       
                       row.names(mash_mat) = mash_dat$V1
                       
                       # cut clones
                       
                       hc <- hclust(as.dist(mash_mat), method = "average")
                       clusters <- cutree(hc, h = mash_threshold)
                       
                       # Get unique cluster IDs
                       unique_clusters <- unique(clusters)
                       
                       # Sample one genome per cluster
                       representatives <- sapply(unique_clusters, function(clust_id) {
                         genomes_in_cluster <- names(clusters)[clusters == clust_id]
                         sample(genomes_in_cluster, 1)
                       })
                       
                       # if(length(representatives) > 5){
                       #   representatives = sample(representatives, 5, replace = FALSE)
                       # }
                       
                       }else{
                         representatives = st_genos
                       }
                       
                       list(data.frame(ST_brd = i, genomes = representatives))
                     }

rep_genos <- rbindlist(rep_genos)

rep_genos[, n:=.N, ST_brd]

fwrite(rep_genos[,.(genomes)],
       paste0(outdir_dat, "/st_no_clones_mash2e-05.csv"))


# 
# unique(rep_genos[,.(ST_brd,n)])[,n] %>% hist(ylab = "Number of STs",
#                                              xlab = "Number of genomes",
#                                              col = "dodgerblue",
#                                              ylim = c(0,75), yaxt = "n",
#                                              breaks = as.numeric(rep_genos[,.(max(n))]),
#                                              main = paste0("Mash cutoff = ",mash_threshold,
#                                                            "\nn genomes = ", nrow(rep_genos)))
# axis(side = 2, at = seq(0,70, 10),
#      labels = seq(0,70, 10), las = 2)
# }

# # Plot HC
# par(mar = c(1,6,4,2))
# plot(hc, cex = 0.7,
#      xlab = "", ylab = "",
#      yaxt = "n",
#      main = "ST512 Cluster Dendrogram")
# axis(side = 2, at = seq(0,6e-4, 1e-4), las = 2,
#      labels = sprintf("%1.4f", seq(0,6e-4, 1e-4)))
# 
# axis(side = 2, at = 0.00025,
#      line = 3.5,
#      tick = FALSE,
#      labels = "Mash Distance")
# abline(h = mash_threshold, col = "red", lty = 2)
# text(15, 0.00002,"< 2e-05", pos = 3, col = "red")


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

all_hc_plot_dat <- merge(all_hc_plot_dat, ST_grps, 
                         by = "geno_id")

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



