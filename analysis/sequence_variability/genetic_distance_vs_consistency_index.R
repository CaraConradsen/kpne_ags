
pirate_alignments <- list.files("./input_data/PIRATE_485_lng_rds_out/feature_sequences/",
                                full.names = TRUE, pattern = ".nucleotide.fasta")

for (i in 1:length(pirate_alignments)) {
  if (nrow(read.dna(pirate_alignments[i], format = "fasta")) < 70 & 
      nrow(read.dna(pirate_alignments[i], format = "fasta")) > 35){
    print(i)
    break
  }
}

val = 441 # g000409_2

focal_align <- read.dna(pirate_alignments[val], format = "fasta")

# JC69 model, with variance estimates
dist_obj <- dist.dna(focal_align, model = "JC69", 
                     variance = TRUE, as.matrix = TRUE)

# Extract the matrices
var_vec <- attr(dist_obj, "variance")
n <- nrow(dist_obj)
labels <- rownames(dist_obj)

# Rebuild full symmetric variance matrix
var_mat <- matrix(0, n, n)
var_mat[lower.tri(var_mat)] <- var_vec
var_mat <- var_mat + t(var_mat)
rownames(var_mat) <- colnames(var_mat) <- labels

mean_dist <- rowMeans(dist_obj, na.rm = TRUE)
var_dist  <- apply(dist_obj, 1, var, na.rm = TRUE)
n_seq <- nrow(dist_obj)
var_of_mean <- var_dist / (n_seq - 1)
se_dist <- sqrt(var_of_mean)


overall_mean <- mean(mean_dist)
z_rel <- (mean_dist - overall_mean) / se_dist


div_df <- data.frame(
  Sequence = rownames(dist_obj),
  Mean_Dist = mean_dist,
  SE_Empirical = se_dist,
  Z_Empirical = z_rel
)

boxplot(div_df$Z_Empirical, pch = 16)

div_df <- div_df[order(div_df$Z_Empirical, decreasing = TRUE), ]

temp <- readDNAStringSet(pirate_alignments[val])

msaR(temp[names(temp) %in% div_df$Sequence[1:15]])
