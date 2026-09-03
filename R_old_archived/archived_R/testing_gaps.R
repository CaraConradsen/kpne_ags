diagnose_outlier <- function(fasta_path,
                             gap_threshold = 0.05,
                             plot_tree = TRUE) {
  
  # Defensive: check file exists and has content
  if (!file.exists(fasta_path)) {
    warning(sprintf("File not found: %s", fasta_path))
    return(NULL)
  }
  
  tryCatch({
    aln <- read.dna(fasta_path, format = "fasta", as.matrix = TRUE)
    n_seq <- nrow(aln)
    aln_len <- ncol(aln)
    
    aln_chr <- as.character(aln)
    gap_prop  <- rowMeans(aln_chr == "-")
    max_gap <- max(gap_prop)
    is_gappy  <- gap_prop > gap_threshold
    col_gap_freq <- colMeans(aln_chr == "-")
    
    d <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE)
    
    if (any(is.na(d))) {
      warning("NA distances present — some pairs share no ungapped sites.")
    }
    
    if (n_seq > 2) {
      tr <- nj(d)
      tr <- ladderize(tr)
      hc       <- hclust(d, method = "average")
      clusters <- cutree(hc, k = 2)
      clusters <- clusters[tr$tip.label]
      is_gappy_aligned <- is_gappy[tr$tip.label]
      tab <- table(cluster = clusters, gappy = is_gappy_aligned)
      ft  <- tryCatch(fisher.test(tab), error = function(e) NULL)
      d_vec <- as.vector(d)
      d_vec <- d_vec[!is.na(d_vec)]
      bimod_coef <- {
        m  <- mean(d_vec); s <- sd(d_vec)
        sk <- mean((d_vec - m)^3) / s^3
        ku <- mean((d_vec - m)^4) / s^4 - 3
        (sk^2 + 1) / (ku + 3 * (length(d_vec) - 1)^2 /
                        ((length(d_vec) - 2) * (length(d_vec) - 3)))
      }
      
      result <- list(
        file              = gsub(".fasta", "", fasta_path),
        n_seq             = n_seq,
        alignment_length  = aln_len,
        n_gappy           = sum(is_gappy),
        prop_gappy        = mean(is_gappy),
        max_gap           = max_gap,
        tree              = tr,
        clusters          = clusters,
        cluster_gap_table = tab,
        fisher_p          = if (!is.null(ft)) ft$p.value else NA_real_,
        mean_dist         = mean(d_vec),
        max_dist          = max(d_vec),
        bimodality_coef   = bimod_coef,
        col_gap_freq      = col_gap_freq
      )
    } else {
      result <- list(
        file              = gsub(".fasta", "", fasta_path),
        n_seq             = n_seq,
        alignment_length  = aln_len,
        n_gappy           = sum(is_gappy),
        prop_gappy        = mean(is_gappy),
        max_gap           = max_gap,
        tree              = NULL,
        clusters          = NULL,
        cluster_gap_table = NULL,
        fisher_p          = NA_real_,
        mean_dist         = NA_real_,
        max_dist          = NA_real_,
        bimodality_coef   = NA_real_,
        col_gap_freq      = col_gap_freq
      )
    }
    
    if (plot_tree && n_seq > 2) {
      op <- par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))
      on.exit(par(op), add = TRUE)
      tip_col <- ifelse(is_gappy_aligned, "tomato", "steelblue")
      plot(tr, tip.color = tip_col, cex = 0.4,
           main = sprintf("NJ tree (n=%d)\nred = gappy seqs", n_seq))
      add.scale.bar()
      d_vec <- as.vector(d)
      d_vec <- d_vec[!is.na(d_vec)]
      hist(d_vec, breaks = 40, col = "grey80", border = "white",
           xlab = "pairwise raw distance",
           main = sprintf("Pairwise distances\nbimod_coef = %.3f", result$bimodality_coef))
      abline(v = mean(d_vec), col = "red", lty = 2)
      plot(col_gap_freq, type = "h", col = "darkgreen", ylim = c(0,1),
           xlab = "alignment column", ylab = "gap frequency",
           main = "Where the gaps fall")
    }
    
    result
    
  }, error = function(e) {
    warning(sprintf("Error processing %s: %s", fasta_path, e$message))
    NULL
  })
}

# ---- Apply to a directory of outlier alignments ----

outlier_dir   <- paste0(outdir_dat, "/ag_fus_alns/")
outlier_files <- list.files(outlier_dir, pattern = ".fasta",
                            full.names = TRUE)

results_list <- lapply(outlier_files, diagnose_outlier, plot_tree = FALSE)

# Keep track of which genes failed
failed_genes <- names(results_list)[sapply(results_list, is.null)]
results_list <- results_list[!sapply(results_list, is.null)]

# Convert to data.table
results_dt <- rbindlist(lapply(results_list, function(x) {
  data.table(
    gene_family = basename(x$file),
    number_genomes = x$n_seq,
    aln_len = x$alignment_length,
    prop_gappy = x$prop_gappy,
    max_gap = x$max_gap,
    mean_dist = x$mean_dist,
    max_dist = x$max_dist,
    bimodality_coef = x$bimodality_coef,
    fisher_p = x$fisher_p,
    cluster1_size = if (!is.null(x$clusters)) sum(x$clusters == 1) else NA_integer_,
    cluster2_size = if (!is.null(x$clusters)) sum(x$clusters == 2) else NA_integer_
  )
}), fill = TRUE)

# Now add back empty genes with NAs
all_genes <- data.table(gene_family = unique(c(results_dt$gene_family, failed_genes)))
summary_dt <- all_genes[results_dt, on = "gene_family"]

# add gene set info

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))[,gene_family]
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))[,gene_family]

summary_dt[, set := fcase(gene_family %in% set_1_ags, 1,
                                 default = 0)]
summary_dt[, set := fcase(gene_family %in% set_2_ags, 2,
                                 default = set)]


pirate <- fread("./input_data/PIRATE_260_hybrid_chr_out/PIRATE.gene_families.ordered.tsv")[,1:21]

summary_dt <- merge(summary_dt, pirate[,.(gene_family,threshold)],
      all.x = TRUE, by = "gene_family")

outliers <- fread("./output/data/set_outliers.csv",
                  select = c("gene_family", "theta_outlier", "w_theta_sel", "V2"))

#remove dups to merge
dupes_set1 <- outliers[, .N, by = gene_family][N > 1, gene_family]
outliers <- outliers[!(V2 == "set_1" & gene_family %in% dupes_set1)]
outliers[,V2:=NULL]

summary_dt <- merge(summary_dt, outliers,
                    all.x = TRUE, by = "gene_family")



# thresholds impacting the clustering
pielou_j <- function(x) {
  p <- x / sum(x)
  H <- -sum(p[p > 0] * log(p[p > 0]))
  H_max <- log(length(x))
  return(H / H_max)
}

summary_dt[, pielou_j := pielou_j(c(cluster1_size, cluster2_size)), by = 1:nrow(summary_dt)]

summary_dt[, cluster_balance := fcase(pielou_j <0.4, "uneven",
                                      default = "even")]

summary_dt[, is_gappy := fcase(max_gap> 0.1, "Y",
                               default = "N")]

barplot_dat <- summary_dt[cluster_balance!="NA" & set != 0,
                          .(n = .N), by = c("threshold", "cluster_balance", "is_gappy")]
setorderv(barplot_dat, cols = c("threshold", "cluster_balance", "is_gappy"))

# Pivot to wide format: each row is threshold + cluster_balance, columns are Y and N counts
barplot_dat_wide <- dcast(barplot_dat, threshold + cluster_balance ~ is_gappy, value.var = "n")

# Create matrix for stacked barplot
mat <- as.matrix(barplot_dat_wide[, .(N, Y)])
mat[is.na(mat)] = 0
rownames(mat) <- paste(barplot_dat_wide$threshold, barplot_dat_wide$cluster_balance, sep = "\n")

# Stacked barplot
barplot(t(mat),
        col = c("lightblue","coral"),
        main = "Number of Outlier AGs by threshold and cluster evenness\nwhen Pielou's J < 0.4 (uneven)",
        ylab = "Number of AGs",
        xlab = "PIRATE cluster threshold",
        ylim = c(0, max(rowSums(mat)) * 1.1),
        legend.text = c("not gappy", "max gap> 0.1"),
        args.legend = list(x = "topleft"))



fwrite(summary_dt, paste0(outdir_dat,"/testing_gaps.csv"))
