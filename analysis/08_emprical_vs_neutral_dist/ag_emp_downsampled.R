
# Get pangraph information ------------------------------------------------

pangraph_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                                  "number_genomes", "start","end", "average_dose",
                                  "strand","ST","ag_type"))


# remove estimates of core, paralogs and singletons
pangraph_anno <- pangraph_anno[ag_type!="core"][number_genomes > 1][average_dose <= 1]


# Import 2-fold and 4-fold ag sites est -----------------------------------

syn_sites_ag_dt <- fread(paste0(outdir_dat, "/ag_aln_sites.csv"))

# subset by set 1
set_1_names <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

# gene length 
ag_m <- unique(syn_sites_ag_dt[, .(gene_family, gene_length)])
ag_m[, set_1 := fcase(gene_family %in% set_1_names , "Y",
                      default = "N")]

with(ag_m,
     boxplot(gene_length~set_1,
             ylim =c(0, 2500),
             boxwex  = 0.25,
             pch =  16, col = "dodgerblue"))

quantile(ag_m[set_1=="Y", gene_length], probs = c(0.1, 0.25))
# 10% 25% 
# 246 435

# set m target --------------------------------------------------------------

writeLines(as.character(quantile(ag_m[set_1=="Y", gene_length], 
                                 probs = c(0.25))), 
           paste0(outdir_dat, "/m_target.txt"))

# inspect distribution ----------------------------------------------------

png(filename = paste0(outdir_fig, "/set_1_gene_len_dist.png"),
    width = 8.5,height = 9,
    units = "cm",res = 300)

par(mar = c(4.5,4,0.5,0.5))
with(ag_m[set_1=="Y",],
     hist(gene_length,
          yaxt = "n",
          main = "",
          xlab = "Gene length (bp)",
          border = "grey30",
          col = "grey30", breaks = 100))
axis(side = 2, at = seq(0,350, 50),
     labels = seq(0,350, 50), las = 2)
abline(v = 435, lty = 2, col = "red")
text(435, 300, expression("Q"[1]~"="~435*"bp"), 
     pos = 4)

dev.off()

syn_sites_ag_dt <- syn_sites_ag_dt[gene_family %chin% set_1_names]

# Calculate AG segregating sites for target gene size ------------------------

ag_genes_hypergeometric <- function(sites_dt, m_target = 435, min_sites = 2) {
  
  # drop genes shorter than m_target
  sites_dt <- sites_dt[gene_length >= m_target]
  
  # per-gene processing 
  sites_dt[, {
    Ls <- .N
    
    Ls_keep <- rhyper(1, Ls, gene_length[1] - Ls, 375)
    
    if (Ls_keep < min_sites) {
      # return an empty row; filtered out below
      .(Ls = Ls, Ls_keep = Ls_keep, pi_s = NA_real_, 
        n_median = NA_real_, seg_s_sites = NA_integer_)
    } else {
      idx      <- sample.int(Ls, size = Ls_keep)
      n_j      <- n[idx]
      counts   <- cbind(A[idx], C[idx], G[idx], T[idx])
      
      p_mat    <- counts / n_j
      h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))
      seg      <- rowSums(counts > 0) >= 2
      
      .(Ls        = as.integer(Ls),
        Ls_keep   = as.integer(Ls_keep),
        pi_s      = as.numeric(mean(h)),
        n_median  = as.numeric(median(n_j)),
        seg_s_sites = as.integer(sum(seg)))
    }
  }, by = .(gene_family, class, gene_length)]
}


ag_seg_sites_dt <- ag_genes_hypergeometric(syn_sites_ag_dt)


# add back frequency
ag_seg_sites_dt <- merge(ag_seg_sites_dt,
                         unique(pangraph_anno[, .(gene_family, number_genomes)]),
                         all.x = TRUE, by = "gene_family")

# Save output
fwrite(ag_seg_sites_dt, paste0(outdir_dat, "/ag_seg_sites_dt.csv"))
