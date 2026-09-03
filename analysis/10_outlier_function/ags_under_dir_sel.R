


anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                          select = c("fus_locus_tag","geno_id", "gene_family", "ST", "threshold", 
                                     "number_genomes", "consensus_gene_name","consensus_product",
                                     "alleles_at_maximum_threshold")))
anno <- anno[gene_family %chin% c("g000371", "g000222", "g001303", "g001257", 
                                  "g000145", "g001254", "g003299", "g004223")]

unique(anno[,.(gene_family, number_genomes, consensus_gene_name, consensus_product, alleles_at_maximum_threshold)])
#     gene_family number_genomes consensus_gene_name                                                              consensus_product alleles_at_maximum_threshold
# 1:     g000371            210                kexD                         multidrug efflux RND transporter permease subunit KexD                            8
# 2:     g001257            254                ramR                                TetR/AcrR family transcriptional regulator RamR                            3
# 3:     g001254            249       AbiEii/AbiGii                     Nucleotidyl transferase AbiEii/AbiGii toxin family protein                           11
# 4:     g000145            251                ompA                                            OmpA-like domain-containing protein                           74
# 5:     g003299            249               AbiEi Transcriptional regulator AbiEi antitoxin N-terminal domain-containing protein                           11
# 6:     g001303            215                acrA                                 Co/Zn/Cd efflux system membrane fusion protein                           25
# 7:     g000222            215                acrB                                    Cobalt-zinc-cadmium resistance protein CzcA                           61
# 8:     g004223             67                tssF                                type VI secretion system baseplate subunit TssF                            7




# Using aa ----------------------------------------------------------------
# accounting for fusions
faa_files <-  list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_hps0.6_out/feature_sequences/",
                         full.names = TRUE, pattern = ".aa.fasta")

faa_files <- faa_files[grepl(paste(unique(anno$gene_family), collapse = "|"), faa_files)]

# lapply(seq_len(nrow(faa_files)), function(i) {
#   
#   prots <- readAAStringSet(faa_files[i])
#   
#   ag = gsub(".aa.fasta","", basename(faa_files[i]))
#   
#   names <- sub("^\\((.*?)[:;].*$", "\\1", anno[gene_family == ag, fus_locus_tag])
#   
#   prots <- prots[names(prots) %in% names]
# }
# 

anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                     select = c("fus_locus_tag","geno_id", "gene_family", "ST", "threshold", 
                                "number_genomes", "consensus_gene_name","consensus_product",
                                "seqnames", "fstart",  "fend", "fstrand",
                                "alleles_at_maximum_threshold")))
anno <- anno[gene_family %chin% c("g000371", "g000222", "g001303", "g001257", 
                                  "g000145", "g001254", "g003299", "g004223")]

setnames(anno, c( "fstart",  "fend", "fstrand"),
         c( "start",  "end", "strand"))

fna_files <-  list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_260_chr_fasta/",
                         full.names = TRUE, pattern = ".fasta")


focal_ags <- unique(anno$gene_family)

# Set up parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Paralleled loop to read sequences and compute consensus
# Per isolate, per msu, concatenate core gene sequences, compute GC
ag_st_allele_counts <- foreach(g = focal_ags, 
                               .combine = "c", 
                               .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                                 
                                 g_dat <- anno[gene_family == g, .(geno_id, seqnames, start, end, strand, ST)]
                                 
                                 seqs <- lapply(seq_len(nrow(g_dat)), function(i) {
                                   
                                   gr <- g_dat[i]
                                   
                                   genome <- readDNAStringSet(fna_files[grepl(gr$geno_id, fna_files)])
                                   
                                   getSeq(genome, dt2gr(gr[, .(seqnames, start, end, strand)]))
                                 })
                                 
                                 seqs <- do.call(c, seqs)
                                 
                                 names(seqs) <- g_dat$geno_id
                                 
                                 seqs <- suppressWarnings(Biostrings::translate(seqs, if.fuzzy.codon = "solve"))
                                 
                                 allele <- match(as.character(seqs), unique(as.character(seqs)))
                                 
                                 g_dat$allele_name = paste0("a_", allele)
                                 g_dat[, n := .N, by = allele]
                                 
                                 list(dcast(g_dat, allele_name ~ ST, fun.aggregate = length, value.var = "geno_id"))
                                 
                               }

# Stop cluster after execution
stopCluster(cl)

names(ag_st_allele_counts) <- focal_ags



# Using nt seqs -----------------------------------------------------------

# # alleles
# alleles_info <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_hps0.6_out/PIRATE.unique_alleles.tsv",
#                       fill = TRUE, select = c("allele_name", "gene_family", "number_genomes", "threshold",
#                                               unique(anno$geno_id)))
# alleles_info <- alleles_info[gene_family %chin% c("g000371", "g000222", "g001303", "g001257", 
#                                   "g000145", "g001254", "g003299", "g004223")][threshold==100]
# 
# alleles_info[,threshold:=NULL]
# 
# alleles_info <- melt(alleles_info, 
#                    id.vars = c("allele_name","gene_family","number_genomes"),
#                    variable.name = "geno_id", 
#                    value.name = "fus_locus_tag"
# )
# 
# alleles_info <- alleles_info[fus_locus_tag!=""][, c("fus_locus_tag", "number_genomes"):= NULL]
# 
# setorderv(alleles_info, c("gene_family","geno_id", "allele_name"), c(1,1,1))
# 
# alleles_info <- alleles_info[, .(allele_name = paste(sort(unique(allele_name)), collapse = ";")),
#                              by = .(gene_family, geno_id)]
# 
# 
# alleles_vs_st <- merge(anno, alleles_info,
#                        all.x = TRUE, by = c("gene_family", "geno_id"))
# 
# allele_st_tabs <- split(alleles_vs_st, by = "gene_family", keep.by = FALSE)
# 
# allele_st_tabs <- lapply(allele_st_tabs, function(d)
#   dcast(d, allele_name ~ ST, value.var = "geno_id",
#         fun.aggregate = length, fill = 0L))
# 
# drop_zero_cols <- function(d) {
#   num <- names(d)[sapply(d, is.numeric)]
#   keep <- num[colSums(d[, ..num]) > 0]
#   d[, c("allele_name", keep), with = FALSE]
# }
# 
# allele_st_tabs <- lapply(allele_st_tabs, drop_zero_cols)

# 
# # Allele diversity plot ---------------------------------------------------
# 
collapse_st <- function(tab, min_n = 10, merge_lv = TRUE) {
  long <- melt(tab, id.vars = "allele_name",
               variable.name = "ST", value.name = "n")
  long[, ST := as.character(ST)]
  if (merge_lv) long[, ST := sub("-[0-9]+LV$", "", ST)]
  long <- long[, .(n = sum(n)), by = .(allele_name, ST)]

  tot  <- long[, .(N = sum(n)), by = ST][order(-N)]
  keep <- tot[N >= min_n, ST]
  long[, ST_grp := fifelse(ST %in% keep, ST, "other")]

  out <- long[, .(n = sum(n)), by = .(allele_name, ST_grp)]
  lev <- c(tot[ST %in% keep][order(-N), ST], "other")
  out[, ST_grp := factor(ST_grp, levels = lev)]
  attr(out, "n_other_st") <- tot[!ST %in% keep, .N]
  out[]
}
# 
# long_g000371 <- collapse_st(allele_st_tabs[["g000371"]], min_n = 2)
# long_g001257 <- collapse_st(allele_st_tabs[["g001257"]], min_n = 2)
# long_g001303 <- collapse_st(allele_st_tabs[["g001303"]], min_n = 2)
# long_g000222 <- collapse_st(allele_st_tabs[["g000222"]], min_n = 2)
# long_g000145 <- collapse_st(allele_st_tabs[["g000145"]], min_n = 2)
# long_g003299 <- collapse_st(allele_st_tabs[["g003299"]], min_n = 2)
# long_g001254 <- collapse_st(allele_st_tabs[["g001254"]], min_n = 2)
# long_g004223 <- collapse_st(allele_st_tabs[["g004223"]], min_n = 2)
# 


# Look at amino acids -----------------------------------------------------
# If the hypothesis is functional or antigenic, peptide alleles
# are the right unit and nucleotide alleles inflate counts with synonymous variation.


long_g000371 <- collapse_st(ag_st_allele_counts[["g000371"]], min_n = 2)
long_g001257 <- collapse_st(ag_st_allele_counts[["g001257"]], min_n = 2)
long_g001303 <- collapse_st(ag_st_allele_counts[["g001303"]], min_n = 2)
long_g000222 <- collapse_st(ag_st_allele_counts[["g000222"]], min_n = 2)
long_g000145 <- collapse_st(ag_st_allele_counts[["g000145"]], min_n = 2)
long_g003299 <- collapse_st(ag_st_allele_counts[["g003299"]], min_n = 2)
long_g001254 <- collapse_st(ag_st_allele_counts[["g001254"]], min_n = 2)
long_g004223 <- collapse_st(ag_st_allele_counts[["g004223"]], min_n = 2)




plot_panelA <- function(long, gname, focal = "ST307", seedy = 21) {
  M <- as.matrix(dcast(long, allele_name ~ ST_grp, value.var = "n", fill = 0),
                 rownames = "allele_name")
  # p1 : p3 freq
  p_n = sort(rowSums(M)/sum(M), decreasing = TRUE)[1:3]
  p_n <- sprintf("%1.2f", p_n)
  
  n_per_st <- colSums(M)
  nb  <- length(n_per_st)
  
  max_shown = nrow(M)
  
  ord   <- order(rowSums(M), decreasing = TRUE)
  top   <- rownames(M)[ord][seq_len(min(max_shown, nrow(M)))]
  Mtop  <- M[top, , drop = FALSE]
  other <- colSums(M) - colSums(Mtop)
  Mp    <- rbind(Mtop, other = other)
  Mp    <- sweep(Mp, 2, colSums(Mp), "/")
  
  pal <- c(hcl.colors(40, "Batlow"),
           hcl.colors(40, "Mako"))
  
  set.seed(seedy)
  cols <- sample(pal, nrow(M))
  
  # cols <- c(palette.colors(nrow(Mtop), "Alphabet"), "#D2D2D2")
  wid  <- sqrt(n_per_st) / max(sqrt(n_per_st))
  
  par(mar = c(3, 4.5, 3, 0.5))
  if (nb < 6) {
    bp <- barplot(Mp, col = cols, border = "white", lwd = 0.4,
                  width = wid, space = 0.35, 
                  xlim = c(0, (sum(wid) + nb * 0.35 * mean(wid)) * 2),
                  ylim = c(0, 1), yaxt = "n", names.arg = rep("", nb))
  } else {
    bp <- barplot(Mp, col = cols, border = "white", width = wid,
                  space = 0.35, ylim = c(0, 1), yaxt = "n",
                  names.arg = rep("", ncol(Mp)))
  }
  axis(2, at = c(0, .5, 1), las = 1)
  mtext("allele frequency", 2, line = 2.8, cex = 0.8)
  text(bp, -0.06, colnames(Mp), srt = 45, adj = 1, xpd = NA, cex = 0.8)
  text(bp, 1.02, paste0("n=", n_per_st), cex = 0.65, col = "grey35", xpd = NA)
  
  pre <- bquote(.(gname) ~ "  |  pan. freq =" ~ .(sum(M)) ~ "  |  " ~ .(nrow(M)) ~ "alleles  |  ")
  p1  <- bquote(italic(p)[1] ~ "=" ~ .(p_n[[1]]) * ",  ")
  p2  <- bquote(italic(p)[2] ~ "=" ~ .(p_n[[2]]) * ",  ")
  p3 <- bquote(italic(p)[3] ~ "=" ~ .(p_n[[3]]))
  
  mtext(bquote(.(pre) * phantom(.(p1) * .(p2) * .(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75)
  mtext(bquote(phantom(.(pre)) * .(p1) * phantom(.(p2) * .(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[1])
  mtext(bquote(phantom(.(pre) * .(p1)) * .(p2) * phantom(.(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[2])
  mtext(bquote(phantom(.(pre) * .(p1) * .(p2)) * .(p3)), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[3])
  
  # mtext(mtext(bquote(.(gname) ~ "  |  pan. freq =" ~ .(sum(M)) ~ "  |  " 
  #                    ~ .(nrow(M)) ~ "alleles  |  " 
  #                    ~ italic(p)[1]~"="~.(p_n[[1]])*","
  #                    ~ italic(p)[2]~"="~.(p_n[[2]])),
  #             side = 3, line = 1.5, adj = -0.1, cex = 0.75),
  #       3, line = 1.5, adj = -0.1, cex = 0.75)
}

png(paste0(outdir_fig,"/outlier_alleles_freqs_2.png"),
    width = 23.85, height = 41.5, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

par(mfrow = c(8,1))
plot_panelA(long_g003299, "AbiEi, Phage defence antitoxin")
plot_panelA(long_g001254, "AbiEii/AbiGii, Phage defence toxin")
plot_panelA(long_g000371, "kexD, Antibiotic efflux pump")
plot_panelA(long_g001257, "ramR, Efflux pump switch")
plot_panelA(long_g001303, "acrA, Antibiotic efflux pump")
plot_panelA(long_g000222, "acrB, Antibiotic efflux pump")
plot_panelA(long_g000145, "ompA, Outer membrane porin")
plot_panelA(long_g004223, "tssF, Bacterial harpoon gun")

dev.off()


# # balancing selection example
# selected_ags <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05 & expected_S < seg_s_sites]
# family_2d <- fread(paste0(outdir_dat, "/family_2d.csv"))
# selected_ags[gene_family %chin% family_2d[hgt_class=="none", gene_family]]




