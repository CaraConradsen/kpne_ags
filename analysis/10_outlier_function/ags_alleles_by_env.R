


anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                          select = c("fus_locus_tag","geno_id", "gene_family", "ST", "threshold", 
                                     "freq", "consensus_gene_name","consensus_product",
                                     "alleles_at_maximum_threshold")))
# anno <- anno[gene_family %chin% c(#"g000371", "g000222", "g001303", "g001257", 
#                                   #"g000145", "g001254", "g003299", "g004223", 
#                                   # new fellas 
#                                   "g001247","g004140","g004984", "g003175", "g000242_2","g002556")]

# Focal AGs: robust outliers (see AG outlier candidates doc, 2026-10-14).
# Labels use corrected gene names, not PIRATE consensus names.
focal_dt <- data.table(
  gene_family = c("g004051", "g003783", "g001147", "g000322",  # O1/O2 O-antigen locus (A308:A309)
                  "g003390", "g003914",                        # T6SS (A164:A165)
                  "g003299", "g001254",                        # AbiE toxin-antitoxin operon (A22:A23)
                  "g002556",                                   # acid resistance
                  "g001058",                                   # chromate efflux
                  "g002655", "g003239",                        # acetoin catabolism operon (A169:A170)
                  "g002664", "g002490",                        # ectoine degradation (A68:A69)
                  "g003377", "g004320",                        # sorbose/sorbitol (A7:A8)
                  "g002126_1",                                 # glycerol-3-phosphate uptake
                  "g015591",                                   # ribosomal protein L18
                  "g006452","g003974"),                                  #
  label = c("wzm, O-antigen export permease",
            "wzt, O-antigen export ATP-binding protein",
            "glf, UDP-galactopyranose mutase (O-antigen)",
            "wbbM, O-antigen glycosyltransferase (DUF4422)",
            "vgrG, T6SS tip protein",
            "tssA, T6SS assembly protein",
            "abiEi, Phage defence antitoxin",
            "abiEii, Phage defence toxin",
            "asr, Acid resistance protein",
            "chrA, Chromate transporter",
            "acoA, Acetoin dehydrogenase E1 alpha",
            "acoB, Acetoin dehydrogenase E1 beta",
            "phyH, Ectoine dioxygenase",
            "phyH, Phytanoyl-CoA dioxygenase",
            "sorF, Sorbose PTS EIIA",
            "sorC, Sorbitol operon regulator",
            "glpT, Glycerol-3-phosphate transporter (fosfomycin uptake)",
            "rplR, 50S ribosomal protein L18",
            "ttdA, Fumarate hydratase",
            "tdh, NAD(P)-dependent alcohol dehydrogenase"))

anno <- anno[gene_family %chin% focal_dt$gene_family]


unique(anno[,.(gene_family, freq, consensus_gene_name, consensus_product, alleles_at_maximum_threshold)])


# Using aa ----------------------------------------------------------------
# accounting for fusions
faa_files <-  list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_260_hybrid_chr_hps0.6_out/feature_sequences/",
                         full.names = TRUE, pattern = ".aa.fasta")

faa_files <- faa_files[grepl(paste(unique(anno$gene_family), collapse = "|"), faa_files)]

anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                     select = c("fus_locus_tag","geno_id", "gene_family", "threshold", 
                                "freq", "consensus_gene_name","consensus_product",
                                "seqnames", "fstart",  "fend", "fstrand",
                                "alleles_at_maximum_threshold")))

anno <- anno[gene_family %chin% focal_dt$gene_family]

# add environmental info
meta <- fread("../SpARK data/spark_metadata.csv", select = c("id", "group_summary"))
setnames(meta, c("id","group_summary"), c("geno_id", "env"))

anno <- merge(anno, meta, 
              all.x = TRUE, by = "geno_id")

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
ag_env_allele_counts <- foreach(g = focal_ags, 
                               .combine = "c", 
                               .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                                 
                                 g_dat <- anno[gene_family == g, .(geno_id, seqnames, start, end, strand, env)]
                                 
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
                                 
                                 list(dcast(g_dat, allele_name ~ env, fun.aggregate = length, value.var = "geno_id"))
                                 
                               }

# Stop cluster after execution
stopCluster(cl)

names(ag_env_allele_counts) <- focal_ags


# Allele diversity plot ---------------------------------------------------
# 
collapse_env <- function(tab, min_n = 10, merge_lv = TRUE) {
  long <- melt(tab, id.vars = "allele_name",
               variable.name = "env", value.name = "n")
  long[, env := as.character(env)]
  if (merge_lv) long[, env := sub("-[0-9]+LV$", "", env)]
  long <- long[, .(n = sum(n)), by = .(allele_name, env)]

  tot  <- long[, .(N = sum(n)), by = env][order(-N)]
  keep <- tot[N >= min_n, env]
  long[, env_grp := fifelse(env %in% keep, env, "other")]

  out <- long[, .(n = sum(n)), by = .(allele_name, env_grp)]
  lev <- c(tot[env %in% keep][order(-N), env], "other")
  out[, env_grp := factor(env_grp, levels = lev)]
  attr(out, "n_other_st") <- tot[!env %in% keep, .N]
  out[]
}

# Look at amino acids -----------------------------------------------------
# If the hypothesis is functional or antigenic, peptide alleles
# are the right unit and nucleotide alleles inflate counts with synonymous variation.

long_focal <- lapply(ag_env_allele_counts[focal_dt$gene_family],
                     collapse_env, min_n = 1)
 

plot_panelA <- function(long, gname, focal = "ST307", seedy = 11) {
  M <- as.matrix(dcast(long, allele_name ~ env_grp, value.var = "n", fill = 0),
                 rownames = "allele_name")
  # p1 : p3 freq
  p_n = sort(rowSums(M)/sum(M), decreasing = TRUE)[1:3]
  p_n <- sprintf("%1.2f", p_n)
  
  n_per_env <- colSums(M)
  nb  <- length(n_per_env)
  
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
  wid  <- sqrt(n_per_env) / max(sqrt(n_per_env))
  
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
  text(bp, 1.02, paste0("n=", n_per_env), cex = 0.65, col = "grey35", xpd = NA)
  
  pre <- bquote(.(gname) ~ "  |  pan. freq =" ~ .(sum(M)) ~ "  |  " ~ .(nrow(M)) ~ "alleles  |  ")
  p1  <- bquote(italic(p)[1] ~ "=" ~ .(p_n[[1]]) * ",  ")
  p2  <- bquote(italic(p)[2] ~ "=" ~ .(p_n[[2]]) * ",  ")
  p3 <- bquote(italic(p)[3] ~ "=" ~ .(p_n[[3]]))
  
  mtext(bquote(.(pre) * phantom(.(p1) * .(p2) * .(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75)
  mtext(bquote(phantom(.(pre)) * .(p1) * phantom(.(p2) * .(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[1])
  mtext(bquote(phantom(.(pre) * .(p1)) * .(p2) * phantom(.(p3))), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[2])
  mtext(bquote(phantom(.(pre) * .(p1) * .(p2)) * .(p3)), 3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[3])
  
}

# # Interesting candidates
# # Interesting candidates: 2 pages x 9 panels -> _1.png, _2.png
png(paste0(outdir_fig, "/outlier_alleles_freqs_vary_env_%d.png"),
    width = 23.85, height = 46.7, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

par(mfrow = c(5,1))
for (i in seq_len(nrow(focal_dt))) {
  plot_panelA(long_focal[[focal_dt$gene_family[i]]], focal_dt$label[i])
}

dev.off()


# ChrA --------------------------------------------------------------------




# # Interesting candidates
# # Interesting candidates: 2 pages x 9 panels -> _1.png, _2.png
png(paste0(outdir_fig, "/ChrA_vary_env.png"),
    width = 23.85, height = 11.675, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

M <- as.matrix(dcast(long_focal[["g001058"]], allele_name ~ env_grp, value.var = "n", fill = 0),
               rownames = "allele_name")

env_levels <- c("human - hospital carriage", "human - hospital disease", "human - community carriage",
                "surface - farm", "water - farm", "water - environment",
                "cow", "pig", "dog", "turtle", "rhizoplane", "rhizosphere")

dropped <- setdiff(colnames(M), env_levels)
if (length(dropped)) message("not plotted: ", paste(dropped, collapse = ", "),
                             " (", sum(M[, dropped]), " genomes)")

M <- M[, intersect(env_levels, colnames(M)), drop = FALSE]
M <- M[rowSums(M) > 0, , drop = FALSE]   # alleles only seen in dropped groups

# p1 : p3 freq
p_n = sort(rowSums(M)/sum(M), decreasing = TRUE)[1:6]
p_n <- sprintf("%1.2f", p_n)

n_per_env <- colSums(M)
nb  <- length(n_per_env)

max_shown = nrow(M)

ord   <- order(rowSums(M), decreasing = TRUE)
top   <- rownames(M)[ord][seq_len(min(max_shown, nrow(M)))]
Mtop  <- M[top, , drop = FALSE]
other <- colSums(M) - colSums(Mtop)
Mp    <- rbind(Mtop, other = other)
Mp    <- sweep(Mp, 2, colSums(Mp), "/")

pal  <- palette.colors(NULL, "Polychrome 36")[-(1:2)]   # drop its dark grey and off-white
cols <- c("#201158",unname(pal[seq_len(nrow(M))]))

# cols <- c(palette.colors(nrow(Mtop), "Alphabet"), "#D2D2D2")
wid  <- sqrt(n_per_env) / max(sqrt(n_per_env))

par(mar = c(8, 4.5, 3, 0.5))
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
text(bp, -0.06, colnames(Mp), srt = 30, adj = 1, xpd = NA, cex = 0.8)
text(bp, 1.02, paste0("n=", n_per_env), cex = 0.65, col = "grey35", xpd = NA)

pre <- bquote(.(focal_dt$label[focal_dt$gene_family =="g001058"]) ~ "  |  pan. freq =" ~ .(sum(M)) ~ "  |  " ~ .(nrow(M)) ~ "alleles  |  ")
p1  <- bquote(italic(p)[1] ~ "=" ~ .(p_n[[1]]) * ",  ")
p2  <- bquote(italic(p)[2] ~ "=" ~ .(p_n[[2]]) * ",  ")
p3  <- bquote(italic(p)[3] ~ "=" ~ .(p_n[[3]]) * ",  ")
p4  <- bquote(italic(p)[4] ~ "=" ~ .(p_n[[4]]) * ",  ")
p5  <- bquote(italic(p)[5] ~ "=" ~ .(p_n[[5]]) * ",  ")
p6 <- bquote(italic(p)[6] ~ "=" ~ .(p_n[[6]]))

mtext(bquote(.(pre) * phantom(.(p1) * .(p2) * .(p3) * .(p4) * .(p5) * .(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75)
mtext(bquote(phantom(.(pre)) * .(p1) * phantom(.(p2) * .(p3) * .(p4) * .(p5) * .(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[1])
mtext(bquote(phantom(.(pre) * .(p1)) * .(p2) * phantom(.(p3) * .(p4) * .(p5) * .(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[2])
mtext(bquote(phantom(.(pre) * .(p1) * .(p2)) * .(p3) * phantom(.(p4) * .(p5) * .(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[3])
mtext(bquote(phantom(.(pre) * .(p1) * .(p2) * .(p3)) * .(p4) * phantom(.(p5) * .(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[4])
mtext(bquote(phantom(.(pre) * .(p1) * .(p2) * .(p3) * .(p4)) * .(p5) * phantom(.(p6))),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[5])
mtext(bquote(phantom(.(pre) * .(p1) * .(p2) * .(p3) * .(p4) * .(p5)) * .(p6)),
      3, line = 1.5, adj = -0.1, cex = 0.75, col = cols[6])

dev.off()

