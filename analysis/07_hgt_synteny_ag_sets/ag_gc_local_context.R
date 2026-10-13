# Local GC context of AGs, compared only within the genomes that carry them
#
# For every AG copy (set 1, syntenic, one insertion site per family):
#   d_local  = GC3(AG) - GC3(flanking core anchors, same genome, same MSU)
#   d_genome = GC3(AG) - GC3(all core anchors, same genome)
#
# Family metrics:
#   offset_local : median d_local -> foreign composition (single, possibly old, acquisition)
#   het_obs      : sd of d_local across carriers -> heterogeneity among copies
#                  (multiple introductions / allelic replacement at the same locus)
#
# Null for het_obs: the same statistic for core anchor genes of similar length,
# computed over exactly the same carrier genomes. Core genes are vertically
# inherited, so this is the between-copy GC3 variation expected from descent
# through the clonal frame of those genomes alone.
#
# Replaces the "Get GC content" and "Calculate AG GC zscores" sections.

# settings ----------------------------------------------------------------
flank_k  <- 5   # core anchors taken either side of the AG insertion
n_len_bins <- 5 # gene-length bins for matching AGs to core null genes

set1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

fna_files <-  list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_260_chr_fasta/",
                         full.names = TRUE, pattern = ".fasta")

# get locations -----------------------------------------------------------
full_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                   select = c("gene_family", "geno_id", "locus_tag", "seqnames",
                              "start", "end", "strand",
                              "anchor", "msu", "ag_type"))[msu != ""]

# one row per locus (guards against duplicated MGE merges)
full_anno <- unique(full_anno, by = c("geno_id", "locus_tag"))

# core anchors: one gene per genome per anchor
core_anch <- full_anno[!grepl(":", anchor)]
core_anch[, a_idx := as.integer(sub("^A", "", anchor))]

# set 1 AGs, with the anchors either side of the insertion
ag_loc <- full_anno[gene_family %chin% set1_ags & grepl(":", anchor)]
ag_loc[, c("a_left", "a_right") := tstrsplit(anchor, ":", fixed = TRUE)]
ag_loc[, c("a_left", "a_right") := lapply(.SD, function(x) as.integer(sub("^A", "", x))),
       .SDcols = c("a_left", "a_right")]
ag_loc[, `:=`(lo = pmin(a_left, a_right) - flank_k + 1L,
              hi = pmax(a_left, a_right) + flank_k - 1L)]


# Step 1, GC and GC3 for every core anchor and AG locus --------------------
seq_dt <- rbind(core_anch[, .(geno_id, locus_tag, seqnames, start, end, strand)],
                ag_loc[,    .(geno_id, locus_tag, seqnames, start, end, strand)])

list_seqs <- split(seq_dt, by = "geno_id")

# Set up parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)

locus_gc <- foreach(f = fna_files,
                    .combine = rbind,
                    .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {

                      genome_id = gsub(".fasta", "", basename(f))

                      # GC at third codon positions (CDS is in frame from start, strand aware)
                      gc3_fun <- function(x) {
                        vapply(seq_along(x), function(k) {
                          s   <- x[[k]]
                          idx <- seq(3L, length(s) - length(s) %% 3L, by = 3L)
                          letterFrequency(s[idx], letters = "GC", as.prob = TRUE)[[1]]
                        }, numeric(1))
                      }

                      genome <- readDNAStringSet(f)

                      temp_dt <- list_seqs[[genome_id]]
                      if (is.null(temp_dt)) return(NULL)

                      genes_gr <- dt2gr(temp_dt[, .(seqnames, start, end, strand)])

                      temp_seqs <- getSeq(genome, genes_gr)

                      data.table(geno_id   = genome_id,
                                 locus_tag = temp_dt$locus_tag,
                                 len       = width(temp_seqs),
                                 gc        = letterFrequency(temp_seqs, letters = "GC", as.prob = TRUE)[, 1],
                                 gc3       = gc3_fun(temp_seqs))
                    }

# Stop cluster after execution
stopCluster(cl)

core_anch <- merge(core_anch, locus_gc, by = c("geno_id", "locus_tag"))
ag_loc    <- merge(ag_loc,    locus_gc, by = c("geno_id", "locus_tag"))


# Step 2, genome-wide core GC per genome -----------------------------------
genome_core <- core_anch[, .(genome_gc  = weighted.mean(gc,  len),
                             genome_gc3 = weighted.mean(gc3, len)), by = geno_id]


# Step 3, flanking core GC per AG copy (same genome, same MSU) -------------
flank_dt <- core_anch[ag_loc,
                      on = .(geno_id, msu, a_idx >= lo, a_idx <= hi),
                      .(locus_tag = i.locus_tag, gc = x.gc, gc3 = x.gc3, len = x.len),
                      allow.cartesian = TRUE, nomatch = NULL]

flank_dt <- flank_dt[, .(flank_gc  = weighted.mean(gc,  len),
                         flank_gc3 = weighted.mean(gc3, len),
                         n_flank   = .N), by = locus_tag]

ag_comp <- merge(ag_loc[, .(gene_family, geno_id, msu, locus_tag, len, gc, gc3)],
                 flank_dt, by = "locus_tag")

ag_comp <- merge(ag_comp, genome_core, by = "geno_id")

ag_comp[, `:=`(d_local         = gc3 - flank_gc3,
               d_genome        = gc3 - genome_gc3,
               flank_vs_genome = flank_gc3 - genome_gc3)]


# Step 4, the same statistic for core anchors (null) -----------------------
# flanking window of each core anchor, excluding itself
core_anch[, `:=`(lo = a_idx - flank_k, hi = a_idx + flank_k)]

core_flank <- core_anch[core_anch,
                        on = .(geno_id, msu, a_idx >= lo, a_idx <= hi),
                        .(locus_tag = i.locus_tag, f_tag = x.locus_tag,
                          gc3 = x.gc3, len = x.len),
                        allow.cartesian = TRUE, nomatch = NULL][locus_tag != f_tag]

core_flank <- core_flank[, .(flank_gc3 = weighted.mean(gc3, len)), by = locus_tag]

core_anch <- merge(core_anch, core_flank, by = "locus_tag")
core_anch[, d_local := gc3 - flank_gc3]

# core family x genome matrix of d_local
core_mat <- dcast(core_anch, gene_family ~ geno_id, value.var = "d_local", fun.aggregate = mean)
core_fams <- core_mat$gene_family
core_mat  <- as.matrix(core_mat[, -1])
core_mat[is.nan(core_mat)] <- NA
rownames(core_mat) <- core_fams

# length bins
core_len <- core_anch[, .(len = as.numeric(median(len))), by = gene_family][match(core_fams, gene_family)]
len_breaks <- quantile(core_len$len, probs = seq(0, 1, length.out = n_len_bins + 1))
len_breaks[c(1, length(len_breaks))] <- c(0, Inf)
core_len[, len_bin := findInterval(len, len_breaks, rightmost.closed = TRUE)]


# Step 5, family metrics and test ------------------------------------------
fam_gc <- ag_comp[, .(freq            = .N,
                      len             = as.numeric(median(len)),
                      offset_local    = median(d_local),
                      offset_genome   = median(d_genome),
                      flank_vs_genome = median(flank_vs_genome),
                      het_obs         = sd(d_local),
                      carriers        = list(geno_id)), by = gene_family]

fam_gc[, len_bin := findInterval(len, len_breaks, rightmost.closed = TRUE)]

het_null <- foreach(i = seq_len(nrow(fam_gc)),
                    .combine = rbind,
                    .packages = c("data.table", "matrixStats")) %do% {

                      cols <- intersect(fam_gc$carriers[[i]], colnames(core_mat))
                      rows <- core_len[len_bin == fam_gc$len_bin[i], gene_family]

                      null_sd <- rowSds(core_mat[rows, cols, drop = FALSE], na.rm = TRUE)
                      null_sd <- null_sd[!is.na(null_sd)]

                      data.table(gene_family  = fam_gc$gene_family[i],
                                 het_null_med = median(null_sd),
                                 n_null       = length(null_sd),
                                 p_het        = (1 + sum(null_sd >= fam_gc$het_obs[i])) /
                                                (length(null_sd) + 1))
                    }

fam_gc <- merge(fam_gc, het_null, by = "gene_family")

fam_gc[, het_ratio := het_obs / het_null_med]
fam_gc[, p_bh_het  := p.adjust(p_het, method = "BH")]
fam_gc[, hgt_het   := p_bh_het <= 0.05]

fam_gc[, carriers := NULL]

table(fam_gc$hgt_het)


# output ------------------------------------------------------------------
fwrite(ag_comp, paste0(outdir_dat, "/ag_gc3_local_per_copy.csv"))
fwrite(fam_gc,  paste0(outdir_dat, "/ag_gc3_local_family.csv"))
# fam_gc  <- fread(paste0(outdir_dat, "/ag_gc3_local_family.csv"))
# ag_comp <- fread(paste0(outdir_dat, "/ag_gc3_local_per_copy.csv"))

col_map <- c("FALSE" = adjustcolor("steelblue", 0.3),
             "TRUE"  = adjustcolor("firebrick", 0.8))
cex_f   <- function(f) sqrt(f / max(fam_gc$freq)) * 2 + 0.3

par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))

# A. observed vs expected heterogeneity -----------------------------------
eps <- 1e-4
with(fam_gc[tested == TRUE],
     plot(het_null_med + eps, het_obs + eps, log = "xy",
          pch = 16, cex = cex_f(freq), col = col_map[as.character(hgt_het)],
          xlab = expression("Core null median SD(" * Delta * GC3[local] * ")"),
          ylab = expression("AG SD(" * Delta * GC3[local] * ")"),
          main = "A. Between-copy GC3 heterogeneity"))
abline(0, 1, lty = 2, col = "grey40")
legend("topleft", bty = "n", pch = 16, cex = 0.8,
       col = col_map[c("TRUE", "FALSE")],
       legend = c(paste0("HGT (", sum(fam_gc$hgt_het), ")"),
                  paste0("not flagged (", sum(fam_gc$tested & !fam_gc$hgt_het), ")")))

# B. foreignness vs heterogeneity -----------------------------------------
with(fam_gc[tested == TRUE],
     plot(offset_local, z_het,
          pch = 16, cex = cex_f(freq), col = col_map[as.character(hgt_het)],
          xlab = expression("Median " * Delta * GC3[local] * " (AG - flanking core)"),
          ylab = expression(z[het]),
          main = "B. Single foreign acquisition vs multiple introductions"))
abline(v = 0, lty = 3, col = "grey40")
abline(h = qnorm(1 - 0.05 / sum(fam_gc$tested)), lty = 2, col = "grey40")  # Bonferroni line for reference

# C. p-value calibration ----------------------------------------------------
hist(fam_gc[tested == TRUE, p_het], breaks = 40,
     col = "grey70", border = "grey40",
     xlab = expression(p[het]), main = "C. p-value distribution (tested families)")
abline(h = sum(fam_gc$tested) / 40, lty = 2, col = "firebrick")   # uniform expectation

# D. per-copy values for the top flagged families -------------------------
top_fams <- fam_gc[hgt_het == TRUE][order(-z_het)][1:min(.N, 12), gene_family]
ref_fams <- fam_gc[tested == TRUE & hgt_het == FALSE][order(abs(z_het))][1:3, gene_family]
plot_fams <- c(top_fams, ref_fams)

strip_dt <- ag_comp[gene_family %chin% plot_fams]
strip_dt[, gene_family := factor(gene_family, levels = plot_fams)]

stripchart(d_local ~ gene_family, data = strip_dt,
           vertical = TRUE, method = "jitter", jitter = 0.2,
           pch = 16, cex = 0.6,
           col = c(rep(adjustcolor("firebrick", 0.6), length(top_fams)),
                   rep(adjustcolor("steelblue", 0.6), length(ref_fams))),
           las = 2, cex.axis = 0.7,
           ylab = expression(Delta * GC3[local] ~ "per copy"),
           main = "D. Copies within flagged (red) vs typical (blue) families")
abline(h = 0, lty = 3, col = "grey40")

par(mfrow = c(1, 1))


abs_d  <- abs(delta$delta)
sign_d <- ifelse(delta$delta < 0, "negative", "positive")

breaks <- seq(0, max(abs_d, na.rm = TRUE) + 0.001, by = 0.001)

counts <- rbind(
  negative = hist(abs_d[sign_d == "negative"], breaks = breaks, plot = FALSE)$counts,
  positive = hist(abs_d[sign_d == "positive"], breaks = breaks, plot = FALSE)$counts
)

bp <- barplot(counts, space = 0, border = NA,
              col = c(adjustcolor("dodgerblue3", 0.8), adjustcolor("firebrick", 0.8)),
              xlab = expression("|" * delta * "| (|AG " * pi[S] * " - flanking core " * pi[S] * "|)"),
              ylab = "Frequency",
              main = expression("Folded distribution of " * delta))

# x-axis on the |delta| scale (barplot positions bars at 0.5, 1.5, ...)
at_vals <- pretty(breaks)
axis(1, at = at_vals / diff(breaks[1:2]), labels = at_vals)

legend("topright", bty = "n", fill = c(adjustcolor("dodgerblue3", 0.8), adjustcolor("firebrick", 0.8)),
       border = NA, legend = c(paste0("AG < core (", sum(sign_d == "negative"), ")"),
                               paste0("AG > core (", sum(sign_d == "positive"), ")")))
