# =========================================================================
# Single vs multiple acquisition filter: AG vs flanking core divergence ceiling
# =========================================================================

library(data.table)
library(Biostrings)
library(DECIPHER)
library(seqinr)
library(gUtils)
library(foreach)
library(doParallel)

set.seed(42)

fasta_dir       <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_260_chr_fasta/"
min_k           <- 5     # below this, report as untestable
n_calib         <- 200   # core genes used to calibrate tol
calib_n_genomes <- 30    # genomes sampled per calibration gene
tol_quantile    <- 0.95  # ~5% false positive rate on vertically inherited genes

# -------------------------------------------------------------------------
# 1. Annotation
# -------------------------------------------------------------------------
anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                     select = c("gene_family", "geno_id", "seqnames", "start",
                                "end", "strand", "anchor", "msu", "ag_type")))[msu != ""]

set1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

# -------------------------------------------------------------------------
# 2. Map each AG to its two flanking core gene families
# -------------------------------------------------------------------------
core_map <- unique(anno[ag_type == "core", .(msu, anchor, core_family = gene_family)])
if (anyDuplicated(core_map, by = c("msu", "anchor")))
  warning("Some core anchors map to more than one gene family; first match used.")
core_map <- core_map[, .SD[1], by = .(msu, anchor)]

# most common anchor per AG
ag_map <- anno[gene_family %chin% set1_ags, .N, by = .(gene_family, msu, anchor)
][order(-N)][, .SD[1], by = gene_family]
ag_map[, c("anchor_L", "anchor_R") := tstrsplit(anchor, ":", fixed = TRUE, keep = 1:2)]

ag_map <- merge(ag_map, core_map[, .(msu, anchor_L = anchor, flank_L = core_family)],
                by = c("msu", "anchor_L"), all.x = TRUE)
ag_map <- merge(ag_map, core_map[, .(msu, anchor_R = anchor, flank_R = core_family)],
                by = c("msu", "anchor_R"), all.x = TRUE)

message(ag_map[is.na(flank_L) | is.na(flank_R), .N], " AGs have an unmapped flank")

# -------------------------------------------------------------------------
# 3. Calibration set: core gene triples (left neighbour, gene, right neighbour)
# -------------------------------------------------------------------------
ref_geno   <- anno[ag_type == "core", .N, by = geno_id][which.max(N), geno_id]
core_order <- anno[ag_type == "core" & geno_id == ref_geno][order(seqnames, start)]
core_order[, `:=`(flank_L = shift(gene_family),
                  flank_R = shift(gene_family, type = "lead")), by = seqnames]

calib_map <- core_order[!is.na(flank_L) & !is.na(flank_R)
][sample(.N, min(.N, n_calib)), .(gene_family, flank_L, flank_R)]

# -------------------------------------------------------------------------
# 4. Coordinates for every family needed (single-copy per genome only)
# -------------------------------------------------------------------------
needed <- unique(na.omit(c(ag_map$gene_family, ag_map$flank_L, ag_map$flank_R,
                           unlist(calib_map))))

coords <- unique(anno[gene_family %chin% needed,
                      .(gene_family, geno_id, seqnames, start, end, strand)])
coords[, n_copy := .N, by = .(gene_family, geno_id)]
message(coords[n_copy > 1, uniqueN(paste(gene_family, geno_id))],
        " family-genome combinations dropped as multi-copy")
coords <- coords[n_copy == 1][, n_copy := NULL]
coords[, id := paste(gene_family, geno_id, sep = "|")]

coords_by_geno <- split(coords, by = "geno_id")

# -------------------------------------------------------------------------
# 5. Extract sequences (strand-aware) from each genome
# -------------------------------------------------------------------------
fna_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)

cl <- makeCluster(num_cores)
registerDoParallel(cl)

seq_list <- foreach(f = fna_files,
                    .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                      gid <- sub("\\.fasta$", "", basename(f))
                      dt  <- coords_by_geno[[gid]]
                      if (is.null(dt)) return(DNAStringSet())
                      genome <- readDNAStringSet(f)
                      s <- getSeq(genome, dt2gr(dt[, .(seqnames, start, end, strand)]))
                      names(s) <- dt$id
                      s
                    }

stopCluster(cl)

all_seqs <- do.call(c, unname(seq_list))
fam      <- sub("\\|.*", "", names(all_seqs))
names(all_seqs) <- sub(".*\\|", "", names(all_seqs))   # names become geno_id
seqs_by_fam <- split(all_seqs, fam)

rm(seq_list, all_seqs); gc()

# -------------------------------------------------------------------------
# 6. Helpers
# -------------------------------------------------------------------------

# codon-aware alignment; drops sequences not in frame
codon_align <- function(s) {
  s <- s[width(s) %% 3 == 0]
  if (length(s) < 2) return(NULL)
  AlignTranslation(s, verbose = FALSE)
}

# remove codons with gaps, ambiguous bases, or stop codons in any sequence
clean_codons <- function(aln) {
  m     <- as.matrix(aln)
  n_cod <- ncol(m) %/% 3
  if (n_cod == 0) return(NULL)
  m     <- m[, seq_len(n_cod * 3), drop = FALSE]
  
  ok_base <- matrix(m %in% c("A", "C", "G", "T"), nrow = nrow(m))
  cod_ok  <- colSums(!ok_base[, c(TRUE, FALSE, FALSE), drop = FALSE] |
                       !ok_base[, c(FALSE, TRUE, FALSE), drop = FALSE] |
                       !ok_base[, c(FALSE, FALSE, TRUE), drop = FALSE]) == 0
  
  cods <- matrix(paste0(m[, c(TRUE, FALSE, FALSE)],
                        m[, c(FALSE, TRUE, FALSE)],
                        m[, c(FALSE, FALSE, TRUE)]), nrow = nrow(m))
  no_stop <- colSums(matrix(cods %in% c("TAA", "TAG", "TGA"), nrow = nrow(m))) == 0
  
  keep <- rep(cod_ok & no_stop, each = 3)
  if (!any(keep)) return(NULL)
  
  out <- apply(m[, keep, drop = FALSE], 1, paste0, collapse = "")
  names(out) <- rownames(m)
  out
}

# pairwise Ks matrix (Li 1993 via seqinr); saturated/undefined -> NA
ks_matrix <- function(seqs) {
  a <- seqinr::as.alignment(nb = length(seqs), nam = names(seqs), seq = tolower(seqs))
  k <- as.matrix(seqinr::kaks(a)$ks)
  k[k >= 10] <- NA
  k[k < 0]   <- 0
  k
}
# ceiling test for one focal gene against its two flanks
run_one <- function(focal, fL, fR, tol = 1, genomes = NULL) {
  
  s_f <- seqs_by_fam[[focal]]; s_L <- seqs_by_fam[[fL]]; s_R <- seqs_by_fam[[fR]]
  if (is.null(s_f) || is.null(s_L) || is.null(s_R)) return(NULL)
  
  g <- Reduce(intersect, list(names(s_f), names(s_L), names(s_R)))
  if (!is.null(genomes)) g <- intersect(g, genomes)
  if (length(g) < 2) return(NULL)
  
  a_f <- codon_align(s_f[g]); a_L <- codon_align(s_L[g]); a_R <- codon_align(s_R[g])
  if (is.null(a_f) || is.null(a_L) || is.null(a_R)) return(NULL)
  
  g <- Reduce(intersect, list(names(a_f), names(a_L), names(a_R)))
  if (length(g) < 2) return(NULL)
  
  c_f <- clean_codons(a_f[g])
  c_L <- clean_codons(a_L[g])
  c_R <- clean_codons(a_R[g])
  if (is.null(c_f) || is.null(c_L) || is.null(c_R)) return(NULL)
  
  c_fl <- paste0(c_L[g], c_R[g]); names(c_fl) <- g
  c_f  <- c_f[g]
  
  d_ag <- ks_matrix(c_f)
  d_fl <- ks_matrix(c_fl)
  ut   <- upper.tri(d_ag)
  
  if (all(is.na(d_fl[ut]))) return(NULL)
  
  # floor of ~one synonymous difference, so clonal carriers don't give a zero ceiling
  floor_ks <- 1 / (0.25 * nchar(c_f[1]))
  max_fl   <- max(d_fl[ut], na.rm = TRUE)
  ref      <- max(max_fl, floor_ks)
  ceil     <- ref * tol
  
  d_ag[is.na(d_ag)] <- Inf          # saturated AG pairs exceed any ceiling
  diag(d_ag) <- 0
  max_ag <- max(d_ag[ut])
  
  d_cl <- d_ag
  d_cl[!is.finite(d_cl)] <- 1e6
  n_acq <- length(unique(cutree(hclust(as.dist(d_cl), method = "single"), h = ceil)))
  
  data.table(gene_family = focal, k = length(g),
             cod_used = nchar(c_f[1]) / 3,
             max_d_ag = max_ag, max_d_fl = max_fl,
             R_max = max_ag / ref, n_acq = n_acq)
}

# -------------------------------------------------------------------------
# 7. Calibrate tol on core genes vs their own neighbours
# -------------------------------------------------------------------------
all_genomes <- unique(coords$geno_id)

calib_res <- rbindlist(lapply(seq_len(nrow(calib_map)), function(i) {
  run_one(calib_map$gene_family[i], calib_map$flank_L[i], calib_map$flank_R[i],
          genomes = sample(all_genomes, calib_n_genomes))
}))

tol <- as.numeric(quantile(calib_res[is.finite(R_max), R_max], tol_quantile))
message("Calibrated tol = ", round(tol, 3), " from ", nrow(calib_res), " core genes")
# Calibrated tol = 2.809 from 200 core genes

hist(calib_res[is.finite(R_max), R_max], breaks = 50,
     main = "Core genes vs own flanks", xlab = expression(R[max]))
abline(v = tol, col = "red", lwd = 2)

# -------------------------------------------------------------------------
# 8. Run on AGs
# -------------------------------------------------------------------------
ag_run <- ag_map[!is.na(flank_L) & !is.na(flank_R)]

cl <- makeCluster(num_cores)
registerDoParallel(cl)

acq_res <- foreach(i = seq_len(nrow(ag_run)), .combine = rbind,
                   .packages = c("data.table", "Biostrings", "DECIPHER", "seqinr"),
                   .export = c("seqs_by_fam", "run_one", "codon_align",
                               "clean_codons", "ks_matrix", "tol")) %dopar% {
                                 run_one(ag_run$gene_family[i], ag_run$flank_L[i], ag_run$flank_R[i], tol = tol)
                               }

stopCluster(cl)

# -------------------------------------------------------------------------
# 9. Classify and save
# -------------------------------------------------------------------------
acq_res <- merge(ag_map[, .(gene_family, msu, anchor, flank_L, flank_R)],
                 acq_res, by = "gene_family", all.x = TRUE)

acq_res[, status := fcase(is.na(flank_L) | is.na(flank_R), "no_flank",
                          is.na(n_acq),                    "failed",
                          k < min_k,                       "untestable",
                          n_acq == 1,                      "single",
                          default =                        "multiple")]

print(acq_res[, .N, by = status])


fwrite(acq_res,   paste0(outdir_dat, "/ag_acquisition_ceiling.csv"))
fwrite(calib_res, paste0(outdir_dat, "/ag_acquisition_calibration.csv"))

single_acq_ags <- acq_res[status == "single", gene_family]
