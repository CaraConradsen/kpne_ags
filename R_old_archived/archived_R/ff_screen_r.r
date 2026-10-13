## =====================================================================
## ff_screen.R
## Screen PIRATE AG alignments (ag_fus_alns) for fission/fusion artefacts
## and test whether flagged families are enriched for segregating sites.
## Dependencies: data.table only (plots in base R).
## =====================================================================

library(data.table)

## ---- settings: EDIT -------------------------------------------------
base_dir  <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data"
aln_dir   <- file.path(base_dir, "ag_fus_alns")
anno_file <- file.path(aln_dir, "kpne_anno_full_3490_ags.csv")
test_file <- file.path(aln_dir, "test_results.csv")
out_dir   <- file.path(base_dir, "ff_screen")
dir.create(out_dir, showWarnings = FALSE)

set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

# anno column names (check names(anno))
col_geno   <- "geno_id"
col_family <- "gene_family"
col_ftag   <- "fus_locus_tag"   # ":" or ";" in this column = fission/fusion called

# test_results column names (check names(test))
col_test_family <- "gene_family"
col_freq    <- "freq"
col_seg     <- "seg_sites"
col_len     <- "gene_length"
col_pbh     <- "p_bh"

len_hi      <- 1.2   # > this x family median ungapped length: possible unflagged fusion
len_lo      <- 0.8   # < this: truncated / fragment
priv_mult   <- 5     # private-site outlier: > priv_mult x family median (min priv_min)
priv_min    <- 5
sentinel    <- "g000145"

stops <- c("TAA", "TAG", "TGA")
bases <- c("A", "C", "G", "T")

## ---- helpers ----------------------------------------------------------
read_fasta <- function(f) {
  x <- sub("\r$", "", readLines(f, warn = FALSE))
  x <- x[nzchar(x)]
  h <- grepl("^>", x)
  id <- sub("^>", "", x[h])
  grp <- factor(cumsum(h)[!h], levels = seq_along(id))
  s <- vapply(split(x[!h], grp), paste, "", collapse = "")
  setNames(toupper(s), id)
}

codons <- function(s) {
  n <- nchar(s) %/% 3
  if (n == 0) return(character(0))
  substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
}

seq_metrics <- function(s) {
  ug   <- gsub("[-.!?]", "", s)
  L    <- nchar(ug)
  cd   <- codons(ug)
  nc   <- length(cd)
  core <- sub("^-+", "", sub("-+$", "", s))
  runs <- nchar(regmatches(core, gregexpr("-+", core))[[1]])
  lead <- nchar(sub("^(-*).*$", "\\1", s))
  ac   <- codons(s)
  ok   <- !grepl("[-.!?]", ac)
  last <- if (any(ok)) max(which(ok)) else 0
  list(
    aln_len      = nchar(s),
    ungap_len    = L,
    gap_frac     = 1 - L / nchar(s),
    len_mod3     = L %% 3,
    lead_gap     = lead,
    n_bad_gap    = sum(runs %% 3 != 0),          # internal gap runs not multiple of 3
    int_stop     = if (nc > 1) sum(cd[-nc] %in% stops) else 0L,
    int_stop_aln = sum(ac[ok & seq_along(ac) < last] %in% stops),
    start_ok     = substr(ug, 1, 3) %in% c("ATG", "GTG", "TTG"),
    stop_ok      = nc > 0 && cd[nc] %in% stops,
    n_amb        = nchar(gsub("[ACGT.!?-]", "", s)),
    n_bang       = nchar(s) - nchar(gsub("!", "", s, fixed = TRUE))
  )
}

# Private (singleton) segregating sites per sequence, split by codon position.
# Artefactual fragments add private sites roughly evenly across positions;
# real variation is concentrated at position 3.
private_sites <- function(seqs) {
  m <- do.call(rbind, strsplit(seqs, ""))
  n <- nrow(m)
  pos <- ((seq_len(ncol(m)) - 1) %% 3) + 1
  nvalid <- colSums(matrix(m %in% bases, n))
  ndist  <- Reduce(`+`, lapply(bases, function(b) colSums(m == b) > 0))
  keep   <- nvalid >= 3 & ndist >= 2
  hits   <- matrix(FALSE, n, ncol(m))
  for (b in bases) {
    mb <- m == b
    single <- colSums(mb) == 1 & keep
    hits <- hits | (mb & rep(single, each = n))
  }
  sapply(1:3, function(p) rowSums(hits[, pos == p, drop = FALSE]))
}

## ---- 1. read alignments -----------------------------------------------
files <- list.files(aln_dir, "\\.fasta$", full.names = TRUE)
message(length(files), " alignment files")

# look at outliers only:
focal_files <- lapply(set_1_ags, function(x){
  paste0("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/ag_fus_alns/",x,".fasta")
})

focal_files <- unlist(focal_files)
files[files %in% focal_files] -> files

res <- rbindlist(lapply(files, function(f) {
  s   <- read_fasta(f)
  fam <- sub("\\.fasta$", "", basename(f))
  aligned <- length(unique(nchar(s))) == 1
  met <- rbindlist(lapply(s, seq_metrics))
  priv <- if (aligned && length(s) >= 3) private_sites(s) else
    matrix(NA_integer_, length(s), 3)
  colnames(priv) <- c("priv_p1", "priv_p2", "priv_p3")
  data.table(gene_family = fam,
             header      = names(s),
             geno_id     = sub("_[^_]+$", "", names(s)),
             pirate_idx  = sub("^.*_", "", names(s)),
             aligned     = aligned,
             n_in_file   = length(s),
             met, priv)
}))

res[, len_ratio  := ungap_len / median(ungap_len), by = gene_family]
res[, priv_total := priv_p1 + priv_p2 + priv_p3]
res[, priv_outlier := priv_total > pmax(priv_min, priv_mult * median(priv_total, na.rm = TRUE)),
    by = gene_family]

## ---- 2. join fission/fusion flags -------------------------------------
anno <- fread(anno_file)
stopifnot(all(c(col_geno, col_family, col_ftag) %in% names(anno)))
v <- as.character(anno[[col_ftag]])
anno[, ff_colon := !is.na(v) & grepl(":", v, fixed = TRUE)]
anno[, ff_semi  := !is.na(v) & grepl(";", v, fixed = TRUE)]
anno[, ff_flag  := ff_colon | ff_semi]
print(anno[, .N, by = .(ff_colon, ff_semi)])   # sanity check: should be a minority
flags <- anno[, .(ff_flag  = any(ff_flag),
                  ff_colon = any(ff_colon),
                  ff_semi  = any(ff_semi)), by = c(col_geno, col_family)]
setnames(flags, c(col_geno, col_family), c("geno_id", "gene_family"))

res <- merge(res, flags, by = c("geno_id", "gene_family"), all.x = TRUE)
n_unmatched <- res[is.na(ff_flag), .N]
if (n_unmatched > 0)
  warning(n_unmatched, " sequences had no geno_id/gene_family match in anno; check family naming")
res[is.na(ff_flag), `:=`(ff_flag = FALSE, ff_colon = FALSE, ff_semi = FALSE)]

## ---- 3. classify each sequence ----------------------------------------
res[, status := fcase(
  len_mod3 != 0 | n_bad_gap > 0 | lead_gap %% 3 != 0, "frame_broken",
  int_stop > 0 | int_stop_aln > 0,                    "internal_stop",
  len_ratio > len_hi,                                 "length_long",
  len_ratio < len_lo,                                 "length_short",
  default = "clean")]

## ---- 4. per-family summary --------------------------------------------
p3frac <- function(p3, tot) if (sum(tot, na.rm = TRUE) > 0)
  sum(p3, na.rm = TRUE) / sum(tot, na.rm = TRUE) else NA_real_

fam <- res[, .(
  n_seq            = .N,
  aligned          = all(aligned),
  n_ff             = sum(ff_flag),
  n_ff_colon       = sum(ff_colon),
  n_ff_semi        = sum(ff_semi),
  n_problem        = sum(status != "clean"),
  n_problem_ff     = sum(status != "clean" & ff_flag),
  n_problem_noff   = sum(status != "clean" & !ff_flag),
  n_frame_broken   = sum(status == "frame_broken"),
  n_internal_stop  = sum(status == "internal_stop"),
  n_length_long    = sum(status == "length_long"),
  n_length_short   = sum(status == "length_short"),
  n_priv_outlier   = sum(priv_outlier, na.rm = TRUE),
  n_priv_out_ff    = sum(priv_outlier & ff_flag, na.rm = TRUE),
  priv_ff          = sum(priv_total[ff_flag], na.rm = TRUE),
  priv_other       = sum(priv_total[!ff_flag], na.rm = TRUE),
  p3frac_ff        = p3frac(priv_p3[ff_flag], priv_total[ff_flag]),
  p3frac_other     = p3frac(priv_p3[!ff_flag], priv_total[!ff_flag]),
  median_len       = as.numeric(median(ungap_len))
), by = gene_family]

## ---- 5. enrichment test (handover 5.2) --------------------------------
test <- fread(test_file)
setnames(test, c(col_test_family, col_freq, col_seg, col_len, col_pbh),
         c("gene_family", "freq", "seg_sites", "gene_length", "p_bh"))

ff_fam <- flags[, .(n_ff_genomes = sum(ff_flag)), by = gene_family]
d <- merge(test, ff_fam, by = "gene_family", all.x = TRUE)
d <- merge(d, fam[, .(gene_family, n_problem_ff, n_problem, n_priv_out_ff)],
           by = "gene_family", all.x = TRUE)
d[is.na(n_ff_genomes), n_ff_genomes := 0]
d[, screened := !is.na(n_problem)]
d[is.na(n_problem_ff), `:=`(n_problem_ff = 0, n_problem = 0, n_priv_out_ff = 0)]
d[, has_ff      := n_ff_genomes > 0]
d[, has_prob_ff := n_problem_ff > 0]
d[, sig := p_bh < 0.05]
dm <- d[freq >= 2 & gene_length > 0]

m1 <- glm(seg_sites ~ has_ff + log(freq) + log(gene_length),
          family = quasipoisson, data = dm)
m2 <- glm(seg_sites ~ has_prob_ff + log(freq) + log(gene_length),
          family = quasipoisson, data = dm)
m3 <- glm(seg_sites ~ n_ff_genomes + log(freq) + log(gene_length),
          family = quasipoisson, data = dm)

sink(file.path(out_dir, "ff_screen_enrichment.txt"))
cat("Families in test_results:", nrow(d), "\n")
cat("Families with any ff flag:", sum(d$has_ff), "\n")
cat("Flagged families not found in ag_fus_alns:", d[has_ff & !screened, .N], "\n")
cat("Sequences unmatched to anno:", n_unmatched, "\n\n")
cat("Sequence status counts (by ff_flag):\n")
print(res[, .N, by = .(status, ff_flag)][order(status, ff_flag)])
cat("\nPrivate sites by codon position, ff vs other:\n")
print(res[, .(p1 = sum(priv_p1, na.rm = TRUE), p2 = sum(priv_p2, na.rm = TRUE),
              p3 = sum(priv_p3, na.rm = TRUE)), by = ff_flag])
cat("\n== m1: seg_sites ~ has_ff ==\n");       print(summary(m1)$coefficients)
cat("\n== m2: seg_sites ~ has_prob_ff ==\n");  print(summary(m2)$coefficients)
cat("\n== m3: seg_sites ~ n_ff_genomes ==\n"); print(summary(m3)$coefficients)
cat("\n== Significance (p_bh < 0.05) by has_ff ==\n")
tb <- table(has_ff = dm$has_ff, sig = dm$sig); print(tb)
if (all(dim(tb) == 2)) print(fisher.test(tb))
cat("\n== Significance by has_prob_ff ==\n")
tb2 <- table(has_prob_ff = dm$has_prob_ff, sig = dm$sig); print(tb2)
if (all(dim(tb2) == 2)) print(fisher.test(tb2))
cat("\n== Sentinel", sentinel, "==\n")
print(fam[gene_family == sentinel])
print(res[gene_family == sentinel & (status != "clean" | ff_flag | priv_outlier %in% TRUE),
          .(header, ff_flag, status, ungap_len, len_ratio, int_stop_aln,
            priv_p1, priv_p2, priv_p3)])
sink()

## ---- 6. write outputs ---------------------------------------------------
fwrite(fam[order(-n_problem_ff, -n_problem)], file.path(out_dir, "ff_screen_families.csv"))
fwrite(res, file.path(out_dir, "ff_screen_sequences.csv"))
fwrite(res[status != "clean" | priv_outlier %in% TRUE,
           .(gene_family, geno_id, header, ff_flag, status, priv_outlier,
             ungap_len, len_ratio, len_mod3, n_bad_gap, int_stop_aln, priv_total)],
       file.path(out_dir, "ff_mask_list.csv"))

## ---- 7. quick base-R plot ------------------------------------------------
png(file.path(out_dir, "ff_seg_vs_freq.png"), 1600, 800, res = 150)
par(mfrow = c(1, 2))
cols <- ifelse(dm$has_ff, "firebrick", adjustcolor("grey40", 0.4))
plot(dm$freq, dm$seg_sites / dm$gene_length, col = cols, pch = 16, cex = 0.7,
     log = "x", xlab = "freq (genomes)", ylab = "seg_sites / gene_length",
     main = "Flagged (red) vs unflagged")
boxplot(seg_sites / gene_length ~ has_ff, data = dm, outline = FALSE,
        names = c("no ff", "ff"), ylab = "seg_sites / gene_length")
dev.off()

message("Done. Outputs in ", out_dir)
