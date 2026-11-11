# Screen fission/fusion AG loci for nonsense mutations, frameshifts and 
# in-frame deletions
# Each fission/fusion locus (fus_locus_tag "(a:b)") is extracted over its full
# span (fstart:fend), aligned to the PIRATE representative of its family, and
# the premature stop (if any) is attributed to a substitution or an indel.

# Parameters --------------------------------------------------------------

hp_min      <- 5    # homopolymer run length that flags an indel as a possible
                    # ONT assembly error
term_codons <- 1    # stops within this many codons of the reference end are 
                    # not premature

# Import data -------------------------------------------------------------

set_1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full_intact_fis_fus.csv"),
                     select = c("fus_locus_tag", "geno_id", "locus_tag", "gene_family",
                                "gene", "seqnames", "start", "end", "strand",
                                "fstart", "fend", "fstrand")))

# fission/fusion loci of set 1 AGs
fis_dt <- anno[gene_family %chin% set_1_ags & grepl("^\\(", fus_locus_tag)]

fis_loci <- fis_dt[, .(n_orfs   = .N,
                       genes    = paste(unique(gene[!is.na(gene) & gene != ""]), collapse = ","),
                       seqnames = unique(seqnames),
                       start    = unique(fstart),
                       end      = unique(fend),
                       strand   = unique(fstrand)),
                   by = .(fus_locus_tag, geno_id, gene_family)]

# loci spanning the origin cannot be extracted as one interval
fis_loci[, ori_span := start > end]

# Representative sequences --------------------------------------------------

rep_seqs <- readDNAStringSet("./input_data/PIRATE_260_hybrid_chr_hps0.6_out/representative_sequences.ffn")

rep_info <- data.table(header = names(rep_seqs))
rep_info[, gene_family := sub(";.*", "", header)]
rep_info[, rep_locus := sub(".*locus_tag=([^;]+).*", "\\1", header)]

names(rep_seqs) <- rep_info$gene_family

# flag families whose representative is itself a fission/fusion locus
rep_info[, rep_is_fis_fus := rep_locus %chin% fis_dt$locus_tag]

fis_loci[rep_info, on = .(gene_family), rep_is_fis_fus := i.rep_is_fis_fus]

# Functions ---------------------------------------------------------------

# pairwiseAlignment moved to pwalign in Bioconductor 3.19
if (requireNamespace("pwalign", quietly = TRUE)) {
  aln_fun <- pwalign::pairwiseAlignment
  aln_pat <- pwalign::alignedPattern
  aln_sub <- pwalign::alignedSubject
  sub_mat <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
} else {
  aln_fun <- Biostrings::pairwiseAlignment
  aln_pat <- Biostrings::alignedPattern
  aln_sub <- Biostrings::alignedSubject
  sub_mat <- Biostrings::nucleotideSubstitutionMatrix(match = 1, mismatch = -3)
}

# longest run of one base covering ref position p (and p + 1 for insertions)
hp_run <- function(ref_chr, p) {
  p <- min(max(p, 1), length(ref_chr))
  b <- ref_chr[p]
  lo <- p; hi <- p
  while (lo > 1 && ref_chr[lo - 1] == b) lo <- lo - 1
  while (hi < length(ref_chr) && ref_chr[hi + 1] == b) hi <- hi + 1
  hi - lo + 1
}

# call indels and premature stops for one locus against its reference
call_lof <- function(loc_seq, ref_seq) {
  
  aln <- aln_fun(pattern = loc_seq, subject = ref_seq, type = "global",
                 substitutionMatrix = sub_mat, 
                 gapOpening = 5, gapExtension = 2)
  
  pa <- strsplit(as.character(aln_pat(aln)), "")[[1]]
  sa <- strsplit(as.character(aln_sub(aln)), "")[[1]]
  
  ref_chr <- strsplit(as.character(ref_seq), "")[[1]]
  
  # alignment column -> ungapped position in locus and reference
  loc_pos <- cumsum(pa != "-")
  ref_pos <- cumsum(sa != "-")
  
  n_col <- length(pa)
  
  # internal gap runs only (terminal gaps = start/end differences)
  gap_runs <- function(gapped, type) {
    r <- rle(gapped)
    ends   <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    keep   <- r$values & starts > 1 & ends < n_col
    if (!any(keep)) return(NULL)
    data.table(type = type, col_start = starts[keep], len = r$lengths[keep])
  }
  
  indels <- rbind(gap_runs(pa == "-", "deletion"),   # bases missing in locus
                  gap_runs(sa == "-", "insertion"))  # bases added in locus
  
  if (!is.null(indels)) {
    indels[, ref_at := ref_pos[col_start]]
    indels[, loc_at := loc_pos[col_start]]
    indels[, in_frame := len %% 3 == 0]
    indels[, hp_len := vapply(ref_at, function(p) hp_run(ref_chr, p), numeric(1))]
    indels[, homopolymer := hp_len >= hp_min]
  }
  
  # first stop in frame 1 of the locus
  loc_len <- length(loc_seq) - length(loc_seq) %% 3
  aa <- as.character(suppressWarnings(
    translate(subseq(loc_seq, 1, loc_len), genetic.code = getGeneticCode("11"),
              if.fuzzy.codon = "X")))
  
  stop_codon <- regexpr("\\*", aa)[1]
  ref_codons <- floor(length(ref_seq) / 3)
  
  premature <- FALSE; stop_ref_frac <- NA_real_; stop_cause <- NA_character_
  
  if (stop_codon > 0) {
    
    stop_nt  <- (stop_codon - 1) * 3 + 1
    stop_col <- which(loc_pos == stop_nt & pa != "-")[1]
    stop_ref <- ref_pos[stop_col]
    
    premature <- stop_ref <= (ref_codons - term_codons) * 3
    stop_ref_frac <- stop_ref / length(ref_seq)
    
    if (premature) {
      
      fs_up <- if (is.null(indels)) FALSE else 
        indels[!in_frame & loc_at < stop_nt, .N > 0]
      
      ref_codon <- paste(sa[stop_col:min(stop_col + 2, n_col)], collapse = "")
      
      stop_cause <- fcase(
        fs_up,                                       "frameshift",
        !ref_codon %chin% c("TAA", "TAG", "TGA"),    "nonsense",
        default =                                    "other")
    }
  }
  
  # summarise
  fs_dt  <- if (is.null(indels)) NULL else indels[in_frame == FALSE]
  del_dt <- if (is.null(indels)) NULL else indels[in_frame & type == "deletion"]
  
  # locus bases on reference coordinates ("-" = deleted in locus; bases
  # inserted in the locus are dropped), used to add the copy to the AG alignment
  loc_on_ref <- paste(pa[sa != "-"], collapse = "")
  
  data.table(
    loc_len           = length(loc_seq),
    ref_len           = length(ref_seq),
    pid               = pid(aln),
    premature_stop    = premature,
    stop_ref_frac     = stop_ref_frac,
    stop_cause        = stop_cause,
    n_frameshift      = if (is.null(fs_dt)) 0L else nrow(fs_dt),
    n_frameshift_hp   = if (is.null(fs_dt)) 0L else fs_dt[, sum(homopolymer)],
    n_inframe_del     = if (is.null(del_dt)) 0L else nrow(del_dt),
    max_inframe_del_codons = if (is.null(del_dt) || nrow(del_dt) == 0) 0L else 
      as.integer(max(del_dt$len) / 3),
    indels            = if (is.null(indels)) "" else 
      indels[, paste0(type, ":", len, "@", ref_at, fifelse(homopolymer, "(hp)", ""), 
                      collapse = ";")],
    loc_on_ref        = loc_on_ref
  )
}

# Run per genome ----------------------------------------------------------

fna_files <- list.files("./input_data/kpne_260_chr_fasta/", 
                        full.names = TRUE, pattern = ".fasta")

focal_genos <- unique(fis_loci[ori_span == FALSE, geno_id])

cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

lof_calls <- foreach(g = focal_genos,
                     .combine = "rbind",
                     .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                       
                       g_dat  <- fis_loci[geno_id == g & ori_span == FALSE]
                       
                       genome <- readDNAStringSet(fna_files[grepl(paste0(g, ".fasta"), fna_files)])
                       names(genome) <- sub(" .*", "", names(genome))
                       
                       loc_seqs <- getSeq(genome, dt2gr(g_dat[, .(seqnames, start, end, strand)]))
                       
                       out <- lapply(seq_len(nrow(g_dat)), function(i) {
                         cbind(g_dat[i, .(fus_locus_tag, geno_id, gene_family)],
                               call_lof(loc_seqs[[i]], rep_seqs[[g_dat$gene_family[i]]]))
                       })
                       
                       rbindlist(out)
                     }

stopCluster(cl)

# Classify ----------------------------------------------------------------

lof_calls <- merge(fis_loci[, .(fus_locus_tag, geno_id, gene_family, n_orfs, genes, 
                                rep_is_fis_fus, ori_span)],
                   lof_calls, by = c("fus_locus_tag", "geno_id", "gene_family"), 
                   all.x = TRUE)

lof_calls[, lof_class := fcase(
  ori_span,                                               "not_assessed_origin",
  premature_stop & stop_cause == "nonsense",              "nonsense",
  premature_stop & stop_cause == "frameshift" & 
    n_frameshift_hp == n_frameshift,                      "frameshift_homopolymer",
  premature_stop & stop_cause == "frameshift",            "frameshift",
  premature_stop,                                         "premature_stop_other",
  n_inframe_del > 0L,                                     "inframe_deletion",
  default =                                               "no_lof_detected")]

# loci joining differently named genes may be true fusions rather than
# a gene split by a stop
lof_calls[, multi_gene := grepl(",", genes)]

lof_calls[, .N, by = .(lof_class, rep_is_fis_fus)][order(-N)]

# Output ------------------------------------------------------------------

fwrite(lof_calls, paste0(outdir_dat, "/fis_fus_lof_calls.csv"))
