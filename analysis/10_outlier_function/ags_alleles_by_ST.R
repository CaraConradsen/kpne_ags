# Peptide alleles of every set 1 AG, tabulated against ST
#
# Each locus is pulled from its genome, translated, and identical peptides
# collapsed into one allele. The result is an allele x ST count matrix per
# gene, written to allele_counts_by_ST.csv so that the extraction is run
# once and every downstream figure reads the cached table
# (see ags_alleles_by_ST_plot.R).
#
# Column sums per gene reproduce its pangenome frequency, so alleles and the
# neutrality test are counted over the same carriers.

allele_file <- paste0(outdir_dat, "/allele_counts_by_ST.csv")

recompute_alleles <- FALSE   # TRUE rebuilds the table from the genomes

if (recompute_alleles || !file.exists(allele_file)) {

  # Using aa ----------------------------------------------------------------
  # accounting for fusions
  anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                       select = c("locus_tag","geno_id", "gene_family", "ST", "threshold",
                                  "freq", "consensus_gene_name","consensus_product",
                                  "seqnames", "start",  "end", "strand",
                                  "alleles_at_maximum_threshold")))

  set_1 <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

  anno <- anno[gene_family %chin% set_1]

  fna_files <-  list.files("./input_data/kpne_260_chr_fasta/",
                           full.names = TRUE, pattern = ".fasta")

  focal_ags <- unique(anno$gene_family)

  # Set up parallel processing
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Paralleled loop to read sequences and translate
  # Per gene, pull every carrier locus, translate, and count peptides by ST
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

  # list must be named (A, B, ...) for idcol to pick up the names
  allele_dt <- rbindlist(ag_st_allele_counts, idcol = "gene_family", fill = TRUE, use.names = TRUE)

  # STs absent from a gene come through as NA -> 0
  st_cols <- setdiff(names(allele_dt), c("gene_family", "allele_name"))
  setnafill(allele_dt, type = "const", fill = 0L, cols = st_cols)

  # optional: order ST columns
  setcolorder(allele_dt, c("gene_family", "allele_name", sort(st_cols)))

  fwrite(allele_dt, allele_file)

} else {

  allele_dt <- fread(allele_file)

}
