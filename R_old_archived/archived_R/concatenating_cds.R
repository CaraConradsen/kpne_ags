
# cl <- makePSOCKcluster(num_cores-6)
# registerDoParallel(cl)

cds_fasta_gff_dt <- foreach(genome = geno_list,
                            .packages = c("Biostrings", "rtracklayer",
                                          "GenomicRanges", "gUtils")) %do% {
                                            # genome = geno_list[26]
                                            focal_genome = genome
                                            
                                            # import fasta
                                            
                                            fasta_file = lng_rd_fasta_files[grepl(focal_genome, lng_rd_fasta_files)]
                                            
                                            coord_fasta = Biostrings::readDNAStringSet(fasta_file)
                                            
                                            # Get PIRATE cords
                                            coord_fasta = coord_fasta[names(coord_fasta) %in% chr_contigs$seqnames]
                                            
                                            #--- Load PIRATE GFF cord info file ---
                                            
                                            cords = fread(paste0(data_dir, "/co-ords/", focal_genome,".co-ords.tab"))
                                            
                                            colnames(cords)[1] = "locus_tag"
                                            
                                            cords = merge(cords, gene_families[geno_id == focal_genome,
                                                                               .(gene_family, locus_tag, allele_name)],
                                                          all.x = TRUE, by = "locus_tag")
                                            
                                            # Convert strand to "+" / "-"
                                            cords$Strand <- ifelse(cords$Strand == "Forward", "+", "-")
                                            
                                            #convert to lower case
                                            colnames(cords) = tolower(colnames(cords))
                                            
                                            # set contig as seqname
                                            colnames(cords)[which(names(cords)=="contig")] = "seqnames"
                                            
                                            # remove missing loci (short proteins)
                                            cords <- cords[!is.na(gene_family)]
                                            
                                            # remove plasmids
                                            cords <- cords[seqnames %chin% chr_contigs$seqnames]
                                            
                                            # # set order
                                            # setorderv(cords, cols = 'start')
                                            # cords[gene %chin% c("dnaA", "dnaN", "recF", "gyrB", "gyrA")]
                                            # 
                                            # Create grange file
                                            gr_cords <- dt2gr(cords)
                                            
                                            # subset cols
                                            mcols(gr_cords) <- mcols(gr_cords)[, "locus_tag", drop = FALSE]
                                            
                                            # determinenew cds positions
                                            cum_lengths <- cumsum(width(gr_cords))
                                            
                                            adj_gff_dt <- data.frame(start = as.integer(c(1, cum_lengths[1:length(cum_lengths)-1]+1)),
                                                                     end = as.integer(cum_lengths),
                                                                     strand = as.character(strand(gr_cords)),
                                                                     locus_tag = gr_cords$locus_tag)
                                            
                                            adj_gff_dt$geno_id = focal_genome
                                            
                                            # Merge overlapping CDS on same strand
                                            gr_cords <- reduce(gr_cords, ignore.strand = TRUE)
                                            
                                            # Extract sequences
                                            cds_seqs <- getSeq(coord_fasta, gr_cords)
                                            
                                            concatenated_cds <- paste(as.character(cds_seqs), collapse = "")
                                            
                                            concatenated_cds_dna <- DNAString(concatenated_cds)
                                            
                                            # Create a DNAStringSet with one "sequence" per genome
                                            final_fasta <- DNAStringSet(concatenated_cds_dna)
                                            tot_len = width(final_fasta)
                                            names(final_fasta) <- focal_genome  # genome identifier
                                            
                                            # Write CDS FASTA
                                            writeXStringSet(final_fasta, paste0(pangraph_dir,focal_genome, ".fasta"))
                                            
                                            cat("Processed", focal_genome, "\n")
                                            
                                            adj_gff_dt$tot_len = tot_len
                                            
                                            setDT(adj_gff_dt)
                                            
                                          }

cds_fasta_gff_dt <- rbindlist(cds_fasta_gff_dt)
# stopCluster(cl)

# fwrite(cds_fasta_gff_dt, paste0(outdir_dat, "/cds_fasta_gff_dt.csv"))
