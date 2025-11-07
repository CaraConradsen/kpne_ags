# find backbone pairs

msu_paths_dt <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"))
path_dt <- fread(paste0(outdir_dat, "/path_dt.csv"))
pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"))
nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))
msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv"))


# MSU maintain core gene rank order ---------------------------------------

cl <- makePSOCKcluster(num_cores-6)
registerDoParallel(cl)

anno_msu_dt<- foreach(i = unique(pangraph_anno$geno_id),
                      .packages = c("data.table", "GenomicRanges",
                                    "gUtils", "rlang"),
                      .export = c("shift"),
                      .combine = "rbind") %dopar% {
                        # annoation into grange
                        focal_anno = pangraph_anno[geno_id == i]
                        
                        focal_anno[, seqnames:= geno_id]
                        
                        # get genome length
                        genome_length = focal_anno[,.(max(end))]
                        
                        focal_anno_gr <- dt2gr(focal_anno)
                        
                        # Get MSU into grange
                        msu_pos_dt <- merge(nodes_dt[geno_id == i], msu_mergers_dt,
                                            all.x = TRUE, by = "block_id")
                        
                        msu_pos_dt[, seqnames:= geno_id]
                        
                        msu_pos_dt <- msu_pos_dt[!is.na(msu_mergers),.(seqnames, start, end, strand, msu_mergers)]
                        
                        wrapped_msu = msu_pos_dt[start > end, msu_mergers]
                        
                        # linearlize
                        msu_pos_dt <- rbind(
                          # keep normal ranges
                          msu_pos_dt[start <= end],
                          
                          # split circular ranges into two
                          msu_pos_dt[start > end][, .(
                            seqnames = c(seqnames, seqnames),
                            start = c(start, 0L),
                            end = as.integer(c(genome_length, end)),
                            strand = c(strand, strand),
                            msu_mergers = c(msu_mergers, msu_mergers)
                          )]
                        )[order(seqnames, start)]
                        
                        if(is_empty(wrapped_msu)){
                          msu_pos_dt <- msu_pos_dt[, .(
                            seqnames = unique(seqnames),
                            start = min(start),
                            end = max(end)
                          ), by = msu_mergers]
                        }else{
                          wrapped_msu_dt <- msu_pos_dt[msu_mergers %chin% wrapped_msu]
                          setorder(wrapped_msu_dt, start)
                          wrapped_msu_dt[, start_next := data.table::shift(start, -1)]
                          wrapped_msu_dt[, start_next := log10(start_next - start)]
                          threshold = wrapped_msu_dt[start_next == max(start_next, na.rm = TRUE), end]
                          
                          wrapped_msu_dt <- rbind(wrapped_msu_dt[end <= threshold, .(
                            seqnames = unique(seqnames),
                            start = min(start),
                            end = max(end)
                          ), by = msu_mergers],
                          wrapped_msu_dt[end > threshold, .(
                            seqnames = unique(seqnames),
                            start = min(start),
                            end = max(end)
                          ), by = msu_mergers])
                          
                          msu_pos_dt <- rbind(
                            msu_pos_dt[!msu_mergers %chin% wrapped_msu, .(
                              seqnames = unique(seqnames),
                              start = min(start),
                              end = max(end)
                            ), by = msu_mergers], 
                            wrapped_msu_dt)
                          
                        }
                        
                        msu_gr <- suppressWarnings(dt2gr(msu_pos_dt))
                        
                        # fix lengths
                        seqlengths(msu_gr) <- msu_pos_dt[,.(max(end))][[1]]
                        seqlengths(focal_anno_gr) <- msu_pos_dt[,.(max(end))][[1]]
                        
                        # find hits
                        hits <- findOverlaps(msu_gr, focal_anno_gr)
                        
                        # # Get overlap widths
                        # overlap_width <- width(pintersect(msu_gr[queryHits(hits)],
                        #                                   focal_anno_gr[subjectHits(hits)]))
                        # 
                        # # Put in dt
                        # overlap_width_dt <- data.table(data.frame(
                        #   msu_mergers = msu_gr$msu_mergers[queryHits(hits)],
                        #   locus_tag = focal_anno_gr$locus_tag[subjectHits(hits)],
                        #   overlap = overlap_width
                        # )
                        # )
                        # 
                        # # For each block_id, keep the row with the maximum overlap
                        # best_hits <- overlap_width_dt[, .SD[which.max(overlap)], by = locus_tag]
                        # 
                        # focal_anno_blk <- merge(focal_anno, best_hits[,.(msu_mergers,locus_tag)],
                        #                         all.x = TRUE, by = "locus_tag")
                        
                        # Extract hits
                        q <- msu_gr[queryHits(hits)]
                        s <- focal_anno_gr[subjectHits(hits)]
                        
                        # Option 1: create a GRanges representing the overlapping region
                        ovl_gr <- pintersect(q, s)
                        
                        # Add metadata from both sides
                        mcols(ovl_gr) <- cbind(mcols(q), mcols(s))
                        
                        gr2dt(ovl_gr)
                        
                      }

stopCluster(cl)

# rank core gene order in a msu
rank_core <- anno_msu_dt[ag_type == "core",
                         .(geno_id, start, gene_family, msu_mergers, strand)]

# Sort properly within groups
rank_core_fwd <- rank_core[strand == "+"]
rank_core_rev <- rank_core[strand == "-"]

setorderv(rank_core_fwd, 
          cols = c("geno_id", "msu_mergers", "start"),
          c(1,1,1))
setorderv(rank_core_rev, 
          cols = c("geno_id", "msu_mergers", "start"),
          c(1,1,-1))

# Create next_gene using shift
rank_core <- rbind(rank_core_fwd, rank_core_rev)

rank_core[, next_gene := shift(gene_family, type = "lead"), 
          by = .(geno_id, msu_mergers)]

rank_core[, core_pair_str := paste(gene_family,next_gene, collapse = ":")]

rank_core <- rank_core[!is.na(next_gene), .(msu_mergers, geno_id, core_pair_str)]

core_pair_counts <- rank_core[, .N, by = .(msu_mergers, core_pair_str)]

backbone_pairs <- core_pair_counts[N==93,]

backbone_pairs[, c("x1", "x2"):= tstrsplit(core_pair_str, ":", fill = TRUE, keep = 1:2)]

backbone_pairs <- melt(backbone_pairs,
                       id.vars= c("msu_mergers","core_pair_str"), 
                       measure.vars = c("x1","x2"),
                       value.name = "gene_family")

backbone_pairs[,variable:= NULL]

backbone_pairs$in_pair = "anchor"

pairs_anno <- pangraph_anno

core_pairs = unique(backbone_pairs$core_pair_str)


cl <- makePSOCKcluster(num_cores - 6)
registerDoParallel(cl)

locus_backbone <- foreach(i = core_pairs,
                          .packages = c("data.table"),
                          .combine = "rbind") %dopar% {
                            # make a local copy of pairs_anno to avoid modifying the global one
                            pa <- copy(pairs_anno)
                            
                            # get the two anchor genes for this pair
                            anchor_genes <- backbone_pairs[core_pair_str == i, gene_family]
                            
                            # get their coordinates
                            anchor_genes_loc <- pa[gene_family %chin% anchor_genes,
                                                   .(geno_id, locus_tag, start, end)]
                            anchor_genes_loc[, in_pair := "anchor"]
                            
                            # find the start and end coordinates of the region
                            anchor_genes_pos <- anchor_genes_loc[, .(
                              interval_start = min(start),
                              interval_end = max(end)
                            ), by = geno_id]
                            
                            setkey(pa, geno_id, start, end)
                            setkey(anchor_genes_pos, geno_id, interval_start, interval_end)
                            
                            # find genes that overlap the interval
                            subset_dt <- foverlaps(
                              pa,
                              anchor_genes_pos,
                              by.x = c("geno_id", "start", "end"),
                              by.y = c("geno_id", "interval_start", "interval_end"),
                              nomatch = 0L
                            )
                            
                            # keep only genes fully within interval
                            subset_dt <- subset_dt[start >= interval_start & end <= interval_end]
                            
                            # remove duplicates (some may overlap boundary edges)
                            subset_dt <- unique(subset_dt[, .(locus_tag)])
                            
                            # label these genes as being within this pair
                            subset_dt[, in_pair := paste0("in_", i)]
                            
                            # combine interior genes with anchors
                            res <- rbind(subset_dt,
                                         anchor_genes_loc[, .(locus_tag, in_pair)],
                                         use.names = TRUE, fill = TRUE)
                            
                            unique(res)
                          }

stopCluster(cl)



pangraph_anno

# count AGs in and out of msus# core_pair_strcount AGs in and out of msus
count_msu_ags <-rbind(anno_msu_dt[number_genomes!= 93 & !is.na(msu_mergers), 
                                  .(n=.N, in_msu = "y"), geno_id],
                      anno_msu_dt[number_genomes!= 93 & is.na(msu_mergers), 
                                  .(n=.N, in_msu = "n"), geno_id])

count_msu_ags <- dcast(count_msu_ags,
                       geno_id~in_msu, 
                       value.var = "n")

count_msu_ags[, prop_in := y / (n + y)]

count_msu_ags[, .(mean_prop = mean(prop_in), se_prop = sd(prop_in)/sqrt(.N))]

