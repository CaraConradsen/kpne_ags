# find backbone pairs

msu_paths_dt <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"))
path_dt <- fread(paste0(outdir_dat, "/path_dt.csv"))
pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"))
nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))

msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv")) # core block ids

# parms
# geno_ids
geno_ids = unique(pangraph_anno$geno_id)

# Get backbone for MSU_0 --------------------------------------------------

msu_0_core_ids = msu_mergers_dt[msu_mergers=="MSU_0"]

# Get graph nodes
msu_0_core_nodes <- nodes_dt[block_id %chin% msu_0_core_ids$block_id]

msu_0_region <- msu_0_core_nodes[, .(ST = unique(ST), 
                                    start = min(start), 
                                    end = max(end)),
                                by = geno_id]

setnames(msu_0_region, c("start","end"), c("region_start","region_end"))

setkey(pangraph_anno, geno_id, ST, start, end)

setkey(msu_0_region, geno_id, ST, region_start, region_end)

res <- foverlaps(
  x = pangraph_anno,
  y = msu_0_region,
  by.x = c("geno_id", "ST", "start", "end"),
  by.y = c("geno_id", "ST", "region_start", "region_end"),
  type = "within",
  nomatch = 0L
)

res[, n := .N, gene_family]

# get anchors in order
anchors = res[n==93,.(geno_id, start, gene_family)]

result <- anchors[
  order(start),                                # sort within each group
  .(gene_families = paste(gene_family, collapse = ", ")), 
  by = geno_id
]

# Count occurrences and return the most frequent gene_families string
mode_string <- result[
  , .N, by = gene_families
][order(-N)][1, gene_families]


ordered_families <- data.table(
  order = seq_len(length(strsplit(mode_string, ", ")[[1]])),
  gene_family = strsplit(mode_string, ", ")[[1]]
)

ordered_families[, anchor:= paste0("A", order)]

syn_blocks <- lapply(seq_len(nrow(ordered_families) - 1), function(i) {
  
  focal_gene_fams <- ordered_families[i:(i+1), gene_family]
  anchor_names    <- paste(ordered_families[i:(i+1), anchor], collapse = ":")
  
  between_genes <- res[
    , {
      idx <- which(gene_family %in% focal_gene_fams)
      
      # Require exactly two anchors per geno_id
      if (length(idx) == 2L) {
        
        first_pos  <- idx[1]
        second_pos <- idx[2]
        forward <- first_pos < second_pos
        
        if (forward) {
          i1 <- first_pos
          i2 <- second_pos
          block <- .SD[(i1+1):(i2-1)]
        } else {
          i1 <- second_pos
          i2 <- first_pos
          block <- .SD[(i1+1):(i2-1)]
          if (nrow(block) > 0L) {
            block <- block[.N:1]    # reverse orientation
          }
        }
        
        block
        
      } else {
        NULL
      }
    },
    by = geno_id
  ]
  
  between_genes[, anchor := anchor_names]
  between_genes[, .(locus_tag, anchor)]
})



syn_blocks <- rbindlist(syn_blocks)

# # fix trails
# first_fam <- ordered_families[1, gene_family]
# last_fam  <- ordered_families[.N, gene_family]
# 
# first_anchor_name <- paste0(":", ordered_families[1, anchor])
# last_anchor_name  <- paste0(ordered_families[.N, anchor], ":")
# 
# # leading block
# leading_block <- res[
#   , {
#     idx <- which(gene_family == first_fam)
#     if (length(idx) == 1L && idx > 1L) {
#       .SD[1:(idx - 1)]
#     } else NULL
#   },
#   by = geno_id
# ][, .(locus_tag)]
# 
# if (nrow(leading_block) > 0){
#   leading_block[, anchor := first_anchor_name]
# }
# 
# # lagging block
# lagging_block <- res[
#   , {
#     idx <- which(gene_family == last_fam)
#     if (length(idx) == 1L && idx < .N) {
#       .SD[(idx + 1):.N]
#     } else NULL
#   },
#   by = geno_id
# ][, .(locus_tag)]
# 
# if (nrow(lagging_block) > 0){
#   lagging_block[, anchor := last_anchor_name]
# }
  
# add anchors

syn_blocks <- rbind(#leading_block,
                    syn_blocks,
                    #lagging_block,
                    merge(ordered_families[,.(gene_family, anchor)], 
                          pangraph_anno[, .(gene_family, locus_tag)],
                          all.x = TRUE, by = "gene_family")[,.(locus_tag, anchor)])

syn_blocks[, n := .N, locus_tag]

msu_anchors = rbind(syn_blocks[n==1][,.(locus_tag, anchor)],
                    syn_blocks[n>1 & !grepl(":", anchor)][,.(locus_tag, anchor)])


msu_anchors_anno <- merge(res[, !c("region_start", "region_end", "n")], 
                          msu_anchors, 
                          all.x = TRUE, by = "locus_tag")

# SPARK_121_C2

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

