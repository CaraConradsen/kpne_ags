# find backbone pairs
# needs fixing for the origin spanning msu - creates a large junction, overlaps all other msus

msu_paths_dt <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"))
path_dt <- fread(paste0(outdir_dat, "/path_dt.csv"))

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))
nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))
msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv")) # core block ids

# parms
# geno_ids
geno_ids = unique(pan_anno$geno_id)
all_msu_nums = unique(msu_mergers_dt$msu_mergers)

# Get backbone for MSUs --------------------------------------------------

# make ori-spanning region continuous
linearize_ori_intervals<- function(msu_dat, id = "node"){
  
  # set general variables
  min_val = min(msu_dat$start, msu_dat$end, na.rm = TRUE)

  col_order = colnames(msu_dat)
  
  # Circular msus are split
  circular_msus <- msu_dat[start > end]
  
  if(id == "node"){  
    # non-circular msus stay as-is
    linear_msus <- msu_dat[!node_id %chin% circular_msus$node_id]
    
    # Split into two intervals
    split_msus <- rbind(
      # Interval 1: region_start → genome_length
      circular_msus[, .(
        geno_id, path_id, node_id, 
        block_id, strand, 
        start, 
        end = tot_len, 
        tot_len, ST 
      )],
      # Interval 2: 1 → region_end
      circular_msus[, .(
        geno_id, path_id, node_id, 
        block_id, strand, 
        start = min_val, 
        end, 
        tot_len, ST 
      )]
    ) 
    
    split_msus <- split_msus[end != min_val]
    
    }else{
      # non-circular msus stay as-is
      # pan_anno is linearised already
      linear_msus <- msu_dat[!fus_locus_tag %chin% circular_msus$fus_locus_tag]
      
      split_msus = NULL
    }
  
  
  # Combine back with linear MSUs
  msu_dat <- rbindlist(list(linear_msus, split_msus), use.names = TRUE)
  
  # define where the gene lies 
  msu_dat[, mid:= (start + end)/2]
  msu_dat[, post_ori := fcase(mid < (tot_len/2), "post",
                                       default = "pre")]
  
  # adjust length to make region continuous
  msu_dat[post_ori=="post", `:=`(start = start + tot_len,
                                          end = end + tot_len)]
  # clean dat
  msu_dat[, c("mid", "post_ori"):= NULL]
  
  if(id == "node"){ 
  # get grouping columns
  grping_cols = col_order[!col_order %in% c("start", "end")]
  
  
  # collapse split genes into one
  msu_dat <- msu_dat[,
                     .(
                       start   = min(start),
                       end     = max(end)
                     ),
                     by = grping_cols
  ]
  }
  
  setcolorder(msu_dat, col_order)
  
  return(msu_dat)
  
}

anchor_msu_regions <- function(
    msu_num,
    msu_mergers_dt,
    pan_anno,
    nodes_dt,
    ori_span = FALSE){
  
  anno = copy(pan_anno)
  
  msu_x_core_ids = msu_mergers_dt[msu_mergers == msu_num]
  
  # Get graph nodes
  msu_x_core_nodes <- nodes_dt[block_id %chin% msu_x_core_ids$block_id]
  
  # adjust python values to +1
  msu_x_core_nodes[, `:=`(start = start + 1,
                          end = end + 1, 
                          tot_len = tot_len + 1)]
  
  # switch circular value
  ori_span = ifelse(nrow(msu_x_core_nodes[start > end]) > 0, 
                    TRUE, FALSE)
  
  # detect origin spanning genes
  if(isTRUE(ori_span)){
    msu_x_core_nodes <- linearize_ori_intervals(msu_x_core_nodes)
  }
  
  # get msu regions
  msu_x_region <- msu_x_core_nodes[end != tot_len, .(ST = unique(ST), 
                                                     start = min(start), 
                                                     end = max(end)),
                                   by = geno_id]
  
  # set keys for core nodes
  setnames(msu_x_core_nodes, c("start","end"), c("region_start","region_end"))
  
  setkey(msu_x_core_nodes, geno_id, ST, region_start, region_end)
  
  # set keys for regions
  setnames(msu_x_region, c("start","end"), c("region_start","region_end"))
  
  # adjust annotations to match msu_nodes continuity 
  if(isTRUE(ori_span)){
    anno <- linearize_ori_intervals(anno, id = "anno")
  }
  
  setkey(anno, geno_id, ST, start, end)
  
  setkey(msu_x_region, geno_id, ST, region_start, region_end)
  
  msu_region <- foverlaps(
    x = anno,
    y = msu_x_region,
    by.x = c("geno_id", "ST", "start", "end"),
    by.y = c("geno_id", "ST", "region_start", "region_end"),
    type = "any",
    nomatch = 0L
  )
  
  msu_anchors <- foverlaps(
    x = anno,
    y = msu_x_core_nodes,
    by.x = c("geno_id", "ST", "start", "end"),
    by.y = c("geno_id", "ST", "region_start", "region_end"),
    type = "within",
    nomatch = 0L
  )
  

  # get anchors in order
  msu_anchors[, n := .N, gene_family]
  
  anchors = msu_anchors[n==length(geno_ids),.(geno_id, start, gene_family)]
  
  anchor_path <- anchors[
    order(start),                                # sort within each group
    .(gene_families = paste(gene_family, collapse = ", ")), 
    by = geno_id
  ]
  
  # Count occurrences and return the most frequent gene_families string
  mode_path <- anchor_path[
    , .N, by = gene_families
  ][order(-N)][1, gene_families]
  
  
  ordered_families <- data.table(
    order = seq_len(length(strsplit(mode_path, ", ")[[1]])),
    gene_family = strsplit(mode_path, ", ")[[1]]
  )
  
  ordered_families[, anchor:= paste0("A", order)]
  
  syn_blocks <- lapply(seq_len(nrow(ordered_families) - 1), function(i) {
    
    focal_gene_fams <- ordered_families[i:(i+1), gene_family]
    anchor_names    <- paste(ordered_families[i:(i+1), anchor], collapse = ":")
    
    between_genes <- msu_region[
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
    between_genes[, .(fus_locus_tag, anchor)]
  })
  
  
  
  syn_blocks <- rbindlist(syn_blocks)
  
  # add anchors
  
  syn_blocks <- rbind(#leading_block,
    syn_blocks,
    #lagging_block,
    merge(ordered_families[,.(gene_family, anchor)], 
          anno[, .(gene_family, fus_locus_tag)],
          all.x = TRUE, by = "gene_family")[,.(fus_locus_tag, anchor)])
  
  syn_blocks[, n := .N, fus_locus_tag]
  
  msu_anchors = rbind(syn_blocks[n==1][,.(fus_locus_tag, anchor)],
                      syn_blocks[n>1 & !grepl(":", anchor)][,.(fus_locus_tag, anchor)])
  
  
  # msu_anchors_anno <- merge(msu_region[, !c("region_start", "region_end")], 
  #                           msu_anchors, 
  #                           all.x = TRUE, by = "locus_tag")
  
  msu_anchors$msu = msu_num
  
  # remove missing values
  msu_anchors <- msu_anchors[!is.na(anchor) & !is.na(fus_locus_tag)]
  
  return(msu_anchors)
  
}

cl <- makePSOCKcluster(num_cores-4)
registerDoParallel(cl)

msu_regions_anchored <- foreach(msu = all_msu_nums,
                                .combine = "rbind",
                                .packages = "data.table") %dopar%{
                                 
                                  # run anchor function
                                    anchor_msu_regions(
                                    msu_num        = msu,
                                    msu_mergers_dt = msu_mergers_dt,
                                    pan_anno  = pan_anno,
                                    nodes_dt       = nodes_dt
                                  )
                                } 
  

stopCluster(cl)

msu_regions_anchored <- merge(pan_anno,
                              msu_regions_anchored,
                              all.x = TRUE, by = "fus_locus_tag")

# id duplicates within and across msu/junctions ---------------------------

# Examine the msu blocks 

msu_loci <- msu_regions_anchored[grepl(":", anchor)| (anchor=="" | is.na(anchor))]

# Are there any AGs across msu?
msu_loci_grp = unique(msu_loci[,.(ag_type, gene_family, msu)])

across_msu_loci = msu_loci_grp[,.(n = .N), by = c("ag_type", "gene_family")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_msu = 0
msu_regions_anchored[gene_family %chin% across_msu_loci, acrs_msu := 1]

# Examine the syntenic junctions ---------------------------------------------------
# take into account null estimates

junction_loci <- msu_regions_anchored[grepl(":", anchor)]

# Are there any AGs across junctions?
junction_loci_grp = unique(junction_loci[,.(ag_type, gene_family, msu, anchor)])

across_junction_loci = junction_loci_grp[,.(n = .N), by = c("ag_type", "gene_family", "msu")][n>1][,gene_family]


# code across junctions
msu_regions_anchored$acrs_jun = 0
msu_regions_anchored[gene_family %chin% across_junction_loci, acrs_jun := 1]



# Output anchored genes ---------------------------------------------------

fwrite(msu_regions_anchored, paste0(outdir_dat, "/msu_regions_anchored.csv"))
