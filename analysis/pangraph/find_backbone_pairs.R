# find backbone pairs

msu_paths_dt <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"))
path_dt <- fread(paste0(outdir_dat, "/path_dt.csv"))

pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"))
nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))
msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv")) # core block ids

# parms
# geno_ids
geno_ids = unique(pangraph_anno$geno_id)
all_msu_nums = unique(msu_mergers_dt$msu_mergers)

# Get backbone for MSUs --------------------------------------------------

anchor_msu_regions <- function(
    msu_num,
    msu_mergers_dt,
    pangraph_anno,
    nodes_dt
){
  
  msu_x_core_ids = msu_mergers_dt[msu_mergers== msu_num]
  
  # Get graph nodes
  msu_x_core_nodes <- nodes_dt[block_id %chin% msu_x_core_ids$block_id]
  
  msu_x_region <- msu_x_core_nodes[, .(ST = unique(ST), 
                                       start = min(start), 
                                       end = max(end)),
                                   by = geno_id]
  
  # set keys for core nodes
  setnames(msu_x_core_nodes, c("start","end"), c("region_start","region_end"))
  
  setkey(msu_x_core_nodes, geno_id, ST, region_start, region_end)
  
  # set keys for regions
  setnames(msu_x_region, c("start","end"), c("region_start","region_end"))
  
  setkey(pangraph_anno, geno_id, ST, start, end)
  
  setkey(msu_x_region, geno_id, ST, region_start, region_end)
  
  msu_region <- foverlaps(
    x = pangraph_anno,
    y = msu_x_region,
    by.x = c("geno_id", "ST", "start", "end"),
    by.y = c("geno_id", "ST", "region_start", "region_end"),
    type = "any",
    nomatch = 0L
  )
  
  msu_anchors <- foverlaps(
    x = pangraph_anno,
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
    between_genes[, .(locus_tag, anchor)]
  })
  
  
  
  syn_blocks <- rbindlist(syn_blocks)
  
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
  
  
  msu_anchors_anno <- merge(msu_region[, !c("region_start", "region_end")], 
                            msu_anchors, 
                            all.x = TRUE, by = "locus_tag")
  
  msu_anchors_anno$msu = msu_num
  
  # remove missing values
  msu_anchors_anno <- msu_anchors_anno[!is.na(anchor)]
  
  # # fix unachored leading and lagging genes
  # msu_anchors_anno <- msu_anchors_anno[
  #   order(geno_id, start)
  # ][
  #   , {
  #     # Identify anchor rows
  #     anchor_rows <- which(!is.na(anchor) & !grepl(":", anchor))
  #     
  #     if (length(anchor_rows) == 0L) {
  #       # No anchors at all → return unchanged
  #       .SD
  #     } else {
  #       
  #       first_row <- anchor_rows[1]
  #       last_row  <- anchor_rows[length(anchor_rows)]
  #       
  #       first_anchor <- anchor[first_row]
  #       last_anchor  <- anchor[last_row]
  #       
  #       # tidy names
  #       first_anchor <-  ifelse(first_anchor == "A1", ":A1", paste0(first_anchor,":"))
  #       last_anchor <-  ifelse(last_anchor == "A1", ":A1", paste0(last_anchor,":"))
  #       
  #       # Copy anchor col (so := doesn't modify .SD inside the expression)
  #       out <- copy(.SD)
  #       
  #       # Leading unanchored
  #       if (first_row > 1L) {
  #         out[1:(first_row-1), anchor := first_anchor]
  #       }
  #       
  #       # Trailing unanchored
  #       if (last_row < .N) {
  #         out[(last_row+1):.N, anchor := last_anchor]
  #       }
  #       
  #       out
  #     }
  #   },
  #   by = geno_id
  # ]
  
  return(msu_anchors_anno)
  
}

cl <- makePSOCKcluster(num_cores-4)
registerDoParallel(cl)

msu_regions_anchored <- foreach(msu = all_msu_nums[-1],
                                .combine = "rbind",
                                .packages = "data.table") %dopar%{
                                 
                                  # run anchor function
                                    anchor_msu_regions(
                                    msu_num        = msu,
                                    msu_mergers_dt = msu_mergers_dt,
                                    pangraph_anno  = pangraph_anno,
                                    nodes_dt       = nodes_dt
                                  )
                                } 
  

stopCluster(cl)


fwrite(msu_regions_anchored, paste0(outdir_tab, "/msu_regions_anchored.csv"))
