
# How does it contrast to core? -------------------------------------------

core_theta_dt <- fread(paste0(outdir_dat, "/core_theta_dt.csv"))

ag_dir_dt <- fread(paste0(outdir_dat, "/set_outliers.csv"))[theta_outlier==1 & w_theta_sel!="balancing"]

with(core_theta_dt,
     boxplot(seg_s_sites,
             ylab = "Synonymous segregating sites",
             xlab = "Core genes", 
             pch = 16,
             col = rgb(0.1,0.1,0.1, alpha = 0.5)))

# What does it physically look like? --------------------------------------


nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))
path_dt <- fread(paste0(outdir_dat, "/path_dt.csv"))

# get path order
nodes_long <- path_dt[, {
  v <- strsplit(nodes, "|", fixed = TRUE)[[1]]
  .(node_id = v, node_order = seq_along(v))
}, by = .(geno_id, path_id)]

nodes_dt <- merge(nodes_dt, nodes_long,
                  all.x = TRUE, 
                  by = c("geno_id","path_id", "node_id"))



# convert start and length to 1 based
nodes_dt[, c("start", "end", "tot_len") := .(start + 1L, end + 1L, tot_len + 1L)]


# localise selected gene --------------------------------------------------
ag_dir_sel <- fread(paste0(outdir_dat, "/set_outliers.csv"))[theta_outlier==1 & p_bh <= 0.05 & w_theta_sel!="balancing"]

ag_dir_sel <- ag_dir_sel$gene_family


ag_dir_sel_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("geno_id", "gene_family",
                                         "start", "end", "strand",  
                                         "anchor", "msu"))[gene_family==ag_dir_sel]

ag_ds_dt <- ag_dir_sel_anno[, .(
  geno_id,
  rstart = start,
  rend   = end,
  strand
)]

setkey(ag_ds_dt, geno_id, rstart, rend, strand)


# blocks that cross the origin
wrap <- nodes_dt[start > end]
lin  <- nodes_dt[start <= end]

nodes_dt_lin <- rbind(
  lin,
  wrap[, .(geno_id, path_id, node_id, block_id, strand,
           start = start, end = tot_len, tot_len, ST, node_order)],
  wrap[, .(geno_id, path_id, node_id, block_id, strand,
           start = 0L,   end = end,     tot_len, ST, node_order)]
)

ag_node_dt <- nodes_dt_lin[, .(
  geno_id,
  block_id,
  nstart = start,
  nend   = end,
  strand
)]

setkey(ag_node_dt, geno_id, nstart, nend, strand)

gene_block <- foverlaps(
  ag_node_dt,
  ag_ds_dt,
  by.x = c("geno_id", "nstart", "nend","strand"),
  by.y = c("geno_id", "rstart", "rend", "strand"),
  type = "any",   
  nomatch = NULL
)
