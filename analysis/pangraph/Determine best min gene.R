# Summarise pangraph syntenic blocks

# Import json -------------------------------------------------------------
# Read from a local file
jsons <- list.files("./input_data/test_pangraph/jsons",
                    full.names = TRUE)

get_pangeno_info <- function(json_file){
  pangraph_data <- fromJSON(json_file)
  
  min_size = gsub("graph_|\\.json","", basename(json_file))
  
  # Extract nodes and align to gene_families --------------------------------
  nodes_dt <- rbindlist(lapply(pangraph_data$nodes, function(x) {
    data.table(
      node_id = as.character(format(x$id, scientific = FALSE)),
      block_id = as.character(format(x$block_id, scientific = FALSE)),
      path_id = x$path_id,
      strand = x$strand,
      start = x$position[1],
      end = x$position[2]
    )
  }), fill = TRUE)
  nodes_dt$min_size = min_size
  
  # extract paths
  # Extract relevant fields from each path
  path_dt <- rbindlist(lapply(pangraph_data$paths, function(x) {
    data.table(
      path_id = x$id,
      geno_id = x$name,
      tot_len = x$tot_len,
      # convert node IDs to character to preserve precision
      nodes = list(as.character(format(x$nodes, scientific = FALSE)))
    )
  }), fill = TRUE)
  path_dt$min_size = min_size
  
  # extract paths
  # Extract relevant fields from each path
  blocks_dt <- rbindlist(lapply(pangraph_data$blocks, function(x) {
    data.table(
      block_id = as.character(format(x$id, scientific = FALSE)),
      block_con_len = nchar(x$consensus)
    )
  }), fill = TRUE)
  blocks_dt$min_size = min_size
  
  return(list(nodes_dt, path_dt, blocks_dt))

}

pan_jsons <- lapply(jsons, function(json) get_pangeno_info(json))

node_info_dts <- rbindlist(lapply(pan_jsons, `[[`, 1))
path_info_dts <- rbindlist(lapply(pan_jsons, `[[`, 2))
block_info_dts <- rbindlist(lapply(pan_jsons, `[[`, 3))


#plot blocks
block_info_dts[, n :=.N, by= min_size]
block_info_dts$size = as.integer(gsub("s", "", block_info_dts$min_size))

library(vioplot)

plot(NULL, xlab = "set min syntenic block size (bp)",
     ylab = "log10 block sizes (bp)", bty = "L", 
     xlim = c(0, 1600),
     yaxt = "n",
     ylim = round(range(log10(block_info_dts$block_con_len))))
text(-100, 5.25, 
     "n syntenic\nblocks", xpd=TRUE)
axis(side = 2, at = seq(2, 5, 0.5),
     sprintf("%1.1f", seq(2, 5, 0.5)), las = 2)

# loop over each group and add violins
for (x in unique(block_info_dts$min_size)) {
  temp_dt <- block_info_dts[min_size == x]
  
  vioplot(log10(temp_dt$block_con_len),
          at = unique(temp_dt$size),
          col = "dodgerblue",
          border = "dodgerblue",
          add = TRUE,
          wex = 90)  # roughly scales violin width
  
  text(unique(temp_dt$size),
       5.15, srt=45, pos = 3,
       temp_dt$n[1],
       xpd = TRUE)
}

# aggregate(subtree_br_len ~ freq, data = res, FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x)))

