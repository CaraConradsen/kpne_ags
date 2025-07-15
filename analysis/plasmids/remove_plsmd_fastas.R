# list_fasta <- list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_gffs")
# 
# 
# list_f_gffs <- list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_p_filtered_gffs")
# 
# list_f_gffs <- gsub("_filtered", "", list_f_gffs)
# 
# 
# gff_list <- list_fasta[!list_fasta %in% list_f_gffs]
# 
# 
# # Define input and output directories (Windows paths, escaped)
# input_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_gffs"
# output_dir <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_p_filtered_gffs"
# 
# # Create full file paths
# source_paths <- file.path(input_dir, gff_list)
# dest_paths <- file.path(output_dir, gff_list)
# 
# # Copy files
# file.copy(from = source_paths, to = dest_paths, overwrite = TRUE)

#dir.create("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/plsm_filtered_gffs", recursive = TRUE)

fasta_dir = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_fasta"
out_fasta_dir = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/plsm_filtered_fasta"

list_fasta <- list.files(fasta_dir)

plasmid_contigs <- fread(paste0(outdir_dat, "/plasmid_contigs.txt"), header = FALSE)
plasmid_contigs$V1 <- gsub(".gff", ".fasta", plasmid_contigs$V1)


# Function to decide if block should be kept
should_keep_block <- function(i) {
  header_idx <- header_lines[i]
  header_name <- sub("^>", "", tmp_fasta[header_idx])
  !(header_name %in% contigs)
}

start.time <- Sys.time()

for(i in list_fasta){
  print(paste0("Processing ....", i))

  # get gff
  tmp_fasta <- readLines(paste(fasta_dir, i, sep="/"))
  
  if(i %in% unique(plasmid_contigs$V1)){
    
    # assign contigs
    contigs <- plasmid_contigs[V1==i, V2]
    
    # clean up fasta
    # Find all header indices
    header_lines <- grep("^>", tmp_fasta)
    header_lines <- c(header_lines, length(tmp_fasta) + 1)
    
    # Apply over all blocks except the sentinel
    keep_blocks <- sapply(seq_along(header_lines[-length(header_lines)]), should_keep_block)
    
    # For each block, generate the sequence of indices for that block
    blocks_indices <- mapply(function(start, end) seq(start, end),
                             header_lines[-length(header_lines)],
                             header_lines[-1] - 1,
                             SIMPLIFY = FALSE)
    
    # Combine indices for blocks we want to keep
    indices_to_keep <- unlist(blocks_indices[keep_blocks])
    
    # Subset the original vector
    result <- tmp_fasta[indices_to_keep]
    
    writeLines(result, con = paste(out_fasta_dir, i, sep = "/"))
    
  } else {
    writeLines(tmp_fasta, con = paste(out_fasta_dir, i, sep = "/"))
  }
}

end.time <- Sys.time()
end.time - start.time# Time difference of 6.649714 mins



