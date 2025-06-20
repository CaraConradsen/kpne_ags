# list_gffs <- list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_gffs")
# 
# 
# list_f_gffs <- list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_p_filtered_gffs")
# 
# list_f_gffs <- gsub("_filtered", "", list_f_gffs)
# 
# 
# gff_list <- list_gffs[!list_gffs %in% list_f_gffs]
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

gff_dir = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_gffs"
out_gff_dir = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/plsm_filtered_gffs"

list_gffs <- list.files(gff_dir)

plasmid_contigs <- fread(paste0(outdir_dat, "/plasmid_contigs.txt"), header = FALSE)

start.time <- Sys.time()

for(i in list_gffs){
  print(paste0("Processing ....", i))
  #locus name
  loc_nm = paste0("=kps", gsub("SPARK_|.gff|_","", i))
  
  # get gff
  tmp_gff <- readLines(paste(gff_dir, i, sep="/"))
  
  # fix locus names
  tmp_gff <- gsub("=k__",loc_nm, tmp_gff)
  
  if(i %in% unique(plasmid_contigs$V1)){
    
    # assign contigs
    contigs <- plasmid_contigs[V1==i, V2]
    
    # remove header lines
    rm_gff_head_lines = paste(paste0("##sequence-region ",contigs), collapse = "|")
    tmp_gff <- tmp_gff[!grepl(rm_gff_head_lines, tmp_gff)]
    
    # Remove gff lines
    rm_gff_lines = paste(paste0("^",contigs), collapse = "|")
    tmp_gff <- tmp_gff[!grepl(rm_gff_lines, tmp_gff)]
    
    # clean up fasta
    fasta_brk = grep("##FASTA", tmp_gff)
    
    fasta_tmp = tmp_gff[fasta_brk:length(tmp_gff)]
    
    # Find all header indices
    header_lines <- grep("^>", fasta_tmp)
    header_lines <- c(header_lines, length(fasta_tmp) + 1)
    
    # Function to decide if block should be kept
    should_keep_block <- function(i) {
      header_idx <- header_lines[i]
      header_name <- sub("^>", "", fasta_tmp[header_idx])
      !(header_name %in% contigs)
    }
    
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
    result <- fasta_tmp[indices_to_keep]
    
    # add everything together
    tmp_gff <- c(tmp_gff[1:fasta_brk], result)
    
    writeLines(tmp_gff, con = paste(out_gff_dir, i, sep = "/"))
    
  } else {
    writeLines(tmp_gff, con = paste(out_gff_dir, i, sep = "/"))
  }
}

end.time <- Sys.time()
end.time - start.time# Time difference of 19.24803 mins



