# library(Biostrings)
# library(data.table)
library(seqinr)

# Load standard codon table (1)
std_genetic_code <- getGeneticCode("1")

# Helper: check if a codon is valid (no gaps, Ns, or ambiguous bases)
is_valid_codon <- function(codon) {
  !any(grepl("-", codon)) && codon %in% names(std_genetic_code)
}

# Helper: is a substitution synonymous?
is_synonymous <- function(codon1, codon2) {
  std_genetic_code[[codon1]] == std_genetic_code[[codon2]]
}

# Function to calculate pairwise piS between two codon sequences
calc_piS_pair <- function(codons1, codons2) {
  syn_sites <- 0
  syn_diffs <- 0
  
  for (i in seq_along(codons1)) {
    c1 <- codons1[i]
    c2 <- codons2[i]
    
    if (!is_valid_codon(c1) || !is_valid_codon(c2)) next
    if (c1 == c2) {
      # Count 1 synonymous site
      syn_sites <- syn_sites + 1
    } else {
      # Count whether substitution is synonymous
      syn_sites <- syn_sites + 1
      if (is_synonymous(c1, c2)) {
        syn_diffs <- syn_diffs + 1
      }
    }
  }
  
  piS <- if (syn_sites > 0) syn_diffs / syn_sites else NA_real_
  return(list(S = syn_sites, Sd = syn_diffs, piS = piS))
}

# Main function for a single FASTA alignment
calc_piS_from_fasta <- function(fasta_file) {
  aa <- readDNAStringSet(fasta_file)
  gene <- tools::file_path_sans_ext(basename(fasta_file))
  
  # Split into codons
  codon_list <- lapply(aa, function(seq) {
    codons <- substring(as.character(seq), seq(1, length(seq) - 2, 3), seq(3, length(seq), 3))
    toupper(codons)
  })
  
  seq_ids <- names(codon_list)
  pairwise <- combn(seq_ids, 2, simplify = FALSE)
  
  # Prepare data.table to store results
  results <- data.table(gene = character(0), seq1 = character(0), seq2 = character(0), S = numeric(0
                                                                                                   ), Sd = numeric(0), piS = numeric(0))
  
  for (pair in pairwise) {
    res <- calc_piS_pair(codon_list[[pair[1]]], codon_list[[pair[2]]])
    
    new_row <- data.table(
      gene = gene,
      seq1 = pair[1],
      seq2 = pair[2],
      S = res$S,
      Sd = res$Sd,
      piS = res$piS
    )
    
    results <- rbind(results, new_row)
  }
  
  return(results)
}

# === Run on all files in directory ===
input_dir <- "./input_data/PIRATE_1695_out/feature_sequences"  # change to your actual folder
fasta_files <- list.files(input_dir, pattern = "nucleotide.fasta", full.names = TRUE)
fasta_files <- as.data.table(fasta_files)

# Subset to AGs
chunk_size <- 100

ag_names <- unique(ag_syn_block_count[n_syn_geno !="1", gene_fam])# remove singletons

chunks <- split(ag_names, ceiling(seq_along(ag_names) / chunk_size))

# Apply each chunk's pattern to subset dt
fasta_files_flt <- lapply(chunks, function(chunk) {
  ag_pattern <- paste(chunk, collapse = "|")
  fasta_files[grepl(ag_pattern, fasta_files, perl = TRUE)]
})

# Combine results into one data.table
fasta_files_flt <- rbindlist(fasta_files_flt, use.names = TRUE, fill = TRUE)

#detect and remove duplicates
fasta_files_flt[, n:=.N, by = fasta_files]

fasta_files_flt <- fasta_files_flt[n <2]

fasta_files_flt$n = NULL




cl <- makeCluster(num_cores-4)
registerDoParallel(cl)

# Use lapply instead of map_dfr
all_piS <- foreach(file = fasta_files[,fasta_files], 
                   .combine = rbind, 
                   .packages = c("Biostrings", "data.table")) %dopar% {
  calc_piS_from_fasta(file)
}

# Stop the parallel backend after the loop is done
stopCluster(cl)

# Save result
fwrite(all_piS, paste0(outdir_dat,"/pairwise_piS_results.csv"))

# View result
print(head(all_piS))