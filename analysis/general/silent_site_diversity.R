calculate_piS <- function(seq1, seq2) {
  
  seq1 = as.character(seq1)
  seq2 = as.character(seq2)
  
  # Ensure input sequences are same length and divisible by 3
  if (nchar(seq1) != nchar(seq2)) stop("Sequences must be same length")
  if (nchar(seq1) %% 3 != 0) stop("Sequences must be codon-aligned")
  
  # Split into codons
  codons1 <- substring(seq1, seq(1, nchar(seq1), 3), seq(3, nchar(seq1), 3))
  codons2 <- substring(seq2, seq(1, nchar(seq2), 3), seq(3, nchar(seq2), 3))
  
  # Initialise counts
  synonymous_differences <- 0
  synonymous_sites <- 0
  
  # Codon table for translation
  std_table <- Biostrings::getGeneticCode("Standard", full.search = TRUE)
  
  for (i in seq_along(codons1)) {
    codon1 <- toupper(codons1[i])
    codon2 <- toupper(codons2[i])
    
    # Skip codons with ambiguous bases
    if (grepl("[^ACGT]", codon1) || grepl("[^ACGT]", codon2)) next
    
    aa1 <- std_table[[codon1]]
    aa2 <- std_table[[codon2]]
    
    # Estimate synonymous sites using Li (1985) approximate method
    syn_sites <- sum(sapply(1:3, function(pos) {
      bases <- c("A", "C", "G", "T")
      alt_codons <- sapply(bases, function(b) {
        tmp <- strsplit(codon1, "")[[1]]
        tmp[pos] <- b
        paste0(tmp, collapse = "")
      })
      sum(sapply(alt_codons, function(ac) {
        ac != codon1 && std_table[[ac]] == aa1
      }))
    })) / 3  # average per codon
    
    synonymous_sites <- synonymous_sites + syn_sites
    
    # Count synonymous difference if codons differ but amino acids are the same
    if (codon1 != codon2 && aa1 == aa2) {
      nuc_diff <- sum(strsplit(codon1, "")[[1]] != strsplit(codon2, "")[[1]])
      synonymous_differences <- synonymous_differences + nuc_diff
    }
  }
  
  piS <- synonymous_differences / synonymous_sites
  return(piS)
}
