# Codon table for translation
std_table <- Biostrings::getGeneticCode("Standard", full.search = TRUE)

s2triplet <- function(s) {
  n <- nchar(s)
  starts <- seq(1, n - 2, by = 3)
  substring(s, starts, starts + 2)
}

seq1 = as.character(focal_ag_aln[1])
seq2 = as.character(focal_ag_aln[3])

# seq1 = "ATGACTCAGATCCAAAAACAGCGCGTAGTCCGTTTCGATGGCAATAAGCAGATCGTTGAAGTTCCCGATCCGGCGCCGGCGGTAATTGGCGCACCCACGACCACTGACTACGGCGGCGTAAAGCTCGGAGCCTCCATTGCCGTCGCTGCGGCAGCTACCGCAACCGCAGATACCGCATCCAGTGCAAGCGATGTGGCTGGTCTCCTTGTTGATCACAACGACCTGGTGACGAAGTACAACGCGCTTCTCAATGATGCCGCGGCTCTTCGTACCACGCTGAATGCTGTCCTTACGCAGCTGAAAGCCAAAACAATTCCTGTTTAA" 
# seq2 = "ATGACTCAGATCCAAAAACAGCGCGTAGTCCGTTTCGATGGCAATAAGCAGATCGTTGAAGTTCCCGATCCGGCGCCGGCGGTAATTGGCGCGCCCACGGCCACTGAATACGGCGGCGTAAAGCTCGGAGCCTCCATTGCCGCCGCTGCGGCAGCTACCGCAACCGCAGATACCGCATCCAGTGCAAGCGATGTGGCTGGTCTCCTTGTTGATCACAACGACCTGGTGACGAAGTACAACGCGCTTCTCAATGATGCCGCGGCTCTTCGTACCACGCTGAATGCTGTCCTTACGCAGCTGAAAGCCAAAACAATTCCTGTTTAA" 

codons <- lapply(c(seq1, seq2), s2triplet)

diff_pos <- which(mapply(`!=`, codons[[1]], codons[[2]]))

codons_diff <- lapply(codons, `[`, diff_pos)

# Logical vector of valid codon pairs
is_good <- mapply(function(x, y) {
  !grepl("[^ATGC]", x) && !grepl("[^ATGC]", y)
}, codons_diff[[1]], codons_diff[[2]])

# Subset both lists by valid positions
codons_diff_clean <- lapply(codons_diff, `[`, is_good)

aa1 <- std_table[codons_diff_clean[[1]]]
aa2 <- std_table[[codon2]]

AA1 <- lapply(codons_diff_clean[[1]], function(x) seqinr::getTrans(s2c(x)))
AA2 <- lapply(codons_diff_clean[[2]], function(x) seqinr::translate(s2c(x)))

count_syn_nonsyn_changes <- mapply(function(x, y) {
  c(ifelse(std_table[[x]] == std_table[[y]], "s", "n"),
  sum(s2c(x) != s2c(y)))
}, codons_diff_clean[[1]], codons_diff_clean[[2]])


res <- mapply(function(x, y) {
  type <- ifelse(std_table[[x]] == std_table[[y]], "s", "n")
  diffs <- sum(s2c(x) != s2c(y))
  list(type = type, diffs = diffs)
}, codons_diff_clean[[1]], codons_diff_clean[[2]], SIMPLIFY = FALSE)

# Extract type and diffs vectors
types <- vapply(res, `[[`, character(1), "type")
diffs <- vapply(res, `[[`, numeric(1), "diffs")

# synonymous_differences / synonymous_sites

names(table(types))
as.vector(tapply(diffs, types, sum)) / as.vector(table(types))

# Tabulate results
summary <- data.frame(
  type = names(table(types)),
  count = as.vector(table(types)),
  sum = as.vector(tapply(diffs, types, sum)),
  row.names = NULL
)







diff_pos <- which(codons1 != codons2)

codon1 = codons1[diff_pos]
codon2 = codons2[diff_pos]

# check for missing nucleotidesS
bad_idx <- unique(c(grep("[^ATGC]", codon1),
             grep("[^ATGC]", codon2)))

if (length(bad_idx) != 0) {
  codon1 = codon1[-bad_idx]
  codon2 = codon2[-bad_idx]
}

AA1 <- getTrans(s2c(paste0(codon1, collapse = "")), numcode = 11)
AA2 <- getTrans(s2c(paste0(codon2, collapse = "")), numcode = 11)

diff_aa <- which(AA1 != AA2)



which(mapply(`!=`, x[[1]], x[[2]]))



set = lapply(list(codon1, codon2), function(codon)
  paste0(codon, collapse = ""))

dna_str <- DNAStringSet(unlist(set))

tln_nt <- Biostrings::translate(dna_str)




dna <- DNAStringSet(c("ATGGCGTTCGAA", "ATGGCATTCAAG"))

# Translate each sequence
aa <- translate(dna)
aa

Ls <- sum(Ls_i)
Ln <- sum(Ln_i)




# Define a pairwise function (example: Hamming distance)
pairwise_fun <- function(s1, s2) {
  sum(strsplit(s1, "")[[1]] != strsplit(s2, "")[[1]])
}

# Generate all unique pairs (combinations of 2)
pairs <- combn(names(seqs), 2, simplify = FALSE)

# Apply your function to each pair using mapply
results <- mapply(function(p1, p2) {
  pairwise_fun(seqs[[p1]], seqs[[p2]])
}, sapply(pairs, `[`, 1), sapply(pairs, `[`, 2), SIMPLIFY = TRUE)

# Label results with pair names
names(results) <- sapply(pairs, function(p) paste(p, collapse = "_"))