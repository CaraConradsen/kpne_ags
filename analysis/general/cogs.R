library(data.table)

COG_dt <- data.frame(
  COG_letter = c(
    "A","B","C","D","E","F","G","H","I","J",
    "K","L","M","N","O","P","Q","R","S","T",
    "U","V","W","X","Y","Z"
  ),
  COG_function = c(
    "RNA processing and modification",
    "Chromatin structure and dynamics",
    "Energy production and conversion",
    "Cell cycle control, cell division, chromosome partitioning",
    "Amino acid transport and metabolism",
    "Nucleotide transport and metabolism",
    "Carbohydrate transport and metabolism",
    "Coenzyme transport and metabolism",
    "Lipid transport and metabolism",
    "Translation, ribosomal structure and biogenesis",
    "Transcription",
    "Replication, recombination and repair",
    "Cell wall/membrane/envelope biogenesis",
    "Cell motility",
    "Post-translational modification, protein turnover, chaperones",
    "Inorganic ion transport and metabolism",
    "Secondary metabolites biosynthesis, transport, and catabolism",
    "General function prediction only",
    "Function unknown",
    "Signal transduction mechanisms",
    "Intracellular trafficking, secretion, and vesicular transport",
    "Defense mechanisms",
    "Extracellular structures",
    "Mobilome",
    "Nuclear structure",
    "Cytoskeleton"
  ),
  COG_num = 1:26,
  stringsAsFactors = FALSE
)

setDT(COG_dt)