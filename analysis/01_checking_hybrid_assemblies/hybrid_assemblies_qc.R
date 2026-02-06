
hybrid_files = list.files("C:/Users/carac/Dropbox/Vos_Lab/Spark_plasmids/1150_SPARK_hybrid_assemblies/1150_SPARK_hybrid_assemblies",
                          pattern = "fasta", full.names = TRUE)


hybrid_contig_info <- foreach(i = hybrid_files,
                              .packages = c("Biostrings"),
                              .combine = "rbind") %do% {
                                
                                fasta_file = i
                                
                                genome = gsub(".fasta","", basename(i))
                                
                                contigs = Biostrings::readDNAStringSet(fasta_file)
                                
                                data.frame(filename = genome,
                                           header = names(contigs),
                                           length = width(contigs))
                              }

# tidy info

setDT(hybrid_contig_info)

hybrid_contig_info[, seqnames := tstrsplit(header, " ", fill = TRUE, keep = 1)]

hybrid_contig_info[!grepl("SPARK", seqnames), 
                   seqnames := paste(filename, seqnames, sep="_")]


dnaA_tsv = list.files("C:/Users/carac/Dropbox/Vos_Lab/Spark_plasmids/1150_SPARK_hybrid_assemblies/1150_SPARK_hybrid_assemblies/blast_results/",
                      pattern = "tsv", full.names = TRUE)

dna_loc <- foreach(t = dnaA_tsv,
                   .packages = "data.table",
                   .combine = "rbind") %do% {
                     
                     genome = gsub("_dnaA.tsv","", basename(t))
                     
                     cbind(genome, fread(t, select = c("V2", "V9", "V10")))
                   }

colnames(dna_loc) = c("filename", "seqnames", "dnaA_start", "dnaA_end")

# tidy names
dna_loc[!grepl("SPARK", seqnames), 
                   seqnames := paste(filename, seqnames, sep="_")]

dna_loc[,filename:=NULL]

hybrid_contig_info <- merge(hybrid_contig_info, dna_loc,
                            all.x=TRUE, by = "seqnames")

hybrid_contig_info[, n_contig:= .N, filename]

hybrid_contig_info[, max_contig := ifelse(.I %in% .I[which.max(length)], "Y", NA), by = filename]


# get kpne

spark_meta = fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv",
                   select = c("id", "species_abbv", "ST"))


colnames(spark_meta)[1] = "filename"

hybrid_contig_info <- merge(hybrid_contig_info, spark_meta,
                            all.x=TRUE, by = "filename")



# Export K.pne circular chr -----------------------------------------------

kpne_chr_circ <- hybrid_contig_info[dnaA_start==1 & species_abbv=="K.pne"]

# select 1 colony from same plate

kpne_chr_circ[, plate := sub("_C[1-9]$", "", filename)]

# max length per plate
kpne_chr_circ[, max_plate := ifelse(.I %in% .I[which.max(length)], "Y", "N"), by = plate]

kpne_chr_circ <- kpne_chr_circ[max_plate == "Y"]# 412

hybrid_dir = "C:/Users/carac/Dropbox/Vos_Lab/Spark_plasmids/1150_SPARK_hybrid_assemblies/1150_SPARK_hybrid_assemblies/"

out_dir = "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_412_chr_fasta/"

for (geno in kpne_chr_circ$filename) {
  cat("\nprocessing...", geno)
  
  keep_contig = kpne_chr_circ[filename==geno, seqnames] 
  
  chr_contig = Biostrings::readDNAStringSet(paste0(hybrid_dir,geno,".fasta"))
  
  chr_contig <- chr_contig[grepl("_1($|\\s)", names(chr_contig))]
  
  if(length(chr_contig)>1){
    cat("; num of contigs = ", length(chr_contig))
    break
  }
  
  if(length(chr_contig)==0){
    
    # re-read fasta
    chr_contig = Biostrings::readDNAStringSet(paste0(hybrid_dir,geno,".fasta"))
    
    # Get index of the largest contig
    max_idx <- which.max(width(chr_contig))
    
    # Subset only the largest contig
    chr_contig <- chr_contig[max_idx]
  }
  
  if(length(chr_contig)==0){
    cat("; num of contigs = ", length(chr_contig))
    break
  }
  
  # clean contig name
  names(chr_contig) = keep_contig
  
  writeXStringSet(chr_contig, 
                  filepath = paste0(out_dir, geno, ".fasta"), 
                  format="fasta")
}


