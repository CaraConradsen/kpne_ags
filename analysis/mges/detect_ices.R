# Investigate ICEs in the pangenome, contrast this against the plasmid list



# Import and sort blastn hits against ICEberg 2.0 -------------------------
iceberg_files = list.files("./input_data/iceberg_kpne", 
                           recursive = TRUE,
                           full.names = TRUE, pattern = ".txt")


tmp_ice <- unique(fread(iceberg_files[i]))

tmp_ice <- tmp_ice[ , .SD[which.max(V12)], by = .(V1)]

tmp_ice[, c("V2", "accession") := tstrsplit(V2, "\\|", fill = TRUE, keep = c(3, 5))]

tmp_ice[, .(V1, V2, V3, V7, V8, accession)]


colnames(tmp_ice) = c("contig_id", "ICE", "pident", "start", "end", 
                      "accession")
