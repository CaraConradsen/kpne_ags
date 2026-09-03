# will estimate the expected number of segregating sites per gene length

ag_gene_lengths <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"),
                         select = c("gene_family", "gene_length"))

# estimates names for set 1 
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

ag_gene_lengths <- ag_gene_lengths[gene_family %chin% set_1_ags]

ag_gene_lengths <- sort(unique(ag_gene_lengths$gene_length))

fwrite(list(ag_gene_lengths), 
       paste0(outdir_dat, "/syntenic_ag_gene_lengths.txt"), 
       col.names = FALSE)


# check spread for hpc
length(ag_gene_lengths)
tail(ag_gene_lengths)
ag_gene_lengths[721:724]


for (t in 2:91) {
     idx <- ((t-1)*8+1):min(t*8, length(ag_gene_lengths))
     cat(t, ":", min(ag_gene_lengths[idx]), "-", max(ag_gene_lengths[idx]), "\n")
}
