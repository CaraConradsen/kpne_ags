# will estimate the expected number of segregating sites per syntenic sites

ag_syn_sites <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"),
                         select = c("gene_family", "syn_sites"))

# estimates names for set 1 
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

ag_syn_sites <- ag_syn_sites[gene_family %chin% set_1_ags]

ag_syn_sites <- sort(unique(ag_syn_sites$syn_sites))

fwrite(list(ag_syn_sites), 
       paste0(outdir_dat, "/syntenic_ag_syn_sites.txt"), 
       col.names = FALSE)


# check spread for hpc
length(ag_syn_sites)
tail(ag_syn_sites)


for (t in 1:86) {
     idx <- ((t-1)*8+1):min(t*8, length(ag_syn_sites))
     cat("\t", t, ":", min(ag_syn_sites[idx]), "-", max(ag_syn_sites[idx]), "\n")
}
