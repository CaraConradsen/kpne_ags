

# import all the gene probability info, from hpc, previously separ --------

ag_m_prob_files <- list.files(paste0(outdir_dat, "/prod_ag_m/"), 
                              pattern = "^pi_prod_lng_m[0-9]+\\.csv$", 
                              full.names = TRUE)

theta_prod_lng <- lapply(ag_m_prob_files, function(x){
  gene_length = as.integer(sub(".*_m([0-9]+)\\.csv$", "\\1", x))
  cbind(fread(x), gene_length)
})

theta_prod_lng <- rbindlist(theta_prod_lng)

# estimates names for syntenic ags
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

ag_age_S_dt <- ag_age_S_dt[gene_family %chin% set_1_ags]



# get piS theta probability distribution of expected seg sites ------------

# calculate p

# add expected S to gene_table
expected_by_freq <- theta_prod_lng[, .(expected_S = sum(seg_sites * p)), by = c("freq", "gene_length")]

ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = .(freq, gene_length)]

# two-sided p via the "minimum tail doubled" 
gene_pvals <- ag_age_S_dt[theta_prod_lng, on = .(freq, gene_length, seg_sites), 
                         nomatch = NULL]

# for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
gene_pvals[, p_lower := cdf]                    # P(X <= observed)
gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]


# Benjamini-Hochberg
gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]

fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_fdr.csv"))









