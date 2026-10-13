
# For Ls: import all the gene probability info, from hpc, previously separated --------

ag_Ls_prob_files <- list.files(paste0(outdir_dat, "/prod_ag_Ls_2.92/"), 
                               pattern = "^pi_prod_lng_m[0-9]+\\.csv$", 
                               full.names = TRUE)

theta_prod_Ls_lng <- lapply(ag_Ls_prob_files, function(x){
  syn_sites = as.integer(sub(".*_m([0-9]+)\\.csv$", "\\1", x))
  cbind(fread(x), syn_sites)
})

theta_prod_Ls_lng <- rbindlist(theta_prod_Ls_lng)

# estimates names for syntenic ags
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))

ag_age_S_dt <- ag_age_S_dt[gene_family %chin% set_1_ags]


# get piS theta probability distribution of expected seg sites ------------
# add expected S to gene_table
expected_by_freq <- theta_prod_Ls_lng[, .(expected_S = sum(seg_sites * p)), by = c("freq", "syn_sites")]

ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = .(freq, syn_sites)]

# left join: keep every gene, even if seg_sites is beyond the truncated range
gene_pvals <- theta_prod_Ls_lng[ag_age_S_dt, on = .(freq, syn_sites, seg_sites)]

# p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
gene_pvals[, p_lower := cdf]
gene_pvals[, p_upper := 1 - cdf + p]

# beyond truncation: P(X >= obs) < tail mass (~1e-10)
gene_pvals[is.na(p) & seg_sites > expected_S, `:=`(p_lower = 1, p_upper = 1e-10)]

gene_pvals[, p_two := pmin(1, 2 * pmin(p_lower, p_upper))]

gene_pvals[is.na(p_two)]

# Benjamini-Hochberg
gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]

fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_fdr.csv"))
