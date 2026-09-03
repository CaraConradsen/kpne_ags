ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

colnames(ag_age_S_dt)[10] = "freq"

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_s_sites),
                               freq = mean(freq)), gene_family]


# get piS theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]

# get expected distribution
theta_prod_lng <- fread(paste0(outdir_dat,"/pi_prod_lng.csv"))

# calculate p

# add expected S to gene_table
expected_by_freq <- theta_prod_lng[, .(expected_S = sum(seg_sites * p)), by = freq]

ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = "freq"]

# two-sided p via the "minimum tail doubled" 
gene_pvals <- ag_age_S_dt[theta_prod_lng, on = .(freq, seg_s_sites = seg_sites), 
                         nomatch = NULL]

# for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
gene_pvals[, p_lower := cdf]                    # P(X <= observed)
gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]


# Benjamini-Hochberg
gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]

fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_fdr.csv"))









