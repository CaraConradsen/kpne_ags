ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

colnames(ag_age_S_dt)[10] = "freq"

ag_age_S_dt <- ag_age_S_dt[, .(seg_s_sites = mean(seg_s_sites),
                               freq = mean(freq)), gene_family]

# get waterson's theta probability distribution of expected seg sites
# round up seg sites
ag_age_S_dt[, seg_s_sites := round(seg_s_sites)]

# get expected distribution
theta_prod_lng_list <- list.files(paste0(outdir_dat, "/theta_pis_out"),
                                  pattern ="prob_pi_lng", 
                                  full.names = TRUE, recursive = TRUE)

theta_prod_lng <- lapply(theta_prod_lng_list, function(lst){
  theta_dat <- fread(lst)
  theta_dat$theta_val <- as.numeric(gsub("prob_pi_lng_", "", 
                                      gsub(".csv", "", basename(lst))))
  
  theta_dat
}) 

theta_prod_lng <- rbindlist(theta_prod_lng)

theta_vals = unique(theta_prod_lng$theta_val)

# calculate p

# add expected S to gene_table
expected_by_freq <- theta_prod_lng[, .(expected_S = sum(seg_sites * p)), 
                                   by = c("freq", "theta_val")]

ag_age_S_lst <- lapply(theta_vals, function(theta){
  temp_exp <- expected_by_freq[theta_val == theta]
  temp_exp[ag_age_S_dt, on = "freq"]
})
  
ag_age_S_dt  <- rbindlist(ag_age_S_lst)

# two-sided p via the "minimum tail doubled" 
gene_pvals <- ag_age_S_dt[theta_prod_lng, on = .(freq, theta_val, seg_s_sites = seg_sites), 
                         nomatch = NULL]

# for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
gene_pvals[, p_lower := cdf]                    # P(X <= observed)
gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]


# Benjamini-Hochberg
gene_pvals[, p_bh := p.adjust(p_two, method = "BH"), theta_val]

fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_mutli_pis_fdr.csv"))









