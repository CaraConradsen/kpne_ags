# 
# 
# # # import all the gene probability info, from hpc, previously separ --------
# # 
# # ag_m_prob_files <- list.files(paste0(outdir_dat, "/prod_ag_m/"), 
# #                               pattern = "^pi_prod_lng_m[0-9]+\\.csv$", 
# #                               full.names = TRUE)
# # 
# # theta_prod_m_lng <- lapply(ag_m_prob_files, function(x){
# #   gene_length = as.integer(sub(".*_m([0-9]+)\\.csv$", "\\1", x))
# #   cbind(fread(x), gene_length)
# # })
# # 
# # theta_prod_m_lng <- rbindlist(theta_prod_m_lng)
# 
# # estimates names for syntenic ags
# set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family
# 
# ag_age_S_dt <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))
# 
# ag_age_S_dt <- ag_age_S_dt[gene_family %chin% set_1_ags]
# 
# 
# 
# # get piS theta probability distribution of expected seg sites ------------
# 
# # calculate p
# 
# # add expected S to gene_table
# expected_by_freq_m <- theta_prod_m_lng[, .(expected_S = sum(seg_sites * p)), by = c("freq", "gene_length")]
# 
# ag_age_S_m_dt <- expected_by_freq_m[ag_age_S_dt, on = .(freq, gene_length)]
# 
# # two-sided p via the "minimum tail doubled" 
# m_gene_pvals <- ag_age_S_m_dt[theta_prod_m_lng, on = .(freq, gene_length, seg_sites), 
#                               nomatch = NULL]
# 
# # for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
# m_gene_pvals[, p_lower := cdf]                    # P(X <= observed)
# m_gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
# m_gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]
# 
# 
# # Benjamini-Hochberg
# m_gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]
# 
# # fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_fdr.csv"))
# 



#Ls
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


# Ls !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# get piS theta probability distribution of expected seg sites ------------

# calculate p

# add expected S to gene_table
expected_by_freq <- theta_prod_Ls_lng[, .(expected_S = sum(seg_sites * p)), by = c("freq", "syn_sites")]

ag_age_S_dt <- expected_by_freq[ag_age_S_dt, on = .(freq, syn_sites)]

# two-sided p via the "minimum tail doubled" 
gene_pvals <- ag_age_S_dt[theta_prod_Ls_lng, on = .(freq, syn_sites, seg_sites), 
                          nomatch = NULL]

# for each gene, p_two_sided = 2 * min(P(X <= obs), P(X >= obs))
gene_pvals[, p_lower := cdf]                    # P(X <= observed)
gene_pvals[, p_upper := 1 - cdf + p]          # P(X >= observed)
gene_pvals[, p_two  := pmin(1, 2 * pmin(p_lower, p_upper))]


# Benjamini-Hochberg
gene_pvals[, p_bh := p.adjust(p_two, method = "BH")]

fwrite(gene_pvals, paste0(outdir_dat, "/gene_pvals_fdr.csv"))




# TO REMOVE, just checking ------------------------------------------------
# setnames(m_gene_pvals, "expected_S", "expected_S_m")
# 
# temp_dat <- merge(gene_pvals, m_gene_pvals[,.(gene_family, expected_S_m)],
#                   all.x = TRUE, by = "gene_family")
# 
# 
# temp_dat[,.(gene_family, freq, gene_length, syn_sites,seg_sites, expected_S, expected_S_m)] %>% head
# 
# with(temp_dat,
#      plot(expected_S, expected_S_m,
#           xlab = "E[s] Ls", ylab = "E[s] m",
#           bty = "L",
#           pch = 16, col = rgb(0.3,0.3,1, alpha = 0.5)))
# 
# abline(a=0,b=1, lty = 2)


