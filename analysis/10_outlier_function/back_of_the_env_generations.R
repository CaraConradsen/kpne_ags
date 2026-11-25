gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))
gene_pvals[, call := "ns"]
gene_pvals[p_bh < alpha_fdr & seg_sites > expected_S, call := "high"]
gene_pvals[p_bh < alpha_fdr & seg_sites < expected_S, call := "low"]
gene_pvals[, call := factor(call, levels = c("ns", "high", "low"))]


high_out_ags <- gene_pvals[call=="high"]

high_out_ags[, .(mean = mean(pi_s), median = median(pi_s),
        min = min(pi_s), max = max(pi_s))]

high_out_ags[which.min(pi_s), .(gene_family, freq, seg_sites, syn_sites, pi_s)]
high_out_ags[which.max(pi_s), .(gene_family, freq, seg_sites, syn_sites, pi_s)]

high_out_ags[, age_years := pi_s / (6.34e-7)]
high_out_ags[, age_gens  := pi_s / (1.74e-9)]

high_out_ags[, .(mean_yr = mean(age_years), median_yr = median(age_years),
        min_yr = min(age_years), max_yr = max(age_years))]
high_out_ags[, .(mean_gen = mean(age_gens), median_gen = median(age_gens),
        min_gen = min(age_gens), max_gen = max(age_gens))]

high_out_ags[, summary(age_gens / age_years)]   # should be ~364 for every gene
