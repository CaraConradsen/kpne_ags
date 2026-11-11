
# COG ENRICHMENT ANALYSIS, USING FISHER'S EXACT

# Import data -------------------------------------------------------------

set1_cogs <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                              select = c("gene_family", "egng_cog", "ag_type")))

# AG set used throughout
set1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

set1_cogs <- set1_cogs[gene_family %chin% set1_ags]

# # remove unknown genes
set1_cogs <- set1_cogs[egng_cog != ""]

#unlist COGs
set1_cogs <- set1_cogs[, .(egng_cog = unlist(strsplit(egng_cog, "")),
                           ag_type = ag_type), by = "gene_family"]

COGs <- COG_dt[COG_letter!="", COG_letter]

# Outlier calls by direction ----------------------------------------------

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites", "syn_sites",
                               "expected_S", "p_bh",
                               "p", "cdf"))  

gene_pvals <- merge(gene_pvals,
                    unique(set1_cogs[, .(gene_family, ag_type)]),
                    all.y = TRUE, by = "gene_family")

gene_pvals[, s_per_site   := seg_sites / syn_sites]
gene_pvals[, log2_obs_exp := log2((seg_sites + 0.5) / (expected_S + 0.5))]

cog_metrics <- c("s_per_site", "log2_obs_exp")

min_band_n <- 10   # smallest category x band cell worth testing

gene_pvals[, call := fcase(p_bh <= 0.05 & seg_sites > expected_S, "high",
                           p_bh <= 0.05 & seg_sites < expected_S, "low",
                           default = "ns")]

high_ags <- gene_pvals[call == "high", gene_family]
low_ags  <- gene_pvals[call == "low",  gene_family]

# COG enrichment among outliers -------------------------------------------

cog_annotated <- unique(set1_cogs$gene_family)

ag_call <- gene_pvals[gene_family %chin% cog_annotated]

min_cog_n <- 5   # smallest category total across the two arms worth testing

fisher_cogs <- function(focal_ags, bg_ags, label) {

  res_dt <- rbindlist(lapply(COGs, function(cog){

    in_cog <- set1_cogs[egng_cog == cog, gene_family]

    a  <- sum(focal_ags %chin% in_cog)      # focal, in this COG
    b  <- length(focal_ags) - a             # focal, not in this COG
    cc <- sum(bg_ags %chin% in_cog)         # background, in this COG
    d  <- length(bg_ags) - cc

    if (a + cc < min_cog_n) return(NULL)

    cog_mat <- matrix(c(a, b, cc, d), ncol = 2, byrow = TRUE,
                      dimnames = list(c("cog", "other"),
                                      c("focal", "background")))

    res <- fisher.test(cog_mat)

    data.table(
      comparison  = label,
      cog         = cog,
      n_focal     = a,
      n_focal_tot = length(focal_ags),
      n_bg        = cc,
      n_bg_tot    = length(bg_ags),
      OR          = unname(res$estimate),
      CI_lower    = round(res$conf.int[1], 4),
      CI_upper    = round(res$conf.int[2], 4),
      pval        = round(res$p.value, 4))
  }))

  if (!nrow(res_dt)) return(NULL)

  # BH within each comparison: each is its own family of tests
  res_dt[, p_bh := round(p.adjust(pval, "BH"), 4)][]
}

run_both_directions <- function(dt, suffix = "") {
  rbind(
    fisher_cogs(dt[call == "high", gene_family],
                dt[call == "ns",   gene_family],
                paste0("high vs non-outlier", suffix)),
    if (dt[call == "low", .N] >= min_cog_n)
      fisher_cogs(dt[call == "low", gene_family],
                  dt[call == "ns",  gene_family],
                  paste0("low vs non-outlier", suffix)),
    fill = TRUE)
}

# overall, then the same comparison within each pangenome frequency band
fisher_out <- rbind(
  run_both_directions(ag_call),
  rbindlist(lapply(unique(ag_call$ag_type), function(b)
    run_both_directions(ag_call[ag_type == b], paste0(", ", b)))),
  fill = TRUE)

fisher_out[, direction := fcase(is.na(p_bh) | p_bh > 0.05, "",
                                OR > 1, "enriched", default = "depleted")]

cat("\nAGs entering the comparisons:\n")
print(ag_call[, .N, by = .(call, ag_type)][order(call, ag_type)])

cat("\nCOG categories at BH 0.05:\n")
print(fisher_out[direction != "",
                 .(comparison, cog, n_focal, n_bg, OR = round(OR, 2),
                   p_bh, direction)][order(comparison, -OR)])

fwrite(fisher_out, paste0(outdir_dat, "/cog_fisher_outliers.csv"))

