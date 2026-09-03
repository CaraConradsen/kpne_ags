
# COG (and Environment) ENRICHMENT ANALYSIS, USING FISHER'S EXACT
# watch out for double counting with the multi-identified COGs,
# exclude AGs with unassigned cogs


# Import data -------------------------------------------------------------

set1_cogs <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("gene_family", "egng_cog")))

# get all AGs in both sets
set1_ags <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/set_1_names.csv")$gene_family
# set2_ags <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/set_2_names.csv")$gene_family

set1_cogs <- set1_cogs[gene_family %chin% set1_ags]

# # remove unknown genes
set1_cogs <- set1_cogs[egng_cog != ""]

#unlist COGs
set1_cogs <- set1_cogs[, .(egng_cog = unlist(strsplit(egng_cog, ""))), by = "gene_family"]

selected_ags <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]$gene_family

COGs <- COG_dt[COG_letter!="", COG_letter]



# loop over cogs ----------------------------------------------------------

fexact_cogs <- lapply(COGs, function(cog){
  
  # set cog groups up
  focal_cog_ags = set1_cogs[egng_cog==cog, gene_family]
  others = unique(set1_cogs[!gene_family %chin% focal_cog_ags, gene_family])
  
  cog_sel <- length(focal_cog_ags[focal_cog_ags %chin% selected_ags])
  cog_unsel <- length(focal_cog_ags[!focal_cog_ags %chin% selected_ags])
  
  othr_sel <- length(others[others %chin% selected_ags])
  othr_unsel <- length(others[!others %chin% selected_ags])
  
  cog_mat = matrix(c(cog_sel,cog_unsel,
                     othr_sel, othr_unsel),
                   ncol = 2, byrow = TRUE,
                   dimnames = list(c("cog", "other"),
                                   c("sel", "unsel")))
  
  res <- fisher.test(cog_mat)
  data.table(
    cog      = cog,
    pval     = round(res$p.value, 4),
    OR       = res$estimate,
    CI_lower = round(res$conf.int[1], 4),
    CI_upper = round(res$conf.int[2], 4)
  )
  
})


fexact_cogs <- rbindlist(fexact_cogs)

fexact_cogs$p_bh <- round(p.adjust(fexact_cogs$pval, "BH"), 4)



# Export csv --------------------------------------------------------------

setorderv(fexact_cogs, cols = c("pval", "cog"),
          order = c(1, 1))

fexact_cogs[, pval := sprintf("%1.3f", pval)]
fexact_cogs[, p_bh  := sprintf("%1.3f", p_bh)]
fexact_cogs[, OR  := sprintf("%1.2f", OR)]
# fix small vals
fexact_cogs[pval == "0.000", pval := "<0.001"]
fexact_cogs[p_bh == "0.000", p_bh := "<0.001"]

fwrite(fexact_cogs, paste0(outdir_dat, "/fexact_cogs.csv"))



