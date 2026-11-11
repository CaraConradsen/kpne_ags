# create and export an xlxs worksheet for examining outlier AGs

outlier_ags <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                              select = c("gene_family", "freq","egng_cog", "ag_type","gene_length",
                                         "consensus_product","consensus_gene_name",
                                         "tidy_name", "phg_state","most_common_phage",
                                         "ice","is_fam"))

cols <- c("phg_state", "most_common_phage", "ice", "is_fam")
outlier_ags[, (cols) := lapply(.SD, function(x) fifelse(x == "", NA_character_, x)), .SDcols = cols]

count_str <- function(x) {
  tb <- sort(table(x[!is.na(x)]), decreasing = TRUE)
  paste(c(if (length(tb)) paste(tb, names(tb))),
        collapse = " | ")
}


# most common non-NA value (NA if none)
mode_chr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  names(which.max(table(x)))
}

phg_levels <- c("intact", "questionable", "incomplete")

outlier_sum <- outlier_ags[, .(
  phg_state = paste(c(paste(vapply(phg_levels, function(l) sum(phg_state == l, na.rm = TRUE), integer(1)),
                            phg_levels)),
                    collapse = " | "),
  most_common_phage = mode_chr(most_common_phage),
  ice               = count_str(ice),
  is_fam            = count_str(is_fam)
), by = .(gene_family, freq, egng_cog, ag_type, gene_length, consensus_product, tidy_name, consensus_gene_name)]

# clean phg_state
outlier_sum[is.na(most_common_phage), phg_state := NA]



# take longest gene
outlier_ags <- unique(outlier_sum[, .SD[which.max(gene_length)], by = gene_family])


gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "seg_sites", "syn_sites", 
                               "expected_S", "pi_s","p_two","p_bh"))

## three-way call: high, low, or not significant
alpha_fdr = 0.05
gene_pvals[, call := "ns"]
gene_pvals[p_bh < alpha_fdr & seg_sites > expected_S, call := "high"]
gene_pvals[p_bh < alpha_fdr & seg_sites < expected_S, call := "low"]
gene_pvals[, call := factor(call, levels = c("ns", "high", "low"))]

outlier_ags <- merge(gene_pvals[call != "ns"],
                     outlier_ags, all.x = TRUE, 
                  by = "gene_family")

setcolorder(outlier_ags, c( "call", "gene_family","pi_s", "freq", 
                            "ag_type","egng_cog", "tidy_name", "consensus_product",  "consensus_gene_name",
                            "gene_length", "syn_sites", "seg_sites", "expected_S",  "p_two","p_bh"))

outlier_ags[, pi_s := sprintf("%1.4f", pi_s)]
outlier_ags[, p_two := sprintf("%1.4f", p_two)]
outlier_ags[, p_bh  := sprintf("%1.4f", p_bh)]
# fix small vals
outlier_ags[p_two == "0.0000", p_two := "<0.0001"]
outlier_ags[p_bh == "0.0000", p_bh := "<0.0001"]


# associations ------------------------------------------------------------
msu_anchor <- unique(fread("./output/data/msu_regions_anchored.csv", 
                           select = c("gene_family", "msu", "anchor")))

outlier_ags <- merge(outlier_ags,
                     msu_anchor, all.x = TRUE, 
                     by = "gene_family")

clean_coin_ags <- fread(paste0(outdir_dat,"/clean_coin_ags.csv"),
                        select = c("gene_family", "clust_id", "score"))

outlier_ags <- merge(outlier_ags,
                     clean_coin_ags, all.x = TRUE, 
                     by = "gene_family")

setorderv(outlier_ags, cols = c("call", "freq"),
          order = c(-1, 1))

setnames(outlier_ags,
         c( "call", "gene_family","pi_s", "freq", 
            "ag_type","egng_cog", "tidy_name", "consensus_product",  
            "gene_length", "syn_sites", "seg_sites", "expected_S",  "p_two","p_bh", 
            "msu", "anchor", "clust_id",  "score", 
            "phg_state","most_common_phage",
            "ice","is_fam"),
         
         c( "Outlier type", "Gene","pi_s", "Gene Freq.", 
            "Pangenome category","COG", "Tidy gene name", "Consensus product",  
            "Longest gene length", "Synonymous sites", "Segregating synonymous sites, S", "E[S]",  "P value","FDR Q",
            "Core-block number", "Flanking core genes", "Co-occurence group",  "Detection method", 
            "Phage Completeness","Phage Type",
            "ICE","Insertion sequence Family"))




# Create a coloured xlxs workbook for analysis ----------------------------

wb <- openxlsx::createWorkbook()

# Define styles
orange_sig    <- createStyle(fgFill = "#D0A855")
blue_sig      <- createStyle(fgFill = "#7FA8D4")
bold_header   <- createStyle(textDecoration = "bold")


# Sheet 1: Data dictionary ------------------------------------------------
dict_outliers <- data.table(
  Column = c(
    "Gene", "Outlier type", "pi_s", "Gene Freq.", "Pangenome category", "COG",
    "Tidy gene name", "Consensus product", "Longest gene length",
    "Synonymous sites", "Segregating synonymous sites, S", "E[S]",
    "P value", "FDR Q",
    "Phage Completeness", "Phage Type", "ICE", "Insertion sequence Family",
    "Core-block number", "Flanking core genes",
    "Co-occurence group", "Detection method"
  ),
  Description = c(
    "PIRATE gene family identifier",
    "Direction of the outlier: high (S > E[S], excess synonymous diversity) or low (S < E[S], deficit)",
    "Synonymous nucleotide diversity (piS) across 4-fold and 2-fold degenerate sites",
    "Number of genomes carrying an intact copy (fission/fusion loci excluded)",
    "Pangenome category from PIRATE frequency: core (>= 99%), soft (95-99%), shell (15-95%), cloud (< 15%)",
    "eggNOG COG functional category letter code(s)",
    "Gene name, or shortened product annotation where no gene name is assigned",
    "Consensus product annotation from PIRATE (Bakta)",
    "Length (bp) of the longest representative gene in the family",
    "Number of synonymous (4-fold and 2-fold) sites assayed",
    "Observed number of segregating synonymous sites",
    "Expected segregating synonymous sites under the neutral coalescent null (msprime, demography fitted to the core genome)",
    "Two-tailed p-value from the selection test",
    "Benjamini-Hochberg corrected p-value",
    "Number of genomes in which the gene lies within a PHASTEST prophage region, by region completeness (intact | questionable | incomplete); NA: not within a prophage in any genome",
    "Most common PHASTEST prophage hit among genomes where the gene lies in a prophage",
    "Number of genomes in which the gene lies within each ICEfinder-predicted ICE",
    "Number of genomes in which the gene lies within each ISEScan insertion sequence family",
    "Syntenic core block (MSU) containing the gene",
    "Pair of core anchor genes flanking the insertion site",
    "Goldfinder co-occurrence cluster identifier",
    "Goldfinder scoring method: simultaneous, subsequent, or both"
  )
)
dict_outliers[, `:=`(Sheet  = "BH corrected selection outliers",
                     Colour = "Blue: high outlier; Orange: low outlier")]
setcolorder(dict_outliers, c("Sheet", "Column", "Description", "Colour"))

addWorksheet(wb, "Data dictionary")
openxlsx::writeData(wb, "Data dictionary", dict_outliers)
addStyle(wb, "Data dictionary", bold_header,
         rows = 1, cols = 1:ncol(dict_outliers), gridExpand = TRUE)

# Sheet 2: BH corrected selection outliers --------------------------------
addWorksheet(wb, "BH corrected selection outliers")
openxlsx::writeData(wb, "BH corrected selection outliers", outlier_ags)

dir_rows <- which(outlier_ags[["Outlier type"]] == "low") + 1
bal_rows <- which(outlier_ags[["Outlier type"]] != "low") + 1

addStyle(wb, "BH corrected selection outliers", orange_sig, rows = dir_rows,
         cols = 1:ncol(outlier_ags), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "BH corrected selection outliers", blue_sig, rows = bal_rows,
         cols = 1:ncol(outlier_ags), gridExpand = TRUE, stack = TRUE)
addStyle(wb, "BH corrected selection outliers", bold_header,
         rows = 1, cols = 1:ncol(outlier_ags), gridExpand = TRUE)

saveWorkbook(wb, paste0(outdir_xlsx, "/outliers_ag.xlsx"), overwrite = TRUE)
