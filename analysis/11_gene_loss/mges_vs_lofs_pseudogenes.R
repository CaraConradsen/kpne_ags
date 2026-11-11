
# Two-panel version of plot_seg_sites_per_site.R
# Panel a: all syntenic AGs (as before)
# Panel b: same, with MGE-associated AGs removed
# BH correction and Spearman's rho are recomputed within each panel, so panel b
# is a standalone test on the reduced set rather than a subset of panel a's calls.

alpha_fdr = 0.05

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

# Identify MGEs -----------------------------------------------------------
# Same COG source as plot_outliers_by_COG.R: egng_cog, renamed COG_funct_cat.

anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                     select = c("gene_family", "egng_cog",
                                "consensus_gene_name", "consensus_product")))

setnames(anno, old = "egng_cog", new = "COG_funct_cat")

# collapse to one row per family; flag any family carrying conflicting annotation
dup_fams <- anno[, .N, by = gene_family][N > 1, gene_family]
if (length(dup_fams))
  message("Families with >1 annotation row, first row taken: ", length(dup_fams))
anno <- unique(anno, by = "gene_family")

# COG X = Mobilome (prophages, transposons). COG strings can hold several
# letters, so search for X rather than matching exactly.
anno[, is_mge := grepl("X", COG_funct_cat, fixed = TRUE)]

# eggNOG assigns X conservatively, so also screen the product text.
# Set to FALSE to use the COG call alone.
use_product_screen <- TRUE

mge_regex <- paste0("transposase|transposon|integrase|excisionase|resolvase|",
                    "recombinase|insertion sequence|IS[0-9]|phage|prophage|",
                    "conjugal|conjugative|relaxase|mobilization|mobile element|",
                    "plasmid|replicase|type IV secretion")

# phage-defence genes are not themselves MGEs; keep them
mge_keep <- "anti-phage|phage resistance|phage defen|abortive infection|restriction"

if (use_product_screen) {
  anno[grepl(mge_regex, consensus_product, ignore.case = TRUE), is_mge := TRUE]
  anno[grepl(mge_keep, consensus_product, ignore.case = TRUE), is_mge := FALSE]
}

gene_pvals <- merge(gene_pvals, anno, by = "gene_family", all.x = TRUE)

# unannotated families are retained in panel b; note this in the legend
gene_pvals[is.na(is_mge), is_mge := FALSE]

cat("AGs tested:", nrow(gene_pvals),
    "| flagged MGE:", gene_pvals[is_mge == TRUE, .N],
    "| non-MGE:", gene_pvals[is_mge == FALSE, .N], "\n")

# add in lofs, collapse per family

lof_calls <- fread(paste0(outdir_dat, "/fis_fus_lof_calls.csv"))[lof_class != "no_lof_detected"]

lof_ags <- data.frame(gene_family = unique(lof_calls[, gene_family]),
              type = "lof")

setDT(lof_ags)

gene_pvals <- merge(gene_pvals, lof_ags,
                    all.x = TRUE, by = "gene_family")


gene_pvals[is.na(type), type := "not_pseudo"]


gene_pvals[, bin := cut(freq,
                        breaks = seq(min(freq), max(freq) + 5, by = 4),
                            right = FALSE)]

lof_plot <- gene_pvals[type == "lof",.(n = .N), "bin"][order(bin)]

lof_vec <- setNames(lof_plot$n, as.character(lof_plot$bin))

mge_plot <- gene_pvals[is_mge==TRUE,.(n = .N), "bin"][order(bin)]

mge_vec <- setNames(mge_plot$n, as.character(mge_plot$bin))


par(mfrow=c(2,1), mar = c(0.5,4,0.5,0.5))

barplot(lof_vec, las = 2, col = c("#0B3C49"),
        border = c("#0B3C49"),
        xaxt = "n")
axis(side = 2, at = 40, "Genes with at least one LOF event",
     tick = FALSE, line = 2)

par(mar = c(4,4,0,0.5))
bp <- barplot(mge_vec,
        ylim = c(max(mge_vec), 0),
        las = 2, col="tomato1",
        border = c("tomato1"),
        xaxt = "n",
        xlab = expression("Pangenome frequency ("*italic(k)*")"))
axis(side = 2, at = 150, "MGEs",
     tick = FALSE, line = 2)

axis(side = 1, at = seq(0, max(bp) + 0.5, length.out = 5),
     labels = round(seq(0,257, length.out = 5), digits = 0))

