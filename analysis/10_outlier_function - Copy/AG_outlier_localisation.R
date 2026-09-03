
cnv_info <- fread(paste0(outdir_dat, "/cnv_info.csv"))

ag_outliers <- fread(paste0(outdir_dat, "/set_outliers.csv"))[theta_outlier==1]

cnv_info <- cnv_info[gene_family %in% ag_outliers$gene_family]

gene_pvals <- fread(paste0(outdir_dat, "/set_outliers.csv"),
                    select = c("gene_family","seg_s_sites", "p_two", "p_bh","w_theta_sel", "theta_outlier"))[theta_outlier==1]

out_cnv <- cnv_info[gene_family %chin% gene_pvals$gene_family & gene_family %chin% cnv_info[cnv_on_plsmd==1, gene_family]]

out_cnv <- merge(out_cnv, gene_pvals, 
                 all.x = TRUE, by = "gene_family")

write.xlsx(out_cnv, paste0(outdir_xlsx, "/AG_outliers.xlsx"), sheetName="cnv", append = TRUE)




# Single gene with a cnv
unique(cnv_info[cnv_on_plsmd>0, .(gene_family)])
# gene_family
# <char>
#   1:     g002647


# Inspect the cnv, g002647 ------------------------------------------------
core_gub_tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_bootstrapped_tree.tre")
core_gub_tree <- midpoint(core_gub_tree)

states <- setNames(rep("absent", length(core_gub_tree$tip.label)), core_gub_tree$tip.label)
states[cnv_info[gene_family == "g002647" & cnv_on_plsmd == 0, geno_id]] <- "no_cnv"
states[cnv_info[gene_family == "g002647" & cnv_on_plsmd != 0, geno_id]] <- "cnv"
states <- states[core_gub_tree$tip.label]

cols <- c(absent = "white", no_cnv = "black", cnv = "red")

# ST aligned to tip order, but NOT written into tip.label
pan_anno  <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))
ST_labels <- unique(pan_anno[geno_id %chin% core_gub_tree$tip.label, .(geno_id, ST)])
ST_labels <- ST_labels[match(core_gub_tree$tip.label, geno_id)]
st_vec    <- ST_labels$ST

# plot without the built-in labels, then add ST text and state symbols ourselves
ape::plot.phylo(core_gub_tree, show.tip.label = FALSE,
                no.margin = TRUE, type = "fan")

tiplabels(pch = 21, bg = cols[as.character(states)], cex = 0.9)          # state ring

# pull tip coordinates from the last plot
lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
n      <- length(core_gub_tree$tip.label)
tip_x  <- lastPP$xx[1:n]
tip_y  <- lastPP$yy[1:n]

# angle of each tip relative to the fan centre, in degrees
angles <- atan2(tip_y, tip_x) * 180 / pi

# draw each label individually with its own rotation
tiplabels(pch = 21, bg = cols[as.character(states)], cex = 0.7)
for (i in seq_len(n)) {
  tiplabels(text = st_vec[i], tip = i, frame = "none", cex = 0.4,
            offset = 1500, srt = angles[i])
}
legend("bottomleft",
       legend = c("AG absent", "AG present, no CNV", "AG present + plasmid CNV"),
       pt.bg = cols, pch = 21, bty = "n")


# Environment? ------------------------------------------------------------
spark_metadata <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv")

spark_metadata <- spark_metadata[id %chin% names(states)]

spark_metadata[id %chin% names(states[states=="cnv"])]
spark_metadata[id %chin% names(states[states=="absent"])]
spark_metadata[id %chin% names(states[states=="no_cnv"])]


# Location? ---------------------------------------------------------------
# most_common_ord_msu <- readLines("most_common_ord_msu.txt")
# ag_outliers[w_theta_sel=="directional" & V2 == "set_2"]

ag_synteny <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))

# genome.blocks <- ag_synteny[msu=="MSU_3" & anchor %chin% c("A143","A143:A144", "A144"),
#                             .(geno_id, gene_family, start, end,
#                               strand, anchor, consensus_product)]
# genome.blocks <- ag_synteny[msu=="MSU_13" & anchor %chin% c("A30","A30:A31", "A31"),
#                             .(geno_id, gene_family, start, end,
#                               strand, anchor, consensus_product)]
# genome.blocks <- ag_synteny[msu=="MSU_6" & anchor %chin% c("A138","A138:A139", "A139"),
#                             .(geno_id, gene_family, start, end, 
#                               strand, anchor, consensus_product)]

# balancing selection
genome.blocks <- ag_synteny[msu=="MSU_2" & anchor %chin% c("A339","A339:A340", "A340"),
                            .(geno_id, gene_family, start, end, 
                              strand, anchor, consensus_product)]
genome.blocks <- ag_synteny[msu=="MSU_1" & anchor %chin% c("A22","A22:A23", "A23"),
                            .(geno_id, gene_family, start, end, 
                              strand, anchor, consensus_product)]

genome.blocks <- ag_synteny[msu=="MSU_0" & anchor %chin% c("A113","A113:A114", "A114"),
                            .(geno_id, gene_family, start, end, 
                              strand, anchor, consensus_product)]

genome.blocks_cols = data.frame(gene_family = unique(genome.blocks$gene_family),
                      cols = colorRampPalette(c("#420D55","blueviolet",
                                                "blue", "#799000",
                                                "yellow","#FF8000",
                                                "brown", "salmon"))(length(unique(genome.blocks$gene_family))))

genome.blocks <- genome.blocks[genome.blocks_cols, on = "gene_family"]

pdf(paste0(outdir_fig,"/gene_sections.pdf"),
    width  = 17,
    height = 13)
plot_gene_sections(genome.blocks, order = "dendro")
dev.off()



# the graveyard -----------------------------------------------------------
states2 <- setNames(rep("absent", length(core_gub_tree$tip.label)), core_gub_tree$tip.label)
states2[ag_synteny[gene_family=="g002959", geno_id]] <- "island"
states2[ag_synteny[gene_family=="g006097", geno_id]] <- "g1"

states2 <- states2[core_gub_tree$tip.label]

cols <- c(absent = "white", g1 = "black", island = "forestgreen")

# ST aligned to tip order, but NOT written into tip.label
pan_anno  <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))
ST_labels <- unique(pan_anno[geno_id %chin% core_gub_tree$tip.label, .(geno_id, ST)])
ST_labels <- ST_labels[match(core_gub_tree$tip.label, geno_id)]
st_vec    <- ST_labels$ST

# plot without the built-in labels, then add ST text and state symbols ourselves
ape::plot.phylo(core_gub_tree, show.tip.label = FALSE,
                no.margin = TRUE, type = "fan")

tiplabels(pch = 21, bg = cols[as.character(states2)], cex = 0.9)          # state ring

# pull tip coordinates from the last plot
lastPP <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
n      <- length(core_gub_tree$tip.label)
tip_x  <- lastPP$xx[1:n]
tip_y  <- lastPP$yy[1:n]

# angle of each tip relative to the fan centre, in degrees
angles <- atan2(tip_y, tip_x) * 180 / pi

# draw each label individually with its own rotation
tiplabels(pch = 21, bg = cols[as.character(states2)], cex = 0.7)
for (i in seq_len(n)) {
  tiplabels(text = st_vec[i], tip = i, frame = "none", cex = 0.4,
            offset = 1500, srt = angles[i])
}
legend("bottomleft",
       legend = c("absent", "Phage graveyard", "Metabolic island"),
       pt.bg = cols, pch = 21, bty = "n")

spark_metadata <- fread("C:/Users/carac/Dropbox/Vos_Lab/SpARK data/spark_metadata.csv")

spark_metadata <- spark_metadata[id %chin% names(states2)]

spark_metadata[id %chin% names(states2[states2=="island"])] %>% View()
spark_metadata[id %chin% names(states2[states2=="absent"])]
spark_metadata[id %chin% names(states2[states2=="g1"])]


# association / dissociation -------------------------------------------------------




