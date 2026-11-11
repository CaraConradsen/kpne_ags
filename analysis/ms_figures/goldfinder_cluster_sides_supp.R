
# Supplementary figure: selection within co-occurrence clusters ----------
#   a  clusters confined to a single syntenic locus
#   b  clusters whose members sit at multiple loci
#   c  example cluster recovered by the simultaneous score
#   d  example cluster recovered by the subsequent score
#
# Panels a and b are built by goldfinder_catepillar_plot.R, sourced with
# make_fig = FALSE so that it returns syn_cnt, non_cnt and draw_panel()
# without opening a device of its own. Keeping one copy of that code means
# the supplementary figure can never drift from the analysis.
#
# heatmap() calls layout() internally and layout() cannot be nested inside
# split.screen(), so c and d are written to PNG first and rasterised back
# in, the same approach used for the previous version of this figure.

# Both examples are drawn from panel a (single syntenic locus) and are
# bolded there, so the reader can find them. Each is scored by one method
# only, so the two panels contrast the mechanisms cleanly.
#
# Cluster 72 is the ydj sugar catabolic operon discussed in the Results:
# all eight members carry fewer synonymous sites than expected and only
# tdh passes correction, which is the partial-detection pattern panel a is
# built to show. Cluster 174 is the opposite case, four fimbrial genes all
# passing correction in the excess direction.
clst_sim <- "72"       # simultaneous score only, 8 genes, 1 outlier (tdh)
clst_sub <- "174"      # subsequent score only, 4 genes, all four outliers

heat_dir <- paste0(outdir_fig, "/clusters_heatmaps")
if (!dir.exists(heat_dir)) dir.create(heat_dir, recursive = TRUE)

# Panels a and b ----------------------------------------------------------

make_fig <- FALSE
highlight_clusters <- c(clst_sim, clst_sub)   # bolded in panel a
source("./analysis/10_outlier_function/goldfinder_catepillar_plot.R")

# Panels c and d ----------------------------------------------------------
# Colours follow the rest of the manuscript: orange for reduced diversity,
# blue for excess, with the lighter shade for genes that departed in that
# direction without passing correction, and grey for genes with no test.

heat_breaks <- c(-1.5, -0.75, -0.25, 0.25, 0.75, 1.5, 2.5)
heat_cols   <- c("#F26C2A", "#F9A876", "white", "#4DB3FF", "#1A7ACC", "grey70")
out_grp_cols <- setNames(heat_cols, c(-1, -0.5, 0, 0.5, 1, 2))

gene_pvals_h <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                      select = c("gene_family", "seg_sites", "expected_S", "p_bh"))

gene_pvals_h[, sel_type := fifelse(seg_sites > expected_S, "balancing",
                                   "directional")]

gene_pvals_h[, heat_val := fcase(
  sel_type == "balancing"   & p_bh <= q_cut,  1,
  sel_type == "balancing"   & p_bh >  q_cut,  0.5,
  sel_type == "directional" & p_bh <= q_cut, -1,
  sel_type == "directional" & p_bh >  q_cut, -0.5,
  default = 2)]

sim_clus_dt <- fread(paste0(outdir_dat, "/simul_clust_fams.csv"))
sub_clus_dt <- fread(paste0(outdir_dat, "/sub_clust_fams.csv"))

pan_anno_h <- fread(paste0(outdir_dat, "/pangenome_anno.csv"),
                    select = c("geno_id", "gene_family"))

# pa_gf_gene_heatmap.csv caches presence/absence already painted with
# heat_val, but it was written before the current outlier calls and stores
# tdh as non-significant. The matrix is therefore rebuilt here from the
# current p-values, for the cluster being drawn only.

cluster_pa <- function(genes, tips) {

  d <- pan_anno_h[gene_family %chin% genes]

  m <- dcast(d, geno_id ~ gene_family, fun.aggregate = length,
             value.var = "gene_family")

  m <- m[data.table(geno_id = tips), on = .(geno_id)]      # pad absent genomes

  for (g in genes) {
    if (!g %chin% names(m)) m[, (g) := 0L]
    hv <- gene_pvals_h[gene_family == g, heat_val]
    if (!length(hv)) hv <- 2                               # no neutrality test
    v <- as.double(m[[g]]); v[is.na(v)] <- 0
    m[, (g) := fifelse(v > 0, hv, 0)]
  }

  m[, c("geno_id", genes), with = FALSE]
}

# Tree as a dendrogram ----------------------------------------------------
# Branch lengths are log-compressed so that the shallow structure within
# sequence types stays visible alongside the deep splits.

tree <- read.tree(
  "./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")

if (!is.rooted(tree)) {
  tree <- root(tree, outgroup = tree$tip.label[1], resolve.root = TRUE)
}

tree$edge.length <- log1p(tree$edge.length)

dend <- as.dendrogram(as.hclust(phytools::force.ultrametric(tree)))

# One cluster heatmap -----------------------------------------------------

plot_gene_heatmap <- function(clst_id, dend, gene_pvals, tree) {

  temp_clst_dt <- clean_coin_ags[as.character(clust_id) == as.character(clst_id)]

  # the strongest outlier in the cluster carries the branch colouring
  focal_out_ag <- temp_clst_dt[p_bh == min(p_bh, na.rm = TRUE), gene_family][1]
  focal_ag_col <- gene_pvals[gene_family == focal_out_ag, heat_val]

  present <- pan_anno_h[gene_family == focal_out_ag, geno_id]

  dend %>% set("by_labels_branches_col",
               value     = present,
               TF_values = c(out_grp_cols[[as.character(focal_ag_col)]], Inf)) -> dend_col

  # simultaneous members first, then the rest of their cluster, then the
  # subsequent members, so related genes sit together along the x axis
  sim_ags <- if (nrow(temp_clst_dt[score == "simultaneous"]) > 0) {
    temp_clst_dt[score == "simultaneous"][order(p_bh), gene_family]
  } else NULL

  sim_clus_num <- if (length(sim_ags)) {
    unique(sim_clus_dt[gene_family %chin% sim_ags, assoc_cluster])
  } else NULL

  othr_sim_ags <- if (length(sim_clus_num)) {
    setdiff(sim_clus_dt[as.character(assoc_cluster) %chin%
                          as.character(sim_clus_num), gene_family], sim_ags)
  } else NULL

  sub_ags <- if (nrow(temp_clst_dt[score != "simultaneous"]) > 0) {
    temp_clst_dt[score != "simultaneous"][order(p_bh), gene_family]
  } else NULL

  sub_clus_num <- if (length(sub_ags)) {
    unique(sub_clus_dt[gene_family %chin% sub_ags, assoc_cluster])
  } else NULL

  othr_sub_ags <- if (length(sub_clus_num)) {
    setdiff(sub_clus_dt[as.character(assoc_cluster) %chin%
                          as.character(sub_clus_num), gene_family], sub_ags)
  } else NULL

  temp_coin_genes <- unique(c(sim_ags, othr_sim_ags,
                              setdiff(sub_ags, c(sim_ags, othr_sim_ags)),
                              setdiff(othr_sub_ags, c(sim_ags, othr_sim_ags))))

  temp_coin_genes <- intersect(temp_coin_genes, pan_anno_h$gene_family)

  temp_mat <- cluster_pa(temp_coin_genes, tree$tip.label)
  temp_mat[, geno_id := factor(geno_id, levels = tree$tip.label)]
  setorder(temp_mat, geno_id)

  mat <- as.matrix(temp_mat[, -1])
  row.names(mat) <- temp_mat$geno_id

  # the panel tag in the composite carries the cluster and score, so the
  # heatmap's own title would only duplicate it
  heatmap(mat, Colv = NA, Rowv = dend_col, scale = "none", labRow = NA,
          col = heat_cols, cexCol = 1, margins = c(4, 2), breaks = heat_breaks,
          main = "")
}

# Rendered at the aspect ratio of the panel it will occupy, so rasterImage
# does not have to letterbox it and leave white bands either side.
write_heatmap <- function(clst_id) {
  f <- paste0(heat_dir, "/cluster_", clst_id, ".png")
  png(f, width = 8, height = 4.6, units = "in", res = 300, pointsize = 15,
      bg = "white")
  plot_gene_heatmap(clst_id, dend, gene_pvals_h, tree)
  dev.off()
  f
}

heat_sim <- write_heatmap(clst_sim)
heat_sub <- write_heatmap(clst_sub)

stopifnot(file.exists(heat_sim), file.exists(heat_sub))

# one raster panel, scaled to fit without distorting the aspect ratio
draw_panel_png <- function(path, tag) {
  img <- readPNG(path)
  par(mar = c(0, 0, 1.6, 0), xaxs = "i", yaxs = "i")
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
  d      <- dim(img)
  ar_img <- d[1] / d[2]
  pin    <- par("pin")
  ar_reg <- pin[2] / pin[1]
  if (ar_img > ar_reg) { h <- 1; w <- ar_reg / ar_img } else { w <- 1; h <- ar_img / ar_reg }
  rasterImage(img, (1 - w) / 2, (1 - h) / 2, (1 + w) / 2, (1 + h) / 2)
  mtext(tag, side = 3, adj = 0, font = 2, cex = 0.95, line = 0.2,
        at = grconvertX(0.01, "nfc", "user"), xpd = NA)
}

# Compose -----------------------------------------------------------------

png(paste0(outdir_fig, "/goldfinder_cluster_sides_supp.png"),
    width = 32, height = 31, units = "cm", res = 300, pointsize = 12,
    type = "cairo")

# a and b are sized to their own row counts rather than a shared one, so
# panel b no longer carries 18 empty rows; the x scale is still shared, so
# bar lengths remain comparable between them.
split.screen(rbind(c(0,   0.5, 0.24, 1),
                   c(0.5, 1,   0.24, 1),
                   c(0,   0.5, 0,    0.24),
                   c(0.5, 1,   0,    0.24)))

screen(1)
draw_panel(syn_cnt, sprintf("a  Single syntenic locus (n = %d)", nrow(syn_cnt)))
legend(x = xmin - pad, y = nmax + 2 * grp_gap + 6.4, horiz = TRUE, bty = "n",
       cex = 0.66, fill = c(fire, "grey78", sky), border = "white",
       legend = c("Reduced diversity", "Not significant", "Excess diversity"),
       xpd = NA)
close.screen(1)

screen(2)
draw_panel(non_cnt, sprintf("b  Multiple loci (n = %d)", nrow(non_cnt)))
close.screen(2)

screen(3)
draw_panel_png(heat_sim, paste0("c  Simultaneous, cluster ", clst_sim))
close.screen(3)

screen(4)
draw_panel_png(heat_sub, paste0("d  Subsequent, cluster ", clst_sub))
close.screen(4)

close.screen(all.screens = TRUE)
dev.off()

cat("wrote", paste0(outdir_fig, "/goldfinder_cluster_sides_supp.png"), "\n")
