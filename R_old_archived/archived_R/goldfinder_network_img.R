clean_coin_ags <- fread(paste0(outdir_dat,"/clean_coin_ags.csv"),
                        select = c("gene_family","p_bh", "w_theta_sel", "set","perc_in_set1",
                                   "n_genes_assoc_cluster", "score", "clust_id"))
# Id the set 2 AG
coin_ag_set2 <- clean_coin_ags[set == 2, gene_family]

clean_coin_ags[, set:=NULL]

# do we use only AGs in set 1?
# clean_coin_ags <- clean_coin_ags[perc_in_set1==100]

# fix overlappers
clean_coin_ags[, score := paste(sort(unique(score)), collapse = ":"), 
   by = .(gene_family, clust_id)]

# Collapse to one row per gene_family/clust_id
clean_coin_ags <- unique(clean_coin_ags, by = c("gene_family", "clust_id"))


# read in pvals and network info ----------------------------------------------

# fdr corrected p-values
# Pool all p-values, one global BH correction
gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

sim_edges <- fread("./input_data/goldfinder/simul_output/cytoscape_input.csv")

sim_edges$score = "simultaneous"

sub_edges <- fread("./input_data/goldfinder/sub_output/cytoscape_input.csv")

sub_edges$score = "subsequent"

edges <- rbind(sim_edges, sub_edges)

# fix edge names
edges[, Node1:= tstrsplit(Node1, "__", fill = TRUE, keep = 2)]
edges[, Node2:= tstrsplit(Node2, "__", fill = TRUE, keep = 2)]

# keep only focal genes
edges <- edges[Node1 %chin% clean_coin_ags$gene_family & 
                 Node2 %chin% clean_coin_ags$gene_family]

edges <- edges[!(Node1 == "g005811" & Node2 == "g002445")]

# remove duplicated edges
edges[, c("n1", "n2") := .(pmin(Node1, Node2), pmax(Node1, Node2))]
edges <- unique(edges, by = c("n1", "n2"))
edges[, c("n1", "n2") := NULL]

# fwrite(edges, paste0(outdir_dat,"/out_edges.csv"))



# set up dt ---------------------------------------------------------------
# colour

clean_coin_ags[, pt_colr:= fcase(w_theta_sel == "balancing" & p_bh <= 0.05, "#1A7ACC",
                                 w_theta_sel == "balancing" & p_bh > 0.05,"#4DB3FF",
                                 w_theta_sel == "directional" & p_bh <= 0.05,"#F26C2A",
                                 default = "#F9A876")]

# add pch
clean_coin_ags[, pchy := 16]

if(length(coin_ag_set2) > 0 ){
  clean_coin_ags[, bg_colr:= pt_colr]
  clean_coin_ags[gene_family %in% coin_ag_set2, bg_colr:= "white"]
  
  # add pch
  clean_coin_ags[gene_family %in% coin_ag_set2, pchy := 21]
}


cex_min <- 0.8
cex_max <- 1.8

n <- clean_coin_ags$n_genes_assoc_cluster
clean_coin_ags[, cex_val := cex_min + (sqrt(n) - sqrt(min(n))) / 
             (sqrt(max(n)) - sqrt(min(n))) * 
             (cex_max - cex_min)]


# put in groups
# x_grp lookup table
border_clusts <- c(54L, 60L, 120L)

xgrp_clust <- data.table(
  clust_id = c(
    # x_grp 1
    64L, 81L, 135L, 1L, 
    # x_grp 2
    142L,
    # x_grp 3
    25L, 27L,   118L, 113L,125L, 12L,
    # x_grp 4
    5L,  19L, 35L, 40L, 112L, 99L, 3L,
    # x_grp 5
    45L, 9L, 15L, 39L, 126L,
    # x_grp 6
    115L, 21L
  ),
  x_grp = c(
    rep(1L, 4),
    rep(2L, 1),
    rep(3L, 6),
    rep(4L, 7),
    rep(5L, 5),
    rep(6L, 2)
  )
)

clean_coin_ags[, x_grp := fcase(
  clust_id %in% border_clusts & score == "simultaneous", 6L,
  clust_id %in% border_clusts & score == "simultaneous:subsequent", 7L,
  score == "simultaneous:subsequent" & !clust_id %in% border_clusts & clust_id != 73, 8L,
  score == "simultaneous:subsequent" & clust_id == 73, 9L,
  score == "subsequent" & clust_id != 143, 10L,
  score == "subsequent" & clust_id == 143, 11L,
  default = xgrp_clust[.SD, on = "clust_id", x_grp]
)]

max_groups = 11
 
# plot --------------------------------------------------------------------
# global values
max_y = ceiling(max(-log10(clean_coin_ags$p_bh)))
x_positions <- 1:11
ylim_val <- c(0, max(-log10(clean_coin_ags$p_bh)) * 1.05)

# Broad category bracket positions
# simultaneous:  x_grp 1-6
# sim:subsequent: x_grp 7-8  (border 7, pure 8)
# subsequent:    x_grp 9
broad_bounds <- list(
  simultaneous            = c(1, 6),
  "simultaneous:subsequent" = c(7, 9),
  subsequent              = c(10, 11)
)


#beeswarm

par(mar = c(2, 4, 2, 10), xpd = FALSE)

plot(NA, xlim = c(0.3, max(x_positions)+0.7), ylim = c(0, max_y),
     xaxt = "n", yaxt = "n",yaxs = "i",
     xlab = "", bty="n", ylab = "")

abline(h = -log10(0.05), lty = 2, col = "grey60", lwd = 0.8)

axis(2, at = 0:max_y, labels = 0:max_y, las = 1)
mtext(expression(-log[10](italic(p))), side = 2, line = 2.5)
mtext("Scoring method", side = 1, line = 2.5)

mtext("a", side = 3, font = 2, cex = 1.1, line = 0.5, adj = -0.05)

# Add beeswarm for each x_grp, capturing xy output
node_coords <- vector("list", max_groups)

for (grp in 1:max_groups) {
  sub <- clean_coin_ags[x_grp == grp]
  if (nrow(sub) == 0) next
  
  bs <- beeswarm(
    -log10(sub$p_bh),
    at        = grp,
    add       = TRUE,
    pwcol     = sub$pt_colr,
    pwpch     = sub$pchy,
    pwbg      = sub$bg_colr,
    cex       = max(sub$cex_val),
    pwcex     = sub$cex_val,
    method    = "swarm",
    spacing = 1.1, 
    do.plot = FALSE
  )
  
  # beeswarm() returns a data.frame with x, y (and original data)
  # attach gene_family so we can match nodes later
  bs$gene_family <- sub$gene_family
  bs$clust_id <- sub$clust_id
  node_coords[[grp]] <- bs
}

# Combine all node positions
nodes <- rbindlist(lapply(node_coords, as.data.table), fill = TRUE)
# nodes now has: x, y, gene_family

# Bracket annotations for broad categories
# -----------------------------------------------------------------------
par(xpd = NA)
brk_y  <- ylim_val[1] - 0.15 * 0.12
tick_y <- ylim_val[1] + 1 * 0.08

# simultaneous
segments(x0 = 0.7,  x1 = 6.3,  y0 = brk_y, y1 = brk_y, lwd = 1.5)
segments(x0 = c(0.7, 6.3), x1 = c(0.7, 6.3), y0 = brk_y, y1 = tick_y, lwd = 1.5)
text(x = 3.5, y = brk_y - 2 *0.06,
     "Simultaneous", cex = 0.8, col = "grey40")

# simultaneous:subsequent
segments(x0 = 6.7,  x1 = 9.3,  y0 = brk_y, y1 = brk_y, lwd = 1.5)
segments(x0 = c(6.7, 9.3), x1 = c(6.7, 9.3), y0 = brk_y, y1 = tick_y, lwd = 1.5)
text(x = 8, y = brk_y - 2 *0.06,
     "Both", cex = 0.8, col = "grey40")

# subsequent
segments(x0 = 9.7, x1 = 11.3, y0 = brk_y, y1 = brk_y, lwd = 1.5)
segments(x0 = c(9.7, 11.3), x1 = c(9.7, 11.3), y0 = brk_y, y1 = tick_y, lwd = 1.5)
text(x = 10.5, y = brk_y - 2 *0.06,
     "Subsequent", cex = 0.8, col = "grey40")

# -----------------------------------------------------------------------
# Size legend (top-right outer margin)
# -----------------------------------------------------------------------
n_vals   <- round(sort(c(range(clean_coin_ags$n_genes_assoc_cluster), 
                   mean(clean_coin_ags$n_genes_assoc_cluster))), digits = 0)
cex_vals <- sort(c(range(clean_coin_ags$cex_val), 
                   mean(clean_coin_ags$cex_val)))

leg1 <- legend(
  x = par("usr")[2] + 0.05,
  y = par("usr")[4],
  legend    = n_vals,
  pt.cex    = cex_vals * max(clean_coin_ags$cex_val) * par("cex"),
  pch       = 16,
  col       = "grey60",
  bty       = "n",
  title     = "Cluster size",
  title.adj = 0,
  y.intersp = 1.1,
  x.intersp = 1.5,
  cex = 0.75,
  xjust     = 0,
  adj       = c(0, 0.5),
  text.width = strwidth("57") * 3
)

leg2 <- legend(
  x = par("usr")[2] + 0.2,
  y = leg1$rect$top - leg1$rect$h - strheight("M") * 0.8,
  legend    = c("Set 1", "Set 2"),
  pch       = c(16, 21),
  col       = "grey60",
  pt.bg     = "white",
  pt.cex    = 1.2 * par("cex"),
  bty       = "n",
  title     = "Set",
  title.adj = 0,
  y.intersp = 1.1,
  x.intersp = 1.5,
  cex = 0.75,
  xjust     = 0,
  adj       = c(0, 0.5),
  text.width = strwidth("57") * 3
)

legend(
  x = par("usr")[2] + 0.2,
  y = leg2$rect$top - leg2$rect$h - strheight("M") * 0.8,
  legend    = c("Balancing, p \u2264 0.05",
                "Balancing, p > 0.05",
                "Directional, p \u2264 0.05",
                "Directional, p > 0.05"),
  pch       = 16,
  col       = c("#1A7ACC", "#4DB3FF", "#F26C2A", "#F9A876"),
  pt.cex    = 1.2 * par("cex"),
  bty       = "n",
  title     = "Selection",
  title.adj = 0,
  y.intersp = 1.1,
  x.intersp = 1.5,
  cex = 0.75,
  xjust     = 0,
  adj       = c(0, 0.5),
  text.width = strwidth("Directional, p \u2264 0.05") * 1.2
)

par(xpd = FALSE)


# Save node coordinates for igraph

g <- graph_from_data_frame(
  d        = edges[, .(Node1, Node2)],
  directed = FALSE,
  vertices = nodes[, .(name = gene_family, x, y)]
)

# Fixed layout from beeswarm coordinates
layout_fixed <- as.matrix(nodes[, .(x, y)])

# Plot network on top of existing beeswarm plot
par(new = TRUE)   # overlay on same plot device
plot(g,
     layout       = layout_fixed,
     rescale      = FALSE,
     vertex.size  = 0,
     vertex.shape = "none",
     vertex.color = "white",
     vertex.label = NA,
     edge.color   = adjustcolor("grey40", alpha.f = 0.5),
     edge.width   = 0.8,
     add          = TRUE)

# add points
with(nodes,
     points(x, y, 
            pch = pch, 
            col = col, 
            cex = cex))

# add labels
for (clst in unique(clean_coin_ags$clust_id)) {
  temp <- nodes[clust_id == clst]
  
  text(x      = max(temp$x), 
       y      = max(temp$y) + (temp[temp$x==max(temp$x)]$cex)/100, 
       labels = clst,
       cex    = 0.5,
       pos    = ifelse(clst %in% c(15,45,39,86,20), 4, 
                       ifelse(clst %in%  c(3,115),2, 3)),
       col    = "grey30")
}
