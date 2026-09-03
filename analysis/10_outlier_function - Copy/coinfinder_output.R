
ass_ags <- fread("./input_data/coinfinder/association_output_pairs.tsv")
colnames(ass_ags) <- tolower(colnames(ass_ags))

#tidy gene names to gene family
ass_ags[, source:= tstrsplit(source, "__", fill = TRUE, keep = 2)]
ass_ags[, target:= tstrsplit(target, "__", fill = TRUE, keep = 2)]

# odds ratio
ass_ags[, OE_ratio := `successes (oa(ij))` / pmax(`expected (ea(ij))`, 1)] 

# Or — much more powerful and what I'd actually recommend — Benjamini-Hochberg
ass_ags[, p_bh := p.adjust(p, method = "BH")]
ass_sig_bh <- ass_ags[p_bh < 0.05]
nrow(ass_sig_bh)

sweep <- rbindlist(lapply(1:50, function(thr) {
  sub <- ass_sig_bh[OE_ratio >= thr]
  if (nrow(sub) == 0) return(NULL)
  g_thr <- simplify(graph_from_data_frame(
    sub[, .(from = source, to = target)], directed = FALSE))
  cc_thr <- components(g_thr)
  data.table(
    OE_thr      = thr,
    n_edges     = ecount(g_thr),
    n_nodes     = vcount(g_thr),
    n_components = cc_thr$no,
    largest_comp = max(cc_thr$csize),
    frac_in_giant = max(cc_thr$csize) / vcount(g_thr)
  )
}))
print(sweep)

set_oe_ratio = sweep[frac_in_giant<0.7,.(min(OE_thr))][[1]]

# build  graph from filtered edges  ---------------------------------------

g <- graph_from_data_frame(
  ass_sig_bh[OE_ratio >= set_oe_ratio, .(from = source, to = target)],
  directed = FALSE
)

# Simplify: collapse any duplicate edges and remove self-loops
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

cat("Graph: ", vcount(g), "nodes,", ecount(g), "edges\n")

# find connected components 
cc <- components(g)

# Per-family cluster assignment
fam_clusters <- data.table(
  gene_family  = names(cc$membership),
  cluster = as.integer(cc$membership)
)

new_name <- paste0("ass_clust_bh_oe", set_oe_ratio)
setnames(fam_clusters, old = "cluster", new = new_name)

fwrite(fam_clusters, paste0(outdir_dat, "/", new_name, ".csv"))
