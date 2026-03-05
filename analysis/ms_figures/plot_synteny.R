
# Visualise AG inter-regions ----------------------------------------------

pirate_full <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"))

# msu_paths_dt <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"))
nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"))
msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv")) # core block ids


# recreate paths dt (the saved list is truncated)
pangraph_data <- yyjsonr::read_json_file("./input_data/pangraph/graph.json",
                                         int64 = "string"
)

# extract paths
# Extract relevant fields from each path
path_dt <- rbindlist(lapply(pangraph_data$paths, function(x) {
  data.table(
    path_id = x$id,
    geno_id = x$name,
    tot_len = x$tot_len,
    # convert node IDs to character to preserve precision
    nodes = list(as.character(x$nodes))
  )
}), fill = TRUE)


# for space remove Json file (not concerned with blocks atm)
rm(pangraph_data)

# Get block ids
focal_msu_blocks <- msu_mergers_dt[msu_mergers == "MSU_1", block_id]

get_region_ag_blocks <- function(node_dat, path_dat, fcore_blocks){
  
  # nodes
  msu_nodes = node_dat[block_id %chin% fcore_blocks, node_id]
  
  path_dat[, subset_nodes := lapply(nodes, function(nds) {
    hits <- match(msu_nodes, nds)
    hits <- hits[!is.na(hits)]
    if (length(hits) == 0) return(NA_character_)
    nds[seq(min(hits), max(hits))]
  })]
  
  # iterate over genomes ----
  node_sub_list <- lapply(seq_len(nrow(path_dat)), function(i) {
    geno <- path_dat$geno_id[i]
    keep_nodes <- path_dat$subset_nodes[[i]]
    if (length(keep_nodes) == 0) return(NULL)
    node_dat[geno_id == geno & node_id %chin% keep_nodes]
  })
  
  node_sub <- rbindlist(node_sub_list, use.names = TRUE, fill = TRUE)
  
  return(node_sub)
}


msu_1_dt <- get_region_ag_blocks(node_dat = nodes_dt, path_dat = path_dt, 
                                 fcore_blocks = focal_msu_blocks) 


# Ci vs parsimony ---------------------------------------------------------
phylo_info <- fread(, paste0(outdir_dat, "/phylo_info.csv"))

box_plots <- phylo_info[!is.na(par_gains),.(gene_family, consistencyindex, par_gains)]

box_plots[, ci_bin:= sprintf("%1.2f",round(consistencyindex / 0.05) * 0.05)]


box_plots[, ci_bin:= factor(ci_bin, levels = sprintf("%1.2f", seq(0,1,0.05)))]


# Split gains by ci_bin
gain_list <- split(box_plots$par_gains, box_plots$ci_bin)

# Remove empty groups
gain_list_clean <- gain_list[sapply(gain_list, length) > 0]

positions <- match(names(gain_list_clean),
                   levels(box_plots$ci_bin))



# PLOT!!!!! ---------------------------------------------------------------

# 
png(paste0(outdir_fig,"/plot_synteny.png"),
    width = 15.9, height = 18.5, units = "cm", res = 300,
    pointsize = 12, type = "cairo")

mat <- matrix(c(1,1,1,2,2,2,3,4,4,5,5,6), nrow = 4, byrow = TRUE)
layout(mat, heights = c(1.5,1,3,2),
       widths = c(1,2.5,3.5))
# layout.show(n=6)

# Core-block schematic -----------------------------------------------------
par(mar = c(2,12,2,10), xpd = NA)
plot(NULL,
     xlim = c(0, 200),          # adjust to your msu range
     ylim = c(0.5, 4.5),
     xlab = "",
     ylab = "",
     yaxt = "n",
     xaxt = "n",
     bty = "n")                 # remove box if desired

axis(2, line = -1, tick = FALSE,
     at = c(2, 3, 4),
     labels = c("Genome 3", "Genome 2", "Genome 1"),
     las = 1)

axis(3, line = -1.7, tick = FALSE,
     at = c(14,65,115,174.5),cex = 0.95,
     labels = c("core 1", "core 2", 
                expression("core"~italic("n")-1), 
                expression("core"~italic("n"))),
     las = 1)

height <- 0.25

polygon(x = c(5, 80, 80, 5), 
        y = c(2 - height, 2 - height, 2 + height, 2 + height),
          col = "grey80", border = NA)
polygon(x = c(0, 90, 90, 0), 
        y = c(3 - height, 3 - height, 3 + height, 3 + height),
        col = "grey80", border = NA)
polygon(x = c(8, 70, 70, 8), 
        y = c(4 - height, 4 - height, 4 + height, 4 + height),
        col = "grey80", border = NA)
polygon(x = c(110, 180, 180, 110), 
        y = c(4 - height, 4 - height, 4 + height, 4 + height),
        col = "grey80", border = NA)

segments(75,4,105,4,
         lwd = 4, lty = 3)

genes <- data.frame(
  genome = c(rep(4,7), rep(3,5), rep(2,4)),
  start  = c(20,26,35,120,130,140,155,
             8,35,46,55,72,
             16, 22, 38, 65),
  end    = c(26,35,60,130, 140,155,169,
             35,46,55,72, 83,
             22, 38, 65, 70),
  col    = c(rep("white",2),rep("grey40",3),rep("white",2),
             "white","grey40","white", rep("red",2),
             rep("white",2),rep("grey40",2)),
  border = c(rep("black",10), rep("firebrick4",2),
             rep("black",4)),
  lwd = c(rep(1.5,16))
)

gene_height <- 0.3

for (i in 1:nrow(genes)) {
  y_center <- genes$genome[i]
  
  rect(xleft   = genes$start[i],
       xright  = genes$end[i],
       ybottom = y_center - gene_height,
       ytop    = y_center + gene_height,
       col     = genes$col[i],
       border  = genes$border[i],
       lwd  = genes$lwd[i])
}

height <- 0.27

polygon(x = c(195, 220, 220, 195), xpd=TRUE,
        y = c(4.25 - height, 4.25 - height, 4.25 + height, 4.25 + height),
        col = "grey80", border = NA)

polygon(x = c(195, 220, 220, 195), xpd=TRUE,
        y = c(3.65 - height, 3.65 - height, 3.65 + height, 3.65 + height),
        col = "white", border = "black", lwd = 1.5)

polygon(x = c(195, 220, 220, 195), xpd=TRUE,
        y = c(3.05 - height, 3.05 - height, 3.05 + height, 3.05 + height),
        col = "grey40", border = "black", lwd = 1.5)

polygon(x = c(195, 220, 220, 195), xpd=TRUE,
        y = c(2.45 - height, 2.45 - height, 2.45 + height, 2.45 + height),
        col = "red", border = "firebrick4", lwd = 1.5)

axis(4, tick = FALSE, line = 2.5, xpd=TRUE,
     at = c(2.45, 3.05, 3.65, 4.25),
     labels = c("cloud < 0.15","shell < 0.95","soft-core < 0.99", "core"),
     las = 1)

# highlight AGs
segments(16, 1.25, 70, 1.25,
         lwd = 1.2)
segments(16, 1.25, 16, 1.35,
         lwd = 1.2)
segments(70, 1.25, 70, 1.35,
         lwd = 1.2)
text(43,1.25,"Accessory\ngenes", 
      pos = 1, cex =0.85)

# highlight AGs
segments(0,0.25, 180, 0.25,
         lwd = 2)
segments(0,0.25, 0, 0.35,
         lwd = 2)
segments(180,0.25, 180, 0.35,
         lwd = 2)
text(90,0.25,"Aggregate collinearly oriented core genes into a 'core block'", 
     pos = 1, xpd = TRUE)

boxed.labels(x = -105, y = 3, cex = 0.8,
             labels = "Synteny at\nlocal level",
             bg = "white", xpad=1.4,ypad=1.4,
             border = "black", lwd = 1.5)

boxed.labels(x = -100, y = 0.25,cex = 0.8,
             labels = "Genome-wide\nsynteny",
             bg = "white", xpad=1.4,ypad=1.4,
             border = "black", lwd = 1.5)

usr <- par("usr")
y_gap = (abs(usr[3])+usr[4])/30
text(x = usr[1] - 80, y = usr[4] + y_gap, labels = "a",
     xpd = TRUE, font = 2,cex = 1.5)
# core orientation --------------------------------------------------------
par(mar = c(0,10,5,8), xpd = NA)

n <- 35
gene_lengths <- runif(n, 0.5, 2)   # variable gene widths
x_start <- cumsum(c(0, head(gene_lengths, -1)))
x_end   <- x_start + gene_lengths

plot(0, 0,
     type = "n",
     xlim = c(0, sum(gene_lengths)),
     ylim = c(0, 1),
     xlab = "",
     main = "core-block orientation",
     ylab = "",
     xaxt = "n",
     yaxt = "n",
     bty = "n")

cols_rain = rainbow(n)

for (i in 1:n) {
  rect(x_start[i], 0.4,
       x_end[i],   0.6,
       col = cols_rain[i],
       border = "black")
}

axis(side = 2, at = 0.5, tick = FALSE,
     labels = "genome A", las = 2)

usr <- par("usr")
y_gap = (abs(usr[3])+usr[4])/30
text(x = usr[1]- 10, y = usr[4] + y_gap, labels = "b",
     xpd = TRUE, font = 2,cex = 1.5)

# synteny plot data  ------------------------------------------------------

test_dat = msu_1_dt[geno_id %chin% c("SPARK_1006_C1", "SPARK_1294_C1","SPARK_2668_C1") 
                    & start > 700000 
                    & end < 2000000, ]

# Part 1. get tree
core_gub_tree <- read.tree("./input_data/bootstrapped_gubbins/tmp8yl_0c9w/RAxML_bestTree.core_genome_aln.iteration_20")

core_gub_tree = keep.tip(core_gub_tree, unique(test_dat$geno_id))
core_gub_tree$node.label <- NULL

# Convert ape phylo → Newick string
nwk <- write.tree(core_gub_tree)

# Convert Newick string → ade4 phylog
tree_phylog <- ade4::newick2phylog(nwk)

# Part 2. get gene annotations
generate_dna_segs <- function(anno_file, synt_dat, geno_ord){
  # get anno min max
  subset_anno_minmax <- synt_dat[, .(
    start = min(start),
    end   = max(end)
  ), by = geno_id]
  
  # get annotations
  subset_anno <- subset_anno_minmax[
    anno_file,
    on = .(geno_id, start <= start, end >= end),
    nomatch = 0
  ]
  
  # set locus_tag as the name column
  colnames(subset_anno) = gsub("locus_tag", "name", colnames(subset_anno))
  
  # fix strand
  subset_anno[, strand:= fcase(strand == "+", 1,
                               default = -1)]
  
  # add plot values
  subset_anno$fill = "grey50"
  subset_anno$col = "grey50"
  subset_anno$gene_type = "headless_arrows"#"arrows"
  subset_anno$lty = 1
  subset_anno$lwd = 1
  subset_anno$pch = 8
  subset_anno$cex = 1
  
  # highlight core
  subset_anno[ag_type == "soft", 
              `:=`(fill = "white", col = "black", gene_type = "blocks")]
  subset_anno[ag_type == "shell", 
              `:=`(fill = "grey20", col = "black", gene_type = "blocks")]
  subset_anno[ag_type == "cloud", 
              `:=`(fill = "brown1", col = "darkred", gene_type = "blocks")]
  
  # order by genome order
  subset_anno[, geno_id := factor(geno_id, levels = geno_ord)]
  
  # order the data.table by that factor order
  setorder(subset_anno, geno_id)
  
  # order columns
  first_cols <- c("name","start", "end", "strand")
  
  # reassign column order
  setcolorder(subset_anno, c(intersect(first_cols, names(subset_anno)), setdiff(names(subset_anno), first_cols)))
  
  
  # Split by ST, drop ST column in each sub-table
  
  seq_list = split(subset_anno, by = "ST")
  
  # seq_list = lapply(seq_list, dna_seg)
  
  return(seq_list)
  
  
}

dna_seqs <- generate_dna_segs(anno_file = pirate_full, geno_ord = core_gub_tree$tip.label,
                              synt_dat = test_dat)


# Part 3. get sequential gene comparisons
# get colours
get_block_colours <- function(geno_ord, synt_dat){
  
  col_dat = synt_dat
  core = unique(synt_dat[count==93 & n_strains == 93, 
                         .(block_id = block_id, col = "grey70")])
  
  col_dat = col_dat[!block_id %chin% core$block_id]
  
  col_dat[, ST := factor(ST, levels = geno_ord)]
  
  setorderv(col_dat, cols = c("ST", "start"))
  
  block_ls_dt = col_dat[, .(block_ls = list(block_id),
                            n = .N), by = geno_id]
  
  tmp_cols <- data.frame(block_id = unlist(block_ls_dt[n == max(n)][1]$block_ls))
  
  tmp_cols$col = colorRampPalette(c("#420D55","darkblue","#008080",
                                    "#799000","#FF8000"))(nrow(tmp_cols))
  
  # light
  # tmp_cols$col = colorRampPalette(c("#895C9E","#5A7FCF","#4FB5B5",
  # "#B7C24D","#FFB366"))(nrow(tmp_cols))
  
  # add in missing colours
  tmp2_cols <- data.frame(block_id = col_dat[!block_id %chin% tmp_cols$block_id, block_id])
  
  if(nrow(tmp2_cols)>0){
    # light
    tmp2_cols$col = colorRampPalette(c("#FFB300", "#D94E1F", "#B30059",
                                       "#6A1B9A", "#283593"))(nrow(tmp2_cols))
    # tmp2_cols$col = colorRampPalette(c("#FFD666","#ED8C66","#E06699",
    #                                    "#A966C9","#6D78C6"))(nrow(tmp2_cols))
    
    tmp_cols <- rbind(tmp_cols, tmp2_cols)
  }
  
  tmp_cols <- rbind(tmp_cols, core)
  
  return(tmp_cols)
}

# start1, end1, start2 and end2, col, direction
sequential_comp <- function(geno_ord, synt_dat){
  
  # get colours
  block_colours <- get_block_colours(geno_ord, synt_dat)
  
  synt_dat[, rep := .N, by = c("geno_id", "block_id", "strand")]
  
  # find gene pairs
  pairs <- lapply(seq_along(geno_ord[-length(geno_ord)]), function(i) geno_ord[i:(i+1)])
  
  comp_list <- lapply(pairs, function(p){
    cln_synt_dat <- synt_dat[geno_id %chin% p & !block_id %chin% synt_dat[rep>1, block_id]]
    
    # Reshape to wide
    synt_dat_wd <- dcast(
      cln_synt_dat[, .(block_id, start, end, geno_id, strand)],
      block_id + strand  ~ geno_id,
      value.var = c("start", "end")
    )
    
    # add colour
    synt_dat_wd <- merge(synt_dat_wd, block_colours, all.x = TRUE, by ="block_id")
    
    # fix pos names
    colnames(synt_dat_wd) = gsub("strand", "direction", names(synt_dat_wd))
    
    synt_dat_wd$direction = ifelse(synt_dat_wd$direction == "+", 1, -1)
    
    # fix strand
    if(grepl(p[1], p[2])){
      colnames(synt_dat_wd) = gsub(paste0("_", p[1]), "1", 
                                   gsub(paste0("_", p[2]), "2", names(synt_dat_wd)))
    } else{
      colnames(synt_dat_wd) = gsub(paste0("_", p[2]), "2", 
                                   gsub(paste0("_", p[1]), "1", names(synt_dat_wd)))
    }
    
    first_cols <- c("start1", "end1", "start2", "end2", "col","direction")
    
    # reassign column order
    setcolorder(synt_dat_wd, c(intersect(first_cols, names(synt_dat_wd)), setdiff(names(synt_dat_wd), first_cols)))
    
    # remove nas
    # comp <- genoPlotR::comparison(synt_dat_wd)
    comp <- synt_dat_wd
    comp[!is.na(comp$direction),]
    
  })
  
  return(comp_list)
}


comparisons <- sequential_comp(geno_ord = core_gub_tree$tip.label,
                               synt_dat = test_dat)

# Annotations
# mids <- apply(dna_seqs[[1]][,c("start", "end")], 1, mean)

# # get gene names for non-core genes
# text <- dna_seqs[[1]][,c("ag_type", "gene")]
# text$gene[text$ag_type == "core"] <- ""
# text <- text$gene
# 
# annot <- annotation(x1=mids, text=text, rot=30)

# names(dna_seqs) <- gsub("-", "_", names(dna_seqs))



# plot_gene_map(dna_segs=dna_seqs,
#               comparisons=comparisons,
#               tree=tree_phylog,
#               dna_seg_line = FALSE,
#               tree_width=1.5,
#               tree_scale = TRUE,
#               annotations=annot)


# Plot synteny plot ---------------------------------------------


# gene tree
sub_tree = keep.tip(core_gub_tree, unique(test_dat$geno_id))

par(mar = c(4.5, 0.05,0.05,0))
plot(sub_tree,
     edge.width = 1,                 
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.65,
     no.margin = FALSE, main="")
# 
# usr <- par("usr")
# y_gap = (abs(usr[3])+usr[4])/40
# text(x = usr[2]+3, y = usr[4] + y_gap, labels = "c",
#      xpd = TRUE, font = 2,cex = 1.5, xpd = TRUE)


# then synteny plot
plot_synteny <- function(gene_ls, gene_comp, height = 0.7, syn_transparency = 0.6) {
  # get y-range
  n_genomes = length(gene_ls)
  
  # get x-range
  x_bp_minmax <- rbindlist(
    lapply(gene_ls, function(g) g[, .(start = min(start), end = max(end))]),
    idcol = "gene_set"
  )[, .(start = min(start), end = max(end))]
  
  par(mar = c(4.5,5,2,0.2))
  plot(
    0, 0, type = "n",
    bty= "n",
    xlim = unlist(x_bp_minmax),
    ylim = c(0.5, n_genomes + 0.5),
    xlab = "Position (bp)",
    ylab = "", xaxs = "i",
    yaxt = "n", yaxs = "i",
    xaxt = "n"
  )
  
  axis(1)
  
  # plot syntenies
  
  invisible(
    lapply(seq_len(length(gene_comp)), function(i){
      comp_tmp = gene_comp[[i]]
      
      invisible(
        lapply(seq_len(nrow(comp_tmp)), function(k) {
          with(comp_tmp[k],
               polygon(
                 x = c(start1, end1, end2, start2),
                 y = c(i + height/3, i + height/3, 
                       i +1 - height/3, i + 1 - height/3),
                 col = adjustcolor(col, alpha.f = syn_transparency),
                 border = col
               ))
        })
      )
      
      
    })
  )
  
  # Then add some genes
  invisible(
    lapply(seq_len(length(gene_ls)), function(i){
      genes_tmp = gene_ls[[i]]
      
      ST_temp = as.character(unique(gene_ls[[i]]$ST))
      
      gene_width <- mean(genes_tmp$end - genes_tmp$start)
      
      invisible(
        lapply(seq_len(nrow(genes_tmp)), function(j) {
          berryFunctions::roundedRect(genes_tmp$start[j],
                                      i - ifelse(genes_tmp$ag_type[j]=="core",
                                                 height/3,height/2),
                                      genes_tmp$end[j],
                                      i + ifelse(genes_tmp$ag_type[j]=="core",
                                                 height/3,height/2), 
                                      col = genes_tmp$fill[j],
                                      rounding = ifelse(genes_tmp$ag_type[j]=="core",
                                                        0, 0.05),
                                      border = genes_tmp$col[j])
        })
      )
      
      axis(side = 2, at = i, las = 2, col.axis = "navyblue", cex.axis = 0.85,
           labels = ST_temp, tick = FALSE)
      
      
    }
    )
  )
  # get gene names for non-core genes
  mids <- apply(gene_ls[[n_genomes]][,c("start", "end")], 1, mean)
  g_text <- gene_ls[[n_genomes]][,c("ag_type", "gene")]
  g_text$gene[g_text$ag_type == "core"] <- ""
  g_text <- g_text$gene
  
  text(mids, rep(n_genomes, n_genomes)+ height, 
       labels = g_text, pos = 4, xpd = TRUE,
       col = "grey20",
       cex = 0.55,
       srt = 45)
  
  
}
plot_synteny(gene_ls = dna_seqs,
             gene_comp = comparisons, 
             height = 0.25, 
             syn_transparency = 0.6)

# add zoom lines
segments(1823886.0, 3.5, 1833886.0, 4.5,
         lty = 2, col= "grey60", outer = TRUE)

segments(1995395.0, 3.5, 1895395.0, 4.5,
         lty = 2, col= "grey60", outer = TRUE)


# violin plot ci vs inferred gains ----------------------------------------

par(mar = c(4, 4,2,1))

with(box_plots,
     boxplot(par_gains  ~ ci_bin, 
             yaxt = "n",
             xlab = "Consistency index bin", 
             ylab = "Number of gains",
             cex= 0.75, bty = "o",
             col = "white",
             xaxt = "n",
             border = "white",
             pch = 16, boxwex = 0.5))

axis(side = 1, at = seq(0.5, 21.5, length.out = 22),
     labels = rep("",22))
axis(side = 1, at = seq(0.5, 20.5, length.out = 5),
     labels = seq(0, 1, 0.25))

abline(h=1, col ="firebrick3", xpd=FALSE)  

vioplot(gain_list_clean, 
        add = TRUE,
        ylim = c(1:21),
        at = positions,
        col = "dodgerblue",
        border = "dodgerblue4")

axis(side = 2, at = seq(0,40, 10), 
     labels = seq(0,40, 10),
     las = 2)

text(16, 1,"Number of gains = 1", cex = 0.75,
     col = "firebrick3", pos = 3)

usr <- par("usr")
y_gap = (abs(usr[3])+usr[4])/30
text(x = usr[1] - 3, y = usr[4] + y_gap, labels = "c",
     xpd = TRUE, font = 2,cex = 1.5)

# MGE proportions ---------------------------------------------------------
par(mar = c(4,5,3,1))
data <- matrix(c(4,10,5,3,12,7,13,6,15,1), nrow=2, byrow = TRUE)
colnames(data) <- c("prophage","ICE/IMEs","integrons","transposons","insertion\nsequences")
rownames(data) <- c("Non-syntenic","Discordant_phylogeny")


# Get the stacked barplot
barplot(data, 
        col= c("steelblue4","lightblue4"),
        border = "white",
        bty = "o",
        ylab = "Number of AGs",
        ylim = c(0,25),
        space=0.3, 
        yaxt = "n",
        font.axis=2, 
        xlab="")
box()

axis(side = 2, at = seq(0,25, 5), 
     labels = seq(0,25, 5),
     las = 2)

legend("topright",
       legend = c("Non-syntenic", "Phylogentic"),
       fill = c("steelblue4", "lightblue4"),
       border = "white",
       bty = "n")


usr <- par("usr")
y_gap = (abs(usr[3])+usr[4])/30
text(x = usr[1]-0.6, y = usr[4]+y_gap, labels = "d",
     xpd = TRUE, font = 2,cex = 1.5)

dev.off()























# # post hoc labels
# 
# 
# # get gene names for non-core genes
# mids <- apply(dna_seqs$gen[ag_type=="cloud"][,c("start", "end")], 1, mean)
# mids <- mids[c(3,22,37)]
# g_text <- c("tpn", "phgprot", "phgtail")
# 
# text(mids, rep(2, 3)+ height/2, 
#      labels = g_text, pos = 4, xpd = TRUE,
#      col = "grey20",
#      cex = 0.55,
#      srt = 45)


# # dotplots ----------------------------------------------------------------
# # consensus "ST23-1LV"
# focal_ST = c("ST23-1LV","ST2703","ST86","ST268","ST3586",
#              "ST37-1LV","ST11","ST258","ST258-1LV","ST405-1LV")
# 
# # run dotplots in python
# # variables
# str_i <- "ST23-1LV"
# str_j <- "ST2703"
# output_dir <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/figures"
# 
# file_dir <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph"
# py_script <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/dotplot.py"
# 
# # Construct WSL command
# cmd_dotplots_py <- paste(
#   "wsl",
#   "cd", file_dir, "&&",
#   "/home/carac/anaconda3/bin/python", py_script,
#   str_i, str_j, output_dir
# )
# 
# # Run it
# system(cmd_dotplots_py)

