# Visualise AG inter-regions ----------------------------------------------

# MSU consensus order
most_common_ord_msu

# let's look at the first two blocks
msu_pair = unlist(tstrsplit(most_common_ord_msu, ",", fill = TRUE, keep = 18))

# Get block ids
focal_msu_blocks <- msu_mergers_dt[msu_mergers %chin% msu_pair, block_id]

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



# synteny plot data ------
test_dat = msu_1_dt[ST %chin% c("ST252","ST22","ST2029-1LV",
                                "ST307", "ST392","ST2447-1LV",      
                                "ST230","ST629","ST3068","ST301",     
                                "ST485", "ST101","ST13"  ), ]
#"ST423","ST3586","ST17","ST45","ST20","ST29","ST219","ST1428",
#"ST512","ST661","ST3345","ST187-1LV","ST37","ST307-1LV","ST1198",

# Part 1. get tree
ST_gene_tree = keep.tip(core_gub_tree, unique(test_dat$ST))
ST_gene_tree$node.label <- NULL

# Convert ape phylo → Newick string
nwk <- write.tree(ST_gene_tree)

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
  subset_anno[, ST := factor(ST, levels = geno_ord)]
  
  # order the data.table by that factor order
  setorder(subset_anno, ST)
  
  # order columns
  first_cols <- c("name","start", "end", "strand")
  
  # reassign column order
  setcolorder(subset_anno, c(intersect(first_cols, names(subset_anno)), setdiff(names(subset_anno), first_cols)))
  
  
  # Split by ST, drop ST column in each sub-table
  
  seq_list = split(subset_anno, by = "ST")
  
  # seq_list = lapply(seq_list, dna_seg)
  
  return(seq_list)
  
  
}

dna_seqs <- generate_dna_segs(anno_file = pangraph_anno, geno_ord = ST_gene_tree$tip.label,
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
  
  synt_dat[, rep := .N, by = c("ST", "block_id", "strand")]
  
  # find gene pairs
  pairs <- lapply(seq_along(geno_ord[-length(geno_ord)]), function(i) geno_ord[i:(i+1)])
  
  comp_list <- lapply(pairs, function(p){
    cln_synt_dat <- synt_dat[ST %chin% p & !block_id %chin% synt_dat[rep>1, block_id]]
    
    # Reshape to wide
    synt_dat_wd <- dcast(
      cln_synt_dat[, .(block_id, start, end, ST, strand)],
      block_id + strand  ~ ST,
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


comparisons <- sequential_comp(geno_ord = ST_gene_tree$tip.label,
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


# Try plotting with own ggenes ---------------------------------------------


# roundedRect: draw a rounded rectangle on base R plot
# xleft= 0.5; ybottom = y0 - height/2; xright = 1.8; ytop = y0 + height/2; r = 0.15; rel = TRUE

# gene tree
sub_tree = keep.tip(core_gub_tree, unique(test_dat$ST))


layout(matrix(c(1,2), nrow=1), widths=c(1,6))  # left = tree, right = other plot
par(mar = c(4.5, 0.05,0.05,0))
plot(sub_tree,
     edge.width = 1,                 
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.65,
     no.margin = FALSE, main="")


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

# post hoc labels


# get gene names for non-core genes
mids <- apply(dna_seqs$ST22[ag_type=="cloud"][,c("start", "end")], 1, mean)
mids <- mids[c(3,22,37)]
g_text <- c("tpn", "phgprot", "phgtail")

text(mids, rep(2, 3)+ height/2, 
     labels = g_text, pos = 4, xpd = TRUE,
     col = "grey20",
     cex = 0.55,
     srt = 45)


# dotplots ----------------------------------------------------------------
# consensus "ST23-1LV"
focal_ST = c("ST23-1LV","ST2703","ST86","ST268","ST3586",
             "ST37-1LV","ST11","ST258","ST258-1LV","ST405-1LV")

# run dotplots in python
# variables
str_i <- "ST23-1LV"
str_j <- "ST2703"
output_dir <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/figures"

file_dir <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph"
py_script <- "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/dotplot.py"

# Construct WSL command
cmd_dotplots_py <- paste(
  "wsl",
  "cd", file_dir, "&&",
  "/home/carac/anaconda3/bin/python", py_script,
  str_i, str_j, output_dir
)

# Run it
system(cmd_dotplots_py)

