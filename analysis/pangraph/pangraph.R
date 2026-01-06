# test pangraph


# Import PIRATE pangenome information -------------------------------------

# pangenome info
data_dir = "./input_data/PIRATE_485_lng_rds_out/"

# map loci
gene_families = fread(paste0(data_dir, "/PIRATE.gene_families.ordered.tsv"))
colnames = colnames(gene_families)[c(1:2, 23:ncol(gene_families))]
gene_families = gene_families[,..colnames]

gene_families = melt(gene_families, id.vars = c("allele_name", "gene_family"),
                     variable.name = "geno_id", value.name = "locus_tag")

setDT(gene_families)

gene_families <- gene_families[locus_tag!=""]


# Retain only correct oriented chromosomes --------------------------------
pan_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))
correct_start <- pan_anno[gene=="dnaA" & strand == "+" & 
                            grepl("_1$", seqnames) & start == 70, 
                          .(geno_id, ST)]

setorderv(correct_start, cols = "geno_id")

# take 2-5 genomes per ST (30 STs)
correct_start[,n := .N, ST]
correct_start_unique <- correct_start[n!=1, .SD[1:5], by = ST][!is.na(geno_id)]# subset taking first 5 STs

geno_list = correct_start_unique[,geno_id]# 115 genomes

# get fasta files

lng_rd_fasta_files = list.files("./input_data/kpne_485_fasta/",
                                full.names = TRUE)# take first 10 for now

pangraph_dir =  "./input_data/pangraph/"

chr_contigs <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
                     select = c("seqnames", "asmbly_type")))

# REMOVE PLASMIDS
chr_contigs <- chr_contigs[asmbly_type!="plasmid"][,asmbly_type:=NULL]

# get minimum gene length for pangenome
min_gene_length <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
                            select = c("geno_id", "start", "end", "strand", "gene_family"))
min_gene_length <- min_gene_length[geno_id %chin% geno_list] # focal genomes
min_gene_length <- min_gene_length[grepl("g", gene_family)]# discard small proteins not used by PIRATE
colnames(min_gene_length)[1]="seqnames"
min_gene_length <- dt2gr(min_gene_length)
min_gene_length = min(width(min_gene_length))

# Create a fasta of concatenated CDS using in PIRATE ----------------------

# cl <- makePSOCKcluster(num_cores-6)
# registerDoParallel(cl)

cds_fasta_gff_dt <- foreach(genome = geno_list,
                            .packages = c("Biostrings", "rtracklayer",
                                          "GenomicRanges", "gUtils")) %do% {
  # genome = geno_list[26]
  focal_genome = genome
  
  # import fasta
  
  fasta_file = lng_rd_fasta_files[grepl(focal_genome, lng_rd_fasta_files)]
  
  coord_fasta = Biostrings::readDNAStringSet(fasta_file)
  
  # Get PIRATE cords
  coord_fasta = coord_fasta[names(coord_fasta) %in% chr_contigs$seqnames]
  
  #--- Load PIRATE GFF cord info file ---
  
  cords = fread(paste0(data_dir, "/co-ords/", focal_genome,".co-ords.tab"))
  
  colnames(cords)[1] = "locus_tag"
  
  cords = merge(cords, gene_families[geno_id == focal_genome,
                                     .(gene_family, locus_tag, allele_name)],
                all.x = TRUE, by = "locus_tag")
  
  # Convert strand to "+" / "-"
  cords$Strand <- ifelse(cords$Strand == "Forward", "+", "-")
  
  #convert to lower case
  colnames(cords) = tolower(colnames(cords))
  
  # set contig as seqname
  colnames(cords)[which(names(cords)=="contig")] = "seqnames"
  
  # remove missing loci (short proteins)
  cords <- cords[!is.na(gene_family)]
  
  # remove plasmids
  cords <- cords[seqnames %chin% chr_contigs$seqnames]
  
  # # set order
  # setorderv(cords, cols = 'start')
  # cords[gene %chin% c("dnaA", "dnaN", "recF", "gyrB", "gyrA")]
  # 
  # Create grange file
  gr_cords <- dt2gr(cords)
  
  # subset cols
  mcols(gr_cords) <- mcols(gr_cords)[, "locus_tag", drop = FALSE]
  
  # determinenew cds positions
  cum_lengths <- cumsum(width(gr_cords))
  
  adj_gff_dt <- data.frame(start = as.integer(c(1, cum_lengths[1:length(cum_lengths)-1]+1)),
                           end = as.integer(cum_lengths),
                           strand = as.character(strand(gr_cords)),
                           locus_tag = gr_cords$locus_tag)
  
  adj_gff_dt$geno_id = focal_genome
  
  # Merge overlapping CDS on same strand
  gr_cords <- reduce(gr_cords, ignore.strand = TRUE)
  
  # Extract sequences
  cds_seqs <- getSeq(coord_fasta, gr_cords)
  
  concatenated_cds <- paste(as.character(cds_seqs), collapse = "")
  
  concatenated_cds_dna <- DNAString(concatenated_cds)
  
  # Create a DNAStringSet with one "sequence" per genome
  final_fasta <- DNAStringSet(concatenated_cds_dna)
  names(final_fasta) <- focal_genome  # genome identifier
  
  # Write CDS FASTA
  writeXStringSet(final_fasta, paste0(pangraph_dir,focal_genome, ".fasta"))
  
  cat("Processed", focal_genome, "\n")
  
  setDT(adj_gff_dt)
  
}

cds_fasta_gff_dt <- rbindlist(cds_fasta_gff_dt)
# stopCluster(cl)

# fwrite(cds_fasta_gff_dt, paste0(outdir_dat, "/cds_fasta_gff_dt.csv"))

# Assign Core vs accessory ------------------------------------------------
# cds_fasta_gff_dt <- fread(paste0(outdir_dat, "/cds_fasta_gff_dt.csv"))

cds_fasta_gff_dt <- merge(cds_fasta_gff_dt, gene_families[,!"allele_name"],
                          all.x = TRUE, by = c("geno_id", "locus_tag"))

n_genomes = length(unique(cds_fasta_gff_dt[, geno_id]))

cds_fasta_gff_dt[, n := .N, gene_family]

cds_fasta_gff_dt[, ag_type := fcase(n < 0.99 * n_genomes & n >=  0.95 * n_genomes, "soft",
                                      n < 0.95 * n_genomes & n >=  0.15 * n_genomes, "shell",
                                      n < 0.15 * n_genomes, "cloud",
                                      default = "core")]

cds_fasta_gff_dt[, n := NULL]

# Run Noll et al.'s Pangraph ----------------------------------------------

# 1. ### Generate pangraph

pangraph_path = "/home/carac/anaconda3/bin/pangraph" # pangraph program dir 
threads = 18 # Number of threads
min_length = round(min_gene_length, digits = -1) # min block size in bp
circular_flag = "--circular"   # or "" if not circular
input_dir = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph"# input folder dir
output_file = file.path(input_dir, paste0("graph.json"))# output file dir
kern_sens = "-k minimap2 -s 20"# alignment sensitivity
beta_diversity = "-b 5"# energy cost for diversity alignment

# pangraph command string 
cmd_pangraph <- sprintf(
  "wsl %s build -j %d -l %d %s %s %s %s/*.fasta -o %s",
  pangraph_path,
  threads,
  min_length,
  circular_flag,
  kern_sens,
  beta_diversity,
  input_dir,
  output_file
)

# check command
cat("Running command:\n", cmd_pangraph, "\n\n")


start <- Sys.time()# Start timer

# run pangraph 
system(cmd_pangraph)

end <- Sys.time()# End timer
print(end - start) # Print runtime
# 10 seqs Time difference of 0.58 mins
# 25 seqs Time difference of 2.177748 mins
# 93 seqs Time difference of 11.95702 mins
# 115 seqs Time difference of 13.96003 mins

# 2. ### Generate clonalframe for trees using pangraph and gubbins

# get core
pangraph_input = file.path(input_dir, "graph.json")
set_ref_strain = "SPARK_1023_C1"
core_output_file = file.path(input_dir, "core_genome_aln.fa")# output file dir

# pangraph core command string 
cmd_core <- sprintf(
  "wsl %s export core-genome %s --guide-strain %s -o %s",
  pangraph_path,
  pangraph_input,
  set_ref_strain,
  core_output_file
)

# check command
cat("Running command:\n", cmd_core, "\n\n")

# run pangraph 
system(cmd_core)


# 3. ### get clean clonal frame with gubbins (in future maybe try STs?)

conda_path = "/home/carac/anaconda3/bin/conda" # path to conda inside WSL
conda_env = "gubbins_env" # conda environment name
threads = 12 # number of threads
tree_method = "fasttree" # tree builder to use
# input alignment file (from pangraph output)
input_alignment = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/core_genome_aln.fa"
# output prefix 
output_prefix = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/gub_graph"

# build gubbins command string 
cmd_gubbins <- sprintf(
  "wsl %s run -n %s run_gubbins.py --threads %d -t %s --prefix %s %s",
  conda_path,
  conda_env,
  threads,
  tree_method,
  output_prefix,
  input_alignment
)

# ---- Print command for verification ----
cat("Running command:\n", cmd_gubbins, "\n\n")

start <- Sys.time()# Start timer

# run pangraph 
system(cmd_gubbins)

end <- Sys.time()# End timer
print(end - start) # Print runtime
# 10 seqs Time difference of 1.921924 mins
# 25 seqs Time difference of 2.211768 mins
# 93 seqs Time difference of 25.72868 mins
# 115 seqs Time difference of 25.72868 mins

# Import json -------------------------------------------------------------
# # Read from a local file

pangraph_data <- yyjsonr::read_json_file("./input_data/pangraph/graph.json",
                                         int64 = "string"
)

# Extract nodes and align to gene_families --------------------------------
nodes_dt <- rbindlist(lapply(pangraph_data$nodes, function(x) {
  data.table(
    node_id = as.character(x$id),
    block_id = as.character(x$block_id),
    path_id = x$path_id,
    strand = x$strand,
    start = x$position[1],
    end = x$position[2]
  )
}), fill = TRUE)

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

nodes_dt <- merge(nodes_dt, path_dt[,.(path_id,geno_id, tot_len)] ,
                  all.x = TRUE, by ="path_id")

nodes_dt <- merge(nodes_dt, unique(pan_anno[,.(geno_id, ST)]),
                   all.x = TRUE, by ="geno_id")

# fwrite(nodes_dt, paste0(outdir_dat, "/nodes_dt.csv"))


# Plot core
nodes_dt[, count := .N, by = block_id]
count_strains = unique(nodes_dt[, .(block_id, geno_id)]) 
count_strains = count_strains[,.(n_strains = .N), by = block_id]
nodes_dt <- merge(nodes_dt, count_strains,
                  all.x = TRUE, by = "block_id")

  
# core blocks
max_genomes = length(unique(nodes_dt$geno_id))
core_dt <- nodes_dt[count == n_strains & n_strains == max_genomes,]

# Plot Core MSU order ---------------------------------------------------------------

# Left panel: Phylogenetic tree 
core_gub_tree <- read.tree("./input_data/pangraph/gub_graph.node_labelled.final_tree.tre") 

ST_labels <- unique(pan_anno[geno_id %chin% core_gub_tree$tip.label,
                            .(geno_id, ST)])

# re-label with ST
ST_labels <- ST_labels[match(core_gub_tree$tip.label, geno_id)]  # reorder

core_gub_tree$tip.label <- ST_labels$ST


# run msu in python

file_dir = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph"
py_script = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/find_msu.py"

# Construct the WSL command
cmd_msu_py <- paste(
  "wsl",                # invoke WSL
  "cd", file_dir, "&&", # change directory in WSL
  "/home/carac/anaconda3/bin/python", py_script   # run Python script
)

# Run it
system(cmd_msu_py)

# Import data
MSU_mergers <- yyjsonr::read_json_file("./input_data/pangraph/MSU_mergers.json",
                                       int64 = "string")
MSU_paths <- yyjsonr::read_json_file("./input_data/pangraph/MSU_paths.json",
                                     int64 = "string")
MSU_len <- yyjsonr::read_json_file("./input_data/pangraph/MSU_len.json",
                                   int64 = "string")
                                                                                                               
# Extract MSU data --------------------------------
msu_mergers_dt <- rbindlist(lapply(MSU_mergers, function(x) {
  data.table(
    msu_mergers = x
  )
}), fill = TRUE)
msu_mergers_dt$block_id = names(MSU_mergers)

# fwrite(msu_mergers_dt, paste0(outdir_dat, "msu_mergers_dt.csv"))

msu_paths_dt <- rbindlist(lapply(MSU_paths, function(x) {
  data.table(
    msu_paths = x
  )
}), fill = TRUE)
msu_paths_dt$geno_id = names(MSU_paths)

msu_len_dt <- rbindlist(lapply(MSU_len, function(x) {
  data.table(
    msu_length = x
  )
}), fill = TRUE)


path_dt <- merge(path_dt, ST_labels, all.x=TRUE, by ="geno_id")

# Get cds_gff
pangraph_anno <- merge(cds_fasta_gff_dt, 
                       pan_anno[,.(geno_id, locus_tag, gene, product, COG_funct_cat,KEGG, ST)],
                       all.x = TRUE, by =c("geno_id", "locus_tag"))


msu_len_dt$msu_mergers = names(MSU_len)

num_msu_blocks = length(names(MSU_len))

msu_paths_dt[, paste0("order_",1:num_msu_blocks):= tstrsplit(msu_paths, "\\]_\\[", 
                                                                     keep = 1:num_msu_blocks)]
msu_paths_dt <- melt(msu_paths_dt[,-1], id.vars = "geno_id",
                     variable.name = "order")

msu_paths_dt$value = gsub("\\]|\\[","", msu_paths_dt$value)
msu_paths_dt$order = as.integer(gsub("order_","", msu_paths_dt$order))
msu_paths_dt[,c("msu_mergers", "strand"):= tstrsplit(value, "\\|", keep =1:2)]

msu_paths_dt[,value:=NULL]

msu_paths_dt[, ord_msu := paste(unlist(msu_mergers), collapse=","), by = geno_id]

# fwrite(msu_paths_dt, paste0(outdir_dat, "/msu_paths_dt.csv"))
# fwrite(path_dt, paste0(outdir_dat, "/path_dt.csv"))
# fwrite(pangraph_anno, paste0(outdir_dat, "/pangraph_anno.csv"))


most_common_ord_msu <- msu_paths_dt[, .N, by = ord_msu][which.max(N), ord_msu]

msu_cols = data.frame(msu_mergers = unlist(strsplit(most_common_ord_msu, ",")),
                      cols = colorRampPalette(c("#420D55","darkblue","#008080",
                                                "#799000","#FF8000"))(num_msu_blocks))
msu_plot <- merge(msu_paths_dt, msu_cols, all.x = TRUE, by = "msu_mergers")

# change to ST and set order
msu_plot <- merge(msu_plot, ST_labels, all.x = TRUE, by = "geno_id")

# Create a factor with levels matching the tree tip order
msu_plot[, ST := factor(ST, levels = core_gub_tree$tip.label)]

# Now use setorderv — it will respect the factor level order
setorderv(msu_plot, c("ST", "order"))

# Setup
yvals <- 1:length(unique(msu_plot$geno_id))      # 24 unique Y values
nblocks <- nrow(msu_cols)       # 5 squares per row
# Define square size
sq_size <- 1       # both width and height = 1

# consensus "ST23-1LV"
focal_ST = c("ST23-1LV","ST2703","ST86","ST268","ST3586",
             "ST37-1LV","ST11","ST258","ST258-1LV","ST405-1LV")

# plot Core units + phylogenetic tree
layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE), 
       widths=c(1,2), heights = c(2,4))  # left = tree, right = other plot
plot.new()

msu_plot_lengths <- merge(msu_plot[ST == "ST23-1LV", 
                                   .(msu_mergers,order,cols)], 
                          msu_len_dt, all.x = TRUE, by = "msu_mergers")
msu_plot_lengths[, labs:= paste0(msu_length, " bp")]

setorderv(msu_plot_lengths,cols = "order")

par(mar = c(6,5,0.2,0.1))
with(msu_plot_lengths,
     barplot(msu_length, col = cols,
             border = "white", 
             space = 0, xaxs = "i",
             ylab = "length (bp)"))
with(msu_plot_lengths, 
     axis(side = 1, las = 2, 
          tick = FALSE, cex.axis = 0.85,
          at = order - 0.5, 
          label = labs)
)

plot(core_gub_tree,
     edge.width = 2,
     show.tip.label = FALSE,
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.65,
     no.margin = TRUE, main="")

# Add labels only for focal STs
focal_tips <- which(core_gub_tree$tip.label %in% focal_ST)
tiplabels(text = core_gub_tree$tip.label[focal_tips],
          tip = focal_tips,
          frame = "none",
          adj = -0.1,              # nudges labels a bit right of tip
          col = "dodgerblue4",
          font = 2,
          cex = 0.65)

# Create empty plot
par(mar = c(0,5,0,0.1))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = c(0, nblocks * sq_size),
  ylim = c(0.5, length(focal_ST)+ 0.5), #max(yvals) 
  xlab = "", xaxs = "i",
  ylab = "",
  yaxt = "n",
  xaxt = "n"
)

t=1
# Draw squares (no gaps)
for (i in unique(msu_plot[ST %chin% focal_ST, ST])) {#unique(msu_plot$geno_id)
  dat <- msu_plot[ST==i]
  
  y <- t
  for (b in 1:max(dat$order)) {
    xleft  <- (b - 1) * sq_size
    xright <- xleft + sq_size
    rect(xleft, y - 0.5, xright, y + 0.5,
         col = dat[order==b,cols],
         border = "white")
    if (dat[order==b,strand]=="-") {
      arrows(xright, t, xleft, t,
             lwd = 3,
             length = 0.05, col = "white")
    }
  }
  label_txt <- if (i != "ST23-1LV") unique(dat$ST) else paste0(unique(dat$ST), "\n(consensus)")
  
  axis(side = 2, tick = FALSE, las = 2,
       cex.axis = 0.85,
       at = y, 
       labels = label_txt)
  t=t+1
}


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


# # ancestral reconstruction ------------------------------------------------
# tree <- read.tree("./input_data/pangraph/gub_graph.final_tree.tre")
# P <- read.table("./input_data/PIRATE_485_lng_rds_out/pres_abs_file.tsv", header=TRUE, row.names=1)
# row.names(P) = sub("_[0-9]+$", "", row.names(P))
# P <- P[tree$tip.label]
# # P = t(P)
# # P <- as.data.frame(P)
# 
# library(phytools)
# 
# # get AGs
# focal_AGs <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
#                   select = c("geno_id","asmbly_type","gene_family", "locus_tag","number_genomes"))
# 
# focal_AGs <- focal_AGs[gene_family!="" & number_genomes < 481 & 
#                          asmbly_type!= "plasmid" &
#                          geno_id %chin% tree$tip.label]
# 
# genes <- unique(focal_AGs$gene_family)
# 
# 
# # Example for one gene
# gene <- as.character(as.numeric(P["g000077", tree$tip.label]))  # use character, not factor
# names(gene) <- tree$tip.label
# 
# # Must be a character matrix with rownames matching tip labels
# gene_mat <- matrix(gene, ncol = 1)
# rownames(gene_mat) <- tree$tip.label
# 
# # Confirm it's a character matrix (not list or factor)
# str(gene_mat)
# # chr [1:50, 1] "0" "0" "0" "1" ...
# 
# # Run parsimony ancestral reconstruction
# anc <- ancestral.pars(tree, gene_mat)
# 

