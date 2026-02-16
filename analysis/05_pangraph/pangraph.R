# run pangraph to get syntenic blocks
# get minimal syntenic units

# Retain only correct oriented chromosomes --------------------------------
pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))

geno_list = unique(pan_anno[,geno_id])

gene_family = unique(pan_anno[,.(gene_family)])

# get fasta files

hybrid_fasta_files = list.files("./input_data/kpne_260_chr_fasta/",
                                full.names = TRUE)

pangraph_dir =  "./input_data/pangraph/"

# get minimum gene length for pangenome
min_gene_length <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
                            select = c("geno_id", "start", "end", "strand", "gene_family"))
min_gene_length <- min_gene_length[geno_id %chin% geno_list] # focal genomes
min_gene_length <- min_gene_length[grepl("g", gene_family)]# discard small proteins not used by PIRATE
colnames(min_gene_length)[1]="seqnames"
min_gene_length <- dt2gr(min_gene_length)
min_gene_length = min(width(min_gene_length))


# Run Noll et al.'s Pangraph ----------------------------------------------

# 1. ### Generate pangraph

# ! Make sure fasta headers are the genome name and not contig #

pangraph_path = "/home/carac/anaconda3/bin/pangraph" # pangraph program dir 
threads = 20 # Number of threads
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
# 214 seqs Time difference of 32.05445 mins
# 260 seqs Time difference of 46.88418 mins

# 2. ### Generate clonalframe for trees using pangraph and gubbins

# get core
pangraph_input = file.path(input_dir, "graph.json")
set_ref_strain = "SPARK_1006_C1"
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
iterations = 20
threads = 22 # number of threads
tree_method = "fasttree" # tree builder to use
# input alignment file (from pangraph output)
input_alignment = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/core_genome_aln.fa"
# output prefix 
output_prefix = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/gub_graph"

# build gubbins command string 
cmd_gubbins <- sprintf(
  "wsl %s run -n %s run_gubbins.py --iterations %d --threads %d -t %s --prefix %s %s",
  conda_path,
  conda_env,
  iterations,
  threads,
  tree_method,
  output_prefix,
  input_alignment
)

# Print command for verification 
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
# 214 seqs Time difference of 2.804536 hours
# 260 seqs Time difference of 5.134598 hours

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
core_gub_tree_geno <- read.tree("./input_data/pangraph/gub_graph.node_labelled.final_tree.tre") 

core_gub_tree <- core_gub_tree_geno

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

most_common_ord_msu <- msu_paths_dt[, .N, by = ord_msu][which.max(N), ord_msu]

msu_cols = data.frame(msu_mergers = unlist(strsplit(most_common_ord_msu, ",")),
                      cols = colorRampPalette(c("#420D55","blueviolet",
                                                "darkblue","#008080",
                                                "forestgreen","#799000",
                                                "yellow","#FF8000",
                                                 "brown"))(num_msu_blocks))

msu_plot <- merge(msu_paths_dt, msu_cols, all.x = TRUE, by = "msu_mergers")

# change to ST and set order
msu_plot <- merge(msu_plot, ST_labels, all.x = TRUE, by = "geno_id")

# Create a factor with levels matching the tree tip order
msu_plot[, geno_id  := factor(geno_id, levels = core_gub_tree_geno$tip.label)]

# Now use setorderv — it will respect the factor level order
setorderv(msu_plot, c("geno_id", "order"))


# All genomes msu plot ----------------------------------------------------

png(filename = paste0(outdir_fig, "/msu_260_genos.png"),
    width = 14.27,height = 23.38,
    units = "in",res = 300)

# Setup
yvals <- 1:length(unique(msu_plot$geno_id))      # 24 unique Y values
nblocks <- nrow(msu_cols)       # 5 squares per row
# Define square size
sq_size <- 1       # both width and height = 1

# plot Core units + phylogenetic tree
layout(matrix(c(1,2), nrow=1, byrow = TRUE), 
       widths=c(1,1))  # left = tree, right = other plot

plot(core_gub_tree_geno,
     edge.width = 0.75,
     show.tip.label = FALSE,
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.15,
     no.margin = TRUE, main="")

# Add labels 
tiplabels(text = core_gub_tree_geno$tip.label,
          # tip = focal_tips,
          frame = "none",
          adj = -0.1,              # nudges labels a bit right of tip
          col = "dodgerblue4",
          font = 2,
          cex = 0.45)

# Create empty plot
par(mar = c(0,4,0,1))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = c(0, nblocks * sq_size),
  ylim = c(0.5, length(geno_list)+ 0.5), #max(yvals) 
  xlab = "", xaxs = "i",
  ylab = "",
  yaxt = "n",
  xaxt = "n"
)

t=1
# Draw squares (no gaps)
for (i in core_gub_tree_geno$tip.label) {#unique(msu_plot$geno_id)
  dat <- msu_plot[geno_id==i]
  
  y <- t
  for (b in 1:max(dat$order)) {
    xleft  <- (b - 1) * sq_size
    xright <- xleft + sq_size
    rect(xleft, y - 0.5, xright, y + 0.5,
         col = dat[order==b,cols],
         border = "white")
    if (dat[order==b,strand]=="-") {
      arrows(xright, t, xleft, t,
             lwd = 0.5,
             length = 0.05, col = "white")
    }
  }
  axis(side = 2, tick = FALSE, las = 2,
       cex.axis = 0.35,
       at = y, 
       labels = i)
  t=t+1
}

dev.off()

# Plot to closer inspect novel genome order ------------------------------------------

# detect genomes that differ in msu order/orientation
unordered_msu_genos <- as.character(
  unique(
    msu_plot[ord_msu != most_common_ord_msu | strand != "+", geno_id]))

# Setup
yvals <- 1:length(unordered_msu_genos)      # 24 unique Y values
nblocks <- nrow(msu_cols)       # 5 squares per row
# Define square size
sq_size <- 1       # both width and height = 1

# consensus "SPARK_1006_C1"
focal_geno = c("SPARK_1332_C1","SPARK_1006_C1",
               unordered_msu_genos, "SPARK_587_C1")

png(filename = paste0(outdir_fig, "/msu_subset260_unordered_msu.png"),
    width = 9.69, height = 10.74, type = "cairo",
    units = "in",res = 300)


# plot Core units + phylogenetic tree
layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE), 
       widths=c(1,2), heights = c(2,4))  # left = tree, right = other plot
plot.new()

msu_plot_lengths <- merge(msu_plot[geno_id == "SPARK_1006_C1", 
                                   .(msu_mergers,order,cols)], 
                          msu_len_dt, all.x = TRUE, by = "msu_mergers")
msu_plot_lengths[, labs:= paste0(msu_length, " bp")]

setorderv(msu_plot_lengths,cols = "order")

par(mar = c(4,6,0.2,1))

barp_y_lim = pretty(range(msu_plot_lengths$msu_length))[-7] 

with(msu_plot_lengths,
     barplot(msu_length, col = cols,
             border = "white", 
             space = 0, xaxs = "i",
             yaxt="n",
             ylab = "length (Kbp)"))
with(msu_plot_lengths, 
     axis(side = 1, las = 2, 
          tick = FALSE, cex.axis = 0.85,
          at = order - 0.5, 
          label = labs))

axis(side = 2, at = barp_y_lim,
     labels = barp_y_lim/1000, las = 2)

# phylogenetic tree
sub_tree <- keep.tip(core_gub_tree_geno, focal_geno)

ntip <- Ntip(sub_tree)

par(xpd = TRUE, mar = c(4,5,0.2,0))
plot(sub_tree,
     edge.width = 2,
     y.lim = c(3.25, ntip + 0.25),
     x.lim = c(5000, max(node.depth.edgelength(sub_tree))+100),
     show.tip.label = FALSE,
     edge.color = "darkgrey",
     tip.color = "dodgerblue4",
     font = 1,
     cex = 0.65,
     # no.margin = TRUE,
     main="")

# Create empty plot
par(mar = c(0,6,0,1))
plot(
  0, 0, type = "n",
  bty= "n",
  xlim = c(0, nblocks * sq_size),
  ylim = c(0.5, length(focal_geno)+ 0.5), #max(yvals) 
  xlab = "", xaxs = "i",
  ylab = "",
  yaxt = "n",
  xaxt = "n"
)

t=1
# Draw squares (no gaps)
for (i in sub_tree$tip.label) {
  dat <- msu_plot[geno_id==i]
  
  y <- t
  for (b in 1:max(dat$order)) {
    xleft  <- (b - 1) * sq_size
    xright <- xleft + sq_size
    rect(xleft, y - 0.5, xright, y + 0.5,
         col = dat[order==b,cols],
         border = "white")
    if (dat[order==b,strand]=="-") {
      arrows(xright, t, xleft, t,
             lwd = 1,
             length = 0.05, col = "white")
    }
  }
  label_txt <- if (!i %in% c("SPARK_1006_C1","SPARK_587_C1","SPARK_1332_C1")) unique(dat$geno_id) else paste0(unique(dat$geno_id), "\n(consensus)")
  
  axis(side = 2, tick = FALSE, las = 2,
       cex.axis = 0.85,
       at = y, 
       labels = label_txt)
  t=t+1
}

dev.off()

