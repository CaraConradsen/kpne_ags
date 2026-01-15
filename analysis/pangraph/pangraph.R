# run pangraph to get syntenic blocks
# get minimal synteic units


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

# # take 2-5 genomes per ST (30 STs)
# correct_start[,n := .N, ST]
# correct_start_unique <- correct_start[n!=1, .SD[1:5], by = ST][!is.na(geno_id)]# subset taking first 5 STs
# 
# geno_list = correct_start_unique[,geno_id]# 115 genomes

# no clones (mash > 2e-05)

# geno_list = fread(paste0(outdir_dat, "/st_no_clones_mash2e-05.csv"))[,genomes]# 214 genomes

geno_list = read.csv(paste0(outdir_dat, "/all_no_clones_mash2e-05.csv"),
                     header = FALSE)$V1 # 214 genomes

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
# 214 seqs Time difference of 32.05445 mins

# 2. ### Generate clonalframe for trees using pangraph and gubbins

# get core
pangraph_input = file.path(input_dir, "graph.json")
set_ref_strain = "SPARK_1004_C1"
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

pdf(paste0(outdir_fig, "/msu_214_genos.pdf"), width = 14.27, height = 23.38)

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
library(ggtree)
library(ggplot2)


# detect genomes that differ in msu order/orientation
unordered_msu_genos <- as.character(
  unique(
    msu_plot[ord_msu != most_common_ord_msu | strand != "+", geno_id]))

# Setup
yvals <- 1:length(unordered_msu_genos)      # 24 unique Y values
nblocks <- nrow(msu_cols)       # 5 squares per row
# Define square size
sq_size <- 1       # both width and height = 1

# consensus "SPARK_1004_C1"
focal_geno = c("SPARK_1004_C1",
               unordered_msu_genos, "SPARK_192_C1")


pdf(paste0(outdir_fig, "/msu_subset214_unordered_msu.pdf"), width = 11.69, height = 8.74)

# plot Core units + phylogenetic tree
layout(matrix(c(1,2,3,4), nrow=2, byrow = TRUE), 
       widths=c(1,2), heights = c(2,4))  # left = tree, right = other plot
plot.new()

msu_plot_lengths <- merge(msu_plot[geno_id == "SPARK_1004_C1", 
                                   .(msu_mergers,order,cols)], 
                          msu_len_dt, all.x = TRUE, by = "msu_mergers")
msu_plot_lengths[, labs:= paste0(msu_length, " bp")]

setorderv(msu_plot_lengths,cols = "order")

par(mar = c(6,6,0.2,1))

barp_y_lim = pretty(range(msu_plot_lengths$msu_length)) 

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
          label = labs)
)
axis(side = 2, at = barp_y_lim,
     labels = barp_y_lim/1000, las = 2)

# phylogenetic tree
sub_tree <- keep.tip(core_gub_tree_geno, focal_geno)

ntip <- Ntip(sub_tree)

par(xpd = TRUE)
plot(sub_tree,
     edge.width = 2,
     y.lim = c(3.85, ntip + 0.5),
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
  label_txt <- if (!i %in% c("SPARK_1004_C1","SPARK_192_C1")) unique(dat$geno_id) else paste0(unique(dat$geno_id), "\n(consensus)")
  
  axis(side = 2, tick = FALSE, las = 2,
       cex.axis = 0.85,
       at = y, 
       labels = label_txt)
  t=t+1
}

dev.off()

