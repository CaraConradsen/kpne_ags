# Use python GFApy to create gene networks

library(reticulate)
reticulate::py_install("gfapy")

# Load gfapy Python library
gfapy <- reticulate::import("gfapy")

# Get gff data from Pirate
# Find syntenic blocks


# Create igraphs for each genome ----------------------------------------------------------

# pangenome info
data_dir = "./input_data/PIRATE_485_lng_rds_out/"

# map gene families and loci from pirate 
gene_families = fread(paste0(data_dir, "/PIRATE.gene_families.ordered.tsv"))
colnames = colnames(gene_families)[c(1:2, 23:ncol(gene_families))]
gene_families = gene_families[,..colnames]

gene_families = melt(gene_families, id.vars = c("allele_name", "gene_family"),
                     variable.name = "geno_id", value.name = "locus_tag")

gene_families <- gene_families[locus_tag!=""]

geno_list = as.character(unique(gene_families[,geno_id]))[1:100]# depth of 10 genomes

# determine core genes
pirate_pangenome <- fread(paste0(data_dir,"/PIRATE.gene_families.ordered.tsv"), 
                          select = c("allele_name", "gene_family", "consensus_product",
                                     "number_genomes"))
core_genes <- pirate_pangenome[number_genomes == 485, gene_family]


# Create import cords for dt ----------------------------------------------------------

pangenomes_dt <- foreach(genome = geno_list,
                         .packages = c("data.table")) %do% {
  focal_genome = genome
  
  #--- Load PIRATE GFF cord info file ---
  
  cords = fread(paste0(data_dir, "/co-ords/",focal_genome,".co-ords.tab"),
                select = c("Name", "Contig", "Start", "End", "Strand"))
  
  #tidy names
  colnames(cords) = tolower(colnames(cords))
  colnames(cords)[1] = "locus_tag"
  
  # add gene families
  cords = merge(cords, gene_families[geno_id == focal_genome,
                                     .(gene_family, locus_tag)],
                all.x = TRUE, by = "locus_tag")
  
  # Convert strand to "+" / "-"
  cords$strand <- ifelse(cords$strand == "Forward", "+", "-")
  
  # add genome id
  cords$geno_id = focal_genome
  
  # remove missing loci (short proteins)
  cords <- cords[!is.na(gene_family)]
  
  cords[, .(geno_id, locus_tag, contig, start, end, strand, gene_family)] #first thousand genes [1:500,]
}

# collapse list
pangenomes_dt <- rbindlist(pangenomes_dt)

# remove plasmid contigs
plasmid_contigs <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"), 
                         select = c("seqnames", "asmbly_type")))[asmbly_type == "chromosome"]

pangenomes_dt <- pangenomes_dt[contig %chin% plasmid_contigs$seqnames]

# Sort genes by genome and genomic position
setorderv(pangenomes_dt, c("geno_id", "contig", "start"))

# Create GFA network ------------------------------------------------------
# Helper for strand formatting
strand_symbol <- function(s) ifelse(s == "-", "-", "+")

# --- Create new GFA object ---
gfa <- gfapy$Gfa()

# --- Step 1: Define segments (orthogroups) ---
gene_fam <- unique(pangenomes_dt$gene_family)
for (fam in gene_fam) {
  color_tag <- ifelse(fam %in% core_genes, "core", "accessory")  # tag core genes
  # Add an optional GFA tag 'CL' (color label) for Cytoscape styling
  gfa$add_line(sprintf("S\t%s\t*\tCL:Z:%s", fam, color_tag))
}

# --- Step 2: Add links (adjacent genes per genome, respecting strand) ---
links <- pangenomes_dt[, .(
  from         = gene_family[-.N],
  # from_strand  = strand_symbol(strand[-.N]),
  to           = gene_family[-1],
  # to_strand    = strand_symbol(strand[-1])
), by = geno_id]

# Deduplicate adjacencies
links_unique <- unique(links[, .(from, from_strand, to, to_strand)])

for (i in seq_len(nrow(links_unique))) {
  gfa$add_line(sprintf(
    "L\t%s\t%s\t%s\t%s\t0M",
    links_unique$from[i],
    links_unique$from_strand[i],
    links_unique$to[i],
    links_unique$to_strand[i]
  ))
}


# --- Step 3: Add genome paths (ordered by coordinates and strand) ---
for (g in unique(pangenomes_dt$geno_id)) {
  sub <- pangenomes_dt[geno_id == g]
  path <- paste0(sub$gene_family, strand_symbol(sub$strand), collapse = ",")
  gfa$add_line(sprintf("P\t%s\t%s\t*", g, path))
}

node_attr <- data.frame(
  gene_family = pangenomes_dt$gene_family,
  color_group = ifelse(pangenomes_dt$gene_family %in% core_genes, "core", "accessory"),
  strand = pangenomes_dt$strand,
  stringsAsFactors = FALSE
)

node_attr$gene_family = paste0(node_attr$gene_family, node_attr$strand)

node_attr = unique(node_attr)
node_attr$strand = NULL

write.table(node_attr, "node_attributes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)



# Get path ----------------------------------------------------------------
# Extract paths
paths <- gfa$paths  # returns a dictionary of path objects in gfapy

# Build adjacency counts
adj_counts <- data.table(from=character(), to=character(), count=integer())

for (g in unique(pangenomes_dt$geno_id)) {
  sub <- pangenomes_dt[geno_id == g]
  segs <- paste0(sub$gene_family, strand_symbol(sub$strand))
  from <- segs[-length(segs)]
  to   <- segs[-1]
  adj_counts <- rbind(adj_counts, data.table(from=from, to=to, count=1))
}

# Sum counts
adj_counts <- adj_counts[, .(count=sum(count)), by=.(from, to)]
setorder(adj_counts, -count)

adj_counts[grepl("g001385", adj_counts$from)]

# Build the consensus path
# Create a weighted directed graph
g <- graph_from_data_frame(adj_counts, directed = TRUE)

# Option 1: use a simple greedy approach: pick the heaviest outgoing edge at each step
start_nodes <- setdiff(V(g)$name, adj_counts$to)  # nodes with no incoming edges
consensus_path <- c()

current <- start_nodes[1]  # pick one start
while (!is.null(current) && !(current %in% consensus_path)) {
  consensus_path <- c(consensus_path, current)
  
  # pick the most frequent outgoing edge
  out_edges <- adj_counts[from == current]
  if (nrow(out_edges) == 0) break
  current <- out_edges$to[which.max(out_edges$count)]
}

consensus_gfa_path <- paste0(consensus_path, collapse = ",")
gfa$add_line(sprintf("P\tconsensus_backbone\t%s\t*", consensus_gfa_path))


# --- Step 4: Save GFA file ---
out_path <- "gene_order_graph.gfa"
gfa$to_file(out_path)
