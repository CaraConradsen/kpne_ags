# Find syntenic blocks


# Create igraphs for each genome ----------------------------------------------------------

# pangenome info
data_dir = "./input_data/PIRATE_485_lng_rds_out/"

pirate_pangenome <- fread(paste0(data_dir,"/PIRATE.gene_families.ordered.tsv"), 
                    select = c("allele_name", "gene_family", "consensus_product",
                               "number_genomes", "cluster", "cluster_order"))

# map loci
gene_families = fread(paste0(data_dir, "/PIRATE.gene_families.ordered.tsv"))
colnames = colnames(gene_families)[c(1:2, 23:ncol(gene_families))]
gene_families = gene_families[,..colnames]

gene_families = melt(gene_families, id.vars = c("allele_name", "gene_family"),
                     variable.name = "geno_id", value.name = "locus_tag")

gene_families <- gene_families[locus_tag!=""]

geno_list = as.character(unique(gene_families[,geno_id]))[1:10]

# identify core gene clusters
core_clusters <- pirate_pangenome[number_genomes == 485]
core_clusters[, cl_length := .N, by = cluster]
core_genes <- as.character(unique(pirate_pangenome[number_genomes == 485, cluster]))

# Create igraphs ----------------------------------------------------------

for(genome in geno_list){
  focal_genome = genome
  
  #--- Load PIRATE GFF cord info file ---
  
  cords = fread(paste0(data_dir, "/co-ords/",focal_genome,".co-ords.tab"))
  
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
  
  ########## REMOVE MIDDLE CORE GENES
  cords <- cords[!gene_family %chin% core_clusters[cluster_order!=1, gene_family]]
  cords[core_clusters, on = "gene_family", cluster := i.cluster]
  cords[!is.na(cluster), gene_family := cluster]
  cords[, cluster:=NULL]
  
  # Create grange file
  gr_cords <- dt2gr(cords)
  
  # Find contig lengths (max end - min start per contig)
  contig_lengths <- tapply(width(ranges(gr_cords)), seqnames(gr_cords), sum)
  
  gr_cords <- gr_cords[seqnames(gr_cords) == names(which.max(contig_lengths))]
  
  # Build igraph ------------------------------------------------------------
  
  # Ensure GRanges is sorted
  gr_cords <- sort(gr_cords)
  
  # Create globally unique node names by combining allele and contig
  node_ids <- paste0(gr_cords$locus_tag)
  
  # Vertex data
  vertex_data <- data.frame(
    name = node_ids,
    gene_family = gr_cords$gene_family,
    locus_tag = gr_cords$locus_tag,
    product = gr_cords$product,
    strand = as.character(strand(gr_cords)),
    contig = as.character(seqnames(gr_cords)),
    stringsAsFactors = FALSE
  )
  
  # Build edge list
  edges <- do.call(rbind, lapply(split(gr_cords, seqnames(gr_cords)), function(contig_genes) {
    if (length(contig_genes) < 2) return(NULL)
    allele_ids <- contig_genes$locus_tag
    from <- allele_ids[-length(allele_ids)]
    to   <- allele_ids[-1]
    data.frame(from = from, to = to, contig = as.character(seqnames(contig_genes))[1], stringsAsFactors = FALSE)
  }))
  
  # Now build the graph
  g <- graph_from_data_frame(edges, vertices = vertex_data, directed = TRUE)
  
  g <- as_directed(g, mode = "mutual")
  
    # Save data
  write_graph(g, file = paste0("./input_data/genome_graphs/", focal_genome, ".graphml"),
              format = "graphml")
  
  print(paste0("Finished ... ", focal_genome))
  
}


# Get core linkers and assembly information -------------------------------
ln_rd_pangenome_info <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"),
                              select = c("geno_id","locus_tag","number_genomes", "asmbly_type"))



# Create directed paths with long reads -----------------------------------
# Get core genome ----------------------------------

# Set up parallel backend
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)

# Parallel graph loading and processing
graph_list <- foreach(i = geno_list, .packages = "igraph") %do% {

  focal_g <- read_graph(file = paste0("./input_data/genome_graphs/", i, ".graphml"),
                            format = "graphml")
  
  p_names <- ln_rd_pangenome_info[geno_id == i & (asmbly_type=="plasmid"), locus_tag]
  
  # remove core and plasmid
  keep_vertices <- V(focal_g)[!name %in% p_names]

  # Create a filtered subgraph
  focal_g <- induced_subgraph(focal_g, vids = keep_vertices)

  V(focal_g)$name <- V(focal_g)$gene_family

  E(focal_g)$weight <- 1

  # Remove unnecessary vertex
  delete_vertices <- c("locus_tag", "gene_family", "id", "product", "strand", "contig")
  focal_g <- Reduce(function(gx, attr) delete_vertex_attr(gx, name = attr), delete_vertices, focal_g)

  focal_g
}

# stopCluster(cl)  # Clean up

names(graph_list) = geno_list

# Extract syntenic blocks
get_syntenic_blocks <- function(g) {
  comps <- components(g)
  split(V(g)$name, comps$membership)
}

syntenic_blocks <- lapply(graph_list, get_syntenic_blocks)

# Flatten all gene sets into a list
syntenic_blocks <- unlist(syntenic_blocks, recursive = FALSE)

syntenic_blocks <- Filter(function(x) length(x) > 20, syntenic_blocks)


jaccard <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}


# Compute pairwise Jaccard distances
n <- length(syntenic_blocks)
dist_mat <- matrix(0, n, n)
rownames(dist_mat) = names(syntenic_blocks)
colnames(dist_mat) = names(syntenic_blocks)

# build a distance matrix
for(i in 1:n) {
  for(j in 1:n) {
    dist_mat[i,j] <- 1 - jaccard(syntenic_blocks[[i]], syntenic_blocks[[j]])  # distance = 1-similarity
  }
}

# Cluster the sets
hc <- hclust(as.dist(dist_mat), method = "average")

# Cutting the tree at a threshold gives clusters of similar gene sets:
clusters <- cutree(hc, h = 0.99)  # adjust h for similarity threshold
group_blocks <- split(names(syntenic_blocks), clusters)
group_blocks <- Filter(function(x) length(x) == 6, group_blocks)


subset_g <- foreach(i = 1:length(syntenic_blocks[c(group_blocks$`1`)]),
                    .packages = "igraph") %do% {
                      geno_id <- sub("\\..*$", "",names(syntenic_blocks[c(group_blocks$`1`)][i]))
                      
                      focal_genes = syntenic_blocks[c(group_blocks$`1`)][i][[1]]
                      
                      keep_vertices <- V(graph_list[geno_id][[1]])[name %in% focal_genes]
                      
                      # Create a filtered subgraph
                      focal_net <- induced_subgraph(graph_list[geno_id][[1]], vids = keep_vertices)
                      
                      focal_net
                    }


