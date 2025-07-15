# Find syntenic blocks


# Create igraphs for each genome ----------------------------------------------------------


# pangenome info
data_dir = "./input_data/PIRATE_i95_noplsm_1695_out"

pirate_pangenome <- fread(paste0(data_dir,"/PIRATE.gene_families.ordered.tsv"), 
                    select = c("allele_name", "gene_family", "consensus_product",
                               "number_genomes"))

# map loci
gene_families = fread(paste0(data_dir, "/PIRATE.gene_families.ordered.tsv"))
colnames = colnames(gene_families)[c(1:2, 23:ncol(gene_families))]
gene_families = gene_families[,..colnames]

gene_families = melt(gene_families, id.vars = c("allele_name", "gene_family"),
                     variable.name = "geno_id", value.name = "locus_tag")

gene_families <- gene_families[locus_tag!=""]

geno_list = as.character(unique(gene_families[,geno_id]))

# Check STs ---------------------------------------------------------------

ST_type <- fread("C:/Users/carac/Dropbox/Vos_Lab/kleb_AG_mobilisation/data/Thorpe_etal_2022_TableS2.csv",
                 select = c("id", "species_abbv","ST"))[species_abbv=="K.pne"][, species_abbv:=NULL]

ST_type <- ST_type[id %chin% geno_list]

freq_ST = ST_type[, .(n = .N), by = "ST"][which.max(n), ST]

# Create igraphs ----------------------------------------------------------

for(genome in geno_list){
  focal_genome = genome
  
  #--- Load PIRATE GFF cord info file ---
  
  cords = fread(paste0(data_dir, "/co-ords/",focal_genome,".co-ords.tab"))
  
  colnames(cords)[1] ="locus_tag"
  
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
  
  # Create grange file
  gr_cords <- dt2gr(cords)
  
  
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
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = vertex_data)
  
  
  # Save data
  write_graph(g, file = paste0("./input_data/genome_graphs/", focal_genome, ".graphml"),
              format = "graphml")
  
  print(paste0("Finished ... ", focal_genome))
  
}


# Remove high degree connections ------------------------------------------
# remove bad actors from .gfa

gfa_dir = "C:/Users/carac/Dropbox/Vos_Lab/pipeline_benchmarking_kpne_ST3/PIRATE_ST3_out/pangenome.gfa"

# Read all lines from the GFA file
gfa_lines <- readLines(gfa_dir)

# Extract links (edges): lines starting with "L"
link_lines <- gfa_lines[grepl("^L", gfa_lines)]
link_df <- read.table(text = link_lines, sep = "\t", stringsAsFactors = FALSE)

# Columns: 1=L, 2=from, 3=orientation1, 4=to, 5=orientation2, 6=cigar
colnames(link_df) <- c("type", "from", "from_orient", "to", "to_orient", "cigar")


# Create edge list
edges <- link_df[, c("from", "to")]

# Build graph
g_pan <- graph_from_data_frame(edges, directed = FALSE)


high_degree_nodes <- V(g_pan)[igraph::degree(g_pan) > 5]


g_pan <- delete_vertices(g_pan, high_degree_nodes)

# remove short chains
# Assume g is your graph
components <- components(g_pan)

# Find component IDs with 3 or more vertices
keep_ids <- which(components$csize >= 3)

# Get vertex IDs that belong to those components
vertices_to_keep <- V(g_pan)[components$membership %in% keep_ids]

# Create filtered graph
g_pan_filtered <- induced_subgraph(g_pan, vertices_to_keep)



# Get the pangenome for most common ST,  ----------------------------------

freq_ST


g_combined <- graph.empty(n = 0, directed = FALSE)

# import data and parse
for(i in geno_list) {
  
  tmp_contigs <- read_graph(file = paste0("./genome_graphs/", i , ".graphml"),
                            format = "graphml")
  
  tmp_contigs <- as_directed(tmp_contigs, mode = "mutual")  # this makes bidirectional
  
  # set focal loci
  tmp_contigs <- set_vertex_attr(tmp_contigs, name = i, value = V(tmp_contigs)$name)
  
  V(tmp_contigs)$name <- V(tmp_contigs)$gene_family
  
  focal_vertices = V(tmp_contigs)$name[V(tmp_contigs)$name %in% V(g_pan_filtered)$name]  
  
  focal_g <- induced_subgraph(tmp_contigs, focal_vertices)
  
  E(focal_g)$weight = 1
  
  if(vcount(g_combined) != 0 && ecount(g_combined) != 0){
    
    g_combined <- igraph::union(g_combined, focal_g)
    
    E(g_combined)$weight <- rowSums(
      cbind(E(g_combined)$weight_1, E(g_combined)$weight_2),
      na.rm = TRUE
    )
    
    delete_vertices <- c("locus_tag", "gene_family", "id", "gene", 
                         "product", "strand", "contig", "allele_name")
    
    g_combined <- Reduce(
      function(g, attr) delete_vertex_attr(g, name = attr),
      delete_vertices,
      init = g_combined
    )
    
    delete_edges <- c("weight_1", "weight_2", "contig")
    
    g_combined <- Reduce(
      function(g, attr) delete_edge_attr(g, name = attr),
      delete_edges,
      init = g_combined
    )
  } else {
    
    g_combined = focal_g
    
    delete_vertices <- c("locus_tag", "gene_family", "id", "gene", 
                         "product", "strand", "contig", "allele_name")
    
    g_combined <- Reduce(
      function(g, attr) delete_vertex_attr(g, name = attr),
      delete_vertices,
      init = g_combined
    )
    
    g_combined <- delete_edge_attr(g_combined, name = "contig")
  }
  
}


# Find the main path ------------------------------------------------------

fas <- feedback_arc_set(g_combined)
g_dag <- delete_edges(g_combined, fas)

# remove spurious edges
g_dag <- delete_edges(g_dag, E(g_dag)[weight <= 20])


