


# Import pirate gff information -------------------------------------------

i95_cords_files <- list.files("./input_data/PIRATE_i95_noplsm_1695_out/co-ords",
                              recursive = TRUE, full.names = TRUE)


cl <- makeCluster(num_cores)
registerDoParallel(cl)

kpne_genes <- foreach(i = 1:length(i95_cords_files),
                   .combine = rbind,
                   .packages = "data.table") %dopar% {
                     
                     fread(i95_cords_files[i])[,TYPE:=NULL]
                   } 

stopCluster(cl)

colnames(kpne_genes) <- tolower(colnames(kpne_genes))
colnames(kpne_genes)[1] = "locus_tag"

# Add PIRATE info ---------------------------------------------------------
pangenome = fread("./input_data/PIRATE_i95_noplsm_1695_out/PIRATE.gene_families.ordered.tsv")

colnames = colnames(pangenome)[c(2,4:5,7, 23:ncol(pangenome))]

pangenome = pangenome[,..colnames]

pangenome = melt(pangenome, id.vars = c("gene_family",
                                        "consensus_product", "threshold", "number_genomes"),
                     variable.name = "geno_id", value.name = "locus_tag")

pangenome = pangenome[locus_tag != ""]

# add allele information
pangenome_alleles = fread("./input_data/PIRATE_i95_noplsm_1695_out/PIRATE.unique_alleles.tsv")


colnames = colnames(pangenome_alleles)[c(1:2,5,21:ncol(pangenome_alleles))]

pangenome_alleles = pangenome_alleles[,..colnames]

pangenome_alleles = melt(pangenome_alleles, 
                         id.vars = c("allele_name","gene_family", "threshold"),
                 variable.name = "geno_id", value.name = "locus_tag")

pangenome_alleles = pangenome_alleles[locus_tag != ""]


pangenome <- merge(pangenome, pangenome_alleles, all.x = TRUE, 
                   by = c("gene_family","geno_id","locus_tag", "threshold"))

# free up memory
rm(pangenome_alleles)

pangenome[, freq:= number_genomes/1695]

pangenome[, pan_grp := fcase(freq >=0.99, "core",
                             freq < 0.99 & freq >= 0.95, "soft-core",
                             freq < 0.95 & freq >= 0.15, "shell",
                                      default = "cloud")]

pangenome <- merge(kpne_genes, pangenome,
                   all.x = TRUE, by = "locus_tag")

# fwrite(pangenome, paste0(outdir_dat, "/pangenome.csv"))



# Add discordant measures -------------------------------------------------
# Load libraries
library(ape)
library(phangorn)
library(phytools)
library(castor)

# Get newick tree
# Newick tree from PathogenWatch
kpne_tree <- read.tree("C:/Users/carac/Dropbox/Vos_Lab/kleb_AG_mobilisation/data/k.pne_spec_tree.nwk")

kpne_tree <- keep.tip(kpne_tree, 
                      as.character(unique(pangenome[!is.na(geno_id), geno_id])))


# 1. Load presence/absence matrix
pa_matrix <- read.table("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/presence_absence/gene_presence_absence.tsv", 
                        header = TRUE, row.names = 1, check.names = FALSE)

# convert rownames to gene_families
rownames(pa_matrix) <- gsub("_[0-9]+$", "", rownames(pa_matrix))

# keep only AGs
pa_matrix <- pa_matrix[rownames(pa_matrix) %in% unique(pangenome[pan_grp!="core",gene_family]), ]

# Transpose if genes are rows
pa_matrix <- t(pa_matrix)

# Calculate parsimony scores per gene
parsimony_scores <- apply(pa_matrix, 2, function(gene_col) {
  gene_col <- as.factor(gene_col)
  names(gene_col) <- rownames(pa_matrix)
  phyDat_ag <- phyDat(as.matrix(gene_col), type = "USER", levels = c("0", "1"))
  parsimony(kpne_tree, phyDat_ag)
})

# Calculate Consistency Index (CI) per gene
CI_scores <- apply(pa_matrix, 2, function(gene_col) {
  gene_col <- as.factor(gene_col)
  names(gene_col) <- rownames(pa_matrix)
  phyDat_ag <- phyDat(as.matrix(gene_col), type = "USER", levels = c("0", "1"))
  CI(kpne_tree, phyDat_ag)
})

# plots
hist(parsimony_scores, main="Parsimony scores per gene", breaks = 100)
hist(CI_scores, main="Parsimony scores per gene", breaks = 100)

for (i in names(parsimony_scores[parsimony_scores <50 & parsimony_scores >1 ])) {
  print(unique(pangenome[gene_family == i, c(gene_family, product)]))
}

for (i in names(CI_scores[CI_scores == 0.5])) {
  print(unique(pangenome[gene_family == i, c(gene_family, product, pan_grp, number_genomes)]))
}


# g004636, g002724_1


# Syntenic networks 1. Get Dups-------------------------------------------------------
gene_strings <- apply(pa_matrix, 2, paste, collapse = "_")

# Group by identical pattern
dup_sets <- split(seq_along(gene_strings), gene_strings)

# get gene names
dup_sets <- lapply(dup_sets, function(x) colnames(pa_matrix)[x])

dup_sets_ls <- lapply(names(dup_sets), function(set_name) {
  set <- dup_sets[[set_name]]
  data.table(
    number_genomes = sum(strsplit(set_name, "_")[[1]] == "1"),
    n_genes = length(set), 
    genes = paste(set, collapse = ",")
  )
})

# Convert to data.frame
dup_df <- rbindlist(dup_sets_ls)


# Syntenic networks 2. create igraphs -------------------------------------
focal_genes = c("g001872","g002450","g002519","g002538","g003816","g003875","g004532","g004534","g005139","g005922","g006148")

cords = pangenome[gene_family %in% focal_genes]

# Convert strand to "+" / "-"
cords$strand <- ifelse(cords$strand == "Forward", "+", "-")

# set contig as seqname
colnames(cords)[which(names(cords)=="contig")] = "seqnames"

# fix missing loci
cords[is.na(gene_family), `:=`(gene_family = paste0(product,gsub(focal_genome,"", locus_tag)))]

# Create grange file
gr_cords <- dt2gr(cords)

# Ensure GRanges is sorted
gr_cords <- sort(gr_cords)

# Create globally unique node names by combining allele and contig
node_ids <- paste0(gr_cords$locus_tag)

# Vertex data
vertex_data <- data.frame(
  name = node_ids,
  gene_family = gr_cords$gene_family,
  locus_tag = gr_cords$locus_tag,
  gene = gr_cords$gene,
  product = gr_cords$product,
  strand = as.character(strand(gr_cords)),
  contig = as.character(seqnames(gr_cords)),
  allele_name = gr_cords$allele_name,
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
g <- graph_from_data_frame(edges, directed = TRUE, vertices = vertex_data)

E(g)$weight = 1
g <- decompose(g, mode = "weak")

for(i in 1: length(g)){
  g_tmp =  g[[i]]

}







