require(ggplot2)
require(cowplot)
require(gggenes)
require(ggdendro)
require(reshape2)
require(vegan) # for Jaccard distance

# args = commandArgs(trailingOnly=TRUE)
# 1 - gfa csv
# 2 - name of block with gene in (for centering)
# 3 - png from bandage
# 4 - output pdf

# Genome blocks approach
genome.blocks <- read.csv("./input_data/test_pangraph/graph.gfa.blocks.csv",
                          header=T, stringsAsFactors = F,
                          colClasses=c("block"="character"))
#colnames(genome.blocks) <- c("genome", "block", "strand", "start", "end", "colour")
genome.blocks$forward <- ifelse(genome.blocks$strand=="+", TRUE, FALSE)
block.counts <- table(genome.blocks$block)
blocks.which.need.colours <- names(block.counts)[which(block.counts>1)]

genome.blocks$block.coloured <- sapply(genome.blocks$block,
                                       function(x) ifelse(x %in% blocks.which.need.colours,
                                                          x,
                                                          "_other"))
block.colours <- unique(genome.blocks$colour)
names(block.colours) <- unique(genome.blocks$block.coloured)

# Get all the genome paths (in terms of blocks)
genome.paths <- sapply(unique(genome.blocks$genome), function(x) paste(genome.blocks[which(genome.blocks$genome==x), "block"], collapse=","))
genome.paths <- sort(genome.paths, decreasing = TRUE)

# Subset to unique paths
# Keep one representative genome for each
genome.path.reps <- names(genome.paths)[!duplicated(genome.paths)]
genome.blocks.unique <- genome.blocks[which(genome.blocks$genome %in% genome.path.reps),]

# And store the number of examples of each
genome.blocks.unique$genome.path <- genome.paths[genome.blocks.unique$genome]
genome.blocks.unique$n.reps <- sapply(genome.blocks.unique$genome.path,
                                      function(x) table(genome.paths)[x])
genome.blocks.unique$genome.path.name <- paste0("Type", as.numeric(as.factor(genome.blocks.unique$genome.path )))
genome.blocks.unique$genome.n <- sapply(genome.blocks.unique$n.reps, function(x) ifelse(x==1,
                                        "", paste0("n=", x)))

# Now use a dendrogram to order the genomes
m <- acast(genome ~ block, data=genome.blocks.unique, fill=0, fun.aggregate=length)
m.dist <- vegdist(m, method="jaccard") # jaccard distances based on block presence/absence
dendro <- as.dendrogram(hclust(m.dist))
dendro_order <- order.dendrogram(dendro)
genome_labels <- dendro_data(dendro)$labels$label
genome.blocks.unique$genome.ordered <- ordered(genome.blocks.unique$genome,
                                             levels=genome_labels)

# Make the linear block plot
p.blocks <- ggplot(genome.blocks.unique, aes(xmin = start, xmax = end, forward = forward, y = genome.ordered, fill = block.coloured)) +
  geom_gene_arrow(arrow_body_height = unit(2, "mm"),
                  arrowhead_height = unit(2, "mm"),
                  arrowhead_width = unit(1, "mm")) +
  theme_genes()+
  scale_fill_manual(values=block.colours)+
  ylab("")+
  theme(legend.position = "none")+
  scale_y_discrete(breaks=genome.blocks.unique$genome.ordered, labels=genome.blocks.unique$genome.n)+
  ggtitle("Linearized blocks")+
  theme(plot.title=element_text(hjust=0.5))+
  xlab("Position (bp)")
# Include the bandage graph plot
p.graph <- ggdraw() + draw_image("./input_data/test_pangraph/graph.gfa.png")+
  ggtitle("Graph representation")+
  theme(plot.title=element_text(hjust=0.5))

# Make the output pdf with both
pdf("./input_data/test_pangraph/graph.gfa_plot.pdf", height=5, width=10)
cowplot::plot_grid(p.blocks, p.graph, nrow=1, align='h')
dev.off()
