# plot pangenome graphs
plot_g = focal_g

V(plot_g)$color <- "dodgerblue"
# V(plot_g)$color[V(plot_g)$name %in% target_genes] <- "green2"
V(plot_g)$color[V(plot_g)$name %in% core_genes] <- "brown3"#core_genes2

# plot(plot_g,
#      layout=layout_nicely,
#      edge.color = "grey20",
#      edge.width = (E(plot_g)$weight),
#      vertex.label = NA,
#      # vertex.label.cex = 0.5,
#      # vertex.label.dist = 0.5,
#      vertex.color = V(plot_g)$color,
#      vertex.size=3, 
#      edge.arrow.size=0.5)
# # vertex.color=my_color, 
# # vertex.label.cex=0.7,
# # vertex.label.color="white",
# # vertex.frame.color="transparent")



# # Visualise graphs --------------------------------------------------------
# cytoscapePing()# check cytoscape
createNetworkFromIgraph(plot_g,"SPARK")
