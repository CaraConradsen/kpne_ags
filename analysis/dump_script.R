

# # Add MGE information -----------------------------------------------------
# # prophage
# phastest_files = list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/phastest_1695_results",
#                             full.names = TRUE)
# phastest_files = phastest_files[!grepl(".zip", phastest_files)]
# 
# cl <- makeCluster(num_cores)
# registerDoParallel(cl)
# 
# phages_dt <- foreach(i = 1:length(phastest_files),
#                       .combine = rbind,
#                       .packages = c("data.table",
#                                     "jsonlite")) %dopar% {
#                         
#                         fromJSON(paste0(phastest_files[i], "/predicted_phage_regions.json"))
#                       } 
# 
# stopCluster(cl)
# 
# setDT(phages_dt)
# 
# # fwrite(phages_dt, paste0(outdir_dat, "/phages_dt.csv"))
# 
# phages_dt <- phages_dt[, c("region", "start", "stop", "GC", "most_common_phage") := NULL]
# 
# colnames(phages_dt) = c("contig", "start", "end", "phage")
# 
# phages_dt$phage = paste0("phg_", phages_dt$phage)
# 
# 
# # Set keys for foverlaps
# setkey(phages_dt, contig, start, end)
# 
# # Create a 'dummy' column in genes so it can be used with foverlaps
# pangenome[, phage := NA_character_]  # placeholder
# 
# # foverlaps requires genes to have start and end
# overlaps <- foverlaps(pangenome, phages_dt, by.x = c("contig", "start", "end"), 
#                       by.y = c("contig", "start", "end"), nomatch = 0L)
# 
# pangenome[,phage := NULL]
# 
# pangenome <- merge(pangenome, overlaps[,.(locus_tag, phage)], 
#                      by = c("locus_tag"), 
#                      all.x = TRUE)


# 
# # set unique loci
# find_loci_g <- focal_graphs[[1]]
# 
# V(find_loci_g)$name <- V(find_loci_g)$allele_name
# 
# # Find weakly connected components (treat direction as undirected)
# comps <- components(find_loci_g, mode = "weak")
# 
# unique_grps <- lapply(1:max(comps$membership), function(i){
#   tmp_comps = names(comps$membership)[comps$membership == i]
#   tmp_comps  <- tmp_comps[tmp_comps %in% target_genes]
#   tmp_comps  <- tmp_comps[1]
# }
# )
# 
# unique_grps <- unlist(unique_grps)
# 
# membership_list <- list()
# # further divide graphs to subgraphs 


# # Canonical path as ordered gene family names
# canon_path <- topo_sort(graph_list[[1]], mode = "out")$gene_family
# 
# canon_edges <- data.frame(
#   from = head(canon_path, -1),
#   to   = tail(canon_path, -1)
# )
# 

# # orient g_combined
# combined_edges <- as_data_frame(g_combined, what = "edges")
# 
# # Check if each edge is in the same direction as canonical
# canon_edge_pairs <- paste(canon_edges$from, canon_edges$to, sep = "→")
# combined_pairs <- paste(combined_edges$from, combined_edges$to, sep = "→")
# 
# # For reverse matches (i.e., edge exists but backward)
# reverse_pairs <- paste(combined_edges$to, combined_edges$from, sep = "→")
# reversed <- reverse_pairs %in% canon_edge_pairs
# 
# # Flip edges that are backwards
# combined_edges_flipped <- combined_edges
# combined_edges_flipped[reversed, c("from", "to")] <- combined_edges_flipped[reversed, c("to", "from")]
# 
# g_oriented <- graph_from_data_frame(combined_edges_flipped, directed = TRUE)
# E(g_oriented)$weight <- combined_edges_flipped$weight
# E(g_oriented)$canonical_oriented <- combined_edges_flipped$canonical_oriented


# focal_subgraphs <- lapply(focal_graphs, function(g_tmp) {
#   
#   V(g_tmp)$name <- V(g_tmp)$allele_name
#   
#   # sub_mem <- components(g_tmp)$membership
#   # 
#   # sub_mem <- sub_mem[names(sub_mem) %in% unique_grps]
#   # 
#   # membership_list <- c(membership_list, list(sub_mem))
#   
#   # decompose networks into separate chunks
#   subgraphs <- decompose(g_tmp, mode = "weak")
#   
#   subgraphs
# })
# 
# 
# similarity_score <- function(g1, g2) {
#   score <- 0
#   score <- score + abs(vcount(g1) - vcount(g2))
#   score <- score + abs(ecount(g1) - ecount(g2))
#   return(score)
# }
# 
# # Get all (i, j) indices
# graph_indices <- list()
# for (i in seq_along(focal_subgraphs)) {
#   for (j in seq_along(focal_subgraphs[[i]])) {
#     graph_indices[[length(graph_indices) + 1]] <- list(i = i, j = j)
#   }
# }
# 
# 
# # Generate similarity matrix
# n <- length(graph_indices)
# similarity_matrix <- matrix(NA, n, n)
# 
# for (a in 1:(n-1)) {
#   for (b in (a+1):n) {
#     i1 <- graph_indices[[a]]$i
#     j1 <- graph_indices[[a]]$j
#     i2 <- graph_indices[[b]]$i
#     j2 <- graph_indices[[b]]$j
#     
#     g1 <- focal_subgraphs[[i1]][[j1]]
#     g2 <- focal_subgraphs[[i2]][[j2]]
#     
#     similarity_matrix[a, b] <- similarity_score(g1, g2)
#     similarity_matrix[b, a] <- similarity_matrix[a, b]
#   }
#   similarity_matrix[a, a] <- 0
# }
# similarity_matrix[n, n] <- 0
# 
# # culster similar graphs
# dist_matrix <- as.dist(similarity_matrix)
# hc <- hclust(dist_matrix, method = "average")
# 
# # Choose a similarity threshold
# threshold <- 5  # You can tune this
# clusters <- cutree(hc, h = threshold)
# 
# clustered_graph_refs <- split(graph_indices, clusters)
# 
# # Optionally, you can get the graphs too:
# clustered_graphs <- lapply(clustered_graph_refs, function(refs) {
#   lapply(refs, function(ref) focal_subgraphs[[ref$i]][[ref$j]])
# })
# 
# merge_graphs_weighted <- function(graph_list) {
#   if (length(graph_list) == 0) return(make_empty_graph())  # empty cluster
#   
#   all_edges <- character()
#   
#   for (g in graph_list) {
#     if (ecount(g) == 0) next  # skip graphs with no edges
#     
#     el <- as_edgelist(g, names = TRUE)
#     edge_keys <- apply(el, 1, function(e) paste(sort(e), collapse = "_"))  # for undirected
#     all_edges <- c(all_edges, edge_keys)
#   }
#   
#   if (length(all_edges) == 0) {
#     # Return a merged graph of isolated nodes if no edges at all
#     all_nodes <- unique(unlist(lapply(graph_list, function(g) V(g)$name)))
#     return(make_empty_graph(n = length(all_nodes), directed = FALSE) %>%
#              set_vertex_attr("name", value = all_nodes))
#   }
#   
#   # Count edge frequencies
#   edge_table <- table(all_edges)
#   edge_df <- data.frame(
#     edge = names(edge_table),
#     weight = as.numeric(edge_table),
#     stringsAsFactors = FALSE
#   )
#   edge_df <- transform(edge_df,
#                        from = sapply(strsplit(edge, "_"), `[`, 1),
#                        to   = sapply(strsplit(edge, "_"), `[`, 2))
#   
#   graph_from_data_frame(edge_df[, c("from", "to", "weight")], directed = FALSE)
# }
# 
# 
# 
# merged_graphs <- lapply(clustered_graphs, merge_graphs_weighted)
# 
# 
# result <- lapply(focal_graphs, function(g_tmp) {
#   V(g_tmp)$name <- V(g_tmp)$allele_name
#   
#   sub_mem <- components(g_tmp)$membership
#   sub_mem <- sub_mem[names(sub_mem) %in% unique_grps]
#   
#   subgraphs <- decompose(g_tmp, mode = "weak")
#   
#   list(
#     membership = sub_mem,
#     subgraphs = subgraphs
#   )
# })
# 
# # Extract all membership vectors
# membership_list <- lapply(result, `[[`, "membership")
# 
# membership_list <-  lapply(unique_grps, function(target){
#   df <- do.call(rbind, lapply(seq_along(membership_list), function(i) {
#     vec <- membership_list[[i]]
#     if (target %in% names(vec)) {
#       data.frame(
#         list_index = i,
#         allele_name = target,
#         sublist = vec[[target]],
#         stringsAsFactors = FALSE
#       )
#     }
#   }))
#   
#   df
# })

# lapply(membership_list, function(x){
#   extracted_subgraphs <- Map(
#     function(list_i, sub_i) {
#       g_tmp = result[[list_i]]$subgraphs[[sub_i]]
#       V(g_tmp)$name <- V(g_tmp)$gene_family
#       g_edge <- igraph::as_data_frame(g_tmp, what = "edges")[, c("from", "to")]
#       g_edge$weight <- 1
#       g_edge
#     },
#     x$list_index,
#     x$sublist
#   )
#   
#   g_edges_weight <- rbindlist(extracted_subgraphs)
#   
#   g_edges_weight <- g_edges_weight[, .(weight = sum(weight)), 
#                                    by = c("from","to")]
#   
#   g_combined <- graph_from_data_frame(g_edges_weight, directed = FALSE)
# }


# 

