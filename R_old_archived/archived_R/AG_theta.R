tree <- read.tree("./input_data/bootstrapped_pangraph_gubbins/boot_gub_pang_100.final_tree.tre")

ntip     <- Ntip(tree)
internal <- (ntip + 1):(ntip + tree$Nnode)   

desc      <- Descendants(tree, internal, type = "tips")  
size      <- lengths(desc)
clades    <- lapply(desc, function(i) tree$tip.label[i])

# group clades by their size
clades_by_size <- split(clades, size)


ag_seg_sites_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

plot_dat <- ag_seg_sites_dt[,.(theta_W = mean(thetaW_s)), by = c("gene_family", "number_genomes")]


set_1_names <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

with(plot_dat[gene_family %chin% set_1_names],
     boxplot(theta_W ~ number_genomes,
             col = "lightgrey", pch = 16, cex  = 0.6))

abline(h=0.02, col = "red")


