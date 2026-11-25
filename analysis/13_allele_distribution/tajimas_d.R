library(ape); library(pegas)

focal_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

ag_fastas <- paste0("./output/data/ag_fus_alns/",
                    focal_ags, ".fasta")

cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

ag_tajima_stats <- foreach::foreach(ag = ag_fastas,
                                    .combine = "rbind",
                                    .packages = c("data.table", "ape", "pegas")) %dopar% {
                                      gene_family = gsub(".fasta", "", basename(ag))
                                      aln <- read.dna(ag,
                                                      format = "fasta", as.matrix = TRUE)   
                                      data.frame(gene_family, tajima.test(aln))
                                    }

stopCluster(cl)

fwrite(ag_tajima_stats, paste0(outdir_dat, "ag_tajima_stats.csv"))

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))
ag_tajima_stats <- fread(paste0(outdir_dat, "ag_tajima_stats.csv"))
trouble_children <- fread(paste0(outdir_dat, "/trouble_children.csv"))
trouble_children <- merge(trouble_children,
                          gene_pvals, 
                          by = "gene_family", all.x = TRUE)

outliers_t <- merge(trouble_children,
                    ag_tajima_stats, 
                    by = "gene_family", all.x = TRUE)
