# Identify and separate out plasmids


# mob_suite predictions ---------------------------------------------------
mob_contig_summary_dir <- list.files("input_data/mob_recon_results",
                                     recursive = TRUE,
                                     full.names = TRUE)
mob_contig_summary_dir <- mob_contig_summary_dir[grepl("contig_report",
                                                       mob_contig_summary_dir)]

cl <- makeCluster(num_cores)
registerDoParallel(cl)

mob_sum <- foreach(i = 1:length(mob_contig_summary_dir),
                   .combine = rbind,
                   .packages = "data.table") %dopar% {
                     
                     fread(mob_contig_summary_dir[i], 
                           select = c("sample_id", "primary_cluster_id",
                             "contig_id", "molecule_type"))
                   } 

stopCluster(cl)

mob_sum[primary_cluster_id!="-", sample_id := paste(sample_id, 
                                                    primary_cluster_id, sep = ":"), by=.I]

mob_sum$primary_cluster_id = NULL


# add in the mod_recon info
mob_mobtyper_dir <- list.files("input_data/mob_recon_results",
                                     recursive = TRUE,
                                     full.names = TRUE)
mob_mobtyper_dir <- mob_mobtyper_dir[grepl("mobtyper_results",
                                           mob_mobtyper_dir)]

cl <- makeCluster(num_cores)
registerDoParallel(cl)

mob_typer <- foreach(i = 1:length(mob_mobtyper_dir),
                   .combine = rbind,
                   .packages = "data.table") %dopar% {
                     
                     fread(mob_mobtyper_dir[i], 
                           select = c("sample_id", "predicted_mobility",
                                      "mash_nearest_neighbor", "mash_neighbor_distance",
                                      "mash_neighbor_identification"))
                   } 

stopCluster(cl)

plasmid_sum <- merge(mob_sum, mob_typer, all.x = TRUE, by="sample_id")




# add machine learning likelihood to data ---------------------------------


# devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids")
library(mlplasmids)


# Change the following object with the system location of your input file
my_path <- list.files("input_data/kpne_1695_fna", 
                      recursive = TRUE,
                      full.names = TRUE)

kpne_prediction <- plasmid_classification(path_input_file = my_path,
                                             full_output = TRUE,
                                             species = 'Klebsiella pneumoniae')

kpne_prediction <- as.data.table(kpne_prediction)

fwrite(kpne_prediction, paste0(outdir_dat, "/kpne_prediction.csv"))
# kpne_prediction <- fread(paste0(outdir_dat, "/kpne_prediction.csv"))

# Quick and dirty ----------------------------------------------------------

temp_dat <- kpne_prediction

temp_dat[, contig_origin := fcase(Prediction=="Plasmid" & Prob_Plasmid <0.9,
                                  "Chromosome",
                                  Prediction=="Plasmid" & Prob_Plasmid >= 0.9,"Plasmid",
                                  default = "Chromosome")
         ]


fwrite(temp_dat[,.(Contig_name, contig_origin)], 
       paste0(outdir_dat, "/kpneu_plasmids.csv"))
