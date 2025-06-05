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

plasmid_sum[, geno_id := tstrsplit(sample_id, ":", fill=TRUE, keep = 1)]

# 10% of genomes do not have plasmids


# add machine learning likelihood to data ---------------------------------


# devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids")
library(mlplasmids)


# Change the following object with the system location of your input file
my_path <- list.files("input_data/kpne_1695_fna", 
                      recursive = TRUE,
                      full.names = TRUE)

kpne_prediction <- plasmid_classification(path_input_file = my_path,
                                          full_output = TRUE,
                                          species = 'Klebsiella pneumoniae',
                                          min_length = 150)

kpne_prediction <- as.data.table(kpne_prediction)
colnames(kpne_prediction)[4] <- "contig_id"
kpne_prediction$Prediction = tolower(kpne_prediction$Prediction)

fwrite(kpne_prediction, paste0(outdir_dat, "/kpne_prediction.csv"))
# kpne_prediction <- fread(paste0(outdir_dat, "/kpne_prediction.csv"))


# Marjorie's plasmid contigs -----------------------------------------------
# 
# marj_plasmid_contig = list.files("C:/Users/carac/Dropbox/Vos_Lab/Spark_plasmids/4100_contigs/4100_contigs")
# 
# marj_plasmid_contig <- marj_plasmid_contig[grepl("SPARK_", marj_plasmid_contig)]

# Read json from Plasmid finder ---------------------------------------------

plasmid_pf_dir = "./input_data/plasmid_finder"
plasmid_pf_files = list.files(plasmid_pf_dir, 
                              recursive = TRUE,
                              full.names = TRUE,
                              pattern = "data.json")

cl <- makeCluster(num_cores)
registerDoParallel(cl)
# foreach loop in parallel
results_json <- foreach(i = 1:length(plasmid_pf_files), 
                        .combine = rbind, 
                        .packages = c("jsonlite", "data.table")) %dopar% {
                          
                          # import json file
                          json_data <- fromJSON(plasmid_pf_files[i])
                          
                          # extract results
                          json_data <- unlist(json_data$plasmidfinder$results)
                          
                          if(any(!grepl("No hit found", json_data))){
                            tmp_data <- json_data[!grepl("No hit found", json_data)]
                            
                            tmp_data <- as.data.table(data.frame(contig_info = names(tmp_data),
                                                                 data = tmp_data))
                            
                            tmp_data[, c("contig","replicon", "var") := tstrsplit(contig_info, 
                                                                       "\\.", keep = c(3,5, 7))]
                            
                            tmp_data$contig_info = NULL
                            
                            tmp_data[, contig := tstrsplit(contig, 
                                                           ":", keep = 1L)]
                            
                            # tmp_data[, replicon := tstrsplit(replicon, 
                            #                                ":", keep = 2L)]
                            
                            # reshape
                            tmp_data <- dcast(unique(tmp_data[var %in% c("plasmid", "accession", "identity")]), # seems to be mutliple duplications of the same motif 
                                              contig + replicon ~ var, value.var = "data")
                            
                            tmp_data <- tmp_data[, .(
                              replicons = paste(unique(plasmid), collapse = ", "),
                              accessions = paste(unique(accession), collapse = ", "),
                              identities = paste(identity, collapse = ", ")
                            ), by = contig]
                            
                            tmp_data
                          }
                        }

# Stop the parallel backend after the loop is done
stopCluster(cl)

colnames(results_json)[1] = "contig_id"

results_json$pf_pred = "plasmid"

# Add plasclass info ------------------------------------------------------

plasclass_dir = "./input_data/plasclass"
plasclass_files = list.files(plasclass_dir, 
                              recursive = TRUE,
                              full.names = TRUE,
                              pattern = ".tsv")

plasclass_dat <- foreach(i = plasclass_files, 
                         .combine = "rbind") %do% {
                           
                           # Import data
                           tmp_plsclas = fread(i)
                           
                           tmp_plsclas[, plcl_type := fcase(V2 > 0.5, 'plasmid',
                                                            V2 <= 0.5, 'chromosome',
                                                            default = NA)]
                           
                           tmp_plsclas
                         }

# tidy data
colnames(plasclass_dat)[1:2] <- c("contig_id", "plcl_prob")

plasclass_dat$plcl_prob = round(plasclass_dat$plcl_prob, digits = 5)

# Combine data ----------------------------------------------------------
all_plas_dat <- merge(plasmid_sum, kpne_prediction[,.(contig_id, Prediction, Prob_Plasmid)], 
                      all.x = TRUE, by = "contig_id")

all_plas_dat <- merge(all_plas_dat, results_json, 
                      all.x = TRUE, by = "contig_id")


all_plas_dat <- merge(all_plas_dat, plasclass_dat, 
                      all.x = TRUE, by = "contig_id")

all_plas_dat[, plasmid_count := rowSums(.SD == "plasmid", na.rm = TRUE), 
             .SDcols = c("molecule_type", "plcl_type", 
                         "pf_pred", "Prediction")]

all_plas_dat[, origin := ifelse(molecule_type == "plasmid" & 
                        (pf_pred == "plasmid" | plcl_type == "plasmid" | Prediction == "plasmid"),
                      "plasmid", "chromosome")]

# fix the one contig that has both mob-suite and pf, but not ml or plclass
all_plas_dat[origin!="plasmid" & molecule_type == "plasmid" & pf_pred == "plasmid",
             origin := "plasmid"]

fwrite(all_plas_dat,
       paste0(outdir_dat, "/kpneu_plasmids.csv"))


# all_plas_dat <- fread(paste0(outdir_dat, "/kpneu_plasmids.csv"))


