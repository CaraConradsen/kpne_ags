# Analyse the output from Goldfinder


# Part 1 Simultaneous estimates of association ----------------------------

ass_clusters <- read.delim("./input_data/goldfinder/simul_output/association_clusters.txt",
                           header = FALSE)

header_idx <- which(startsWith(ass_clusters$V1, ">"))

# end of each cluster = one before the next header; last cluster runs to end of df
starts <- header_idx + 1
ends   <- c(header_idx[-1] - 1, nrow(ass_clusters))

# list of row indices per cluster
member_rows <- Map(seq, starts, ends)

a_clus_dt <- lapply(seq_along(member_rows), function(g) {
  assoc_cluster = gsub(">","", ass_clusters$V1[header_idx[g]])
  gene_family = ass_clusters$V1[member_rows[[g]]]
  
  data.table(cbind(assoc_cluster,gene_family))
})

a_clus_dt <- rbindlist(a_clus_dt)

a_clus_dt[,c("assoc_cluster", "n_genes_assoc_cluster"):= tstrsplit(assoc_cluster, ",", fill = TRUE)]
a_clus_dt[, c("assoc_cluster", "n_genes_assoc_cluster") := lapply(.SD, as.numeric), .SDcols = c("assoc_cluster", "n_genes_assoc_cluster")]
a_clus_dt[, gene_family:= gsub(",", "", gene_family)]
a_clus_dt[, gene_family:= tstrsplit(gene_family, "__", fill = TRUE, keep = 2)]

a_clus_dt$assc_type = "associate"

# assign whether dis
dis_pairs <- fread("./input_data/goldfinder/simul_output/simultaneous_dissociation_significant_pairs.csv")
dis_pairs[, Gene_1:= tstrsplit(Gene_1, "__", fill = TRUE, keep = 2)]
dis_pairs[, Gene_2:= tstrsplit(Gene_2, "__", fill = TRUE, keep = 2)]

all_clus_dt <- rbind(a_clus_dt, data.table(assoc_cluster = max(a_clus_dt$assoc_cluster)+1,
                                           gene_family = c(dis_pairs$Gene_1, dis_pairs$Gene_2),
                                           n_genes_assoc_cluster = 2,
                                           assc_type = "disassociate"))


all_clus_dt$score = "simultaneous"

fwrite(all_clus_dt, paste0(outdir_dat, "/simul_clust_fams.csv"))


# Part 2 Subsequent estimates of association ----------------------------
ass_sub_clusters <- read.delim("./input_data/goldfinder/sub_output/association_clusters.txt",
                           header = FALSE)

header_idx <- which(startsWith(ass_sub_clusters$V1, ">"))

# end of each cluster = one before the next header; last cluster runs to end of df
starts <- header_idx + 1
ends   <- c(header_idx[-1] - 1, nrow(ass_sub_clusters))

# list of row indices per cluster
member_rows <- Map(seq, starts, ends)

a_clus_sub_dt <- lapply(seq_along(member_rows), function(g) {
  assoc_cluster = gsub(">","", ass_sub_clusters$V1[header_idx[g]])
  gene_family = ass_sub_clusters$V1[member_rows[[g]]]
  
  data.table(cbind(assoc_cluster,gene_family))
})

a_clus_sub_dt <- rbindlist(a_clus_sub_dt)

a_clus_sub_dt[,c("assoc_cluster", "n_genes_assoc_cluster"):= tstrsplit(assoc_cluster, ",", fill = TRUE)]
a_clus_sub_dt[, c("assoc_cluster", "n_genes_assoc_cluster") := lapply(.SD, as.numeric), .SDcols = c("assoc_cluster", "n_genes_assoc_cluster")]
a_clus_sub_dt[, gene_family:= gsub(",", "", gene_family)]
a_clus_sub_dt[, gene_family:= tstrsplit(gene_family, "__", fill = TRUE, keep = 2)]

a_clus_sub_dt$assc_type = "associate"

a_clus_sub_dt$score = "subsequent"

fwrite(a_clus_sub_dt, paste0(outdir_dat, "/sub_clust_fams.csv"))



# combine networks and gene outliers ----------------------------------------------
coin_ag <- rbind(all_clus_dt, a_clus_sub_dt)

# create parameter for ags that didn't pass set 1 filter
set1_ags <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/set_1_names.csv")$gene_family

coin_ag[, n_in_set1 := uniqueN(gene_family[gene_family %chin% set1_ags]),
        by = .(assc_type, score, assoc_cluster)]

coin_ag[, perc_in_set1 := round((n_in_set1/n_genes_assoc_cluster)*100, digits = 0)]

# Find gene outlier associations 
ag_outliers <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]

coin_ag <- merge(ag_outliers, coin_ag,
                     all.x = TRUE, by = "gene_family")

associated_outliers <- coin_ag[!is.na(assoc_cluster)][order(assoc_cluster)]

associated_outliers[, c("expected_S", "p",
                        "cdf","p_lower","p_upper") := NULL]

setnames(associated_outliers, "freq", "number_genomes")


# Summary of findings -----------------------------------------------------

sum_dat <-associated_outliers[p_bh <= 0.05, .(gene_family, score)]

sum_dat <- sum_dat[, .(score = paste(sort(unique(score)), collapse = ":")), by = .(gene_family)]

sum_dat[,.(n = .N), by = c("score")]

# w_theta_sel                   score     n
# 1: directional            simultaneous     1
# 2:   balancing            simultaneous    25
# 3:   balancing simultaneous:subsequent    10
# 4:   balancing              subsequent     6

sum_dat <-associated_outliers[p_bh <= 0.05 & perc_in_set1 == 100, .(gene_family, score)]

sum_dat <- sum_dat[, .(score = paste(sort(unique(score)), collapse = ":")), by = .(gene_family)]

sum_dat[,.(n = .N), by = c("score")]

# w_theta_sel                   score     n
# 1: directional            simultaneous     1
# 2:   balancing            simultaneous    10
# 3:   balancing simultaneous:subsequent    6
# 4:   balancing              subsequent    8

# Clean data for network image --------------------------------------------
# take set 2 over set 1 (only 1 dw); simultaneous over subsequent

clean_coin_ags <- associated_outliers[score=="simultaneous"]

# new clus id
clean_coin_ags$clust_id = clean_coin_ags$assoc_cluster

max_clus = max(clean_coin_ags$assoc_cluster)

subs <-associated_outliers[score!="simultaneous"]

sort_sim_clust <- clean_coin_ags$assoc_cluster

#loop over and re-cluster sim and subs
for (cl in sort_sim_clust) {
  temp_dat <- clean_coin_ags[assoc_cluster == cl]
  
  temp_subs <- rbind(subs[gene_family %chin% temp_dat$gene_family], NULL)
  
  if(!is.null(temp_subs)){
    temp_subs$clust_id <- unique(temp_dat$clust_id)
    
    clean_coin_ags <- rbind(clean_coin_ags, temp_subs)
    
    subs <- subs[!gene_family %chin% temp_subs$gene_family]
  }

}

# fix remaining subs
unique_clusters <- subs[, unique(assoc_cluster)]
cluster_map <- data.table(
  assoc_cluster = unique_clusters,
  clust_id = seq(max_clus + 1, max_clus + length(unique_clusters))
)
subs <- cluster_map[subs, on = "assoc_cluster"]


setcolorder(subs, neworder = colnames(clean_coin_ags))

# put everything together
clean_coin_ags <- rbind(clean_coin_ags, subs)

fwrite(clean_coin_ags, paste0(outdir_dat,"/clean_coin_ags.csv"))

clean_coin_ags <- fread(paste0(outdir_dat,"/clean_coin_ags.csv"))

# Fisher's exact test on clusters in Set 1 ------------------------------------
pan_anno <- unique(fread("./output/data/pangenome_anno.csv", 
                  select = c("gene_family", "number_genomes")))

# remove genes less than 5%
pan_anno <- pan_anno[number_genomes >= (260*0.05)]

# add co-occurence & genealogically ordered cluster id
coincident_genes <- unique(rbind(all_clus_dt, a_clus_sub_dt)$gene_family)

pan_anno <- pan_anno[, gf_hits := fcase(gene_family %chin% coincident_genes, "Y",
                                        default = "N")]

# add selection group info 
selected_ags <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05, gene_family]

pan_anno <- pan_anno[, sel_grp := fcase(gene_family %chin% selected_ags, "sel", 
                                        default = "unsel")]

# get set info
set1_ags <- fread("./output/data/set_1_names.csv")$gene_family

# retain set 1 genes
pan_anno_set1 <- pan_anno[gene_family %chin% set1_ags]

gene_sets <- list(pan_anno_set1)

fexact_coincident_genes  <- lapply(seq_along(gene_sets), function(lst_num){
  
  ag_lst <- gene_sets[[lst_num]]
  
  tab <- table(ag_lst[,.(gf_hits,sel_grp)])
  
  res <- fisher.test(tab)
  
  data.table(
    set      = lst_num,
    pval     = round(res$p.value, 4),
    OR       = res$estimate,
    CI_lower = round(res$conf.int[1], 4),
    CI_upper = round(res$conf.int[2], 4)
  )
  
})


fexact_coincident_genes <- rbindlist(fexact_coincident_genes)

fexact_coincident_genes$p_bh <- round(p.adjust(fexact_coincident_genes$pval, "BH"), 4)

print(fexact_coincident_genes)

