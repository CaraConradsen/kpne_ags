# Retain only correct oriented chromosomes --------------------------------
pan_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_cogs.csv"))[,c(1:11, 13,17)]


# subset anno by pangraph anno
pangraph_anno <- fread(paste0(outdir_dat, "/pangraph_anno.csv"),
                       select = c("locus_tag", "ag_type"))# 214 genomes


pan_anno <- merge(pangraph_anno, pan_anno,
                  all.x = TRUE, by ="locus_tag")

# add syntenic identifiers
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                              select = c("locus_tag", "anchor",
                                         "msu", "acrs_msu","acrs_jun"))# accounts for 99.68% loci 

pan_anno <- merge(pan_anno,msu_regions_anchored,
                  all.x = TRUE, by ="locus_tag")

# tidy binary cats
pan_anno[is.na(acrs_msu), acrs_msu := NA]
pan_anno[is.na(acrs_jun), acrs_jun := NA]

# add msu order
msu_paths <- fread(paste0(outdir_dat, "/msu_paths_dt.csv"), 
                   select = c("geno_id", "order", "msu_mergers"))

# use ref order
msu_paths <- msu_paths[geno_id == "SPARK_1004_C1",]
colnames(msu_paths)[2:3] = c("msu_order", "msu")

pan_anno <- merge(pan_anno,msu_paths[,.(msu_order, msu)],
                  all.x = TRUE, by = "msu")

# add msu length
MSU_len <- yyjsonr::read_json_file("./input_data/pangraph/MSU_len.json",
                                   int64 = "string")

msu_len_dt <- rbindlist(lapply(MSU_len, function(x) {
  data.table(
    msu_length = x
  )
}), fill = TRUE)

msu_len_dt[, msu := names(MSU_len)]

pan_anno <- merge(pan_anno,msu_len_dt,
                  all.x = TRUE, by = "msu")



# Explore regions of high hgt ---------------------------------------------

multi_hgt <- unique(pan_anno[acrs_msu==1 | acrs_jun==1,.(gene_family, msu, anchor,
                                                  acrs_msu, acrs_jun, msu_order, msu_length)])

multi_hgt[, num_hgt := .N, by = c("msu")]

multi_hgt[, rel_num:= num_hgt/msu_length]
multi_hgt[, rel_num:= rel_num*1000]

rel_counts = unique(multi_hgt[,.(rel_num, msu, msu_order)])

temp <- msu_paths[!msu %chin% rel_counts$msu,]
temp[,rel_num := 0]

rel_counts = rbind(rel_counts, temp[,.(rel_num, msu,msu_order)])

setorderv(rel_counts, cols="msu_order", order = 1L)

nodes_dt <- fread(paste0(outdir_dat, "/nodes_dt.csv"), 
                  select = c("block_id", "start"))
msu_mergers_dt <- fread(paste0(outdir_dat, "msu_mergers_dt.csv")) # core block ids


msu_bp <- merge(nodes_dt, msu_mergers_dt, by = "block_id")

colnames(msu_bp)[3] = "msu"

msu_bp <- rbind(
  msu_bp[msu != "MSU_6",.(bp = median(start)), by = msu],
  msu_bp[msu == "MSU_6" & start <18552,.(bp = median(start)), by = msu])

rel_counts <- merge(rel_counts, msu_bp,
                    all.x = TRUE, by ="msu")

setorderv(rel_counts, cols="bp", order = 1L)

with(rel_counts,
     plot(bp/1e03, rel_num, type = "b",
          bty = "L",yaxt = "n",
          xlab = "Genome position (Kbp)",
          ylab = "Relative HGT enrichment (x1000)",
          pch = 16, col = "dodgerblue"))
axis(side = 2, at = seq(0,8,2),
     labels = seq(0,8,2), las = 2)


# COGs --------------------------------------------------------------------
multi_cogs <- unique(pan_anno[!is.na(COG_funct_cat) & gene_family %chin% multi_hgt$gene_family,
                              .(gene_family, COG_funct_cat)])

multi_cogs <- multi_cogs[COG_funct_cat!=""]

multi_cogs <- multi_cogs[, (n=.N), COG_funct_cat]
