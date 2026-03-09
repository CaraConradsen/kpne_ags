
# Import ASR and CI -------------------------------------------------------
ci <- fread(paste0(outdir_dat,"/consistencyindex_syn.csv"))
gl <- fread(paste0(outdir_dat,"/maximum_parsimony_gain_loss_results.csv"))
colnames(gl)[3:5] <- c("par_gains","par_losses","par_anc_state")
# ace <- fread(paste0(outdir_dat, "/gain_loss_results_syn.csv"))
# colnames(ace)[2:4] <- c("ace_gains","ace_losses","ace_anc_state")
# # import boot and get summary
# gain_boot <- fread(paste0(outdir_dat, "/gain_loss_boot_results.csv"))
# 
# gain_boot_sum <- gain_boot[, .(
#   median_gain = as.numeric(median(gains, na.rm = TRUE)),
#   gain_ci_lo = as.numeric(quantile(gains, 0.025, na.rm = TRUE)),
#   gain_ci_up = as.numeric(quantile(gains, 0.975, na.rm = TRUE)),
#   median_loss = as.numeric(median(losses, na.rm = TRUE)),
#   loss_ci_lo = as.numeric(quantile(losses, 0.025, na.rm = TRUE)),
#   loss_ci_up = as.numeric(quantile(losses, 0.975, na.rm = TRUE)),
#   reps = max(rep)
# ), by = gene_family]

# get ags
msu_regions_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"))
ag_homo_dat <- msu_regions_anchored[ag_type != "core"][number_genomes != 1][!grepl("_",gene_family)]# [acrs_msu!=1][acrs_jun != 1]

# remove anchor dups
ag_homo_dat <- unique(ag_homo_dat[,.(gene_family, number_genomes, anchor, acrs_msu, acrs_jun)])

nonsyntenic_loci <- unique(ag_homo_dat[acrs_msu == 1 | acrs_jun == 1 , gene_family])

syntenic_key <- ag_homo_dat[!gene_family %chin% nonsyntenic_loci, .(
  syn_jun = ifelse(all(anchor != ""), 1, 0)
), by = gene_family]

# # Key
# 0 = at_least_one_unknown
# 1 = all_known
# 2 = non-syntenic

# put data.frames together

phylo_info = merge(ci, gl, by = "gene_family")

# phylo_info = merge(phylo_info, ace, by = "gene_family")

# phylo_info = merge(phylo_info, gain_boot_sum, 
#                    all.x = TRUE, by = "gene_family")

# add gene frequency
phylo_info = merge(phylo_info, unique(ag_homo_dat[,.(gene_family, number_genomes)]),
                   all.x = TRUE, by = "gene_family")

# add synteny info
phylo_info = merge(phylo_info, syntenic_key,
                   all.x = TRUE, by = "gene_family")

# fix nonsytenic
phylo_info[is.na(syn_jun), syn_jun := 2]

# Retain genes for set 2

# fwrite(phylo_info, paste0(outdir_dat, "/phylo_info.csv"))
# fwrite(phylo_info[syn_jun!=2, .(gene_family)], paste0(outdir_dat, "/set_1_names.csv"))

# Plot homoplasy vs parsimony ---------------------------------------------

png(paste0(outdir_fig,"/ci_vs_par_gains.png"),
    width = 16, height = 11, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

mat <- matrix(c(1:6), nrow = 2, byrow = TRUE)
layout(mat, widths = c(1.25,1,1))

#syntenic

col_map <- c(
  "0" = rgb(0.1,0.1, 0.8, alpha = 0.35),
  "1" = rgb(0.8,0.1, 0.1, alpha = 0.35),
  "2" = rgb(0.1, 0.1, 0.1, alpha = 0.35)
)

gene_set <- c(
  "0" = "Unknown synteny",
  "1" = "Within syntenic block",
  "2" = "Nonsyntenic"
)

for (i in c(2,1,0)) {
  colr <- col_map[as.character(i)]
  
  if(i==2){
    par(mar = c(0, 4,4,0))
  }else{
    par(mar = c(0, 0,4,0))
  }
  
  with(phylo_info[syn_jun==i],
       plot(number_genomes, consistencyindex,
            pch = 19,
            main = "",
            cex = 0.5, yaxt ="n",
            col = colr,
            xaxt = "n",
            xlab = "", 
            ylab = ifelse(i==2,"Consistency index", "")))
  
  if(i==2){
    axis(side=2, at = seq(0,1, 0.2),
         las = 2,
         labels = sprintf("%.1f", seq(0,1, 0.2)))
  }
  
  axis(side = 3, at = 125, tick = FALSE,
       labels = bquote(
         .(gene_set[as.character(i)])~" (" * italic(n) ~ "=" ~ .(nrow(phylo_info[syn_jun == i])) * ")"
       ))
  
}


for (i in c(2,1,0)) {
  colr <- col_map[as.character(i)]
  
  if(i==2){
    par(mar = c(4, 4,0,0))
  }else{
    par(mar = c(4, 0,0,0))
  }
  
  with(phylo_info[syn_jun==i],
       plot(number_genomes, par_gains,
            pch = 19, #bty = "L",
            cex = 0.5, yaxt ="n",
            col = colr,
            xlab = ifelse(syn_jun==1, "Number of genomes",""), 
            ylab = "Inferred parsimony gains"))
  
  if(i==2){
    axis(side=2, at = seq(0,40, 10),
         las = 2,
         labels = seq(0,40, 10))
  }
  
}

dev.off()



# Set 2 strict equal or less than 1-------------------------------------------------------

set2_ags <- rbind(phylo_info[syn_jun != 2 & par_gains == 1 & par_anc_state == 0],
                  phylo_info[syn_jun != 2 & par_gains == 0 & par_anc_state == 1])


# fwrite(set2_ags[, .(gene_family)], paste0(outdir_dat, "/set_2_names.csv"))

# png(paste0(outdir_fig,"/ci_vs_par_gains.png"),
#     width = 16, height = 11, units = "cm", res = 300,
#     pointsize = 14, type = "cairo")

mat <- matrix(c(1:2), nrow = 1, byrow = TRUE)
layout(mat, widths = c(1.25,1))

#syntenic

col_map <- c(
  "0" = rgb(0.1,0.1, 0.8, alpha = 0.35),
  "1" = rgb(0.8,0.1, 0.1, alpha = 0.35)
)

gene_set <- c(
  "0" = "Unknown synteny",
  "1" = "Within syntenic block"
)

for (i in c(1,0)) {
  colr <- col_map[as.character(i)]
  
  if(i==1){
    par(mar = c(4, 4,4,0))
  }else{
    par(mar = c(4, 0,4,0))
  }
  
  with(set2_ags[syn_jun==i],
       plot(number_genomes, consistencyindex,
            pch = 19,
            main = "",
            cex = 0.5, yaxt ="n",
            col = colr,
            xlab = ifelse(syn_jun==1, "Number of genomes",""), 
            ylab = ifelse(i==1,"Consistency index", "")))
  
  if(i==1){
    axis(side=2, at = seq(0,1, 0.2),
         las = 2,
         labels = sprintf("%.1f", seq(0,1, 0.2)))
  }
  
  axis(side = 3, at = 125, tick = FALSE,
       labels = bquote(
         .(gene_set[as.character(i)])~" (" * italic(n) ~ "=" ~ .(nrow(set2_ags[syn_jun == i])) * ")"
       ))
  
}


