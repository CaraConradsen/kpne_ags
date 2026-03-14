
pangraph_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                                  "number_genomes", "start","end",
                                  "strand","ST","ag_type", "egng_cog", 
                                  "consensus_gene_name", "consensus_product"))

# rename eggnog_cog to COG_funct_cat
setnames(pangraph_anno, old = "egng_cog", new = "COG_funct_cat")

ag_age_S_dt <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))


# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
neutral_summary <- fread(paste0(outdir_dat, "/neutral_summary.csv"))

# estimates for set 1 and 2
set_1_ags = fread(paste0(outdir_dat, "/set_1_names.csv"))
set_2_ags = fread(paste0(outdir_dat, "/set_2_names.csv"))

# Anchor genealogy estimates by empirical means ---------------------------
ag_age_S_dt_set1 <- ag_age_S_dt[gene_family %chin% set_1_ags$gene_family]
ag_age_S_dt_set2 <- ag_age_S_dt[gene_family %chin% set_2_ags$gene_family]


# Join neutral quantiles and define outliers
ag_age_S_dt_set1[neutral_summary[set==1], on = "freq",
          outlier := fifelse(
            mean_ks < pi_lower | mean_ks > pi_upper, 1L, 0L
          )
]

ag_age_S_dt_set2[neutral_summary[set==2], on = "freq",
                 outlier := fifelse(
                   mean_ks < pi_lower | mean_ks > pi_upper, 1L, 0L
                 )
]

plot_dat_list <- list(ag_age_S_dt_set1, ag_age_S_dt_set2)

# Plot graph --------------------------------------------------------------
# nice light blue: 0.400, 0.714, 0.922

png(paste0(outdir_fig,"/emp_vs_neut_outliers.png"),
    width = 32, height = 16, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

par(mfrow = c(1,2), mar = c(2,3,1,0.5),
    oma = c(1.5,1,0.5,0))

for (i in seq_along(plot_dat_list)) {
    plot(NULL, pch = 16,
       xlim = c(0, 260),  
       ylim = c(0, 1),
       cex = 0.45,
       yaxt = "n", bty="L",
       ylab = "",
       xlab = "")
  
  line_vals <- neutral_summary[set==i][order(freq)]
  
  polygon(
    c(line_vals$freq, rev(line_vals$freq)),
    c(line_vals$pi_upper, rev(line_vals$pi_lower)),
    col = "grey80",
    border = "grey80"
  )
  
  with(plot_dat_list[[i]], 
       points(freq, mean_ks, pch = 16,
              ylim = c(0, 1),
              cex = 0.45,
              yaxt = "n", bty="L",
              ylab = "",
              xlab = "",
              col = ifelse(outlier==1, rgb(0.2, 0.6, 1, alpha = 0.5),
                           rgb(0.1,0.1,0.1, alpha = 0.5))))

  axis(side = 2, at = seq(0,1,0.25),
       labels = sprintf("%1.2f", seq(0,1,0.25)), las = 2)
  
  mtext(paste0("Set ", i), side = 3)
}

mtext(expression(pi[S]), side = 2, outer = TRUE)
mtext("Frequency in the pangenome", side = 1, outer = TRUE)

legend("topleft", #horiz = TRUE,
       legend = c("Outliers", "Non-outlier", "Neutral 95% CI"),
       col = c(rgb(0.2, 0.6, 1, alpha = 1), "black", "grey80"),
       pch = c(16, 16,15),
       bty = "n")

dev.off()

# Deviation values ------------------------------------------------------
# Join neutral quantiles and define outliers
ag_age_S_dt_set1[neutral_summary[set==1], on = "freq",
          D := (mean_ks - pi_median)/iqr
]

ag_age_S_dt_set2[neutral_summary[set==2], on = "freq",
                 D := (mean_ks - pi_median)/iqr
]

# add cogs
ag_age_S_dt_set1[pangraph_anno[,.(gene_family, COG_funct_cat)], on = "gene_family",
          COG := COG_funct_cat
]

ag_age_S_dt_set2[pangraph_anno[,.(gene_family, COG_funct_cat)], on = "gene_family",
                 COG := COG_funct_cat
]

set1_cogs <- ag_age_S_dt_set1[outlier==1 & COG != "" & mean_ks <= 1, .(D, COG, gene_family)]
set2_cogs <- ag_age_S_dt_set2[outlier==1 & COG != "" & mean_ks <= 1, .(D, COG, gene_family)] 

nrow(set1_cogs); nrow(set2_cogs)

set1_cogs <- set1_cogs[, .(COG_letter = unlist(strsplit(COG, ""))), by = c("gene_family","D")]
set2_cogs <- set2_cogs[, .(COG_letter = unlist(strsplit(COG, ""))), by = c("gene_family","D")]

# add cog annotation from COG_dt
set1_cogs <- merge(set1_cogs, COG_dt,
                   all.x = TRUE, by = "COG_letter")
set2_cogs <- merge(set2_cogs, COG_dt,
                   all.x = TRUE, by = "COG_letter")


# plot set 1 --------------------------------------------------------------
# 
png(paste0(outdir_fig,"/cog_D_set_1.png"),
    width = 32, height = 16, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

par(mar = c(7,17,1,1))                   
plot(NULL,
     bty = "L",yaxt = "n",
     ylab = "",
     xlab = "Robust deviation score",
     ylim = c(26, 1),
     xlim = c(-10, 30))

with(set1_cogs,
     axis(side = 2, at = broad_ord, cex.axis = 0.5,
          col = "grey30",
          labels = paste0(COG_letter,", ",COG_function), las = 2))

abline(v = 0, col = "grey30", lwd = 1)
abline(v = 1.5, col = "grey30", lty = 2)
abline(v = -1.5, col = "grey30", lty = 2)

for (i in 1:26) {
  if(nrow(set1_cogs[broad_ord == i]) >= 5){
    with(set1_cogs[broad_ord == i],
         vioplot(D, at = i, 
                 col = "salmon1",
                 border = NA,
                 horizontal = TRUE,
                 drawRect = FALSE,
                 wex = 0.75,
                 add = TRUE))
  }
}

for (i in 1:26) {
  with(set1_cogs[broad_ord == i],
       points(D, broad_ord, pch = 16, cex = 0.65,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  if(nrow(set1_cogs[broad_ord == i]) >= 10){
    with(set1_cogs[broad_ord == i],
         points(median(D), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}

# add broad categories
for (i in unique(set1_cogs$desc)) {
  if(i != ""){
    axis_grp_dat <- set1_cogs[desc==i]
    
    n_cogs = length(unique(axis_grp_dat$COG_letter))
    
    axis_grp_dat <- axis_grp_dat[, .(min = min(broad_ord),
                                     max = max(broad_ord))]
    axis(side = 4,
         line = -50,
         at = c(axis_grp_dat$min,
                axis_grp_dat$max),
         labels = rep("", 2))
    
    axis(side = 4,
         line = -52.5,
         tick = FALSE,
         cex.axis = 0.5,
         at = mean(c(axis_grp_dat$min,
                     axis_grp_dat$max)),
         labels = paste(strwrap(i, width = 20), 
                        collapse = "\n"))
  }
}

legend("topright", legend = expression(bold("Set 1")),
       bty = "n", cex = 1.1)

# fewer segregating sites than expected → candidate positive selection

# more segregating sites than expected → candidate balancing selection

arrows(x0 = -10, y0 = 34,      # start of arrow
       x1 = -1.5, y1 = 34,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey20",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(-5, 35.5, cex = 0.95, col = "grey30",
     "positive selection candidate", xpd = TRUE)

arrows(x0 = 30, y0 = 34,      # start of arrow
       x1 = 1.5, y1 = 34,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey30",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(15, 35.5, cex = 0.95, col = "grey30",
     "balancing selection candidate", xpd = TRUE)

dev.off()


# plot set 2 --------------------------------------------------------------
# 
png(paste0(outdir_fig,"/cog_D_set_2.png"),
    width = 32, height = 16, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

par(mar = c(7,17,1,1))                   
plot(NULL,
     bty = "L",yaxt = "n",
     ylab = "",
     xlab = "Robust deviation score",
     ylim = c(26, 1),
     xlim = c(-10, 10))

with(set2_cogs,
     axis(side = 2, at = broad_ord, cex.axis = 0.5,
          col = "grey30",
          labels = paste0(COG_letter,", ",COG_function), las = 2))

abline(v = 0, col = "grey30", lwd = 1)
abline(v = 1.5, col = "grey30", lty = 2)
abline(v = -1.5, col = "grey30", lty = 2)

for (i in 1:26) {
  if(nrow(set2_cogs[broad_ord == i]) >= 5){
    with(set2_cogs[broad_ord == i],
         vioplot(D, at = i, 
                 col = "salmon1",
                 border = NA,
                 horizontal = TRUE,
                 drawRect = FALSE,
                 wex = 0.75,
                 add = TRUE))
  }
}

for (i in 1:26) {
  with(set2_cogs[broad_ord == i],
       points(D, broad_ord, pch = 16, cex = 0.65,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  if(nrow(set2_cogs[broad_ord == i]) >= 10){
    with(set2_cogs[broad_ord == i],
         points(median(D), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}

# add broad categories
for (i in unique(set2_cogs$desc)) {
  if(i != ""){
    axis_grp_dat <- set2_cogs[desc==i]
    
    n_cogs = length(unique(axis_grp_dat$COG_letter))
    
    axis_grp_dat <- axis_grp_dat[, .(min = min(broad_ord),
                                     max = max(broad_ord))]
    axis(side = 4,
         line = -50,
         at = c(axis_grp_dat$min,
                axis_grp_dat$max),
         labels = rep("", 2))
    
    axis(side = 4,
         line = -52.5,
         tick = FALSE,
         cex.axis = 0.5,
         at = mean(c(axis_grp_dat$min,
                     axis_grp_dat$max)),
         labels = paste(strwrap(i, width = 20), 
                        collapse = "\n"))
  }
}

legend("topright", legend = expression(bold("Set 2")),
       bty = "n", cex = 1.1)

# fewer segregating sites than expected → candidate positive selection

# more segregating sites than expected → candidate balancing selection

arrows(x0 = -10, y0 = 34,      # start of arrow
       x1 = -1.5, y1 = 34,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey20",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(-5, 35.5, cex = 0.95, col = "grey30",
     "positive selection candidate", xpd = TRUE)

arrows(x0 = 10, y0 = 34,      # start of arrow
       x1 = 1.5, y1 = 34,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey30",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(5, 35.5, cex = 0.95, col = "grey30",
     "balancing selection candidate", xpd = TRUE)

# label extremes

extreme_labs <- set2_cogs[, .SD[which.max(abs(D))], by = COG_letter][,.(gene_family, D, broad_ord)]

extreme_labs <- merge(extreme_labs, unique(pangraph_anno[,.(gene_family, consensus_gene_name, consensus_product)]),
                      all.x = TRUE, by = "gene_family")

fwrite(extreme_labs, paste0(outdir_dat, "/extreme_labs.csv"))

with(extreme_labs,
     text(D, broad_ord,
          labels = consensus_gene_name, 
          cex = 0.45,
          pos = ifelse(D>0,4,2),
          offset = 0.25))

dev.off()



