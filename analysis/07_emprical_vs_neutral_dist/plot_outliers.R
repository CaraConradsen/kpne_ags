
ag_age_S_dt <-fread(paste0(outdir_dat, "/ag_age_S_dt.csv"))


# Examine outliers --------------------------------------------------------
# get 95% quantiles for the neutral distribution from rcoal.R
neutral_summary <- fread(paste0(outdir_dat, "/neutral_summary.csv"))

# Join neutral quantiles and define outliers
ag_age_dt[neutral_summary, on = "freq",
          outlier := fifelse(
            max_ks < pi_lower | max_ks > pi_upper, 1L, 0L
          )
]



# Plot graph --------------------------------------------------------------

pdf(paste0(outdir_fig, "/emp_vs_neut_outliers.pdf"), width = 8.27, height = 5.845)

par(mar = c(4,4,1,1))
with(ag_age_dt, 
     plot(freq, max_ks, pch = 16,
          ylim = c(0, 2.5),
          cex = 0.85,
          yaxt = "n", bty="L",
          ylab = expression(pi[S]),
          xlab = "Frequency in the pangenome",
          col = ifelse(outlier==1, rgb(0.8,0.1,0.1, alpha = 0.75),
                       rgb(0.1,0.1,0.8, alpha = 0.75))))

axis(side = 2, at = seq(0,2.5,0.5),
     labels = sprintf("%1.2f", seq(0,2.5,0.5)), las = 2)

legend("top", horiz = TRUE,
       legend = c("Outliers", "Non-outlier"),
       col = c("red", "blue"),
       pch = 16,
       bty = "n")

dev.off()

# Deviation values ------------------------------------------------------
# Join neutral quantiles and define outliers
ag_age_dt[neutral_summary, on = "freq",
          D := (max_ks - pi_median)/iqr
]


# add cogs
cogs <- pangraph_anno[COG_funct_cat!="", .(
  COG_mode = names(which.max(table(COG_funct_cat)))
), by = gene_family]


ag_age_dt[cogs, on = "gene_family",
          COG_letter := COG_mode
]

ag_age_dt <- merge(ag_age_dt,COG_dt,
                   all.x = TRUE,
                   by = "COG_letter")


cog_bxplt_dt = ag_age_dt[!is.na(COG_num) & D <= 5 & D >= -5]


pdf(paste0(outdir_fig, "/cog_D.pdf"), width = 8.27, height = 11.69)

par(mar = c(7,14,1,1))                   
plot(NULL,
     bty = "L",yaxt = "n",
     ylab = "",
     xlab = "Robust deviation score",
     ylim = c(26, 1),
     xlim = c(-5, 5))

with(COG_dt,
     axis(side = 2, at = COG_num, cex.axis = 0.5,
          col = "grey30",
          labels = paste0(COG_letter,", ",COG_function), las = 2))

abline(v = 0, col = "grey30", lwd = 1)
abline(v = 1.5, col = "grey30", lty = 2)
abline(v = -1.5, col = "grey30", lty = 2)

for (i in 1:26) {
  if(nrow(cog_bxplt_dt[COG_num == i]) >= 10){
    with(cog_bxplt_dt[COG_num == i],
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
  with(cog_bxplt_dt[COG_num == i],
       points(D, COG_num, factor = 0.5, pch = 16, cex = 0.65,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  if(nrow(cog_bxplt_dt[COG_num == i]) >= 10){
    with(cog_bxplt_dt[COG_num == i],
         points(median(D), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}

# fewer segregating sites than expected → candidate positive selection

# more segregating sites than expected → candidate balancing selection

arrows(x0 = -5, y0 = 29.75,      # start of arrow
       x1 = -1.5, y1 = 29.75,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey20",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(-3.25, 30.25, cex = 0.95, col = "grey30",
     "positive selection candidate", xpd = TRUE)

arrows(x0 = 5, y0 = 29.75,      # start of arrow
       x1 = 1.5, y1 = 29.75,    # end of arrow
       code = 1,             # arrows at both ends
       angle = 20,           # arrowhead angle
       length = 0.1,         # arrowhead length
       col = "grey30",
       lwd = 2,
       xpd = NA)             # allow drawing outside plot region

text(3.25, 30.25, cex = 0.95, col = "grey30",
     "balancing selection candidate", xpd = TRUE)

dev.off()

