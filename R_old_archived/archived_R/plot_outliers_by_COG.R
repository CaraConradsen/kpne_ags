# don't double count genes!!

pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", 
                                  "fus_locus_tag","egng_cog", 
                                  "consensus_gene_name", "consensus_product")))

# rename eggnog_cog to COG_funct_cat
setnames(pangraph_anno, old = "egng_cog", new = "COG_funct_cat")

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]


# Calculate signed log10(p) -----------------------------------------------
# Join neutral quantiles and define outliers
# Poisson-Based Approaches
# Since segregating sites are count data, Poisson tests are statistically appropriate:
# Calculate the mean from neutral genes (λ)
# Use Poisson probability: P(X ≥ observed | λ)
# Convert to a score: -log10(p-value)

# signed magnitude: positive = excess, negative = deficit
gene_pvals[, direction := fifelse(seg_s_sites > expected_S, 1, -1)]
gene_pvals[, signed_logp := direction * -log10(p_two)]

# Plot by functional category -----------------------------------------
cog_set <- gene_pvals[,.(gene_family, signed_logp)]

# add cogs
cog_set[pangraph_anno[,.(gene_family, COG_funct_cat)], on = "gene_family",
        COG := COG_funct_cat
]

# How many dup COGs are there?
cat("Number of genes with dup COGS:", nrow(cog_set[nchar(COG) > 1]))
# Number of genes with dup COGS: 75 (new piS analysis)

cog_set <- rbind(cog_set[, .(COG_letter = unlist(strsplit(COG, ""))), by = c("gene_family","signed_logp")],
                 cog_set[COG=="",.(gene_family, signed_logp, COG_letter= COG)])

# add cog annotation from COG_dt
cog_set <- merge(cog_set, COG_dt,
                   all.x = TRUE, by = "COG_letter")

# cogs not represented
COG_dt[!broad_ord %chin% cog_set$broad_ord]

# tidy data
new_cols <- unique(cog_set[, .(broad_ord)])[order(broad_ord)]
new_cols[, new_y := .I]
mx_y = max(new_cols$new_y)
cog_set <- cog_set[new_cols, on = "broad_ord"]

fwrite(cog_set, paste0(outdir_dat, "/bh_p_cog_set.csv"))

# How many COGs categories are missing?
cat("Missing COGs:\n")
print(COG_dt[!COG_letter %chin% cog_set$COG_letter, .(COG_letter, COG_function, desc)])

# plot set 1 --------------------------------------------------------------
# 


png(paste0(outdir_fig,"/piS_theta_cog_signed_logp_set_1.png"),
    width = 32, height = 17, units = "cm", res = 300,
    pointsize = 14, type = "cairo")


# par(mar = c(7,17,1,1))
par(mar = c(4,17,1,1))  
plot(NULL,
     bty = "L",yaxt = "n",
     ylab = "",
     xlab = expression("Signed log"[10]~italic(p)),
     ylim = c(mx_y + 1, 0),
     xlim = c(0,12))

with(cog_set,
     axis(side = 2, at = new_y, cex.axis = 0.5,
          col = "grey30",
          labels = ifelse(COG_function == "Unassigned", 
                          COG_function,
                          paste0(COG_letter,", ",COG_function)), 
          las = 2))

abline(v = 0, col = "grey30", lwd = 1, lty = 2)


for (i in 1:mx_y) {
  if(nrow(cog_set[new_y == i]) >= 5){
    with(cog_set[new_y == i],
         vioplot(signed_logp, at = i, 
                 col = "grey80",
                 border = NA,
                 horizontal = TRUE,
                 drawRect = FALSE,
                 wex = 0.75,
                 add = TRUE))
  }
}

for (i in 1:mx_y) {
  with(cog_set[new_y == i],
       points(signed_logp, new_y, pch = 16, cex = 0.65,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  # # overlay set 2
  # with(cog_set[set==2 & new_y == i],
  #      points(signed_logp, new_y, pch = 21, cex = 0.7,
  #             col = "darkgoldenrod3",
  #             bg = "darkgoldenrod1"))
  
  if(nrow(cog_set[new_y == i]) >= 10){
    with(cog_set[new_y == i],
         points(median(signed_logp), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}

# add legend
legend("topright",c("AG outliers", "Median"),#c("Set 1", "Set 2"),
       pch = c(16,23),#,21), 
       col = c(rgb(0.1,0.1,0.1, alpha = 0.75), "black"),
               # "darkgoldenrod3"),
       pt.bg = c(NULL,"white"), 
       # pt.bg = c(NULL,"darkgoldenrod1"),
       pt.cex = c(0.65,1.2),
       horiz = TRUE,
       bty = "n", cex = 0.95)

# add broad categories
for (i in unique(cog_set$desc)) {
  if(i != ""){
    axis_grp_dat <- cog_set[desc==i]
    
    n_cogs = length(unique(axis_grp_dat$COG_letter))
    
    axis_grp_dat <- axis_grp_dat[, .(min = min(new_y),
                                     max = max(new_y))]
    axis(side = 4,
         line = -50.5,
         at = c(axis_grp_dat$min,
                axis_grp_dat$max),
         labels = rep("", 2))
    
    axis(side = 4,
         line = -53,
         tick = FALSE,
         cex.axis = 0.5,
         at = mean(c(axis_grp_dat$min,
                     axis_grp_dat$max)),
         labels = paste(strwrap(i, width = 20), 
                        collapse = "\n"))
  }
}

# legend("topright", legend = expression(bold("Set 1")),
#        bty = "n", cex = 1.1)

# fewer segregating sites than expected → candidate positive selection

# more segregating sites than expected → candidate balancing selection
# 
# arrows(x0 = 8, y0 = mx_y+8,      # start of arrow
#        x1 = 1, y1 = mx_y+8,    # end of arrow
#        code = 1,             # arrows at both ends
#        angle = 20,           # arrowhead angle
#        length = 0.1,         # arrowhead length
#        col = "grey30",
#        lwd = 2,
#        xpd = NA)             # allow drawing outside plot region
# 
# text(4,mx_y+9.5, cex = 0.95, col = "grey30",
#      "balancing selection", xpd = TRUE)
# 
# arrows(x0 = -3, y0 = mx_y+8,      # start of arrow
#        x1 = -1, y1 = mx_y+8,    # end of arrow
#        code = 1,             # arrows at both ends
#        angle = 20,           # arrowhead angle
#        length = 0.1,         # arrowhead length
#        col = "grey30",
#        lwd = 2,
#        xpd = NA)             # allow drawing outside plot region
# 
# text(-2,mx_y+9.5, cex = 0.95, col = "grey30",
#      "directional selection", xpd = TRUE)


extreme_labs <- cog_set[, .SD[which.max(signed_logp)], by = COG_letter][,.(gene_family, signed_logp, new_y)]

extreme_lo_labs <- cog_set[signed_logp< 0, .SD[which.min(signed_logp)], by = COG_letter][,.(gene_family, signed_logp, new_y)]

extreme_labs <- rbind(extreme_labs, extreme_lo_labs)

extreme_labs <- merge(extreme_labs, unique(pangraph_anno[,.(gene_family, consensus_gene_name, consensus_product)]),
                      all.x = TRUE, by = "gene_family")

extreme_labs[consensus_gene_name == "mdaB", consensus_product:= "Putative NADPH-quinone reductase mdaB"]
extreme_labs[consensus_gene_name == "wbdD", consensus_product:= "methyltransferase/kinase WbdD"]
extreme_labs[gene_family == "g003299", consensus_product:= "Transcriptional regulator AbiEi antitoxin"]
extreme_labs[consensus_gene_name == "wcaJ", consensus_product:= "Transferase wcaJ"]
# fwrite(extreme_labs, paste0(outdir_dat, "/extreme_labs.csv"))

with(extreme_labs,
     text(signed_logp, new_y,xpd = TRUE,
          labels = consensus_product, 
          cex = 0.45,
          pos = 3,
          offset = 0.25))


dev.off()



