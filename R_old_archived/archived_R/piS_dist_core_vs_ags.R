
# Import pis values -------------------------------------------------------

ag_piS_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

core_piS_dt <- fread(paste0(outdir_dat, "/core_theta_dt.csv"))

# sum per gene

ag_piS_dt <- ag_piS_dt[,.(
  pi_s = mean(pi_s, na.rm = TRUE),
  freq = mean(number_genomes)
), gene_family]

core_piS_dt <- core_piS_dt[,.(
  pi_s = mean(pi_s, na.rm = TRUE)
), gene_family]



# plot --------------------------------------------------------------------

par(mfrow=c(2,1), mar = c(3,3,0.5,1),
    oma = c(1.5,1.5,0,0))

plot(NULL,
     bty = "n",
     xaxt = "n", 
     yaxt = "n",
     ylab = "", 
     xlab = "",
     ylim = c(0, 3),
     xlim = c(0, 0.18))

axis(side = 2, at = seq(0,3, 0.5),
     labels = sprintf("%1.1f",seq(0,3, 0.5)), las = 2)

axis(side = 1, at = seq(0,0.18, 0.02),
     labels = sprintf("%1.2f",seq(0,0.18, 0.02)))

core_hist <- with(core_piS_dt,
                  hist(pi_s, breaks = 80, plot = FALSE))

rect(xleft = core_hist$breaks[-length(core_hist$breaks)],
     xright = core_hist$breaks[-1],
     ybottom = 0,
     ytop = log10(core_hist$counts+1),
     col = "tomato", border = F)

mtext("Core genes", side = 3, adj = 0.9, line = -2, col = "tomato")


plot(NULL,
     bty = "n",
     xaxt = "n", 
     yaxt = "n",
     ylab = "", 
     xlab = "",
     ylim = c(0, 3.5),
     xlim = c(0, 0.18))

axis(side = 2, at = seq(0,3.5, 0.5),
     labels = sprintf("%1.1f",seq(0,3.5, 0.5)), las = 2)

axis(side = 1, at = seq(0,0.18, 0.02),
     labels = sprintf("%1.2f",seq(0,0.18, 0.02)))


ag_hist <- with(ag_piS_dt,
                hist(pi_s, breaks = core_hist$breaks, plot = FALSE))

rect(xleft = ag_hist$breaks[-length(ag_hist$breaks)],
     xright = ag_hist$breaks[-1],
     ybottom = 0,
     ytop = log10(ag_hist$counts+1),
     col = "steelblue", border = F)

mtext("Accessory genes", side = 3, adj = 0.9, line = -2, col = "steelblue")



mtext(expression("Log"[10]~"(number of genes + 1)"),
      side  = 2 , outer = TRUE)

mtext(expression("Synonymous diversity ("*pi[S]*")"),
      side  = 1 , outer = TRUE)




# Non-logged --------------------------------------------------------------

layout(matrix(c(1,2), nrow = 2), 
       heights = c(0.33, 0.66))

par(mar = c(3,5,0.5,1),
    oma = c(2.5,0,0,0))

plot(NULL,
     bty = "n",
     xaxt = "n", 
     yaxt = "n",
     xlab = "", 
     ylab = "Number of genes",
     ylim = c(0, 500),
     xlim = c(0, 0.18))

axis(side = 2, at = seq(0,500, 100),
     labels = sprintf("%3.0f",seq(0,500, 100)), las = 2)

axis(side = 1, at = seq(0,0.18, 0.02),
     labels = sprintf("%1.2f",seq(0,0.18, 0.02)))

core_hist <- with(core_piS_dt,
                  hist(pi_s, breaks = 80, plot = FALSE))

rect(xleft = core_hist$breaks[-length(core_hist$breaks)],
     xright = core_hist$breaks[-1],
     ybottom = 0,
     ytop = core_hist$counts,
     col = "tomato", border = "tomato4")

mtext("Core genes", side = 3, adj = 0.9, line = -2, col = "tomato4")


plot(NULL,
     bty = "n",
     xaxt = "n", 
     yaxt = "n",
     ylab = "Number of genes", 
     xlab = "",
     ylim = c(0, 1500),
     xlim = c(0, 0.18))

axis(side = 2, at = seq(0,1500, 500),
     labels = sprintf("%3.0f",seq(0,1500, 500)), las = 2)

axis(side = 1, at = seq(0,0.18, 0.02),
     labels = sprintf("%1.2f",seq(0,0.18, 0.02)))


ag_hist <- with(ag_piS_dt,
                hist(pi_s, breaks = core_hist$breaks, plot = FALSE))

rect(xleft = ag_hist$breaks[-length(ag_hist$breaks)],
     xright = ag_hist$breaks[-1],
     ybottom = 0,
     ytop = ag_hist$counts,
     col = "steelblue", border = "steelblue4")

mtext("Accessory genes", side = 3, adj = 0.9, line = -2, col = "steelblue4")


mtext(expression("Synonymous diversity ("*pi[S]*")"),
      side  = 1 , line =1, outer = TRUE)



dev.off()



# by freq
with(ag_piS_dt,
     boxplot(pi_s~freq, col = "dodgerblue",
             pch = 16, cex = 0.65,
             xlab = "Pangenome frequency",
             ylab = expression("Synonymous diversity ("*pi[S]*")")))

abline(h = mean(core_piS_dt$pi_s),
       lty = 2, col = "red")


