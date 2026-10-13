png(paste0(outdir_fig,"/3d_outliers_seg_sites_pi.png"),
    width = 32, height = 24, units = "cm", res = 300,
    pointsize = 14, type = "cairo")

prob_pi <- fread(paste0(outdir_dat, "/prob_pi_norm.csv"))

# ---- build the surface from your neutral expectation table ----
z <- as.matrix(prob_pi[, !"freq"])                        # rows = freq, cols = s_0..s_n
x <- prob_pi$freq                                          # pangenome frequency
y <- as.numeric(sub("s_", "", colnames(z)))           # segregating sites (numeric)

# ---- draw the surface, capture the 3D->2D projection matrix ----
pmat <- persp(x, y, z,
              theta = 60, phi = 30, expand = 0.75,
              xlab = "\nPangenome frequency (k)",
              ylab = "\nSegregating synonymous sites",
              zlab = "\nProbability",
              ticktype = "detailed",
              shade = 0.05, 
              col = "white", border = "grey15", lwd = 0.4)

# ---- overlay observed AGs ----
# obs is a data.frame with columns: freq, seg_s_sites, sel (sel = "balancing"/"directional")
# decide what z should be for each point — here, the neutral probability at (freq, s):

obs <- fread(paste0(outdir_dat, "/set_outliers.csv"))
colnames(obs)[8] = "set"

# look at pi only, set 1
obs <- obs[set=="set_1" & !is.na(pi_sel)]

get_p   <- function(f, s) z[match(f, x), match(s, y)]
obs$p   <- mapply(get_p, obs$freq, round(obs$seg_s_sites))

pts  <- trans3d(obs$freq, obs$seg_s_sites, obs$p, pmat)
base <- trans3d(obs$freq, obs$seg_s_sites, 0,     pmat)

segments(base$x, base$y, pts$x, pts$y, col = "grey30", lwd = 0.6)
points(pts, pch = 16, cex = 1, 
       col = ifelse(obs$pi_sel == "balancing", "#3FA9F5", "#F26C2A"))

dev.off()