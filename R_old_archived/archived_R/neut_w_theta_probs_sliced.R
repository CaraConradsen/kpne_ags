
prob_theta <- fread(paste0(outdir_dat, "/prob_theta_norm.csv"))

# slice_freq <- round(quantile(prob_theta$freq, probs = c(0,0.25,0.5,0.75,1)), digits = 0)

slice_freq =c(6, 32, 119, 199, 249)

s_range = range(as.integer(gsub("s_","", names(prob_theta[,-1]))))

# slice the probability of the outlier

# ---- build the surface from your neutral expectation table ----
z <- as.matrix(prob_theta[, !"freq"])                        # rows = freq, cols = s_0..s_n
x <- prob_theta$freq                                          # pangenome frequency
y <- as.numeric(sub("s_", "", colnames(z)))           # segregating sites (numeric)

obs <- fread(paste0(outdir_dat, "/set_outliers.csv"))

# look at waterson's theta only, set 1
obs <- obs[set==1 & !is.na(w_theta_sel)]

get_p   <- function(f, s) z[match(f, x), match(s, y)]
obs$p   <- mapply(get_p, obs$freq, round(obs$seg_s_sites))

# plot

for (i in 1:5) {
  fdat = t(prob_theta[freq==slice_freq[i]][,-1])
  
  fdat = cbind(fdat, s_range[1]:s_range[2])
  
  fdat = as.data.frame(fdat)
  
  colnames(fdat) = c("P", "S")
  
  if(i==1){
    par(mar=c(4,4.5,3,0.5))
  }else{
    par(mar=c(4,0.5,3,0.5))
  }
  
  # plot
  plot(NULL, pch = 16,
       xlim = c(s_range[1], 20),  
       ylim = 0:1,
       yaxt = "n",
       ylab = "",
       xlab = "")

  polygon(
    c(fdat$S, rev(fdat$S)),
    c(fdat$P, rep(0,length(fdat$P))),
    col = "grey40",
    border = "grey40"
  )
  
  if(i==1){
    axis(side = 2, at = seq(0,1, by=0.2),
         labels = sprintf("%1.1f", seq(0,1, by=0.2)), las = 2)
    axis(side = 2, at =0.5, tick = FALSE,line = 2,
         labels = "Neutral probability")
  }
  
  if(i==3){
    axis(side = 1, at = 10, tick = FALSE,line = 1.5,
         labels = "Number of synonymous segregating sites", xpd = TRUE)
  }
  
  legend("top", legend = bquote(italic(k)~"="~.(slice_freq[i])),
         bty = "n", cex = 1.1)
  
  #add outliers
  with(obs[freq==slice_freq[i]],
       points(seg_s_sites, p, pch = 16, cex = 1, 
              col = ifelse(w_theta_sel == "balancing" &  p_bh <= 0.05, "#3FA9F5", 
                           ifelse(w_theta_sel == "directional" &  p_bh <= 0.05, "#F26C2A", 
                                  rgb(0.1,0.1,0.1, alpha = 0.5))))
       )
}




