# plot the pangenome 

pan_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"))


# Assign Core vs accessory ------------------------------------------------

n_genomes = length(unique(pan_anno[, geno_id]))


pan_anno[, ag_type := fcase(number_genomes < 0.99 * n_genomes & number_genomes >=  0.95 * n_genomes, "soft",
                            number_genomes < 0.95 * n_genomes & number_genomes >=  0.15 * n_genomes, "shell",
                            number_genomes < 0.15 * n_genomes, "cloud",
                            default = "core")]


# plot distributions ------------------------------------------------------

# pie dat
pie_dat = unique(pan_anno[,.(gene_family, ag_type)])


pie_dat = pie_dat[,.(n = .N), by=ag_type]

pie_dat[, lab := paste0(ag_type, " (", 
                        formatC(n, format = "d", big.mark = ","), 
                        ")")]

# frequncy data
freq_dat = unique(pan_anno[,.(gene_family, number_genomes)])

freq_dat = freq_dat[,.(n = .N), number_genomes]

freq_dat[, log10_n := log10(n)]

freq_dat = merge(as.data.table(
  data.frame(number_genomes = 1:260)
), freq_dat,
all.x = TRUE, by = "number_genomes")

freq_dat[is.na(log10_n), `:=`(n = 0, log10_n = 0)]


png(paste0(outdir_fig,"/core_pangenome_stats.png"),
    width = 17, height = 21.7, units = "cm", res = 300,
    pointsize = 12, type = "cairo")


# plots
mat <- matrix(c(1,1,2,rep(3,3)), byrow =TRUE, ncol = 3)


layout(mat, widths = c(0.75,0.75,1),
       heights = c(1,2))

# gene distribution
par(mar = c(6,6,1,0.5))
with(freq_dat,
     barplot(log10_n~number_genomes,
             border = NA,
             yaxt = "n",
             col = "grey40",
             ylab = expression("log"[10]~"Number of genes"),
             xlab = "Number of genomes"))

abline(h=0)
axis(side=2, at = seq(0,3.5,0.5),
     labels = sprintf("%.1f", seq(0,3.5,0.5)),
     las = 2)

usr <- par("usr")

text(x = usr[1]- 55, y = usr[4], labels = "a", 
     xpd = TRUE, font = 2,cex = 1.5)

# Green-themed colours
blue_palette <- c("grey80", "lightskyblue3",  "steelblue", "dodgerblue")  # dark to light greens

par(mar = c(0.1,0.1,1,6))
order_idx <- c(1,2,3,4)
# Create pie chart
with(pie_dat,
     pie(n[order_idx],
         labels = c(lab[order_idx][1:3], ""),
         col = blue_palette,
         main = "",
         border = "#fafaf8",
         clockwise = TRUE,
         init.angle = 90
     )
)


# adjust label
angle_radians <- 2 * pi * 14258/2 / sum(pie_dat$n)


# Radius at which to place labels (1 is on the circumference)
r <- -0.9
x_pos <- r * cos(angle_radians)
y_pos <- r * sin(angle_radians)


# Add text at custom positions
text(x_pos-0.05, y_pos, xpd=TRUE, pos =4,
     labels = pie_dat[ag_type=="cloud", lab])

r <- -0.8
x_pos1 <- r * cos(angle_radians)
y_pos1 <- r * sin(angle_radians)

r <- -0.85
x_pos2 <- r * cos(angle_radians)
y_pos2 <- r * sin(angle_radians)

segments(x_pos1,y_pos1,x_pos2,y_pos2)

usr <- par("usr")

text(x = usr[1]+.1, y = usr[4], labels = "b", 
     xpd = TRUE, font = 2,cex = 1.5)

# Add phylogenetic tree
par(mar = c(2,1,3,1))
tree <- read.tree("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/gub_graph.node_labelled.final_tree.tre") 

ST_groups <- unique(pan_anno[,.(geno_id, ST)])

# get colours
ST_cols = ST_groups[,.(n=.N), ST][n>1, .(ST)]
setorderv(ST_cols, "ST",1)

ST_cols[, cols := colorRampPalette(c("#420D55","blueviolet",
                          "darkblue","#008080",
                          "forestgreen","#799000",
                          "yellow","#FF8000",
                          "brown"))(nrow(ST_cols))]

ST_groups <- merge(ST_groups, ST_cols,
                   all.x = TRUE, by ="ST")
ST_groups[is.na(cols), cols := "black" ]

# re-label with ST
ST_groups <- ST_groups[match(tree$tip.label, geno_id)]  # reorder

tree$tip.label <- paste0("  ", tree$tip.label, "  ")

ST_seg = ST_groups$ST
names(ST_seg) = tree$tip.label

ST_cols = ST_groups$cols
names(ST_cols) = tree$tip.label

# Plot circular tree
lim <- c(-1.2, 1.2)
plot(tree,
     edge.width = 0.8,
     type = "fan",
     cex = 0.4,
     tip.color = "grey20",
     use.edge.length = FALSE,
     show.tip.label = TRUE,
     no.margin = TRUE, x.lim = lim, y.lim = lim)


lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
ntips <- Ntip(tree)

# Blue filled tip points
points(lp$xx[1:ntips],
       lp$yy[1:ntips],
       pch = 16,
       cex = 0.75,
       col = ST_cols)


# Calculate angles of tips in radians
tip_angles <- atan2(lp$yy[1:ntips], lp$xx[1:ntips])


# Find consecutive runs
split_consecutive_angle <- function(indices, angle_order){
  pos <- match(indices, angle_order)
  pos <- sort(pos)
  breaks <- c(TRUE, diff(pos) != 1)
  split(indices[order(pos)], cumsum(breaks))
}

draw_clade_arc <- function(tip_indices, col="grey30", lwd=2, radius=1.17, label_radius=1.18,
                           label_text=NULL){
  
  if(length(tip_indices) == 0) return()  # skip empty
  
  # angles for these tips
  angles <- tip_angles[tip_indices]
  
  # unwrap if boundary crossing
  if(diff(range(angles)) > pi){
    angles[angles < mean(angles)] <- angles[angles < mean(angles)] + 2*pi
  }
  
  # arc coordinates
  theta <- seq(min(angles), max(angles), length.out=100)
  theta <- ifelse(theta > pi, theta - 2*pi, theta)
  
  x <- radius * cos(theta)
  y <- radius * sin(theta)
  lines(x, y, col=col, lwd=lwd, lend=2)
  
  # label
  if(!is.null(label_text)){
    label_angle <- mean(angles)
    if(label_angle > pi) label_angle <- label_angle - 2*pi
    
    xpos <- label_radius * cos(label_angle)
    ypos <- label_radius * sin(label_angle)
    
    angle_deg <- label_angle * 180/pi
    on_left <- xpos < 0
    if(on_left){
      angle_deg <- angle_deg + 180
      adj_val <- c(1,0.5)
    } else {
      adj_val <- c(0,0.5)
    }
    
    text(xpos, ypos,
         labels = label_text,
         srt = angle_deg,
         adj = adj_val,
         cex = 0.5,
         font = 2)
  }
}

# loop through and ad arcs
for(i in unique(ST_groups$ST)){
  idx <- which(ST_seg == i)
  
  # optional: skip empty
  if(length(idx)==0) next
  
  # determine blocks in angular order
  blocks <- split_consecutive_angle(idx, order(tip_angles))
  
  for(block in blocks){
    draw_clade_arc(block, col="grey30", lwd=2, 
                   label_text=ifelse(i!="ST512", i, ""))
  }
  
  if(i=="ST512"){
    split_idx <- idx[ceiling(length(idx)/2)]
    xpos <- lp$xx[split_idx] * 1.18
    ypos <- lp$yy[split_idx]
    
    label_angle <- atan2(ypos, xpos)
    angle_deg <- label_angle * 180 / pi
    
    # Flip text if on left side
    on_left <- xpos < 0
    if(on_left){
      angle_deg <- angle_deg + 180
      adj_val <- c(1, 0.5)  # right-justified
    } else {
      adj_val <- c(0, 0.5)  # left-justified
    }
    
    text(x = xpos, y = ypos,
         labels = "ST512",
         srt = angle_deg,
         adj = adj_val,
         cex = 0.5,
         font = 2)
    
    segments(xpos+0.01,ypos-0.02,xpos+0.01,ypos+0.03,
             col = "grey30",
             lwd=2, lend=2)
    
  }
  
}

# ---- Scale bar (bottom-right) ----
usr <- par("usr")

add.scale.bar()

text(x = usr[1]+0.2, y = usr[4]-0.3, labels = "c", 
     xpd = TRUE, font = 2,cex = 1.5)

dev.off()


