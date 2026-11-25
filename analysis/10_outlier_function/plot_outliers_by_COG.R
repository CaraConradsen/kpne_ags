# don't double count genes!!

pangraph_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                       select = c("geno_id", "gene_family", 
                                  "fus_locus_tag","egng_cog", 
                                  "consensus_gene_name", "consensus_product",
                                  "tidy_name")))

# rename eggnog_cog to COG_funct_cat
setnames(pangraph_anno, old = "egng_cog", new = "COG_funct_cat")

cog_x <- "signed_logp"

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))

if (cog_x == "signed_logp") gene_pvals <- gene_pvals[p_bh <= 0.05]

gene_pvals[, s_per_site   := seg_sites / syn_sites]
gene_pvals[, log2_obs_exp := log2((seg_sites + 0.5) / (expected_S + 0.5))]

# Pangenome frequency ag_types ------------------------------------------------

cog_ag_types <- TRUE          # FALSE draws the pooled figure

n_genomes <- 260

# cloud <= 15%, shell > 15% and < 95%, soft-core >= 95%
gene_pvals[, ag_type := fcase(freq * 100 <= 15 * n_genomes, "Cloud (≤15%)",
                              freq * 100 <  95 * n_genomes, "Shell (>15% and <95%)",
                              default = "Soft-core (≥95%)")]

# keep as an ordered factor: ag_type_cols and b_off are built from levels()
gene_pvals[, ag_type := factor(ag_type, levels = c("Cloud (≤15%)",
                                                   "Shell (>15% and <95%)",
                                                   "Soft-core (≥95%)"))]

ag_type_cols <- setNames(c("#3FB0E8", "#B22222", "#00285B"),
                      levels(gene_pvals$ag_type))

print(gene_pvals[, .(n = .N, median_k = as.numeric(median(freq)),
                     pct_detectable = round(100 * mean(expected_S > 6))),
                 by = ag_type][order(ag_type)])


# Calculate signed log10(p) -----------------------------------------------

# signed magnitude: positive = excess, negative = deficit
gene_pvals[, direction := fifelse(seg_sites > expected_S, 1, -1)]
gene_pvals[, signed_logp := direction * -log10(p_two)]

gene_pvals[, cog_val := get(cog_x)]

cat(cog_x, ":", sprintf("%.2f to %.2f, median %.2f over %d genes",
                        min(gene_pvals$cog_val), max(gene_pvals$cog_val),
                        median(gene_pvals$cog_val), nrow(gene_pvals)), "\n")

# Plot by functional category -----------------------------------------
cog_set <- gene_pvals[,.(gene_family, signed_logp = cog_val, ag_type)]

# add cogs
cog_set[pangraph_anno[,.(gene_family, COG_funct_cat)], on = "gene_family",
        COG := COG_funct_cat
]

# How many dup COGs are there?
cat("Number of genes with dup COGS:", nrow(cog_set[nchar(COG) > 1]))
# Number of genes with dup COGS: 75 (new piS analysis)

cog_set <- rbind(cog_set[, .(COG_letter = unlist(strsplit(COG, ""))), by = c("gene_family","signed_logp","ag_type")],
                 cog_set[COG=="",.(gene_family, signed_logp, ag_type, COG_letter= COG)])

# add cog annotation from COG_dt
cog_set <- merge(cog_set, COG_dt,
                   all.x = TRUE, by = "COG_letter")

# cogs not represented
COG_dt[!broad_ord %chin% cog_set$broad_ord]

# tidy data
grp_ord  <- c("Information storage and processing", "Metabolism",
              "Cellular processes and signalling", "")
new_cols <- unique(cog_set[, .(broad_ord, desc)])[order(match(desc, grp_ord), broad_ord)][, .(broad_ord)]
new_cols[, new_y := .I]
mx_y = max(new_cols$new_y)
cog_set <- cog_set[new_cols, on = "broad_ord"]

# fwrite(cog_set, paste0(outdir_dat, "/bh_p_cog_set.csv"))

# x axis scale ------------------------------------------------------------
log_mod_x <- cog_x == "signed_logp"

lm_trans <- function(x) if (log_mod_x) sign(x) * log10(1 + abs(x)) else x

x_ticks <- if (log_mod_x) c(-5, -2, 0, 2, 5, 10, 25, 50, 100) else
  pretty(range(cog_set$signed_logp), 8)
x_ticks <- x_ticks[x_ticks >= min(cog_set$signed_logp) &
                     x_ticks <= max(cog_set$signed_logp)]

x_lab <- switch(cog_x,
                signed_logp  = expression("Signed" ~ -log[10] ~ italic(P)),
                s_per_site   = "Segregating sites per synonymous site",
                log2_obs_exp = expression(log[2]~"(observed / neutral expectation)"))

cog_set[, x_plot := lm_trans(signed_logp)]

# Sign of the departure
# violins and medians are drawn per side, only when side >4
cog_set[, side := fifelse(signed_logp < 0, -1L, 1L)]

# How many COGs categories are missing?
cat("Missing COGs:\n")
print(COG_dt[!COG_letter %chin% cog_set$COG_letter, .(COG_letter, COG_function, desc)])

# draw_cog_plot --------------------------------------------------------------

draw_cog_plot <- function(){

or_strip_w    <- 3.4    # width of one strip
or_strip_gap  <- 1    # gap between strips
or_n_strips   <- 4      # pooled, and one split by AG type
or_strip_in   <- 0.8    # gap between the innermost strip and the panel

# inner and outer edge of strip i, counting outward from the panel
or_strip_at <- function(i) {
  inn <- or_strip_in + (i - 1) * (or_strip_w + or_strip_gap)
  c(inn, inn + or_strip_w)
}

lab_line <- or_strip_at(or_n_strips)[2] + 0.4

par(mar = if (cog_ag_types) c(5, 19.1, 2, lab_line + 1.4) else c(4, 17, 1, 7))
plot(NULL,
     bty = "L", yaxt = "n", xaxt = "n",
     ylab = "",
     xlab = x_lab,
     ylim = c(mx_y + 1, 0),
     xlim = range(cog_set$x_plot) +
       c(-0.17, if (cog_ag_types) 0.17 else 0.05) * diff(range(cog_set$x_plot)))

# width of n margin lines, expressed in user x units. Defined here because
# the row shading needs it before the strips are drawn.
ln_x <- function(n) n * par("csi") * diff(par("usr")[1:2]) / par("pin")[1]

# Row shading --------------------------------------------------------------
# Alternating bands run from the panel straight through the odds-ratio
# strips, so a bar can be traced back to its category without counting rows
# across the gap. Drawn before anything else, so the violins, points and
# bars all sit on top of it; box() then restores the axis lines the shading
# would otherwise cover.
if (cog_ag_types) {
  shade_y <- seq(1, mx_y, by = 2)
  rect(par("usr")[2], shade_y - 0.5,
       par("usr")[2] + ln_x(or_strip_at(or_n_strips)[2] + 0.3), shade_y + 0.5,
       col = "grey95", border = NA, xpd = NA)
  box(bty = "L")
}

axis(side = 1, at = lm_trans(x_ticks), labels = x_ticks, col = "grey30")

# Category labels
y_lab_dt <- unique(cog_set[, .(new_y, COG_letter, COG_function)])[order(new_y)]
y_lab_dt[, lab := ifelse(COG_function == "Unassigned", COG_function,
                         paste0(COG_letter, ", ", COG_function))]

cog_abbr <- c(D = "Cell cycle con., cell div., chrom. part.",
              O = "Post-transl. mod., prot. turnover, chaperones",
              U = "Intra. traf., secret. and ves. trans.",
              Q = "Sec. metabol. biosyn., trans. and cata.")
y_lab_dt[COG_letter %chin% names(cog_abbr),
         lab := paste0(COG_letter, ", ", cog_abbr[COG_letter])]

axis(side = 2, at = y_lab_dt$new_y, labels = FALSE, col = "grey30")
axis(side = 2, at = y_lab_dt$new_y, labels = y_lab_dt$lab, tick = FALSE,
     line = 0.3, cex.axis = 0.72, col.axis = "black", las = 2)

# second y-axis
axis(side = 4, at = y_lab_dt$new_y, labels = y_lab_dt$COG_letter, tick = FALSE,
     line = lab_line - 0.4, cex.axis = 1, col.axis = "grey30", las = 2)

abline(v = 0, col = "grey30", lwd = 1, lty = 2)

# vertical offsets of the three ag_type violins within a COG row
b_off <- setNames(c(-0.24, 0, 0.24), levels(cog_set$ag_type))

# top of a row's drawn content, where the gene labels are anchored
row_top <- 0.24 + 0.10

# add text label
mtext("a", side = 3, adj = -0.05, font = 2, line = -1, cex = 8/7)

# Gene label text ----------------------------------------------------------
# Labels come from tidy_name
max_lab <- 24  

lab_short <- c(
  "Transcriptional regulator AbiEi antitoxin N-terminal domain" =
    "AbiEi antitoxin N-term domain")


# Smallest cell drawn as a violin
min_vio <- 3

if (cog_ag_types) {

  for (i in 1:mx_y) {
    for (b in levels(cog_set$ag_type)) {
      for (sd in c(-1L, 1L)) {
        sub <- cog_set[new_y == i & ag_type == b & side == sd]
        if (nrow(sub) >= min_vio) {
          vioplot(sub$x_plot, at = i + b_off[[b]],
                  col = adjustcolor(ag_type_cols[[b]], 0.55),
                  border = NA, horizontal = TRUE, drawRect = FALSE,
                  wex = 0.20, add = TRUE)
        }
      }
    }
  }

} else {

for (i in 1:mx_y) {
  if(nrow(cog_set[new_y == i]) >= 5){
    with(cog_set[new_y == i],
         vioplot(x_plot, at = i, 
                 col = "grey80",
                 border = NA,
                 horizontal = TRUE,
                 drawRect = FALSE,
                 wex = 0.75,
                 add = TRUE))
  }
}

}

for (i in 1:mx_y) {
  if (cog_ag_types) {
    for (b in levels(cog_set$ag_type)) {
      sub <- cog_set[new_y == i & ag_type == b]
      if (!nrow(sub)) next
      points(sub$x_plot, sub$new_y + b_off[[b]], pch = 16, cex = 0.5,
             col = adjustcolor(ag_type_cols[[b]], 0.85))
      # one median per side, matching the violins
      for (sd in c(-1L, 1L)) {
        ss <- sub[side == sd]
        if (nrow(ss) >= min_vio)
          points(median(ss$x_plot), i + b_off[[b]], pch = 23,
                 bg = "white", col = "grey20", cex = 0.75, lwd = 0.6)
      }
    }
    next
  }
  with(cog_set[new_y == i],
       points(x_plot, new_y, pch = 16, cex = 0.55,
              col = rgb(0.1,0.1,0.1, alpha = 0.75)))
  
  if(nrow(cog_set[new_y == i]) >= 10){
    with(cog_set[new_y == i],
         points(median(x_plot), i, pch = 23,
                bg = "white", cex = 1.2)
    )
  }
}


# odds-ratio strips ---------------------------------------------------------

or_file      <- paste0(outdir_dat, "/cog_fisher_outliers.csv")
or_clip      <- 2                       # log2 OR clipped to +/- this
or_ticks     <- c(-1, 0, 1)             # axis ticks, in log2 OR
or_tick_labs <- c("0.5", "1", "2")

# suffixes cog_enrichment_fisher.R appends when it stratifies
ag_type_key <- setNames(c("cloud", "shell", "soft"), levels(cog_set$ag_type))

# the four strips, in the order they read left to right
or_spec <- list(
  list(cmp = "low vs non-outlier",  by_type = FALSE, title = "all AGs"),
  list(cmp = "low vs non-outlier",  by_type = TRUE,  title = "by AG type"),
  list(cmp = "high vs non-outlier", by_type = FALSE, title = "all AGs"),
  list(cmp = "high vs non-outlier", by_type = TRUE,  title = "by AG type")
)

# the two strips of each direction are bracketed under one label
or_groups <- list(list(a = 1L, b = 2L, title = "Low diversity"),
                  list(a = 3L, b = 4L, title = "High diversity"))

if (cog_ag_types && file.exists(or_file)) {

  or_all <- fread(or_file)

  or_get <- function(cmp) {
    d <- or_all[comparison == cmp, .(COG_letter = cog, OR, p_bh)]
    d <- d[y_lab_dt, on = "COG_letter"][!is.na(OR) & OR > 0]
    d[, l2   := pmax(-or_clip, pmin(or_clip, log2(OR)))]
    d[, star := fcase(p_bh <= 0.001, "***",
                      p_bh <= 0.01,  "**",
                      p_bh <= 0.05,  "*",
                      default = "")]
    d[]
  }

  u <- par("usr")

  # strip i, counting outward from the right edge of the panel
  strip_geom <- function(i) {
    list(mid = u[2] + ln_x(mean(or_strip_at(i))),
         hw  = ln_x(or_strip_w) / 2)
  }

  or_xat <- function(g, l2) g$mid + g$hw * l2 / or_clip

  or_bars <- function(g, d, bar_col, yoff, half_h) {
    if (!nrow(d)) return(0L)
    for (sig in c(FALSE, TRUE)) {
      dd <- d[(star != "") == sig]
      if (!nrow(dd)) next
      rect(g$mid, dd$new_y + yoff - half_h,
           or_xat(g, dd$l2), dd$new_y + yoff + half_h,
           col = adjustcolor(bar_col, if (sig) 0.95 else 0.30),
           border = NA, xpd = NA)
    }
    d[star != "", .N]
  }

  # BH-significant categories are starred as well as filled, so the by-type
  # strips read without having to judge a thin bar by its fill alone.
  or_stars <- function(g, d, yoff, cx) {
    dd <- d[star != ""]
    if (!nrow(dd)) return(invisible(NULL))
    text(or_xat(g, dd$l2), dd$new_y + yoff, dd$star,
         pos = fifelse(dd$l2 >= 0, 4, 2), offset = 0.08,
         cex = cx, col = "grey20", xpd = NA)
  }

  or_axis <- function(g, title) {
    y0 <- mx_y + 1          # ylim is reversed, so larger y is further down
    segments(or_xat(g, -or_clip), y0, or_xat(g, or_clip), y0,
             col = "grey30", lwd = 0.8, xpd = NA)
    segments(or_xat(g, or_ticks), y0, or_xat(g, or_ticks), y0 + 0.18,
             col = "grey30", lwd = 0.8, xpd = NA)
    text(or_xat(g, or_ticks), y0 + 0.70, or_tick_labs,
         cex = 0.72, col = "grey20", xpd = NA)
    text(g$mid, y0 + 1.55, title, cex = 0.72, col = "grey20", xpd = NA)
    segments(g$mid, 0.4, g$mid, mx_y + 0.6, col = "grey70", lwd = 0.6,
             xpd = NA)
  }

  draw_or_strip <- function(i, sp) {

    g <- strip_geom(i)
    or_axis(g, sp$title)

    # pooled: one bar per category
    if (!sp$by_type) {
      d <- or_get(sp$cmp)
      n <- or_bars(g, d, "grey35", 0, 0.30)
      or_stars(g, d, 0, 1)
      return(n)
    }

    # by type: one bar per AG type, at its violin's offset and colour
    sum(vapply(levels(cog_set$ag_type), function(b) {
      d <- or_get(paste0(sp$cmp, ", ", ag_type_key[[b]]))
      n <- or_bars(g, d, ag_type_cols[[b]], b_off[[b]], 0.15)
      or_stars(g, d, b_off[[b]], 1)
      n
    }, numeric(1)))
  }

  # bracket spanning the pair of strips that share a comparison
  or_group <- function(gr) {
    ga <- strip_geom(gr$a)
    gb <- strip_geom(gr$b)
    x0 <- min(ga$mid - ga$hw, gb$mid - gb$hw)
    x1 <- max(ga$mid + ga$hw, gb$mid + gb$hw)
    y0 <- mx_y + 0.75
    segments(x0, y0 + 2.30, x1, y0 + 2.30, col = "grey30", lwd = 0.8,
             xpd = NA)
    segments(c(x0, x1), y0 + 2.30, c(x0, x1), y0 + 2.12, col = "grey30",
             lwd = 0.8, xpd = NA)
    text((x0 + x1) / 2, y0 + 2.95, gr$title, cex = 0.72, col = "grey20",
         xpd = NA)
  }

  # strip 1 sits against the panel's right edge, so the spec order is the
  # order the strips appear
  for (i in seq_along(or_spec)) {
    sp <- or_spec[[i]]
    n  <- draw_or_strip(i, sp)
    cat(sprintf("OR strip %-14s %-22s %d at BH 0.05\n",
                sp$title, sp$cmp, n))
  }

  invisible(lapply(or_groups, or_group))
  
}

mtext("Odds ratio", side = 1, line = 3.5,
      at = par("usr")[2] + ln_x(mean(c(or_strip_at(1)[1], or_strip_at(or_n_strips)[2]))))
mtext("b", side = 3, adj =1.05, font = 2, line = -1)

# add legend
pt_lab <- if (cog_x == "signed_logp") "AG outliers" else "Accessory genes"

if (cog_ag_types) {

  # above the plot box, outside it. ylim is reversed, so the top of the
  # panel is usr[4] and moving further that way needs a negative step.
  u <- par("usr")

  l_txt <- c(names(ag_type_cols), "Median")
  l_y   <- u[4] + (u[4] - u[3]) * 0.05

  # shrink until the row fits the device, then left-align it on the box
  l_x   <- mean(u[1:2])
  l_cex <- 1
  repeat {
    lw <- legend(x = l_x, y = l_y, plot = FALSE, legend = l_txt, xjust = 0.5,
                 pch = c(22, 22, 22, 23), horiz = TRUE, bty = "n",
                 cex = l_cex, x.intersp = 0.7,
                 text.width = strwidth(l_txt, cex = l_cex))$rect$w
    if (lw <= diff(grconvertX(c(0.02, 0.98), "ndc", "user")) || l_cex < 0.72) break
    l_cex <- l_cex - 0.02
  }

  legend(x = l_x, y = l_y, xpd = NA, xjust = 0.5,
         legend = l_txt,
         pch = c(22, 22, 22, 23),
         pt.bg = c(adjustcolor(ag_type_cols, 0.55), "white"),
         col = c(rep(NA, 3), "grey20"),
         pt.cex = c(1.1, 1.1, 1.1, 0.75),
         horiz = TRUE, bty = "n", cex = l_cex, x.intersp = 0.7,
         text.width = strwidth(l_txt, cex = l_cex))

} else {

legend("top",c(pt_lab, "Median"),#c("Set 1", "Set 2"),
       pch = c(16,23),#,21), 
       col = c(rgb(0.1,0.1,0.1, alpha = 0.75), "black"),
               # "darkgoldenrod3"),
       pt.bg = c(NULL,"white"), 
       # pt.bg = c(NULL,"darkgoldenrod1"),
       pt.cex = c(0.72,1.2),
       horiz = TRUE,
       bty = "n", cex = 0.95)

}



# add broad categories
# anchored off the current margins, so that they track the device width
# rather than the hard-coded offsets of the narrower figure
# side 4 lines run outward from the right edge of the plot region, so the
# offset is the plot width in lines plus the depth into the left margin
brd_tick <- -(par("pin")[1] / par("csi") + 15.9)
brd_lab  <- brd_tick - 2.5

for (i in unique(cog_set$desc)) {
  if(i != ""){
    axis_grp_dat <- cog_set[desc==i]
    
    n_cogs = length(unique(axis_grp_dat$COG_letter))
    
    axis_grp_dat <- axis_grp_dat[, .(min = min(new_y),
                                     max = max(new_y))]
    axis(side = 4,
         line = brd_tick,
         at = c(axis_grp_dat$min,
                axis_grp_dat$max),
         labels = rep("", 2))
    
    axis(side = 4,
         line = brd_lab,
         tick = FALSE,
         cex.axis = 1,
         at = mean(c(axis_grp_dat$min,
                     axis_grp_dat$max)),
         labels = paste(strwrap(i, width = 20), 
                        collapse = "\n"))
  }
}


if (cog_ag_types) {

  # the most extreme gene each way per category. Labels are anchored above
  # the whole row rather than against their own ag_type, so that they sit in
  # the gap between rows instead of over a neighbouring violin; the x
  # position still identifies which gene each one belongs to.
  extreme_labs <- rbind(
    cog_set[, .SD[which.max(signed_logp)], by = COG_letter],
    cog_set[signed_logp < 0, .SD[which.min(signed_logp)], by = COG_letter])

  extreme_labs <- extreme_labs[, .(gene_family, COG_letter, signed_logp,
                                   x_plot, new_y = new_y - row_top)]

} else {

extreme_labs <- cog_set[, .SD[which.max(signed_logp)], by = COG_letter][,.(gene_family, signed_logp, x_plot, new_y)]

extreme_lo_labs <- cog_set[signed_logp< 0, .SD[which.min(signed_logp)], by = COG_letter][,.(gene_family, signed_logp, x_plot, new_y)]

extreme_labs <- rbind(extreme_labs, extreme_lo_labs)

}

extreme_labs <- merge(extreme_labs,
                      unique(pangraph_anno[,.(gene_family, tidy_name,
                                              consensus_gene_name, consensus_product)]),
                      all.x = TRUE, by = "gene_family")

# mdaB, wbdD and wcaJ no longer need a hand-written product here: wbdD is
# now named from Kaptive and the other two are corroborated by their own
# product, so tidy_name resolves all three.
# fwrite(extreme_labs, paste0(outdir_dat, "/extreme_labs.csv"))

if (cog_ag_types) {

  extreme_labs[, lab := tidy_name]
  
  # gene symbols e.g. mobH to italic; set before lab_short and truncation change lab
  extreme_labs[, is_gene_symbol := grepl("^[a-z]{3}[A-Z]?[0-9]*$", lab)]
  

  hit <- match(extreme_labs$lab, names(lab_short))
  extreme_labs[!is.na(hit), lab := unname(lab_short[hit[!is.na(hit)]])]

  extreme_labs[, hw := strwidth(lab, cex = 0.72) / 2]
  extreme_labs[, x_lab := fifelse(signed_logp > 0,
                                  pmin(pmax(x_plot, lm_trans(0) + ln_x(0.3) + hw),
                                       par("usr")[2] - ln_x(0.3) - hw),
                                  pmin(pmax(x_plot, par("usr")[1] + ln_x(0.3) + hw),
                                       min(cog_set[signed_logp > 0, x_plot]) - ln_x(0.3) - hw))]

  with(extreme_labs,
       text(x_lab, new_y, xpd = NA,
            labels = lab,
            font = fifelse(is_gene_symbol, 3, 1),
            cex = 0.72,
            pos = 3,
            offset = 0.1))

} else {
  with(extreme_labs,
       text(x_plot, new_y,xpd = TRUE,
            labels = consensus_product, 
            cex = 0.72,
            pos = 3,
            offset = 0.25))
}
}



# plot outputs ------------------------------------------------------------

png(paste0(outdir_fig, "/piS_theta_cog_by_ag_type_set_1.png"),
    width = 180, height = 125,
    units = "mm", res = 300,
    pointsize = 7, type = "cairo")
draw_cog_plot(); dev.off()


cairo_ps(paste0(outdir_fig, "/ms_eps/piS_theta_cog_by_ag_type_set_1.eps"),
         width = 180 / 25.4, height = 125 / 25.4,
         pointsize = 7, family = "Arial")

draw_cog_plot(); dev.off()



