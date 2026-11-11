# ST composition by isolation source ----------------------------------------
# Stacked bar of the 260 sequenced genomes, STs on the x-axis and SpARK
# isolation source stacked within each bar. Source categories and colours
# follow Fig. 1 of Thorpe et al. 2022 (Nat Microbiol 7:2054-2067).
#
# NOTE: the hex values below were read off the published legend, not supplied
# by the authors, so they are approximations. They are kept in one named
# vector so they can be corrected in a single place.

spark_meta  <- "../SpARK data/spark_metadata.csv"   # relative to the project root
min_n       <- 2                                    # STs below this are pooled as "Other"
collapse_lv <- TRUE                                 # treat ST307-1LV as ST307


# Source levels and colours, in the order Thorpe et al. Fig. 1 lists them -----

src_cols <- c(
  "human - community carriage" = "#A6CEE3",   # pale blue
  "human - hospital carriage"  = "#3B8EA5",   # mid teal-blue
  "human - hospital disease"   = "#1F3B73",   # dark navy
  "dog"                        = "#FBB4B9",   # pale pink
  "cat"                        = "#F4A582",   # light salmon
  "cow"                        = "#E31A1C",   # red
  "pig"                        = "#FF7F00",   # orange
  "sheep"                      = "#FFFFB3",   # cream
  "other mammals"              = "#CAB2D6",   # lavender
  "poultry"                    = "#6A3D9A",   # dark purple
  "other birds"                = "#9ECAE1",   # periwinkle
  "turtle"                     = "#33A02C",   # bright green
  "invertebrates"              = "#FFFF33",   # yellow
  "water - environment"        = "#1F78B4",   # bright blue
  "water - river"              = "#8B1A1A",   # dark red
  "water - farm"               = "#08306B",   # deep blue
  "surface - environment"      = "#9E9E6B",   # olive
  "surface - hospital"         = "#CCCCCC",   # light grey
  "surface - farm"             = "#666666",   # dark grey
  "soil"                       = "#8B4513",   # brown
  "epiphytes"                  = "#B2DF8A",   # light green
  "endophytes"                 = "#66BD63",   # mid green
  "rhizosphere"                = "#00CED1",   # cyan
  "rhizoplane"                 = "#000000"    # black
)


# Load -----------------------------------------------------------------------

# the 260 sequenced genomes
genomes <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                        select = "geno_id")$geno_id)

meta <- fread(spark_meta, select = c("id", "ST", "group_summary"))
setnames(meta, "id", "geno_id")

dt <- meta[geno_id %chin% genomes]

if (nrow(dt) != length(genomes))
  warning(length(genomes) - nrow(dt), " genomes not found in the SpARK metadata")


# Tidy ST and source ----------------------------------------------------------

if (collapse_lv) dt[, ST := sub("-[0-9]+LV$", "", ST)]

dt[, group_summary := trimws(tolower(group_summary))]

unknown <- setdiff(unique(dt$group_summary), names(src_cols))
if (length(unknown))
  warning("source levels with no colour assigned: ", paste(unknown, collapse = ", "))

# pool rare STs, keep the common ones ordered by total count
st_n   <- dt[, .N, by = ST][order(-N)]
keep   <- st_n[N >= min_n, ST]
dt[, ST_grp := fifelse(ST %chin% keep, ST, "Singletons")]
dt[, ST_grp := factor(ST_grp, levels = c(keep, "Singletons"))]

dt[, group_summary := factor(group_summary,
                             levels = rev(names(src_cols)[names(src_cols) %chin%
                                                          unique(dt$group_summary)]))]

cat("genomes:        ", nrow(dt), "\n")
cat("STs shown:      ", length(keep), " (n >= ", min_n, ")\n", sep = "")
cat("pooled as Other:", st_n[N < min_n, sum(N)], "genomes in",
    st_n[N < min_n, .N], "STs\n")


# Matrix for the stacked bars -------------------------------------------------

st_src <- dcast(dt[, .N, by = .(group_summary, ST_grp)],
                group_summary ~ ST_grp, value.var = "N", fill = 0, drop = FALSE)

st_mat <- as.matrix(st_src[, -1])
rownames(st_mat) <- as.character(st_src$group_summary)

bar_cols <- src_cols[rownames(st_mat)]


# Plot -------------------------------------------------------------------------
# 
# png(paste0(outdir_fig, "/st_by_environment.png"),
#     width = 180, height = 95, units = "mm", res = 600,
#     pointsize = 7, type = "cairo", family = "sans")
par(mar = c(7,5, 0.5, 9.5), xpd = TRUE)

bp <- barplot(st_mat,
              col = bar_cols, border = "white", space = 0.25,
              las = 2, xaxt = "n", yaxt = "n",
              ylim = c(0, ceiling(max(colSums(st_mat)) / 10) * 10),
              xlab = "", ylab = "")

axis(side = 1, at = bp, labels = FALSE, lwd = 0, lwd.ticks = 0.75)

text(x = bp, y = -10,
     labels = colnames(st_mat), srt = 45, adj = 1,
     xpd = NA, cex = 0.75)

axis(side = 1, at = mean(bp), tick = FALSE, line = 2.5,
     labels = "Sequence type")

axis(side = 2, at = pretty(c(0, max(colSums(st_mat))-10)), las = 1)

axis(side = 2, at = max(colSums(st_mat)) / 2, tick = FALSE, line = 2,
     labels = "Number of genomes")

legend(34,82, inset = FALSE,
       legend = rev(rownames(st_mat)),
       fill   = rev(bar_cols),
       border = "white",
       bty = "n", cex = 0.75, ncol = 1)

# dev.off()
