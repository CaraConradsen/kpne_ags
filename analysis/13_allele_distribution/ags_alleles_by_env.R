# Peptide alleles of one accessory gene, by ST and by isolation source
#
# One column, four rows: ST on top, isolation source beneath. Rows 1 and 3
# are counts, so the sampling imbalance is on the page; rows 2 and 4 are the
# allele-frequency panels of allelic_diversity.R, bars scaled to 1 with width
# proportional to sqrt(carriers) and the three major allele frequencies named
# in the header, in the colour of the allele each refers to.
#
# Every source level is shown. Rare STs are still pooled as "other", as in
# 10_outlier_function/ags_alleles_by_ST_plot.R, because there are 74 of them.
#
# Source colours follow Fig. 1 of Thorpe et al. 2022 (Nat Microbiol
# 7:2054-2067) via plot_st_by_environment.R. The bars themselves are coloured
# by allele, since the same alleles must be comparable across all four
# panels, so the source colours are carried on the source labels instead.
#
# Description only. Nothing here tests association between allele and source.

make_fig <- FALSE                      # the standalone ST figure is not wanted
source("./analysis/13_allele_distribution/ags_alleles_by_ST_plot.R")

# Source colours ------------------------------------------------------------
# Duplicated from analysis/04_assemblying_datasets/plot_st_by_environment.R
# rather than sourced, because that script draws a figure and reads the
# metadata on source, which would cut across the device opened below. If the
# hex values are corrected there, correct them here too.

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

unknown <- setdiff(unique(al$env), names(src_cols))
if (length(unknown))
  warning("source levels with no colour assigned: ", paste(unknown, collapse = ", "))

# allele x source, every level kept, in the order Thorpe et al. Fig. 1 uses
env_lev <- names(src_cols)[names(src_cols) %chin% unique(al$env)]
al[, env_g := factor(env, levels = c(env_lev, setdiff(unique(al$env), env_lev)))]

M_env <- as.matrix(dcast(al[, .N, by = .(allele, env_g)], allele ~ env_g,
                         value.var = "N", fill = 0, drop = FALSE),
                   rownames = "allele")
M_env <- M_env[rowSums(M_env) > 0, colSums(M_env) > 0, drop = FALSE]

print(data.table(source = colnames(M_env), carriers = colSums(M_env),
                 alleles = colSums(M_env > 0)))

# Figure ---------------------------------------------------------------------

png(paste0(outdir_fig, "/alleles_by_ST_env_", fam, ".png"),
    width = 18, height = 30, units = "cm", res = 600,
    pointsize = 11, type = "cairo")

layout(matrix(1:4, ncol = 1), heights = c(1, 1.1, 1.3, 1.35))

# a, b: ST on the x-axis, alleles stacked
par(mar = c(3.6, 4.8, 4.8, 1.4), bty = "l")
plot_allele_grp(M_st, prop = FALSE, hdr = paste0(gene_lab, " (", fam, "), by ST"))
panel_tag("a")

par(mar = c(5.2, 4.8, 4.0, 1.4))
plot_allele_grp(M_st, prop = TRUE)
panel_tag("b", dy = 0.22)

# c, d: isolation source on the x-axis, alleles stacked as in a and b. The
# source categories carry the Thorpe et al. colours on their labels.
par(mar = c(7.8, 4.8, 4.8, 1.4))
plot_allele_grp(M_env, prop = FALSE, lab_cols = src_cols, hdr_stats = FALSE,
                hdr = paste0(gene_lab, " (", fam, "), by isolation source"))
panel_tag("c")

par(mar = c(9.5, 4.8, 4.0, 1.4))
plot_allele_grp(M_env, prop = TRUE, lab_cols = src_cols)
panel_tag("d", dy = 0.22)

dev.off()
cat("wrote", paste0(outdir_fig, "/alleles_by_ST_env_", fam, ".png"), "\n")
