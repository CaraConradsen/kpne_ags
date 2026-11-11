
# Where does selection sit within a co-occurrence cluster? ----------------
#
# Goldfinder clusters are only interpretable as cassettes if their members
# are physically linked. Clusters confined to one syntenic location are
# therefore drawn separately (left panel) from those whose members are
# scattered across the chromosome (right panel), because the two carry
# different claims: the first is co-selection of a linked unit, the second
# co-occurrence between unlinked loci.
#
# Within a cluster every member is counted on the side of the neutral
# expectation it falls on: fewer synonymous sites than expected to the left,
# more than expected to the right. Colour marks whether the departure passed
# the FDR threshold. A cassette under uniform selection gives a one-sided
# row mixing grey with one colour; genuinely variable selection within a
# cassette gives a row that straddles the line.
#
# Genes with no neutrality test carry no direction and are reported in the
# margin rather than forced onto one side.

# Set make_fig <- FALSE before sourcing to build syn_cnt, non_cnt and
# draw_panel() without opening a device, so the panels can be composed into
# a larger figure elsewhere (see analysis/ms_figures/).
if (!exists("make_fig")) make_fig <- TRUE

# Cluster IDs to print in bold, for clusters shown in detail elsewhere in
# the figure. Set before sourcing.
if (!exists("highlight_clusters")) highlight_clusters <- character(0)

q_cut         <- 0.05
min_outliers  <- 1          # 0 draws every cluster, not just those with an outlier
min_label_pc  <- 25         # % of annotated members needed to name a cluster
row_shade     <- TRUE
grp_gap       <- 1.3        # blank rows opened at each pangenome-class break
show_untested <- "margin"   # "margin" or "bar"
lab_trunc     <- 38

# Genes highlighted in the manuscript. A cluster containing any of these is
# annotated with them; otherwise the first few named outlier members are
# used, so every cluster gets whatever handle its annotation allows.
highlight <- c("tdh", "abiEi", "abiEii", "chrA", "wzm", "wbbM", "wbbM1",
               "rfbA", "rfbB", "cas1e", "cas2e", "cas3", "cas5e", "casA",
               "casB", "cse1", "cse2", "vgrG", "tssA", "tssF", "tssG",
               "vasK", "ccmA", "fimA", "fimD", "ureD", "traD", "traC",
               "trbL", "repA", "emrA", "emrB", "mrkC")
max_genes     <- 3

# BAKTA leaves 760 of the 1,204 outliers unnamed and labels 14 of them
# "unitary", so a name has to be recovered from the product string or set
# by hand. These are the families whose manuscript name differs from, or is
# absent from, the annotation.
gene_alias <- c(g003299 = "abiEi",  g001254 = "abiEii", g004051 = "wzm",
                g000322 = "wbbM",   g003974 = "tdh",    g001058 = "chrA")

junk_names <- c("", "NA", "unitary", "unknown")

sky  <- "#3FA9F5"           # excess diversity, as used throughout the ms
fire <- "#F26C2A"           # reduced diversity

# Inputs ------------------------------------------------------------------

clean_coin_ags <- fread(paste0(outdir_dat, "/clean_coin_ags.csv"))

# Read every row rather than unique(), because the modal location can only
# be found with the per-genome counts. A family annotated at more than one
# junction would otherwise duplicate in the merge and inflate the counts.
pan_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                  select = c("gene_family", "ag_type", "msu", "anchor",
                             "consensus_gene_name", "consensus_product"))

pan_anno[, syn_loc := paste(msu, anchor, sep = "_")]

loc_n <- pan_anno[, .N, by = .(gene_family, ag_type, msu, anchor, syn_loc,
                               consensus_gene_name, consensus_product)]
setorder(loc_n, gene_family, -N)

cat(sprintf("%d gene families annotated at >1 syntenic location; taking the modal location\n",
            loc_n[, .N, by = gene_family][N > 1, .N]))

pan_anno <- loc_n[, .SD[1], by = gene_family,
                  .SDcols = c("ag_type", "msu", "anchor", "syn_loc",
                              "consensus_gene_name", "consensus_product")]

clean_coin_ags <- merge(clean_coin_ags, pan_anno, all.x = TRUE, by = "gene_family")

clean_coin_ags[, anno_txt := paste(fifelse(is.na(consensus_gene_name), "",
                                           consensus_gene_name),
                                   consensus_product)]

# Syntenic or not ---------------------------------------------------------

clust_syn <- clean_coin_ags[, .(n_unique = uniqueN(syn_loc)), by = clust_id]

cat(nrow(clust_syn[n_unique > 1, ]), "; clusters with >1 syntenic loci\n",
    nrow(clust_syn[n_unique == 1, ]), "; clusters with 1 syntenic loci\n",
    "out of ", nrow(clust_syn), " clusters\n")

clean_coin_ags[clust_syn, on = .(clust_id),
               panel := fifelse(i.n_unique == 1, "syntenic", "non-syntenic")]

clean_coin_ags[!call %chin% c("high", "ns", "low"), call := "NA"]

# What is the cluster? ----------------------------------------------------
# A keyword table rather than the modal product, because members of a real
# cassette usually carry different products that point at one process. The
# rules are deliberately visible and editable; a cluster is named only when
# at least min_label_pc of its annotated members match, otherwise it stays
# "Mixed". Patterns are matched case-insensitively, so single-letter gene
# suffixes are written as explicit character classes rather than [A-Z].

rules <- data.table(
  label = c("CRISPR-Cas", "T6SS", "Prophage", "O-antigen",
            "Fimbriae", "Conjugation", "Defence/TA",
            "Siderophore", "Sugar metab.", "Transport",
            "Secretion", "Regulation"),
  pattern = c(
    "CRISPR|\\bcas[0-9]\\b|\\bcse[0-9]\\b",
    "type VI|T6SS|VgrG|\\bHcp\\b|PAAR|\\btss[a-z]\\b|VasK|IcmF",
    "phage|capsid|portal|baseplate|terminase|\\btail\\b|tape measure|holin|\\blysin",
    paste0("O-antigen|O[0-9]+ family|\\brfb[a-z]?\\b|\\bwb[a-z]{1,2}\\b|",
           "lipopolysacchar|glycosyltransferase|mannosyltransferase|undecaprenyl|\\bwz[a-z]\\b"),
    "fimbri|pilus|pilin|adhesin|usher|\\bmrk[a-z]\\b",
    "conjug|\\btra[bcdegijklmnuvwy]\\b|\\btrb[a-l]\\b|mobilis|relaxase|\\brepa\\b|oriT|type IV secretion",
    "restriction|abortive|\\babi[a-z]{1,2}\\b|toxin-antitoxin|antitoxin|immunity|retron",
    "enterobactin|siderophore|ferric|ferrous|\\bent[a-z]\\b|yersiniabactin|aerobactin",
    paste0("\\bPTS\\b|sugar|galactitol|mannose|fructose|glucos|transketolase|aldolase|",
           "arabinose|xylose|rhamnose|sorbitol|galactos|alcohol dehydrogenase|phosphomannomutase"),
    "efflux|major facilitator|\\bMFS\\b|ABC transporter|permease|\\btransporter\\b",
    "secretion|lipoprotein|exported",
    "transcriptional regulator|LysR|GntR|AraC|helix-turn-helix|sigma factor|two-component"))

# Genes to print beside the label: manuscript highlights first, then any
# other named outlier in the cluster. Deficit outliers are flagged with *.
cluster_genes <- function(fam, nm, prod, call) {

  keep <- call %chin% c("high", "low")
  fam <- fam[keep]; nm <- trimws(nm[keep]); prod <- prod[keep]; call <- call[keep]
  if (!length(fam)) return("")

  # a manuscript alias wins; then a highlight term appearing in the product
  # (which is how abiEi and abiEii are recoverable at all); then the
  # annotation's own name, if it is a name rather than a placeholder
  disp <- vapply(seq_along(fam), function(i) {
    if (fam[i] %chin% names(gene_alias)) return(unname(gene_alias[fam[i]]))
    j <- which(vapply(highlight,
                      function(hh) grepl(paste0("\\b", hh, "\\b"), prod[i],
                                         ignore.case = TRUE), logical(1)))
    if (length(j)) return(highlight[j[1]])
    if (tolower(nm[i]) %chin% junk_names) return(NA_character_)
    nm[i]
  }, character(1))

  ok <- !is.na(disp)
  if (!any(ok)) return("")
  disp <- disp[ok]; call <- call[ok]

  tag <- paste0(disp, fifelse(call == "low", "*", ""))
  hit <- toupper(disp) %chin% toupper(c(highlight, names(gene_alias)))
  tag <- unique(c(tag[hit], tag[!hit]))

  paste0(" (", paste(head(tag, max_genes), collapse = ", "),
         if (length(tag) > max_genes) ", ..." else "", ")")
}

# Cutting a label at a fixed character count leaves half-words and orphan
# gene names, so cut back to the last complete token and close the bracket.
trunc_lab <- function(x, n) {
  vapply(x, function(z) {
    if (nchar(z) <= n) return(z)
    cut <- substr(z, 1, n)
    sp  <- regexpr("[^ ,(]*$", cut)
    if (sp > 1L) cut <- substr(cut, 1, sp - 1L)
    cut <- sub("[ ,]+$", "", cut)
    if (grepl("\\(", cut) && !grepl("\\)", cut)) paste0(cut, ", ...)") else
                                                     paste0(cut, "...")
  }, character(1), USE.NAMES = FALSE)
}

classify_cluster <- function(txt) {
  txt <- txt[nzchar(trimws(txt)) & !grepl("ypothetical", txt)]
  if (!length(txt)) return("Mixed")
  hits <- vapply(rules$pattern,
                 function(p) sum(grepl(p, txt, ignore.case = TRUE)), integer(1))
  if (max(hits) == 0L) return("Mixed")
  j  <- which.max(hits)
  pc <- 100 * hits[j] / length(txt)
  if (pc >= min_label_pc) rules$label[j] else "Mixed"
}

# Per-cluster counts ------------------------------------------------------

build_counts <- function(dt) {

  dt <- copy(dt)

  # "ns" genes still have a direction, and that direction is what makes the
  # diverging layout meaningful, so they are split rather than pooled.
  dt[, side := fcase(
    call == "low",                           "low",
    call == "high",                          "high",
    call == "ns" & seg_sites <  expected_S,  "ns_low",
    call == "ns" & seg_sites >= expected_S,  "ns_high",
    default =                                "untested")]

  cnt <- dcast(dt, clust_id ~ side, fun.aggregate = length, value.var = "gene_family")
  for (s in c("low", "ns_low", "ns_high", "high", "untested")) {
    if (!s %chin% names(cnt)) cnt[, (s) := 0L]
  }

  # Members of a cluster co-occur by construction, so their frequencies are
  # nearly identical and most clusters carry a single ag_type; a median k
  # and a modal class are therefore well defined rather than an average over
  # a mixture.
  meta <- dt[, .(k_med = as.numeric(median(freq, na.rm = TRUE)),
                 grp   = names(sort(table(ag_type), decreasing = TRUE))[1],
                 what  = paste0(classify_cluster(anno_txt),
                                cluster_genes(gene_family, consensus_gene_name,
                                              consensus_product, call))),
             by = clust_id]

  cnt <- meta[cnt, on = .(clust_id)]

  cnt[, `:=`(n_gene    = low + ns_low + ns_high + high + untested,
             n_outlier = low + high,
             left      = low + ns_low,
             right     = ns_high + high + untested,
             net       = (ns_high + high) - (ns_low + low))]

  cnt[, grp := factor(grp, levels = c("cloud", "shell", "soft"))]

  # cloud first, and within each class from the rows favouring the left to
  # those favouring the right
  cnt[n_outlier >= min_outliers][order(grp, net, n_gene)]
}

syn_cnt <- build_counts(clean_coin_ags[panel == "syntenic"])
non_cnt <- build_counts(clean_coin_ags[panel == "non-syntenic"])

cat(sprintf("drawing %d syntenic and %d non-syntenic clusters with >= %d outlier(s)\n",
            nrow(syn_cnt), nrow(non_cnt), min_outliers))

fwrite(rbind(syn_cnt[, panel := "syntenic"], non_cnt[, panel := "non-syntenic"]),
       paste0(outdir_dat, "/goldfinder_cluster_side_counts.csv"))

# Figure ------------------------------------------------------------------

grp_lab <- c(cloud = "Cloud", shell = "Shell", soft = "Soft core")

seg_left  <- c("ns_low", "low")                 # inner -> outer
seg_right <- if (show_untested == "bar") c("ns_high", "high", "untested") else
                                         c("ns_high", "high")

fill <- c(ns_low = "grey78", low = fire, ns_high = "grey78", high = sky,
          untested = "grey93")
brd  <- c(ns_low = "white", low = "white", ns_high = "white", high = "white",
          untested = "grey60")

# one scale for both panels so the widths are comparable
xr <- function(d) if (show_untested == "bar") max(d$right) else max(d$ns_high + d$high)
xmax <- max(xr(syn_cnt), xr(non_cnt))
xmin <- -max(syn_cnt$left, non_cnt$left)
pad  <- 0.06 * (xmax - xmin)
nmax <- max(nrow(syn_cnt), nrow(non_cnt))

draw_panel <- function(cnt, main, n_slots = nmax) {

  n  <- nrow(cnt)
  bh <- 0.38

  # A gap is opened at each class break and the class name sits in it, so
  # the labels read horizontally on the plot rather than rotated in the
  # margin. Row positions are counted down from a top that is shared across
  # panels, so bar heights stay identical between a and b.
  brk  <- c(0L, as.integer(head(as.character(cnt$grp), -1) !=
                           tail(as.character(cnt$grp), -1)))
  y_of <- (n_slots + 2 * grp_gap) - (seq_len(n) - 1) - cumsum(brk) * grp_gap

  col_k <- xmax + 1.6 * pad
  col_n <- xmax + 2.4 * pad
  col_w <- xmax + 5.2 * pad

  par(mar = c(4.2, 3.8, 6.2, 14), bty = "n", xaxs = "i")

  plot(NA, xlim = c(xmin - pad, xmax + pad),
       ylim = c(0.4, n_slots + 2 * grp_gap + 1),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n")

  if (row_shade) {
    for (i in seq_len(n)) {
      if (i %% 2 == 0) next
      y <- y_of[i]
      rect(par("usr")[1], y - 0.5, par("usr")[2], y + 0.5, col = "grey97",
           border = NA)
    }
  }

  grp_run <- rle(as.character(cnt$grp))
  g_end   <- cumsum(grp_run$lengths)
  g_start <- g_end - grp_run$lengths + 1

  for (j in seq_along(grp_run$values)) {
    y_top <- y_of[g_start[j]]
    if (j > 1) {
      segments(par("usr")[1], y_top + grp_gap * 0.78,
               par("usr")[2], y_top + grp_gap * 0.78, col = "grey65", lwd = 0.8)
    }
    text(xmin - pad, y_top + 0.60, labels = grp_lab[grp_run$values[j]],
         adj = c(0, 0), font = 2, cex = 0.7, xpd = NA)
  }

  abline(v = 0, lty = 2, col = "grey30")

  for (i in seq_len(n)) {

    y <- y_of[i]

    x <- 0
    for (s in seg_left) {
      w <- cnt[[s]][i]
      if (w > 0) {
        rect(x - w, y - bh, x, y + bh, col = fill[s], border = brd[s], lwd = 0.6)
        x <- x - w
      }
    }

    x <- 0
    for (s in seg_right) {
      w <- cnt[[s]][i]
      if (w > 0) {
        rect(x, y - bh, x + w, y + bh, col = fill[s], border = brd[s], lwd = 0.6)
        x <- x + w
      }
    }

    nlab <- if (show_untested == "margin" && cnt$untested[i] > 0)
              sprintf("%d (%d)", cnt$n_gene[i], cnt$untested[i]) else
              as.character(cnt$n_gene[i])

    text(col_k, y, round(cnt$k_med[i]), adj = 1, cex = 0.55, col = "grey35", xpd = NA)
    text(col_n, y, nlab, adj = 0, cex = 0.55, col = "grey35", xpd = NA)
    text(col_w, y, trunc_lab(cnt$what[i], lab_trunc), adj = 0, cex = 0.55,
         col = "grey20", xpd = NA)
  }

  at <- pretty(c(xmin, xmax))
  axis(1, at = at, labels = abs(at), cex.axis = 0.78)
  # drawn with text() rather than axis() so that the clusters shown in
  # detail in panels c and d can be picked out in bold
  hl <- as.character(cnt$clust_id) %chin% as.character(highlight_clusters)
  text(xmin - pad - 0.03 * (xmax - xmin), y_of,
       labels = cnt$clust_id, adj = 1, cex = 0.6,
       font = fifelse(hl, 2L, 1L), xpd = NA)

  mtext("Genes per cluster", side = 1, line = 3.3, cex = 0.9)
  mtext("fewer sites than expected", side = 1, line = 2.1, adj = 0, cex = 0.62,
        col = "grey40")
  mtext("more than expected", side = 1, line = 2.1, adj = 1, cex = 0.62,
        col = "grey40")
  mtext("Co-occurrence cluster", side = 2, line = 2.4, cex = 0.9)

  mtext(expression(italic(k)), side = 3, line = -0.2, at = col_k, adj = 1,
        cex = 0.55, col = "grey35")
  mtext(if (show_untested == "margin") "n (unt.)" else "n", side = 3,
        line = -0.2, at = col_n, adj = 0, cex = 0.55, col = "grey35")
  mtext("function (genes; * reduced)", side = 3,
        line = -0.2, at = col_w, adj = 0, cex = 0.55, col = "grey35")

  # anchored to the left edge of the figure region rather than to the plot
  # box, so the tag lines up with the panels below, which have different
  # margins
  mtext(main, side = 3, line = 4.4, adj = 0, font = 2, cex = 0.95,
        at = grconvertX(0.01, "nfc", "user"), xpd = NA)
}

if (isTRUE(make_fig)) {

  png(paste0(outdir_fig, "/goldfinder_cluster_sides.png"),
      width = 32, height = 3.6 + 0.40 * nmax, units = "cm", res = 300,
      pointsize = 12, type = "cairo")

  layout(matrix(1:2, nrow = 1))

  draw_panel(syn_cnt, sprintf("a  Single syntenic locus (n = %d)", nrow(syn_cnt)))

  legend(x = xmin - pad, y = nmax + 2 * grp_gap + 6.4, horiz = TRUE, bty = "n",
         fill = c(fire, "grey78", sky), border = "white", cex = 0.66,
         legend = c(paste0("Reduced diversity (q < ", q_cut, ")"),
                    "Not significant",
                    paste0("Excess diversity (q < ", q_cut, ")")),
         xpd = NA)

  draw_panel(non_cnt, sprintf("b  Multiple loci (n = %d)", nrow(non_cnt)))

  dev.off()

  cat("wrote", paste0(outdir_fig, "/goldfinder_cluster_sides.png"), "\n")

}

# Summary -----------------------------------------------------------------

for (p in c("syntenic", "non-syntenic")) {
  d <- if (p == "syntenic") syn_cnt else non_cnt
  cat(sprintf("\n%s: %d clusters, %d one-sided, %d straddle, %d hold both tails\n",
              p, nrow(d),
              d[left == 0 | (right - untested) == 0, .N],
              d[left > 0 & (right - untested) > 0, .N],
              d[low > 0 & high > 0, .N]))
  print(d[, .(clusters = .N, median_k = as.numeric(median(k_med))), by = .(grp)])
  print(d[, .N, by = what][order(-N)])
}
