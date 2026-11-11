## operon_ribbon_plots.R
## Summarising gene-ribbon (synteny Sankey) plots for selected accessory-gene operons.
## Base R only.
##
## Synteny comes from the pangraph MSU / anchor columns:
##   anchor "A238"       = the 238th core gene of that MSU
##   anchor "A238:A239"  = accessory junction between core genes 238 and 239
## For each operon the region is: core genes A[L-k+1..L], the junction A[L]:A[R]
## (and any junctions between the flanking core genes), core genes A[R..R+k-1],
## all on the same MSU and contig. Gene order within a genome is taken from
## start/end, and each genome is oriented so the left core gene (A[L]) comes first.
##
## Each column is a gene position aligned across genomes; each box is a gene
## family at that position (height = number of genomes); ribbons join
## consecutive genes within a genome (width = number of genomes with that
## adjacency). Boxes are coloured by the per-gene diversity test in pval_file.

## ---------------------------------------------------------------- settings
anno_file <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/all_pirate_anno_full.csv"
pval_file <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/gene_pvals_fdr.csv"
out_dir   <- "C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/figures"

q_alpha   <- 0.05   # BH threshold on p_bh
min_frac  <- 0.02   # families in < this fraction of genomes at a position are pooled as "other"
label_min <- 0.05   # only label boxes holding >= this fraction of genomes
mge_min   <- 0.5    # mark a box as MGE if >= this fraction of its genes carry an MGE call
max_ins_run <- 6    # insertion runs longer than this are drawn as one "+n genes" block

## MGE calls come from is_fam, integron, ice and phg_state (any non-empty value)
mge_tag <- c("Prophage" = "P", "IS element" = "IS", "ICE" = "ICE", "Integron" = "Int")
mge_col <- c("Prophage" = "#7B3FA0", "IS element" = "#1E8449", "ICE" = "#B7950B", "Integron" = "#C2185B")

## core_flank: number of core genes shown either side of the junction
## ag_flank:   genes kept either side of the target genes (NA = keep whole region);
##             give c(left, right) for different flanks on each side;
##             use for long junctions such as the capsule/O-antigen locus
## labels:     optional display names for gene families
operons <- list(
  list(name       = "g004645 locus",
       targets    = "g004645",
       core_flank = 2, ag_flank = NA,
       labels     = c(g004645 = "ompK36")),
  list(name       = "O1/O2 O-antigen operon",
       targets    = c("g004051", "g003783", "g000322", "g001147"),
       core_flank = 1, ag_flank = 4,
       labels     = c(g004051 = "wzm", g003783 = "wzt", g000322 = "wbbM", g001147 = "glf")),
  list(name       = "AbiE operon",
       targets    = c("g003299", "g001254"),
       core_flank = 2, ag_flank = NA,
       labels     = c(g003299 = "abiEi", g001254 = "abiEii")),
  list(name       = "T6SS (vgrG-tssA)",
       targets    = c("g003390", "g003914"),
       core_flank = 1, ag_flank = c(4, 32),   # downstream prophage in ~38 genomes
       labels     = c(g003390 = "vgrG", g003914 = "tssA")),
  list(name       = "Acetoin dehydrogenase (acoA-acoB)",
       targets    = c("g002655", "g003239"),
       core_flank = 1, ag_flank = 4,
       labels     = c(g002655 = "acoA", g003239 = "acoB",
                      g000332 = "acoC")),   # annotated "accB" but product is the E2 dihydrolipoamide acetyltransferase
  list(name       = "Sorbose PTS (sorF-sorC)",
       targets    = c("g003377", "g004320"),
       core_flank = 2, ag_flank = NA,
       labels     = c(g003377 = "sorF", g004320 = "sorC")),
  list(name       = "Chromate efflux (chrA)",
       targets    = "g001058",
       core_flank = 1, ag_flank = 6,
       labels     = c(g001058 = "chrA")),
  ## Two type I-E CRISPR-Cas systems, at different chromosomal positions and
  ## carried by mutually exclusive sets of genomes (30 and 15, no overlap).
  ## Locus A sits in a compact junction of 13 families with no MGE calls;
  ## locus B sits in a 283-family hotspot carrying an intact prophage in 50
  ## genomes, so its flanks are trimmed.
  list(name       = "CRISPR-Cas type I-E, locus A",
       targets    = c("g002579", "g003424", "g005291", "g004669",
                      "g004793", "g005762", "g006512", "g005961"),
       core_flank = 2, ag_flank = NA,
       labels     = c(g002579 = "cas3",  g003424 = "cse1/casA",
                      g005291 = "cse2/casB", g004669 = "cas7",
                      g004793 = "cas5",  g005762 = "cas6e",
                      g006512 = "cas1",  g005961 = "cas2")),
  list(name       = "CRISPR-Cas type I-E, locus B",
       targets    = c("g003179", "g004638", "g014756", "g005653",
                      "g005811", "g004752", "g005162", "g006892"),
       core_flank = 1, ag_flank = 5,
       labels     = c(g003179 = "cas3",  g004638 = "cse1/casA",
                      g014756 = "cse2/casB", g005653 = "cas7",
                      g005811 = "cas5",  g004752 = "cas6e",
                      g005162 = "cas1",  g006892 = "cas2"))
)

## ---------------------------------------------------------------- data
keep <- c("geno_id", "gene_family", "seqnames", "start", "end", "strand",
          "gene", "consensus_gene_name", "kap_gene", "ag_type", "anchor", "msu",
          "phg_state", "integron", "ice", "is_fam")
hdr <- names(read.csv(anno_file, nrows = 1, check.names = FALSE))
cc  <- rep("NULL", length(hdr))
cc[match(keep, hdr)] <- "character"
cc[match(c("start", "end"), hdr)] <- "integer"
anno <- read.csv(anno_file, colClasses = cc, na.strings = c("", "NA"))

## numeric anchor position: "A238" -> 238, "A238:A239" -> 238.5
anchor_pos <- function(a) {
  v <- regmatches(a, gregexpr("[0-9]+", a))
  vapply(v, function(x) {
    x <- as.numeric(x)
    if (length(x) == 1) x
    else if (length(x) == 2 && abs(diff(x)) == 1) mean(x)
    else NA_real_                       # wrap-around / malformed junctions
  }, numeric(1))
}
ua <- unique(anno$anchor)
anno$apos <- anchor_pos(ua)[match(anno$anchor, ua)]

pv <- read.csv(pval_file)
pv$sel <- ifelse(is.na(pv$p_bh), "Not in test set",
          ifelse(pv$p_bh >= q_alpha, "Not significant",
          ifelse(pv$p_upper < pv$p_lower, "Excess diversity", "Deficit diversity")))
sel_col <- c("Excess diversity"  = "#3FA9F5",
             "Deficit diversity" = "#F26C2A",
             "Not significant"   = "#C9CDD2",
             "Not in test set"        = "#FFFFFF",
             "Other (pooled)"    = "#EDEDED")

## display label per family: user label > kap_gene > gene > consensus name > id
first_mode <- function(x) { x <- x[!is.na(x)]; if (length(x)) names(sort(table(x), decreasing = TRUE))[1] else NA }
fam_label <- function(fams, d, user) {
  out <- vapply(fams, function(f) {
    if (f %in% names(user)) return(user[[f]])
    s <- d[d$gene_family == f, ]
    for (col in c("kap_gene", "gene", "consensus_gene_name")) {
      m <- first_mode(s[[col]]); if (!is.na(m)) return(m)
    }
    f
  }, character(1))
  ifelse(nchar(out) > 12, paste0(substr(out, 1, 11), "."), out)
}

## ---------------------------------------------------------------- build one operon
build_operon <- function(op) {
  tg  <- anno[anno$gene_family %in% op$targets, ]
  tab <- sort(table(paste(tg$msu, tg$anchor, sep = "|")), decreasing = TRUE)
  key <- strsplit(names(tab)[1], "|", fixed = TRUE)[[1]]
  msu <- key[1]; anc <- key[2]
  if (length(tab) > 1)
    message(op$name, ": targets also found at ", paste(names(tab)[-1], tab[-1], sep = " (n=", collapse = "), "), ")")
  jpos <- anchor_pos(anc)
  L <- floor(jpos); R <- ceiling(jpos)
  k <- op$core_flank
  lo <- if (k > 0) L - (k - 1) else jpos
  hi <- if (k > 0) R + (k - 1) else jpos

  reg <- anno[anno$msu %in% msu & !is.na(anno$apos) & anno$apos >= lo & anno$apos <= hi, ]

  ## per-genome oriented gene order
  gl <- lapply(split(reg, reg$geno_id), function(d) {
    d <- d[d$seqnames == names(which.max(table(d$seqnames))), ]
    d <- d[order(d$start, d$end), ]
    l <- d$start[d$apos < jpos]; r <- d$start[d$apos > jpos]
    flip <- length(l) && length(r) && median(l) > median(r)
    if (flip) {
      d <- d[nrow(d):1, ]
      d$strand <- ifelse(d$strand == "+", "-", "+")
    }
    ## MGE call per gene (most specific first)
    mge <- ifelse(!is.na(d$is_fam), "IS element",
           ifelse(!is.na(d$integron), "Integron",
           ifelse(!is.na(d$ice), "ICE",
           ifelse(!is.na(d$phg_state), "Prophage", NA))))
    data.frame(geno = d$geno_id[1], fam = d$gene_family, strand = d$strand,
               idx = seq_len(nrow(d)), mge = mge, stringsAsFactors = FALSE)
  })

  ## reference = the most common complete gene-family order among genomes
  ## carrying the most targets; reference column 0 = its first target gene
  ntg  <- vapply(gl, function(g) sum(unique(g$fam) %in% op$targets), numeric(1))
  sigs <- vapply(gl[ntg == max(ntg)], function(g) paste(g$fam, collapse = ";"), character(1))
  ref  <- strsplit(names(which.max(table(sigs))), ";", fixed = TRUE)[[1]]
  rcol <- seq_along(ref) - min(which(ref %in% op$targets))
  tcols <- rcol[ref %in% op$targets]

  ## global (Needleman-Wunsch) alignment of each genome's gene-family order to
  ## the reference. A mismatch puts an alternative family in the reference
  ## slot; an insertion is placed between reference slots.
  nw <- function(a, b, ms = 2, mm = -1, gp = -1) {
    n <- length(a); m <- length(b)
    S <- matrix(0, n + 1, m + 1); S[, 1] <- gp * 0:n; S[1, ] <- gp * 0:m
    for (i in seq_len(n)) for (j in seq_len(m))
      S[i + 1, j + 1] <- max(S[i, j] + if (a[i] == b[j]) ms else mm,
                             S[i, j + 1] + gp, S[i + 1, j] + gp)
    map <- rep(NA_integer_, n); i <- n; j <- m
    while (i > 0 && j > 0) {
      if (S[i + 1, j + 1] == S[i, j] + if (a[i] == b[j]) ms else mm) { map[i] <- j; i <- i - 1; j <- j - 1 }
      else if (S[i + 1, j + 1] == S[i, j + 1] + gp) i <- i - 1
      else j <- j - 1
    }
    map
  }
  gl <- lapply(gl, function(g) {
    ## align on shared families only (mismatch costlier than an indel), then
    ## treat a run of unmatched genes as a replacement of the reference slots
    ## between the same two anchors only when the run lengths are equal
    ## (e.g. ompC for ompK36); otherwise they stay insertions
    map <- nw(g$fam, ref, mm = -3)
    i <- 1; n <- length(map)
    while (i <= n) {
      if (is.na(map[i])) {
        i2 <- i; while (i2 < n && is.na(map[i2 + 1])) i2 <- i2 + 1
        pj <- if (i > 1) map[i - 1] else 0
        nj <- if (i2 < n) map[i2 + 1] else length(ref) + 1
        if (nj - pj - 1 == i2 - i + 1) map[i:i2] <- (pj + 1):(nj - 1)
        i <- i2 + 1
      } else i <- i + 1
    }
    ## slot = reference column; insertions get the column of the preceding
    ## aligned gene (-Inf before the first) plus an insertion rank
    prev <- -Inf; k <- 0
    g$rc <- NA_real_; g$ins <- 0L
    for (i in seq_len(nrow(g))) {
      if (!is.na(map[i])) { prev <- rcol[map[i]]; k <- 0; g$rc[i] <- prev }
      else { k <- k + 1; g$rc[i] <- prev; g$ins[i] <- k }
    }
    g$rc[is.infinite(g$rc)] <- min(rcol) - 1
    ## collapse insertion runs longer than max_ins_run into one block
    g$run_len <- 1L
    if (any(g$ins > 0)) {
      rid <- cumsum(g$ins == 1 | g$ins == 0)       # new id at each run start / aligned gene
      for (r in unique(rid[g$ins > 0])) {
        w <- which(rid == r & g$ins > 0)
        if (length(w) > max_ins_run) {
          m <- g$mge[w]
          g$fam[w[1]] <- "ins_run"
          g$run_len[w[1]] <- length(w)
          g$mge[w[1]] <- if (mean(!is.na(m)) >= mge_min) names(which.max(table(m))) else NA
          g <- g[-w[-1], ]; rid <- rid[-w[-1]]
        }
      }
    }
    g
  })
  al <- do.call(rbind, gl)

  if (!all(is.na(op$ag_flank))) {
    fl  <- rep(op$ag_flank, length.out = 2)
    wlo <- min(tcols) - fl[1]
    whi <- max(tcols) + fl[2]
    al  <- al[al$rc >= wlo & (al$rc < whi | (al$rc == whi & al$ins == 0)), ]
  }
  ## integer plotting columns: reference slots plus insertion slots, in order
  slots <- unique(al[, c("rc", "ins")])
  slots <- slots[order(slots$rc, slots$ins), ]
  slots$col <- seq_len(nrow(slots))
  al$col <- slots$col[match(paste(al$rc, al$ins), paste(slots$rc, slots$ins))]
  col_lab <- setNames(ifelse(slots$ins == 0, as.character(slots$rc), ""), slots$col)
  N <- length(gl)

  ## pool rare families per column
  cnt <- table(paste(al$col, al$fam))
  rare <- cnt[paste(al$col, al$fam)] < max(2, ceiling(min_frac * N)) &
          !al$fam %in% c(op$targets, "ins_run")
  al$node_fam <- ifelse(rare, "other", al$fam)
  al <- al[order(al$geno, al$col), ]

  list(op = op, msu = msu, anchor = anc, L = L, R = R, N = N,
       n_with = sum(ntg > 0), al = al, reg = reg, col_lab = col_lab)
}

## ---------------------------------------------------------------- layout + draw
plot_ribbon <- function(b, box_w = 0.62, gap = 0.025) {
  al <- b$al; N <- b$N; op <- b$op
  al$node <- paste(al$col, al$node_fam, sep = "|")

  ## nodes
  nd <- aggregate(geno ~ node + col + node_fam, data = al, FUN = length)
  names(nd)[names(nd) == "geno"] <- "n"
  nd$plus <- tapply(al$strand == "+", al$node, mean)[nd$node]
  ## node MGE class: the most common MGE call, if >= mge_min of its genes carry one
  mg <- tapply(al$mge, al$node, function(m) {
    if (mean(!is.na(m)) < mge_min) NA_character_ else names(which.max(table(m)))
  })
  nd$mge <- mg[nd$node]
  nd$h <- nd$n / N

  ## edges between consecutive genes in a genome
  sp <- split(al, al$geno)
  ed <- do.call(rbind, lapply(sp, function(g) {
    if (nrow(g) < 2) return(NULL)
    data.frame(from = g$node[-nrow(g)], to = g$node[-1], stringsAsFactors = FALSE)
  }))
  ed <- aggregate(list(n = rep(1, nrow(ed))), by = ed, FUN = sum)

  cols <- sort(unique(nd$col))
  ## initial order: targets, then by count, "other" last
  nd$ord <- ifelse(nd$node_fam %in% op$targets, 0, ifelse(nd$node_fam == "other", 2, 1))
  place <- function(nd) {
    nd$y0 <- NA; nd$y1 <- NA
    for (cc in cols) {
      i <- which(nd$col == cc); i <- i[order(nd$key[i])]
      tot <- sum(nd$h[i]) + gap * (length(i) - 1)
      top <- 0.5 + tot / 2
      for (j in i) { nd$y1[j] <- top; nd$y0[j] <- top - nd$h[j]; top <- nd$y0[j] - gap }
    }
    nd
  }
  nd$key <- nd$ord * 10 - nd$h
  nd <- place(nd)
  ## barycentre sweeps to reduce ribbon crossings (targets kept on top)
  for (it in 1:4) {
    yc <- setNames((nd$y0 + nd$y1) / 2, nd$node)
    for (cc in (if (it %% 2) cols else rev(cols))) {
      i <- which(nd$col == cc)
      nb <- vapply(nd$node[i], function(v) {
        e <- if (it %% 2) ed[ed$to == v, ] else ed[ed$from == v, ]
        o <- if (it %% 2) e$from else e$to
        if (nrow(e)) sum(yc[o] * e$n) / sum(e$n) else 0.5
      }, numeric(1))
      nd$key[i] <- ifelse(nd$node_fam[i] %in% op$targets, -10 - nb,
                   ifelse(nd$node_fam[i] == "other", 10, -nb))
    }
    nd <- place(nd)
  }
  rownames(nd) <- nd$node

  ## ribbon slots
  ed$ymid_to   <- (nd[ed$to, "y0"] + nd[ed$to, "y1"]) / 2
  ed$ymid_from <- (nd[ed$from, "y0"] + nd[ed$from, "y1"]) / 2
  ed$s1 <- ed$s0 <- ed$t1 <- ed$t0 <- NA
  for (v in nd$node) {
    o <- which(ed$from == v); o <- o[order(-ed$ymid_to[o])]
    top <- nd[v, "y1"]
    for (j in o) { ed$s1[j] <- top; ed$s0[j] <- top - ed$n[j] / N; top <- ed$s0[j] }
    i <- which(ed$to == v); i <- i[order(-ed$ymid_from[i])]
    top <- nd[v, "y1"]
    for (j in i) { ed$t1[j] <- top; ed$t0[j] <- top - ed$n[j] / N; top <- ed$t0[j] }
  }

  ## selection category per node
  fam_sel <- pv$sel[match(nd$node_fam, pv$gene_family)]
  fam_sel[is.na(fam_sel)] <- "Not in test set"
  fam_sel[nd$node_fam %in% c("other", "ins_run")] <- "Other (pooled)"
  nd$fill <- sel_col[fam_sel]
  rl <- tapply(al$run_len, al$node, function(x)
    if (min(x) == max(x)) as.character(x[1]) else paste0(min(x), "-", max(x)))
  nd$run_lab <- rl[nd$node]

  real <- unique(nd$node_fam[!nd$node_fam %in% c("other", "ins_run")])
  lab  <- setNames(fam_label(real, b$reg, op$labels), real)

  ## canvas
  xr <- range(cols) + c(-0.8, 0.8)
  yr <- range(c(nd$y0, nd$y1)) + c(-0.06, 0.04)
  plot.new(); plot.window(xlim = xr, ylim = yr, xaxs = "i", yaxs = "i")

  ## ribbons
  tt <- seq(0, 1, length.out = 40); sm <- 3 * tt^2 - 2 * tt^3
  for (j in seq_len(nrow(ed))) {
    x0 <- nd[ed$from[j], "col"] + box_w / 2; x1 <- nd[ed$to[j], "col"] - box_w / 2
    xs <- x0 + (x1 - x0) * tt
    tgt <- all(c(nd[ed$from[j], "node_fam"], nd[ed$to[j], "node_fam"]) %in% op$targets)
    polygon(c(xs, rev(xs)),
            c(ed$s1[j] + (ed$t1[j] - ed$s1[j]) * sm, rev(ed$s0[j] + (ed$t0[j] - ed$s0[j]) * sm)),
            col = if (tgt) adjustcolor("#5B6770", 0.55) else adjustcolor("#9AA3AD", 0.35),
            border = NA)
  }

  ## gene boxes as arrows (direction = majority strand after orientation)
  head <- 0.14
  for (j in seq_len(nrow(nd))) {
    x <- nd$col[j]; y0 <- nd$y0[j]; y1 <- nd$y1[j]
    xl <- x - box_w / 2; xr2 <- x + box_w / 2; ym <- (y0 + y1) / 2
    if (nd$node_fam[j] %in% c("other", "ins_run")) {
      px <- c(xl, xr2, xr2, xl); py <- c(y0, y0, y1, y1)
    } else if (nd$plus[j] >= 0.5) {
      px <- c(xl, xr2 - head, xr2, xr2 - head, xl); py <- c(y0, y0, ym, y1, y1)
    } else {
      px <- c(xr2, xl + head, xl, xl + head, xr2); py <- c(y0, y0, ym, y1, y1)
    }
    is_t <- nd$node_fam[j] %in% op$targets
    polygon(px, py, col = nd$fill[j], border = if (is_t) "black" else "grey35",
            lwd = if (is_t) 1.6 else 0.6, lty = if (nd$node_fam[j] == "ins_run") 3 else 1)
    if (!is.na(nd$mge[j])) {
      mc <- mge_col[[nd$mge[j]]]
      polygon(px, py, col = NA, border = mc, lwd = 1.4, lty = 2)
      tw <- strwidth(mge_tag[[nd$mge[j]]], cex = 0.45) * 1.3
      th <- strheight("I", cex = 0.45) * 1.6
      rect(x - tw / 2, y1 - th / 2, x + tw / 2, y1 + th / 2, col = mc, border = NA)
      text(x, y1, mge_tag[[nd$mge[j]]], cex = 0.45, col = "white", font = 2)
    }
    if (nd$h[j] >= label_min && nd$node_fam[j] != "other") {
      txt <- if (nd$node_fam[j] == "ins_run") paste0("+", nd$run_lab[j], " genes")
             else lab[[nd$node_fam[j]]]
      dark <- fam_sel[j] %in% c("Excess diversity", "Deficit diversity")
      text(x, ym, paste0(txt, "\n", nd$n[j]), cex = 0.55,
           font = if (is_t) 2 else 1, col = if (dark) "white" else "black")
    }
  }

  ## annotation
  title(main = op$name, adj = 0, line = 1.6, cex.main = 1)
  mtext(sprintf("%s  |  junction %s  |  targets present in %d / %d genomes",
                b$msu, b$anchor, b$n_with, N), side = 3, adj = 0, line = 0.4, cex = 0.7)
  axis(1, at = cols, labels = b$col_lab[as.character(cols)], cex.axis = 0.6, tck = -0.02, mgp = c(2, 0.3, 0))
  mtext("Reference gene position relative to first target (unlabelled ticks = insertions)",
        side = 1, line = 1.4, cex = 0.65)
  invisible(nd)
}

draw_legend <- function() {
  plot.new()
  legend("center", horiz = TRUE, bty = "n", cex = 0.75, x.intersp = 0.6,
         legend = c(names(sel_col), "Target gene"),
         fill = c(sel_col, NA), border = c(rep("grey35", length(sel_col)), "black"),
         title = sprintf("Box height = genomes; ribbons = adjacency; BH q < %.2f", q_alpha))
  legend("bottom", horiz = TRUE, bty = "n", cex = 0.75, x.intersp = 0.6,
         legend = paste0(mge_tag, " = ", names(mge_tag)),
         text.col = mge_col[names(mge_tag)], lty = 2, lwd = 1.4, col = mge_col[names(mge_tag)],
         seg.len = 1.5)
}

## ---------------------------------------------------------------- run
## set options(ribbon_functions_only = TRUE) before sourcing to load the
## operon definitions and drawing functions without generating the figures
if (!isTRUE(getOption("ribbon_functions_only"))) {

built <- lapply(operons, build_operon)

ncols <- vapply(built, function(b) length(unique(b$al$col)), numeric(1))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## one panel + legend per page; width scales with the number of gene positions
make_fig <- function(b) {
  layout(matrix(1:2, ncol = 1), heights = c(1, 0.3))
  par(mar = c(2.6, 1, 3, 1))
  plot_ribbon(b)
  par(mar = c(0, 0, 0, 0))
  draw_legend()
}
fig_w <- function(nc) max(9, 0.62 * nc + 1.5)
fig_h <- 4.4

## multi-page PDF (each page sized to its own panel)
pdf(file.path(outdir_fig, "/operon_ribbon_plots.pdf"), width = fig_w(max(ncols)), height = fig_h)
for (b in built) make_fig(b)
dev.off()
## one PNG per operon
for (i in seq_along(built)) {
  png(file.path(outdir_fig, paste0("/operon_ribbon_", gsub("[^A-Za-z0-9]+", "_", built[[i]]$op$name), ".png")),
      width = fig_w(ncols[i]), height = fig_h, units = "in", res = 300)
  make_fig(built[[i]]); dev.off()
}

## per-operon adjacency summaries, for checking against the figure
for (b in built) {
  s <- aggregate(geno ~ col + node_fam, data = b$al, FUN = length)
  names(s)[3] <- "n_genomes"
  m <- aggregate(list(mge_frac = !is.na(b$al$mge)), by = b$al[, c("col", "node_fam")], FUN = mean)
  s <- merge(s, m, by = c("col", "node_fam"))
  s <- s[order(s$col, -s$n_genomes), ]
  write.csv(s, file.path(out_dir, paste0("ribbon_nodes_", gsub("[^A-Za-z0-9]+", "_", b$op$name), ".csv")),
            row.names = FALSE)
}

}  ## end of figure-generating block
