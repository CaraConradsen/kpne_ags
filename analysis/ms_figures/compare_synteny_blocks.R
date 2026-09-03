## Compare syntenic accessory-gene sections across genomes
## ----------------------------------------------------------------------------
## Each genome contributes a slice between two consecutive core genes (the
## "anchors", e.g. A143 .. A144); the accessory genes (AGs) lie between them.
## Sections may run in reverse and may be translocated, so raw bp coordinates
## are not comparable across genomes. This function orients every genome to a
## reference, then draws one row per genome in one of two layouts:
##
##   layout = "aligned" (default)
##       Genes are placed in ordinal COLUMNS built by a progressive alignment
##       of gene-family order. Homologous genes (same gene_family) line up
##       vertically; a genome missing a gene shows a gap in that column.
##       x is ordinal, not bp -- arrow widths are nominal.
##
##   layout = "coords"
##       True bp retained; each genome is flipped (if reverse) and translated
##       so the start anchor sits at x = 0. Arrow widths are real gene lengths.
##
## Input: data.frame / data.table with columns
##   geno_id, gene_family, start, end, strand, anchor, cols
## 'anchor' carries the bare start/end anchor names on the two core genes
## (e.g. "A143", "A144") and "A143:A144" on the accessory genes between them.
## Orientation of a genome is decided by the strand of its start-anchor gene
## relative to the reference's start-anchor gene.
## ----------------------------------------------------------------------------

plot_gene_sections <- function(genes,
                               start_anchor = NULL,   # auto-detected if NULL
                               end_anchor   = NULL,
                               order     = c("input", "dendro"),
                               ref_geno  = NULL,       # defaults to first genome
                               body_half = 0.32,       # arrow half-height (y-units)
                               label     = TRUE,       # label columns with gene_family
                               label_cex = 0.55,
                               mark_anchors = TRUE) {   # shade the two anchor columns
  
  order  <- match.arg(order)
  genes  <- as.data.frame(genes)
  
  ## --- identify the two anchors -------------------------------------------
  bare <- unique(genes$anchor[!grepl(":", genes$anchor)])
  if (is.null(start_anchor) || is.null(end_anchor)) {
    if (length(bare) != 2L)
      stop("Expected exactly 2 colon-free anchor values, found: ",
           paste(bare, collapse = ", "))
    bare <- bare[order(as.numeric(gsub("\\D", "", bare)))]  # by numeric suffix
    start_anchor <- bare[1]; end_anchor <- bare[2]
  }
  
  ## --- split into per-genome tables (input order), sorted by start ---------
  geno_ids0 <- unique(genes$geno_id)
  gl <- split(genes, factor(genes$geno_id, levels = geno_ids0))
  gl <- lapply(gl, function(d) d[order(d$start), ])
  geno_ids <- geno_ids0
  
  if (is.null(ref_geno)) ref_geno <- geno_ids[1]
  if (!ref_geno %in% geno_ids) stop("ref_geno not found among geno_id values.")
  
  anchor_row <- function(d, a) {
    r <- which(d$anchor == a)
    if (length(r) == 0L)
      stop("Genome ", d$geno_id[1], " lacks anchor '", a, "'.")
    r[1]
  }
  
  ## reference orientation = strand of the reference's start-anchor gene
  ref_strand <- gl[[ref_geno]]$strand[anchor_row(gl[[ref_geno]], start_anchor)]
  
  ## --- orient every genome to the reference --------------------------------
  flip_swap <- function(s) ifelse(s == "+", "-", "+")
  oriented <- lapply(gl, function(d) {
    flipped <- d$strand[anchor_row(d, start_anchor)] != ref_strand
    if (flipped) {                       # mirror the gene order, swap strands
      d <- d[rev(seq_len(nrow(d))), ]
      d$strand <- flip_swap(d$strand)
    }
    d$flipped <- flipped
    d
  })
  
  anchor_fams <- unique(genes$gene_family[genes$anchor %in%
                                            c(start_anchor, end_anchor)])
  
  ## ---- arrow helper (fixed-width, for the aligned layout) -----------------
  draw_arrow <- function(xc, y, fwd, fill, w, hw, bh) {
    lo <- xc - w / 2; hi <- xc + w / 2
    if (isTRUE(fwd)) xs <- c(lo, hi - hw, hi, hi - hw, lo)
    else             xs <- c(hi, lo + hw, lo, lo + hw, hi)
    polygon(xs, y + c(-bh, -bh, 0, bh, bh),
            col = fill, border = "grey20", lwd = 0.5)
  }
  
  
   ## ============================ ALIGNED ==================================
    lcs_merge <- function(columns, s2) {
      n <- length(columns); m <- length(s2)
      if (n == 0L) return(s2[!duplicated(s2)])
      if (m == 0L) return(columns)
      L <- matrix(0L, n + 1L, m + 1L)
      for (i in seq_len(n)) for (j in seq_len(m))
        L[i + 1L, j + 1L] <- if (columns[i] == s2[j]) L[i, j] + 1L
      else max(L[i, j + 1L], L[i + 1L, j])
      i <- n; j <- m; merged <- character(0)
      while (i > 0L || j > 0L) {
        if (i > 0L && j > 0L && columns[i] == s2[j] &&
            L[i + 1L, j + 1L] == L[i, j] + 1L) {
          merged <- c(columns[i], merged); i <- i - 1L; j <- j - 1L
        } else if (j > 0L && (i == 0L || L[i + 1L, j] >= L[i, j + 1L])) {
          merged <- c(s2[j], merged); j <- j - 1L           # genome insertion
        } else {
          merged <- c(columns[i], merged); i <- i - 1L      # gap in genome
        }
      }
      merged[!duplicated(merged)]
    }
    
    seqs <- lapply(oriented, function(d) d$gene_family)
    columns <- seqs[[ref_geno]]
    for (g in setdiff(geno_ids, ref_geno))
      columns <- lcs_merge(columns, seqs[[g]])
    
    # gene label
    gnames <- lapply(oriented, function(d) d$consensus_product)
    gene_names <- gnames[[ref_geno]]
    for (g in setdiff(geno_ids, ref_geno))
      gene_names <- lcs_merge(gene_names, gnames[[g]])
    
    ## genome (row) ordering
    if (order == "dendro" && length(geno_ids) >= 3L) {
      pa <- t(sapply(seqs, function(s) as.integer(columns %in% s)))
      rownames(pa) <- geno_ids
      ord <- order.dendrogram(as.dendrogram(
        hclust(dist(pa, method = "binary"))))   # Jaccard
      geno_ids <- geno_ids[ord]
    }
    
    K <- length(columns)
    n <- length(geno_ids)
    xcol <- setNames(seq_len(K) - 1L, columns)          # start anchor at x = 0
    arrow_w <- 0.8; head_w <- arrow_w * 0.35
    
    par(mar = c(7, 5, 4, 0.5))
    plot(NA, xlim = c(-0.7, K - 0.3), ylim = c(0.5, n + 0.5),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    
    ## shade the two anchor columns
    if (mark_anchors)
      for (af in intersect(anchor_fams, columns)) {
        xc <- xcol[af]
        rect(xc - 0.5, 0.45, xc + 0.5, n + 0.55, col = "grey92", border = NA)
      }
    
    ## one row per genome (first genome at top)
    for (r in seq_len(n)) {
      g <- geno_ids[r]; y <- n - r + 1
      d <- oriented[[g]]
      fam_row <- setNames(seq_len(nrow(d)), d$gene_family)
      segments(-0.5, y, K - 0.5, y, col = "grey90", lwd = 0.5)   # baseline
      for (k in seq_len(K)) {
        fam <- columns[k]; xc <- xcol[fam]
        i <- fam_row[fam]
        if (!is.na(i)) {
          draw_arrow(xc, y, d$strand[i] == "+", d$cols[i],
                     arrow_w, head_w, body_half)
        } else {
          segments(xc - 0.28, y, xc + 0.28, y, col = "grey75", lwd = 0.7)
        }
      }
    }
    
    if (label){
      text(seq_len(K) - 1L, n + 0.65, labels = columns, srt = 90,
           adj = c(0, 0.5), cex = label_cex, xpd = NA)
      text(seq_len(K) - 0.5, 0 - 0.65, labels = gene_names, srt = 90,
           cex = label_cex*0.95, pos= 2, xpd = NA)
      
    }
    
    axis(2, at = n:1, labels = geno_ids, las = 1, tick = FALSE,
         cex.axis = max(0.4, min(0.85, 16 / n)))
    mtext(sprintf("Accessory-gene section  %s -> %s  (aligned columns)",
                  start_anchor, end_anchor), side = 3, line = 2, cex = 0.95)

    
    return(invisible(list(anchors = c(start = start_anchor, end = end_anchor),
                          columns = columns, order = geno_ids,
                          flipped = vapply(oriented[geno_ids],
                                           function(d) d$flipped[1], logical(1)))))
  
}

## ---- example ----------------------------------------------------------------
# library(data.table)
# genes <- fread("sections.csv")              # geno_id, gene_family, start,
#                                             # end, strand, anchor, cols
# n_geno <- length(unique(genes$geno_id))
# n_fam  <- length(unique(genes$gene_family))
# pdf("gene_sections.pdf",
#     width  = max(7, 0.32 * n_fam  + 3),
#     height = max(4, 0.32 * n_geno + 2))
# info <- plot_gene_sections(genes, order = "dendro")   # aligned by default
# dev.off()
# info$flipped        # which genomes were reverse-oriented
# info$columns        # the aligned gene-family column order
