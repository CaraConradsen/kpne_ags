col_steelblue  <- "steelblue4"#"#2EAA5A"
col_blue   <- "dodgerblue2"# "#2ABFBF"
col_skyblue   <- "deepskyblue1"# "blue"
col_lteal <- "#2ECC9A"#"#B85CC8"
col_dteal    <- "#117864"#"#D95030"
col_red <- "#D95030"#"#C8B820"
col_grey   <- "#BBBBBB"
col_dgrey   <- "grey55"
col_white   <- "white"

# ── Block row definitions (list of list(width, colour)) ──────────────────────

# Panel 1 & 4: accessory blocks (with grey gaps)
acc_rows <- list(
  list(list(2,col_grey), list(4,col_steelblue), list(5,col_blue), list(2,col_dgrey), list(5,col_skyblue), list(1,col_grey), list(2,col_red), list(3,col_lteal), list(4,col_dteal), list(1,col_grey)),
  list(list(4,col_steelblue), list(1,col_grey), list(5,col_blue), list(2,col_dgrey), list(5,col_skyblue), list(2,col_grey),  list(2,col_red), list(1,col_grey), list(3,col_lteal), list(4,col_dteal)),
  list(list(4,col_steelblue), list(2,col_dgrey), list(1,col_grey), list(5,col_blue), list(5,col_skyblue), list(2,col_red), list(2,col_grey), list(3,col_lteal), list(1,col_grey), list(4,col_dteal)),
  list(list(2,col_red), list(1,col_grey), list(4,col_steelblue), list(5,col_blue), list(1,col_grey),list(5,col_skyblue), list(1,col_grey), list(4,col_dteal),list(2,col_dgrey), list(3,col_lteal), list(1,col_grey))
)

# Panel 2: core blocks (no grey)
core_rows <- list(
  list(list(3,col_white), list(4,col_steelblue), list(5,col_blue), list(5,col_skyblue), list(2,col_red), list(3,col_lteal), list(4,col_dteal),list(3,col_white)),
  list(list(3,col_white), list(4,col_steelblue), list(5,col_blue), list(5,col_skyblue), list(2,col_red), list(3,col_lteal), list(4,col_dteal), list(3,col_white)),
  list(list(3,col_white), list(4,col_steelblue), list(5,col_blue), list(5,col_skyblue), list(2,col_red), list(3,col_lteal),  list(4,col_dteal), list(3,col_white)),
  list(list(3,col_white), list(2,col_red), list(4,col_steelblue), list(5,col_blue), list(5,col_skyblue), list(4,col_dteal), list(3,col_lteal), list(3,col_white))
)

# ── draw blocks, return x positions for each block ───────────────────
draw_blocks <- function(row, y0, y1, gap = 0.01) {
  total  <- sum(sapply(row, function(b) b[[1]]))
  usable <- 1 - gap * (length(row) - 1)
  x <- 0
  positions <- list()
  for (i in seq_along(row)) {
    b <- row[[i]]
    w <- (b[[1]] / total) * usable
    rect(x, y0, x + w, y1, col = b[[2]], border = NA)
    positions[[i]] <- c(x, x + w)
    x <- x + w + gap
  }
  invisible(positions)
}

# ── draw all 4 rows in a panel, spaced evenly ────────────────────────
draw_panel_rows <- function(rows, gap = 0, 
                            sep_4 = FALSE, gap_4 = 0.1,
                            sep_top = FALSE) {
  # y positions for 4 rows within ylim c(0,1)
  row_h  <- ifelse(sep_top == FALSE, 0.15, 0.145)
  row_gap <- 0.05
  # stack from top
  y_tops <- 1 - 0.05 - seq(0, 3) * (row_h + row_gap)
  # sep out row 4
  if(isTRUE(sep_4)){
    y_tops[4] <- y_tops[4] - gap_4
  }
  
  if(isTRUE(sep_top)){
    y_tops[1] = 0.95
    y_tops[2] = 0.7
    y_tops[3] = 0.45
    y_tops[4] = 0.2
  }
  
  pos_list <- list()
  for (i in seq_along(rows)) {
    pos_list[[i]] <- draw_blocks(rows[[i]], y_tops[i] - row_h, y_tops[i], gap)
  }
  invisible(list(positions = pos_list, y_tops = y_tops, row_h = row_h))
}

# ── draw a ribbon between two rows ───────────────────────────────────
draw_ribbon <- function(x1a, x1b, y_top, x2a, x2b, y_bot, fill, alpha = 0.5) {
  # use bezier-like polygon with control points
  n  <- 20
  t  <- seq(0, 1, length.out = n)
  # top edge: straight
  # sides: linear interpolation
  x_left  <- x1a + (x2a - x1a) * t
  x_right <- x1b + (x2b - x1b) * t
  y_mid   <- y_top + (y_bot - y_top) * t

  xs <- c(x_left, rev(x_right))
  ys <- c(y_mid,  rev(y_mid))

  polygon(xs, ys,
          col    = adjustcolor(fill, alpha.f = alpha),
          border = adjustcolor(fill, alpha.f = 0.8),
          lwd    = 0.8)
}

# ── merge coriented blocks ───────────────────────────────────
merge_blocks <- function(row_dat, pos_info, grad_num = 100,
                         block_grps = list(list(2:4, 6:7), 
                                           list(2:4, 6:7),
                                           list(2:4, 6:7),
                                           list(3:5, 6:7))){
  y_tops   <- pos_info$y_tops
  row_h    <- pos_info$row_h

  for (i in seq_along(row_dat)){
    for (j in 1:length(block_grps[[i]])){
      
      block_cols = sapply(row_dat[[i]][block_grps[[i]][[j]]], `[[`, 2)
      
      colfunc <- colorRampPalette(block_cols)
      
      block_x_pos = unlist(pos_info$positions[[i]][block_grps[[i]][[j]]])
      
      x_seq <- seq(min(block_x_pos), max(block_x_pos), length.out = grad_num)
      
      rect(
        x_seq[-grad_num],        # xleft: all but last
        y_tops[i],
        x_seq[-1],         # xright: all but first
        y_tops[i] - row_h,
        col = colfunc(grad_num-1), # 49 rectangles now
        border = NA
      )
    }
  }
  
}

# ── highlight focal AG ───────────────────────────────────
outline_grey_blocks <- function(rows, positions, y_tops, row_h, 
                                outline_idx, border_col = "grey40", lwd = 2.25) {
  for (i in seq_along(outline_idx)) {
    for (j in outline_idx[[i]]) {
      pos <- positions[[i]][[j]]
      rect(pos[1], y_tops[i] - row_h, pos[2], y_tops[i],
           col = NA, border = border_col, lwd = lwd)
    }
  }
}

# ════════════════════════════════════════════════════════════════════════════
# par(mfrow = c(4,1))

# ── Panel 1: Accessory blocks ─────────────────────────────────────────────
screen(2)               # first sub-screen (screens 1-4 already taken)
par(mar = c(0.75, 2, 0.25, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
draw_panel_rows(acc_rows)
mtext("b", side = 3, adj= -0.1)
mtext("remove accessory blocks (gray)", side = 1, line = -1, cex = 0.7, col = "grey40", font = 2)
close.screen(2)

# ── Panel 2: Core blocks ──────────────────────────────────────────────────
screen(3)
par(mar = c(0.75, 2, 0.25, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
draw_panel_rows(core_rows)
mtext("merge co-oriented core blocks", side = 1, line = -1, cex = 0.7, col = "grey40", font = 2)
close.screen(3)

# ── Panel 3: Rearrangements ───────────────────────────────────────────────
screen(4)
par(mar = c(1, 2, 0, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
info <- draw_panel_rows(core_rows,
                        sep_4 = TRUE, gap_4 = 0.1)

# Extract positions and y coords for ribbon connections
pos_r3   <- info$positions[[3]]   # row 3 block positions
pos_r4   <- info$positions[[4]]   # row 4 block positions
y_tops   <- info$y_tops
row_h    <- info$row_h

y_top <- y_tops[3] - row_h   # bottom of row 3
y_bot <- y_tops[4]            # top of row 4

# merge colours

merge_blocks(core_rows,info,grad_num = 100,
             block_grps = list(list(2:4, 6:7),
                               list(2:4, 6:7),
                               list(2:4, 6:7),
                               list(3:5, 6:7)))


draw_ribbon(pos_r3[[7]][2], pos_r3[[6]][1], y_top,
            pos_r4[[6]][1], pos_r4[[7]][2], y_bot,
            fill = col_lteal)

draw_ribbon(pos_r3[[5]][1], pos_r3[[5]][2], y_top,
            pos_r4[[2]][1], pos_r4[[2]][2], y_bot,
            fill = col_red)

# Arrow
arrows(pos_r4[[8]][1], y_tops[4] - row_h - 0.07,
       pos_r4[[6]][1], y_tops[4] - row_h  - 0.07,
       length = 0.1, col = "grey30", lwd = 1.1)

mtext("account for rearrangements", side = 1, line = -0.25, cex = 0.7, col = "grey40", font = 2)
close.screen(4)

# ── Panel 4: remap AGS ─────────────────────────────────────────────
screen(5)
par(mar = c(0.7, 2, 0.3, 2))
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
info_last <- draw_panel_rows(acc_rows, sep_top = TRUE)

merge_blocks(acc_rows,info_last,grad_num = 100,
             block_grps = list(list(c(2,3,5), 8:9),
                               list(c(1,3,5), 9:10),
                               list(c(1,4,5), c(8,10)),
                               list(c(3,4,6), c(8,10))))

mtext("re-map accessory genes", side = 1, line = -0.5, cex = 0.7, col = "grey40", font = 2)

acc_blocks <- lapply(acc_rows, function(row) {
  which(sapply(row, function(x) x[[2]] == col_grey | x[[2]] == col_dgrey))
})

acc_rows_grey <- lapply(acc_rows, function(row) {
  lapply(row, function(b) {
    if (b[[2]] == col_grey | b[[2]] == col_dgrey) b else list(b[[1]], NA)
  })
})

# add back AGS
draw_panel_rows(acc_rows_grey, sep_top = TRUE)

# attach syntenic loci

# Extract positions and y coords for ribbon connections
pos_r3   <- info_last$positions[[1]]   # row 3 block positions
pos_r4   <- info_last$positions[[2]]   # row 4 block positions
y_tops   <- info_last$y_tops
row_h    <- info_last$row_h

y_top <- y_tops[1] - row_h   # bottom of row 3
y_bot <- y_tops[2]            # top of row 4


draw_ribbon(pos_r3[[4]][1], pos_r3[[4]][2], y_top,
            pos_r4[[4]][1], pos_r4[[4]][2], y_bot,
            fill = col_dgrey)

outline_grey_blocks(
  rows      = acc_rows,
  positions = info_last$positions,
  y_tops    = info_last$y_tops,
  row_h     = info_last$row_h,
  outline_idx = list(
    c(4),   # row 1
    c(4),   # row 2
    c(2),   # row 3
    c(9)    # row 4
  )
)

text(mean(info_last$positions[[1]][[4]]),
     info_last$y_tops[1] + 0.05,
     cex = 0.55,
     "focal accessory gene")

text(mean(info_last$positions[[3]][[2]]),
     info_last$y_tops[3] + 0.05,
     cex = 0.55,
     "upstream non-syntenic")

text(mean(info_last$positions[[4]][[9]]),
     info_last$y_tops[4] + 0.05,
     cex = 0.55,
     "different core block")

# # Redraw row 4 on top of ribbons
# draw_blocks(rear_rows[[4]],
#             y_tops[4] - row_h, y_tops[4])

close.screen(5)

