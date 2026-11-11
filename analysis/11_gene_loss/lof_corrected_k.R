# Correct k and S for pseudogenised copies, and re-test against the null
# pangenome_full_annotation.R drops fission/fusion loci, so these copies are
# missing from both k and S. Here each loss-of-function copy (from 
# fis_fus_loss_of_function.R) is added back: k' = k + copies, and S' is 
# recounted with the copies placed into the AG alignment by homology 
# (via the PIRATE representative, which is a row of that alignment).
# Synonymous sites are defined from the intact copies only, as in 
# obs_ag_vals.R, so L is unchanged and the existing null tables apply.

# Parameters --------------------------------------------------------------

add_classes  <- c("nonsense", "frameshift", "frameshift_homopolymer",
                  "premature_stop_other")
drop_multi   <- TRUE   # loci joining differently named genes (e.g. benB,benC)
                       # are likely the ancestral split state, not a pseudogene
fdr_cut      <- 0.05
max_k_null   <- 259    # largest k in the neutral tables (prod_ag_Ls_2.89)
k_breaks     <- c(1, 5, 10, 20, 50, 100, 150, 200, 240, 259)

fourfold_prefixes           <- c("CT","GT","TC","CC","AC","GC","CG","GG")
twofold_purine_prefixes     <- c("AA","AG","CA","GA","TT")
twofold_pyrimidine_prefixes <- c("AA","AG","CA","GA","TT","TA","TG")

# Import data -------------------------------------------------------------

gene_pvals <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"),
                    select = c("gene_family", "freq", "seg_sites", "syn_sites", 
                               "expected_S", "p_two"))

lof_calls <- fread(paste0(outdir_dat, "/fis_fus_lof_calls.csv"))

lof_add <- lof_calls[lof_class %chin% add_classes & rep_is_fis_fus == FALSE &
                       gene_family %chin% gene_pvals$gene_family]

if (drop_multi) lof_add <- lof_add[multi_gene == FALSE]

rep_headers <- names(readDNAStringSet("./input_data/PIRATE_260_hybrid_chr_hps0.6_out/representative_sequences.ffn"))
rep_info <- data.table(gene_family = sub(";.*", "", rep_headers),
                       rep_locus   = sub(".*locus_tag=([^;]+).*", "\\1", rep_headers))

# Functions ---------------------------------------------------------------

# synonymous sites and segregating sites, before and after adding copies.
# Site class is set from intact rows; an added copy contributes at a site
# only if it is ungapped there, shares the 1st/2nd positions and keeps the
# 2-fold base class, otherwise it is treated as missing at that site.
syn_counts <- function(mat, n_intact) {
  
  intact <- seq_len(n_intact)
  added  <- setdiff(seq_len(nrow(mat)), intact)
  n_cod  <- ncol(mat) %/% 3
  
  out <- vector("list", n_cod)
  
  for (i in seq_len(n_cod)) {
    
    sub <- mat[, (3*i - 2):(3*i), drop = FALSE]
    ok  <- rowSums(sub == "N" | sub == "-") == 0
    
    si <- sub[intact[ok[intact]], , drop = FALSE]
    if (nrow(si) < 2) next
    
    pre <- paste0(si[, 1], si[, 2])
    if (length(unique(pre)) > 1) next
    p <- pre[1]
    b <- si[, 3]
    
    allowed <- if (p %in% fourfold_prefixes) {
      c("A", "C", "G", "T")
    } else if (all(b %in% c("A", "G")) && p %in% twofold_purine_prefixes) {
      c("A", "G")
    } else if (all(b %in% c("C", "T")) && p %in% twofold_pyrimidine_prefixes) {
      c("C", "T")
    } else {
      NULL
    }
    if (is.null(allowed)) next
    
    b_add <- character(0)
    if (length(added)) {
      sa <- sub[added[ok[added]], , drop = FALSE]
      keep <- paste0(sa[, 1], sa[, 2]) == p & sa[, 3] %in% allowed
      b_add <- sa[keep, 3]
    }
    
    out[[i]] <- c(seg_before = length(unique(b)) >= 2,
                  seg_after  = length(unique(c(b, b_add))) >= 2)
  }
  
  out <- do.call(rbind, out)
  c(syn_sites = nrow(out), 
    S_before  = sum(out[, "seg_before"]), 
    S_after   = sum(out[, "seg_after"]))
}

# Recount S with pseudogenised copies -----------------------------------------

focal_fams <- unique(lof_add$gene_family)

recount <- rbindlist(lapply(focal_fams, function(f) {
  
  aln <- readDNAStringSet(paste0(outdir_dat, "/ag_fus_alns/", f, ".fasta"))
  mat <- as.matrix(aln)
  mat[!mat %in% c("A", "C", "G", "T", "-")] <- "N"
  
  rep_row <- match(rep_info[gene_family == f, rep_locus], rownames(mat))
  if (is.na(rep_row)) return(data.table(gene_family = f, note = "rep not in alignment"))
  
  # one added copy per genome, not already a carrier
  aln_genos <- sub("_[0-9]+$", "", rownames(mat))
  f_add <- lof_add[gene_family == f & !geno_id %chin% aln_genos]
  f_add <- f_add[!duplicated(geno_id)]
  if (nrow(f_add) == 0) return(data.table(gene_family = f, note = "no new genomes"))
  
  # reference position of each alignment column
  rep_chr <- mat[rep_row, ]
  ref_pos <- cumsum(rep_chr != "-")
  
  add_mat <- t(vapply(f_add$loc_on_ref, function(s) {
    s_chr <- strsplit(s, "")[[1]]
    row <- rep("-", ncol(mat))
    row[rep_chr != "-"] <- s_chr[ref_pos[rep_chr != "-"]]
    row
  }, character(ncol(mat)), USE.NAMES = FALSE))
  add_mat[!add_mat %in% c("A", "C", "G", "T", "-")] <- "N"
  
  sc <- syn_counts(rbind(mat, add_mat), n_intact = nrow(mat))
  
  data.table(gene_family = f, n_added = nrow(f_add),
             n_nonsense   = f_add[lof_class == "nonsense", .N],
             n_frameshift = f_add[lof_class %chin% c("frameshift", "frameshift_homopolymer"), .N],
             syn_sites_check = sc[["syn_sites"]], 
             S_before_check  = sc[["S_before"]],
             seg_sites_new   = sc[["S_after"]], note = "")
}), fill = TRUE)

recount[note != "", .N, by = note]

gene_pvals <- merge(gene_pvals, recount[note == ""], all.x = TRUE, by = "gene_family")

# the recount must reproduce obs_ag_vals.R on the intact copies
gene_pvals[!is.na(n_added), .(L_match = mean(syn_sites_check == syn_sites),
                              S_match = mean(S_before_check == seg_sites))]

gene_pvals[is.na(n_added), `:=`(n_added = 0L, seg_sites_new = seg_sites)]
gene_pvals[, k_new := pmin(freq + n_added, max_k_null)]

# Re-test at k' ---------------------------------------------------------------

need <- unique(gene_pvals[n_added > 0L, .(syn_sites, k_new)])

theta_prod_new <- rbindlist(lapply(unique(need$syn_sites), function(L) {
  d <- fread(paste0(outdir_dat, "/prod_ag_Ls_2.89/pi_prod_lng_m", L, ".csv"))
  cbind(d[freq %in% need[syn_sites == L, k_new]], syn_sites = L)
}))

expected_new <- theta_prod_new[, .(expected_S_new = sum(seg_sites * p)), 
                               by = .(freq, syn_sites)]

gene_pvals[expected_new, on = .(k_new = freq, syn_sites), 
           expected_S_new := i.expected_S_new]

gene_pvals[theta_prod_new, on = .(k_new = freq, syn_sites, seg_sites_new = seg_sites),
           `:=`(p_new = i.p, cdf_new = i.cdf)]

# p values as in fdr_correction.R
gene_pvals[, p_lower_new := cdf_new]
gene_pvals[, p_upper_new := 1 - cdf_new + p_new]
gene_pvals[n_added > 0L & is.na(p_new) & seg_sites_new > expected_S_new, 
           `:=`(p_lower_new = 1, p_upper_new = 1e-10)]
gene_pvals[, p_two_new := pmin(1, 2 * pmin(p_lower_new, p_upper_new))]

gene_pvals[n_added == 0L, `:=`(p_two_new = p_two, expected_S_new = expected_S)]

gene_pvals[is.na(p_two_new)]

# BH over the same gene set before and after
gene_pvals[, p_bh     := p.adjust(p_two, method = "BH")]
gene_pvals[, p_bh_new := p.adjust(p_two_new, method = "BH")]

gene_pvals[, call := fcase(p_bh < fdr_cut & seg_sites > expected_S, "high",
                           p_bh < fdr_cut & seg_sites < expected_S, "low",
                           default = "ns")]

gene_pvals[, call_new := fcase(p_bh_new < fdr_cut & seg_sites_new > expected_S_new, "high",
                               p_bh_new < fdr_cut & seg_sites_new < expected_S_new, "low",
                               default = "ns")]

gene_pvals[, transition := paste(call, "->", call_new)]

gene_pvals[, k_band := cut(freq, breaks = k_breaks, include.lowest = TRUE,
                           labels = paste0(head(k_breaks, -1) + 1, "-", k_breaks[-1]))]

# summary of affected AGs
gene_pvals[n_added > 0L, .N, by = transition][order(-N)]

# Plot -----------------------------------------------------------------

call_cols <- c(high = "#1E6FB8", ns = "grey45", low = "#D9531E")

trans_lv <- c("high -> high", "high -> ns", "high -> low",
              "ns -> high", "ns -> low",
              "low -> high", "low -> ns", "low -> low")

trans_cols <- c("#1E6FB8", "#9CC9EF", "#7A4FA3",
                "#4C9A5E", "#F6B48A",
                "#4B2E6B", "#F9D5BD", "#D9531E")

aff <- gene_pvals[n_added > 0L]

png(paste0(outdir_fig, "/lof_corrected_k.png"),
    width = 12, height = 5.5, units = "in", res = 600, pointsize = 12)

par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(5, 5, 2, 9), bty = "o", xpd = FALSE)

## a: S per site vs k, AGs with pseudogenised copies moved to (k', S') ------

plot(gene_pvals$freq, gene_pvals$seg_sites / gene_pvals$syn_sites, 
     type = "n", las = 1, xlim = c(0, 260),
     ylim = c(0, max(c(gene_pvals$seg_sites, gene_pvals$seg_sites_new) / gene_pvals$syn_sites)),
     xlab = expression("Pangenome frequency (" * italic(k) * ")"),
     ylab = "Segregating sites per site")

points(gene_pvals[n_added == 0L, freq], gene_pvals[n_added == 0L, seg_sites / syn_sites],
       pch = 16, cex = 0.5, col = adjustcolor("grey75", 0.4))

segments(x0 = aff$freq, x1 = aff$k_new,
         y0 = aff$seg_sites / aff$syn_sites, y1 = aff$seg_sites_new / aff$syn_sites,
         col = adjustcolor("grey30", 0.5), lwd = 0.6)

points(aff$freq, aff$seg_sites / aff$syn_sites, pch = 1, cex = 0.6, 
       col = adjustcolor(call_cols[aff$call], 0.8))

points(aff$k_new, aff$seg_sites_new / aff$syn_sites, pch = 16, cex = 0.7, 
       col = adjustcolor(call_cols[aff$call_new], 0.8))

par(xpd = NA)
legend(x = 270, y = par("usr")[4], bty = "n", cex = 0.75,
       legend = c("No LoF copies", "Before (k, S)", "After (k', S')",
                  "High outlier", "ns", "Low outlier"),
       pch = c(16, 1, 16, 15, 15, 15),
       col = c("grey75", "grey30", "grey30", call_cols))
par(xpd = FALSE)

mtext("a", side = 3, adj = -0.15, line = 0.5, font = 2, cex = 1.2)

## b: call transitions of affected AGs by k band -------------------------
# ns -> ns not plotted

trans_tab <- dcast(aff, transition ~ k_band, value.var = "gene_family", 
                   fun.aggregate = length, drop = FALSE)

keep    <- trans_lv %in% trans_tab$transition
bar_mat <- as.matrix(trans_tab[match(trans_lv[keep], transition), -1])

n_aff  <- aff[, .N, by = k_band][match(colnames(bar_mat), k_band), N]
n_aff[is.na(n_aff)] <- 0L
n_band <- gene_pvals[, .N, by = k_band][match(colnames(bar_mat), k_band), N]

bp <- barplot(bar_mat, col = trans_cols[keep], border = NA, las = 1,
              ylab = "AGs with pseudogenised copies", xlab = "", 
              names.arg = rep(NA, ncol(bar_mat)),
              ylim = c(0, max(colSums(bar_mat), 1) * 1.15))

# angled band labels (horizontal ones overlap)
par(xpd = NA)
text(bp, par("usr")[3] - 0.03 * diff(par("usr")[3:4]), colnames(bar_mat),
     srt = 45, adj = c(1, 1), cex = 0.7)
par(xpd = FALSE)
mtext(expression("Pangenome frequency (" * italic(k) * ")"), side = 1, line = 3.5)

# affected / total AGs in band
text(bp, colSums(bar_mat), paste0(n_aff, "/", n_band), pos = 3, cex = 0.65)

par(xpd = NA)
legend(x = max(bp) + 1, y = par("usr")[4], bty = "n", cex = 0.75,
       legend = trans_lv[keep], fill = trans_cols[keep], border = NA,
       title = "Call (k -> k')")
par(xpd = FALSE)

mtext("b", side = 3, adj = -0.15, line = 0.5, font = 2, cex = 1.2)

dev.off()

# Output ------------------------------------------------------------------

fwrite(gene_pvals, paste0(outdir_dat, "/lof_corrected_k.csv"))
