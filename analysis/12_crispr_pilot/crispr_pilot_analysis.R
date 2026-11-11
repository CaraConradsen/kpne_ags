
# CRISPR-Cas pilot: are the two type I-E loci doing anything detectable? ---
#
# Three questions, each with a pre-specified test:
#   1. Are the two loci mutually exclusive, or is zero overlap what chance
#      gives when 30 and 15 genomes are drawn from 260?
#   2. Do they target different mobile elements? If they partition phage
#      space, that is a mechanism for exclusion and for balancing selection.
#   3. Do carriers carry fewer prophages? The functional test.
#
# Inputs come from extract_spacers.py and the blastn run in out/.

q_cut <- 0.05
pal   <- c(A = "#3FB0E8", B = "#B22222", none = "grey75")

dir_out <- "analysis/12_crispr_pilot/out"

arr  <- fread(file.path(dir_out, "array_summary.csv"))
load <- fread(file.path(dir_out, "prophage_load.csv"))

# 1. Mutual exclusivity ---------------------------------------------------
# 30 genomes carry locus A and 15 carry locus B out of 260. Under
# independence the expected overlap is 30 * 15 / 260 = 1.7 genomes, so an
# observed overlap of zero is not by itself evidence of exclusion.

nA <- load[crispr == "A", .N]; nB <- load[crispr == "B", .N]; N <- nrow(load)
tab <- matrix(c(0L, nB, nA, N - nA - nB), 2, 2, byrow = TRUE,
              dimnames = list(c("B+", "B-"), c("A+", "A-")))
ex  <- fisher.test(tab)

cat("\n1. MUTUAL EXCLUSIVITY\n")
print(tab)
cat(sprintf("   expected overlap under independence: %.2f genomes\n", nA * nB / N))
cat(sprintf("   Fisher's exact P = %.3f\n", ex$p.value))

# 2. Spacer targets -------------------------------------------------------
# A protospacer is called at >= 90%% identity over >= 90%% of the spacer.
# Targets are the 944 prophage regions and 731 plasmid contigs of this same
# 260-genome set, so "targeting" means targeting elements circulating in
# this population rather than a reference database.

read_hits <- function(L) {
  f <- file.path(dir_out, paste0("hits_", L, ".tsv"))
  h <- fread(f, col.names = c("spacer", "target", "pident", "alen",
                              "mm", "qlen", "evalue", "bits"))
  h <- h[pident >= 90 & alen >= 0.9 * qlen]
  h[, `:=`(locus      = L,
           src_genome = tstrsplit(spacer, "\\|", keep = 2)[[1]],
           tgt_genome = tstrsplit(target, "\\|", keep = 1)[[1]],
           tgt_kind   = fifelse(grepl("\\|reg", target), "prophage", "plasmid"))]
  h[, phage := fifelse(tgt_kind == "prophage",
                       tstrsplit(target, "\\|", keep = 4)[[1]], NA_character_)]
  h[]
}

hits <- rbindlist(lapply(c("A", "B"), read_hits))

n_sp <- arr[, .(n_spacers = sum(n_spacers)), by = locus]
hit_sp <- hits[, .(hit = uniqueN(spacer)), by = locus]

cat("\n2. SPACER TARGETS\n")
print(merge(n_sp, hit_sp, by = "locus")[, .(locus, n_spacers, hit,
                                            pct = round(100 * hit / n_spacers))])
print(hits[, .N, by = .(locus, tgt_kind)][order(locus, tgt_kind)])

# do the two loci target the same phages?
ph <- hits[!is.na(phage), .(locus, phage)]
A_set <- ph[locus == "A", unique(phage)]
B_set <- ph[locus == "B", unique(phage)]
jac   <- length(intersect(A_set, B_set)) / length(union(A_set, B_set))
cat(sprintf("   phage types: A = %d, B = %d, shared = %d, Jaccard = %.2f\n",
            length(A_set), length(B_set), length(intersect(A_set, B_set)), jac))

# cross-targeting: does each locus target elements in the other's carriers?
cA <- arr[locus == "A" & n_spacers > 0, geno_id]
cB <- arr[locus == "B" & n_spacers > 0, geno_id]
cat(sprintf("   A spacers hit elements in %d/%d B carriers; B spacers in %d/%d A carriers\n",
            uniqueN(hits[locus == "A" & tgt_genome %chin% cB, tgt_genome]), length(cB),
            uniqueN(hits[locus == "B" & tgt_genome %chin% cA, tgt_genome]), length(cA)))

# self-targeting is selected against, so it should be rare; a check on the
# protospacer criterion rather than a result
cat(sprintf("   self-targeting hits: %d of %d (%.1f%%)\n",
            hits[src_genome == tgt_genome, .N], nrow(hits),
            100 * hits[src_genome == tgt_genome, .N] / nrow(hits)))

# 3. Prophage load --------------------------------------------------------

cat("\n3. PROPHAGE LOAD\n")
print(load[, .(n = .N, median_intact = as.numeric(median(n_intact)),
               mean_intact = round(mean(n_intact), 2)), by = crispr][order(crispr)])
load[, grp := fifelse(crispr == "none", "none", "crispr")]
cat(sprintf("   any CRISPR vs none : Wilcoxon P = %.3f\n",
            wilcox.test(n_intact ~ grp, data = load)$p.value))
cat(sprintf("   locus B vs none    : Wilcoxon P = %.3f\n",
            wilcox.test(load[crispr == "B", n_intact],
                        load[crispr == "none", n_intact])$p.value))

fwrite(hits, file.path(dir_out, "spacer_hits_filtered.csv"))

# Figure ------------------------------------------------------------------

png(file.path(outdir_fig, "/crispr_pilot.png"), width = 26, height = 7.5,
    units = "cm", res = 300, pointsize = 10, type = "cairo")
par(mfrow = c(1, 4), mar = c(4.2, 4.4, 3, 1), bty = "L")

# a. array size
bx <- split(arr[n_spacers > 0, n_spacers], arr[n_spacers > 0, locus])
boxplot(bx, col = pal[names(bx)], border = "grey30", ylab = "Spacers per genome",
        xlab = "Locus", outline = FALSE)
stripchart(bx, vertical = TRUE, method = "jitter", jitter = 0.12, pch = 16,
           cex = 0.6, col = "grey25", add = TRUE)
mtext("a  Array size", side = 3, line = 0.8, adj = 0, font = 2)

# b. what they hit
tk <- dcast(hits[, .N, by = .(locus, tgt_kind)], tgt_kind ~ locus, value.var = "N")
m  <- as.matrix(tk[, -1]); rownames(m) <- tk$tgt_kind
barplot(prop.table(m, 2) * 100, col = c("#00285B", "grey70"), border = NA,
        ylab = "Protospacer hits (%)", xlab = "Locus")
legend("topright", bty = "n", fill = c("#00285B", "grey70"),
       legend = rownames(m), cex = 0.85)
mtext("b  Target type", side = 3, line = 0.8, adj = 0, font = 2)

# c. phage overlap
sh <- length(intersect(A_set, B_set))
barplot(c(`A only` = length(setdiff(A_set, B_set)), Shared = sh,
          `B only` = length(setdiff(B_set, A_set))),
        col = c(pal[["A"]], "grey60", pal[["B"]]), border = NA,
        ylab = "Phage types targeted")
mtext(sprintf("c  Overlap (Jaccard %.2f)", jac), side = 3, line = 0.8, adj = 0, font = 2)

# d. prophage load
bl <- split(load$n_intact, factor(load$crispr, levels = c("none", "A", "B")))
boxplot(bl, col = pal[names(bl)], border = "grey30", outline = FALSE,
        ylab = "Intact prophages", xlab = "CRISPR status")
stripchart(bl, vertical = TRUE, method = "jitter", jitter = 0.12, pch = 16,
           cex = 0.4, col = "grey35", add = TRUE)
mtext("d  Prophage load", side = 3, line = 0.8, adj = 0, font = 2)

dev.off()
cat("\nwrote", file.path(outdir_fig, "crispr_pilot.png"), "\n")
