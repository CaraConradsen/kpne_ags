# Do our outlier AGs have more diversity then the strains that contain them?

set.seed(42)

# Import data -------------------------------------------------------------

piS_outlier_genes <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]$gene_family

# get piS values
ag_piS_dt <- fread(paste0(outdir_dat, "/ag_seg_sites_dt.csv"))

ag_piS_dt <- ag_piS_dt[gene_family %chin% piS_outlier_genes]

ag_piS_dt <- ag_piS_dt[,.(
  ag_pi_s = mean(pi_s, na.rm = TRUE),
  freq = mean(number_genomes)
), gene_family]

# set up genome sets
geno_anno <- unique(fread(paste0(outdir_dat, "/pangenome_anno.csv"), 
                         select = c("gene_family", "geno_id")))[gene_family %chin% piS_outlier_genes]

AG_genomes <- split(geno_anno$geno_id, geno_anno$gene_family)


# Get core genes  ---------------------------------------------------------

core_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                   select = c("geno_id", "gene_family", "locus_tag", "fus_locus_tag", 
                              "number_genomes", "start","end", "average_dose",
                              "strand","ST","ag_type"))


# remove get only core genes and remove paralogs
core_anno <- core_anno[ag_type=="core"][average_dose <= 1]

list_unique_cores = unique(core_anno$gene_family)


# Prefix tables -----------------------------------------------------------

fourfold_prefixes           <- c("CT","GT","TC","CC","AC","GC","CG","GG")
twofold_purine_prefixes     <- c("AA","AG","CA","GA","TT")              # A/G synonymous
twofold_pyrimidine_prefixes <- c("AA","AG","CA","GA","TT","TA","TG")    # C/T synonymous

# codon-walker ------------------------------------------------------------

codon_walker <- function(aln, mode = c("4fold","2fold")) {
  mode  <- match.arg(mode)
  L     <- ncol(aln)
  n_cod <- L / 3L
  out   <- vector("list", n_cod)
  
  for (i in seq_len(n_cod)) {
    sub  <- aln[, (3*i-2):(3*i)]
    keep <- rowSums(sub == "N" | sub == "-") == 0
    sub  <- sub[keep, , drop = FALSE]
    if (nrow(sub) < 2) next
    
    prefix <- paste0(sub[,1], sub[,2])
    if (length(unique(prefix)) > 1) next       # mixed 1st/2nd position -> skip
    p     <- prefix[1]
    third <- sub[, 3]
    
    if (mode == "4fold") {
      if (!(p %in% fourfold_prefixes)) next
    } else {                                   # 2fold
      bases  <- unique(third)
      is_pur <- all(bases %in% c("A","G")) && (p %in% twofold_purine_prefixes)
      is_pyr <- all(bases %in% c("C","T")) && (p %in% twofold_pyrimidine_prefixes)
      if (!is_pur && !is_pyr) next
    }
    
    counts   <- tabulate(match(third, c("A","C","G","T")), 4)
    out[[i]] <- c(n = nrow(sub), counts)
  }
  do.call(rbind, out)
}

# gene walker -- reads each gene ONCE, slices across all frequencies -------

walk_ag_genomes <- function(fasta_dir, unique_cores, AG_genomes,
                                mode = "4fold",
                                n_cores = parallel::detectCores() - 4){
  
  ag_gene <- names(AG_genomes)
  
  per_gene <- function(focal_gene){
    
    temp_string <- tryCatch(
      readDNAStringSet(paste0(fasta_dir, focal_gene, ".fasta")),
      error = function(e) return(NULL)
    )
    if (is.null(temp_string)) return(NULL)   # per-gene skip, not whole run
    
    gene_length <- unique(width(temp_string))
    
    # convert once
    gene_mat <- as.matrix(temp_string)
    gene_mat[!gene_mat %in% c("A","C","G","T","-")] <- "N"
    
    # strip trailing _<digits> suffix so headers match clade members exactly
    core_names <- sub("_[0-9]+$", "", rownames(gene_mat))
    
    # loop frequencies, then clades within each frequency
    gene_out <- list()
    for (k in ag_gene){
      set <- AG_genomes[[k]]                          # the full vector of carrier genomes
      sub_mat <- gene_mat[core_names %in% set, , drop = FALSE]
      if (nrow(sub_mat) < 2) next
      cc <- codon_walker(sub_mat, mode = mode)
      if (is.null(cc)) next
      gene_out[[k]] <- cbind(cc,
                             gene_length   = gene_length,
                             core_gene     = focal_gene,
                             gene_family = k)
    }
    
    do.call(rbind, gene_out)
  }
  
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(Biostrings))
    parallel::clusterExport(cl,
                            c("codon_walker", "fourfold_prefixes",
                              "twofold_purine_prefixes", "twofold_pyrimidine_prefixes",
                              "AG_genomes", "ag_gene", "mode", "fasta_dir"),
                            envir = environment())
    results <- parallel::parLapply(cl, unique_cores, per_gene)
  } else {
    results <- parallel::mclapply(unique_cores, per_gene, mc.cores = n_cores)
  }
  
  do.call(rbind, results)
}

# Run ---------------------------------------------------------------------

fasta_dir <- paste0(outdir_dat, "/core_fus_alns/")   # core fusion alignments must exist

start.time <- Sys.time()
sites_4 <- walk_ag_genomes(fasta_dir, list_unique_cores, AG_genomes,
                               mode = "4fold")
sites_2 <- walk_ag_genomes(fasta_dir, list_unique_cores, AG_genomes,
                               mode = "2fold")
end.time <- Sys.time()
end.time - start.time # Time difference of 41.7497 mins

sites_4_dt <- as.data.table(sites_4); setnames(sites_4_dt, c("n","A","C","G","T", "gene_length", "core_gene", "gene_family"))
cols = c("n","A","C","G","T", "gene_length")
sites_4_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_4_dt[, class := "4fold"]

sites_2_dt <- as.data.table(sites_2); setnames(sites_2_dt, c("n","A","C","G","T", "gene_length", "core_gene", "gene_family"))
cols = c("n","A","C","G","T", "gene_length")
sites_2_dt[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
sites_2_dt[, class := "2fold"]

sites_combined_dt <- rbind(sites_4_dt, sites_2_dt)



fwrite(sites_combined_dt, paste0(outdir_dat,"/ag_genos_core_theta.csv"))
# sites_combined_dt <- fread(paste0(outdir_dat,"/ag_genos_core_theta.csv"))


# Get piS ------------------------------------------------

per_ag <- sites_combined_dt[, {
  Ls <- .N
  n_j      <- n
  counts   <- cbind(A, C, G, T)
  
  p_mat    <- counts / n_j
  h        <- (n_j / (n_j - 1)) * (1 - rowSums(p_mat^2))

  .(pi_s      = as.numeric(mean(h)))
  
}, by = .(core_gene, gene_family, class)]


combined_pi_dt <- per_ag[, .(
  pi_s     = mean(pi_s)
), by = gene_family]


core_theta_per_ag_dt <- combined_pi_dt[,.(
  core_pi_s = mean(pi_s)),
  by = "gene_family"]


all_pi_s_dt <- merge(ag_piS_dt, core_theta_per_ag_dt,
                     all.x = TRUE, by = "gene_family")

fwrite(all_pi_s_dt, paste0(outdir_dat,"/ag_vs_local_core_piS.csv"))
# all_pi_s_dt <- fread(paste0(outdir_dat,"/ag_vs_local_core_piS.csv"))


# add mobilome cogs ----------------------------------------------------------------

cogs_outliers <- fread(paste0(outdir_dat, "/bh_p_cog_set.csv"))

cogs_outliers <- cogs_outliers[grepl("X", COG_letter)]

# plot --------------------------------------------------------------------

png(paste0(outdir_fig,"/ag_vs_local_core_piS_with_Xcog.png"),
    width = 16, height = 16, units = "cm", res = 300,
    pointsize = 11, type = "cairo")

with(all_pi_s_dt,
     plot(core_pi_s, ag_pi_s,
          pch = 16, col = "white",
          bty = "L",
          # xlim = c(0,0.12),
          # ylim = c(0,0.12),
          main = "Synonymous divergence of outlier AGs versus\n core-genome divergence of their carrier strains",
          xlab = expression("Core"~pi[S]),
          ylab = expression("Accessory"~pi[S])))

xmax <- max(all_pi_s_dt$core_pi_s)

# Boundaries in SYNONYMOUS units. Do NOT paste the all-sites literature
# values here directly. Replace each with a value computed in your own
# synonymous scale (see note below).
b_within  <- xmax   # within-phylogroup  — VERIFY in synonymous units
# b_species <- 0.04    # between-KpSC-species — PLACEHOLDER, all-sites value; convert
b_ani     <- 0.07    # ANI species boundary — PLACEHOLDER; convert

band_col <- function(col, alpha) adjustcolor(col, alpha.f = alpha)

# band 1: within-phylogroup up to between-species
polygon(c(0, xmax, xmax, 0),
        c(0, 0 ,b_within, b_within),
        col = band_col("grey70", 0.15), border = NA)

# # band 2: between-species up to ANI boundary
# polygon(c(0, xmax, xmax, 0),
#         c(b_species, b_species, b_ani, b_ani),
#         col = band_col("orange", 0.15), border = NA)

# band 3: above ANI boundary (open-topped — use plot's upper limit)
ytop <- par("usr")[4]
polygon(c(0, xmax, xmax, 0),
        c(b_within, b_within, ytop, ytop),
        col = band_col("orange", 0.15), border = NA)#band_col("red", 0.15)

# boundary lines + labels
abline(h = c(b_within), lty = 2, col = "grey40")#, b_species, b_ani
text(0, b_within+b_within/10,  expression("max core"~pi[S]~"(most divergent carrier set)"),
     pos = 4, cex = 0.8, adj = 1, xpd = TRUE)
# text(xmax, b_species, "between-KpSC species",   pos = 3, cex = 0.7, adj = 1, xpd = TRUE)
# text(xmax, b_ani,     "ANI species boundary",   pos = 3, cex = 0.7, adj = 1, xpd = TRUE)

with(all_pi_s_dt,
     points(core_pi_s, ag_pi_s,
          pch = ifelse(gene_family %in% cogs_outliers$gene_family, 17, 16), 
          col = ifelse(gene_family %in% cogs_outliers$gene_family, rgb(0.1,0.1, 0.1, alpha = 0.95), rgb(0.2,0.2, 0.8, alpha = 0.5)), 
          bty = "L",
          # xlim = c(0,0.12),
          # ylim = c(0,0.12),
          main = "Synonymous divergence of outlier AGs versus\n core-genome divergence of their carrier strains",
          xlab = expression("Core"~pi[S]),
          ylab = expression("Accessory"~pi[S])))

legend(
  0, 0.13,
  legend = c("Other", "Mobilome (X)"),
  pch = c(16, 17),
  horiz = TRUE,
  col = c(
    rgb(0.2, 0.2, 0.8, alpha = 0.5),
    rgb(0.1, 0.1, 0.1, alpha = 0.95)
  ),
  bty = "n"
)

abline(a=0, b=1, lwd = 2, col = "red")
text(0.003,0.003, "1:1",col = "red", pos = 3, cex = 0.8,font =2, adj = 1, xpd = TRUE)

dev.off()


# Cluster -----------------------------------------------------------------
pan_anno <- unique(fread("./output/data/all_pirate_anno_full.csv",
                         select = c("gene_family", "consensus_product", "ice", "is_fam")))


pan_anno[gene_family %chin%  all_pi_s_dt[core_pi_s<0.005 & ag_pi_s<0.015 & freq==32, gene_family]]



# by frequency ------------------------------------------------------------
all_pi_s_dt[, mag_scale := ag_pi_s/core_pi_s]

with(all_pi_s_dt,
     boxplot(ag_pi_s/core_pi_s~freq, 
          main = "Synonymous divergence of outlier AGs scaled by \n core-genome divergence of their carrier strains, plotted for frequency",
          ylab = expression("Accessory"~pi[S]~"/"~"Core"~pi[S]),
          pch = 16, ylim = c(0,20),
          xlab = "Pangenome frequency"))


all_pi_s_dt[core_pi_s<0.005 & ag_pi_s<0.015 & freq==32]
