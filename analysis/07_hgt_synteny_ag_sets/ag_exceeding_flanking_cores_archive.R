# do wee see signs of HGT when we contrast an AG to it's cores?

# get gene and core block info
full_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                  select = c("gene_family","consensus_product", # "number_genomes", 
                             "anchor", "msu")))[msu!=""]

# get core pis
core_genes <- full_anno[!grepl(":", anchor)]

# get core pis
core_piS <- fread(paste0(outdir_dat, "/core_age_S_dt.csv"),
                       select = c("gene_family", "pi_s"))

# add to main dt
core_genes <- merge(core_genes, core_piS,
                    all.x = TRUE,
                    by = "gene_family")

# Now for the AGs
ag_genes <- full_anno[grepl(":", anchor)]

# get pis
ag_piS <- fread(paste0(outdir_dat, "/ag_age_S_dt.csv"),
                       select = c("gene_family", "pi_s"))

# add to main dt
ag_genes <- merge(ag_genes, ag_piS,
                    by = "gene_family")

# keep only the AGs analysed
set1_ags <- fread("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/set_1_names.csv")$gene_family

ag_genes <- ag_genes[gene_family %chin% set1_ags]

# split AG anchors
n <- lengths(strsplit(ag_genes$anchor, ":", fixed = TRUE))

ag_long <- ag_genes[rep(seq_len(.N), n)][
  , anchor := unlist(strsplit(ag_genes$anchor, ":", fixed = TRUE))]

# merge the two together
core_genes_sub <- copy(core_genes)

core_genes_sub[,c("gene_family", "consensus_product") := NULL]

setnames(core_genes_sub, "pi_s", "core_pi_s")

ag_long <- merge(ag_long, core_genes_sub,
                  all.x = TRUE, by = c("msu", "anchor"))


# Calculate difference ----------------------------------------------------

delta = ag_long[,.(
  pi_s = mean(pi_s),
  core_pi_s = mean(core_pi_s)
), gene_family][, delta := pi_s - core_pi_s]



# Get GC content ----------------------------------------------------------
# get locations
full_anno <- unique(fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"),
                          select = c("gene_family","geno_id", "seqnames","fstart", "fend", "fstrand",  
                                     "anchor", "msu", "ag_type")))[msu!=""]
setnames(full_anno, c("fstart","fend","fstrand"), c("start","end","strand"))


fna_files <-  list.files("C:/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_260_chr_fasta/",
                         full.names = TRUE, pattern = ".fasta")

msus <- unique(full_anno[order(msu),msu])

# step 1, get core GC per genome per msu, and account for msu variation
core_seqs <- full_anno[ag_type == "core", .(geno_id, msu, seqnames,start, end, strand)]

list_cores <- split(core_seqs, by = c("geno_id", "msu"))

# Set up parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Paralleled loop to read sequences and compute consensus
# Per isolate, per msu, concatenate core gene sequences, compute GC
core_msu_gc <- foreach(f = fna_files, 
                       .combine = rbind, 
                       .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                         
                         genome_id = gsub(".fasta", "", basename(f))
                         
                         # Read aligned sequences
                         genome <- readDNAStringSet(f)
                         
                         # get core list names
                         msu_lst_names <- names(list_cores)[grepl(genome_id, names(list_cores))]
                         
                         geno_gc_perc <- lapply(msu_lst_names, function(msu_x){
                           temp_core_dt <- list_cores[[msu_x]]
                           
                           genes_gr <- dt2gr(temp_core_dt[,.(seqnames, start, end, 
                                                             strand)])
                           
                           temp_core_seqs <- getSeq(genome, genes_gr)
                           
                           gc_perc = sum(letterFrequency(temp_core_seqs, letters = "GC"))/
                             sum(width(temp_core_seqs))
                           
                           data.table(id = msu_x, gc_perc)
                         })
                         
                         rbindlist(geno_gc_perc)
                      

                    }

# Stop cluster after execution
stopCluster(cl)

core_msu_gc[, c("geno_id", "msu") := tstrsplit(id, ".", fixed=TRUE, keep = 1:2)]
core_msu_gc[,id:=NULL]

# Step 2, get CG perc per AG
ag_seqs <- full_anno[gene_family %chin% set1_ags,
                     .(gene_family,geno_id,seqnames,start,end,strand,msu)]

list_ags <- split(ag_seqs, by = c("geno_id"))

# Set up parallel processing
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Paralleled loop to read sequences and compute consensus
ag_gc <- foreach(f = fna_files, 
                       .combine = rbind, 
                       .packages = c("Biostrings", "BSgenome", "data.table", "gUtils")) %dopar% {
                         
                         genome_id = gsub(".fasta", "", basename(f))
                         
                         # Read aligned sequences
                         genome <- readDNAStringSet(f)
                         
                         temp_ag_dt <- list_ags[[genome_id]]
                         
                         genes_gr <- dt2gr(temp_ag_dt[,.(seqnames, start, end, 
                                                         strand)])
                         
                         temp_ag_seqs <- getSeq(genome, genes_gr)
                           
                         gc_perc = letterFrequency(temp_ag_seqs, letters = "GC", as.prob = TRUE)
                           
                         data.table(temp_ag_dt[,.(geno_id, msu, gene_family)], gc_perc)
                       }

# Stop cluster after execution
stopCluster(cl)

setnames(ag_gc,"G|C", "gc_perc")

# Calculate AG GC zscores -------------------------------------------------
# Global core GC distribution across all MSUs
global_ref <- core_msu_gc[, .(
  msu_gc = mean(gc_perc)
), by = msu]

global_mean <- mean(global_ref$msu_gc)
global_sd   <- sd(global_ref$msu_gc)

# Z-score: how far is this AG from its local MSU mean,
# scaled by genome-wide positional GC variation
ag_comp <- merge(ag_gc, 
                 core_msu_gc[, .(core_msu_gc = mean(gc_perc)), by = msu],
                 by = "msu")
ag_comp[, gc_zscore := (gc_perc - core_msu_gc) / global_sd]# use globals msu varaince for sensible varaince, within core msu has no variance


# Compare the two metrics -------------------------------------------------

family_gc <- ag_comp[, .(
  gc_range  = max(gc_zscore) - min(gc_zscore),
  freq = .N
), by = gene_family]

# Merge
family_2d <- merge(family_gc, delta[,.(gene_family, delta)], by = "gene_family")

# Expected range under normality scales as ~2*qnorm(1 - 1/(2n))*sigma
# Studentised range corrects for this
# family_2d[, expected_range_factor := 2 * qnorm(1 - 1 / (2 * freq))]
# family_2d[, gc_range_adj := gc_range / expected_range_factor]

# gc_thresh <- quantile(family_2d$gc_range_adj, 0.95)
gc_thresh <- quantile(family_2d$gc_range, 0.95)
delta_thresh <- quantile(family_2d$delta, 0.95)

family_2d[, outlier_gc := gc_range > gc_thresh]
family_2d[, outlier_delta := delta > delta_thresh]
family_2d[, outlier_both := outlier_gc & outlier_delta]
family_2d[, outlier_either := outlier_gc | outlier_delta]

family_2d[, hgt_class := fcase(
  outlier_gc & outlier_delta, "strong",
  outlier_gc,                 "gc_only",
  outlier_delta,              "delta_only",
  default =                   "none"
)]

table(family_2d$hgt_class)

col_map <- c(
  "none"       = adjustcolor("steelblue", 0.3),
  "gc_only"    = adjustcolor("darkorange", 0.7),
  "delta_only" = adjustcolor("purple4", 0.7),
  "strong"     = adjustcolor("firebrick", 0.8)
)



# output ------------------------------------------------------------------
fwrite(family_2d, paste0(outdir_dat, "/family_2d.csv"))


# Plot --------------------------------------------------------------------

png(paste0(outdir_fig,"/find_hgt_outliers.png"),
    width = 22, height = 17,  units = "cm", res = 300,
    pointsize = 11, type = "cairo")

with(family_2d,
     plot(delta, gc_range,
          pch = 16,
          cex = sqrt(freq / max(freq)) * 3 + 0.3,
          col = col_map[hgt_class],
          xlab = expression(delta ~ (AG ~ pi[S] - "flanking core" ~ pi[S])),
          ylab = expression("GC deviation range ("*z[AG] ~ "vs MSU core)"),
          main = "Multi-introduction candidates (95th percentile thresholds)"))
abline(h = gc_thresh, lty = 2, col = "grey40")
abline(v = delta_thresh, lty = 2, col = "grey40")
size_vals <- c(10, 50, 150, 257)
size_cex <- sqrt(size_vals / max(family_2d$freq)) * 3 + 0.3
legend("topleft",
       legend = c(
         paste0("Strong (", sum(family_2d$hgt_class == "strong"), ")"),
         paste0("GC only (", sum(family_2d$hgt_class == "gc_only"), ")"),
         paste0("Delta only (", sum(family_2d$hgt_class == "delta_only"), ")"),
         paste0("None (", sum(family_2d$hgt_class == "none"), ")"),
         "", "Allele count:",
         paste0("n = ", size_vals)
       ),
       pch = c(rep(16, 4), NA, NA, rep(16, length(size_vals))),
       pt.cex = c(rep(1.5, 4), 0, 0, size_cex),
       col = c(col_map[c("strong", "gc_only", "delta_only", "none")],
               NA, NA,
               rep(adjustcolor("grey40", 0.6), length(size_vals))),
       bty = "n",
       cex = 0.8,
       y.intersp = 1.3)


# how does this affect our project outliers?
selected_ags <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))[p_bh <= 0.05]$gene_family
sel_cands_dt <- fread(paste0(outdir_dat, "/gene_pvals_fdr.csv"))
sel_cands_dt[, sel := fcase(seg_s_sites>expected_S & p_bh <=0.05, "bal",
                            seg_s_sites < expected_S & p_bh <=0.05, "dir",
                            default = "unsel")]
sel_cands_dt <- sel_cands_dt[,.(gene_family, sel)]
sel_cands_dt <- merge(sel_cands_dt, family_2d[,.(gene_family, hgt_class)],
                      all.x = TRUE, by = "gene_family")
sel_cands_dt[is.na(hgt_class), hgt_class := "none"]
lost_sel_outs <- dcast(sel_cands_dt, sel ~ hgt_class, fun.aggregate = length)
lost_sel_outs[, removed := sum(delta_only,gc_only,strong)/sum(delta_only,gc_only,none,strong), .I]
lost_sel_outs$removed = paste0(sprintf("%2.2f",(lost_sel_outs$removed*100)),"%")
tbl <- as.data.frame(lost_sel_outs[, .(sel, delta_only, gc_only, none, strong, removed)])
addtable2plot(
  x = par("usr")[2],
  y = par("usr")[4],
  table = tbl,
  bty = "n",
  display.rownames = FALSE,
  hlines = TRUE,
  vlines = FALSE,
  cex = 0.65,
  xjust = 1.05,
  yjust = -0.5
)

# key outliers
# Key examples per quadrant
labels <- data.table(
  gene_family = c(
    # GC only
    "g006249","g006368", "g000747", "g000990", "g001058",
    # Strong
    "g000516_3", "g006375", "g007767","g002399_3","g000140_2",
    # Delta only
    "g006323", "g002075", "g001847", "g002467",
    # None
    "g000145", "g000264", "g000378", "g003288"
  ),
  label = c(
    "cdiA toxin(3)", "phage protein (3)", "lpdA (34)", "levR (29)", "chrA (27)",
    "yfmG (11)", "integrase (3)", "xerC (2)","phage-related lysozyme (2)","Ig-like protein (2)",
    "fimB (3)", "yitT (18)", "fimD (17)", "integrase (15)",
    "ompA (74)", "wcaG (56)", "tssK (47)", "hpxX metabolic (11)"
  )
)

# Merge coordinates
labels <- merge(labels, family_2d[, .(gene_family, delta, gc_range, hgt_class)],
                by = "gene_family")

# Highlight labelled points
text(labels$delta, labels$gc_range,
     labels = labels$label,xpd=TRUE,
     pos = 4, cex = 0.6, offset = 0.7,
     col = "grey20")

dev.off()


