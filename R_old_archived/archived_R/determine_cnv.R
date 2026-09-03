# isolates re-run with plasmids, to determine if genes have cnv
# annolyingly the annotations slipped, so have to do a round-about merge

# PART 1. Import gene annotations -------------------------------------------------
plasmid_cnv_gffs_files <- list.files("./input_data/PIRATE_260_plasmid_cnv_hps0.6_out/modified_gffs/",
                                full.names = TRUE, recursive = TRUE, 
                                pattern = ".gff")


cl <- makePSOCKcluster(num_cores)
registerDoParallel(cl)

plasmid_cnv_anno <- foreach(i = 1:length(plasmid_cnv_gffs_files), 
                       .packages = c("rtracklayer", "gUtils"),
                       .combine = "rbind") %dopar% {
                         
                         geno_id = gsub(".gff","", basename(plasmid_cnv_gffs_files[i]))
                         
                         # get pirate anno
                         anno <- rtracklayer::import.gff(plasmid_cnv_gffs_files[i])
                         
                         anno_dt <- gr2dt(anno)
                         
                         anno_dt$geno_id = geno_id
                         
                         anno_dt$gene_length = width(anno)
                         
                         # return dat
                         anno_dt[type == "CDS", .(geno_id, seqnames, start,end, strand,gene_length, 
                                                  locus_tag,gene, product)]
                         
                         
                       }

stopCluster(cl)

# add pirate & core syntenic groups -----------------------------------------
focal_genome_names = unique(plasmid_cnv_anno$geno_id)
plasmid_cnv_lng <- fread("./input_data/PIRATE_260_plasmid_cnv_hps0.6_out/PIRATE.gene_families.ordered.tsv",
                    select = c("gene_family", "consensus_gene_name",
                               "consensus_product", focal_genome_names))

plasmid_cnv_lng <- melt(plasmid_cnv_lng, 
                   id.vars = c("gene_family", 
                               "consensus_gene_name",
                               "consensus_product"),
                   variable.name = "geno_id", 
                   value.name = "fus_locus_tag"
)

plasmid_cnv_lng <- plasmid_cnv_lng[fus_locus_tag!=""]

plasmid_cnv_lng$locus_tag = plasmid_cnv_lng$fus_locus_tag

# tidy fusion brackets
plasmid_cnv_lng$locus_tag = gsub("[()]", "", plasmid_cnv_lng$locus_tag)

#tidy fission pirate loci
plasmid_cnv_lng <- plasmid_cnv_lng[
  ,
  .(locus_tag = unlist(strsplit(locus_tag, ";"))),
  by = setdiff(names(plasmid_cnv_lng), "locus_tag")
]

#tidy fusion pirate loci
plasmid_cnv_lng <- plasmid_cnv_lng[
  ,
  .(locus_tag = unlist(strsplit(locus_tag, ":"))),
  by = setdiff(names(plasmid_cnv_lng), "locus_tag")
]

# add in annotation, removing all short genes not analysed in pirate
plasmid_cnv_anno <- merge(plasmid_cnv_lng, plasmid_cnv_anno, 
                     all.x = TRUE, by = c("geno_id","locus_tag"))


# fix gene_family name

setnames(plasmid_cnv_anno, "gene_family", "cnv_gene_family")

# Add in our chromosomal information 

pirate_anno <- fread(paste0(outdir_dat, "/all_pirate_anno_full.csv"), 
                            select = c("seqnames","strand", "start","end", "gene_family"))


# chr contigs
chr_contig = unique(pirate_anno$seqnames)

# id plasmid contigs
plasmid_cnv_anno[, replicon := fifelse(seqnames %in% chr_contig,
                                       "chromosome", "plasmid")]

# get cnv info per locus
cnv_info = plasmid_cnv_anno[replicon=="chromosome", .(locus_tag, cnv_gene_family, geno_id, consensus_product)]

# use fus_log, don't count twice

plsm_counts = unique(plasmid_cnv_anno[replicon != "chromosome", 
                                      .(geno_id, fus_locus_tag, cnv_gene_family)])

plsm_counts = plsm_counts[,.(cnv = .N), by = c("geno_id", "cnv_gene_family")]


cnv_info = merge(cnv_info, plsm_counts,
                 all.x = TRUE, by =c("geno_id", "cnv_gene_family"))

# set na to zero
cnv_info[is.na(cnv), cnv:= 0]# by locus, so cnv entry is duplicated for fusion



# id focal chr genes ------------------------------------------------------

plasmid_cnv_anno_mrg_dat <- plasmid_cnv_anno[,.(seqnames,start,end,strand,locus_tag)]

# merge messy
plas_set <- copy(plasmid_cnv_anno_mrg_dat)
chr_set <- copy(pirate_anno)


plas_set <- plas_set[, .(
  locus_tag,
  pseqnames =seqnames,
  pstrand = strand,
  pstart = start,
  pend   = end

)]

setkey(plas_set, pseqnames, pstrand, pstart, pend)

chr_set <- chr_set[, .(
  gene_family,
  cseqnames = seqnames,
  cstrand = strand,
  cstart = start,
  cend   = end
)]

setkey(chr_set, cseqnames, cstrand, cstart, cend)

cnv_matches <- foverlaps(
  plas_set,
  chr_set,
  by.x = c("pseqnames", "pstrand", "pstart", "pend"),
  by.y = c("cseqnames", "cstrand", "cstart", "cend"),
  type = "any",  
  nomatch = NA
)

cnv_matches = cnv_matches[!is.na(gene_family),.(locus_tag, gene_family)]


# put everything together -------------------------------------------------

cnv_info <- merge(cnv_matches, cnv_info,
                  all.x = TRUE, by =c("locus_tag"))

cnv_info[, locus_tag:= NULL]

# add back to out chromosomal data

chr_anno <- fread(paste0(outdir_dat, "/pangenome_anno.csv"),
                  select = c("gene_family", "geno_id", 
                             "number_genomes", "consensus_product"))

chr_anno <- merge(chr_anno, cnv_info,
                  all = TRUE, by = c("gene_family", "geno_id"))


cnv_info <- cnv_info[, .(cnv = mean(cnv)), by = c("gene_family", "geno_id") ]

setnames(cnv_info, "cnv", "cnv_on_plsmd")

cnv_info[, cnv_on_plsmd := as.integer(cnv_on_plsmd > 0)]

fwrite(cnv_info, paste0(outdir_dat, "/cnv_info.csv"))



