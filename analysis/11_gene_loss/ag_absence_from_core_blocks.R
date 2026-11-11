# Infer AG absence at its syntenic locus from pangraph core-block anchors
# For each set 1 AG (n = 5270) and each genome, find the core-block junction
# (Ai:Aj) that holds the AG in carriers, check the junction is intact in the
# genome, and record what sits there when the AG is absent.
# NOTE: an empty site alone does not separate loss from never gained. Map the
# status onto the core tree to make that call (absent lineages nested within
# carriers = loss).

# Parameters --------------------------------------------------------------

rearr_mult   <- 5      # junction > rearr_mult x max carrier junction = rearranged
rearr_min_bp <- 20000  # ... and at least this long
nbr_prop     <- 0.5    # neighbour = family in the junction of >= 50% of carriers
remnant_prop <- 0.5    # junction shortfall < 50% of AG length = possible remnant

# Import data -------------------------------------------------------------

set_1_ags <- fread(paste0(outdir_dat, "/set_1_names.csv"))$gene_family

msu_anchored <- fread(paste0(outdir_dat, "/msu_regions_anchored.csv"),
                      select = c("fus_locus_tag", "gene_family", "geno_id", 
                                 "ag_type", "freq", "start", "end", "tot_len", 
                                 "anchor", "msu"))

geno_ids <- unique(msu_anchored$geno_id)

# AG locus: modal msu and junction among carriers --------------------------

ag_loci <- msu_anchored[gene_family %chin% set_1_ags & grepl(":", anchor),
                        .(n_car_locus = .N), by = .(gene_family, msu, anchor)]

ag_loci <- ag_loci[order(-n_car_locus)][!duplicated(gene_family)]

ag_loci[, c("anc_left", "anc_right") := tstrsplit(anchor, ":", fixed = TRUE)]

# AGs without a resolved junction in any carrier. These sit outside every MSU
# (msu is NA), e.g. ahead of the first or after the last anchor of an MSU,
# where no Ai:Aj junction is defined. They cannot be assessed here.
no_locus <- setdiff(set_1_ags, ag_loci$gene_family)
length(no_locus)

# median AG length among carriers
ag_len <- msu_anchored[gene_family %chin% set_1_ags, 
                       .(ag_len = median(abs(end - start) + 1)), by = gene_family]

# Anchor (core gene) coordinates per genome ---------------------------------

anchor_coords <- msu_anchored[grepl("^A[0-9]+$", anchor),
                              .(geno_id, msu, anchor, 
                                a_start = start, a_end = end, tot_len)]

# anchors are single copy and present in all genomes
stopifnot(anchor_coords[, .N, by = .(geno_id, msu, anchor)][, all(N == 1)])

# Junction contents per genome ----------------------------------------------
# junctions are only populated where both anchors were found; an intact,
# empty junction has no rows here

jun_members <- unique(msu_anchored[grepl(":", anchor),
                                   .(geno_id, msu, anchor, member = gene_family,
                                     ag_type, len = abs(end - start) + 1)])

jun_contents <- jun_members[, .(n_genes  = .N,
                                annot_bp = sum(len),
                                families = paste(sort(unique(member)), collapse = ",")),
                            by = .(geno_id, msu, anchor)]

# Grid: every AG x every genome ---------------------------------------------

ag_grid <- ag_loci[, .(geno_id = geno_ids), 
                   by = .(gene_family, msu, anchor, anc_left, anc_right)]

# presence of the AG (anywhere) and where it was placed
ag_presence <- unique(msu_anchored[gene_family %chin% set_1_ags, 
                                   .(gene_family, geno_id, loc_anchor = anchor)])
ag_presence <- ag_presence[!duplicated(ag_presence[, .(gene_family, geno_id)])]

ag_grid[, present := FALSE]
ag_grid[ag_presence, on = .(gene_family, geno_id), 
        `:=`(present = TRUE, loc_anchor = i.loc_anchor)]

# left and right anchor coordinates
ag_grid[anchor_coords, on = .(geno_id, msu, anc_left = anchor),
        `:=`(l_start = i.a_start, l_end = i.a_end, tot_len = i.tot_len)]

ag_grid[anchor_coords, on = .(geno_id, msu, anc_right = anchor),
        `:=`(r_start = i.a_start, r_end = i.a_end)]

# inner distance between anchors, allowing for junctions across the origin
ag_grid[, jun_len := pmax(l_start, r_start) - pmin(l_end, r_end) - 1]

ag_grid[!is.na(jun_len) & jun_len > tot_len/2, 
        jun_len := tot_len - (pmax(l_end, r_end) - pmin(l_start, r_start) + 1)]

ag_grid[jun_len < 0, jun_len := 0]

# junction contents
ag_grid[jun_contents, on = .(geno_id, msu, anchor),
        `:=`(n_genes = i.n_genes, annot_bp = i.annot_bp, families = i.families)]

ag_grid[is.na(n_genes), `:=`(n_genes = 0L, annot_bp = 0L, families = "")]

ag_grid[ag_len, on = .(gene_family), ag_len := i.ag_len]

# Carrier reference values ------------------------------------------------

ag_grid[, at_locus := present & !is.na(loc_anchor) & loc_anchor == anchor]

carrier_ref <- ag_grid[at_locus == TRUE, 
                       .(car_jun_med = median(jun_len, na.rm = TRUE),
                         car_jun_max = max(jun_len, na.rm = TRUE)), 
                       by = gene_family]

ag_grid[carrier_ref, on = .(gene_family), 
        `:=`(car_jun_med = i.car_jun_med, car_jun_max = i.car_jun_max)]

# Neighbours: accessory families sharing the junction with the AG in carriers
# junctions also hold non-anchor core genes, so these are excluded

acc_members <- unique(jun_members[ag_type != "core", .(geno_id, msu, anchor, member)])

car_members <- merge(ag_grid[at_locus == TRUE, .(gene_family, geno_id, msu, anchor)],
                     acc_members,
                     by = c("geno_id", "msu", "anchor"), allow.cartesian = TRUE)

n_car <- ag_grid[at_locus == TRUE, .(n_car = .N), by = gene_family]

nbr_dt <- car_members[member != gene_family, .(n = .N), by = .(gene_family, member)]
nbr_dt[n_car, on = .(gene_family), prop := n / i.n_car]
nbr_dt <- nbr_dt[prop >= nbr_prop, .(gene_family, member)]

# which neighbours remain in genomes lacking the AG
absent_nbrs <- merge(ag_grid[present == FALSE, .(gene_family, geno_id, msu, anchor)],
                     nbr_dt, by = "gene_family", allow.cartesian = TRUE)

absent_nbrs[, nbr_present := FALSE]
absent_nbrs[acc_members, on = .(geno_id, msu, anchor, member), nbr_present := TRUE]

nbr_summary <- absent_nbrs[, .(n_nbrs = .N, n_nbrs_present = sum(nbr_present)),
                           by = .(gene_family, geno_id)]

ag_grid[nbr_summary, on = .(gene_family, geno_id),
        `:=`(n_nbrs = i.n_nbrs, n_nbrs_present = i.n_nbrs_present)]

ag_grid[is.na(n_nbrs), `:=`(n_nbrs = 0L, n_nbrs_present = 0L)]

# accessory genes at the junction that are not the AG's usual neighbours
# (replacement by other content)
other_acc <- merge(ag_grid[present == FALSE, .(gene_family, geno_id, msu, anchor)],
                   acc_members, by = c("geno_id", "msu", "anchor"),
                   allow.cartesian = TRUE)

other_acc <- other_acc[!nbr_dt, on = .(gene_family, member)]

ag_grid[other_acc[, .(n_other_acc = .N), by = .(gene_family, geno_id)],
        on = .(gene_family, geno_id), n_other_acc := i.n_other_acc]

ag_grid[is.na(n_other_acc), n_other_acc := 0L]

# Assign status -----------------------------------------------------------

ag_grid[, unannot_bp := pmax(jun_len - annot_bp, 0)]

ag_grid[, rearr_cut := pmax(rearr_mult * car_jun_max, rearr_min_bp)]

ag_grid[, status := fcase(
  at_locus,                                         "present_at_locus",
  present & (is.na(loc_anchor) | loc_anchor == ""), "present_locus_unresolved",
  present,                                          "present_other_locus",
  is.na(l_start) | is.na(r_start),                  "anchor_missing",
  jun_len > rearr_cut,                              "junction_rearranged",
  n_genes == 0L,                                    "empty_site",
  n_nbrs_present == n_nbrs,                         "ag_only_missing",
  n_nbrs_present == 0L,                             "ag_block_missing",
  default =                                         "partial_block_missing")]

# replacement: other accessory content at the site the AG occupies in carriers
ag_grid[, replaced := !present & n_other_acc > 0L]

# bp missing relative to carriers. Only interpretable where the neighbours are
# intact (ag_only_missing): a shortfall close to ag_len = clean deletion; a
# shortfall well below ag_len = sequence still there but not called as the AG
# (candidate pseudogene / unannotated remnant, check with BLAST)
ag_grid[, len_shortfall := car_jun_med - jun_len]

ag_grid[, possible_remnant := status == "ag_only_missing" &
          len_shortfall < remnant_prop * ag_len]

ag_grid[, .N, by = status][order(-N)]

# Per AG summary ----------------------------------------------------------

ag_absence_summary <- dcast(ag_grid, gene_family ~ status, 
                            value.var = "geno_id", fun.aggregate = length)

ag_absence_summary[ag_grid[, .(n_remnant  = sum(possible_remnant),
                               n_replaced = sum(replaced)), by = gene_family],
                   on = .(gene_family),
                   `:=`(n_remnant = i.n_remnant, n_replaced = i.n_replaced)]

# AGs that could not be placed at a junction
ag_absence_summary <- rbind(ag_absence_summary,
                            data.table(gene_family = no_locus), fill = TRUE)

ag_absence_summary[, locus_resolved := !gene_family %chin% no_locus]

ag_absence_summary[unique(msu_anchored[, .(gene_family, freq, ag_type)]), 
                   on = .(gene_family), `:=`(freq = i.freq, ag_type = i.ag_type)]

ag_absence_summary[ag_loci, on = .(gene_family), `:=`(msu = i.msu, anchor = i.anchor)]

# Output ------------------------------------------------------------------

fwrite(ag_grid[, !c("rearr_cut")], paste0(outdir_dat, "/ag_absence_by_genome.csv"))
fwrite(ag_absence_summary, paste0(outdir_dat, "/ag_absence_summary.csv"))
