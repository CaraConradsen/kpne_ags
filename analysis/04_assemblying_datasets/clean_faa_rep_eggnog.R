# export clean faa for eggnog mapper

rep_faa_dir <- "./input_data/PIRATE_260_hybrid_chr_out/representative_sequences.faa"

rep_faa <- readAAStringSet(rep_faa_dir)

# clean names

names(rep_faa) <- sub(";.*$", "", names(rep_faa))

writeXStringSet(rep_faa, "./input_data/eggnog/klebsiella_proteins.faa")