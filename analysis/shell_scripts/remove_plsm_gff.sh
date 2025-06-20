#!/bin/bash

# Path to plasmid contig list
PLASMID_LIST="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/plasmid_contigs.txt"

# Directory containing the original GFF files
GFFDIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_gffs"

# Output directory for filtered GFFs
OUTDIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_p_filtered_gffs"
mkdir -p "$OUTDIR"

# Loop through unique GFF files in the list
for GFF in $(cut -f1 "$PLASMID_LIST" | sort -u); do
    echo "Processing $GFF"

    # Full path to input GFF
    INFILE="$GFFDIR/$GFF"
    BASENAME=$(basename "$GFF" .gff)
    OUTFILE="$OUTDIR/${BASENAME}_filtered.gff"

    # Get list of contigs to remove for this GFF
    awk -v gff="$GFF" '$1 == gff {print $2}' "$PLASMID_LIST" > tmp_contigs.txt

    # Filter GFF: keep headers and non-plasmid contigs
    awk '
        BEGIN {
            while ((getline < "tmp_contigs.txt") > 0) contigs[$1] = 1
        }
        /^#/ { print; next }
        !($1 in contigs) { print }
    ' "$INFILE" > "$OUTFILE"

    rm tmp_contigs.txt
done
