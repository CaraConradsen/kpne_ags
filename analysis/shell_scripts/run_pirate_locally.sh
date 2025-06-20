#!/bin/bash 

# Set GFF input and output directories
GFF_DIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/plsm_filtered_gffs"
OUTPUT_DIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/PIRATE_i95_noplsm_1695_out"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Activate the Conda environment
source ~/anaconda3/etc/profile.d/conda.sh  # needed for conda activation in scripts
conda activate pirate_env

# Run PIRATE with specified options
PIRATE -i "$GFF_DIR" -s "95,96,97,98,99,100" -k "--cd-low 100 --cd-mem 35840" -t 6 -o "$OUTPUT_DIR"
