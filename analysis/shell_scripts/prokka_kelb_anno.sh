#!/bin/bash

# Define the directory where your genome files are located
GENOME_DIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kleb_AG_mobilisation/data/temp_kleb_data"
OUTPUT_DIR="/mnt/c/Users/carac/Dropbox/Vos_Lab/kleb_AG_mobilisation/output/data/prokka_anno"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through each genome file in the specified directory
for GENOME in $GENOME_DIR/*.fasta; do
    # Extract the filename without the directory and extension to use as a prefix
    GENOME_NAME=$(basename "$GENOME" .fasta)
    
    # Create a directory for Prokka output for this genome
    GENOME_OUTPUT_DIR="$OUTPUT_DIR/$GENOME_NAME"
    mkdir -p "$GENOME_OUTPUT_DIR"
    
    # Run Prokka with specified arguments and save output in the created directory
    # -p <prefix>: Use the genome name as the prefix
    # --cpus 7: Use 7 CPU cores in parallel
    # --outdir: Specify the output directory
    # --force: Overwrite output if it already exists
    # --addgenes: Add extra genes (for example, extra virulence or antibiotic resistance genes)
    # --metagenome: Useful if you're working with metagenomic data
    prokka \
        --outdir "$GENOME_OUTPUT_DIR" \
        --prefix "$GENOME_NAME" \
        --cpus 7 \
        --force \
        --addgenes \
	--locustag "K_" \
        --kingdom Bacteria \
        --genus Klebsiella \
        "$GENOME" &
    
    # Ensure parallel execution (run each genome analysis in the background)
done

# Wait for all background processes to finish before exiting the script
wait

echo "Prokka annotation completed for all genomes."
