#!/bin/bash

# Define paths
src_dir="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/kpne_1695_fasta"
dst_base="/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/split_fastas"

# Number of parts
n_parts=8

# Make 8 subdirectories
for i in $(seq 1 $n_parts); do
  mkdir -p "$dst_base/part_$i"
done

# Get sorted list of FASTA files
files=($(find "$src_dir" -maxdepth 1 -type f -name "*.fasta" | sort))

# Total number of files
total=${#files[@]}
per_dir=$(( (total + n_parts - 1) / n_parts ))  # round up

# Distribute files
for i in "${!files[@]}"; do
  dir_index=$(( i / per_dir + 1 ))
  cp "${files[$i]}" "$dst_base/part_$dir_index/"
done
