# get core alignment and tree


# Get ska2 alignment ------------------------------------------------------
# code to list all .fasta files
# for f in *.fasta; do
# base=$(basename "$f" .fasta)
# echo -e "${base}\t${f}"
# done > input.list

# Detect and remove areas of recombination --------------------------------

system("wsl conda run -n gubbins_env run_gubbins.py /mnt/c/path/to/alignment.fasta")
