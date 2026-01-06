# Will need:
#   1. the recombination free core genome phylogeny
#   2. a binary presence / absence file


# 1. Convert gubbins final tre file to newick -----------------------------
core_gub_tree <- read.tree("./input_data/test_pangraph/gub_graph.node_labelled.final_tree.tre") 



conda_path = "/home/carac/anaconda3/bin/conda" # path to conda inside WSL
conda_env = "gubbins_env" # conda environment name
threads = 12 # number of threads
tree_method = "fasttree" # tree builder to use
# input alignment file (from pangraph output)
input_alignment = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/test_pangraph/core_genome_aln.fa"
# output prefix 
output_prefix = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/test_pangraph/gub_graph"

# build gubbins command string 
cmd_gubbins <- sprintf(
  "wsl %s run -n %s run_gubbins.py --threads %d -t %s --prefix %s %s",
  conda_path,
  conda_env,
  threads,
  tree_method,
  output_prefix,
  input_alignment
)

# ---- Print command for verification ----
cat("Running command:\n", cmd_gubbins, "\n\n")

start <- Sys.time()# Start timer

# run pangraph 
system(cmd_gubbins)

end <- Sys.time()# End timer
print(end - start) # Print runtime
# 10 seqs Time difference of 1.921924 mins
# 25 seqs Time difference of 2.211768 mins
# 93 seqs Time difference of 25.72868 mins

aln_file = res_dict["aln_f"]
aln_L = res_dict["aln_L"]

# instantiate treetime
myTree = treetime.TreeAnc(
  gtr="Jukes-Cantor", tree=tree, aln=aln_file, verbose=0, seq_len=aln_L
)
myTree.tree.root.branch_length = 0.0
myTree.infer_ancestral_sequences(prune_short=True)