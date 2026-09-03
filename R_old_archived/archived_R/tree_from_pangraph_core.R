
# 2. ### Generate clonalframe for trees using pangraph and gubbins

# get core
pangraph_input = file.path(input_dir, "graph.json")
set_ref_strain = "SPARK_1006_C1"
core_output_file = file.path(input_dir, "core_genome_aln.fa")# output file dir

# pangraph core command string 
cmd_core <- sprintf(
  "wsl %s export core-genome %s --guide-strain %s -o %s",
  pangraph_path,
  pangraph_input,
  set_ref_strain,
  core_output_file
)

# check command
cat("Running command:\n", cmd_core, "\n\n")

# run pangraph 
system(cmd_core)


# 3. ### get clean clonal frame with gubbins (in future maybe try STs?)

conda_path = "/home/carac/anaconda3/bin/conda" # path to conda inside WSL
conda_env = "gubbins_env" # conda environment name
iterations = 20
threads = 22 # number of threads
tree_method = "fasttree" # tree builder to use
# input alignment file (from pangraph output)
input_alignment = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/core_genome_aln.fa"
# output prefix 
output_prefix = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/input_data/pangraph/gub_graph"

# build gubbins command string 
cmd_gubbins <- sprintf(
  "wsl %s run -n %s run_gubbins.py --iterations %d --threads %d -t %s --prefix %s %s",
  conda_path,
  conda_env,
  iterations,
  threads,
  tree_method,
  output_prefix,
  input_alignment
)

# Print command for verification 
cat("Running command:\n", cmd_gubbins, "\n\n")

start <- Sys.time()# Start timer

# run pangraph 
system(cmd_gubbins)

end <- Sys.time()# End timer
print(end - start) # Print runtime
# 10 seqs Time difference of 1.921924 mins
# 25 seqs Time difference of 2.211768 mins
# 93 seqs Time difference of 25.72868 mins
# 115 seqs Time difference of 25.72868 mins
# 214 seqs Time difference of 2.804536 hours
# 260 seqs Time difference of 5.134598 hours
