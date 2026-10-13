# library(reticulate)
# # conda_create("msprime", packages = c("msprime", "tskit"), channel = "conda-forge")
# use_condaenv("msprime", required = TRUE)
# py_config()

# set up theta sim --------------------------------------------------------
msprime_env = "/home/carac/anaconda3/envs/msprime/bin/python"
py_file = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/determine_theta.py"
alpha_hat = fread(paste0(outdir_dat, "/alpha_hat.csv"))$alpha_hat
pi_obs = fread(paste0(outdir_dat, "/alpha_hat.csv"))$pi_obs
thetaw_obs = fread(paste0(outdir_dat, "/alpha_hat.csv"))$thetaw_obs 
reps = 10000
out_file = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/theta_pi_curve.csv"
theta_out = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/theta_hat.csv"

# theta msprmine command string
cmd_theta_msprime <- sprintf(
  "wsl %s %s --alpha %.6g --pi-obs %.6g --thetaw-obs %.6g --reps %d --out %s --theta-out %s",
  msprime_env,
  py_file,
  alpha_hat,
  pi_obs,
  thetaw_obs,
  reps,
  out_file,
  theta_out
)

# check command
cat("Running command:\n", cmd_theta_msprime, "\n\n")


start <- Sys.time()# Start timer

# run msprime alpha sim
system(cmd_theta_msprime)

end <- Sys.time()
print(end - start) # Time difference of 2.366954 mins
