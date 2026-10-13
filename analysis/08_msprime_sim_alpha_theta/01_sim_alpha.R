# library(reticulate)
# # conda_create("msprime", packages = c("msprime", "tskit"), channel = "conda-forge")
# use_condaenv("msprime", required = TRUE)
# py_config()

# set up alpha sim --------------------------------------------------------
mr <- fread(paste0(outdir_dat, "/core_mut_rate.csv"))[class == "unweighted_mean"]

msprime_env = "/home/carac/anaconda3/envs/msprime/bin/python"
py_file = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/analysis/python/determine_alpha.py"
sample_n = 260
pi_obs = mr$mean_pi_s
thetaw_obs = mr$mean_thetaW_s
reps = 10000
out_file = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/alpha_theta_curve_test.csv"
alpha_out = "/mnt/c/Users/carac/Dropbox/Vos_Lab/kpne_ags/output/data/alpha_hat.csv"

# alpha msprime command string 
cmd_alpha_msprime <- sprintf(
  "wsl %s %s --n %d --pi-obs %.6g --thetaw-obs %.6g --reps %d --out %s --alpha-out %s",
  msprime_env,
  py_file,
  sample_n,
  pi_obs,
  thetaw_obs,
  reps, 
  out_file,
  alpha_out
)

# check command
cat("Running command:\n", cmd_alpha_msprime, "\n\n")


start <- Sys.time()# Start timer

# run msprime alpha sim 
system(cmd_alpha_msprime)

end <- Sys.time()
print(end - start) # Time difference of 54.30717 secs
