#!/bin/bash
#SBATCH --job-name=run_rcoal_sim
#SBATCH --nodes=1
#SBATCH --partition=hmq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4       # adjust as needed
#SBATCH --mem=120G               # adjust as needed
#SBATCH --time=7:00:00         # adjust as needed
#SBATCH --account=m.vos 
#SBATCH --output=run_rcoal_sim_%j.out
#SBATCH --error=run_rcoal_sim_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cc1376@exeter.ac.uk


cd /nobackup/beegfs/workspace/cc1376/run_rcoal || exit 1
export OUTDIR_DAT=/nobackup/beegfs/workspace/cc1376/run_rcoal

# Load R
module load R/4.4.0
Rscript -e 'stopifnot(getRversion() >= "4.4.0")' || exit 1

# Run the script
Rscript simulate_genealogies_rcoal_revised.R
