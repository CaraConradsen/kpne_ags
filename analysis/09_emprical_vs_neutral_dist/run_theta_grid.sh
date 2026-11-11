#!/bin/bash
#SBATCH --job-name=ag_theta_grid
#SBATCH --partition=defq
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=09:00:00
#SBATCH --account=m.vos
#SBATCH --output=logs/ag_theta_grid_%A_%a.out
#SBATCH --error=logs/ag_theta_grid_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cc1376@exeter.ac.uk

WORKDIR=/nobackup/beegfs/workspace/cc1376/neutral_expectation
cd "$WORKDIR" || exit 1
export OUTDIR_DAT="$WORKDIR"

# one array task per theta; keep --array in step with the length of this list
export THETA_GRID=0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10

# median synonymous sites over the syntenic AGs
export M_TARGET=400

mkdir -p logs "$WORKDIR/theta_grid_out"

# Load R
module load R/4.4.0

Rscript neutral_expectation_hpc_theta_grid.R
