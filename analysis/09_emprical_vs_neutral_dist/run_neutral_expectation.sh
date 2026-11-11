#!/bin/bash
#SBATCH --job-name=ag_neutral
#SBATCH --partition=defq
#SBATCH --array=1-79
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=09:00:00
#SBATCH --account=m.vos
#SBATCH --output=logs/ag_neutral_%A_%a.out
#SBATCH --error=logs/ag_neutral_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cc1376@exeter.ac.uk

WORKDIR=/nobackup/beegfs/workspace/cc1376/neutral_expectation
cd "$WORKDIR" || exit 1
export OUTDIR_DAT="$WORKDIR"
export THETA_HAT=0.031
export GENES_PER_TASK=10

# Load R
module load R/4.4.0

Rscript neutral_expectation_hpc.R
