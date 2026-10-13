#!/bin/bash
#SBATCH --job-name=run_msprime_sweep
#SBATCH --nodes=1
#SBATCH --partition=hmq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=90G
#SBATCH --time=09:00:00
#SBATCH --account=m.vos
#SBATCH --array=0-0
#SBATCH --output=run_msprime_sweep_%A_%a.out
#SBATCH --error=run_msprime_sweep_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=cc1376@exeter.ac.uk

WORKDIR=/nobackup/beegfs/workspace/cc1376/run_msprime
cd "$WORKDIR" || exit 1
export OUTDIR_DAT="$WORKDIR"

# --- growth rates, in coalescent units -----------------------------------
# index 0 : generic expansion   -- placeholder, choose your own value
# index 1 : calibrated          -- from calibrate_growth_rate.py --target 0.65
ALPHAS=(2.88)
ALPHA=${ALPHAS[$SLURM_ARRAY_TASK_ID]}

if [ "$ALPHA" = "999" ]; then
  echo "ERROR: calibrated alpha is still the placeholder 999." >&2
  echo "Run calibrate_growth_rate.py first, then edit ALPHAS." >&2
  exit 1
fi

# --- environment ---------------------------------------------------------
module load cc1376/msprime/0.4.0

# Fail fast if the compute node resolves the module differently to the
# login node, or if the msprime API is 0.x rather than 1.x.
python -c 'import msprime, sys; print("msprime", msprime.__version__, "|", sys.executable)'
python -c 'import msprime; assert hasattr(msprime, "sim_ancestry"), "msprime 1.x required"'

# --- run -----------------------------------------------------------------
python simulate_genealogies_msprime.py \
  --n-genomes 260 \
  --n-trees 500000 \
  --growth-rate "$ALPHA" \
  --chunk-size 10000 \
  --seed $((42 + SLURM_ARRAY_TASK_ID)) \
  --outdir "$OUTDIR_DAT"
