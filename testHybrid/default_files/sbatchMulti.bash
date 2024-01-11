#!/bin/bash
#SBATCH -A m4189
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 3:00:00
#SBATCH -n 1
#SBATCH -c 256

export OMP_NUM_THREADS=128
export SLURM_CPU_BIND="cores"
srun ./ColliderMultiCore.x /global/u2/l/lpkolanz/SpaceLab/jobs/smallerDt2/ 30 2>sim_err.log 1>sim_out.log