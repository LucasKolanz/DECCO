#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:05:00
#SBATCH -J newGPUtests
#SBATCH -N 1
#SBATCH -G 1
export SLURM_CPU_BIND="cores"
srun -n 1 -c 4 --cpu-bind=cores numactl --interleave=all ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/DECCOmain/SpaceLab/GPU_testing/SquareMatrix/jobs/newGPUtests1/ 2>sim_err.log 1>sim_out.log