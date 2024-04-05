#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 10:00:00
#SBATCH -J fullCompSqMa
#SBATCH -N 1
#SBATCH -G 1
export SLURM_CPU_BIND="cores"
srun -n 1 -c 2 --cpu-bind=cores numactl --interleave=all ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/fullCompSqMa1/node_1/ 2>sim_err.log 1>sim_out.log
