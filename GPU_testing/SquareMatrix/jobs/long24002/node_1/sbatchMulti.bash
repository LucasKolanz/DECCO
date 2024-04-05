#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 5:00:00
#SBATCH -J long2400
#SBATCH -N 1
#SBATCH -G 1
export SLURM_CPU_BIND="cores"
srun -n 1 -c 64 --cpu-bind=cores numactl --interleave=all nsys profile -o prof ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/long24002/node_1/ 2>sim_err.log 1>sim_out.log