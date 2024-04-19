#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:30:00
#SBATCH -J computeSqMat
#SBATCH -N 1
#SBATCH -G 1
export SLURM_CPU_BIND="cores"
dcgmi profile --pause
srun -n 1 -c 2 --cpu-bind=cores numactl --interleave=all ncu -o prof --set full ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/computeSqMat2/node_1/ 2>sim_err.log 1>sim_out.log
