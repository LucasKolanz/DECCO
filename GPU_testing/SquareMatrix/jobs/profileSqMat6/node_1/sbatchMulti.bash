#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 1:30:00
#SBATCH -J profileSqMat
#SBATCH -N 1
#SBATCH -G 1
export SLURM_CPU_BIND="cores"
srun -n 1 -c 2 --cpu-bind=cores numactl --interleave=all nsys profile -o SqMaprof6 ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/profileSqMat6/node_1/ 2>sim_err.log 1>sim_out.log
