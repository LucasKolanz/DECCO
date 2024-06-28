#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:30:00
#SBATCH -J strongScaling
#SBATCH -N 2
#SBATCH -G 2
export SLURM_CPU_BIND="cores"
srun -n 2 -c 64 --cpu-bind=cores numactl --interleave=all ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/strongScaling2/node_2/ 2>sim_err.log 1>sim_out.log