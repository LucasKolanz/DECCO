#!/bin/bash
#SBATCH -A m2651
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 0:30:00
#SBATCH -J strongScaling
#SBATCH -N 4
#SBATCH -G 4
export SLURM_CPU_BIND="cores"
srun -n 4 -c 64 --cpu-bind=cores numactl --interleave=all ./ColliderMultiCore.x /pscratch/sd/l/lpkolanz/SpaceLab/testSqMat/jobs/strongScaling2/node_4/ 2>sim_err.log 1>sim_out.log