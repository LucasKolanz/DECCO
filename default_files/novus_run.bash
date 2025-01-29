#!/bin/bash
#SBATCH -J analyze
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 1

# export OMP_NUM_THREADS=1
# module load gnu12/12.3.0
# module load openmpi4/4.1.6
python3 ../gen_data.py
