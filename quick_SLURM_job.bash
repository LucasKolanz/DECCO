#!/bin/bash
#SBATCH -J check_BAPA
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1

srun python3 -u verify_aggregation.py >> tmp/BAPA_out.txt