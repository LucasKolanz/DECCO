#!/bin/bash
#SBATCH -J gen_data
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --account=lazzati
#SBATCH --partition=lazzati.q
#SBATCH --mem=16G

srun python3 gen_data.py >> quick_output.txt