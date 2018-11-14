#!/bin/bash
#
#SBATCH --job-name=test
#SBATCH --output=out.txt
#
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=2048
#SBATCH --nodes=1


srun -N 1 -n 1 -c 1 python /scratch/mape1416/local/untitled0.py &
