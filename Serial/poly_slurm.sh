#!/bin/bash
#SBATCH --job-name=polynomial
#SBATCH --partition=a100
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:10 # time (D-HH:MM)
#SBATCH --output=init.stdout
#SBATCH --error=init.stderr

./run 1000 out 10 11 12 13 14 15 16
