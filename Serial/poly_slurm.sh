#!/bin/bash
#SBATCH --job-name=polynomial
#SBATCH --partition=cosc3500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0-01:00 # time (D-HH:MM)
#SBATCH --output=init.stdout
#SBATCH --error=init.stderr

./run 1000 out 10 11 12 13 14 15 16
