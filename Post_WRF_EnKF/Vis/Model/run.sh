#!/bin/bash

#####header for stampede######
#SBATCH -J ws_dt 
#SBATCH -N 1
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 00:10:00
#SBATCH -o out_UV10
#SBATCH -e error_UV10

source ~/.bashrc
# run python script
srun python Compare_model_USAF.py &> log_model_compare 
