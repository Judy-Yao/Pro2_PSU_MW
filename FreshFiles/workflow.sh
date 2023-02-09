#!/bin/bash
# Script to refresh the time stamp of each files under a certain directory to prevent it being purged 
# Author: Zhu (Judy) Yao. Jan 2, 2023

#####header for stampede######
#SBATCH -J Fresh
#SBATCH -N 1
#SBATCH --ntasks-per-node 48
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -o out_fresh
#SBATCH -e error_fresh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com


# Fill in the storm name and experiment name
Storm=HARVEY
Exp_interest=("J_DA+J_WRF+J_init")
#Exp_interest=("J_DA+J_WRF+J_init" "IR-J_DA+J_WRF+J_init" "J_DA+Y_WRF+J_init-IR+MW" "J_DA+Y_WRF+J_init-IR" "JerryRun") #IR+MW-J_DA+J_WRF+J_init-SP-intel19

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/

# Iterate thru Experiments
for Exper in ${Exp_interest[@]}; do
  echo $Exper
  $(python3 operation.py ${Storm} ${Big_dir} ${Exper})

done # End looping over Experiments

# At this point (Jan 2, 2023), Judy hasn't figured out how to redirect output from python program to a log file
