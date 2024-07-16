#!/bin/bash -x
#SBATCH -J HARVEY
#SBATCH -p icx
#SBATCH -n 48 -N 1
#SBATCH -t 5:00:00
#SBATCH -o HARVEY.batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com
date
module restore intel
python3 Analyses_vars_HARVEY.py > HARVEY.log 
date
