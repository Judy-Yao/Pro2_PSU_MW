#!/bin/bash -x
#SBATCH -J IR_WSM6
#SBATCH -p icx
#SBATCH -n 48 -N 1
#SBATCH -t 0:20:00
#SBATCH -o HARVEY.batch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com
date
module restore intel
python3 IR_WSM6.py > IR_WSM6.log 
date
