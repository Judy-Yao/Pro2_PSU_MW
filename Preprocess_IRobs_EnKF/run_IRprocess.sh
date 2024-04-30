#!/bin/bash
#SBATCH --job-name="IR"
#SBATCH --partition=icx-normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH -t 02:00:00
#SBATCH -o out_IR
#SBATCH -e error_IR
#SBATCH --mail-user=yao.zhu.91@gmail.com
#SBATCH --mail-type=ALL

module load matlab/2022a
matlab -batch main_process_IRobs > log_run_IRprocess
