#!/bin/bash
#SBATCH --job-name="IR"
#SBATCH --partition=development
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH -t 01:00:00
#SBATCH -o out_IR
#SBATCH -e error_IR
#SBATCH --mail-user=yao.zhu.91@gmail.com
#SBATCH --mail-type=ALL

module restore intel17
source ~/.bashrc 
alias matlab19a='/work2/06191/tg854905/stampede2/opt/install_Matlab/R2019a/bin/matlab'
matlab19a  -batch main_IR_vis > log_jose
