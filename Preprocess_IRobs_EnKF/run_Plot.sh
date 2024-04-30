#!/bin/bash
#SBATCH --job-name="IR Plot"
#SBATCH --partition=skx-dev
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --export=ALL
#SBATCH -t 01:00:00
#SBATCH -o out_plot
#SBATCH -e error_plot
#SBATCH --mail-user=yao.zhu.91@gmail.com
#SBATCH --mail-type=ALL

#alias matlab19a='/work2/06191/tg854905/stampede2/opt/install_Matlab/R2019a/bin/matlab'
/work2/06191/tg854905/stampede2/opt/install_Matlab/R2019a/bin/matlab -batch vis_toEnKFobs_IR > log_plot
