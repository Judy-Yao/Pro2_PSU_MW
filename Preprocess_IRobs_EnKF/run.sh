#!/bin/bash
#SBATCH --job-name="validation"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --account=pen116
#SBATCH --export=ALL
#SBATCH -t 00:10:00
#SBATCH --mail-user=yao.zhu.91@gmail.com
#SBATCH --mail-type=ALL


matlab -batch vis_validation_SmallArea > validation.log

 
