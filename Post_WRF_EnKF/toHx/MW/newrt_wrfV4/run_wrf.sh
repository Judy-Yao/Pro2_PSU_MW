#!/bin/bash --login

#####header for stampede######
#SBATCH -J harvey_robert
#SBATCH -N 4
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 2:00:00
#SBATCH -o out_enkf
#SBATCH -e error_enkf

ibrun -n 48 ./test_wrf.exe  ../test/test_crtm_wrf.nml
