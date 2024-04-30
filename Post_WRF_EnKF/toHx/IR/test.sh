#!/bin/bash

#SBATCH -J myjob           # Job name
#SBATCH -o myjob.o%j       # Name of stdout output file
#SBATCH -e myjob.e%j       # Name of stderr error file
#SBATCH -p skx-dev      # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
#SBATCH -n 32              # Total # of mpi tasks
#SBATCH -t 01:30:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=zuy121@psu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
##SBATCH -A myproject       # Allocation name (req'd if you have more than 1)

echo 'what?'
