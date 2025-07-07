#!/bin/bash
# Script to calculate MW Tb from model state post-EnKF

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 1
#SBATCH --ntasks-per-node 8
#SBATCH -p skx-dev
#SBATCH -t 00:30:00
#SBATCH -o out_ir
#SBATCH -e error_ir

source util.sh
module restore intel

# Fill in the storm name and experiment name
Storm=IRMA
Exper=CONV+IR_WSM6

# Parent paths
Big_dir=/scratch/06191/tg854905/Clean_Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/IR

############ User control parameters
max_num_of_crtm=1   # Max number of CRTM.exe to run concurrently 
                    # (make same as # of nodes requested)
cores_per_crtm=8   # Number of cores given to each crtm.exe 
                    # (make same as # of cores per node)
date_st=201709030000        # Start date  
date_ed=201709030000        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues)
time_int=60         # Time interval btwn cycles in minutes
nE=1               # Number of ens members
dom=3
state=( "input" "output") # Input or output or both


##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do

  # if the DAtime directory does not exist (no wrf file of any kind), terminate the loop  
  if [[ ! -d ${Big_dir}/${Storm}/${Exper}/fc/${DAtime} ]]; then
    break
  fi
  
  # check if output dir exists
  outdir=${Big_dir}/${Storm}/${Exper}/Obs_Hx/IR/${DAtime}
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

  # figure out the time for this ensemble
  year=${DAtime:0:4}
  month=${DAtime:4:2}
  day=${DAtime:6:2}
  hour=${DAtime:8:2}
  minute=${DAtime:10:2}

  # Iterate thru states 
  for istate in "${state[@]}"; do
    echo 'Calculating ' $istate '.......'

    wrffile=${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_"$istate"_d0"$dom"_mean
    outfile=${outdir}/TB_GOES_CRTM_"$istate"_mem_mean_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute".bin

    echo $wrffile $outfile

    # Run the CRTM.exe
    ibrun -n $cores_per_crtm -o $(($num_crtm*$cores_per_crtm)) XbtoIR_crtm.exe $wrffile $outfile >& log_XbtoIR &
    # Increment number of active CRTM programs
    num_crtm=$(($num_crtm + 1))

    # If already has max_num_of_crtm CRTM programs active, wait for all processes to clear
    if [[ $num_crtm -ge $max_num_of_crtm ]]; then
      wait
      num_crtm=0
    fi

  done # End looping over states

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done
wait


