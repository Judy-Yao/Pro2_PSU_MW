#!/bin/bash
# Script to forward calculate model state to IR Tb
# Mean of H(x)

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 4
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 02:00:00
#SBATCH -o out_ir
#SBATCH -e error_ir

. util.sh

# Fill in the storm name and experiment name
Storm=HARVEY
Exper=newWRF_MW_THO

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/IR

############ User control parameters
max_num_of_crtm=4   # Max number of CRTM.exe to run concurrently 
                    # (make same as # of nodes requested)
cores_per_crtm=48   # Number of cores given to each crtm.exe 
                    # (make same as # of cores per node)
date_st=201708221200        # Start date  
date_ed=201708231600        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues)
time_int=60         # Time interval btwn cycles in minutes
nE=60               # Number of ens members
dom=3                           # Domain you are running it on 
#state=("output")
state=("input" "output") # Input or output or both


##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

# Initialize date variable
DAtime=$date_st

# Iterate thru states
for istate in "${state[@]}"; do

  # Iterate thru dates
  while [[ $DAtime -le $date_ed ]]; do

    # if the fc directory at DAtime does not exist, terminate the loop
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

    # Iterate thru ens
    for mem in `seq -f "%03g" 1 $nE`; do
      wrffile=${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_"$istate"_d0"$dom"_$mem 
      outfile=${outdir}/TB_GOES_CRTM_"$istate"_mem"$mem"_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute".bin

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

    done # End looping over ens
  
    # Increment date
    DAtime=`advance_time $DAtime $time_int`

  done # End looping over dates

done # End looping over states

wait # Trick: the background process will leave the compute node once the process goes over the shell script. If at that moment the crtm processes are not finished yet, these crtm processes will crash.











