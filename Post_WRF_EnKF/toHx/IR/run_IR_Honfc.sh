#!/bin/bash
# Script to forward calculate model state to IR Tb
# Mean of H(x)

#SBATCH -J IR           # Job name
#SBATCH -o IR.o%j       # Name of stdout output file
#SBATCH -e IR.e%j       # Name of stderr error file
#SBATCH -p skx-dev      # Queue (partition) name
#SBATCH -N 4               # Total # of nodes 
#SBATCH --ntasks-per-node 48
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=zuy121@psu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
##SBATCH -A myproject       # Allocation name (req'd if you have more than 1)

module restore intel
source util.sh

# Fill in the storm name and experiment name
Storm=HARVEY
#Exper=IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi300/wrf_df/201708240600/
Exper=JerryRun/IR_WSM6/fc/201708220000/

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/IR

############ User control parameters
max_num_of_crtm=4  # Max number of CRTM.exe to run concurrently 
                    # (make same as # of nodes requested)
cores_per_crtm=48   # Number of cores given to each crtm.exe 
                    # (make same as # of cores per node)
date_st=201708240900        # Start date  
date_ed=201708240900        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues)
time_int=60         # Time interval btwn cycles in minutes
dom=3                           # Domain you are running it on 

##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do

  # figure out the time for this ensemble
  year=${DAtime:0:4}
  month=${DAtime:4:2}
  day=${DAtime:6:2} 
  hour=${DAtime:8:2}
  minute=${DAtime:10:2}  

  wrffile=${Big_dir}/${Storm}/${Exper}/wrfout_d0${dom}_${year}-${month}-${day}_${hour}:${minute}:00
  outfile=${Big_dir}/${Storm}/${Exper}/CRTM_wrfout_d0${dom}_${year}-${month}-${day}_${hour}:${minute}:00.bin

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

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done # End looping over dates

wait # Trick: the background process will leave the compute node once the process goes over the shell script. If at that moment the crtm processes are not finished yet, these crtm processes will crash.











