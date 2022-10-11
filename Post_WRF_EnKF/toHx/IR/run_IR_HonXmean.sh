#!/bin/bash
# Script to calculate MW Tb from model state post-EnKF

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 3
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 00:30:00
#SBATCH -o out_ir
#SBATCH -e error_ir

. util.sh

# Fill in the storm name and experiment name
Storm=MARIA
Exper=newWRF_IR_only

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/IR

############ User control parameters
max_num_of_crtm=2   # Max number of CRTM.exe to run concurrently 
                    # (make same as # of nodes requested)
cores_per_crtm=48   # Number of cores given to each crtm.exe 
                    # (make same as # of cores per node)
date_st=201709171700        # Start date  
date_ed=201709171800        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues)
time_int=60         # Time interval btwn cycles in minutes
nE=1               # Number of ens members

# DA times 
#Obs_files_str=radiance_d03_201709031100_so
#Obs_files_str=$(ls ${Small_dir}/${Storm}/Obs_y/IR/)
#Obs_files_arr=($Obs_files_str) # String-to-array conversion gives you access to individual element
#DAtimes=()
#for obs_file in ${Obs_files_arr[@]}; do
#  DAtimes+=($(echo $obs_file | cut -c14-25))
#done


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
  # if DAtime dir exists but fort.80071/90071 hasn't been renamed
  if [[ -f ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/fort.80071 ]]; then
    mv ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/fort.80071 ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_input_d03_mean
    mv ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/fort.90071 ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_output_d03_mean
  fi
  # if DAtime dir exists but fort.80071/90071 has been renamed
  if [[ -f ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_input_d03_mean ]]; then

    outdir=${Big_dir}/${Storm}/${Exper}/Obs_Hx/IR/${DAtime}
    if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

    # Loop over each model file
    for xfile in $(ls ${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf*); do
     
      #ln -s ${Code_dir}/XbtoIR_crtm.f90.exe .
      onlyfile=$(basename "$xfile")
      outfile=${outdir}/${onlyfile}_${DAtime}_tb_g16_crtm.bin

      echo $xfile $DAtime $outfile
      
      # Run the CRTM.exe
      #ibrun -n 96 XbtoIR_crtm.exe $xfile $outfile >& XbtoIR.log &
      ibrun -n $cores_per_crtm -o $(($num_crtm*$cores_per_crtm)) XbtoIR_crtm.exe $xfile $outfile >& log_XbtoIR &
      #ibrun -n 144 XbtoIR_crtm.exe $xfile $outfile >& XbtoIR.log &
      # Increment number of active CRTM programs
      num_crtm=$(($num_crtm + 1))

      # If already has max_num_of_crtm CRTM programs active, wait for all processes to clear
      if [[ $num_crtm -ge $max_num_of_crtm ]]; then
        wait
        num_crtm=0
      fi

    done
  fi

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done
wait











