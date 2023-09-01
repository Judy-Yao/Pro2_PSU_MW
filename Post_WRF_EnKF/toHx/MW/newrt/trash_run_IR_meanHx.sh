#!/bin/bash
# Script to calculate IR Tb from model state post-EnKF
# Run the CRTM on every member

# Credit to Yinghui LV and Yunji Zhang
# Author: Zhu (Judy) Yao. July 27 - 28, 2022

#####header for stampede######
#SBATCH -J IR
#SBATCH -N 3
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 02:00:00
#SBATCH -o out_mw
#SBATCH -e error_mw
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com

module purge
module load intel/18.0.2
module load impi/18.0.2
module load parallel-netcdf/4.3.3.1
module load phdf5/1.8.16
module load libfabric/1.7.0
module load python3/3.7.0 

# Fill in the storm name and experiment name
Storm=HARVEY
Exper=J_DA+J_WRF+J_init

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/Post_WRF_EnKF/toHx/IRorMW/newrt #/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/MW/newrt

############ User control parameters
date_st=201708221200        # Start date  
date_ed=201708221200        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues
nE=60               # Number of ens members
dom=3                           # Domain you are running it on 
state=("input" "output") # Input or output or both
sensor='abi_gr'

# ------ DA times where MW obs exists ----------------------
# ----------------------------------------------------------
Obs_files_str=$(ls ${Small_dir}/${Storm}/Obs_y/IR/)
Obs_files_arr=($Obs_files_str) # String-to-array conversion gives you access to individual element
DAtimes=()
for obs_file in ${Obs_files_arr[@]}; do
  DAtimes+=($(echo $obs_file | cut -c14-25))
done

# Iterate thru DAtimes
for DAtime in ${DAtimes[@]}; do
  
  echo $DAtime

  # Determine if DAtime is of interest to this calculation
  if [[ $DAtime -lt $date_st ]] || [[ $DAtime -gt $date_ed ]]; then
    continue
  fi

  # if the fc directory at DAtime does not exist, terminate the loop
  if [[ ! -d ${Big_dir}/${Storm}/${Exper}/fc/${DAtime} ]]; then
    break
  fi
 
  # check if output dir exists
  outdir=${Big_dir}/${Storm}/${Exper}/Obs_Hx/IR/${DAtime}/
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

  cd ${Big_dir}/${Storm}/${Exper}/Obs_Hx/IR/$DAtime
  ln -sf ${Code_dir}/paper_scott_2020.mpi .
  ln -sf ${Code_dir}/getSensorInfo.py .
  ln -sf /work2/06191/tg854905/stampede2/opt/CRTM/PSU_EnKF_CRTM/coefficients coefficients
 
  # Iterate thru states
  for istate in "${state[@]}"; do

    # Iterate thru ens
    for mem in `seq -f "%03g" 1 $nE`; do
      xfile=${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_"${istate}"_d0"$dom"_$mem 
      #onlyfile=$(basename "$xfile") 
      #outfile=${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/${onlyfile} 
      outfile=${outdir}/"${istate}"_mem"$mem"_d0"$dom"_${DAtime}

      #echo ${xfile}
      ln -sf ${xfile} wrffile

      cat > test_crtm_wrf.nml << EOF
&rt_settings
nml_s_rt_program = 'crtm'
nml_l_include_land = T,
nml_s_sensor_id = ${sensor},
nml_a_channels(:,1) = 8,
nml_s_reff_method = 'mp_physics',
nml_i_nicpu = 12,
nml_i_njcpu = 12,
nml_s_crtm_rainLUT='Thompson08_RainLUT_-109z-1.bin',
nml_s_crtm_snowLUT='Thompson08_SnowLUT_-109z-1.bin',
nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1.bin',
/
    
\$rt_input
nml_s_filename_input = '${xfile}'
nml_s_filename_obs='${Small_dir}/${Storm}/Obs_y/IR/radiance_d03_${DAtime}_so' 
/ 
            
\$rt_output
nml_s_filename_output = '${outfile}.tb'
nml_l_output_reff = F,
/
EOF

      # Run the CRTM
      if [[ ! -f ${outfile}.tb.${sensor_Chs[0]}.crtm.nc ]]; then
        ibrun -n 144 ./paper_scott_2020.mpi test_crtm_wrf.nml
      fi
    
    done # End looping over ens
  done #  End looping over states
done # End looping over dates          



#wait


