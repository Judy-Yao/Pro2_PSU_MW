#!/bin/bash
# Script to calculate MW Tb from model state post-EnKF
# Run the CRTM on every member

# Credit to Yinghui LV and Yunji Zhang
# Author: Zhu (Judy) Yao. July 27 - 28, 2022

#####header for stampede######
#SBATCH -J MW
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
Exper=J_DA+Y_WRF+J_init-IR+MW
Exper_obs=J_DA+Y_WRF+J_init-IR+MW


# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/Post_WRF_EnKF/toHx/MW/newrt_wrfV4 #/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/MW/newrt

############ User control parameters
date_st=201708221200        # Start date  
date_ed=201708221300        # End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues
nE=60               # Number of ens members
dom=3                           # Domain you are running it on 
state=("input" "output") # Input or output or both

# ------ DA times where MW obs exists ----------------------
# ----------------------------------------------------------
#Obs_files_str=microwave_d03_201708221200_so
#Obs_files_str=$(ls ${Small_dir}/${Storm}/Obs_y/MW/Processed_2nd_time/$Exper_obs)
Obs_files_str=$(ls ${Small_dir}/${Storm}/Obs_y/MW/Processed_1st_time)
Obs_files_arr=($Obs_files_str) # String-to-array conversion gives you access to individual element
DAtimes=()
for obs_file in ${Obs_files_arr[@]}; do
  DAtimes+=($(echo $obs_file | cut -c15-26))
#  DAtimes+=($(echo $obs_file | cut -c11-22))
done

# Iterate thru DAtimes
for DAtime in ${DAtimes[@]}; do
  
  echo $DAtime
  #if ! ls -l nofolder; then
  #  echo "Folder doesn't exist"
  #  exit 1
  #fi

  # Determine if DAtime is of interest to this calculation
  if [[ $DAtime -lt $date_st ]] || [[ $DAtime -gt $date_ed ]]; then
    continue
  fi

  # if the fc directory at DAtime does not exist, terminate the loop
  if [[ ! -d ${Big_dir}/${Storm}/${Exper}/fc/${DAtime} ]]; then
    break
  fi
 
  # check if output dir exists
  outdir=${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

  # At the DAtime, get sensor information based on microwave_SO (obs)
  if [[ ! -f ${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/${DAtime}_sensorCh ]]; then
    $(python3 getSensorInfo.py ${Storm} ${Exper} ${Exper_obs} ${DAtime})
  fi
  Sensor_Info=${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/${DAtime}_sensorCh 


  cd ${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/$DAtime
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

      # Check how many unique sensor+channel are available at this DA time
      # Num_Ch=$(wc -l < ${Sensor_Info})

      # Read sensor and channel information from file ${DAtime}_sensorCh
      while IFS= read -r line
        do
          sensor_Chs=($line)

          cat > test_crtm_wrf.nml << EOF
&rt_settings
nml_s_rt_program = 'crtm'
nml_l_include_land = T,
nml_s_sensor_id = ${sensor_Chs[0]},
nml_a_channels(:,1) = ${sensor_Chs[@]:1},
nml_s_reff_method = 'mp_physics',
nml_i_nicpu = 12,
nml_i_njcpu = 12,
nml_s_crtm_rainLUT='Thompson08_RainLUT_-109z-1.bin',
nml_s_crtm_snowLUT='Thompson08_SnowLUT_-109z-1.bin',
nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1.bin',
/
    
\$rt_input
nml_s_filename_input = '${xfile}'
nml_s_filename_obs='${Small_dir}/${Storm}/Obs_y/MW/Processed_1st_time/microwave_d03_${DAtime}_so'
!nml_s_filename_obs='${Small_dir}/${Storm}/Obs_y/MW/Processed_2nd_time/${Exper_obs}/microwave_${DAtime}_so' 
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
      
      done < "${Sensor_Info}" # End looping over overpass sensors
    done # End looping over ens
  done #  End looping over states
done # End looping over dates
            

#cp run_MW_meanHx.sh ${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/   


#wait


