#!/bin/bash
# Script to calculate MW Tb from model state post-EnKF
# Run the CRTM on every member

# Credit to Yinghui LV and Yunji Zhang
# Author: Zhu (Judy) Yao. July 27 - 28, 2022

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 4
#SBATCH --ntasks-per-node 48
#SBATCH -p development
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
Exper=newWRF_IR_only

# Parent paths
Big_dir=/scratch/06191/tg854905/Pro2_PSU_MW/
Small_dir=/work2/06191/tg854905/stampede2/Pro2_PSU_MW
Code_dir=/home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/MW/newrt

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
state=("input" "output") # Input or output or both

##### Initialize counting variable to keep track of number of active CRTM.exe
num_crtm=0 

# ------ DA times where MW obs exists ---------------
Obs_files_str=microwave_d03_201708221200_so
#Obs_files_str=$(ls ${Small_dir}/${Storm}/Obs_y/MW/)
Obs_files_arr=($Obs_files_str) # String-to-array conversion gives you access to individual element
DAtimes=()
for obs_file in ${Obs_files_arr[@]}; do
  DAtimes+=($(echo $obs_file | cut -c15-26))
done

# Iterate thru DAtimes
for DAtime in ${DAtimes[@]}; do

  # if the fc directory at DAtime does not exist, terminate the loop
  if [[ ! -d ${Big_dir}/${Storm}/${Exper}/fc/${DAtime} ]]; then
    break
  fi
 
  # check if output dir exists
  outdir=${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi

  # At the DAtime, get sensor information based on microwave_SO (obs)
  if [[ ! -f ${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/${DAtime}_sensorCh ]]; then
    $(python3 getSensorInfo.py ${Storm} ${Exper} ${DAtime})
  fi
  Sensor_Info=${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/${DAtime}/${DAtime}_sensorCh 


  # figure out the time for this ensemble
  year=${DAtime:0:4}
  month=${DAtime:4:2}
  day=${DAtime:6:2}
  hour=${DAtime:8:2} 
  minute=${DAtime:10:2}

  cd ${Big_dir}/${Storm}/${Exper}/Obs_Hx/MW/$DAtime
  ln -sf ${Code_dir}/paper_scott_2020.mpi .
  ln -sf ${Code_dir}/getSensorInfo.py .
  ln -sf /work2/06191/tg854905/stampede2/opt/CRTM/PSU_EnKF_CRTM/coefficients coefficients
 
  # Iterate thru states
  for istate in "${state[@]}"; do

    # Iterate thru ens
    for mem in `seq -f "%03g" 1 $nE`; do
      xfile=${Big_dir}/${Storm}/${Exper}/fc/${DAtime}/wrf_enkf_"$state"_d0"$dom"_$mem 
      outfile=${outdir}/"$state"_mem"$mem"_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute"

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
nml_s_filename_obs='${Small_dir}/${Storm}/Obs_y/MW/microwave_d03_${DAtime}_so' 
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
            

wait


