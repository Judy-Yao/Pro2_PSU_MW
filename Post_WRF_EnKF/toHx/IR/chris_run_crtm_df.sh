#!/bin/bash --login
#####header for stampede######
#SBATCH -J df
#SBATCH -N 4
#SBATCH -n 320
#SBATCH -p icx-normal
#SBATCH -t 6:00:00
#SBATCH -o output_crtm
#SBATCH -e error_crtm
#SBATCH --mail-user=cxh416@psu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

. util.sh

############ User control parameters ##################################
export max_num_of_crtm=4   # Max number of CRTM.exe to run concurrently 
                           # (make same as # of nodes requested)
export cores_per_crtm=80   # Number of cores given to each crtm.exe 
                           # (make same as # of cores per node)
export storm_name=FLORENCE
export exp_names="gts_no_amv_update_w_blend100 gts_hourly_update_w_blend100 gts+ir5_hourly_no_bc_update_w_blend100 gts+ir5_hourly_no_bc_update_w_blend100_blended_with_gts+amv"
export analyses="201808291200 201808291300 201808291400 201808291500 201808291600 201808291700 201808291800 201808291900 201808292000 201808292100 201808292200 201808292300"

export forecast_length=120   # Forecast length in # of hours
export time_int=60			# Time interval btwn cycles in minutes
export dom=1                           # Domain you are running it on

########################################################################

##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

for exp_name in $exp_names
do
  export exp=$exp_name
  for analysis_time in $analyses
  do
    export analysis=$analysis_time
    if [[ $exp != "gts+ir5_hourly_no_bc_update_w_blend100_blended_with_gts+amv" ]]; then
      export rundir=$SCRATCH/"$storm_name"/"$exp"/run/$analysis/deterministic_forecast 	# Directory of experiment
    elif [[ $exp == "gts+ir5_hourly_no_bc_update_w_blend100_blended_with_gts+amv" ]]; then
      export rundir=$SCRATCH/"$storm_name"/gts+ir5_hourly_no_bc_update_w_blend100/run/$analysis/deterministic_forecast_relaxed_to_gts+amv 	# Directory of experiment
    fi 
    export outdir=$rundir/radiance                # Directory for output

    # Initialize date variable
    export date_now=$analysis
    export date_ed=`advance_time $date_now $(($time_int*$forecast_length))`
    
    # Iterate thru dates
    while [[ $date_now -le $date_ed ]]; do
    
      # Figure out input and output file names
      year=${date_now:0:4} 
      month=${date_now:4:2}
      day=${date_now:6:2}
      hour=${date_now:8:2}
      minute=${date_now:10:2}
    
      wrffile=$rundir/wrfout_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute":00
      outfile=$outdir/Radiance_df_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute".bin
    
      echo $exp $analysis $date_now $date_ed $num_crtm
      echo $wrffile
      echo $outfile

      # Run the CRTM.exe
      ibrun -n $cores_per_crtm -o $(($num_crtm*$cores_per_crtm)) crtm.exe $wrffile $outfile >& log.crtm_"$date_now" &
    
      # Increment number of active CRTM programs
      num_crtm=$(($num_crtm + 1))
    
      # If already has max_num_of_crtm CRTM programs active, wait for all processes to clear
      if [[ $num_crtm -ge $max_num_of_crtm ]]; then
        wait
        num_crtm=0
      fi
    
      # Increment date
      date_now=`advance_time $date_now $time_int`
    
    done # -------- End of loop over dates
  done  # --------- End of loop over analyses
done  # ----------- End of loop over experiments
wait
