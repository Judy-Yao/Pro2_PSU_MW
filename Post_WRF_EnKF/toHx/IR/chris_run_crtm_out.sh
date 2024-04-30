#!/bin/bash --login
#####header for stampede######
#SBATCH -J out
#SBATCH -N 4
#SBATCH -n 320
#SBATCH -p icx-normal
#SBATCH -t 2:00:00
#SBATCH -o output_crtm_out
#SBATCH -e error_crtm_out
#SBATCH --mail-user=cxh416@psu.edu
#SBATCH --mail-type=all    # Send email at begin and end of job

. util.sh

############ User control parameters
max_num_of_crtm=4   # Max number of CRTM.exe to run concurrently 
                    # (make same as # of nodes requested)
cores_per_crtm=48   # Number of cores given to each crtm.exe 
                    # (make same as # of cores per node)
date_st=201808291200 		# Start date  
date_ed=201808292300		# End date (24 forecast hrs can be done in < 2 hr w/4 nodes on skx queues)
time_int=60			# Time interval btwn cycles in minutes
nE=60				# Number of ens members
dom=1                           # Domain you are running it on
state=output                     # Input or output
fcdir=$SCRATCH/FLORENCE/gts_no_amv_update_w_blend100/fc 	# FC directory of experiment
outdir=$SCRATCH/FLORENCE/gts_no_amv_update_w_blend100/radiance                # Directory for output

##### Initialize counting variable to keep track of number of active CRTM.exe 
num_crtm=0

# Initialize date variable
date_now=$date_st

# Iterate thru dates
while [[ $date_now -le $date_ed ]]; do

  # Iterate thru ens
  for mem in `seq -f "%03g" 1 $nE`; do

    # Figure out the wrf file name
    wrffile=$fcdir/$date_now/wrf_enkf_"$state"_d0"$dom"_$mem

    # Figure out output file name
    year=${date_now:0:4} 
    month=${date_now:4:2}
    day=${date_now:6:2}
    hour=${date_now:8:2}
    minute=${date_now:10:2}
    outfile=$outdir/Radiance_"$state"_en"$mem"_d0"$dom"_"$year"-"$month"-"$day"_"$hour":"$minute".bin

    echo $date_now $mem

    # Run the CRTM.exe
    ibrun -n $cores_per_crtm -o $(($num_crtm*$cores_per_crtm)) crtm.exe $wrffile $outfile >& log.crtm_"$state"_"$date_now"_"$mem" &

    # Increment number of active CRTM programs
    num_crtm=$(($num_crtm + 1))

    # If already has max_num_of_crtm CRTM programs active, wait for all processes to clear
    if [[ $num_crtm -ge $max_num_of_crtm ]]; then
      wait
      num_crtm=0
    fi

  done #----------- End of loop over members
  
  # Increment date
  date_now=`advance_time $date_now $time_int`

done # -------- End of loop over dates
wait
