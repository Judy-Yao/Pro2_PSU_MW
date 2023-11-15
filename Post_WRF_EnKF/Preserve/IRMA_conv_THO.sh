#!/bin/bash
# Script to create a compressed archive

#####header for stampede######
#SBATCH -J conv_THO
#SBATCH -N 1
#SBATCH --ntasks-per-node 48
#SBATCH -p icx-normal
#SBATCH -t 3:00:00
#SBATCH -o out_fc_IRMA_conv_THO
#SBATCH -e error_fc_IRMA_conv_THO
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zuy121@psu.edu

function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}
export -f advance_time

# Storm
Storm=IRMA
# Experiment
DA=conv
MP=THO
Exper=J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900
# Specify time range
spinup_st=201709021200
date_st=201709021200
date_ed=201709021200

# Navigate to the directory where the file is located
cd /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}/fc
echo 'Enter' /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}/fc/

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do
  if [[ ! -d ${DAtime} ]]; then
    continue
  fi
  # Enter the time
  cd $DAtime
  # Special treatment to spin up start time
  if [[ ${DAtime} == ${spinup_st} ]]; then
    echo 'Dealing with spinup start'
    # Run the tar command to create the compressed archive
    # wrfbdy_d01*
    bdy_files=($(find . -maxdepth 1 -type f -name 'wrfbdy*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfbdy_d01.tar.gz "${bdy_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfbdy* takes '$runtime's'

    # wrfinput_d01_0*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d01_0*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d01_0.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d01_0* takes '$runtime's'
    # wrfinput_d02_0*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d02_0*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d02_0.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d02_0* takes '$runtime's'
    # wrfinput_d03_0*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d03_0*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d03_0.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d03_0* takes '$runtime's'

    # wrfinput_d01_2017*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d01_2017*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d01_2017.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d01_2017* takes '$runtime's'
    # wrfinput_d02_2017*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d02_2017*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d02_2017.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d02_2017* takes '$runtime's'
    # wrfinput_d03_2017*
    wrf_files=($(find . -maxdepth 1 -type f -name 'wrfinput_d03_2017*'))
    start=`date +%s.%N`
    tar -czvf ${DAtime}_wrfinput_d03_2017.tar.gz "${wrf_files[@]}"
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'wrfinput_d03_2017* takes '$runtime's'

  fi
  # Increment date
  DAtime=`advance_time $DAtime $time_int`

  # Return to the parent directory
  cd ..
done

# Check the exit status
echo $?

#Note: with test, running tar only takes longer time than tar+compression
# For 9 files
#Only tar takes 32.047625489s
#Tar with compression takes 9.842801481s
