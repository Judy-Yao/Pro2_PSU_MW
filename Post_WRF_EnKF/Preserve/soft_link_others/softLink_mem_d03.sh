#! /bin/bash

function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}
export -f advance_time


path_from=/scratch/02191/yuz31/HARVEY/IR_WSM6/fc/
path_to=/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_WSM6/fc/

date_st=201708221200
date_ed=201708241200
time_int=60         # Time interval btwn cycles in minutes
nE=60
dom=3                           # Domain you are running it on 

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do
   
  cd $path_from/${DAtime}
  echo 'Enter ' ${path_from}/${DAtime}

  # figure out the time for this ensemble
  year=${DAtime:0:4}
  month=${DAtime:4:2}
  day=${DAtime:6:2}
  hour=${DAtime:8:2}
  minute=${DAtime:10:2}

  outdir=${path_to}/${DAtime}/
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi 
  for mem in `seq -f "%03g" 1 $nE`; do  
    ln -s  ${path_from}/${DAtime}/wrf_enkf_input_d03_${mem}  ${path_to}/${DAtime}/wrf_enkf_input_d03_${mem}
    ln -s  ${path_from}/${DAtime}/wrf_enkf_output_d03_${mem}  ${path_to}/${DAtime}/wrf_enkf_output_d03_${mem}
  done

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done
