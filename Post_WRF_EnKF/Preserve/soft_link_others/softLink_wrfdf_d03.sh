#! /bin/bash

function advance_time {
  ccyymmdd=`echo $1 |cut -c1-8`
  hh=`echo $1 |cut -c9-10`
  mm=`echo $1 |cut -c11-12`
  inc=$2
  date -u -d $inc' minutes '$ccyymmdd' '$hh':'$mm +%Y%m%d%H%M
}
export -f advance_time


path_from=/scratch/02191/yuz31/HARVEY/IR_THO/output/
path_to=/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/wrf_df/

date_st=201708221800
date_ed=201708241200
time_int=360         # Time interval btwn cycles in minutes
nE=60
dom=3                           # Domain you are running it on 

# Initialize date variable
DAtime=$date_st

# Iterate thru dates
while [[ $DAtime -le $date_ed ]]; do

  cd $path_from/${DAtime}
  echo 'Enter ' ${path_from}/${DAtime}

  outdir=${path_to}/${DAtime}/
  if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi 
  ln -s  ${path_from}/${DAtime}/ATCF  ${path_to}/${DAtime}/ATCF_rsl.error.0000

  # Increment date
  DAtime=`advance_time $DAtime $time_int`

done
