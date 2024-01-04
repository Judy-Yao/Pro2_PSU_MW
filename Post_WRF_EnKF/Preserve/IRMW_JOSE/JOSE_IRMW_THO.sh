#!/bin/bash
# Script to create a compressed archive

#####header for stampede######
#SBATCH -J JOSEirmwTHO
#SBATCH -N 1
#SBATCH --ntasks-per-node 48
#SBATCH -p icx-normal
#SBATCH -t 20:00:00
#SBATCH -o out_JOSE_irmw_THO
#SBATCH -e error_JOSE_irmw_THO
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

function clean_icbc {
 mkdir dir_log
 mv *.0000 dir_log
 mv *.log dir_log
 rm *.nc
 rm *.log*
 rm rsl*
 rm FILE*
 rm wrfbdy_d01 wrfinput_d01 wrfinput_d02 wrfinput_d03
 mv dir_log/* .
}
export -f clean_icbc

function clean_perturb_ic {
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    rm  ${id}/wrfinput_d01 ${id}/wrfinput_d02 ${id}/wrfinput_d03
  done
}
export -f clean_perturb_ic

function clean_obsproc {
  rm obs.raw
}
export -f clean_obsproc

function clean_enkf {
  for file in $(ls */fort.800*); do unlink $file; done
  for file in $(ls */fort.900*); do unlink $file; done
  for file in $(ls */fort.600*); do rm $file; done
}
export -f clean_enkf
 
function clean_forecast {
  mv rsl.error.0000  lala
  rm rsl*
  rm wrfout_d01*
  rm wrfout_d02_*
  mv lala rsl.error.0000
}
export -f clean_forecast


# Storm
Storm=JOSE
# Experiment
DA=IR+MW
MP=THO
Exper=IR+MW-J_DA+J_WRF+J_init-SP-intel17-THO-30hr-hroi900
NUM_ENS=60
# Specify time range
spinup_st=201709041200
date_st=201709041200
date_ed=201709070000
# Specify which dir to process
at_run=true
at_wrf_df=true

# run directory
if $at_run; then
  # Navigate to the directory where the file is located
  cd /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}/run/
  echo 'Enter' /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}/run/
  # Initialize date variable
  DAtime=$date_st
  # Iterate thru dates
  while [[ $DAtime -le $date_ed ]]; do
    echo $DAtime
    # skip the time if it does not exist
    if [[ ! -d ${DAtime} ]]; then
      DAtime=`advance_time $DAtime 60`
      continue  
    else
      # Enter the time directory
      cd $DAtime; echo 'Entering /run/' $DAtime 
    fi
    # LBDATE
    export minute_off=$(echo "(${DAtime:8:2}*60+${DAtime:10:2})%360" |bc)
    if [[ $DAtime == $(advance_time "$DAtime" -"$minute_off") ]]; then
      LBDATE=$(advance_time "$DAtime" -"$minute_off")
    fi
    # Special treatment to spin up start time
    if [[ ${DAtime} == ${spinup_st} ]]; then
      echo 'Dealing with spinup start'
      # icbc
      cd icbc; echo 'clean icbc...'
      $(clean_icbc)
      cd ..
      # perturb_ic
      cd perturb_ic; echo 'clean perturb_ic...'
      $(clean_perturb_ic)
      cd ..
    # Synoptic times 
    elif [ ${DAtime} == ${LBDATE} ]; then 
      # icbc
      cd icbc; echo 'clean icbc...'
      $(clean_icbc)
      cd ..
      # perturb_ic
      cd perturb_ic; echo 'clean perturb_ic...'
      $(clean_perturb_ic)
      cd .. 
      # obsproc
      cd obsproc; echo 'clean obsproc...'
      $(clean_obsproc)
      cd ..
      # enkf
      cd enkf; echo 'clean enkf...'
      $(clean_enkf)
      cd ..
    # Other regular times
    else
      # obsproc
      cd obsproc; echo 'clean obsproc...'
      $(clean_obsproc)
      cd ..
      # enkf
      cd enkf; echo 'clean enkf...'
      $(clean_enkf)
      cd ..
    fi
  
    # Increment date
    DAtime=`advance_time $DAtime 60`
    # Return to the parent directory
    cd ..

  done
  # return to the main directory
  cd /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}
  # tar
  start=`date +%s.%N`
  tar -czvf ${Storm}_${DA}_${MP}_run.tar.gz run 
  end=`date +%s.%N`
  runtime=$( echo "$end - $start" | bc -l )
  echo '/run takes '$runtime's' 
  # Check the exit status
  echo $?
fi

# wrf_df
if ${at_wrf_df}; then
  # main directory
  cd /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}
  # list initialization times
  times=$(ls wrf_df)
  # enter wrf_df directory
  cd wrf_df; echo 'Enter wrf_df'
  # compress each forecast
  for it in ${times}; do
    # clean
    cd ${it}; echo 'clean '${it}
    $(clean_forecast)
    cd ..
    # compress
    start=`date +%s.%N`
    tar -czvf ${it}.tar.gz ${it}
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo '/run takes '$runtime's'
    echo $?
  done
fi


#Note: with test, running tar only takes longer time than tar+compression
# For 9 files
#Only tar takes 32.047625489s
#Tar with compression takes 9.842801481s
