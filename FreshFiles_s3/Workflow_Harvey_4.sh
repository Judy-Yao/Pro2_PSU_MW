#!/bin/bash
# Script to refresh the time stamp of each files under a certain directory to prevent it being purged 
# Author: Zhu (Judy) Yao. Jan 2, 2023

#####header for stampede######
#SBATCH -J Fresh
#SBATCH -N 1
#SBATCH --ntasks-per-node 48
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -o out_fresh
#SBATCH -e error_fresh
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com


source util.sh

#---------------------------------------------------------
# Configuration
#---------------------------------------------------------
Storm=HARVEY
Exp_interest=("CONV+IR_THO")
#Exp_interest=("J_DA+J_WRF+J_init" "IR-J_DA+J_WRF+J_init" "J_DA+Y_WRF+J_init-IR+MW" "J_DA+Y_WRF+J_init-IR" "JerryRun") #IR+MW-J_DA+J_WRF+J_init-SP-intel19
StartDate_spinup=201708220000
StartDate_cycling=201708221200
EndDate_cycling=201708241200
NUM_ENS=60
MAX_DOM=3

if_refresh_fc=false
if_refresh_all=true

# Parent paths
Big_dir=/scratch/06191/tg854905/Clean_Pro2_PSU_MW/


#---------------------------------------------------------
# Refresh files under /fc
#---------------------------------------------------------
if $if_refresh_fc; then
  # Iterate thru Experiments
  for Exper in ${Exp_interest[@]}; do
    cd $Big_dir/$Storm/$Exper/fc
    echo 'cd' $Big_dir/$Storm/$Exper/ 'fc'
    for dir in $(ls $pwd); do 
      cd $dir
      for file in $(ls wrf*); do 
        ncdump -h $file
      done
      cd .. 
    done
    #$(python3 operation.py ${Storm} ${Big_dir} ${Exper})
  done # End looping over Experiments
fi

#---------------------------------------------------------
# Refresh other files
#---------------------------------------------------------
if $if_refresh_all; then
  # Iterate thru Experiments
  for Exper in ${Exp_interest[@]}; do
    cd $Big_dir/$Storm/$Exper/

    #---------------
    # /fc
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper'/fc'
    start=`date +%s.%N`
    find $Big_dir/$Storm/$Exper/fc -exec touch {} \;
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

    #---------------
    # /obs
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper'/obs'
    start=`date +%s.%N`
    find $Big_dir/$Storm/$Exper/obs -exec touch {} \;
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

    #---------------
    # /Obs_Hx
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper'/Obs_Hx'
    start=`date +%s.%N`
    find $Big_dir/$Storm/$Exper/Obs_Hx -exec touch {} \;
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

    #---------------
    # /rc
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper/'rc'
    start=`date +%s.%N`
    find $Big_dir/$Storm/$Exper/rc -exec touch {} \;
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

    #---------------
    # /wrf_df
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper/'wrf_df'
    start=`date +%s.%N`
    find $Big_dir/$Storm/$Exper/wrf_df -type l -exec touch --no-dereference {} \;
    find $Big_dir/$Storm/$Exper/wrf_df -exec touch {} \;
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

    #---------------
    # /run
    #---------------
    echo 'At' $Big_dir/$Storm/$Exper/'run'
    start=`date +%s.%N`
    # refresh spinup
    cd $Big_dir/$Storm/$Exper/run/$StartDate_spinup
    echo 'At '$StartDate_spinup
    find . -type l -exec touch --no-dereference {} \;
    find . -exec touch {} \;
    # refresh wrf-enkf cyclings
    DATE=$StartDate_cycling
    NEXTDATE=$StartDate_cycling
    while [[ $NEXTDATE -le $EndDate_cycling ]]; do 
      cd $Big_dir/$Storm/$Exper/run/$DATE
      echo 'At '$DATE
      # delete broken links
      echo 'Broken links are (if any)...'
      find ./enkf -xtype l
      find ./enkf -xtype l -exec rm {} \;
      # refresh time stamp of symbolic links
      find . -type l -exec touch --no-dereference {} \;  
      find . -exec touch {} \;
      NEXTDATE=`advance_time $DATE 60`
      DATE=$NEXTDATE
    done
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo 'It takes '$runtime's'

  done # End looping over Experiments
fi












# At this point (Jan 2, 2023), Judy hasn't figured out how to redirect output from python program to a log file
