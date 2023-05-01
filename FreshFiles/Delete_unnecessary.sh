#!/bin/bash
# Delete files that are not necessary for each experiment
# Author: Zhu (Judy) Yao. April30, 2023

source util.sh

#---------------------------------------------------------
# Configuration
#---------------------------------------------------------
Storm=MARIA
Exper=IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900
StartDate_cycling=201709160000
EndDate_cycling=201709180000
NUM_ENS=60
MAX_DOM=3

if_clean_enkf=false
if_clean_wrf_ens=true

# Set up directories
source_dir=/scratch/06191/tg854905/Pro2_PSU_MW/


if $if_clean_enkf; then
  DATE=$StartDate_cycling
  NEXTDATE=$StartDate_cycling
  while [[ $NEXTDATE -le $EndDate_cycling ]]; do
    for n in `seq 1 $MAX_DOM`; do
      dm=d`expr $n + 100 |cut -c2-`
      for NE in `seq 1 $NUM_ENS`; do
        input_id=fort.`expr 80010 + $NE`
        rm $source_dir/$Storm/$Exper/run/$DATE/enkf/$dm/$input_id
        output_id=fort.`expr 90010 + $NE`
        rm $source_dir/$Storm/$Exper/run/$DATE/enkf/$dm/$output_id
      done    
     done
   NEXTDATE=`advance_time $DATE 60`
   DATE=$NEXTDATE
  done 
fi

if $if_delete_wrf_ens; then
  DATE=$StartDate_cycling   
  NEXTDATE=$StartDate_cycling
  while [[ $NEXTDATE -le $EndDate_cycling ]]; do 
    rm -r $source_dir/$Storm/$Exper/run/$DATE/wrf_ens
    NEXTDATE=`advance_time $DATE 60`
    DATE=$NEXTDATE
  done
fi
