#!/bin/bash --login
#####header for jet######
#PBS -A hfip-psu
#PBS -N run_cycle
#PBS -l walltime=8:00:00
#PBS -q batch
#PBS -l partition=sjet:vjet:xjet
#PBS -l procs=1024
#PBS -j oe
#PBS -d .

#####header for stampede######
##SBATCH -J run_cycle
##SBATCH -n 1024
##SBATCH -n 64
##SBATCH -p normal
##SBATCH -p development
##SBATCH -t 12:00:00
##SBATCH -t 2:00:00
##SBATCH -o out_enkf
##SBATCH -e error_enkf

#source ~/.bashrc

#load configuration files, functions, parameters
cd $WORK/Nystrom_Satellite
export CONFIG_FILE=$WORK/Nystrom_Satellite/config_fcst
. $CONFIG_FILE
. util.sh

if [[ ! -d $WORK_DIR ]]; then mkdir -p $WORK_DIR; fi
cd $WORK_DIR

####total_ntasks####
if [ $JOB_SUBMIT_MODE == 1 ]; then
  if [[ $HOSTTYPE == "stampede" ]]; then
    export total_ntasks=$SLURM_NTASKS
  fi
  if [[ $HOSTTYPE == "jet" ]]; then
    export total_ntasks=$PBS_NP
  fi
else
  export total_ntasks=9999999
fi

#################################
#date

export DATE=$DATE_START
export PREVDATE=$DATE_START
export NEXTDATE=$DATE_START

while [[ $NEXTDATE -le $DATE_CYCLE_END ]]; do

  export OWMAX=`min ${OBS_WIN_MAX[@]}`
  export OWMIN=`min ${OBS_WIN_MIN[@]}`
  export DT=$TIME_STEP
  export MPS=`min ${MINUTE_PER_SLOT[@]}`
  export FCSTM=`min ${FORECAST_MINUTES[@]}`
  export CP=`min ${CYCLE_PERIOD[@]}`
  if [[ $DATE == $DATE_START ]]; then
    export CP=`diff_time $DATE $DATE_CYCLE_START`
  fi

#calculate start_date and run_minutes, used by namelist_wrf.sh to generate correct
#time in namelist.input
  if $RUN_4DVAR; then
    export start_date_cycle=$DATE
    export run_minutes_cycle=`echo $CP+$OWMAX |bc`
    if [[ $DATE -ge $DATE_CYCLE_START ]]; then
      export start_date_cycle=`advance_time $start_date_cycle $OWMIN`
      export run_minutes_cycle=`echo $run_minutes+$OWMIN |bc`
    fi
  else
    export start_date_cycle=$DATE
    export run_minutes_cycle=$CP
  fi
  if $RUN_DETERMINISTIC; then
    export run_minutes_forecast=`diff_time $DATE $DATE_END`
  else
    export run_minutes_forecast=`max $CP $FCSTM`
  fi

#LBDATE
  #export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$LBC_INTERVAL" |bc`
  export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$BC_INTERVAL" |bc`
  if [[ $DATE == `advance_time $start_date_cycle -$minute_off` ]]; then
    export LBDATE=`advance_time $start_date_cycle -$minute_off`
    #export LBDATE=$DATE_START
  fi
#FCSTDATE
  export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$FCST_INTERVAL" |bc`
  export FCSTDATE=`advance_time $start_date_cycle -$minute_off`
#ENS_FCSTDATE
  export minute_off=`echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$ENS_FCST_INTERVAL" |bc`
  export ENSFCSTDATE=`advance_time $start_date_cycle -$minute_off`

  export NEXTDATE=`advance_time $DATE $CP`
  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: `wrf_time_string $DATE` => `wrf_time_string $NEXTDATE`"

  mkdir -p {run,rc,fc,output,obs}/$DATE

#clear error tags
  for d in `ls run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` != "complete" ]]; then
      echo waiting > run/$DATE/$d/stat
    fi
  done
#run components
  #$SCRIPT_DIR/module_wps.sh &
  #$SCRIPT_DIR/module_real.sh &
  $SCRIPT_DIR/module_icbc.sh &
  if $RUN_ENKF && [ $DATE == $DATE_START ]; then
    $SCRIPT_DIR/module_perturb_ic.sh &
    $SCRIPT_DIR/module_update_bc.sh &
    $SCRIPT_DIR/module_wrf_ens.sh &
  fi
  if [ $DATE == $LBDATE ]; then
    $SCRIPT_DIR/module_perturb_ic.sh &
  fi
  if [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ]; then
    if $RUN_ENKF || $RUN_4DVAR; then
      $SCRIPT_DIR/module_obsproc.sh &
    fi
    if $RUN_ENKF; then
      $SCRIPT_DIR/module_enkf.sh &
    fi
    if $RUN_4DVAR; then
      $SCRIPT_DIR/module_4dvar.sh &
    fi
    if $RUN_ENKF || $RUN_4DVAR; then
      $SCRIPT_DIR/module_update_bc.sh &
    fi
    if $RUN_ENKF; then
      $SCRIPT_DIR/module_wrf_ens.sh &
    fi
  fi
  if [ $DATE == $FCSTDATE ] && [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ] ; then
    $SCRIPT_DIR/module_wrf.sh &
    #$SCRIPT_DIR/module_wrf_3km.sh &
    #$SCRIPT_DIR/module_wrf_ISFTCFLX.sh &
    if $RUN_ENS_FCST; then 
        $SCRIPT_DIR/module_wrf_ens_fcst.sh &
    fi
  fi
  if $RUN_ENS_FCST && [ $DATE == $ENSFCSTDATE ] && [ $DATE -ge $DATE_CYCLE_START ] && [ $DATE -le $DATE_CYCLE_END ] ; then
    $SCRIPT_DIR/module_wrf_ens_fcst.sh &
  fi

  wait
  date

#check errors  
  for d in `ls -t run/$DATE/`; do
    if [[ `cat run/$DATE/$d/stat` == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: $d
      exit 1
    fi
  done

#advance to next cycle
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE
echo bottom $MODULEPATH_ROOT
