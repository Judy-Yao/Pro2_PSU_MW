#!/bin/bash -x

module restore intel
export CONFIG_FILE=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/WRF_EnKF/IRMA_run/IR-THO_onlyTCvitals/config_IR_THO.IRMA.S3

source "$CONFIG_FILE"
source "$SCRIPT_DIR"/util.sh

if [[ ! -d "$WORK_DIR" ]]; then mkdir -p "$WORK_DIR"; fi
cd "$WORK_DIR" || exit

#################################
#date

export DATE=$DATE_START
export PREVDATE=$DATE_START
export NEXTDATE=$DATE_START

while [[ $NEXTDATE -le $DATE_CYCLE_END ]]; do
 
  export CP=$(min ${CYCLE_PERIOD[@]})
  if [[ $DATE == $DATE_START ]]; then
    export CP=$(diff_time "$DATE" "$DATE_CYCLE_START")
  fi

  # -----------------------------------------------------------------
  # calculate start_date and run_minutes, used by namelist_wrf.sh to &
  # generate correct time in namelist.input
  # -----------------------------------------------------------------
  export start_date_cycle=$DATE
  export run_minutes_cycle=$CP

  # ------------------
  # LBDATE
  # ------------------
  export minute_off=$(echo "(${start_date_cycle:8:2}*60+${start_date_cycle:10:2})%$BC_INTERVAL" |bc)
  if [[ $DATE == $(advance_time "$start_date_cycle" -"$minute_off") ]]; then
    export LBDATE=$(advance_time "$start_date_cycle" -"$minute_off")
  fi

  # ------------------
  # FCSTDATE
  # ------------------
  export NEXTDATE=$(advance_time "$DATE" "$CP")
  echo "----------------------------------------------------------------------"
  echo "CURRENT CYCLE: $(wrf_time_string "$DATE") => $(wrf_time_string "$NEXTDATE")"

  mkdir -p $WORK_DIR/{run,fc}/"$DATE"

  # ------------------
  # clear error tags
  # ------------------
  for d in run/"$DATE"/*; do
    if [[ $(cat "$d"/stat) != "complete" ]]; then
      echo waiting > "$d"/stat
    fi
  done

  # ------------------
  # run components
  # ------------------
  if [ "$DATE" -ge "$DATE_CYCLE_START" ] && [ "$DATE" -le "$DATE_CYCLE_END" ]; then
    echo "CURRENT CYCLE: $(wrf_time_string "$DATE")"
    mkdir -p $WORK_DIR/run/"$DATE"/enkf/
    "$SCRIPT_DIR"/module_enkf_part1.sh $DATE $PREVDATE
    #"$SCRIPT_DIR"/module_enkf_part2.sh $DATE $PREVDATE
  fi

  wait
  date

  # ------------------
  # check errors  
  # ------------------
  #for d in $(ls -t run/"$DATE"/); do
  for d in run/"$DATE"/*; do
    if [[ $(cat run/"$DATE"/"$d"/stat) == "error" ]]; then
      echo CYCLING STOP DUE TO FAILED COMPONENT: "$d"
      exit 1
    fi
  done

  cp $CONFIG_FILE ${WORK_DIR}/run/$DATE

  # ------------------
  #advance to next cycle
  # ------------------
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE


