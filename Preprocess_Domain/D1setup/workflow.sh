#!/bin/bash 


#source ~/.bashrc
source /home1/06191/tg854905/.bashrc

#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!
PRED_SCRIPT_DIR=$(pwd)

#load configuration files, functions, parameters
if [ -f config.MARIA_test ]; then
  export CONFIG_FILE=$(pwd)/config.MARIA_test
  source "$CONFIG_FILE"
fi
if [ -f util.sh ]; then
  source "$(pwd)/util.sh"
fi
PRED_DIR=$(pwd)
#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!



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

  # ------------------
  # FCSTDATE
  # ------------------
  export NEXTDATE=$(advance_time "$DATE" "$CP")
  echo "----------------------------------------------------------------------"
  echo $DATE 

  rundir=$PRED_DIR/$DATE
  if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

  cd $rundir
  if [[ `cat stat` == "complete" ]]; then exit; fi

  # ------------------
  # run components
  # ------------------
  if $FOLLOW_STORM; then
    echo "  Calculating IJ_PARENT_START using linear interpolation..."
    ${PRED_SCRIPT_DIR}/calc_ij_parent_start.sh $DATE ${rundir}/ij_parent_start >& find_ij_parent_start.log
    watch_file $rundir/ij_parent_start 1 $rundir    

  wait
  date
        
    # ------------------
    # check errors  
    # ------------------
    if $CLEAN; then rm -f *log.????; fi
    echo complete > stat

  fi

  # ------------------
  #advance to next cycle
  # ------------------
  export PREVDATE=$DATE
  export DATE=$NEXTDATE

done
echo CYCLING COMPLETE
echo bottom "$MODULEPATH_ROOT"

