#!/bin/bash 


#source ~/.bashrc
source /home1/06191/tg854905/.bashrc
#module restore intel17  # Intel Compiler Version 19
module restore intel

#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!
PRED_SCRIPT_DIR=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/Preprocess_Domain/FindD3Mesh

#load configuration files, functions, parameters
if [ -f config.HARVEY ]; then
  export CONFIG_FILE=$(pwd)/config.HARVEY
  source "$CONFIG_FILE"
fi
if [ -f util.sh ]; then
  source "$(pwd)/util.sh"
fi
PRED_DIR=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Preprocess_Domain/${STORM_ID}
#!!!!!!!!!!!!!!!!!! Attention !!!!!!!!!!!!!!!!!!!!!!!!!!!!

if [[ ! -d "$PRED_DIR" ]]; then mkdir -p "$PRED_DIR"; fi
cd "$PRED_DIR" || exit


#################################
#date

export DATE=$DATE_CYCLE_START
export PREVDATE=$DATE_CYCLE_START
export NEXTDATE=$DATE_CYCLE_START

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
  echo "Generating simulation domain for: " $DATE 

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
    ${PRED_SCRIPT_DIR}/namelist_wps.sh > namelist.wps

    touch geogrid.log
    
    if [[ ! `tail -n2 geogrid.log |grep Successful` ]]; then
       echo "  Running geogrid.exe..."
       ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
       #./geogrid.exe >& geogrid.log  &
       cat > run_wps_real.sh << EOF
#!/bin/bash -x
#SBATCH -J GetDomain
#SBATCH -p skx-dev
#SBATCH -N 1
#SBATCH -n 48 
#SBATCH -t 02:00:00
#SBATCH -o run_wps.batch

date
cd $rundir
ibrun ./geogrid.exe >& geogrid.log 
date
EOF
      sbatch run_wps_real.sh &>> job_submit_wps_real.log
      watch_log geogrid.log Successful 1 $rundir
    fi

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

cp config.HARVEY /work2/06191/tg854905/stampede2/Pro2_PSU_MW/Preprocess_Domain/${STORM_ID}/ 
