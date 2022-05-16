#!/bin/bash --login
. $CONFIG_FILE
rundir=$WORK_DIR/run/$DATE/wrf_ens_fcst
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
wait_for_module ../update_bc ../icbc 

#Setup for wrf run
echo running > stat
echo "  Running WRF ensemble forecast..."

tid=0  #does not start from 0, because the wrf forecast runs with ens at the same time.
nt=$((total_ntasks/$wrf_ntasks))
for NE in `seq 1 $NUM_ENS_FCST`; do
  id=`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  touch $id/rsl.error.0000
  if [[ `tail -n1 $id/rsl.error.0000 |grep SUCCESS` ]]; then continue; fi

  cd $id
  ln -fs $WRF_DIR_FCST/run/* .
  rm -f namelist.*

  for n in `seq 1 $MAX_DOM`; do
    dm=d`expr $n + 100 |cut -c2-`
    if $FIRST_CYCLE; then
    	ln -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfinput_$dm
    	ncl $SCRIPT_DIR/util_change_nc_time.ncl 'ncfile="wrfinput_d01"' 'time="'`wrf_time_string $DATE`'"'
    else
    	ln -fs -fs $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id wrfrst_${dm}_`wrf_time_string $DATE` 
    fi
  done
  ln -fs $WORK_DIR/fc/$DATE/wrfbdy_d01_$id wrfbdy_d01

  if $FOLLOW_STORM; then
    cp $WORK_DIR/rc/$DATE/ij_parent_start .
    cp $WORK_DIR/rc/$DATE/domain_moves .
  fi
 
  if [[ $SST_UPDATE == 1 ]]; then
    ln -fs $WORK_DIR/rc/$LBDATE/wrflowinp_d?? .
  fi

  export start_date=$start_date_cycle
  if $LAST_CYCLE && $RUN_ENS_FCST ; then
    export run_minutes=`max $run_minutes_cycle $run_minutes_forecast`
  else
    export run_minutes=`diff_time $DATE $DATE_END`
    #export run_minutes=$run_minutes_cycle
  fi
  export wrf_for=forecast
  export inputout_interval=$run_minutes
  export inputout_begin=0
  export inputout_end=$run_minutes
  export GET_PHYS_FROM_FILE=false
  #$SCRIPT_DIR/namelist_wrf_realtime.sh wrfw $RUN_DOMAIN > namelist.input
  $SCRIPT_DIR/namelist_wrf_1km_ISFTCFLX.sh wrfw > namelist.input
  $SCRIPT_DIR/job_submit.sh $wrf_ntasks_fcst $((tid*$wrf_ntasks)) $HOSTPPN ./wrf.exe >& wrf.log &
  tid=$((tid+1))
  if [[ $tid == $nt ]]; then
    tid=0
    wait
  fi
  cd ..
done
wait

for NE in `seq 1 $NUM_ENS_FCST`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/rsl.error.0000 SUCCESS 1 $rundir
  rm $id/rsl.*  
  echo SUCCESS > $id/rsl.error.0000

  #outfile=$id/wrfinput_d01_`wrf_time_string $NEXTDATE`
  #outfile=$id/wrfrst_d01_`wrf_time_string $NEXTDATE`
  #watch_file $outfile 1 $rundir
  #mv $outfile $WORK_DIR/fc/$DATE/wrfinput_d01_`wrf_time_string $NEXTDATE`_$id
  #if [ $MAX_DOM -gt 1 ]; then
  #  for n in `seq 2 $MAX_DOM`; do
  #    dm=d`expr $n + 100 |cut -c2-`
  #    #outfile=$id/wrfout_${dm}_`wrf_time_string $NEXTDATE`
  #    outfile=$id/wrfrst_${dm}_`wrf_time_string $NEXTDATE`
  #    watch_file $outfile 1 $rundir
  #    mv $outfile $WORK_DIR/output/$DATE/wrfinput_${dm}_`wrf_time_string $NEXTDATE`_$id
  #  done
  #fi
done

if $CLEAN; then rm $rundir/$id/wrfout* ; fi
echo complete > stat
