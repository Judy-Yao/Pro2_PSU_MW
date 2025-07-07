#!/bin/bash --login
. "$CONFIG_FILE"

module restore intel

rundir=$WORK_DIR/run/$DATE/enkf
domlist=`seq 3 $MAX_DOM`

if [[ ! -d $rundir ]]; then mkdir -p "$rundir"; echo waiting > "$rundir"/stat; fi
if [[ $(cat stat) == "complete" ]]; then exit; fi
cd "$rundir" || exit

DATE=$1 
PREVDATE=$2


### replacing outside with first guess (GFS/FNL) reanalysis
if $REPLACE_ENVIRONMENT; then
 nt=$((total_ntasks/$HOSTPPN))
 if [ $DATE == $LBDATE ]; then
  tid=0
  for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    if [[ ! -d replace_environment ]]; then mkdir -p replace_environment; fi
    cd replace_environment
    echo "  Replacing environment with GFS for domain $dm"
    for NE in `seq 1 $((NUM_ENS+1))`; do
      id=`expr $NE + 1000 |cut -c2-`
      if [[ ! -d $id ]]; then mkdir $id; fi
      if [ -e ${id}/replace_environment.log ]; then if [[ `tail -n2 ${id}/replace_environment.log |grep Successful` ]]; then continue; fi; fi
      cd $id

      ln -sf $ENKF_DIR/replace_environment_by_gfs.exe .
      ln -sf $TCVITALS_DIR/${DATE:0:4}/${DATE}.${STORM_ID}-tcvitals.dat tcvitals.dat
      #mv ../../fort.`expr 90010 + $NE` wrfinput
      if  [[ $NE == `expr $((NUM_ENS+1))` ]]; then
        ln -sf $CONTROL_DIR/rc/$DATE/wrfinput_${dm} wrfinput_gfs
      else
        ln -sf $CONTROL_DIR/rc/$DATE/wrfinput_${dm}_$id wrfinput_gfs
      fi
      ./replace_environment_by_gfs.exe >& replace_environment.log
      tid=$((tid+1))
      if [[ $tid == $nt ]]; then
        tid=0
        wait
      fi
      cd ..
    done
    cd ../..
  done

  for n in $domlist; do
#  for n in `seq 1 $((MAX_DOM-1))`; do
    dm=d`expr $n + 100 |cut -c2-`
    cd $dm
    for NE in `seq 1 $((NUM_ENS+1))`; do
      id=`expr $NE + 1000 |cut -c2-`
      watch_log replace_environment/$id/replace_environment.log Successful 1 $rundir
      mv replace_environment/$id/wrfinput fort.`expr 90010 + $NE` 
    done
    cd ..
  done
 fi
fi
###

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  #cp $dm/fort.10000 $WORK_DIR/obs/$DATE/assimilated_obs_${dm}_fort.10000
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    mv $dm/fort.`expr 80010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    mv $dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id $dm/fort.`expr 80010 + $NE`
    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $dm/fort.`expr 90010 + $NE`
    ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
  done
  cp $dm/fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  cp $dm/fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean
  ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean $WORK_DIR/fc/$DATE/wrfinput_${dm}
done


echo complete > stat

