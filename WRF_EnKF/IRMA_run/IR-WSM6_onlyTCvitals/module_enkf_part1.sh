#!/bin/bash --login
. "$CONFIG_FILE"
module restore intel

rundir=$WORK_DIR/run/$DATE/enkf
if [[ ! -d $rundir ]]; then mkdir -p "$rundir"; echo waiting > "$rundir"/stat; fi
cd "$rundir" || exit
if [[ $(cat stat) == "complete" ]]; then exit; fi

DATE=$1 
PREVDATE=$2

#Run EnKF
echo running > stat
echo "  Running EnKF..."

domlist=`seq 3 $MAX_DOM`
#link files
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  if [[ ! -d $dm ]]; then mkdir -p $dm; fi
  if [ -f $dm/${DATE}.finish_flag ]; then continue; fi
  cd $dm
  lfs setstripe -c 1 $rundir/$dm

  echo "    Linking files for domain $dm"
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    if $USE_ESTIMATE_INF_RADIANCE; then
      cp -L $CONTROL_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    elif $USE_INF_BT_RADIANCE; then
      cp -L $CONTROL_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    else
      ln -fs $CONTROL_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_$id fort.`expr 80010 + $NE`
    fi
    cp -L fort.`expr 80010 + $NE` fort.`expr 90010 + $NE` >> link.log 2>&1 &
  done
  wait

  #cp -L $WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE` fort.`expr 80011 + $NUM_ENS`
  cp -L fort.80011 fort.`expr 80011 + $NUM_ENS`
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 90011 + $NUM_ENS` 
  cp -L fort.`expr 80011 + $NUM_ENS` fort.`expr 60011 + $NUM_ENS`

  ln -fs $WRF_DIR/run/LANDUSE.TBL .
  # coefficients for CRTM
  # ln -fs $CRTM_DIR/crtm_wrf/coefficients .
  ln -fs $CRTM_COEFF_DIR ./coefficients
  ln -sf $TCVITALS_DIR/${DATE:0:4}/${DATE}.${STORM_ID}-tcvitals.dat tcvitals.dat
  #Observations
  #LITTLE_R format from obsproc
  if ( $INCLUDE_LITTLE_R || $INCLUDE_BUFR || $INCLUDE_MADIS || $INCLUDE_TCV); then
    ln -fs $CONTROL_DIR/obs/$DATE/obs_gts_`wrf_time_string $DATE`.3DVAR obs_3dvar_${DATE}00
  fi
  #AIRBORNE superobs
  if ( $INCLUDE_AIRBORNE ); then
    ln -fs $CONTROL_DIR/obs/$DATE/airborne_`echo $DATE |cut -c1-12`_so .
  fi
  #HURRICANE PI
  if ( $USE_HURRICANE_PI ); then
    ln -sf $TCVITALS_DIR/${DATE:0:4}/${DATE}.${STORM_ID}-tcvitals.dat hurricane_best_track
  fi
  #Radiance
  if ( $USE_RADIANCE ); then
    ln -fs $DATA_DIR/radiance/radiance_${dm}_${DATE}_so radiance_${DATE}_so  # for when $DATA_DIR includes $STORM_ID
    #ln -fs $DATA_DIR/radiance/${STORM_ID}/radiance_${dm}_${DATE}_so radiance_${DATE}_so
  fi
  #Microwave
  if ( $USE_MICROWAVE ); then
    ln -fs $DATA_DIR/microwave/microwave_${dm}_${DATE}_so microwave_${DATE}_so
    #ln -fs $DATA_DIR/microwave/${STORM_ID}/microwave_${dm}_${DATE}_so microwave_${DATE}_so
  fi

  # Lu YH
  # If Radiance and Microwave both exist for one assimilation time, only assimilate IR
  # for clear pixels and only assimilate MW for cloudy pixels
  if ( $USE_RADIANCE && $USE_MICROWAVE ); then
    if [ -e radiance_${DATE}_so  ] && [ -e microwave_${DATE}_so ] ; then
      
      if [ "$MERGE_IRMW_METHOD" != "" ]; then
        merge_script=merge_${MERGE_IRMW_METHOD}.py
        if [ -e $SCRIPT_DIR/${merge_script} ]; then
          mv radiance_${DATE}_so radiance_${DATE}_so_orig # radiance data will NOT be assimilated
          mv microwave_${DATE}_so microwave_${DATE}_so_orig

          ln -sf $SCRIPT_DIR/${merge_script} .
          ./${merge_script} radiance_${DATE}_so_orig microwave_${DATE}_so_orig fort.80011 radiance_${DATE}_so microwave_${DATE}_so
        else
          echo "MERGE_IR_MW_METHOD not supported"
        fi
      fi
    fi
  fi
  export CRTM_OUT_DIR=${WORK_DIR}/run/${DATE}/enkf/${dm}/crtm_out_${CLOUDCOEFF1}
  echo $CRTM_OUT_DIR
  if [[ ! -d ${CRTM_OUT_DIR} ]]; then mkdir -p ${CRTM_OUT_DIR}; touch ${CRTM_OUT_DIR}/test; fi # SBS


  # updating non-Q variables every 1-hour
  if [[ ${DATE:10:2} == '00'  ]]; then
    ln -fs $ENKF_DIR/enkf.mpi .
  else
    ln -fs $ENKF_DIR/enkf_hydro.mpi enkf.mpi
  fi

  # multiplicative inflation
  if $USE_ESTIMATE_INF_RADIANCE; then
    ln -fs $ENKF_DIR/cal_inflate.mpi .
    ln -sf $CONTROL_DIR/run/$PREVDATE/enkf/$dm/parameters_update${PREVDATE} parameters_update
  fi
  if $USE_INF_BT_RADIANCE; then
    ln -sf $CONTROL_DIR/run/$PREVDATE/enkf/$dm/parameters_update${PREVDATE} parameters_update
  fi

  cp $SCRIPT_DIR/namelist.enkf .  
  #$SCRIPT_DIR/namelist_enkf.sh $n > namelist.enkf
  #cd ..

echo "    Running enkf.mpi"
cat > run_enkf.sh << EOF
#!/bin/bash -x
#SBATCH -J enkf
#SBATCH -p skx-dev
#SBATCH -n 192 -N 4
#SBATCH -t 0:30:00
#SBATCH -o enkf.batch

##SBATCH -J enkf
##SBATCH -A pen116
##SBATCH -p compute
##SBATCH --nodes=2
##SBATCH --ntasks-per-node=128
##SBATCH --cpus-per-task=1
##SBATCH -t 0:10:00
##SBATCH -o enkf.batch

date
for n in \`seq 3 $MAX_DOM\`; do
  dm=d\`expr \$n + 100 |cut -c2-\`
  cd $rundir/\$dm
  ibrun -n 192 ./enkf.mpi >& enkf.log
  date
done
EOF
sbatch run_enkf.sh &>> job_submit.log

done

#Check output
for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  watch_log $rundir/$dm/enkf.log Successful 3 $rundir
  #watch_log $dm/${DATE}.finish_flag _ 5 $rundir
done


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
      mv ../../fort.`expr 90010 + $NE` wrfinput
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
    cd $rundir/$dm
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




