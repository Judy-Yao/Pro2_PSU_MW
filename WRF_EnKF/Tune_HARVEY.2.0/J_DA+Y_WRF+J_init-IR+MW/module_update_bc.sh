#!/bin/bash --login
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/update_bc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi
cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#Check dependency
#----------------
#wait_for_module ../real
if $RUN_ENKF; then
  if [[ $DATE == $DATE_START ]]; then
    wait_for_module ../perturb_ic
  else
    wait_for_module ../enkf
  fi
fi

echo running > stat
echo "  Running UpdateBC..."

cat > parame.in << EOF
&control_param
 da_file               = 'wrfinput_d01_update'
 da_file_02            = 'ana02'
 wrf_bdy_file          = 'wrfbdy_d01_update'
 wrf_input             = 'wrfinput_d01'
 domain_id             = 1
 debug                 = .false. 
 update_lateral_bdy    = .true. 
 update_low_bdy        = .true.
 update_lsm            = .false.
 var4d_lbc             = .false.
/
EOF

#members
#tid=0
#nt=$total_ntasks
for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  if [[ ! -d $id ]]; then mkdir $id; fi
  cd $id
  ln -fs ../parame.in .
  ln -fs $WRFDA_DIR/var/da/da_update_bc.exe .

  # CM - modified boundary file to always use 
  # boundary condition file from previous run or
  # start date if first cycle. Also, updated to
  # randomly perturb all boundary times

  if [[ $DATE == $DATE_START ]]; then
    # determine total number of boundary times to update
    LBC_FORECAST_MIN=`diff_time $DATE_START $DATE_END`
    LBC_UPDATE_LOOP=$((${LBC_FORECAST_MIN}/${LBC_INTERVAL}))
    ln -fs $WORK_DIR/rc/$DATE/wrfbdy_d01 .
  elif [[ $DATE == $LBDATE ]]; then
    LBC_FORECAST_MIN=`diff_time $DATE_START $DATE_END`
    LBC_UPDATE_LOOP=$((${LBC_FORECAST_MIN}/${LBC_INTERVAL}))
    ln -fs $WORK_DIR/rc/$DATE/wrfbdy_d01 .
  else
    #echo $id
    LBC_FORECAST_MIN=`diff_time $DATE_START $DATE`
    LBC_UPDATE_LOOP=$((${LBC_FORECAST_MIN}/${LBC_INTERVAL} + 1))
    ln -fs $WORK_DIR/fc/$PREVDATE/wrfbdy_d01_$id wrfbdy_d01
    #ln -fs $WORK_DIR/fc/$PREVDATE/wrfbdy_d01 wrfbdy_d01
  fi

  ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01_$id wrfinput_d01
  cp -L wrfbdy_d01 wrfbdy_d01_update
  rm -f wrfinput_d01_update
  ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01_$id wrfinput_d01_update  

  ./da_update_bc.exe >& update_bc.log
#  tid=$((tid+1))
#  if [[ $tid == $nt ]]; then
#    tid=0
#    wait
#  fi
  cd ..
done
wait

for NE in `seq 1 $NUM_ENS`; do
  id=`expr $NE + 1000 |cut -c2-`
  watch_log $id/update_bc.log successfully 1 $rundir
  mv $id/wrfbdy_d01_update $WORK_DIR/fc/$DATE/wrfbdy_d01_$id
done

#ensemble mean (analysis for deterministic run)
if [[ ! -d mean ]]; then mkdir mean; fi
cd mean
ln -fs ../parame.in .
ln -fs $WRFDA_DIR/var/da/da_update_bc.exe .
#ln -fs $WORK_DIR/rc/$LBDATE/wrfbdy_d01 .
#ln -fs $WORK_DIR/rc/$LBDATE/wrfinput_d01 .
ln -fs $WORK_DIR/fc/$PREVDATE/wrfbdy_d01 wrfbdy_d01
cp -L wrfbdy_d01 wrfbdy_d01_update
ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01 wrfinput_d01
rm -f wrfinput_d01_update
ln -fs $WORK_DIR/fc/$DATE/wrfinput_d01 wrfinput_d01_update
./da_update_bc.exe >& update_bc.log
cd ..
watch_log mean/update_bc.log successfully 1 $rundir
mv mean/wrfbdy_d01_update $WORK_DIR/fc/$DATE/wrfbdy_d01

echo complete > stat

