#!/bin/bash --login
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/icbc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#no dependency
if [[ $DATE -gt $DATE_START ]]; then
  wait_for_module ../../$DATE_START/icbc
fi

echo running > stat

#0. Calculate nested domain locations (centered over storm) and domain move steps
if $FOLLOW_STORM; then
  echo "  Calculating preset nesting..."
  #Nested domain location i,j: calculate from tcvatils if first cycle, otherwise get from previous cycle outputs
  if [ $DATE == $DATE_START ]; then
    $SCRIPT_DIR/calc_ij_parent_start.sh $DATE $WORK_DIR/rc/$DATE/ij_parent_start >& follow_storm.log
    watch_file $WORK_DIR/rc/$DATE/ij_parent_start 1 $rundir
    cp $WORK_DIR/rc/$DATE/ij_parent_start $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
  else
    if $RUN_ENKF; then
      i_parent_start="1 "
      j_parent_start="1 "
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_`wrf_time_string $DATE`_001
        #outfile=$WORK_DIR/fc/$DATE/wrfinput_${dm}_001
        watch_file $outfile 1 $rundir
        i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
        j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
      done
      echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start
      echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start
    fi
    if $RUN_4DVAR; then
      i_parent_start="1 "
      j_parent_start="1 "
      for n in `seq 2 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        if $RUN_ENKF; then
          outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)_mean
        else
          outfile=$WORK_DIR/fc/$PREVDATE/wrfinput_${dm}_$(wrf_time_string `advance_time $DATE $OBS_WIN_MIN`)
        fi
        watch_file $outfile 1 $rundir
        i_parent_start="$i_parent_start $(ncdump -h $outfile |grep :I_PARENT_START |awk '{print $3}')"
        j_parent_start="$j_parent_start $(ncdump -h $outfile |grep :J_PARENT_START |awk '{print $3}')"
      done
      echo $i_parent_start > $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
      echo $j_parent_start >> $WORK_DIR/rc/$DATE/ij_parent_start_4dvar
    fi
  fi

  #Domain move steps
  $SCRIPT_DIR/calc_domain_moves.sh $DATE $NEXTDATE $WORK_DIR/rc/$DATE/domain_moves >& follow_storm.log
  #$SCRIPT_DIR/calc_domain_moves.sh $DATE $DATE_END $WORK_DIR/rc/$DATE/domain_moves >& follow_storm.log
  watch_file $WORK_DIR/rc/$DATE/domain_moves 1 $rundir
  if $RUN_4DVAR; then
    if [ $DATE == $DATE_START ]; then
      $SCRIPT_DIR/calc_domain_moves.sh $DATE `advance_time $NEXTDATE $OBS_WIN_MIN` $WORK_DIR/rc/$DATE/domain_moves_4dvar >& follow_storm.log
    else
      $SCRIPT_DIR/calc_domain_moves.sh `advance_time $DATE $OBS_WIN_MIN` `advance_time $NEXTDATE $OBS_WIN_MIN` $WORK_DIR/rc/$DATE/domain_moves_4dvar >& follow_storm.log
    fi
    watch_file $WORK_DIR/rc/$DATE/domain_moves_4dvar 1 $rundir
  fi
  ln -fs $WORK_DIR/rc/$DATE/domain_moves .
  ln -fs $WORK_DIR/rc/$DATE/ij_parent_start .
fi


#if CP < LBC_INTERVAL, cannot generate wrfinput and wrfbdy from LBC data
#instead, we will fetch wrfbdy from the previous cycle where LBC is available
#and wrfinput will be from the previous cycle wrf run.

##################################################################
if [[ $LBDATE != $DATE ]]; then echo complete > stat; exit; fi
##################################################################

#######
#if $USE_PERTS; then
#  if [[ $LBDATE == $DATE ]]; then
#    cp $PERTS_DIR/rc/$DATE/wrfinput_d0? $WORK_DIR/rc/$DATE/. 
#    cp $PERTS_DIR/rc/$DATE/wrfbdy_d01 $WORK_DIR/rc/$DATE/. 
#  fi
#  echo complete > stat; exit #USING provided wps data and Perturbations!!
#fi
#######

export start_date=$DATE
#export run_minutes=`diff_time $start_date $DATE_END`
#export run_minutes=$((LBC_INTERVAL*2))
if [ $DATE == $DATE_START ]; then
  export run_minutes=`diff_time $DATE_START $NEXTDATE`
else
  #export run_minutes=$((LBC_INTERVAL*2))
  #export run_minutes=$FORECAST_MINUTES
  export run_minutes=`diff_time $DATE $DATE_END`
fi

#export run_minutes=`max $run_minutes $run_minutes_forecast`

$SCRIPT_DIR/namelist_wps.sh > namelist.wps

#1. geogrid.exe --------------------------------------------------------------------
#if [[ $DATE == $DATE_START ]]; then
touch geogrid.log
#if [[ ! `tail -n2 geogrid.log |grep Successful` ]]; then
#  echo "  Running geogrid.exe..."
  ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
#  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log 
#  #./geogrid.exe >& geogrid.log
#  watch_log geogrid.log Successful 1 $rundir
#fi
#  #mv geo_em.d0?.nc $WORK_DIR/rc/.
##fi
##ln -fs ../../../rc/geo_em.d0?.nc .

#2. ungrib.exe --------------------------------------------------------------------
fcst_hr=1000
fgdate=$start_date
touch ungrib.log
if [[ ! `tail -n2 ungrib.log |grep Successful` ]]; then
  echo "  Running ungrib.exe..."
  #Link first guess files (FNL, GFS or ECWMF-interim)
  fgdate=$start_date
  echo $fgdate
  gribfile=""
  while [[ $fgdate -le `advance_time $start_date $run_minutes` ]]; do
    #GFS
    ccyymm=`echo $start_date |cut -c1-6`
    dd=`echo $start_date |cut -c7-8`
    hh=`echo $start_date |cut -c9-10`
    file_ext=`echo $fcst_hr |cut -c2-4`
    file=$FG_DIR_GFS/gfs.0p25.${ccyymm}${dd}${hh}.f${file_ext}.grib2 #GFS Forecast
    echo $file
    #ccyymm=`echo $fgdate |cut -c1-6`
    #dd=`echo $fgdate |cut -c7-8`
    #hh=`echo $fgdate |cut -c9-10`
    #file=$FG_DIR/gfs.0p25.${ccyymm}${dd}${hh}.f000.grib2 #GFS Analysis
    #FNL
    #ccyymm=`echo $fgdate |cut -c1-6`
    #dd=`echo $fgdate |cut -c7-8`
    #hh=`echo $fgdate |cut -c9-10`
    #file="$FG_DIR/gfs.$ccyymm$dd$hh/`date -u -d $ccyymm$dd' '$hh':00' +%y%j%H`000000" #GFS
    #file="$FNL_DIR/fnl_${ccyymm}${dd}_${hh}_00.grib2"      #FNL
    #file="$FG_DIR/gdas1.fnl0p25.${ccyymm}${dd}${hh}.f00.grib2"      #FNL
    #file="$FNL_DIR/ERA-Interim_pl_${ccyymm}${dd}_${hh}_00.grib"      #ERA-Interim
    #echo $file
    if [ -e $file ]; then 
      gribfile="$gribfile $file"
    fi
    ## for ERA-Interim
    #file="$FNL_DIR/ERA-Interim_sfc_${ccyymm}${dd}_${hh}_00.grib"      #ERA-Interim
    #if [ -e $file ]; then
    #  gribfile="$gribfile $file"
    #fi
    fgdate=`advance_time $fgdate $LBC_INTERVAL`
    fcst_hr=$((${fcst_hr}+6))
    #fcst_hr=$((${fcst_hr}+3))
  done
  echo $gribfile
  $WPS_DIR/link_grib.csh $gribfile
  ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFS Vtable
  #ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.ERA-interim.pl Vtable
  ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
  ./ungrib.exe >& ungrib.log
  watch_log ungrib.log Successful 2 $rundir
fi

#3. metgrid.exe --------------------------------------------------------------------
#touch metgrid.log
#if [[ ! `tail -n2 metgrid.log |grep Successful` ]]; then
#  echo "  Running metgrid.exe..."
  ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
  ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
#  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log
#  #./metgrid.exe >& metgrid.log
#  watch_log metgrid.log Successful 2 $rundir
#  #mv met_em* $WORK_DIR/rc/$DATE/.
#fi

#4. real.exe ----------------------------------------------------------------------
touch rsl.error.0000
if [[ ! `tail -n2 rsl.error.0000 |grep SUCCESS` ]]; then
  echo "  Running geogrid.exe, metgrid.exe, and real.exe..."
  export NUM_METGRID_LEVELS=32
  export NUM_METGRID_SOIL_LEVELS=4
  export GET_PHYS_FROM_FILE=false
  #$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
  $SCRIPT_DIR/namelist_real.sh real 1 > namelist.input
  #ln -fs ../../../rc/$DATE/met_em* .
  ln -fs $WRF_DIR/main/real.exe .
#  $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./real.exe >& real.log

  cat > run_wps_real.sh << EOF
#!/bin/bash -x
#SBATCH -A pen116
#SBATCH -J wps
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 48 
#SBATCH -t 0:30:00
#SBATCH -o run_wps.batch
#SBATCH --cpus-per-task=1

date
cd $rundir
ibrun ./geogrid.exe >& geogrid.log 
ibrun ./metgrid.exe >& metgrid.log
ibrun ./real.exe >& real.log 
date
EOF
  sbatch run_wps_real.sh &>> job_submit_wps_real.log
  #./real.exe >& real.log
  watch_log rsl.error.0000 SUCCESS 2 $rundir
  if [ $SST_UPDATE == 1 ]; then
    if [ $CYCLE_PERIOD -lt $LBC_INTERVAL ]; then
      for n in `seq 1 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
		module restore default
		ncl $SCRIPT_DIR/util_linint_nc_time.ncl dmin=$CYCLE_PERIOD 'infile="wrflowinp_'$dm'"' >> lowinp.log 2>&1
		module restore intel
		mv tmp.nc $WORK_DIR/rc/$DATE/wrflowinp_$dm
      done
    else
      cp wrflowinp_d?? $WORK_DIR/rc/$DATE
    fi
  fi
fi

cp wrfinput_d?? $WORK_DIR/rc/$DATE/.
cp wrfbdy_d01 $WORK_DIR/rc/$DATE/.
cp namelist.input $WORK_DIR/rc/$DATE/.
cp namelist.wps $WORK_DIR/rc/$DATE/.
if [[ $DATE == $DATE_START ]]; then
  cp wrfinput_d?? $WORK_DIR/fc/$DATE/.
  cp wrfbdy_d01 $WORK_DIR/fc/$DATE/.
  if $RUN_4DVAR; then
    cp $WORK_DIR/fc/wrfbdy_d01 $WORK_DIR/fc/wrfbdy_d01_window
  fi
fi

if $CLEAN; then rm -f *log.???? met_em* *FILE*; fi
echo complete > stat
