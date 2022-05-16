#!/bin/bash --login
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/icbc
if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi

#no dependency

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

if [[ $LBDATE != $DATE ]]; then echo complete > stat; exit; fi
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

export run_minutes=`diff_time $DATE $DATE_END`

$SCRIPT_DIR/namelist_wps.sh > namelist.wps

#1. geogrid.exe --------------------------------------------------------------------
#if [[ $DATE == $DATE_START ]]; then
touch geogrid.log
if [[ ! `tail -n2 geogrid.log |grep Successful` ]]; then
echo "  Running geogrid.exe..."
ln -sf $WPS_DIR/geogrid/src/geogrid.exe .
$SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./geogrid.exe >& geogrid.log 
#./geogrid.exe >& geogrid.log
watch_log geogrid.log Successful 1 $rundir
fi

mcount=1000
for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    if [[ ! -d $id ]]; then mkdir $id; fi
    cd ${id}

    cp ../namelist.wps .
    
    if [[ $NE -lt 21 ]]; then 
        init_date=2015102012
        FG_DIR2=${FG_DIR}/${init_date}
        fcst_hr=1012
    elif [[ $NE -lt 41 ]]; then 
        init_date=2015102018
        FG_DIR2=${FG_DIR}/${init_date}
        if [[ $NE -eq 21 ]]; then
            mcount=1000
        fi
        fcst_hr=1006
    else 
        init_date=2015102100
        FG_DIR2=${FG_DIR}/${init_date}
        if [[ $NE -eq 41 ]]; then
            mcount=1000
        fi
        fcst_hr=1000
    fi
    mcount=`expr $mcount + 1`
    #2. ungrib.exe --------------------------------------------------------------------
    fgdate=$start_date
    touch ungrib.log
    if [[ ! `tail -n2 ungrib.log |grep Successful` ]]; then
    echo "  Running ungrib.exe..."
    #Link first guess files (FNL, GFS or ECWMF-interim)
    fgdate=$start_date
    gribfile=""
    while [[ $fgdate -le `advance_time $start_date $run_minutes` ]]; do
        ccyymmdd=`echo $init_date |cut -c1-8`
        hh=`echo $init_date |cut -c9-10`
        file_ext=`echo $fcst_hr |cut -c2-4`
        ##GEFS
        mem=`echo ${mcount} |cut -c3-`
        file1=$FG_DIR2/gens-a_3_${ccyymmdd}_${hh}00_${file_ext}_${mem}.grb2
        file2=$FG_DIR2/gens-b_3_${ccyymmdd}_${hh}00_${file_ext}_${mem}.grb2
        echo $file1
        echo $file2
        if [ -e $file1 ]; then 
            gribfile="$gribfile $file1"
        fi
        if [ -e $file2 ]; then 
            gribfile="$gribfile $file2"
        fi
        fgdate=`advance_time $fgdate $LBC_INTERVAL`
        fcst_hr=$((${fcst_hr}+6))
    done
    echo $gribfile
    $WPS_DIR/link_grib.csh $gribfile
    ln -sf $WPS_DIR/ungrib/Variable_Tables/Vtable.GFSENS Vtable
    ln -fs $WPS_DIR/ungrib/src/ungrib.exe .
    ./ungrib.exe >& ungrib.log
    watch_log ungrib.log Successful 2 $rundir
    fi
    cd ..
done

#3. metgrid.exe --------------------------------------------------------------------
for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    cd ${id}
    ln -fs ../geo_em.d0?.nc .
    touch metgrid.log
    if [[ ! `tail -n2 metgrid.log |grep Successful` ]]; then
    echo "  Running metgrid.exe..."
    ln -fs $WPS_DIR/metgrid/METGRID.TBL.ARW METGRID.TBL
    ln -fs $WPS_DIR/metgrid/src/metgrid.exe .
    $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./metgrid.exe >& metgrid.log &
    #./metgrid.exe >& metgrid.log
    #watch_log metgrid.log Successful 2 $rundir
    #mv met_em* $WORK_DIR/rc/$DATE/.
    fi
    cd ..
done

#4. real.exe ----------------------------------------------------------------------
for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    cd ${id}
    #watch_log metgrid.log Successful 2 $rundir
    touch rsl.error.0000
    if [[ ! `tail -n2 rsl.error.0000 |grep SUCCESS` ]]; then
       echo "  Running real.exe..."
        export NUM_METGRID_LEVELS=27
        export NUM_METGRID_SOIL_LEVELS=4
        export GET_PHYS_FROM_FILE=false
        ln -fs ../ij_parent_start .
        ln -fs ../domain_moves .
        #$SCRIPT_DIR/namelist_wrf.sh real > namelist.input
        $SCRIPT_DIR/namelist_real.sh real 1 > namelist.input
        #ln -fs ../../../rc/$DATE/met_em* .
        ln -fs $WRF_DIR/main/real.exe .
        $SCRIPT_DIR/job_submit.sh $real_ntasks 0 $HOSTPPN ./real.exe >& real.log &
        #./real.exe >& real.log
        #watch_log rsl.error.0000 SUCCESS 2 $rundir

        if [ $SST_UPDATE == 1 ]; then
        if [ $CYCLE_PERIOD -lt $LBC_INTERVAL ]; then
            for n in `seq 1 $MAX_DOM`; do
                dm=d`expr $n + 100 |cut -c2-`
                ncl $SCRIPT_DIR/util_linint_nc_time.ncl dmin=$CYCLE_PERIOD 'infile="wrflowinp_'$dm'"' >> lowinp.log 2>&1
                mv tmp.nc $WORK_DIR/rc/$DATE/wrflowinp_$dm
            done
        else
            cp wrflowinp_d?? $WORK_DIR/rc/$DATE
        fi
        fi
    fi
    cd ..
done

for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    cd ${id}
    watch_log rsl.error.0000 SUCCESS 2 $rundir
    # Move files
    cp namelist.wps $WORK_DIR/rc/$DATE/namelist.wps_${id}
    cp namelist.input $WORK_DIR/rc/$DATE/namelist.input_${id}
    cp wrfbdy_d01 $WORK_DIR/rc/$DATE/wrfbdy_d01_${id}
    for n in `seq 1 $MAX_DOM`; do
        dm=d`expr $n + 100 |cut -c2-`
        cp wrfinput_${dm} $WORK_DIR/rc/$DATE/wrfinput_${dm}_${id}
        if [[ $DATE == $DATE_START ]]; then
            cp wrfinput_${dm} $WORK_DIR/fc/$DATE/wrfinput_${dm}_${id}
            if [[ $n -eq 1 ]]; then
                cp wrfbdy_d01 $WORK_DIR/fc/$DATE/wrfbdy_d01_${id}
            fi
        if $RUN_4DVAR; then
            cp $WORK_DIR/fc/wrfbdy_d01 $WORK_DIR/fc/wrfbdy_d01_window
        fi
        fi
    done
    cd ..
done

if $CLEAN; then rm -f *log.???? met_em* *FILE*; fi
echo complete > stat
