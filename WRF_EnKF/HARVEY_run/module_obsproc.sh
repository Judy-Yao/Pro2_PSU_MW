#!/bin/bash --login
. $CONFIG_FILE

rundir=$WORK_DIR/run/$DATE/obsproc

if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi

cd $rundir
if [[ `cat stat` == "complete" ]]; then exit; fi
echo running > stat
echo "  Running ObsProc..."

ln -fs $WRFDA_DIR/var/obsproc/obsproc.exe .
ln -fs $WRFDA_DIR/var/obsproc/obserr.txt .
echo > obs.raw

##### include NCAR_LITTLE_R (3-hourly) #####
if $INCLUDE_LITTLE_R; then
  #cp $DATA_DIR/sfcobs/SURFACE_OBS:${DATE:0:10} $DATA_DIR/littler/little_r:${DATE:0:4}-${DATE:4:2}-${DATE:6:2}_${DATE:8:2}:${DATE:10:2}
  #cat $DATA_DIR/upaobs/OBS:${DATE:0:10} >> $DATA_DIR/littler/little_r:${DATE:0:4}-${DATE:4:2}-${DATE:6:2}_${DATE:8:2}:${DATE:10:2}
  rm -f datelist
  time_lag=$conv_frequency # +/- hours to grab observations from other files (should be ~ obs_interval)
  obs_interval=$conv_frequency #observation availability interval
  for offset in `seq $((OWMIN/60-$time_lag)) $obs_interval $((OWMAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done

  for d in `cat datelist |sort |uniq`; do
    #cat $DATA_DIR/littler/little_r:${d:0:4}-${d:4:2}-${d:6:2}_${d:8:2}:${d:10:2} >> obs.raw
    #cat $DATA_DIR/upaobs/OBS:${DATE:0:10} >> obs.raw
    #cat $DATA_DIR/sfcobs/SURFACE_OBS:${DATE:0:10} >> obs.raw
    #if [ -e $DATA_DIR/conventional/obs.${d:0:10} ]; then
    #  cat $DATA_DIR/conventional/obs.${d:0:10} >> obs.raw
    #fi
    if [ -e $DATA_DIR/conventional/SURFACE_OBS:${d:0:10} ]; then
      cat $DATA_DIR/conventional/SURFACE_OBS:${d:0:10} >> obs.raw
    fi
    if [ -e $DATA_DIR/conventional/OBS:${d:0:10} ]; then
      cat $DATA_DIR/conventional/OBS:${d:0:10} >> obs.raw
    fi
  done

fi

##### include MADIS data (hourly) data ######
if $INCLUDE_MADIS; then
  rm -f datelist
  time_lag=1
  obs_interval=1
  for offset in `seq $((OWMIN/60-$time_lag)) $obs_interval $((OWMAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done
  for d in `cat datelist |sort |uniq`; do
    cat $DATA_DIR/madis/`echo $d |cut -c1-6`/madis_`echo $d |cut -c1-10`_littler >> obs.raw
  done
fi

##### include BUFR ADP Surface and Upperair (6-hourly) data #####
if $INCLUDE_BUFR; then
  bufr_sfc_decode_dir=$CODE_DIR/BUFR/bufr_decode_ADPsfc_littler/exe
  bufr_upa_decode_dir=$CODE_DIR/BUFR/bufr_decode_ADPupa_littler/exe
  bufr_dir=$DATA_DIR/bufr
  rm -f datelist
  time_lag=6
  obs_interval=6
  for offset in `seq $((OWMIN/60-$time_lag)) $obs_interval $((OWMAX/60+$time_lag))`; do
    obsdate=`advance_time $DATE $((offset*60))`
    hh=`echo $obsdate |cut -c9-10`
    inc=`echo $hh%$obs_interval*60 |bc`
    if [[ $inc -lt $((obs_interval*60/2)) ]]; then
      obsdate=`advance_time $obsdate -$inc`
    else
      obsdate=`advance_time $obsdate $((obs_interval*60-inc))`
    fi
    echo $obsdate >> datelist
  done
  for d in `cat datelist |sort |uniq`; do
    ccyymmdd=`echo $d |cut -c1-8`
    hh=`echo $d |cut -c9-10`
    $bufr_sfc_decode_dir/bufr_sfc2ob.x $bufr_dir/sfcobs.$ccyymmdd/gdas.adpsfc.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_sfc_decode_dir/bufr_ship2ob.x $bufr_dir/sfcobs.$ccyymmdd/gdas.sfcshp.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_upa2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.adpupa.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_aircar2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.aircar.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_craft2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.aircft.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    $bufr_upa_decode_dir/bufr_sat2ob.x $bufr_dir/upaobs.$ccyymmdd/gdas.satwnd.t${hh}z.$ccyymmdd.bufr $ccyymmdd$hh >> decode_bufr.log 2>&1
    rm -f files.txt
    for type in Airca Aircraft Satob Ship Surface Upper; do 
      echo $type$ccyymmdd$hh.obs >> files.txt
    done
    $bufr_sfc_decode_dir/runob2lit_imd_obs.x files.txt $ccyymmdd$hh >> decode_bufr.log 2>&1
    cat OBS:$ccyymmdd$hh >> obs.raw
    cat SURFACE_OBS:$ccyymmdd$hh >> obs.raw
  done
  rm -f OBS:* SURFACE_OBS:* *.obs
fi

#AIRBORNE data
if ( $INCLUDE_AIRBORNE_RV ); then
  RADAR_WINDOW=30 #minutes
  echo > $WORK_DIR/obs/$DATE/airborne_`echo $DATE |cut -c1-12`_so
  for offset in `seq -$RADAR_WINDOW $RADAR_WINDOW`; do
    obsdate=`advance_time $DATE $offset`
    obsfile=$SO_DIR/`echo $obsdate |cut -c1-12`_all.so_ass
    #echo $obsfile
    if [ -f $obsfile ]; then
      cat $obsfile >> $WORK_DIR/obs/$DATE/airborne_`echo $DATE |cut -c1-12`_so .
    fi
  done
fi

#COPLANE data
if ( $USE_COPLANE_UV ); then
  RADAR_WINDOW=30 #minutes
  for offset in `seq -$RADAR_WINDOW $RADAR_WINDOW`; do
    obsdate=`advance_time $DATE $offset`
    obsfile=$DATA_DIR/coplane/`echo $DATE |cut -c1-4`/coplane_`echo $obsdate |cut -c1-12`
    if [ -f $obsfile ]; then
      cat $obsfile >> $WORK_DIR/obs/$DATE/coplane_`echo $DATE |cut -c1-12`
    fi
  done
fi

#####  START OBSPROC #####
for var_type in 3DVAR 4DVAR; do
  case $var_type in 
    3DVAR)
      if ! $RUN_ENKF; then continue; fi
    ;;
    4DVAR)
      if ! $RUN_4DVAR; then continue; fi
    ;;
  esac
  echo > obsproc.log
  export use_for=$var_type
  $SCRIPT_DIR/namelist_obsproc.sh > namelist.obsproc
  ./obsproc.exe >& obsproc.log

  watch_log obsproc.log 99999 1 $rundir
  if $INCLUDE_TCV; then
    ln -sf $TCV_DIR/HPI_obs_gts_`wrf_time_string $DATE`.3DVAR TCV_3dvar_${DATE}00
    #ln -sf $TCV_DIR/HPI_EnKF_fmt_`wrf_time_string $DATE` TCV_3dvar_${DATE}00
    if [ -f obs_gts_*.$var_type ]; then
      cp obs_gts_*.$var_type Conventional_3dvar_${DATE}00
      Conv_tot=$(head -n 1 Conventional_3dvar_${DATE}00 | cut -c 8-14)
    else
     Conv_tot=0
    fi
    TCV_tot=$(head -n 1 TCV_3dvar_${DATE}00 | cut -c 8-14)
    tot_obs_num=$(echo "$TCV_tot+$Conv_tot" | bc)
    tot_obs_num=`printf "%.*i" 5 $tot_obs_num`
    $SCRIPT_DIR/gts_header_template.sh $tot_obs_num > obs_gts_`wrf_time_string $DATE`.$var_type
    cat TCV_3dvar_${DATE}00 | tail -n+22 >> obs_gts_*.$var_type
    if [ -f Conventional_3dvar_${DATE}00 ]; then
      cat Conventional_3dvar_${DATE}00 | tail -n+22 >> obs_gts_*.$var_type
    fi
    cp obs_gts_*.$var_type $WORK_DIR/obs/$DATE/.
  else
    cp obs_gts_*.$var_type $WORK_DIR/obs/$DATE/.
  fi
done

if $CLEAN; then rm obs.raw; fi

echo complete > stat

