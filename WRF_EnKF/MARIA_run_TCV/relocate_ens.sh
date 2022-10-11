#!/bin/bash --login
rundir=/lfs3/projects/hfip-psu/rnystrom/scratch/Patricia/201510210000_3km_RV_vRTPP0.5/run/201510212100/reloc_Best2
ENKF_DIR=$WORK/WRF_DA/EnKF_fz/src
TCVITALS_DIR=$WORK/data/TCV/raw_interp
ENS_CENTER_DIR=$WORK/data/TCV/Patricia_Ens
ensdir=/lfs3/projects/hfip-psu/rnystrom/scratch/Patricia/201510210000_3km_RV_vRTPP0.5/fc/201510212100
STORM_ID=PATRICIA
DATE=201510212100
MAX_DOM=3
NUM_ENS=60


if [[ ! -d $rundir ]]; then mkdir -p $rundir; echo waiting > $rundir/stat; fi
cd $rundir

### relocate ens members and relax env to whatever
domlist=`seq 1 $MAX_DOM`
for n in $domlist; do
    dm=d`expr $n + 100 |cut -c2-`
    if [[ ! -d $dm ]]; then mkdir -p $dm; fi
    cd $dm
    echo "  Relocating & Replacing environment for domain $dm"
    for NE in `seq 1 $((NUM_ENS))`; do
      id=`expr $NE + 1000 |cut -c2-`
      if [[ ! -d $id ]]; then mkdir $id; fi
      if [ -e ${id}/reloc.log ]; then if [[ `tail -n2 ${id}/reloc.log |grep Successful` ]]; then continue; fi; fi
      cd $id

      ln -sf $ENKF_DIR/relocate_vortex_replace_environment_by_gfs.exe .
      ln -sf $TCVITALS_DIR/${DATE:0:4}/${DATE}.${STORM_ID}-tcvitalsB.dat tcvitals.dat
      ln -sf $ENS_CENTER_DIR/${DATE}_${id}_ens.dat ens.dat
      cp $ensdir/wrf_enkf_output_${dm}_${id} wrfinput
      # link env to relax to 
      ln -sf $ensdir/wrf_enkf_output_${dm}_${id} wrfinput_gfs
      #ln -sf $WORK_DIR/rc/$DATE/wrfinput_${dm} wrfinput_gfs

      ./relocate_vortex_replace_environment_by_gfs.exe >& reloc.log
      wait
      cd ..
    done
    cd ..
done

### 
