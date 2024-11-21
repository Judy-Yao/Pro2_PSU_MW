#!/bin/bash --login

module restore intel
export CONFIG_FILE=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/WRF_EnKF/IRMA_run/IR-WSM6_onlyTCvitals/config_IR_WSM6.IRMA.S3
. "$CONFIG_FILE"
module restore intel

domlist=`seq 3 3`
DATE=$1

for n in $domlist; do
  dm=d`expr $n + 100 |cut -c2-`
  for NE in `seq 1 $NUM_ENS`; do
    id=`expr $NE + 1000 |cut -c2-`
    mv $WORK_DIR/run/$DATE/enkf/$dm/fort.`expr 80010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id
    mv $WORK_DIR/run/$DATE/enkf/$dm/fort.`expr 90010 + $NE` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id
    #ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_$id $dm/fort.`expr 80010 + $NE`
    #ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $dm/fort.`expr 90010 + $NE`
    #ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_$id $WORK_DIR/fc/$DATE/wrfinput_${dm}_$id
  done
  cp $WORK_DIR/run/$DATE/enkf/$dm/fort.`expr 80011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_input_${dm}_mean
  cp $WORK_DIR/run/$DATE/enkf/$dm/fort.`expr 90011 + $NUM_ENS` $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean
  #ln -fs $WORK_DIR/fc/$DATE/wrf_enkf_output_${dm}_mean $WORK_DIR/fc/$DATE/wrfinput_${dm}
done
