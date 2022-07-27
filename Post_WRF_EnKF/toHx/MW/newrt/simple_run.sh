#!/bin/bash

if (( $# == 2 )); then
fcst_start_time=$1
fcst_evalu_time=$2
else
fcst_start_time=201708251200
fcst_evalu_time=201708252200
fi
echo $fcst_start_time
echo $fcst_evalu_time

nml='/home1/05012/tg843115/src/forward_rt/test/scott_2020.nml'
exe='/home1/05012/tg843115/src/forward_rt/newrt/paper_scott_2020.mpi'

#dirname=/scratch/05012/tg843115/Harvey_run_scott_nograup/Harvey_conv_IR_2212_abei
dirname=/scratch/05012/tg843115/EnKF_crtm_new_run/IR_MW_noShips_2212_DM1_graupelNoLowNoUpdate


dir2=$dirname/output/$fcst_start_time/unthinned/$fcst_evalu_time
echo $dir2

#if [[ -f $dir2/crtm.log && `tail -n 3 $dir2/crtm.log | grep -c SUCCESS` -gt 0 ]]; then
#  echo done
#else
  ibrun -n 48 $exe $nml "nml_s_filename_input=$dir2/fort.80071" "nml_s_filename_output=$dir2/tb.80071" "nml_s_filename_obs=$dir2/microwave_${fcst_evalu_time}_so" > $dir2/crtm.log
#fi

dirname=/scratch/05012/tg843115/Harvey_run_scott_nograup/Harvey_conv_GPM-03x36x24-13x36x24_IR_2212_abei
dir2=$dirname/output/$fcst_start_time/unthinned/$fcst_evalu_time
echo $dir2
#if [[ -f $dir2/crtm.log && `tail -n 3 $dir2/crtm.log | grep -c SUCCESS` -gt 0 ]]; then
#  echo done
#else
#  ibrun -n 48 $exe $nml "nml_s_filename_input=$dir2/fort.80071" "nml_s_filename_output=$dir2/tb.80071" "nml_s_filename_obs=$dir2/microwave_${fcst_evalu_time}_so" > $dir2/crtm.log
#fi


