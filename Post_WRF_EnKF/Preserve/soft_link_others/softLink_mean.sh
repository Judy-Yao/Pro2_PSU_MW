#! /bin/bash

path_from=/scratch/02191/yuz31/HARVEY/IR_THO/fc/
path_to=/scratch/06191/tg854905/Pro2_PSU_MW/HARVEY/JerryRun/IR_THO/fc/
 
cd $path_from
echo 'Enter ' ${path_from}
for file in $(ls .); do
  echo 'Work at ' $file
  lastwo=${file: -2}
  echo $lastwo
  if [ $lastwo = '00' ]; then
    outdir=${path_to}/${file}/
    if [[ ! -d $outdir ]]; then mkdir -p $outdir; fi 
    ln -s  ${path_from}/${file}/wrf_enkf_input_d03_mean  ${path_to}/${file}/wrf_enkf_input_d03_mean
    ln -s  ${path_from}/${file}/wrf_enkf_output_d03_mean  ${path_to}/${file}/wrf_enkf_output_d03_mean
  else
    echo 'Not a directory!'
  fi 
done
