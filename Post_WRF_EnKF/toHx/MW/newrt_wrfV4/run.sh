#!/bin/bash --login

#####header for stampede######
#SBATCH -J harvey_robert
#SBATCH -N 4
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 2:00:00
#SBATCH -o out_enkf
#SBATCH -e error_enkf

#ibrun -n 48 ./test.exe  ../test/test_rttov_mie_hrrr.nml
#ibrun -n 48 ./test.exe  ../test/test_rttov_liu09_hrrr.nml
#ibrun -n 48 ./test.exe  ../test/test_rttov_kuoall_hrrr.nml
#ibrun -n 48 ./test.exe  ../test/test_rttov_kuofit_hrrr.nml
#ibrun -n 48 ./test.exe  ../test/test_rttov_kuohyper_hrrr.nml

YYYY=2020
mm=04


offset=0
#for d in {12..13}; do
for d in {12..14}; do
  dd=`echo $((d+1000)) | cut -c 3-4`

  for H in {0..23}; do 
    HH=`echo $((H+1000)) | cut -c 3-4`

    output_dir=/scratch/05012/tg843115/paper_mw_scatt/$YYYY/$mm/$dd

    if [ ! -d $output_dir ]; then
      mkdir -p $output_dir
    fi

    for LUT in liu09 mie kuoall kuofit kuohyper; do
      #for fcsthour in 0 ; do
      for fcsthour in 0 1 6; do
        FH=`echo $((fcsthour+1000)) | cut -c 3-4`

        args1="nml_s_filename_input=/home1/05012/tg843115/wave/data2/RAP/rap/$YYYY$mm$dd/rap.t${HH}z.wrfnatf${FH}.grib2"

        FHTIME=`date -d "$fcsthour hours $YYYY$mm$dd ${HH}:00" +%Y%m%d%H%M`

        YYYY2=`echo $FHTIME | cut -c1-4`
        mm2=`echo $FHTIME | cut -c5-6`
        dd2=`echo $FHTIME | cut -c7-8`
        HH2=`echo $FHTIME | cut -c9-10`
        args2="nml_s_filename_obs=/home1/05012/tg843115/wave/data2/GPM/GMI/ascii/$YYYY2/$mm2/test_$YYYY2$mm2$dd2${HH2}00_gmi_ch01x063x063_ch03x036x036_ch05x036x036_ch06x036x036_ch08x036x036_ch10x018x018_ch12x018x018_ch13x018x018.txt"

        output_prefix=${output_dir}/tb.rap.$YYYY$mm${dd}${HH}00.f${FH}.${LUT}
        args3="nml_s_filename_output=$output_prefix"
        args4="nml_s_rttov_mietable=/work/05012/tg843115/stampede2/opt/RTTOV/rttov_test_coef/rtcoef_rttov12/test_mietable/${LUT}/mietable_gpm_gmi.dat"
        outlog=${output_prefix}.out
        errlog=${output_prefix}.err
        #echo $outlog
        if [[ -f $outlog && `tail -n 2 $outlog | grep -c SUCCESS` == 1 ]]; then
          continue
        else
          echo $outlog $offset
          echo $errlog $offset
          ibrun -o $offset -n 48 ./paper_mw_scatt.mpi ../test/paper_mw_scatt.nml $args1 $args2 $args3 $args4 > ${output_prefix}.out 2>${output_prefix}.err &
          offset=$((offset + 48))
        fi
        #sleep 1

        if [[ $offset -ge 192 ]]; then
          wait
          offset=0
        fi
      done # fcsthour
      #wait
    done # LUT
  done
done
