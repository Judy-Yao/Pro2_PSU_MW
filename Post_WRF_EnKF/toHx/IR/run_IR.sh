#!/bin/bash --login

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 16
#SBATCH --ntasks-per-node 64
#SBATCH -p development
#SBATCH -t 2:00:00
#SBATCH -o out_enkf
#SBATCH -e error_enkf

module restore crtm_mw
for file in `ls /scratch/02191/yuz31/HARVEY/MW_THO/output/*/ATCF`; do
  cd `dirname $file`
  ln -sf /work/02191/yuz31/stampede2/code/CRTM/forward_rt_simple/newrt/paper_scott_2020.mpi .
  ln -sf /scratch/02191/yuz31/coefficients_l coefficients
  ln -sf wrfout_d03_2017-08-25_12:00:00 wrffile
  cat > test_crtm_wrf.nml << EOF
&rt_settings
  nml_s_rt_program = 'crtm'
  nml_l_include_land = T,
  nml_s_sensor_id='abi_gr'
  nml_a_channels(:,1) = 8,9,10
!  nml_a_channels(:,2) = 13,
!  nml_s_reff_method='default',
  nml_s_reff_method='mp_physics',
!  nml_i_x_calc_beg = 140,
!  nml_i_x_calc_end = 160,
!  nml_i_y_calc_beg = 140,
!  nml_i_y_calc_end = 160,
  nml_i_nicpu = 32,
  nml_i_njcpu = 32,
  nml_s_crtm_rainLUT='Thompson08_RainLUT_-109z-1.bin',
  nml_s_crtm_snowLUT='Thompson08_SnowLUT_-109z-1.bin',
  nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1.bin',
!  nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1_nolowfreq.bin',
/

\$rt_input
  nml_s_filename_input = 'wrffile'
  nml_s_filename_obs='/work/02191/yuz31/stampede2/scripts/HARVEY/HARVEY/radiance/radiance_d03_201708251200_so'
/

\$rt_output
  nml_s_filename_output = 'wrffile.tb'
  nml_l_output_reff = F,
/
EOF
  ibrun -n 1024 ./paper_scott_2020.mpi test_crtm_wrf.nml
done
