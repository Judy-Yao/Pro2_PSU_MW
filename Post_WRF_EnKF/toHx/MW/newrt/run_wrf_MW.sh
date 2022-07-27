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

for dir in 202008100800 202008100830 202008100900 202008101100 202008101215 202008101245 202008101345; do
  cd /scratch/02191/yuz31/20200810/MW-5MIN/icbcs_$dir
  ln -sf /work/02191/yuz31/stampede2/code/CRTM/forward_rt_simple/newrt/paper_scott_2020.mpi .
  ln -sf /scratch/02191/yuz31/coefficients_l coefficients
  sensor=`head -1 ../../MW/enkf_${dir}/microwave_${dir}_so | cut -c13-24`
  ch=`head -1 ../../MW/enkf_${dir}/microwave_${dir}_so | cut -c25-36`
  #if [[ $dir == '202008100830' ]]; then sensor='gmi_gpm'; fi
  for ne in {11..50}; do
    file=wrfinput_d01_$ne
    echo $dir/$file
    cat > test_crtm_wrf.nml << EOF
&rt_settings
  nml_s_rt_program = 'crtm'
  nml_l_include_land = T,
  nml_s_sensor_id='`echo $sensor`',
  nml_a_channels(:,1) = `echo $ch`,
!  nml_a_channels(:,2) = 10,
!  nml_s_reff_method='default',
  nml_s_reff_method='mp_physics',
!  nml_i_x_calc_beg = 140,
!  nml_i_x_calc_end = 160,
!  nml_i_y_calc_beg = 140,
!  nml_i_y_calc_end = 160,
  nml_i_nicpu = 64,
  nml_i_njcpu = 16,
  nml_s_crtm_rainLUT='Thompson08_RainLUT_-109z-1.bin',
  nml_s_crtm_snowLUT='Thompson08_SnowLUT_-109z-1.bin',
  nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1.bin',
!  nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1_nolowfreq.bin',
/

\$rt_input
  nml_s_filename_input = '$file'
  nml_s_filename_obs='../../MW/enkf_${dir}/microwave_${dir}_so'
/

\$rt_output
  nml_s_filename_output = '$file.tb'
  nml_l_output_reff = F,
/
EOF
    if [[ ! -f $file.tb.`echo $sensor`.crtm.nc ]]; then # || [[ ! -f $file.tb.abi_g16.crtm.nc ]]; then
      ibrun -n 1024 ./paper_scott_2020.mpi test_crtm_wrf.nml
    fi
  done
done
