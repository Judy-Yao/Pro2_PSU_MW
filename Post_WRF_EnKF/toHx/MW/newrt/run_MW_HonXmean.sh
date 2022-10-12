#!/bin/bash

#####header for stampede######
#SBATCH -J crtm
#SBATCH -N 3
#SBATCH --ntasks-per-node 48
#SBATCH -p skx-dev
#SBATCH -t 02:00:00
#SBATCH -o out_mw
#SBATCH -e error_mw
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yao.zhu.91@gmail.com

File_SS="/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/fc/201708221200/d03/201708221200_sensorCh"

module load intel/18.0.2
module load impi/18.0.2
module load parallel-netcdf/4.3.3.1
module load phdf5/1.8.16
module load libfabric/1.7.0
module load python3/3.7.0 

for DAtime in 201708221200; do
    for file in $(ls /work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/fc/${DAtime}/d03/wrf_enkf_input_d03_001); do

        cd /work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/Obs_Hx/MW/${DAtime}
        ln -sf /home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/MW/newrt/paper_scott_2020.mpi .
        ln -sf /home1/06191/tg854905/Pro2_PSU_MW/Post_WRF_EnKF/toHx/MW/newrt/getSensorInfo.py .
        ln -sf /work2/06191/tg854905/stampede2/opt/CRTM/PSU_EnKF_CRTM/coefficients coefficients
        ln -sf $file wrffile

        # At the DAtime, get sensor and channel
        $(python3 getSensorInfo.py 201708221200)
        File_SS=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/MW_THO/fc/201708221200/d03/${DAtime}_sensorCh
        while IFS= read -r line
        do
            sensor_Chs=($line)

            cat > test_crtm_wrf.nml << EOF
&rt_settings
nml_s_rt_program = 'crtm'
nml_l_include_land = T,
nml_s_sensor_id='`echo ${sensor_Chs[0]}`',
nml_a_channels(:,1) = `echo ${sensor_Chs[@]:1}`,
nml_s_reff_method='mp_physics',
nml_i_nicpu = 12,
nml_i_njcpu = 12,
nml_s_crtm_rainLUT='Thompson08_RainLUT_-109z-1.bin',
nml_s_crtm_snowLUT='Thompson08_SnowLUT_-109z-1.bin',
nml_s_crtm_graupelLUT='Thompson08_GraupelLUT_-109z-1.bin',
/
    
\$rt_input
nml_s_filename_input = 'wrffile'
nml_s_filename_obs='/work2/06191/tg854905/stampede2/Pro2_PSU_MW/HARVEY/Obs_y/MW/microwave_d03_201708221200_so' 
/ 
            
\$rt_output
nml_s_filename_output = 'wrffile.tb'
nml_l_output_reff = F,
/
EOF
            
            ibrun -n 144 ./paper_scott_2020.mpi test_crtm_wrf.nml

        done < "$File_SS"

    done
done





