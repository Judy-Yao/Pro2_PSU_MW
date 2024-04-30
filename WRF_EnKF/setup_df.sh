# Storm
Storm=JOSE
# Experiment
DA=IR
MP=TuneWSM6
Exper=IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-30hr-hroi900
NUM_ENS=60
# Specify time range
spinup_st=201709041200
date_st=201709041200
date_ed=201709060000
# Specify which dir to process
at_run=false
at_wrf_df=true

# wrf_df
if ${at_wrf_df}; then
  # main directory
  cd /scratch/06191/tg854905/Pro2_PSU_MW/${Storm}/${Exper}
  # list initialization times
  times=$(ls wrf_df)
  # enter wrf_df directory
  cd wrf_df; echo 'Enter wrf_df'
  for it in ${times}; do
    cd ${it}; echo 'linking '${it}
    # prepare ic/bc
    ln -s ../../fc/${it}/wrf_enkf_output_d01_mean wrfinput_d01
    ln -s ../../fc/${it}/wrf_enkf_output_d02_mean wrfinput_d02
    ln -s ../../fc/${it}/wrf_enkf_output_d03_mean wrfinput_d03
    ln -s ../../rc/${it}/wrfbdy_d01 .
    # copy namelist.input
    cp ../../run/${it}/wrf_ens/001/namelist.input .
    cd ..
  done
fi


