#!/bin/bash

# Storm
Storm=MARIA
Exper=("J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900" "J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900" "IR-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900" "IR-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900" "IR-TuneWSM6-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900" "IR+MW-J_DA+J_WRF+J_init-SP-intel17-WSM6-24hr-hroi900"  "IR+MW-J_DA+J_WRF+J_init-SP-intel17-THO-24hr-hroi900")

echo 'Starting at ' $(date)

for iexper in "${Exper[@]}"; do
  echo 'Experiment: ' $iexper
  if [[ ! -d /work2/06191/tg854905/stampede2/Pro2_PSU_MW/${Storm}/${iexper}/Data_analyze/run ]]; then mkdir -p /work2/06191/tg854905/stampede2/Pro2_PSU_MW/${Storm}/${iexper}/Data_analyze/run ; fi


  for time in $(ls ${SCRATCH_S2}/Pro2_PSU_MW/${Storm}/${iexper}/run); do
    mkdir /work2/06191/tg854905/stampede2/Pro2_PSU_MW/${Storm}/${iexper}/Data_analyze/run/${time}
    cp -r $SCRATCH_S2/Pro2_PSU_MW/${Storm}/${iexper}/run/${time}/enkf /work2/06191/tg854905/stampede2/Pro2_PSU_MW/${Storm}/${iexper}/Data_analyze/run/${time}/
  done

done
