#!/bin/bash
#This function generates domain_moves for namelist.input so that the nested domains follow
#  storm center (from tcvitals data). Storm center lat/lon are calculated at each time (including
#  start/end time for wrf and any lbc (0 6 12 18Z) time in between. The storm is assumed to move
#  linearly between any two given location points.
#      
#  Only move the smallest domain, and use corral_dist to force move its parent domains
#
#Input: date1 = wrf run start date
#       date2 = wrf run end date
#Main Output: a text file that contains namelist options related to domain moves.

. $CONFIG_FILE

date1=$1 
echo "Current cycle is at" $date1
date2=$2
echo "Next cycle is at" $date2
outfile=$3
echo "Will output the domain_moves to" $outfile

echo " "
echo "--------------------------------------------------------------------------------------"
echo "Find TCvitals data that covers the period between the current cycle and the next cycle"
echo "--------------------------------------------------------------------------------------"
echo "Minutes between two consecutive TCvital records are" $TCV_INTERVAL
echo "Using the date of current cycle, next cycle, and the interval of TCvital records to calculate the off time......"
off1=`echo "(${date1:8:2}*60+${date1:10:2})%$TCV_INTERVAL" |bc` # string:P:L (P indicates the starting index, L is the length of the substring)
echo "Current cycle is" $off1 "minutes off (off1) from the nearest TCvital record"
off2=`echo "(${date2:8:2}*60+${date2:10:2})%$TCV_INTERVAL" |bc`
echo "Next cycle is" $off2 "minutes off (off2) from the nearest TCvital record"
fgdate1=`advance_time $date1 -$off1`
echo "The start TCvital record (fgdate1) is at" $fgdate1
if [ `echo "${date2:8:2}%6" |bc` -eq "0" ]; then
    fgdate2=$date2
else
    fgdate2=`advance_time $date2 $((TCV_INTERVAL-$off2))`
fi
echo "The end TCvital record (fgdate2) is at" $fgdate2
n=$(echo "`diff_time $fgdate1 $fgdate2`/$TCV_INTERVAL" |bc)
echo "The period between the start and the end TCvital record is" $n "*(360) minutes"
off3=`diff_time $fgdate1 $date2`
echo "Next cycle is" $off3 "minutes off (off3) from the start TCvital record"

echo " "
echo "--------------------------------------------------------------------------------------"
echo "Read lat and lon from TCvitals "
echo "--------------------------------------------------------------------------------------"
# Ge_time $date1 $((TCV_INTERVAL-$off1))` lat lon at each time node
ndate=()
nlat=()
nlon=()
for i in `seq 0 $n`; do
  tcdate=`advance_time $fgdate1 $((i*$TCV_INTERVAL))`
  ndate+=($tcdate)
  echo "Read" $(echo "$i+1" |bc) "th TCvital record at" $tcdate
  tcvitals_data=$TCVITALS_DIR/${tcdate:0:4}/${tcdate}.${STORM_ID}-tcvitals.dat
  if [ ! -f $tcvitals_data ]; then echo "$tcvitals_data not found!"; exit; fi
  latstr=`head -n1 $tcvitals_data |awk '{print $6}'`
  lonstr=`head -n1 $tcvitals_data |awk '{print $7}'`
  if [ ${latstr:0-1} == "N" ]; then
    nlat+=($(echo $(echo ${latstr} | awk '{ print substr( $0, 1, length($0)-1 ) }')/10 |bc -l))
  else
    nlat+=($(echo $(echo -${latstr} | awk '{ print substr( $0, 1, length($0)-1 ) }')/10 |bc -l))
  fi
  if [ ${lonstr:0-1} == "E" ]; then
    nlon+=($(echo $(echo ${lonstr} | awk '{ print substr( $0, 1, length($0)-1 ) }')/10 |bc -l))
  else
    nlon+=($(echo $(echo -${lonstr} | awk '{ print substr( $0, 1, length($0)-1 ) }')/10 |bc -l))
  fi
done
echo "Times are at ${ndate[@]}"
echo "Lats are at" ${nlat[@]}
echo "Lons are at" ${nlon[@]}

echo " "
echo "--------------------------------------------------------------------------------------"
echo "Convert the geolocation of TCvitals to index space (ij) in domain 1 (D1)"
echo "--------------------------------------------------------------------------------------"
len_data=${#nlat[@]}
end_index=$(echo $(($len_data-1)))

# Get domain d01 i/j at storm center each time
for i in $(seq 0 $end_index); do
  echo "Process" $(echo "$i+1" |bc) "th TCvital record..." 
  cat << EOF > ll2ij.ncl
begin
   opt = True
   opt@MAP_PROJ          = 3
   opt@TRUELAT1          = $TRUELAT1
   opt@TRUELAT2          = $TRUELAT2
   opt@STAND_LON         = $STAND_LON
   opt@REF_LAT           = $REF_LAT
   opt@REF_LON           = $REF_LON
   opt@KNOWNI            = `echo "${E_WE[0]}/2" |bc -l`
   opt@KNOWNJ            = `echo "${E_SN[0]}/2" |bc -l`
   opt@DX                = ${DX[0]} 
   opt@DY                = ${DY[0]}
   loc=wrf_ll_to_ij(${nlon[$i]},${nlat[$i]},opt)
   asciiwrite("ij",round(loc,3))
end
EOF
  module restore default
  ncl ll2ij.ncl
  ni[$i]=`cat ij |head -n1`
  nj[$i]=`cat ij |tail -n1`
  echo "ni->`cat ij |head -n1`"
  echo "nj->`cat ij |tail -n1`"
  rm -f ll2ij.ncl ij
  module restore intel
done

echo " "
echo "--------------------------------------------------------------------------------------"
echo "Generate domain_moves options"
echo "--------------------------------------------------------------------------------------"
num_moves=0
dt_from1stVital=0
move_id=""
move_interval=""
move_cd_x=""
move_cd_y=""
echo "Current cycle is" $off1 "minutes off (off1) from the start TCvital record"
echo "Next cycle is" $off3 "minutes off (off3) from the start TCvital record"
echo " "
echo "Loop through each needed TCvitals..."
for i in `seq 0 $((end_index-1))`; do  
  echo '      '
  echo "Evolution of TCvital locations in index space in D1 is:" "${ndate[$i]} (${ni[$i]} ${nj[$i]}) -> ${ndate[$i+1]} (${ni[$i+1]} ${nj[$i+1]})"
  dmin_vitals=`diff_time ${ndate[$i]} ${ndate[$i+1]}`
  echo "The next TCvital record is" $dmin_vitals "minutes apart."
  echo "From this time to next time..."
  di_vitals=`echo "(${ni[$i+1]} - ${ni[$i]})" |bc`
  echo "i direction: two locations are off" $di_vitals "grid points in D1"
  dj_vitals=`echo "(${nj[$i+1]} - ${nj[$i]})" |bc`
  echo "j direction: two locations are off" $dj_vitals "grid points in D1"
  for i in `seq 2 $((MAX_DOM-1))`; do #Convert from D01 to smallest domain's parent domain for move steps
    echo "For D$i..."
    di_vitals=`echo "$di_vitals*${GRID_RATIO[$i-1]}" |bc`
    echo "i direction: two locations are off" $di_vitals "grid points in D$i" 
    dj_vitals=`echo "$dj_vitals*${GRID_RATIO[$i-1]}" |bc`
    echo "j direction: two locations are off" $dj_vitals "grid points in D$i"
  done
  d_ij_vitals=$(echo "`abs $di_vitals` + `abs $dj_vitals`" |bc)
  echo "Two locations are off in absoluate value:" $d_ij_vitals "grid points in D$i" 
  
  echo '-------------------------------------------------'
  if [ $d_ij_vitals -gt 0 ]; then  # If displacement is non-zero, move the domain
    echo "The displacement is not zero, moving the domain......"
    echo "For" $dmin_vitals "minutes, the storm moves" $d_ij_vitals "steps. Every step takes" $(echo $(($dmin_vitals/$d_ij_vitals))) "minutes"
    echo "Note: in the current algorithm, each step represents ONE grid point!! " 
    for step in `seq 1 $d_ij_vitals`; do
      echo '       '
      echo "Setting up the" $step "th step..."
      if [ $step -lt $d_ij_vitals ]; then #break up move steps
        dt_from1stVital=`echo "$dt_from1stVital + $dmin_vitals/$d_ij_vitals" |bc`
        echo $dt_from1stVital "minutes ahead from the start TCvital record"
      else
        dt_from1stVital=`echo "$dt_from1stVital + $dmin_vitals - ($d_ij_vitals-1)*$dmin_vitals/$d_ij_vitals" |bc`
        echo $dt_from1stVital "minutes ahead from the start TCvital record"
      fi

      if [ $off1 -lt $dt_from1stVital ] && [ $dt_from1stVital -le $off3 ]; then
        echo "The time associated with the step is ahead of the current cycle and behind the next cycle. Calculate the preset move....."
        num_moves=$((num_moves+1))
        move_id="$move_id $MAX_DOM,"
        dummy_step=`echo "($TIME_STEP + 59)/60" |bc` #Adjust for timestep (Division only keeps the integer part eg., 10/3 = 3)
        if [ $dt_from1stVital -ge `expr $off3 - $dummy_step` ]; then
         echo "The step is very close to the next cycle! Customizing the move interval..." 
          interval=$(echo $(($dt_from1stVital - $off1 - $dummy_step)))
          echo "The step is $interval minutes in advance from the current cycle"
          move_interval="$move_interval $interval,"
        else
          interval=$(echo $(($dt_from1stVital - $off1)))
          echo "The step is $interval minutes in advance from the current cycle"
          move_interval="$move_interval $interval,"
        fi

        if [ `abs $di_vitals` -gt 0 ] && [ `abs $dj_vitals` -gt 0 ]; then
          echo "For this step, the move happens in both i and j direction between the two TCvital records"
          if [ `abs $di_vitals` -ge `abs $dj_vitals` ]; then
            echo "The moving distance in i direction is larger in j direction"
            ds=$(echo "$d_ij_vitals/`abs $dj_vitals`" |bc)
            if [ `echo "$step % $ds" |bc` == 0 ] && [ $step -le $(echo "`abs $dj_vitals` * $ds" |bc) ]; then
              move_cd_x="$move_cd_x 0,"
              move_cd_y="$move_cd_y $(echo "$dj_vitals/`abs $dj_vitals`" |bc),"
              echo 'Move one step/grid in j direction.'
            else
              move_cd_x="$move_cd_x $(echo "$di_vitals/`abs $di_vitals`" |bc),"
              move_cd_y="$move_cd_y 0,"
              echo 'Move one step/grid in i direction.'
            fi
          else
            echo "The moving distance is further in j direction thatn in i direction"
            ds=$(echo "$d_ij_vitals/`abs $di_vitals`" |bc)
            if [ `echo "$step % $ds" |bc` == 0 ] && [ $step -le $(echo "`abs $di_vitals` * $ds" |bc) ]; then
              move_cd_x="$move_cd_x $(echo "$di_vitals/`abs $di_vitals`" |bc),"
              move_cd_y="$move_cd_y 0,"
              echo 'Move one step/grid in i direction.'
            else
              move_cd_x="$move_cd_x 0,"
              move_cd_y="$move_cd_y $(echo "$dj_vitals/`abs $dj_vitals`" |bc),"
              echo 'Move one step/grid in j direction.'
            fi
          fi
        elif [ `abs $di_vitals` -gt 0 ] && [ `abs $dj_vitals` -eq 0 ]; then
          echo "The move happens only in i direction between the two TCvital records"
          move_cd_x="$move_cd_x $(echo "$di_vitals/`abs $di_vitals`" |bc),"
          move_cd_y="$move_cd_y 0,"
          echo 'Move one step/grid in i direction.'
        else
          echo "The move happens only in j direction between the two TCvital records"
          move_cd_x="$move_cd_x 0,"
          move_cd_y="$move_cd_y $(echo "$dj_vitals/`abs $dj_vitals`" |bc),"
          echo 'Move one step/grid in j direction.'
        fi
      fi
    done
  else
    dt_from1stVital=`echo "$dt_from1stVital + $dmin_vitals" |bc`
  fi
done

# Sanity Check 
echo " "
echo "--------------------------------------------------------------------------------------"
echo "Sanity check the final location"
echo "--------------------------------------------------------------------------------------"
# split string and store it into an array
IFS=', ' read -ra move_x_arr <<< "$move_cd_x"
IFS=', ' read -ra move_y_arr <<< "$move_cd_y"
# accumulate the steps
total_xmove_D2=$(IFS=+; echo "$((${move_x_arr[*]}))")
total_ymove_D2=$(IFS=+; echo "$((${move_y_arr[*]}))")
total_xmove_D1=$(echo "$total_xmove_D2/${GRID_RATIO[1]}" |bc)
total_ymove_D1=$(echo "$total_ymove_D2/${GRID_RATIO[1]}" |bc)
# locate the final location in ij space in D1
i_final_D1=$(echo "${ni[0]}+$total_xmove_D1" | bc)
j_final_D1=$(echo "${nj[0]}+$total_ymove_D1" | bc)
echo "Final calculated location in the index space for D1 is" $i_final_D1 $j_final_D1
# Get lat lon of the final location 
cat << EOF > ij2ll_final.ncl
begin
   opt = True
   opt@MAP_PROJ          = 3
   opt@TRUELAT1          = $TRUELAT1
   opt@TRUELAT2          = $TRUELAT2
   opt@STAND_LON         = $STAND_LON
   opt@REF_LAT           = $REF_LAT
   opt@REF_LON           = $REF_LON
   opt@KNOWNI            = `echo "${E_WE[0]}/2" |bc -l`
   opt@KNOWNJ            = `echo "${E_SN[0]}/2" |bc -l`
   opt@DX                = ${DX[0]} 
   opt@DY                = ${DY[0]}
   loc=wrf_ij_to_ll(${i_final_D1},${j_final_D1},opt)
   asciiwrite("ll_final",round(loc,3))
end
EOF
module restore default
ncl ij2ll_final.ncl
module restore intel
echo "Final lon->`cat ll_final |head -n1`"
echo "Final lat->`cat ll_final |tail -n1`"

echo " "
echo "--------------------------------------------------------------------------------------"
echo "Writing output..."
echo "--------------------------------------------------------------------------------------"
corral_dist="8, 8,"
for i in `seq 3 $MAX_DOM`; do
  corral_dist="$corral_dist $(echo "`min ${I_PARENT_START[$i-1]} ${J_PARENT_START[$i-1]}`-2" |bc),"
  echo "The distance that the moving nest (D3) is allowed to get near the mother domain (D2)  boundary: "  ${corral_dist[@]}
done

# Output
echo "num_moves     = $num_moves," > $outfile
if [ $num_moves -gt 0 ]; then
  echo "move_id       = $move_id" >> $outfile
  echo "move_interval = $move_interval" >> $outfile
  echo "move_cd_x     = $move_cd_x" >> $outfile
  echo "move_cd_y     = $move_cd_y" >> $outfile
fi
echo "corral_dist   = $corral_dist" >> $outfile
