#!/bin/bash
#This function generates domain_moves for namelist.input so that the nested domains follows 
#  storm center (from tcvitals data). Storm center lat/lon are calculated at each time (including
#  start/end time for wrf and any lbc (0 6 12 18Z) time in between. The storm is assumed to move
#  linearly between any two given location points.
#      
#  Only move the smallest domain, and use corral_dist to force move its parent domains
#
#Input: date1 = wrf run start date
#       date2 = wrf run end date
#Output: a text file that contains namelist options related to domain moves.

. util.sh
. $CONFIG_FILE

date1=201708221300
date2=201708221400
#date2=$DATE_END
outfile=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/SourceCode/WRF_EnKF/MARIA_run/IR-J_DA+J_WRF+J_init-SP-intel17/debug/test

echo "---------------------------------------------------"
echo "Find TCvitals data that covers the period"
echo "---------------------------------------------------"
# Find tcvitals data that covers the period
echo "Current cycle is at" $date1
echo "Next cycle is at" $date2
echo "Minutes between two consecutive TCvital records are" $TCV_INTERVAL
off1=`echo "(${date1:8:2}*60+${date1:10:2})%$TCV_INTERVAL" |bc` # string:P:L (P indicates the starting index, L is the length of the substring)
echo "Current cycle is" $off1 "minutes off (off1) from the nearest TCvital record"
off2=`echo "(${date2:8:2}*60+${date2:10:2})%$TCV_INTERVAL" |bc`
echo "Next cycle is" $off2 "minutes off (off2) from the nearest TCvital record"
echo "off1->$off1"
echo "off2->$off2"
fgdate1=`advance_time $date1 -$off1`
echo "The start TCvital record (fgdate1) is at" $fgdate1
fgdate2=`advance_time $date2 $((TCV_INTERVAL-$off2))`
echo "The end TCvital record (fgdate2) is at" $fgdate2
echo "fgdate1->$fgdate1"
echo "fgdate2->$fgdate2"
n=$(echo "`diff_time $fgdate1 $fgdate2`/$TCV_INTERVAL" |bc)
echo "The period between the start and end TCvital records is" $n "(360) minutes"
echo "n->$n"
off3=`diff_time $fgdate1 $date2`
echo "Next cycle is" $off3 "minutes off (off3) from the start TCvital record"
echo "off3->$off3"

echo "---------------------------------------------------"
echo "Read lat and lon from TCvitals "
echo "---------------------------------------------------"
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
echo "Times are ${ndate[@]}"
echo "Lats are" ${nlat[@]}
echo "Lons are" ${nlon[@]}

echo "---------------------------------------------------"
echo "Convert the geolocation of TCvitals to index space (ij) in domain 1 (D1)"
echo "---------------------------------------------------"
# Interpolate if time node is not at lbc time
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
  ncl ll2ij.ncl
  ni[$i]=`cat ij |head -n1`
  nj[$i]=`cat ij |tail -n1`
  echo "ni->`cat ij |head -n1`"
  echo "nj->`cat ij |tail -n1`"
  rm -f ll2ij.ncl ij
done


echo "---------------------------------------------------"
echo "Generate domain_moves options"
echo "---------------------------------------------------"
num_moves=0
dt=0
move_id=""
move_interval=""
move_cd_x=""
move_cd_y=""
echo "Current cycle is" $off1 "minutes off (off1) from the start TCvital record"
echo "Next cycle is" $off3 "minutes off (off3) from the start TCvital record"

for i in `seq 0 $((end_index-1))`; do  #For each time
  echo "Location of storm (from TCvitals) in index space in D1 is:" "${ndate[$i]} (${ni[$i]} ${nj[$i]}) -> ${ndate[$i+1]} (${ni[$i+1]} ${nj[$i+1]})"
  dmin=`diff_time ${ndate[$i]} ${ndate[$i+1]}`
  echo "minutes between two records are" $dmin
  di=`echo "(${ni[$i+1]} - ${ni[$i]})" |bc`
  echo "i direction: two locations are off" $di "grid points in D1"
  dj=`echo "(${nj[$i+1]} - ${nj[$i]})" |bc`
  echo "j direction: two locations are off" $dj "grid points in D1"
  for i in `seq 2 $((MAX_DOM-1))`; do #Convert from D01 to smallest domain's parent domain for move steps
    echo "For D" $i "..."
    di=`echo "$di*${GRID_RATIO[$i-1]}" |bc`
    echo "i direction: two locations are off" $di "grid points in D$i" 
    dj=`echo "$dj*${GRID_RATIO[$i-1]}" |bc`
    echo "j direction: two locations are off" $dj "grid points in D$i"
  done
  s=$(echo "`abs $di` + `abs $dj`" |bc)
  echo "Two locations are off in absoluate value:" $s "grid points in D$i" 

  if [ $s -gt 0 ]; then  # If displacement is non-zero, move the domain
    echo "The displacement is not zero, moving the domain......"
    echo "For" $dmin "minutes, the storm moves" $s "steps. Every step takes" $(echo $(($dmin/$s))) "minutes"
    for t in `seq 1 $s`; do
      echo "Moving" $t "th step..."
      if [ $t -lt $s ]; then #break up move steps
        dt=`echo "$dt + $dmin/$s" |bc`
        echo "Timing:" $dt "minutes (from the start storm location/TCvital records)"
      else
        dt=`echo "$dt + $dmin - ($s-1)*$dmin/$s" |bc`
        echo "Timing:" $dt "minutes (from the start storm location/TCvital records)"
      fi
      if [ $off1 -lt $dt ] && [ $dt -le $off3 ]; then
        echo "The storm (from TCvitals) is moving between the current cycle and the next cycle..."
        num_moves=$((num_moves+1))
        move_id="$move_id $MAX_DOM,"
        dummy_step=`echo "($TIME_STEP + 59)/60" |bc` #Adjust for timestep (Division only keeps the integer part eg., 10/3 = 3)
        if [ $dt -ge `expr $off3 - $dummy_step` ]; then #if domain moves after last time step move earlier 
         echo "dummy_step->$dummy_step"
          interval=$(echo $(($dt - $off1 - $dummy_step)))
          echo "Storm is $interval minutes in advance from the current cycle"
          move_interval="$move_interval $interval,"
        else
          interval=$(echo $(($dt - $off1)))
          echo "Storm is $interval minutes in advance from the current cycle"
          move_interval="$move_interval $interval,"
        fi
        if [ `abs $di` -gt 0 ] && [ `abs $dj` -gt 0 ]; then
          echo "The storm moves both in i and j direction between the two TCvital records"
          if [ `abs $di` -ge `abs $dj` ]; then
            echo "The storm moves further in i direction thatn in j direction"
            ds=$(echo "$s/`abs $dj`" |bc)
            if [ `echo "$t % $ds" |bc` == 0 ] && [ $t -le $(echo "`abs $dj` * $ds" |bc) ]; then
              move_cd_x="$move_cd_x 0,"
              move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
            else
              move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
              move_cd_y="$move_cd_y 0,"
            fi
          else
            echo "The storm moves further in j direction thatn in i direction"
            ds=$(echo "$s/`abs $di`" |bc)
            if [ `echo "$t % $ds" |bc` == 0 ] && [ $t -le $(echo "`abs $di` * $ds" |bc) ]; then
              move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
              move_cd_y="$move_cd_y 0,"
            else
              move_cd_x="$move_cd_x 0,"
              move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
            fi
          fi
        elif [ `abs $di` -gt 0 ]; then
          echo "The storm only moves in i direction between the two TCvital records"
          move_cd_x="$move_cd_x $(echo "$di/`abs $di`" |bc),"
          move_cd_y="$move_cd_y 0,"
        else
          echo "The storm does not move between the two TCvital records"
          move_cd_x="$move_cd_x 0,"
          move_cd_y="$move_cd_y $(echo "$dj/`abs $dj`" |bc),"
        fi
      fi
    done
  else
    dt=`echo "$dt + $dmin" |bc`
  fi
done

corral_dist="8, 8,"
for i in `seq 3 $MAX_DOM`; do
  corral_dist="$corral_dist $(echo "`min ${I_PARENT_START[$i-1]} ${J_PARENT_START[$i-1]}`-2" |bc),"
  echo "The distance that the moving nest is allowed to get near the mother domain boundary: "  ${corral_dist[@]}
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
