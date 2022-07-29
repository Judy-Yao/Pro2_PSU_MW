# Ongoing work: to automatically check if ij_parent_start in one file for all of members are the same
if [[ $DATE == $DATE_START ]]; then
  spinup_end=$(advance_time $DATE $CP) 
  Enddate=$(wrf_time_string $spinup_end)
	
  domlist=$(seq 1 $MAX_DOM)
  for n in $domlist; do
    dm=`expr $n + 100 |cut -c2-`
	for file in $(ls ???/wrfout_${dm}*${Enddate}); do
	  # Dump the line where I_PARENT_START is
	  ips=$(ncdump -h $file | grep I_PARENT_START)
	  # Remove any (beginning, tailing, intermediate) white spaces. Element in array is separated by white space
	  ips_nospace=$(echo -e "${ips}" | tr -d '[:space:]')
	  # Add new element at the end of the array
	  arrVar+=($ips_nospace)
	done	
	# Count the unique item
	n_uniq=$(echo $(tr ' ' '\n' <<<"${arrVar[@]}" | awk '!u[$0]++') | wc -l)
		# use  tr to translate the spaces between  items into newlines (\n):break the array apart into separate lines
		# pipe this to awk, and have it only quote unique lines
		# count the number of lines
	# 
	if [ $n_uniq -eq 1 ] 
	  # continue	
	else
	  #
	fi



  done

else


fi

