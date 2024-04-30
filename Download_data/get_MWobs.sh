#!/bin/bash
# Script to download the Microwave Tbs files on Linux for ONE storm 
# Credit to Joseph Chan and Junchang Ju 
# Author: Zhu (Judy) Yao. June 14 - 16, 2022

# ===================================================================
# Configure: carefully evaluate each item and tailor it to your need!
# ===================================================================

# Fill in the storm name
Storm=IRMA

# Specify the project directory
Project_Path=/work2/06191/tg854905/stampede2/Pro2_PSU_MW/Preprocess_Obs/

# Date range of case study (yyyymmddHH)
Timebeg=2017090300
Timend=2017090500

# Evaluate if GPM-calibrated level 1C is used
GPM_1C="true"
## GPM-calibrated sensors
GPM_sensor=("AMSR2" "ATMS" "GMI" "MHS" "SAPHIR" "SSMIS") 
## GPM-calibrated level 1C version
GPM_1C_ver=07 #05
## Path to all GPM-calibrated level 1C data
GPM_datapath="https://gpm1.gesdisc.eosdis.nasa.gov/data/GPM_L1C/"

# Evaluate if SSMI data is used
SSMI="true"
## Path to SSMI Tb data
SSMI_datapath="https://www.ncei.noaa.gov/data/ssmis-brightness-temperature-csu/access/FCDR/"
## ICDR: interim climate data record is used here because the more accurate one 
## & FCDR is not available for this period. FCDR is better.

# ===============================================================
# Utility
# ===============================================================

# Identify download path
DownloadTo=${Project_Path}/raw_Obs/Microwave/${Storm}/

# Identify the year 
Year=${Timebeg:0:4}

# For GPM data
if [[ ${GPM_1C} == "true" ]]
then
	# If the end time is the beginning of N+1 day (e.g., 2017082300),
	# & only download the files from the beginning day to the Nth day.
	Timend_hh=${Timend:8:10}
	if [[ ${Timend_hh} == '00' ]]
	then
    	Timend=$(date -u -d -1' days '${Timend:0:8}' '${Timend_hh} +%Y%m%d%H)
	fi

	# Convert time string to the day of the year (e.g., 20170823 is 235th day of the year)
	# The data is stored based on the day of the year on EARTHDATA website
	date_Tbeg=$(date -d "${Timebeg:0:8}" +'%Y-%m-%d')
	date_Tbeg_ofYear=$(date +%j --date="${date_Tbeg}")
	date_Tend=$(date -d "${Timend:0:8}" +'%Y-%m-%d')
	date_Tend_ofYear=$(date +%j --date="${date_Tend}")
fi

# For SSMI data
if [[ ${SSMI} == "true" ]]
then
	Daybeg=${Timebeg:0:8}
	Dayend=${Timend:0:8}

    # If the end time is the beginning of N+1 day (e.g., 2017082300),
    # & only download the files from the beginning day to the Nth day.
    Timend_hh=${Timend:8:10}
    if [[ ${Timend_hh} == '00' ]]
    then
        Dayend=$(date -u -d -1' days '"${Dayend}" +%Y%m%d)
    fi
fi

# ===============================================================
# Download GPM-calibrated level 1C Tbs
# ===============================================================
# Data will be download from EARTHDATA website.
# ! Make sure you have set up your Earthdata account. 
# Sequentially check each sensor.
# -------------------------------------------------------------

if [[ ${GPM_1C} == "false" ]]
then
 
	echo "GPM-calibrated 1C data is not being downloaded!"

else
	# Prepare before the download
	## Check if the netrc file exists in your home directy
	if [ ! -f "$HOME"/.netrc ]
	then
		echo "$HOME/.netrc file unavailable" >&2
		echo "Search the web for how to set up .netrc" >&2
		exit 1
	else 
		if ! grep urs.earthdata.nasa.gov "$HOME"/.netrc -q
		then
			echo "urs.earthdata.nasa.gov entry not found in $HOME/.netrc" >&2
			exit 1
		fi
	fi
	## Create a cookie file. This file will be used to persist sessions across calls to wget command
	cd ~ || exit
	touch .urs_cookies 
	
	# Function to download
	Download_GPM_granule () {
		cd "$1" || exit
		dayOfYear=${date_Tbeg_ofYear}
		while (( $((10#${dayOfYear})) <= $((10#${date_Tend_ofYear})) )); do

			echo "Day of the year: ${dayOfYear} in ${Year}"
			dataset="$2/${Year}/${dayOfYear}/"
			echo ${dataset}
			wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A HDF5 --content-disposition "${dataset}"
			dayOfYear=$(($((10#${dayOfYear}))+1))

		done
		cd ..
	}

	# Loop each sensor in GPM mission 
	cd ${DownloadTo} || exit
	
	## Echo date and some msg
	date
	echo "Starting to download GPM data"

	for iss in ${!GPM_sensor[@]}; do
		sensor=${GPM_sensor[iss]}	
		
		case ${sensor} in
			AMSR2)
				cd ${DownloadTo} || exit
				if [ ! -d AMSR2 ]
				then
					mkdir AMSR2
				fi

			    if [ ! -d AMSR2/GCOMW1 ]
                then
                    mkdir AMSR2/GCOMW1
					dataset="${GPM_datapath}GPM_1CGCOMW1AMSR2.${GPM_1C_ver}/"
					Download_GPM_granule AMSR2/GCOMW1 ${dataset}
					cd ${DownloadTo} || exit
                fi	
				;;
			
			ATMS)
				if [ ! -d ATMS ]
				then
					mkdir ATMS
				fi

                if [ ! -d ATMS/NPP ]
                then
                    mkdir ATMS/NPP
					dataset="${GPM_datapath}GPM_1CNPPATMS.${GPM_1C_ver}/"
					Download_GPM_granule ATMS/NPP ${dataset}
					cd ${DownloadTo} || exit
                fi	
				;;

			GMI)
				if [ ! -d GMI ]
				then
					mkdir GMI
				fi

                if [ ! -d GMI/GPM ]
                then
                    mkdir GMI/GPM
					dataset="${GPM_datapath}GPM_1CGPMGMI.${GPM_1C_ver}/"
					Download_GPM_granule GMI/GPM ${dataset} 
					cd ${DownloadTo} || exit
                fi
				;;

			MHS)
				cd ${DownloadTo} || exit
                if [ ! -d MHS ]
                then
                    mkdir MHS
                fi

				if [ ! -d MHS/METOPA ]
				then
					mkdir MHS/METOPA
					dataset="${GPM_datapath}GPM_1CMETOPAMHS.${GPM_1C_ver}/"
					Download_GPM_granule MHS/METOPA ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d MHS/METOPB ]
				then
					mkdir MHS/METOPB
					dataset="${GPM_datapath}GPM_1CMETOPBMHS.${GPM_1C_ver}/"
					Download_GPM_granule MHS/METOPB ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d MHS/NOAA18 ]
				then
					mkdir MHS/NOAA18
					dataset="${GPM_datapath}GPM_1CNOAA18MHS.${GPM_1C_ver}/"
					Download_GPM_granule MHS/NOAA18 ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d MHS/NOAA19 ]
				then
					mkdir MHS/NOAA19
					dataset="${GPM_datapath}GPM_1CNOAA19MHS.${GPM_1C_ver}/"
					Download_GPM_granule MHS/NOAA19 ${dataset}
					cd ${DownloadTo} || exit
				fi
				;;

			SAPHIR)
				cd ${DownloadTo} || exit
                if [ ! -d SAPHIR ]
                then
                    mkdir SAPHIR
                fi

				if [ ! -d SAPHIR/MT1 ]
				then
					mkdir SAPHIR/MT1
					dataset="${GPM_datapath}GPM_1CMT1SAPHIR.${GPM_1C_ver}/"
					Download_GPM_granule SAPHIR/MT1 ${dataset}
					cd ${DownloadTo} || exit
				fi
				;;

			SSMIS)
				cd ${DownloadTo} || exit
                if [ ! -d SSMIS ]
                then
                    mkdir SSMIS
                fi

				if [ ! -d SSMIS/F16 ]
				then
					mkdir SSMIS/F16
					dataset="${GPM_datapath}GPM_1CF16SSMIS.${GPM_1C_ver}/"
					Download_GPM_granule SSMIS/F16 ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d SSMIS/F17 ]
				then
					mkdir SSMIS/F17
					dataset="${GPM_datapath}GPM_1CF17SSMIS.${GPM_1C_ver}/"
					Download_GPM_granule SSMIS/F17 ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d SSMIS/F18 ]
				then
					mkdir SSMIS/F18
					dataset="${GPM_datapath}GPM_1CF18SSMIS.${GPM_1C_ver}/"
					Download_GPM_granule SSMIS/F18 ${dataset}
					cd ${DownloadTo} || exit
				fi

				if [ ! -d SSMIS/F19 ]
				then
					mkdir SSMIS/F19
					dataset="${GPM_datapath}GPM_1CF19SSMIS.${GPM_1C_ver}/"
					Download_GPM_granule SSMIS/F19 ${dataset}
					cd ${DownloadTo} || exit	
				fi
				;;
				
		esac	
	done #---- End of sensor loop---
fi


# ===============================================================
# Download not-GMI-calibrated sensor products: SSMI
# ===============================================================
# 
# -------------------------------------------------------------
cd ${DownloadTo} || exit

if [[ ${SSMI} == "false" ]]
then

    echo "SSMI data is not being downloaded!"

else

	## Echo date and some msg
    date
    echo "Starting to download SSMI data"

	if [ ! -d SSMI ]
	then
		mkdir SSMI
	fi

    if [ ! -d SSMI/F15 ]
    then
        mkdir SSMI/F15
    fi

	cd SSMI/F15 || exit
	currentDay=${Daybeg}
	while (( ${currentDay} <= ${Dayend} )); do

		echo "Date: ${currentDay}"
		dataset="${SSMI_datapath}/${Year}/"
		echo ${dataset}
	    wget --no-verbose --no-parent --recursive --level=1 --no-directories -A "CSU_SSMI_FCDR_V02R00_F15_D${currentDay}_*.nc" --content-disposition "${dataset}"
    	currentDay=$(date -u -d 1' days '${currentDay} +%Y%m%d)

	done

fi


exit 0
