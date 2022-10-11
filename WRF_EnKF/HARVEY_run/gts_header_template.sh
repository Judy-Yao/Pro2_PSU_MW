#!/bin/bash
tot_obs_num=$1

cat << EOF
Total =$tot_obs_num, MISS. =-888888.,
SYNOP =      0, METAR =      0, SHIP  =      0, BUOY  =      0, BOGUS =      0, TEMP  =      0,
AMDAR =      0, AIREP =      0, TAMDAR=      0, PILOT =      0, SATEM =      0, SATOB =      0,
GPSPW =      0, GPSZD =      0, GPSRF =      0, GPSEP =      0, SSMT1 =      0, SSMT2 =      0,
TOVS  =      0, QSCAT =      0, PROFL =      0, AIRSR =      0, OTHER =      0,
PHIC  =   0.19, XLONC =  75.47, TRUE1 =   0.00, TRUE2 = -30.00, XIM11 =   1.00, XJM11 =   1.00,
base_temp= 290.00, base_lapse=  50.00, PTOP  = 2000., base_pres=100000., base_tropo_pres= 20000., base_strat_temp=   215.,
IXC   =    223, JXC   =    334, IPROJ =      3, IDD   =      1, MAXNES=      1,
NESTIX=    223,
NESTJX=    334,
NUMC  =      1,
DIS   =   9.00,
NESTI =      1,
NESTJ =      1,
INFO  = PLATFORM, DATE, NAME, LEVELS, LATITUDE, LONGITUDE, ELEVATION, ID.
SRFC  = SLP, PW (DATA,QC,ERROR).
EACH  = PRES, SPEED, DIR, HEIGHT, TEMP, DEW PT, HUMID (DATA,QC,ERROR)*LEVELS.
INFO_FMT = (A12,1X,A19,1X,A40,1X,I6,3(F12.3,11X),6X,A40)
SRFC_FMT = (F12.3,I4,F7.2,F12.3,I4,F7.3)
EACH_FMT = (3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2),11X,3(F12.3,I4,F7.2))
#------------------------------------------------------------------------------#
EOF
