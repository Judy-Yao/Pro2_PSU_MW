#use_for = "perturb"  Use WRF 3DVar randomcv mode to generate ensemble perturbations
#use_for = "ndown"    Running ndown.exe on nested domain $idom
#use_for = "wrfw"     Running wrf.exe for domain $RUN_DOMAIN across the cycle window, 
#                     wrfinput files will be generated for next cycle

. $CONFIG_FILE
end_date=`advance_time $start_date $run_minutes`
use_for=$1
idom=$2
#echo $idom

domlist=`seq 1 $MAX_DOM`
LBINT=$LBC_INTERVAL
if [ $use_for == "perturb" ]; then
  MAX_DOM=1
  domlist=$idom
fi
if [ $use_for == "ndown" ]; then
  MAX_DOM=2
  domlist="${PARENT_ID[$idom-1]} $idom"
  LBINT=$((WRFOUT_INTERVAL[${PARENT_ID[$idom-1]}-1]))
fi

if [ -f ij_parent_start ]; then
  #echo using calculated ij_parent_start
  i_parent_start=(`cat ij_parent_start |head -n1`)
  j_parent_start=(`cat ij_parent_start |tail -n1`)
else
  for n in `seq 1 $MAX_DOM`; do
    i_parent_start[$n-1]=${I_PARENT_START[$n-1]}
    j_parent_start[$n-1]=${J_PARENT_START[$n-1]}
  done
fi

#echo ${i_parent_start[0]}
#echo ${i_parent_start[1]}
#echo ${i_parent_start[2]}

#=============TIME CONTROL PART=============
echo "&time_control"
#cat << EOF
#run_minutes = $FORECAST_MINUTES, 
#EOF
if $RESTARTING; then
  echo start_year         = `for i in $domlist; do printf ${restart_date:0:4}, ; done`
  echo start_month        = `for i in $domlist; do printf ${restart_date:4:2}, ; done`
  echo start_day          = `for i in $domlist; do printf ${restart_date:6:2}, ; done`
  echo start_hour         = `for i in $domlist; do printf ${restart_date:8:2}, ; done`
  echo start_minute       = `for i in $domlist; do printf ${restart_date:10:2}, ; done`
  echo start_second       = `for i in $domlist; do printf 00, ; done`
else
  echo start_year         = `for i in $domlist; do printf ${start_date:0:4}, ; done`
  echo start_month        = `for i in $domlist; do printf ${start_date:4:2}, ; done`
  echo start_day          = `for i in $domlist; do printf ${start_date:6:2}, ; done`
  echo start_hour         = `for i in $domlist; do printf ${start_date:8:2}, ; done`
  echo start_minute       = `for i in $domlist; do printf ${start_date:10:2}, ; done`
  echo start_second       = `for i in $domlist; do printf 00, ; done`
fi
cat << EOF
end_year           = `for i in $domlist; do printf ${end_date:0:4}, ; done`
end_month          = `for i in $domlist; do printf ${end_date:4:2}, ; done`
end_day            = `for i in $domlist; do printf ${end_date:6:2}, ; done`
end_hour           = `for i in $domlist; do printf ${end_date:8:2}, ; done`
end_minute         = `for i in $domlist; do printf ${end_date:10:2}, ; done`
end_second         = `for i in $domlist; do printf 00, ; done`
input_from_file    = .true.,.true.,.true.,.false.,
fine_input_stream  = 0, 0, 0, 0,                                                         
input_from_hires   = .false.,.false.,.false.,.false.,                                                 
rsmas_data_path    = '.',
interval_seconds   = $((LBINT*60)),
history_interval   = `for i in $domlist; do printf ${WRFOUT_INTERVAL[$i-1]}, ; done`
frames_per_outfile = `for i in $domlist; do printf 1, ; done`
EOF
if $RESTARTING; then
  echo restart                             = .true.,
  echo restart_interval                              = 360,
else
  echo restart                             = .false.,
  echo restart_interval                              = 7300,
fi
cat << EOF
debug_level        = 500,
EOF

#input_from_file    = `for i in $domlist; do printf .true., ; done`

if [[ $use_for == "wrfw" ]]; then
cat << EOF
input_outname="wrfinput_d<domain>_<date>",
write_input=true,
inputout_interval=$inputout_interval,
inputout_begin_m=$inputout_begin,
inputout_end_m=$inputout_end,
EOF
fi

if [[ $use_for == "ndown" ]]; then
  echo io_form_auxinput2=2,
fi

if [[ $SST_UPDATE == 1 ]]; then
  dmin=`min ${CYCLE_PERIOD[@]}`
  echo auxinput4_inname="wrflowinp_d<domain>",
  echo auxinput4_interval=180,180,180,
  echo auxinput4_end       = 0,
  echo io_form_auxinput4=2,
  echo inputout_interval   = 180,
  echo inputout_begin_h    = 0,
  echo inputout_end_h      = 3,
fi

echo "/"

#=============DOMAIN PART=============
echo "&domains"
cat << EOF
time_step=$DT,
time_step_fract_num                 = 0,
time_step_fract_den                 = 1,
time_step_dfi                       = 90,
min_time_step                       = -1, -1, -1, -1,
max_time_step                       = -1, -1, -1, -1,
target_cfl                          = 1.2, 1.2, 1.2, 1.2,
target_hcfl                         = 0.84, 0.84, 0.84, 0.84,
max_step_increase_pct               = 5, 5, 5, 5,
starting_time_step                  = -1, -1, -1, -1,
step_to_output_time                 = .true.,
adaptation_domain                   = 1,
use_adaptive_time_step              = .false.,
use_adaptive_time_step_dfi          = .false.,
max_dom=$MAX_DOM,
lats_to_mic                         = 0,
s_we                                = 1,     1,     1,     1,
e_we       = `for i in $domlist; do printf ${E_WE[$i-1]}, ; done`
e_sn       = `for i in $domlist; do printf ${E_SN[$i-1]}, ; done`
s_sn                                = 1,     1,     1,     1,
e_vert     = `for i in $domlist; do printf ${E_VERT[$i-1]}, ; done`
s_vert                              = 1,     1,     1,     1,
dx         = `for i in $domlist; do printf ${DX[$i-1]}, ; done`
dy         = `for i in $domlist; do printf ${DY[$i-1]}, ; done`
EOF

if [[ $use_for == "ndown" ]] && [[ $idom == 2 ]]; then
cat << EOF
grid_id    = 1,2,
parent_id  = 0,1,
parent_grid_ratio = 1,${GRID_RATIO[$idom-1]},
i_parent_start = 1,${i_parent_start[$idom-1]},
j_parent_start = 1,${j_parent_start[$idom-1]},
EOF
elif [[ $use_for == "ndown" ]] && [[ $idom == 3 ]]; then
cat << EOF
grid_id    = 2,3,
parent_id  = 1,2,
parent_grid_ratio = ${GRID_RATIO[$idom-2]},${GRID_RATIO[$idom-1]},
i_parent_start = ${i_parent_start[$idom-2]},${i_parent_start[$idom-1]},
j_parent_start = ${j_parent_start[$idom-2]},${j_parent_start[$idom-1]},
EOF
else
cat << EOF
grid_id    = `for i in $domlist; do printf $i, ; done`
parent_id  = 0,1,2,3
parent_grid_ratio = 1,`for i in $(seq 2 $MAX_DOM); do printf ${GRID_RATIO[$i-1]}, ; done`
parent_time_step_ratio = 1,`for i in $(seq 2 $MAX_DOM); do printf ${TIME_STEP_RATIO[$i-1]}, ; done`
i_parent_start = 1,`for i in $(seq 2 $MAX_DOM); do printf ${i_parent_start[$i-1]}, ; done`100,
j_parent_start = 1,`for i in $(seq 2 $MAX_DOM); do printf ${j_parent_start[$i-1]}, ; done`100,
EOF
fi

if $TWO_WAY_NESTING; then
  echo "feedback=1,"
else
  echo "feedback=0,"
fi

#Preset moves for following TC
echo time_to_move = 0,
echo corral_dist  = 8, 8, 98, 98,
echo track_level = 70000,

cat << EOF
smooth_option=0,
num_metgrid_levels=27,
p_top_requested=$P_TOP,
num_metgrid_soil_levels=4,
interp_theta                        = .false.,
interp_type                         = 2,
vert_refine_fact                    = 1,
extrap_type                         = 2,
t_extrap_type                       = 2,
hypsometric_opt                     = 2,
lowest_lev_from_sfc                 =.false.,
use_levels_below_ground             =.true.,
use_tavg_for_tsk                    =.false.,
use_surface                         =.true.,
lagrange_order                      = 2,
force_sfc_in_vinterp                = 1,
zap_close_levels                    = 500,
sfcp_to_sfcp                        = .false.,
adjust_heights                      = .false.,
smooth_cg_topo                      = .false.,
nest_interp_coord                   = 0,
aggregate_lu                        = .false.,
rh2qv_wrt_liquid                    = .true.,
rh2qv_method                        = 1,
qv_max_p_safe                       = 10000,
qv_max_flag                         = 1.E-5,
qv_max_value                        = 3.E-6,
qv_min_p_safe                       = 110000,
qv_min_flag                         = 1.E-6,
qv_min_value                        = 1.E-6,
eta_levels                          = 1.0,.9919699,.9827400,.9721600,.9600599,.9462600,.9306099,.9129300,.8930600,.8708600,.8462000,.8190300,.7893100,.7570800,.7224600,.6856500,.6469100,.6066099,.5651600,.5230500,.4807700,.4388600,.3978000,.3580500,.3200099,.2840100,.2502900,.2190100,.1902600,.1640600,.1403600,.1190600,.1000500,.0831600,.0682400,.0551200,.0436200,.0335700,.0248200,.0172200,.0106300,.0049200,.0000000,
EOF

echo "/"

#=============PHYSICS PART=============
echo "&physics"
if $GET_PHYS_FROM_FILE; then
cat << EOF
mp_physics         = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :MP_PHYSICS |awk '{print $3}'), ; done`
ra_lw_physics      = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :RA_LW_PHYSICS |awk '{print $3}'), ; done`
ra_sw_physics      = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :RA_SW_PHYSICS |awk '{print $3}'), ; done`
radt               = `for i in $domlist; do printf ${RADT[$i-1]}, ; done`
sf_sfclay_physics  = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :SF_SFCLAY_PHYSICS |awk '{print $3}'), ; done`
sf_surface_physics = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :SF_SURFACE_PHYSICS |awk '{print $3}'), ; done`
bl_pbl_physics     = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :BL_PBL_PHYSICS |awk '{print $3}'), ; done`
bldt               = `for i in $domlist; do printf ${BLDT[$i-1]}, ; done`
cu_physics         = `for i in $domlist; do printf $(ncdump -h wrfinput_d$(expr $i + 100 |cut -c2-) |grep :CU_PHYSICS |awk '{print $3}'), ; done`
cudt               = `for i in $domlist; do printf ${CUDT[$i-1]}, ; done`
EOF
else
cat << EOF
mp_physics         = `for i in $domlist; do printf ${MP_PHYSICS[$i-1]}, ; done`
ra_lw_physics      = `for i in $domlist; do printf ${RA_LW_PHYSICS[$i-1]}, ; done`
ra_sw_physics      = `for i in $domlist; do printf ${RA_SW_PHYSICS[$i-1]}, ; done`
radt               = `for i in $domlist; do printf ${RADT[$i-1]}, ; done`
sf_sfclay_physics  = `for i in $domlist; do printf ${SF_SFCLAY_PHYSICS[$i-1]}, ; done`
sf_surface_physics = `for i in $domlist; do printf ${SF_SURFACE_PHYSICS[$i-1]}, ; done`
bl_pbl_physics     = `for i in $domlist; do printf ${BL_PBL_PHYSICS[$i-1]}, ; done`
bldt               = `for i in $domlist; do printf ${BLDT[$i-1]}, ; done`
grav_settling      = `for i in $domlist; do printf 0, ; done`
cu_physics         = `for i in $domlist; do printf ${CU_PHYSICS[$i-1]}, ; done`
cudt               = `for i in $domlist; do printf ${CUDT[$i-1]}, ; done`
EOF
fi

cat << EOF
isfflx             = 1,
ifsnow             = 0,
icloud             = 1,
surface_input_source=1,
num_soil_layers    = 4,
maxiens            = 1,
maxens             = 3,
maxens2            = 3,
maxens3            = 16,
ensdim             = 144,
seaice_threshold   = 271,
sst_update         = $SST_UPDATE,
sst_skin           = 1,
sf_ocean_physics   = 3
isftcflx           = 4,
EOF

#extra physics options
#cat << EOF
#levsiz = 59,
#paerlev = 29,
#cam_abs_dim1 = 4,
#cam_abs_dim2 = 45,
#EOF

echo "/"


#=============DYNAMICS PART=============
echo "&dynamics"
cat << EOF
w_damping            = 0,
use_input_w          = .false.,
diff_opt             = 0,
km_opt               = 1,
diff_6th_opt         = `for i in $domlist; do printf 0, ; done`
diff_6th_factor      = `for i in $domlist; do printf 0.12, ; done`
base_temp            = 290.,
damp_opt             = 3,
zdamp                = `for i in $domlist; do printf 5000., ; done`
dampcoef             = `for i in $domlist; do printf 0.2, ; done`
khdif                = `for i in $domlist; do printf 0, ; done`
kvdif                = `for i in $domlist; do printf 0, ; done`
smdiv                = `for i in $domlist; do printf 0.1, ; done`
emdiv                = `for i in $domlist; do printf 0.01, ; done`
epssm                = `for i in $domlist; do printf 0.1, ; done`
non_hydrostatic      = `for i in $domlist; do printf .true., ; done`
h_mom_adv_order      = `for i in $domlist; do printf 5, ; done`
v_mom_adv_order      = `for i in $domlist; do printf 3, ; done`
h_sca_adv_order      = `for i in $domlist; do printf 5, ; done`
v_sca_adv_order      = `for i in $domlist; do printf 3, ; done`
use_baseparam_fr_nml = .true.

EOF
echo "/"

#=============BODY CONTROL PART=============
echo "&bdy_control"
cat << EOF
spec_bdy_width     = 5,
spec_zone          = 1,
relax_zone         = 4,
specified          = .true., .false.,.false.,.false.,
periodic_x         = .false.,.false.,.false.,.false.,
symmetric_xs       = .false.,.false.,.false.,.false.,
symmetric_xe       = .false.,.false.,.false.,.false.,
open_xs            = .false.,.false.,.false.,.false.,
open_xe            = .false.,.false.,.false.,.false.,
periodic_y         = .false.,.false.,.false.,.false.,
symmetric_ys       = .false.,.false.,.false.,.false.,
symmetric_ye       = .false.,.false.,.false.,.false.,
open_ys            = .false.,.false.,.false.,.false.,
open_ye            = .false.,.false.,.false.,.false.,
nested             = .false., .true., .true., .true.,
EOF
echo "/"

#=============OTHERS=============
cat << EOF
&noah_mp
/
&fdda
/
&scm
/
&grib2
/
&fire
/
&diags
/
&namelist_quilt
/
&tc
/
&logging
/
&dfi_control
/
&namelist_quilt 
/
EOF
