#!/bin/bash
# This namelist_enkf.sh is modified mainly based on the version used for the TMeCSR and the enkf code on the GitHub. A few parameters are modified to be consistent with the version that Judy used on Harvey, Irma, Jose, and Maria.

# Fortran reads the namelist based on variable names, not order.
# Any missing variables retain their default values (if not assigned in the input file)

. $CONFIG_FILE
domain_id=$1
dx=$((DX[$domain_id-1]/1000))

# =================================================================
# 
#   LIST OF VARIABLES FOR ENKF TO LOAD AND UPDATE
#
# =================================================================


# Variables to load, domain controls and ensemble size
# -------------------------------------------------------
cat << EOF
/
&enkf_parameter
numbers_en   = $NUM_ENS, 
expername    = '$EXP_NAME', 
enkfvar      = 'T         ', 'U         ', 'V         ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', 'QICE      ', 'QSNOW     ', 'QGRAUP    ', 'PH        ', 'MU        ', 'PSFC      ', 'P         ', 'PHB       ', 'PB        ', 'MUB       ', 'U10       ', 'V10       ', 'TSK       ','W         ',
updatevar    = 'T         ', 'U         ', 'V         ', 'QVAPOR    ', 'QCLOUD    ', 'QRAIN     ', 'QICE      ', 'QSNOW     ', 'QGRAUP    ', 'PH        ', 'MU        ', 'PSFC      ', 'P         ', 'U10       ', 'V10       ', 'TSK       ', 
update_is    = 1,
update_ie    = ${E_WE[$domain_id-1]},
update_js    = 1,
update_je    = ${E_SN[$domain_id-1]},
update_ks    = 1,
update_ke    = ${E_VERT[$domain_id-1]},
inflate      = $INFLATION_COEF,
mixing       = $RELAXATION_COEF,
use_sfc_theta= true,
use_sfc_td   = true
random_order = .false.,
print_detail = 0,
use_vortex_RTPP = .$USE_VORTEX_RTPP.,
mixing_core  = $RELAXATION_COEF_CORE,
use_nonlinear_enkf=.$USE_NONLINEAR_ENKF.,
use_gmm_enkf = .$USE_GMM_ENKF.,
max_kernel_num = $MAX_KERNEL_NUM,
min_expanding_kernel_prior_size= $MIN_EXPANDING_KERNEL_PRIOR_SIZE,
/
# -------------------------------------------------------
# - vortex_RTPP: use_vortex_RTPP,mixing_core 
# -------------------------------------------------------
# - Various version of EnKF filters:
# use_nonlinear_enkf,use_gmm_enkf,max_kernel_num,min_expanding_kernel_prior_size
# -------------------------------------------------------
EOF

# Writing enkf parallelization and observation controls
# -----------------------------------------------------
cat << EOF
/

&parallel
manual_parallel = .true.,
nmcpu  = $NMCPU,
nicpu  = $NICPU,
njcpu  = $NJCPU,
/

&osse
use_ideal_obs    = .false.,
gridobs_is   = 0,
gridobs_ie   = 0,
gridobs_js   = 0,
gridobs_je   = 0,
gridobs_ks   = 0,
gridobs_ke   = 0,
gridobs_int_x= 0,
gridobs_int_k= 0,
use_simulated= .false.,
/

&hurricane_PI 
use_hurricane_PI  = .false.,
hroi_hurricane_PI = 60,
vroi_hurricane_PI = 35,
/

&surface_obs
use_surface      = .$USE_SURFOBS.,
datathin_surface = $THIN_SURFACE,
hroi_surface     = $((HROI_SFC/$dx)),
vroi_surface     = $VROI,
/

&sounding_obs
use_sounding      = .$USE_SOUNDOBS.,
datathin_sounding = $THIN_SOUNDING,
hroi_sounding     = $((HROI_UPPER/$dx)),
vroi_sounding     = $VROI,
/

&profiler_obs
use_profiler      = .$USE_PROFILEROBS.,
datathin_profiler = $THIN_PROFILER,
hroi_profiler     = $((HROI_UPPER/$dx)),
vroi_profiler     = $VROI,
/

&aircft_obs
use_aircft      = .$USE_AIREPOBS.,
datathin_aircft = $THIN_AIRCFT,
hroi_aircft     = $((HROI_UPPER/$dx)),
vroi_aircft     = $VROI,
/

&metar_obs
use_metar      = .$USE_METAROBS.,
datathin_metar = $THIN_METAR,
hroi_metar     = $((HROI_SFC/$dx)),
vroi_metar     = 999,
/

&sfcshp_obs
use_sfcshp      = .$USE_SHIPSOBS.,
datathin_sfcshp = $THIN_SFCSHP,
hroi_sfcshp     = $((HROI_SFC/$dx)),
vroi_sfcshp     = $VROI,
/

&spssmi_obs
use_spssmi      = .$USE_SSMIOBS.,
datathin_spssmi = $THIN_SPSSMI,
hroi_spssmi     = $((HROI_UPPER/$dx)),
vroi_spssmi     = $VROI,
/

&atovs_obs
use_atovs      = .$USE_ATOVS.,
datathin_atovs = $THIN_ATOVS,
hroi_atovs     = $((HROI_UPPER/$dx)),
vroi_atovs     = $VROI,
/

&satwnd_obs
use_satwnd      = .$USE_GEOAMVOBS.,
datathin_satwnd = $THIN_SATWND,
hroi_satwnd     = $((HROI_UPPER/$dx)),
vroi_satwnd     = $VROI,
/

&gpspw_obs
use_gpspw      = .$USE_GPSPWOBS.,
datathin_gpspw = $THIN_GPSPW,
hroi_gpspw     = $((HROI_SFC/$dx)),
vroi_gpspw     = $VROI,
/

&radar_obs
radar_number   = 1,
use_radar_rf   = .$USE_RADAR_RF.,
use_radar_rv   = .$USE_RADAR_RV., 
datathin_radar = $THIN_RADAR,
hroi_radar     = $((HROI_RADAR/$dx)),
vroi_radar     = $VROI_RADAR,
/

&airborne_radar   
use_airborne_rf   = .$USE_AIRBORNE_RF.,
use_airborne_rv   = .$USE_AIRBORNE_RV.,
datathin_airborne = $THIN_RADAR,
hroi_airborne     = $((HROI_RADAR/$dx)),
vroi_airborne     = $VROI_RADAR,
/

&radiance 
use_radiance  = .$USE_RADIANCE.,
hroi_radiance = $((HROI_RADIANCE/$dx)),
vroi_radiance = $VROI_RADIANCE,
datathin_radiance = $THIN_RADIANCE,
use_vroi_radiance_halfsfc = .$USE_VROI_RADIANCE_HALFSFC.,
use_aoei      = .$USE_AOEI_RADIANCE.,
use_elf       = .$USE_ELF_RADIANCE.,
use_abei      = .$USE_INF_BT_RADIANCE.,
/
# -------------------------------------------------------
# "hroi_radiance" and "vroi_radiance" are actually specified 
# in preprocessed radiance_*_so file. Radiance obs are only
# assimilated in the domain 3.
# -------------------------------------------------------

&microwave
use_microwave         = .$USE_MICROWAVE.,
hroi_microwave        = $((HROI_MICROWAVE/$dx)),
vroi_microwave        = $VROI_MICROWAVE,
datathin_microwave    = $THIN_MICROWAVE,
aoei_microwave        = .$USE_AOEI_MICROWAVE.,
use_slant_path        = .$USE_SLANT_PATH.,
EOF
if [[ $domain_id == 3 ]]; then
cat << EOF
crtm_mp             = $CRTM_MP_SCHEME,
crtm_rainLUT        = $CRTM_RAINLUT
crtm_snowLUT        = $CRTM_SNOWLUT
crtm_graupelLUT     = $CRTM_GRAUPELLUT
EOF
fi
cat << EOF
/
# -------------------------------------------------------
# "hroi_microwave" and "vroi_microwave" are actually specified
# in preprocessed microwave_*_so file. Microwave obs are only
# assimilated in the domain 3.
# -------------------------------------------------------
# crtm_mp:  wsm6(0), gfdlfv3(1), thompson08(2), nssl(3)
# -------------------------------------------------------

EOF
