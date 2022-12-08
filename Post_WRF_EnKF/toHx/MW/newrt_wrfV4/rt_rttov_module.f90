module rt_rttov_module
  
  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  use rt_state_module
  use rt_result_module

  IMPLICIT NONE

! IR and vis
#include "rttov_direct.interface"
#include "rttov_parallel_direct.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"

! MW
#include "rttov_scatt.interface"
#include "rttov_parallel_scatt.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_print_cld_profile.interface"

! both
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

  PRIVATE

  public :: rt_rttov_n_channels
  public :: rt_rttov_init
  public :: state_to_rttov
  public :: rt_rttov_destroy
  public :: rt_rttov_main
  public :: rttov_to_result
  
  integer, parameter :: nprof = 1
  integer, parameter :: iprof = 1 

  ! RTTOV variables/structures
  !====================
  ! for VIS, IR, MW
  TYPE(rttov_options)              :: opts                     ! Options structure
  TYPE(rttov_coefs)                :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)    => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:) => NULL() ! Input/output surface BRDF
  TYPE(rttov_profile),     POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_transmission)         :: transmission             ! Output transmittances
  TYPE(rttov_radiance)             :: radiance                 ! Output radiances

  ! only used for MW sensors
  TYPE(rttov_options_scatt)          :: opts_scatt               ! RTTOV-SCATT options structure
  TYPE(rttov_scatt_coef)             :: coefs_scatt              ! RTTOV-SCATT coefficients structure
  INTEGER(KIND=jpim),        POINTER :: frequencies(:) => NULL() ! Channel indexes for Mietable lookup
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)  => NULL() ! Flags to specify channels to simulate
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)=> NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles

  INTEGER(KIND=jpim)               :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status

  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim) :: nlevels

  integer, save :: nchannels = -1
  integer, save :: channels(MAX_CHANNELS) = 0
  logical, save :: l_initialized = .false.
  logical, save :: l_isMW = .false.
  character(len=LENGTH_SENSOR_NAME), save :: sensor_id_rt(MAX_SENSORS)

contains
  function rt_rttov_n_channels()
    integer :: rt_rttov_n_channels
    rt_rttov_n_channels = nchannels
  end function rt_rttov_n_channels

  subroutine rt_rttov_init(state, sensor_name, channel_list)
    type(model_state),  intent(in) :: state
    character(len=*),   intent(in) :: sensor_name
    INTEGER(KIND=jpim), intent(in) :: channel_list(:)
    
    character(len=256) :: coef_filename
    character(len=256) :: cld_coef_filename
    character(len=256) :: mietable_filename
    
!    INTEGER(KIND=jpim) :: nprof
!    INTEGER(KIND=jpim) :: nchannels
    ! loop variables
    INTEGER(KIND=jpim) :: j, jch, lay
    INTEGER(KIND=jpim) :: nch
    INTEGER(KIND=jpim) :: joff
    INTEGER            :: ios
    
    nchannels = size(channel_list)
    channels(1:nchannels) = channel_list
    nlevels   = state%nz
   
    if (my_proc_id == 0) then
      print *, 'rt_rttov_init start'
      print *, 'n_channels ', nchannels
      print *, 'n_levels', nlevels, state%nz!, state%xbeg, state%xend
    end if

    if (sensor_name == 'abi_gr') then
      coef_filename     = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/rttov9pred54L/rtcoef_goes_16_abi.dat'
      cld_coef_filename = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/cldaer_visir/sccldcoef_goes_16_abi.dat'
      l_isMW = .false.
    else if (sensor_name == 'gmi_gpm_lf' .or. sensor_name == 'gmi_gpm_hf') then
      coef_filename     = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/rttov7pred54L/rtcoef_gpm_1_gmi.dat'
      mietable_filename = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/mietable/mietable_gpm_gmi.dat'
      l_isMW = .true.
    else if (sensor_name == 'ssmis_f17') then
      coef_filename     = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/rttov7pred54L/rtcoef_dmsp_17_ssmis.dat'
      mietable_filename = '/work/05012/tg843115/stampede2/opt/RTTOV/rttov123_single_thread/rtcoef_rttov12/mietable/mietable_dmsp_ssmis.dat'
      l_isMW = .true.
    end if

    if (nml_s_rttov_mietable .ne. '') then
      mietable_filename = nml_s_rttov_mietable
    end if

    if (.not. l_isMW) then
      opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
      opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects
      opts % rt_ir % ir_scatt_model      = 2       ! Scattering model for emission source term:
                                                   !   1 => DOM; 2 => Chou-scaling
      opts % rt_ir % vis_scatt_model     = 1       ! Scattering model for solar source term:
                                                   !   1 => DOM; 2 => single-scattering; 3 => MFASIS
      opts % rt_ir % dom_nstreams        = 8       ! Number of streams for Discrete Ordinates (DOM)

      opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
      opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
      opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
      opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
      opts % rt_ir % co_data             = .FALSE. !
      opts % rt_ir % so2_data            = .FALSE. !
    end if
    
    if (nml_l_include_cloud) then
      opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
      if (.not. l_isMW) then
        opts % rt_ir % addclouds           = .TRUE.  ! Include cloud effects
        opts % rt_ir % grid_box_avg_cloud  = .TRUE.  ! Cloud concentrations are grid box averages
      else
        opts_scatt % interp_mode = 1
        opts_scatt % config % apply_reg_limits = .TRUE.
        opts_scatt % config % verbose = .TRUE.
      end if
    else
      opts % rt_all % addrefrac          = .FALSE.  ! Include refraction in path calc
      if (.not. l_isMW) then
        opts % rt_ir % addclouds           = .FALSE.  ! Include cloud effects
        opts % rt_ir % grid_box_avg_cloud  = .FALSE.  ! Cloud concentrations are grid box averages
      end if
    end if
    
    opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
    opts % interpolation % interp_mode = 1       ! Set interpolation method

    opts % rt_all % use_q2m             = .FALSE.
    opts % rt_mw % clw_data            = .FALSE. !

    opts % config % apply_reg_limits   = .TRUE.  ! input profiles outside the limits specified in 
                                                 !   the coefficient files are reset to the min/max
    opts % config % verbose            = .TRUE.  ! Enable printing of warnings
    if (.not. l_isMW .and. nml_l_include_cloud) then
      CALL rttov_read_coefs(errorstatus, coefs, opts, &
                            file_coef=coef_filename, file_sccld=cld_coef_filename)
    else
      CALL rttov_read_coefs(errorstatus, coefs, opts, &
                            file_coef=coef_filename)
    end if

    if (l_isMW .and. nml_l_include_cloud) then
      CALL rttov_read_scattcoeffs(errorstatus, opts_scatt, coefs, coefs_scatt,&
                            file_coef=mietable_filename)
    end if
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'fatal error reading coefficients'
      CALL rttov_exit(errorstatus)
    ENDIF

    ! Ensure input number of channels is not higher than number stored in coefficient file
    IF (nchannels > coefs % coef % fmv_chn) THEN
      nchannels = coefs % coef % fmv_chn
    ENDIF

    ! Ensure the options and coefficients are consistent
    CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'error in rttov options'
      CALL rttov_exit(errorstatus)
    ENDIF

    ! --------------------------------------------------------------------------
    ! 3. Allocate RTTOV input and output structures
    ! --------------------------------------------------------------------------

    ! Determine the total number of radiances to simulate (nchanprof).
    ! In this example we simulate all specified channels for each profile, but
    ! in general one can simulate a different number of channels for each profile.

    nchanprof = nchannels * nprof

    ! Allocate structures for rttov_direct
    CALL rttov_alloc_direct( &
          errorstatus,             &
          1_jpim,                  &  ! 1 => allocate
          nprof,                   &
          nchanprof,               &
          nlevels,                 &
          chanprof,                &
          opts,                    &
          profiles,                &
          coefs,                   &
          transmission,            &
          radiance,                &
          calcemis=calcemis,       &
          emissivity=emissivity,   &
          calcrefl=calcrefl,       &
          reflectance=reflectance, &
          init=.TRUE._jplm)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'allocation error for rttov_direct structures'
      CALL rttov_exit(errorstatus)
    ENDIF

    if (l_isMW .and. nml_l_include_cloud) then
      
      ! Allocate the RTTOV-SCATT cloud profiles structure
      ALLOCATE(cld_profiles(nprof), stat=alloc_status)
      IF (alloc_status /= 0) THEN
        WRITE(*,*) 'allocation error for cld_profiles array'
        errorstatus = errorstatus_fatal
        CALL rttov_exit(errorstatus)
      ENDIF

      CALL rttov_alloc_scatt_prof(   &
            errorstatus,             &
            nprof,                   &
            cld_profiles,            &
            nlevels,                 &
            .false.,            &    ! false => separate ciw and snow; true => totalice
            1_jpim,                  &    ! 1 => allocate
            init = .TRUE._jplm,      &
            mmr_snowrain = .true.)  ! snow/rain input units: false => kg/m2/s; true => kg/kg
      IF (errorstatus /= errorstatus_success) THEN
        WRITE(*,*) 'allocation error for cld_profiles structure'
        CALL rttov_exit(errorstatus)
      ENDIF
      
      ! use_chan array is dimensioned by the total number of instrument channels
      ALLOCATE(use_chan(nprof,coefs%coef%fmv_chn), &
               frequencies(nchanprof))

      ! Set use_chan to .TRUE. only for required channels
      use_chan(:,:) = .FALSE._jplm
      DO j = 1, nprof
        use_chan(j,channel_list(1:nchannels)) = .TRUE._jplm
      ENDDO
      ! Populate chanprof and frequencies arrays
      CALL rttov_scatt_setupindex ( &
            nprof,              &
            coefs%coef%fmv_chn, &
            coefs,              &
            nchanprof,          &
            chanprof,           &
            frequencies,        &
            use_chan)
    end if ! l_isMW .and. use_cloud

    ! --------------------------------------------------------------------------
    ! 4. Build the list of profile/channel indices in chanprof
    ! --------------------------------------------------------------------------

    nch = 0_jpim
    DO j = 1, nprof
      DO jch = 1, nchannels
        nch = nch + 1_jpim
        chanprof(nch)%prof = j
        chanprof(nch)%chan = channel_list(jch)
      ENDDO
    ENDDO
    sensor_id_rt(1) = sensor_name
    l_initialized = .true.
  end subroutine rt_rttov_init  

  ! --------------------------------------------------------------------------
  subroutine state_to_rttov(state, x, y, res)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x, y
    type(rt_result), intent(inout) :: res

    integer :: zmax
    zmax = state%nz


    profiles(iprof) % gas_units = 1 ! kg/kg over moist air
    
    ! pressure, temperature, water_vapor

    profiles(iprof) % p(1:zmax) = state%pressure   (x,y,zmax:1:-1) / 100.
    profiles(iprof) % t(1:zmax) = state%temperature(x,y,zmax:1:-1) 
    profiles(iprof) % q(1:zmax) = state%qvapor(x,y,zmax:1:-1) / (1 + state%qvapor(x,y,zmax:1:-1) ) ! convert wapor/dry to wapor/wet



    ! 2-meter air variable
    profiles(iprof) % s2m % t = state%temperature(x,y,1)
    profiles(iprof) % s2m % q = state%qvapor(x,y,1)
    profiles(iprof) % s2m % p = state%Level_Pressure(x,y,0) / 100.
    profiles(iprof) % s2m % u = state%u10(x,y)
    profiles(iprof) % s2m % v = state%v10(x,y)

    ! skin variables
    profiles(iprof) % skin % t = state%tsk(x,y)
    if (state%landmask(x,y) .eq. 1.0) then
      profiles(iprof) % skin % surftype = 0
      profiles(iprof) % skin % fastem(1:5) = (/3.0, 5.0, 15.0, 0.1, 0.3/)
    else
      if (state%icecover(x,y) .gt. 0.0) then
        profiles(iprof) % skin % surftype = 2 ! sea ice
      else
        profiles(iprof) % skin % surftype = 1 ! ocean
      end if
    end if
    
    ! Elevation, latitude and longitude
    profiles(iprof) % elevation = state%hgt(x,y) / 1000.
    profiles(iprof) % latitude  = state%lat(x,y)
    profiles(iprof) % longitude = state%lon(x,y)
    if ( profiles(iprof) % longitude < 0.) profiles(iprof) % longitude = profiles(iprof) % longitude + 360.

    ! Satellite and solar angles
    profiles(iprof) % zenangle = res%zenith_angle(x,y)
!    if (res%azimuth_angle(x,y) == 999.9) then
!      profiles(iprof) % azangle = 0._jprb
!    else
    profiles(iprof) % azangle  = res%azimuth_angle(x,y)
!    end if

    ! -----------------------------------
    ! Specify the cloud input profiles
    ! -----------------------------------

    ! The cloud coefficient files contain optical parameters for 5 pre-defined OPAC liquid
    ! water cloud types, a single liquid water type parameterised in terms of effective diameter,
    ! and for the SSEC/Baum ice cloud dataset. RTTOV also provides a parameterisation of the
    ! Baran ice optical property database.
    !
    ! The cloud inputs are as follows:
    !
    ! Cloud liquid water:
    !
    ! For clw_scheme = 1 you must specify vertical profiles of cloud concentration for one or
    ! more of the 5 OPAC cloud types in profiles(:)%cloud(1:5,1:nlayers).
    !
    ! For clw_scheme = 2 you must specify a vertical profile of CLW concentration (in any of
    ! indices 1-5 in the profiles(:)%cloud(:,1:nlayers) array) and a vertical profile of cloud
    ! liquid water effective diameter in profiles(:)%clwde(1:nlayers).
    !
    ! Ice cloud:
    !
    ! Specify the vertical profile of ice cloud concentration in profiles(:)%cloud(6,1:nlayers).
    !
    ! For ice_scheme = 1 you must either specify a valid ice effective diameter parameterisation
    ! in profiles(:)%idg (1-4) or provide a vertical profile of ice effective diameter in
    ! profiles(:)%icede(1:nlayers).
    !
    ! For ice_scheme = 2 there is no explicit size-dependence so no other inputs are required.
    !
    ! Finally the vertical profile of total layer cloud fraction (for all cloud in the layer)
    ! is specified in profiles(:)%cfrac(1:nlayers).


    ! Select the CLW and ice cloud properties:
    if (nml_l_include_cloud) then
      if (l_isMW) then
        cld_profiles(iprof) % ph  (1:zmax+1) = state%Level_Pressure(x,y,zmax:0:-1) / 100.
   
        cld_profiles(iprof) % cc  (1:zmax) = 1. ! cloud cover (0-1)
        cld_profiles(iprof) % clw (1:zmax) = state%qcloud(x,y,zmax:0:-1)! liquid water (kg/kg)
        cld_profiles(iprof) % ciw (1:zmax) = state%qice  (x,y,zmax:0:-1)! ice water (kg/kg)
        cld_profiles(iprof) % rain(1:zmax) = state%qrain (x,y,zmax:0:-1)! rain (kg/kg)
        cld_profiles(iprof) % sp  (1:zmax) = state%qsnow (x,y,zmax:0:-1) + state%qgraup(x,y,zmax:0:-1)! frozen precip. (kg/kg)
      else
        profiles(iprof) % clw_scheme = 2
        profiles(iprof) % ice_scheme = 1

        ! Set the ice Deff parameterisation to a suitable value
        profiles(:) % idg = 4

        ! Read flag to indicate cloud units (T => kg/kg; F => g/m^3)
        profiles(iprof) % mmr_cldaer = .true.
        ! Read cloud liquid and ice water content, Deff and cloud fraction profiles
        profiles(iprof) % cloud(1,1:zmax) = state%qcloud(x,y,zmax:1:-1)
        profiles(iprof) % cfrac(:) = 1
        profiles(iprof) % clwde(:) = 16.8
        profiles(iprof) % cloud(6,1:zmax) = state%qice  (x,y,zmax:1:-1)
        profiles(iprof) % icede(:) = 25.0
      end if
    end if
    ! --------------------------------------------------------------------------
    ! 6. Specify surface emissivity and reflectance
    ! --------------------------------------------------------------------------

    ! In this example we have no values for input emissivities
    emissivity(:) % emis_in = 0._jprb

    ! Calculate emissivity within RTTOV where the input emissivity value is
    ! zero or less (all channels in this case)
    calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

    ! In this example we have no values for input reflectances
    reflectance(:) % refl_in = 0._jprb

    ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
    ! (all channels in this case)
    calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

    ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
    reflectance(:) % refl_cloud_top = 0._jprb

  end subroutine state_to_rttov

  !--------------------------

  subroutine rt_rttov_destroy(error_status)
    integer, intent(out) :: error_status


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------

  ! Deallocate structures for rttov_direct
  CALL rttov_alloc_direct( &
        errorstatus,             &
        0_jpim,                  &  ! 0 => deallocate
        nprof,                   &
        nchanprof,               &
        nlevels,                 &
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        transmission,            &
        radiance,                &
        calcemis=calcemis,       &
        emissivity=emissivity,   &
        calcrefl=calcrefl,       &
        reflectance=reflectance)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF

  if (l_isMW .and. nml_l_include_cloud) then

    CALL rttov_alloc_scatt_prof(   &
          errorstatus,             &
          nprof,                   &
          cld_profiles,            &
          nlevels,                 &
          .false.,            &  ! must match value used for allocation
          0_jpim)                     ! 0 => deallocate
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'deallocation error for cld_profiles structure'
      CALL rttov_exit(errorstatus)
    ENDIF

    DEALLOCATE(cld_profiles, stat=alloc_status)
    IF (alloc_status /= 0) THEN
      WRITE(*,*) 'dellocation error for cld_profiles array'
    ENDIF

    CALL rttov_dealloc_scattcoeffs(coefs_scatt)
  end if

  l_initialized = .false.
  end subroutine rt_rttov_destroy

  !--------------------------
  subroutine rt_rttov_main(state, res, error_status)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    integer,           intent(out) :: error_status

    integer :: x, y

    if (.not. l_initialized) then
      write(*,*) "rt_rttov not initialized"
      error_status = 1
      stop
    end if
    
    res%nch = nchannels
    res%channels(1:res%nch) = channels(1:res%nch)
    
    do y = state%ybeg, state%yend
      
!      print *, 'y= ',y 
      do x = state%xbeg, state%xend
        ! skip for land columns if land is not calculated
        if ( ( state%landmask(x,y) == 1.0 .and. .not. nml_l_include_land ) .or. &
             ( state%icecover(x,y) == 1.0 .and. .not. nml_l_include_icecover )  ) then
          cycle
        else
          call state_to_rttov(state, x, y, res)
          if (l_isMW .and. nml_l_include_cloud) then
            CALL rttov_scatt ( &
                    errorstatus,         &! out   error flag
                    opts_scatt,          &! in    RTTOV-SCATT options structure
                    nlevels,             &! in    number of profile levels
                    chanprof,            &! in    channel and profile index structure
                    frequencies,         &! in    channel indexes for Mietable lookup
                    profiles,            &! in    profile array
                    cld_profiles,        &! in    cloud/hydrometeor profile array
                    coefs,               &! in    coefficients structure
                    coefs_scatt,         &! in    Mietable structure
                    calcemis,            &! in    flag for internal emissivity calcs
                    emissivity,          &! inout input/output emissivities per channel
                    radiance)             ! inout computed radiances
          else
            CALL rttov_direct(                &
                    errorstatus,              &! out   error flag
                    chanprof,                 &! in    channel and profile index structure
                    opts,                     &! in    options structure
                    profiles,                 &! in    profile array
                    coefs,                    &! in    coefficients structure
                    transmission,             &! inout computed transmittances
                    radiance,                 &! inout computed radiances
                    calcemis    = calcemis,   &! in    flag for internal emissivity calcs
                    emissivity  = emissivity, &! inout input/output emissivities per channel
                    calcrefl    = calcrefl,   &! in    flag for internal BRDF calcs
                    reflectance = reflectance) ! inout input/output BRDFs per channel
          end if 
          
          IF (errorstatus /= errorstatus_success) THEN
            WRITE(*,*) 'error for rttov_direct'
            CALL rttov_exit(errorstatus)
          ENDIF
                
          call rttov_to_result(x, y, res)
        end if
      end do
    end do
    res%sensor_id = sensor_id_rt(1)
    error_status = errorstatus      
  end subroutine rt_rttov_main

  !--------------------------
  subroutine rttov_to_result(x, y, res)

    integer,                          intent(in) :: x,y
    type(rt_result),                  intent(inout) :: res

    res%TB(x,y,:) =   radiance%bt (:) 
    res%emis(x,y,:) = emissivity(:) % emis_out
    
    res%zenith_angle (x,y) = profiles(iprof) % zenangle
    res%azimuth_angle(x,y) = profiles(iprof) % azangle

  end subroutine rttov_to_result
end module rt_rttov_module
