module rt_crtm_module
  use rt_namelist_module
  use rt_state_module
  use rt_result_module
  use rt_reff_module

  use CRTM_Module
  
  implicit none
  PRIVATE
  
  PUBLIC :: state_to_crtm
  PUBLIC :: crtm_to_rt_result
  PUBLIC :: rt_crtm_main
  PUBLIC :: rt_crtm_init
  PUBLIC :: rt_crtm_destroy
  PUBLIC :: rt_crtm_n_channels
  PUBLIC :: rt_crtm_n_sensors
  PUBLIC :: rt_crtm_sensor_id

  INTEGER, PARAMETER :: N_PROFILES  = 1  ! 11934=117*102
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
  INTEGER, PARAMETER :: N_AEROSOLS  = 0

  integer, save :: n_sensors = 0
  character(len=LENGTH_SENSOR_NAME), save :: sensor_id_rt(MAX_SENSORS)
  
  TYPE(CRTM_ChannelInfo_type), save :: chinfo(MAX_SENSORS)
  TYPE(CRTM_Geometry_type)   , save :: geo(N_PROFILES)
  TYPE(CRTM_Atmosphere_type) , save :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)    , save :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), save, allocatable  :: rts(:,:)
  TYPE(CRTM_Options_type), save,    allocatable  :: opt(:,:)

  integer, save :: N_CLOUDS  = -1
  logical, save :: l_initialized = .false.
  real, save, allocatable :: reff_cloud(:)
  real, save, allocatable :: reff_ice(:) 
  real, save, allocatable :: reff_rain(:) 
  real, save, allocatable :: reff_snow(:) 
  real, save, allocatable :: reff_graup(:)

  real, save, allocatable :: rho_air(:)
contains
  function rt_crtm_n_channels( i_sensor )
    integer, intent(in) :: i_sensor
    integer :: rt_crtm_n_channels
    rt_crtm_n_Channels = CRTM_ChannelInfo_n_Channels(chinfo(i_sensor))
  end function rt_crtm_n_channels

  ! ---------
  function rt_crtm_n_sensors()
    integer :: rt_crtm_n_sensors
    rt_crtm_n_sensors = n_sensors
  end function rt_crtm_n_sensors
  ! ---------
  function rt_crtm_sensor_id(i_sensor)
    integer, intent(in) :: i_sensor
    character(len=LENGTH_SENSOR_NAME) :: rt_crtm_sensor_id
    rt_crtm_sensor_id = sensor_id_rt(i_sensor)
  end function rt_crtm_sensor_id

  ! ---------
  subroutine state_to_crtm(state, x, y, res)
    type(model_state), intent(in) :: state
    integer,           intent(in) :: x, y
    type(rt_result), intent(inout) :: res

    integer :: zmax
    integer :: zc ! z in CRTM  (top-to-bottom)
    integer :: zm ! z in model (bottom-to-top)
    integer :: icl

    real :: reff_ratio
    real :: q2c  ! factor to convert Q to water_content

    real(fp) :: zenith_angle, azimuth_angle

    reff_ratio = 1.0

    zmax = state%nz

    atm(1)%Climatology         = TROPICAL
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)

    atm(1)%Level_Pressure(0:zmax)   = state%level_pressure(x,y,zmax:0:-1) / 100.0
    atm(1)%Pressure      (1:zmax)   = state%pressure      (x,y,zmax:1:-1) / 100.0
    atm(1)%Temperature   (1:zmax)   = state%temperature   (x,y,zmax:1:-1)
    atm(1)%Absorber      (1:zmax,1) = state%qvapor        (x,y,zmax:1:-1) * 1000.0
    atm(1)%Absorber      (1:zmax,2) = 0. !5.0E-2
    
    ! set cloud back to zero
    if (nml_l_include_cloud) then
      do icl = 1, N_CLOUDS
        atm(1)%Cloud(icl)%Water_Content = 0.0
        atm(1)%Cloud(icl)%Effective_Radius = 0.0
      end do

      ! set cloud hydrometer information
      rho_air = state%pressure(x,y,:)/287.2/(state%temperature(x,y,:)+0.61*(state%qvapor(x,y,:)/(1+state%qvapor(x,y,:))))
      call calc_reff_driver(nml_s_reff_method, state, x, y, rho_air, reff_cloud, reff_ice, reff_rain, reff_snow, reff_graup)
      ! remove negative effective radius
      where(reff_cloud .lt. 0.0)  reff_cloud =0.0
      where(reff_ice .lt. 0.0)    reff_ice =0.0
      where(reff_rain .lt. 0.0)   reff_rain =0.0
      where(reff_snow .lt. 0.0)   reff_snow =0.0
      where(reff_graup .lt. 0.0)  reff_graup =0.0

      !print *, 'rho_air', rho_air
!      print *, 'pressure', state%pressure(x,y,:)
!      print *, 'temperature', state%temperature(x,y,:)
!      print *, 'delz', state%delz(x,y,:)
      do zc = 1, zmax
        zm = zmax - zc + 1
        q2c = rho_air(zm)*state%delz(x,y,zm)
        icl = 1
        atm(1)%Cloud(icl)%Type = WATER_CLOUD
        if (state%qcloud(x,y,zm) .gt. 0.0) then
          atm(1)%Cloud(icl)%Effective_Radius(zc) = reff_cloud(zm)
          atm(1)%Cloud(icl)%Water_Content(zc)    = state%qcloud(x,y,zm) * q2c 
        end if
        
        icl = 2
        atm(1)%Cloud(icl)%Type = RAIN_CLOUD
        if (state%qrain(x,y,zm) .gt. 0.0) then
          atm(1)%Cloud(icl)%Effective_Radius(zc) = reff_rain(zm)
          atm(1)%Cloud(icl)%Water_Content(zc)    = state%qrain(x,y,zm) * q2c 
        end if
        
        icl = 3
        atm(1)%Cloud(icl)%Type = ICE_CLOUD
        if (state%qice(x,y,zm) .gt. 0.0) then
          atm(1)%Cloud(icl)%Effective_Radius(zc) = reff_ice(zm)
          atm(1)%Cloud(icl)%Water_Content(zc)    = state%qice(x,y,zm) * q2c 
        end if
        
        icl = 4
        atm(1)%Cloud(icl)%Type = SNOW_CLOUD
        if (state%qsnow(x,y,zm) .gt. 0.0) then
          atm(1)%Cloud(icl)%Effective_Radius(zc) = reff_snow(zm)
          atm(1)%Cloud(icl)%Water_Content(zc)    = state%qsnow(x,y,zm) * q2c 
        end if
        
        icl = 5
        atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
        if (state%qgraup(x,y,zm) .gt. 0.0) then
          atm(1)%Cloud(icl)%Effective_Radius(zc) = reff_graup(zm)
          atm(1)%Cloud(icl)%Water_Content(zc)    = state%qgraup(x,y,zm) * q2c 
        end if

      end do

      ! save the reff values if necessary
      if (nml_l_output_reff) then
        res%reff_cloud(x,y,1:zmax) = atm(1)%Cloud(1)%Effective_radius(zmax:1:-1)
        res%reff_rain (x,y,1:zmax) = atm(1)%Cloud(2)%Effective_radius(zmax:1:-1)
        res%reff_ice  (x,y,1:zmax) = atm(1)%Cloud(3)%Effective_radius(zmax:1:-1)
        res%reff_snow (x,y,1:zmax) = atm(1)%Cloud(4)%Effective_radius(zmax:1:-1)
        res%reff_graup(x,y,1:zmax) = atm(1)%Cloud(5)%Effective_radius(zmax:1:-1)
      end if
    end if 
    ! surface
    if(state%landmask(x,y).eq.1.0) then
      sfc(1)%Water_Coverage = 0.0_fp
      sfc(1)%Land_Coverage  = 1.0_fp
      sfc(1)%Ice_Coverage   = 0.0_fp
      sfc(1)%Land_Temperature = state%tsk(x,y)
      sfc(1)%Soil_Temperature = state%tsk(x,y)
    else if (state%icecover(x,y) .gt. 0.0) then
      sfc(1)%Water_Coverage = 0.0_fp
      sfc(1)%Land_Coverage  = 0.0_fp
      sfc(1)%Ice_Coverage   = 1.0_fp
      sfc(1)%Ice_Type = 1  ! Sea ice
      sfc(1)%Ice_Temperature = state%tsk(x,y)
    else
      sfc(1)%Water_Coverage = 1.0_fp
      sfc(1)%Land_Coverage  = 0.0_fp
      sfc(1)%Ice_Coverage   = 0.0_fp
      sfc(1)%Water_Type = 1  ! Sea water
      sfc(1)%Water_Temperature = state%tsk(x,y)
    endif
  
    ! geometry 
    zenith_angle = res%zenith_angle(x,y)
    azimuth_angle = res%azimuth_angle(x,y)
    call CRTM_Geometry_SetValue( geo, &
                sensor_zenith_angle  = zenith_angle, &
                sensor_azimuth_angle = azimuth_angle )
 

!    call CRTM_Atmosphere_Inspect(atm(1))
!    call crtm_surface_inspect(sfc(1)) 
!    print *, 'qgraup', state%qgraup(x,y,:)

!    stop
  end subroutine state_to_crtm

  ! --------------
  subroutine crtm_to_rt_result(res, x, y)
    integer,                    intent(in) :: x,y
    type(rt_result), intent(inout) :: res

    res%TB  (x,y,1:res%nch)        = rts(1:res%nch,1)%Brightness_Temperature
    res%emis(x,y,1:res%nch)        = rts(1:res%nch,1)%surface_Emissivity
!    res%zenith_angle (x,y) = geo(  1)%Sensor_Zenith_Angle
!    res%azimuth_angle(x,y) = geo(  1)%Sensor_Azimuth_Angle
  end subroutine crtm_to_rt_result
  
  ! --------------

  subroutine rt_crtm_init(sensor_id, n_layers, error_status, channels)
    character(len=*),  intent(in) :: sensor_id(:)
    integer,           intent(in) :: n_layers
    integer,           intent(out):: error_status
    integer,           intent(in), optional :: channels(:,:) ! channel x sensor

    integer :: i_sensor, i, max_n_channels, n_channels
    character(len=LENGTH_SENSOR_NAME) :: sensor_id_crtm(MAX_SENSORS)
    integer :: channels_used(MAX_CHANNELS)
    
    TYPE(CRTM_ChannelInfo_type), save :: chinfo2(MAX_SENSORS)
    sensor_id_crtm(:) = ''
    n_sensors = 0
    do i = 1, MAX_SENSORS
      if (sensor_id(i) == '' ) then
        cycle
      else if ( sensor_id(i) == 'gmi_gpm_hf' .or. sensor_id(i) == 'gmi_gpm_lf' ) then
        n_sensors = n_sensors + 1
        sensor_id_crtm(n_sensors) = 'gmi_gpm'
        sensor_id_rt  (n_sensors) = sensor_id(i)
      else
        n_sensors = n_sensors + 1
        sensor_id_crtm(n_sensors) = sensor_id(i)
        sensor_id_rt  (n_sensors) = sensor_id(i)
      end if
    end do

    if (my_proc_id==0) print *, 'n_sensor=',n_sensors,'sensor_id=',sensor_id_crtm
    allocate( opt(n_sensors, n_profiles), &
              stat=error_status )
    call MPI_Barrier( comm, ierr )
    Error_Status = CRTM_Init( Sensor_Id_crtm(1:n_sensors), &  ! Input... must be an array, hencethe (/../)
                              chinfo(1:n_sensors) , &  ! Output
                              IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                              IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                              !MWwaterCoeff_File = 'FASTEM6.MWwater.EmisCoeff.bin', &
                              File_Path='coefficients/', &
                              CloudCoeff_File_rain  = nml_s_crtm_rainLUT,&
                              CloudCoeff_File_snow  = nml_s_crtm_snowLUT,&
                              CloudCoeff_File_graup = nml_s_crtm_graupelLUT,&
                              quiet=.true. )
!    if (nml_i_nchannels > 0) then                       
!      Error_Status = CRTM_ChannelInfo_Subset( chinfo(1),Channel_Subset=nml_a_channels(1:nml_i_nchannels) )
!    end if
    !print *, 'my_proc_id', my_proc_id, Error_status
    call MPI_Barrier( comm, ierr )
    
    if ( present(channels)) then
      do i_sensor = 1, n_sensors
        channels_used(:) = 0
        n_channels = 0
        do i = 1, MAX_CHANNELS
          if ( channels(i,i_sensor) > 0 ) then
            n_channels = n_channels + 1
            channels_used(n_channels) = channels(i, i_sensor)
          end if
        end do
        if (n_channels > 0) then
          !print *, 'my_proc_id', my_proc_id, channels_used(1:n_channels)
          Error_Status = CRTM_ChannelInfo_Subset( chinfo(i_sensor),Channel_Subset=channels_used(1:n_channels) )
        end if
      end do
    else
      do i_sensor = 1, n_sensors
        if (sensor_id_rt(i_sensor) == 'gmi_gpm_hf') then
          Error_Status = CRTM_ChannelInfo_Subset( chinfo(i_sensor),Channel_Subset=(/10,11,12,13/) )
        else if (sensor_id_rt(i_sensor) == 'gmi_gpm_lf') then
          Error_Status = CRTM_ChannelInfo_Subset( chinfo(i_sensor),Channel_Subset=(/1,2,3,4,5,6,7,8,9/) )
        end if
      end do
    end if
    
    max_n_Channels = MAXVAL(CRTM_ChannelInfo_n_Channels(chinfo(1:n_sensors)))
    
!    nml_i_nchannels = n_Channels
!    nml_a_channels(1:nml_i_nchannels) = CRTM_ChannelInfo_Channels( chinfo(i_sensor) )

    if (nml_l_include_cloud) then
      N_CLOUDS = 5
    else
      N_CLOUDS = 0
    end if

    ! allocate CRTM structures
    ALLOCATE( rts( max_n_Channels, N_PROFILES ), STAT=error_Status ) 
    call CRTM_RTSolution_create( rts, n_layers)
    CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS,  N_AEROSOLS)
    do i_sensor = 1, n_sensors
      n_channels = CRTM_ChannelInfo_n_Channels(chinfo(i_sensor))
      call CRTM_Options_Create( opt(i_sensor,:), N_Channels )
    end do

    allocate( reff_cloud(n_layers),&
              reff_ice  (n_layers),&
              reff_rain (n_layers),&
              reff_snow (n_layers),&
              reff_graup(n_layers),&
              rho_air   (n_layers),&
              STAT=error_status)
    l_initialized = .true.
    opt%RT_Algorithm_ID = RT_SOI
    
    if (nml_i_nstream > 0) then
      opt%Use_N_Streams = .TRUE.
      opt%n_Streams = nml_i_nstream
    else if (.not. nml_s_crtm_rainLUT=='' .or. .not. nml_s_crtm_snowLUT=='' .or. .not. nml_s_crtm_graupelLUT=='' ) then
      opt%Use_N_Streams = .TRUE.
      opt%n_Streams = 16
    end if
  end subroutine rt_crtm_init
  ! ---------------
  subroutine rt_crtm_main(i_sensor, state, res, error_status)
    integer,           intent(in) :: i_sensor
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    integer,           intent(out) :: error_status

    integer :: x, y, n_channels

    if (.not. l_initialized) then
      write(*,*) "rt_crtm not initialized"
      error_status = 1
      stop
    end if
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(i_sensor))
  
    res%sensor_id = sensor_id_rt(i_sensor)
    res%channels = 0
    res%channels(1:res%nch) = CRTM_ChannelInfo_Channels(chinfo(i_sensor) )
    !print *, CRTM_ChannelInfo_Channels(chinfo(i_sensor) )
    !print *, res%channels(1:res%nch) 
!    call rt_result_alloc_from_state( res, state, n_channels, error_status, .true.)
!    call rt_result_GEO_zenith(res, state, -89.5, error_status)

    !do y = 100, 200
    do y = state%ybeg, state%yend
      if (should_print(10) .and. mod(y, 10) == 0) print *, 'y: ', y
      !do x = 120, 180
      do x = state%xbeg, state%xend
        
        if ( ( state%landmask(x,y) == 1.0 .and. .not. nml_l_include_land ) .or. &
             ( state%icecover(x,y) == 1.0 .and. .not. nml_l_include_icecover )  ) then
          cycle
        end if
      
        call CRTM_Atmosphere_zero( atm )
        call CRTM_Surface_Zero( sfc )
        call CRTM_RTSolution_Zero( rts)
        
        ! put values into CRTM structures
        call state_to_crtm(state, x, y, res)

        Error_Status = CRTM_Forward( atm        , &
                                     sfc        , &
                                     geo, &
                                     chinfo(i_sensor:i_sensor), &
                                     rts(1:n_channels,:), &
                                     Options = opt(i_sensor,:) )
        ! put result in rt_result
        call crtm_to_rt_result(res, x, y) 
      end do ! x
    end do  ! y
  end subroutine rt_crtm_main

  ! --------------

  subroutine rt_crtm_destroy(error_status)
    integer, intent(out) :: error_status

    error_status = CRTM_Destroy( chinfo )
    deallocate( reff_cloud, &
                reff_ice,   &
                reff_rain,  &
                reff_snow,  &
                reff_graup, &
                rho_air, &
                STAT=error_status)
    l_initialized = .false.


  end subroutine rt_crtm_destroy




end module rt_crtm_module
