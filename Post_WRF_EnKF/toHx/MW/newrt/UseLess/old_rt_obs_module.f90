module rt_obs_module
  use rt_state_module
  use mapinfo_define
  use obs_define
  use rt_constant_module
  
  implicit none

  integer, parameter :: default_nsensors = MAX_SENSORS
  integer, parameter :: default_nchannels = MAX_CHANNELS

  type obs_info_type
    integer :: nsensors = -1
    character(len=LENGTH_SENSOR_NAME), allocatable :: sensors(:)
    integer, allocatable :: nchannels(:)
    integer, allocatable :: channels(:,:)
    
    logical :: l_initialized = .false.
    integer :: max_nsensors  = -1
    integer :: max_nchannels = -1
  end type obs_info_type

contains
  subroutine obs_info_print(obs_info)
    type(obs_info_type), intent(in) :: obs_info

    integer :: i
    print *,'nsensors', obs_info%nsensors
    print *,'max_nsensors',  obs_info%max_nsensors
    print *,'max_nchannels', obs_info%max_nchannels
    print *,'sensors    ', obs_info%sensors
    print *,'nchannels', obs_info%nchannels
    do i = 1, obs_info%nsensors
    print *,'channels',  obs_info%channels(:,i)
    end do
    print *, '-------------------------'
  end subroutine obs_info_print

  subroutine obs_info_init(obs_info)
    type(obs_info_type), intent(inout) :: obs_info
    integer :: isensor
    
    if (.not. obs_info%l_initialized) then
      obs_info%max_nsensors = default_nsensors
      obs_info%max_nchannels= default_nchannels
      allocate( obs_info%sensors  ( obs_info%max_nsensors ) )
      allocate( obs_info%nchannels( obs_info%max_nsensors ) )
      allocate( obs_info%channels ( obs_info%max_nchannels, obs_info%max_nsensors ) )
    end if
    obs_info%nsensors = 0
    do isensor = 1, obs_info%max_nsensors
      obs_info%sensors(isensor) = ''
    end do
    obs_info%nchannels = 0
    obs_info%channels  = 0
  end subroutine obs_info_init
  ! ---------------------
  subroutine obs_info_grow(obs_info, new_max_nsensors, new_max_nchannels)
    type(obs_info_type), intent(inout) :: obs_info
    integer, intent(in) :: new_max_nsensors
    integer, intent(in) :: new_max_nchannels

    character(len=LENGTH_SENSOR_NAME), allocatable :: tmp_sensors(:)
    integer, allocatable :: tmp_nchannels(:)
    integer, allocatable :: tmp_channels(:,:) 
    integer :: old_max_nsensors
    integer :: old_max_nchannels

    if (new_max_nsensors < obs_info%max_nsensors .or. new_max_nchannels < obs_info%max_nchannels) then
      print *, 'sizes should be greater than original'
      stop
    end if

    old_max_nsensors  = obs_info%max_nsensors
    old_max_nchannels = obs_info%max_nchannels
    allocate( tmp_sensors  ( old_max_nsensors ) )
    allocate( tmp_nchannels( old_max_nsensors ) )
    allocate( tmp_channels ( old_max_nsensors, old_max_nchannels) )
    
    tmp_sensors   = obs_info%sensors
    tmp_nchannels = obs_info%nchannels
    tmp_channels  = obs_info%channels
    
    ! resize the arrays of obs_info
    deallocate( obs_info%sensors, &
                obs_info%nchannels, &
                obs_info%channels )
    
    obs_info%max_nsensors = new_max_nsensors
    obs_info%max_nchannels = new_max_nchannels
    allocate( obs_info%sensors  ( new_max_nsensors ) )
    allocate( obs_info%nchannels( new_max_nsensors ) )
    allocate( obs_info%channels ( new_max_nchannels, new_max_nsensors) )
    obs_info%sensors = ''
    obs_info%nchannels = 0
    obs_info%channels = 0

    obs_info%sensors  (1:old_max_nsensors) = tmp_sensors(1:old_max_nsensors)
    obs_info%nchannels(1:old_max_nsensors) = tmp_nchannels(1:old_max_nsensors)
    obs_info%channels (1:old_max_nchannels,1:old_max_nsensors)= tmp_channels(1:old_max_nchannels,1:old_max_nsensors)

  end subroutine obs_info_grow
  ! ---------------------
  subroutine obs_info_add(obs_info, sensor, channel)
    type(obs_info_type), intent(inout) :: obs_info
    character(len=*), intent(in) :: sensor
    integer, intent(in) :: channel

    logical :: found
    integer :: i, j
    ! exit if sensor or channel is actually empty
    if (sensor == '' .or. channel <=0) return
    ! look for sensor
    i = 1
    do while( i <= obs_info%nsensors .and. sensor > obs_info%sensors(i))
      i = i + 1
    end do
    print *, 'i', i
    ! insert sensor if not found
    if (i > obs_info%nsensors .or. sensor < obs_info%sensors(i) ) then
      print *, 'add new sensor "', sensor, '" at location ', i
      
      if (i >= obs_info%max_nsensors ) then
        call obs_info_grow( obs_info, obs_info%max_nsensors+5, obs_info%max_nchannels)
      end if
      
      if (i <= obs_info%nsensors) then
        obs_info%sensors  (i+1:obs_info%nsensors+1)    = obs_info%sensors  (i:obs_info%nsensors)
        obs_info%nchannels(i+1:obs_info%nsensors+1)    = obs_info%nchannels(i:obs_info%nsensors)
        obs_info%channels (:,i+1:obs_info%nsensors+1) = obs_info%channels (:,i:obs_info%nsensors)
      end if
      obs_info%sensors  (i)   = sensor
      obs_info%nchannels(i)   = 1
      obs_info%channels (1,i) = channel
      obs_info%nsensors = obs_info%nsensors + 1
      print *, 'obs_info%channels',obs_info%channels
    end if
    ! look for channel
    j = 1
    do while( j <= obs_info%nchannels(i) .and. channel > obs_info%channels(j,i))
      j = j + 1
    end do
    print *, 'j', j
    ! insert channel if not found
    if (j > obs_info%nchannels(i) .or. channel < obs_info%channels(j,i) ) then

      if (j >= obs_info%max_nchannels) then
        call obs_info_grow( obs_info, obs_info%max_nsensors, obs_info%max_nchannels+5)
      end if
      
      if (j <= obs_info%nchannels(i) ) then
        obs_info%channels (j+1:obs_info%nchannels(i)+1,i) = obs_info%channels (j:obs_info%nchannels(i),i)
      end if
      obs_info%channels(j,i) = channel
      print *, 'j,i',j,i
      print *, 'obs_info%channels(j,i)',obs_info%channels(j,i)
      obs_info%nchannels(i) = obs_info%nchannels(i) + 1
    end if
    

  end subroutine obs_info_add

  ! ---------------------
  subroutine get_obs_info(obs_info, obs_mw)
    type(obs_info_type), intent(inout) :: obs_info
    type(microwave_data_type), intent(in) :: obs_mw

    integer :: isensor
    integer :: iobs


    do iobs = 1, obs_mw%num
      ! search 
!      if (mod(iobs, 100) == 0) then
!        call obs_info_print(obs_info)
!      end if
      call obs_info_add(obs_info, trim(adjustl(obs_mw%platform(iobs))), int(obs_mw%ch(iobs)))
    end do

  end subroutine get_obs_info

  subroutine get_obs_info_fromfile(obs_info, filename, error_status)
    type(obs_info_type), intent(inout) :: obs_info
    character(len=*), intent(in) :: filename
    integer, intent(out) :: error_status
    
    character (len=12)  :: so_time
    character (len=16)  :: sat_id
    integer             :: ch_info
    integer             :: iost
    open (10, file=filename, status='old', form='formatted', iostat=iost)

!...... get the data number
    do_get_raw_data_loop_read : do
      read(10, '(a12,a16,i12)', iostat = iost ) so_time, sat_id, ch_info
      if( iost .ne. 0 ) exit
      call obs_info_add(obs_info, trim(adjustl(sat_id)), ch_info)
    end do do_get_raw_data_loop_read

  end subroutine get_obs_info_fromfile

  subroutine get_microwave(filename, state, microwave, error_status)
    character(len=*), intent(in) :: filename
    type(model_state), intent(in) :: state
    type(microwave_data_type), intent(inout) :: microwave 
    integer, intent(out) :: error_status
    
    character (len=12)                  :: so_time
    character (len=16)                  :: sat_id
    integer                             :: i, n, iost, num
    integer                             :: ch_info,hroi_rad,hroi_drad
    real                                :: lat, lon, tb, err
    real                                :: s, h, rx, ry, ir1, jr1, is, js
    real                                :: efov_aScan, efov_cScan, scan_angle, zenith_angle, azimuth_angle
    real                                :: sat_lat, sat_lon, sat_alt
    open (10, file=filename, status='old', form='formatted', iostat=iost)

!...... get the data number
    num = 0
    do_get_raw_data_loop_read : do
      read(10, '(a12,a16,i12,3f12.3)', iostat = iost ) so_time, sat_id, ch_info, lat, lon, tb
      if( iost .ne. 0 ) exit
      num = num + 1
    end do do_get_raw_data_loop_read

!...... allocate
    allocate( microwave%lat( num ) )
    allocate( microwave%lon( num ) )
    allocate( microwave%platform( num ) )
    allocate( microwave%ch( num ) )
    allocate( microwave%ii( num ) )
    allocate( microwave%jj( num ) )
    allocate( microwave%tb( num ) )
    allocate( microwave%hroi( num ) )
    allocate( microwave%hroi_d( num ) )
    allocate( microwave%err( num ) )
    allocate( microwave%efov_aScan( num ) )
    allocate( microwave%efov_cScan( num ) )
    allocate( microwave%scan_angle( num ) )
    allocate( microwave%zenith_angle( num ) )
    allocate( microwave%azimuth_angle( num ) )
    allocate( microwave%sat_lat( num ) )
    allocate( microwave%sat_lon( num ) )
    allocate( microwave%sat_alt( num ) )

!...... get data
    rewind(10)
    num = 0
    do_get_raw_data_loop : do

!......... Enkf with odd data, and verify with even data
    read(10, '(a12,a16,i12,3f12.3,2i12,9f12.3)', iostat = iost ) &
         so_time, sat_id, ch_info, lat, lon, tb, hroi_rad, hroi_drad, err, &
         efov_aScan, efov_cScan, scan_angle, zenith_angle, azimuth_angle, sat_lat, sat_lon, sat_alt
    if( iost .ne. 0 ) exit
!......... calculate radar center's position according to wrf domain grid
    call latlon_to_ij( state%proj, lat, lon, is, js )
!......... evaluate
    if ( is > state%xbeg0 .and. is < state%xend0 .and. js > state%ybeg0 .and. js < state%yend0) then
      num = num + 1
      microwave%lat(num) = lat
      microwave%lon(num) = lon
      microwave%platform(num) = sat_id
      microwave%ch(num) = ch_info
      microwave%ii(num) = is
      microwave%jj(num) = js
      microwave%tb(num) = tb
      microwave%hroi(num) = (hroi_rad*1000)/state%proj%dx
      microwave%hroi_d(num) = (hroi_drad*1000)/state%proj%dx
      microwave%err(num) = err
      microwave%efov_aScan(num) = efov_aScan
      microwave%efov_cScan(num) = efov_cScan
      microwave%scan_angle(num) = scan_angle
      microwave%zenith_angle(num) = zenith_angle
      microwave%azimuth_angle(num) = azimuth_angle
      microwave%sat_lat(num) = sat_lat
      microwave%sat_lon(num) = sat_lon
      microwave%sat_alt(num) = sat_alt
    else
    endif

    end do do_get_raw_data_loop
    microwave%num = num
    close (10)




  end subroutine get_microwave


end module rt_obs_module
