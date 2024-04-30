module rt_zenith_module
  use rt_state_module
  use mapinfo_define
  use obs_define
  use rt_result_module
  use rt_constant_module
  use rt_util_module
  implicit none

contains
  ! --------------------
  ! subroutines to calculate zenith angle
  ! --------------------
  subroutine calc_zenith_GEO(state, res, sensor_id, error_status)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    character(len=*),  intent(in)    :: sensor_id 
    integer, intent(out) :: error_status

    real, parameter :: Re     = 6378000.0
    real, parameter :: sat_h = 35780000.0

    real, allocatable :: sat_dis(:,:)
    real, allocatable :: scan_angle(:,:)
    real :: sat_lon

    !sat_lon = sat_lon_deg *deg2rad

    if (sensor_id == 'abi_gr') then
      sat_lon = -89.5 * deg2rad
    else if (sensor_id == 'abi_g16') then
      sat_lon = -75.196 * deg2rad
    else if (sensor_id == 'ahi_himawari8') then
      sat_lon = 140. * deg2rad
    else
      write(*,*) "Error in calc_zenith_GEO()"
      write(*,*) "sensor_id ", sensor_id, " not recognized"
      error_status = FAILURE
      return
    end if 
    
    if ( state%l_allocated == .false. ) then
      write(*,*) "Error in calc_zenith_GEO()"
      write(*,*) "state is not initialized"
      error_status = FAILURE
      return
    end if

    allocate( sat_dis(    state%xbeg:state%xend, state%ybeg:state%yend), &
              scan_angle( state%xbeg:state%xend, state%ybeg:state%yend), &
              stat=error_status)
    sat_dis = sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(state%lon*deg2rad-sat_lon)*cos(state%lat*deg2rad))
    SCAN_ANGLE = rad2deg*asin(Re/sat_dis*sqrt(1-(cos(state%lon*deg2rad-sat_lon)*cos(state%lat*deg2rad))**2))
    res%ZENITH_ANGLE = SCAN_ANGLE+rad2deg*acos(cos(state%lon*deg2rad-sat_lon)*cos(state%lat*deg2rad))
    if (nml_s_rt_program == 'crtm') then
      res%AZIMUTH_ANGLE = DEFAULT_AZIMUTH_CRTM
    else if (nml_s_rt_program == 'rttov') then
      res%AZIMUTH_ANGLE = DEFAULT_AZIMUTH_RTTOV
    end if
    deallocate( sat_dis, scan_angle )
  end subroutine calc_zenith_GEO 
  ! --------------

  !------------------------------------------------------
  ! 
  !------------------------------------------------------ 
  subroutine calc_zenith_obs(state, res, obs, sensor, error_status, method)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    type(obs_type), intent(in) :: obs
    character(len=*), intent(in) :: sensor
    integer,           intent(out)   :: error_status
    character(len=*), intent(in), optional :: method

    logical, allocatable :: idx_used(:)
    integer :: i

    allocate( idx_used(obs%num), stat=error_status )
    do i = 1, obs%num
    if ( obs%type(i) .eq. 'Microwave ' .and. trim(adjustl(obs%sat(i))) .eq. trim(adjustl(sensor))) then
        idx_used(i) = .true.
      else
        idx_used(i) = .false.
      end if
    end do
    if (.not. present(method) .or. method == 'gaussian') then
      call calc_zenith_obs_gaussian(state, res, sensor, obs%num, idx_used, obs%position(:,1), obs%position(:,2), obs%zenith_angle, obs%azimuth_angle, &
                        obs%efov_aScan, obs%efov_cScan, error_status)
    else if ( method == 'boxsearch' ) then
      call calc_zenith_obs_boxsearch(state, res, sensor, obs%num, idx_used, obs%position(:,1), obs%position(:,2), obs%zenith_angle, obs%azimuth_angle, &
                        obs%efov_aScan, obs%efov_cScan, error_status)
    else
      print *, "zenith calculation method not recognized:", method
      stop
    end if
    deallocate( idx_used)
  end subroutine calc_zenith_obs
  ! --------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------ 
  subroutine calc_zenith_obs_microwave(state, res, obs_mw, sensor, error_status, method)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    type(microwave_data_type), intent(in) :: obs_mw
    character(len=*), intent(in) :: sensor
    integer,           intent(out)   :: error_status
    character(len=*), intent(in), optional :: method

    logical, allocatable :: idx_used(:)
    integer :: i

    allocate( idx_used(obs_mw%num), stat=error_status )
    do i = 1, obs_mw%num
      if ( trim(adjustl(obs_mw%platform(i))) .eq. trim(adjustl(sensor))) then
        idx_used(i) = .true.
      else
        idx_used(i) = .false.
      end if
    end do
    if (.not. present(method) .or. method == 'gaussian') then
      call calc_zenith_obs_gaussian(state, res, sensor, obs_mw%num, idx_used, obs_mw%ii, obs_mw%jj, obs_mw%zenith_angle, obs_mw%azimuth_angle, &
                        obs_mw%efov_aScan, obs_mw%efov_cScan, error_status)
    else if ( method == 'boxsearch' ) then
      call calc_zenith_obs_boxsearch(state, res, sensor, obs_mw%num, idx_used, obs_mw%ii, obs_mw%jj, obs_mw%zenith_angle, obs_mw%azimuth_angle, &
                        obs_mw%efov_aScan, obs_mw%efov_cScan, error_status)
    else
      print *, "zenith calculation method not recognized:", method
      stop
    end if
    deallocate( idx_used)
  end subroutine calc_zenith_obs_microwave
  ! --------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------ 
  subroutine calc_zenith_obs_gaussian(state, res, sensor, n, idx_used, obs_x, obs_y, obs_zen, obs_azi, efov_a, efov_c, error_status)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(inout) :: res
    character(len=*),  intent(in)    :: sensor
    integer,           intent(in)    :: n
    logical,           intent(in)    :: idx_used(n)
    real,              intent(in)    :: obs_x(n)
    real,              intent(in)    :: obs_y(n)
    real,              intent(in)    :: obs_zen(n)
    real,              intent(in)    :: obs_azi(n)
    real,              intent(in)    :: efov_a(n) ! effective FOV along-track in units of model grids (e.g., 2.5 for 10km FOV in 4-km grid)
    real,              intent(in)    :: efov_c(n) ! ~ cross-track
    integer,           intent(out)   :: error_status

    integer :: iobs
    real, allocatable :: w(:,:) ! weight
    real, allocatable :: total_w(:,:) ! weight
    real, allocatable :: sindazi(:,:), cosdazi(:,:) ! used to handle angle-wrap around 0
    integer :: xmin, xmax, ymin, ymax ! region that needs to calculate weights
    real :: max_distance 
    real :: obs_ii, obs_jj

    allocate( w      ( state%xbeg:state%xend, state%ybeg:state%yend), &
              total_w( state%xbeg:state%xend, state%ybeg:state%yend), &
              sindazi( state%xbeg:state%xend, state%ybeg:state%yend), &
              cosdazi( state%xbeg:state%xend, state%ybeg:state%yend), &
              stat=error_status)
    total_w = 0.
    sindazi = 0.
    cosdazi = 0.
    res%zenith_angle = 0.
    res%azimuth_angle= 0.

    ! loop through observations
    do iobs = 1, n
      if (idx_used(iobs) ) then
        obs_ii = obs_x(iobs)
        obs_jj = obs_y(iobs)
        max_distance = (max(efov_a(iobs), efov_c(iobs))/state%dx * 2)

        xmin = int(max(state%xbeg,   floor(obs_ii) - ceiling(max_distance)))
        xmax = int(min(state%xend, ceiling(obs_ii) + ceiling(max_distance)))
        ymin = int(max(state%ybeg,   floor(obs_jj) - ceiling(max_distance)))
        ymax = int(min(state%yend, ceiling(obs_jj) + ceiling(max_distance)))
        if (xmin .le. xmax .and. ymin .le. ymax) then
          w = 0
          call gaussian_ellipse_weight(obs_ii, obs_jj, xmax-xmin+1, ymax-ymin+1, &
                    state%x(xmin:xmax,ymin:ymax), state%y(xmin:xmax,ymin:ymax), &
                    efov_a(iobs), efov_c(iobs), state%dx, obs_azi(iobs), &
                    w(xmin:xmax,ymin:ymax) )
          res%zenith_angle (xmin:xmax,ymin:ymax) =   res%zenith_angle (xmin:xmax,ymin:ymax) + obs_zen(iobs) * w(xmin:xmax,ymin:ymax)
          sindazi            (xmin:xmax,ymin:ymax) = sindazi          (xmin:xmax,ymin:ymax) + sind(obs_azi(iobs)) * w(xmin:xmax,ymin:ymax)
          cosdazi            (xmin:xmax,ymin:ymax) = cosdazi          (xmin:xmax,ymin:ymax) + cosd(obs_azi(iobs)) * w(xmin:xmax,ymin:ymax)
          total_w            (xmin:xmax,ymin:ymax) = total_w          (xmin:xmax,ymin:ymax) +                 w(xmin:xmax,ymin:ymax)
        end if
      end if
    end do
    where( total_w .ne. 0.) 
      res%zenith_angle  = res%zenith_angle / total_w
      sindazi = sindazi / total_w
      cosdazi = cosdazi / total_w
      res%azimuth_angle = atan2d(sindazi, cosdazi)
    end where
    where (res%azimuth_angle < 0.) 
      res%azimuth_angle = res%azimuth_angle + 360.
    end where
    deallocate(w, total_w, sindazi, cosdazi)
  end subroutine calc_zenith_obs_gaussian
  ! --------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------ 
  subroutine calc_zenith_obs_boxsearch(state, res, sensor, n, idx_used, obs_x, obs_y, obs_zen, obs_azi, efov_a, efov_c, error_status)
    type(model_state), intent(in) :: state
    type(rt_result), intent(inout) :: res 
    character(len=*),  intent(in)    :: sensor
    integer,           intent(in)    :: n
    logical,           intent(in)    :: idx_used(n)
    real,              intent(in)    :: obs_x(n)
    real,              intent(in)    :: obs_y(n)
    real,              intent(in)    :: obs_zen(n)
    real,              intent(in)    :: obs_azi(n)
    real,              intent(in)    :: efov_a(n) ! effective FOV along-track in units of model grids (e.g., 2.5 for 10km FOV in 4-km grid)
    real,              intent(in)    :: efov_c(n) ! ~ cross-track
    integer,           intent(out)   :: error_status


    real, allocatable :: all_zenith_angle   (:,:)
    real, allocatable :: all_azimuth_sind   (:,:)
    real, allocatable :: all_azimuth_cosd   (:,:)
    real, allocatable :: all_angle_sum_count(:,:)
!    real, allocatable :: my_azimuth_sind   (:,:)
!    real, allocatable :: my_azimuth_cosd   (:,:)

    integer :: iobs, obs_beg, obs_end
    real :: sigma, search_radius
    integer :: search_radius_int
    integer :: search_xmin
    integer :: search_xmax
    integer :: search_ymin
    integer :: search_ymax
    integer :: ix, iy
    integer :: mycount
    real    :: mycount_sum

    real :: my_zenith_angle_sum
    real :: my_azimuth_sind_sum
    real :: my_azimuth_cosd_sum
    
    ! need to operate on the whole domain
    allocate( all_zenith_angle   (state%xbeg0:state%xend0, state%ybeg0:state%yend0), &
              all_azimuth_sind   (state%xbeg0:state%xend0, state%ybeg0:state%yend0), &
              all_azimuth_cosd   (state%xbeg0:state%xend0, state%ybeg0:state%yend0), &
              all_angle_sum_count(state%xbeg0:state%xend0, state%ybeg0:state%yend0), &
!              my_azimuth_sind    (state%xbeg :state%xend , state%yend :state%yend ), &
!              my_azimuth_cosd    (state%xbeg :state%xend , state%yend :state%yend ), &
              stat = error_status )
    all_zenith_angle    = 0.0
    all_azimuth_sind    = 0.0
    all_azimuth_cosd    = 0.0
    all_angle_sum_count = 0.0

    ! loop through observations
    call calc_proc_helper( nprocs, my_proc_id, 1, n, obs_beg, obs_end)
    do iobs = obs_beg, obs_end
      if (idx_used(iobs)) then
        sigma = (0.5 * ( ( ( ( efov_a(iobs) + efov_c(iobs) ) /2 ) / 1.18) ) )
        search_radius = 1 + (2 * sigma)

        search_xmin = max(state%xbeg0, nint(obs_x(iobs)) - ceiling(1+(search_radius / state%dx)))
        search_xmax = min(state%xend0, nint(obs_x(iobs)) + ceiling(1+(search_radius / state%dx)))
        search_ymin = max(state%ybeg0, nint(obs_y(iobs)) - ceiling(1+(search_radius / state%dx)))
        search_ymax = min(state%yend0, nint(obs_y(iobs)) + ceiling(1+(search_radius / state%dx)))

        all_zenith_angle(    search_xmin:search_xmax, search_ymin:search_ymax) = all_zenith_angle(    search_xmin:search_xmax, search_ymin:search_ymax) + abs(obs_zen(iobs) )
        all_azimuth_sind(    search_xmin:search_xmax, search_ymin:search_ymax) = all_azimuth_sind(    search_xmin:search_xmax, search_ymin:search_ymax) + sind(obs_azi(iobs))
        all_azimuth_cosd(    search_xmin:search_xmax, search_ymin:search_ymax) = all_azimuth_cosd(    search_xmin:search_xmax, search_ymin:search_ymax) + cosd(obs_azi(iobs))
        all_angle_sum_count( search_xmin:search_xmax, search_ymin:search_ymax) = all_angle_sum_count( search_xmin:search_xmax, search_ymin:search_ymax) + 1
      end if
    end do ! iobs
    ! send and recv from other processors
    call MPI_Allreduce(MPI_IN_PLACE, all_zenith_angle,    (state%xend0-state%xbeg0+1)*(state%yend0-state%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
    call MPI_Allreduce(MPI_IN_PLACE, all_azimuth_sind,    (state%xend0-state%xbeg0+1)*(state%yend0-state%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
    call MPI_Allreduce(MPI_IN_PLACE, all_azimuth_cosd,    (state%xend0-state%xbeg0+1)*(state%yend0-state%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
    call MPI_Allreduce(MPI_IN_PLACE, all_angle_sum_count, (state%xend0-state%xbeg0+1)*(state%yend0-state%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)

    ! take care of my own scan and zenith angle
    do ix = state%xbeg, state%xend
      do iy = state%ybeg, state%yend
        if (all_angle_sum_count(ix,iy) > 0) then
          res%zenith_angle(ix,iy) = all_zenith_angle(ix,iy) / all_angle_sum_count(ix,iy)
          res%azimuth_angle(ix,iy) = atan2d( all_azimuth_sind(ix,iy)/all_angle_sum_count(ix,iy), all_azimuth_cosd(ix,iy)/all_angle_sum_count(ix,iy) )
        else ! need to increase search radius
          search_radius_int = 0
          mycount = 0
          do while(mycount < 2 .and. state%dx*search_radius_int <= 200) !maximum search radius 200 km; maybe better to use instrument-specific value
            search_radius_int = search_radius_int + 1
            search_xmin = max(state%xbeg0, ix-search_radius_int)
            search_xmax = min(state%xend0, ix+search_radius_int)
            search_ymin = max(state%ybeg0, iy-search_radius_int)
            search_ymax = min(state%yend0, iy+search_radius_int)
            mycount = count( all_angle_sum_count(search_xmin:search_xmax, search_ymin:search_ymax) > 0 )
          end do
          mycount_sum                        = sum(all_angle_sum_count(search_xmin:search_xmax, search_ymin:search_ymax))
          if (mycount_sum > 0) then
            res%zenith_angle(ix,iy) =          sum(all_zenith_angle(search_xmin:search_xmax, search_ymin:search_ymax)) / mycount_sum
            res%azimuth_angle(ix,iy) = atan2d( sum(all_azimuth_sind(search_xmin:search_xmax, search_ymin:search_ymax)) / mycount_sum, &
                                                 sum(all_azimuth_cosd(search_xmin:search_xmax, search_ymin:search_ymax)) / mycount_sum )
          else
            res%zenith_angle(ix,iy) =  get_default_sensor_zenith(sensor)
            res%azimuth_angle(ix,iy) = 0.
          end if
        end if
      end do
    end do
    ! angle wrap
    where (res%azimuth_angle < 0.) 
      res%azimuth_angle = res%azimuth_angle + 360.
    end where

  end subroutine calc_zenith_obs_boxsearch
  ! ----------


  !------------------------------------------------------
  ! calculate 'true' model columns when viewing in a slant direction 
  !------------------------------------------------------ 
  subroutine calc_slant_path(state, res, slant, error_status)
    type(model_state), intent(in) :: state
    type(rt_result),   intent(in) :: res
    type(model_state), intent(inout) :: slant
    integer, intent(out) :: error_status
  
    real :: xoff, yoff, weight
    integer :: x1, y1   ! profile loop
    integer :: x2, y2 ! inner loop
    integer :: x_min, x_max, y_min, y_max
    integer :: z1, z2
    real, allocatable :: alt1(:), alt2(:), pdiff(:)
    real :: altitude

    allocate( alt1(0:state%nz), alt2(0:state%nz), pdiff(state%nz-1), stat=error_status )

    
    if (slant%l_allocated) call model_state_dealloc( slant, error_status ) 
    ! allocate without buffer
    call model_state_alloc(slant, state%xbeg0, state%xend0, state%ybeg0, state%yend0, &
                           state%xbeg, state%xend, state%ybeg, state%yend, state%nz, error_status)
    ! copy 2-d fields
    slant%lat     (state%xbeg:state%xend,state%ybeg:state%yend) = state%lat     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%lon     (state%xbeg:state%xend,state%ybeg:state%yend) = state%lon     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%x       (state%xbeg:state%xend,state%ybeg:state%yend) = state%x       (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%y       (state%xbeg:state%xend,state%ybeg:state%yend) = state%y       (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%tsk     (state%xbeg:state%xend,state%ybeg:state%yend) = state%tsk     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%u10     (state%xbeg:state%xend,state%ybeg:state%yend) = state%u10     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%v10     (state%xbeg:state%xend,state%ybeg:state%yend) = state%v10     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%landmask(state%xbeg:state%xend,state%ybeg:state%yend) = state%landmask(state%xbeg:state%xend,state%ybeg:state%yend)
    slant%icecover(state%xbeg:state%xend,state%ybeg:state%yend) = state%icecover(state%xbeg:state%xend,state%ybeg:state%yend)
    slant%lakemask(state%xbeg:state%xend,state%ybeg:state%yend) = state%lakemask(state%xbeg:state%xend,state%ybeg:state%yend)
    slant%soiltype(state%xbeg:state%xend,state%ybeg:state%yend) = state%soiltype(state%xbeg:state%xend,state%ybeg:state%yend)
    slant%hgt     (state%xbeg:state%xend,state%ybeg:state%yend) = state%hgt     (state%xbeg:state%xend,state%ybeg:state%yend)
    slant%level_pressure(state%xbeg:state%xend,state%ybeg:state%yend,0) = state%level_pressure(state%xbeg:state%xend,state%ybeg:state%yend,0)
    
    ! copy other properties
    slant%proj = state%proj
    slant%dx   = state%dx
    slant%mp_physics_name = state%mp_physics_name
  
    ! deal with 3-d fields
    do x1 = state%xbeg, state%xend
      do y1 = state%ybeg, state%yend
        alt1(0) = state%hgt(x1,y1)
        do z1 = 1, state%nz
          alt1(z1) = alt1(z1-1) + state%delz(x1,y1,z1)
          ! altitude of layer z1 at FOV is 
          altitude = ( alt1(z1-1) + alt1(z1) )/2.
          
          yoff = (alt1(z1)+alt1(z1-1))/2 * tand( res%zenith_angle(x1,y1) ) * cosd( res%azimuth_angle(x1,y1)) / state%dx / 1000
          xoff = (alt1(z1)+alt1(z1-1))/2 * tand( res%zenith_angle(x1,y1) ) * sind( res%azimuth_angle(x1,y1)) / state%dx / 1000

          x_min = x1 + int(floor( xoff ))
          x_max = x1 + int(floor( xoff )) + 1
          y_min = y1 + int(floor( yoff ))
          y_max = y1 + int(floor( yoff )) + 1
          
          ! Since we divide whole domain into multiple subdomains
          ! need to check if the buffer size is large enough
          if (x_min<state%xbeg1 ) then
            if (state%xbeg1 == state%xbeg0) then
              x_min = state%xbeg1 ! first sub-domain
            else
              print *, "buffer size not large enough for slant calculation (1)"
              print *, x_min, state%xbeg1, state%xbeg0
              stop
            end if
          end if
          if (x_max<state%xbeg1 ) then
            if (state%xbeg1 == state%xbeg0) then
              x_max = state%xbeg1 ! first sub-domain
            else
              print *, "buffer size not large enough for slant calculation (2)"
              print *, x_max, state%xbeg1, state%xbeg0
              stop
            end if
          end if
          if (x_min>state%xend1 ) then
            if (state%xend1 == state%xend0) then
              x_min = state%xend1 ! last sub-domain
            else
              print *, "buffer size not large enough for slant calculation (3)"
              print *, x_min, state%xend1, state%xend0
              stop
            end if
          end if
          if (x_max>state%xend1 ) then
            if (state%xend1 == state%xend0) then
              x_max = state%xend1 ! last sub-domain
            else
              print *, "buffer size not large enough for slant calculation (4)"
              print *, x_max, state%xend1, state%xend0
              stop
            end if
          end if
          if (y_min<state%ybeg1 ) then
            if (state%ybeg1 == state%ybeg0) then
              y_min = state%ybeg1 ! first sub-domain
            else
              print *, "buffer size not large enough for slant calculation (5)"
              print *, y_min, state%ybeg1, state%ybeg0
              stop
            end if
          end if
          if (y_max<state%ybeg1 ) then
            if (state%ybeg1 == state%ybeg0) then
              y_max = state%ybeg1 ! first sub-domain
            else
              print *, "buffer size not large enough for slant calculation (6)"
              print *, y_max, state%ybeg1, state%ybeg0
              stop
            end if
          end if
          if (y_min>state%yend1 ) then
            if (state%yend1 == state%yend0) then
              y_min = state%yend1 ! last sub-domain
            else
              print *, "buffer size not large enough for slant calculation (7)"
              print *, y_min, state%yend1, state%yend0
              stop
            end if
          end if
          if (y_max>state%yend1 ) then
            if (state%yend1 == state%yend0) then
              y_max = state%yend1 ! last sub-domain
            else
              print *, "buffer size not large enough for slant calculation (8)"
              print *, y_max, state%yend1, state%yend0, state%yend
              print *, res%zenith_angle(x1,y1), res%azimuth_angle(x1,y1), yoff
              print *, z1, alt1(z1)
              stop
            end if
          end if
          
          ! calculate averaged fields along slant path
          do x2 = x_min, x_max
            do y2 = y_min, y_max
              weight = 1
              if (x_min .ne. x_max) weight = weight * abs( x1+xoff -x2 )
              if (y_min .ne. y_max) weight = weight * abs( y1+yoff -y2 )
              
              alt2(0) = state%hgt(x2,y2)
              do z2 = 1, state%nz
                alt2(z2) = alt2(z2-1) + state%delz(x2,y2,z2)
              end do ! z2
        
              ! if altitude exceeds min/max alt2 at (x2,y2), use profile at (x1,y1)
              ! for example, close to steep terrain slope
              if (altitude < alt2(0) .or. altitude > alt2(state%nz) ) then
                call model_state_add_weight(state, x1, y1, z1, slant, x1, y1, z1, weight)

              else 
                ! search for correct layer at (x2, y2)
                do z2 = 1, state%nz
                  if( altitude < alt2(z2) .and. altitude > alt2(z2-1) ) then
                    exit
                  end if
                end do
                call model_state_add_weight(state, x2, y2, z2, slant, x1, y1, z1, weight)
              end if
            end do ! y2
          end do ! x2
        end do ! z1
        ! re-calculate level_pressure from pressure
        ! top layer
        slant%level_Pressure(x1,y1,slant%nz) = state%pressure(x1,y1,slant%nz)*1.5 - state%pressure(x1,y1,slant%nz-1)*0.5 
        ! lowest layer
        ! we do not need to calculate lowest layer. In principle, it should be
        ! the same as state%level_pressure(x1,y1,0)
        ! slant%level_pressure(x1,y1,0) =  max(state%level_pressure(x1,y1,0), state%pressure(x1,y1,1)*1.5 - state%pressure(x1,y1,2)*0.5)
        do z1 = 2, slant%nz
          slant%level_pressure(x1,y1,z1-1) = (state%pressure(x1,y1,z1) + state%pressure(x1,y1,z1-1) ) * 0.5
        end do ! z

        ! check again to see if pressure levels make sense
        pdiff(1:slant%nz-1) = slant%pressure(x1,y1,1:slant%nz-1) - slant%pressure(x1,y1,2:slant%nz)
        if (ANY(pdiff <= 0.99) ) then
          slant%delz          (x1,y1,:) = state%delz          (x1,y1,:)
          slant%level_pressure(x1,y1,:) = state%level_pressure(x1,y1,:)
          slant%pressure      (x1,y1,:) = state%pressure      (x1,y1,:)
          slant%temperature   (x1,y1,:) = state%temperature   (x1,y1,:)
          slant%qvapor        (x1,y1,:) = state%qvapor        (x1,y1,:)

          if (nml_l_include_cloud) then
            slant%qcloud(x1,y1,:) =  state%qcloud(x1,y1,:) 
            slant%qrain (x1,y1,:) =  state%qrain (x1,y1,:) 
            slant%qice  (x1,y1,:) =  state%qice  (x1,y1,:) 
            slant%qsnow (x1,y1,:) =  state%qsnow (x1,y1,:) 
            slant%qgraup(x1,y1,:) =  state%qgraup(x1,y1,:) 
            slant%nrain (x1,y1,:) =  state%nrain (x1,y1,:) 
          end if
        end if ! end check
      end do ! y1
    end do ! x1

    deallocate( alt1, alt2, pdiff, stat=error_status )
  end subroutine calc_slant_path






end module rt_zenith_module
