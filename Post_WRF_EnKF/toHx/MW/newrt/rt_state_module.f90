module rt_state_module
  use rt_constant_module
  use netcdf
  use mpi_module
  use rt_namelist_module
  use rt_util_module
  use mapinfo_define
  use wrf_tools
  
  use wgrib2api
  implicit none

  type model_state
    ! my state
    logical :: l_allocated = .false.
    type(proj_info) :: proj
    ! dimensions
    integer :: xbeg0 = -1 ! whole domain
    integer :: xend0 = -1
    integer :: ybeg0 = -1
    integer :: yend0 = -1
    integer :: xbeg1 = -1 ! only this processor + buffer
    integer :: xend1 = -1 ! used for slant path
    integer :: ybeg1 = -1
    integer :: yend1 = -1
    integer :: xbeg = -1  ! only this processor without buffer
    integer :: xend = -1  ! radiative transfer performed within these ranges
    integer :: ybeg = -1
    integer :: yend = -1
    integer :: nz = -1
    real    :: dx = -1 ! grid spacing in km at truelats
    
    ! 2-D var
    real, allocatable :: lat(:,:)
    real, allocatable :: lon(:,:)
    real, allocatable :: x(:,:)
    real, allocatable :: y(:,:)
    
    real, allocatable :: tsk(:,:)
    real, allocatable :: u10(:,:)
    real, allocatable :: v10(:,:)
    real, allocatable :: landmask(:,:)
    real, allocatable :: icecover(:,:)
    real, allocatable :: lakemask(:,:)
    real, allocatable :: soiltype(:,:)
    real, allocatable :: hgt(:,:)

    ! 3-D var, z-direction from bottom to top
    real, allocatable :: delz  (:,:,:) ! nx*ny*nz
    real, allocatable :: level_pressure(:,:,:) ! nx*ny*(nz+1)
    real, allocatable :: pressure(:,:,:) ! nx*ny*nz
    real, allocatable :: temperature   (:,:,:) ! nx*ny*nz
    real, allocatable :: qvapor(:,:,:)
    real, allocatable :: qcloud(:,:,:)
    real, allocatable :: qrain (:,:,:)
    real, allocatable :: qice  (:,:,:)
    real, allocatable :: qsnow (:,:,:)
    real, allocatable :: qgraup(:,:,:)
    real, allocatable :: nrain (:,:,:)
    real, allocatable :: nsnow (:,:,:)
    real, allocatable :: ngraup (:,:,:)

    character(len=20) :: mp_physics_name = ''
  end type

contains

  !------------------------------------------------------
  ! Allocate model state (state: structure type)
  !------------------------------------------------------
  subroutine model_state_alloc(state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, nz, error_status, buffer_in)
    type(model_state), intent(inout) :: state
    integer, intent(in) :: xbeg0, xend0, ybeg0, yend0, nz
    integer, intent(in) :: xbeg, xend, ybeg, yend
    integer, intent(out) :: error_status
    integer, intent(in), optional :: buffer_in

    integer :: xbeg1, xend1, ybeg1, yend1
    integer :: buffer
    if (present(buffer_in)) then
      buffer = buffer_in
    else
      buffer = 0
    end if


    if (should_print(1)) then
      write(*,*) '--- Initializing state ---'
      write(*,*) 'xbeg, xend: ', xbeg, xend
      write(*,*) 'ybeg, yend: ', ybeg, yend
      write(*,*) 'nz: ', nz
      write(*,*) '--- Initializing state finish ---'
    end if

    if (state%l_allocated .eq. .true.) then
      write(*,*) 'Error in state_module.f90 model_state_alloc():'
      write(*,*) 'structure already allocated'
      error_status = FAILURE
      return
    end if

    state%xbeg0 = xbeg0
    state%xend0 = xend0
    state%ybeg0 = ybeg0
    state%yend0 = yend0
    state%xbeg = xbeg
    state%xend = xend
    state%ybeg = ybeg
    state%yend = yend
    
    state%xbeg1 = xbeg - buffer
    state%xend1 = xend + buffer
    state%ybeg1 = ybeg - buffer
    state%yend1 = yend + buffer

    if (state%xbeg1 < xbeg0) state%xbeg1 = xbeg0
    if (state%xend1 > xend0) state%xend1 = xend0
    if (state%ybeg1 < ybeg0) state%ybeg1 = ybeg0
    if (state%yend1 > yend0) state%yend1 = yend0

    xbeg1 = state%xbeg1
    xend1 = state%xend1
    ybeg1 = state%ybeg1
    yend1 = state%yend1
    
    state%nz   = nz
    allocate( &
              state%lat(xbeg1:xend1, ybeg1:yend1), &
              state%lon(xbeg1:xend1, ybeg1:yend1), &
              state%x  (xbeg1:xend1, ybeg1:yend1), &
              state%y  (xbeg1:xend1, ybeg1:yend1), &
              state%tsk(xbeg1:xend1, ybeg1:yend1), &
              state%u10(xbeg1:xend1, ybeg1:yend1), &
              state%v10(xbeg1:xend1, ybeg1:yend1), &
              state%landmask(xbeg1:xend1, ybeg1:yend1), &
              state%icecover(xbeg1:xend1, ybeg1:yend1), &
              state%lakemask(xbeg1:xend1, ybeg1:yend1), &
              state%soiltype(xbeg1:xend1, ybeg1:yend1), &
              state%hgt     (xbeg1:xend1, ybeg1:yend1), &
              state%delz    (xbeg1:xend1,ybeg1:yend1,nz),&
              state%level_pressure(xbeg1:xend1,ybeg1:yend1,0:nz),&
              state%pressure      (xbeg1:xend1,ybeg1:yend1,nz),&
              state%temperature   (xbeg1:xend1,ybeg1:yend1,nz),&
              state%qvapor        (xbeg1:xend1,ybeg1:yend1,nz),&
              stat=error_status )

    state%lat(:, :) = 0.0
    state%lon(:, :) = 0.0
    state%x  (:, :) = 0.0
    state%y  (:, :) = 0.0
    state%tsk(:, :) = 0.0
    state%u10(:, :) = 0.0
    state%v10(:, :) = 0.0
    state%landmask(:, :) = 0.0
    state%icecover(:, :) = 0.0
    state%lakemask(:, :) = 0.0
    state%soiltype(:, :) = 0.0
    state%hgt     (:, :) = 0.0
    state%delz    (:,:,:) = 0.0
    state%level_pressure(:,:,:) = 0.0
    state%pressure      (:,:,:) = 0.0
    state%temperature   (:,:,:) = 0.0
    state%qvapor        (:,:,:) = 0.0
   
    ! cloudy sky
    if (nml_l_include_cloud) then
      allocate( &
              state%qcloud(xbeg1:xend1,ybeg1:yend1,nz),&
              state%qrain (xbeg1:xend1,ybeg1:yend1,nz),&
              state%qice  (xbeg1:xend1,ybeg1:yend1,nz),&
              state%qsnow (xbeg1:xend1,ybeg1:yend1,nz),&
              state%qgraup(xbeg1:xend1,ybeg1:yend1,nz),&
              state%nrain (xbeg1:xend1,ybeg1:yend1,nz),&
              state%nsnow (xbeg1:xend1,ybeg1:yend1,nz),&
              state%ngraup (xbeg1:xend1,ybeg1:yend1,nz),&
              stat=error_status )
      state%qcloud(:,:,:) = 0.0
      state%qrain (:,:,:) = 0.0
      state%qice  (:,:,:) = 0.0
      state%qsnow (:,:,:) = 0.0
      state%qgraup(:,:,:) = 0.0
      state%nrain (:,:,:) = 0.0
      state%nsnow (:,:,:) = 0.0
      state%ngraup (:,:,:) = 0.0
    end if
  
  state%l_allocated = .true.        
  end subroutine model_state_alloc
  !---------------


  !------------------------------------------------------
  ! Deallocate model state (state: structure type)
  !------------------------------------------------------
  subroutine model_state_dealloc(state, error_status)
    type(model_state), intent(inout) :: state
    integer, intent(out) :: error_status
    
    if (state%l_allocated .eq. .false.) then
      write(*,*) 'Error in state_module.f90 model_state_dealloc():'
      write(*,*) 'structure not allocated'
      error_status = FAILURE
      return
    end if

    deallocate( &
              state%lat, &
              state%lon, &
              state%x,   &
              state%y,   &
              state%tsk, &
              state%u10, &
              state%v10, &
              state%landmask, &
              state%icecover, &
              state%lakemask, &
              state%soiltype, &
              state%hgt,      &
              state%delz    ,&
              state%level_pressure,&
              state%pressure,&
              state%temperature   ,&
              state%qvapor,&
              stat=error_status )

    ! cloudy sky              
    if (nml_l_include_cloud) then
      deallocate( &
              state%qcloud,&
              state%qrain ,&
              state%qice  ,&
              state%qsnow ,&
              state%qgraup,&
              state%nrain ,&
              state%nsnow ,&
              state%ngraup ,&
              stat=error_status )
    end if

    state%l_allocated = .false.
    state%xbeg0 = -1
    state%xend0 = -1
    state%ybeg0 = -1
    state%yend0 = -1
    state%xbeg1 = -1
    state%xend1 = -1
    state%ybeg1 = -1
    state%yend1 = -1
    state%xbeg = -1
    state%xend = -1
    state%ybeg = -1
    state%yend = -1
    state%nz   = -1

  end subroutine model_state_dealloc
  !---------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------
  subroutine calc_proc_xyrange_simple(xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend )
    integer, intent(in) :: xbeg0, xend0, ybeg0, yend0
    integer, intent(out) :: xbeg, xend, ybeg, yend
      
    call calc_proc_helper(nml_i_nicpu, global_my_iid, xbeg0, xend0, xbeg, xend)
    call calc_proc_helper(nml_i_njcpu, global_my_jid, ybeg0, yend0, ybeg, yend)

  end subroutine calc_proc_xyrange_simple
  !---------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------
  subroutine calc_proc_helper(nycpu, my_id, ybeg0, yend0, ybeg, yend)
    integer, intent(in) :: nycpu, my_id
    integer, intent(in) :: ybeg0, yend0
    integer, intent(out) :: ybeg, yend
    integer :: ny, modyi, nyi

    ! 1 -- y_calc_beg  -- p_y_beg -- p_y_end -- y_calc_end -- ymax
    ny = yend0-ybeg0+1
    modyi = mod(ny, nycpu)
    if (my_id < modyi) then
        nyi = ny/nycpu + 1
        ybeg = ybeg0 + my_id * nyi
        yend = ybeg + nyi -1
    else
        nyi = ny/nycpu
        ybeg = ybeg0 + my_id * nyi + modyi
        yend = ybeg + nyi -1
    end if
  end subroutine calc_proc_helper
  !---------------


  !------------------------------------------------------
  ! Read WRF model output into memory 
  !------------------------------------------------------ 
  subroutine model_state_read_wrf(state, filename_wrf, time_idx, error_status, buffer)
    type(model_state), intent(inout) :: state
    character(len=*), intent(in) :: filename_wrf
    integer, intent(in) :: time_idx
    integer, intent(out) :: error_status
    integer, intent(in), optional :: buffer
    
    integer :: xmax, ymax, zmax
    integer :: xbeg0, xend0, ybeg0, yend0
    integer :: xbeg1, xend1, ybeg1, yend1
    integer :: xbeg, xend, ybeg, yend
    integer :: fid

    REAL, PARAMETER :: P1000MB=100000.D0
    REAL, PARAMETER :: R_D=287.D0
    REAL, PARAMETER :: CP=7.D0*R_D/2.D0
    ! temporary variables
    real, allocatable :: p   (:,:,:)
    real, allocatable :: pb  (:,:,:)
    real, allocatable :: ph  (:,:,:)
    real, allocatable :: phb (:,:,:)
    real, allocatable :: t   (:,:,:)
    real, allocatable :: psfc(:,:)

    integer :: z
    integer :: mp_physics
    integer :: ix, iy

    call getxyzmax_wrfout(filename_wrf, xmax, ymax, zmax)
    if (should_print(1)) then
      write(*,*) 'xmax,ymax,zmax: ', xmax, ymax, zmax
    end if

    ! calculate x,y ranges according to namelist
    nml_i_x_calc_beg = max(1,    nml_i_x_calc_beg)
    nml_i_x_calc_end = min(xmax, nml_i_x_calc_end)
    nml_i_y_calc_beg = max(1,    nml_i_y_calc_beg)
    nml_i_y_calc_end = min(ymax, nml_i_y_calc_end)

    xbeg0 = nml_i_x_calc_beg
    xend0 = nml_i_x_calc_end
    ybeg0 = nml_i_y_calc_beg
    yend0 = nml_i_y_calc_end

    if (should_print(1)) then
      write(*,*) 'xbeg0,xend0,ybeg0,yend0: ', xbeg0,xend0,ybeg0,yend0
    end if
    
    ! 
    call calc_proc_xyrange_simple(xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend )
    if (should_print(10)) then
      write(*,*) 'xbeg,xend,ybeg,yend: ', xbeg,xend,ybeg,yend
    end if

    ! initialize state
    if (present(buffer)) then
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status, buffer)
    else
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status)
    end if

    ! allocate temporary arrays
    xbeg1 = state%xbeg1
    xend1 = state%xend1
    ybeg1 = state%ybeg1
    yend1 = state%yend1
    allocate( p   (xbeg1:xend1, ybeg1:yend1, zmax), &
              pb  (xbeg1:xend1, ybeg1:yend1, zmax), &
              ph  (xbeg1:xend1, ybeg1:yend1, zmax+1), &
              phb (xbeg1:xend1, ybeg1:yend1, zmax+1), &
              t   (xbeg1:xend1, ybeg1:yend1, zmax), &
              psfc(xbeg1:xend1, ybeg1:yend1), &
              stat=error_status)

    ! read model output 
    call open_file( filename_wrf, nf_nowrite, fid)
    call get_variable2d_local(fid,'XLAT',    xbeg1, xend1, ybeg1, yend1,           time_idx, state%lat)
    call get_variable2d_local(fid,'XLONG',   xbeg1, xend1, ybeg1, yend1,           time_idx, state%lon)
    call get_variable3d_local(fid,'P',       xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, p)
    call get_variable3d_local(fid,'PB',      xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, pb)
    call get_variable3d_local(fid,'PH',      xbeg1, xend1, ybeg1, yend1, 1, zmax+1,time_idx, ph)
    call get_variable3d_local(fid,'PHB',     xbeg1, xend1, ybeg1, yend1, 1, zmax+1,time_idx, phb)
    call get_variable3d_local(fid,'T',       xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, t)
    call get_variable3d_local(fid,'QVAPOR',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qvapor)
    call get_variable2d_local(fid,'PSFC',    xbeg1, xend1, ybeg1, yend1,           time_idx, psfc)
    call get_variable2d_local(fid,'TSK',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%tsk)
    call get_variable2d_local(fid,'HGT',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%hgt)
    call get_variable2d_local(fid,'LANDMASK',xbeg1, xend1, ybeg1, yend1,           time_idx, state%landmask)
    call get_variable2d_local(fid,'LAKEMASK',xbeg1, xend1, ybeg1, yend1,           time_idx, state%lakemask)
    call get_variable2d_local(fid,'U10',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%u10)
    call get_variable2d_local(fid,'V10',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%v10)
   
    ! read the microphysics scheme used in WRF
    error_status = nf_get_att_int(fid, nf_global, 'MP_PHYSICS', mp_physics)
    if (mp_physics == 6) then
      state%mp_physics_name = 'WSM6'
    else if (mp_physics == 7) then
      state%mp_physics_name = 'Goddard'
    else if (mp_physics == 8) then
      state%mp_physics_name = 'Thompson08'
    else if (mp_physics == 10) then
      state%mp_physics_name = 'Morrison'
    else if (mp_physics == 17) then
      state%mp_physics_name = 'NSSL'
    else if (mp_physics == 28) then
      state%mp_physics_name = 'Thompson08'
    end if

    ! cloudy sky
    if (nml_l_include_cloud) then
      call get_variable3d_local(fid,'QCLOUD',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qcloud)
      call get_variable3d_local(fid,'QRAIN',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qrain)
      call get_variable3d_local(fid,'QICE',    xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qice)
      call get_variable3d_local(fid,'QSNOW',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qsnow)
      call get_variable3d_local(fid,'QGRAUP',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qgraup)
      if (state%mp_physics_name == 'Thompson08') then
        call get_variable3d_local(fid,'QNRAIN',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%nrain)
        !print *, state%nrain
      else if (state%mp_physics_name == 'NSSL') then
        call get_variable3d_local(fid,'QNRAIN',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%nrain)
        call get_variable3d_local(fid,'QNSNOW',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%nsnow)
        call get_variable3d_local(fid,'QNGRAUPEL',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%ngraup)
      end if
    end if

    ! close the file
    call close_file(fid)

    ! Remove negative hydro
    where(state%qvapor .lt. 1.0e-10) state%qvapor =1.0e-10
    if (nml_l_include_cloud) then
      where(state%qcloud .lt. 0.0) state%qcloud =0.0
      where(state%qice   .lt. 0.0) state%qice   =0.0
      where(state%qrain  .lt. 0.0) state%qrain  =0.0
      where(state%qsnow  .lt. 0.0) state%qsnow  =0.0
      where(state%qgraup .lt. 0.0) state%qgraup =0.0
    end if

    ! ***************************
    ! calculate other variables not directly saved in WRF output
    ! ***************************
    state%pressure = P + PB
    state%temperature = (T + 300.0) * ( (state%pressure / P1000MB) ** (R_D/CP) )
    ! top layer
    state%level_Pressure(:,:,zmax) = state%pressure(:,:,zmax)*1.5 - state%pressure(:,:,zmax-1)*0.5 
    ! lowest layer
    state%delz          (:,:,1) = (PH(:,:,2) + PHB(:,:,2)) / 9.806 - state%hgt(:,:)
    state%level_pressure(:,:,0) =  max(psfc, state%pressure(:,:,1)*1.5 - state%pressure(:,:,2)*0.5)
    do z = 2, zmax
      state%delz(:,:,z)             = ((PH(:,:,z+1) + PHB(:,:,z+1))-(PH(:,:,z) + PHB(:,:,z)))/9.806
      state%level_pressure(:,:,z-1) = (state%pressure(:,:,z) + state%pressure(:,:,z-1) ) * 0.5
    end do ! z
    ! *******
 
    ! set projection
    call set_domain_proj(filename_wrf, state%proj)
    state%dx = state%proj%dx / 1000.

    do ix = state%xbeg1, state%xend1
      do iy = state%ybeg1, state%yend1
        state%x(ix,iy) = ix
        state%y(ix,iy) = iy
      end do
    end do

  end subroutine model_state_read_wrf
  !---------------


  !------------------------------------------------------
  ! Read FV3 restart model output into memory 
  !------------------------------------------------------ 
  subroutine model_state_read_fv3restart(state, &
                                         filename_f3r_core, &
                                         filename_f3r_ak, &
                                         filename_f3r_tracer, &
                                         filename_f3r_phy, &
                                         filename_f3r_sfc, &
                                         filename_f3r_grid, &
                                         time_idx, error_status, buffer)
    !-rw-r--r-- 1 yuz31 G-801549 6.4G Apr 12 08:18 fv_core.res.nest02.tile7.nc
    !-rw-r--r-- 1 yuz31 G-801549  35M Apr 12 08:18 fv_srf_wnd.res.nest02.tile7.nc
    !-rw-r--r-- 1 yuz31 G-801549 8.5G Apr 12 08:18 fv_tracer.res.nest02.tile7.nc
    !-rw-r--r-- 1 yuz31 G-801549 8.7G Apr 12 08:18 phy_data.nest02.tile7.nc
    !-rw-r--r-- 1 yuz31 G-801549 1.5G Apr 12 08:18 sfc_data.nest02.tile7.nc
    
    type(model_state), intent(inout) :: state
    character(len=*), intent(in) :: filename_f3r_core
    character(len=*), intent(in) :: filename_f3r_ak
    character(len=*), intent(in) :: filename_f3r_tracer
    character(len=*), intent(in) :: filename_f3r_phy
    character(len=*), intent(in) :: filename_f3r_sfc
    character(len=*), intent(in) :: filename_f3r_grid
    integer, intent(in) :: time_idx
    integer, intent(out) :: error_status
    integer, intent(in), optional :: buffer
    
    integer :: xmax, ymax, zmax
    integer :: xbeg0, xend0, ybeg0, yend0
    integer :: xbeg1, xend1, ybeg1, yend1
    integer :: xbeg, xend, ybeg, yend
    integer :: fid

    REAL, PARAMETER :: P1000MB=100000.D0
    REAL, PARAMETER :: R_D=287.D0
    REAL, PARAMETER :: CP=7.D0*R_D/2.D0
    ! temporary variables
    real, allocatable :: delz(:,:,:)
    real, allocatable :: delp(:,:,:)
    real, allocatable :: psfc(:,:)
    real, allocatable :: slmsk(:,:) !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real, allocatable :: f10m(:,:)  !< fm at 10m - Ratio of sigma level 1 wind and 10m wind
    real, allocatable :: u1D(:,:,:)
    real, allocatable :: v1D(:,:,:)
    real, allocatable :: u1(:,:)
    real, allocatable :: v1(:,:)
    real, allocatable :: ak(:)

    real :: ptop
    integer :: z
    integer :: mp_physics
    integer :: ix, iy

    
    call getxyzmax_fv3restart(filename_f3r_phy, xmax, ymax, zmax)
    if (should_print(1)) then
      write(*,*) 'xmax,ymax,zmax: ', xmax, ymax, zmax
    end if
    ! calculate x,y ranges according to namelist
    nml_i_x_calc_beg = max(1,    nml_i_x_calc_beg)
    nml_i_x_calc_end = min(xmax, nml_i_x_calc_end)
    nml_i_y_calc_beg = max(1,    nml_i_y_calc_beg)
    nml_i_y_calc_end = min(ymax, nml_i_y_calc_end)

    xbeg0 = nml_i_x_calc_beg
    xend0 = nml_i_x_calc_end
    ybeg0 = nml_i_y_calc_beg
    yend0 = nml_i_y_calc_end

    if (should_print(1)) then
      write(*,*) 'xbeg0,xend0,ybeg0,yend0: ', xbeg0,xend0,ybeg0,yend0
    end if
    
    call calc_proc_xyrange_simple(xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend )
    if (should_print(10)) then
      write(*,*) 'xbeg,xend,ybeg,yend: ', xbeg,xend,ybeg,yend
    end if
    ! initialize state
    if (present(buffer)) then
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status, buffer)
    else
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status)
    end if
    ! allocate temporary arrays
    xbeg1 = state%xbeg1
    xend1 = state%xend1
    ybeg1 = state%ybeg1
    yend1 = state%yend1
    
    
    allocate( delz(xbeg1:xend1, ybeg1:yend1, zmax), &
              delp(xbeg1:xend1, ybeg1:yend1, zmax), &
              psfc(xbeg1:xend1, ybeg1:yend1), &
              slmsk(xbeg1:xend1, ybeg1:yend1), &
              f10m(xbeg1:xend1, ybeg1:yend1), &
              stat=error_status)
    allocate( u1D(xbeg1:xend1,     ybeg1:(yend1+1), 1 ), & ! D-grid
              v1D(xbeg1:(xend1+1), ybeg1:yend1    , 1 ),    &
              u1 (xbeg1:xend1, ybeg1:yend1),    & ! C-grid
              v1 (xbeg1:xend1, ybeg1:yend1),    &
              stat=error_status)

    allocate( ak(zmax+1), stat=error_status )

    ! check the following file for defination
    ! /work/05012/tg843115/stampede2/FV3/Harvey/FV3_combine/FV3GFS_SRC/fv3_gfsphysics/GFS_layer/GFS_typedefs.F90
    ! /work/05012/tg843115/stampede2/FV3/Harvey/FV3_combine/FV3GFS_SRC/fv3_gfsphysics/FV3GFS/FV3GFS_io.F90
  

    ! 3-d dynamical fields
    ! /work/05012/tg843115/stampede2/FV3/Harvey/FV3_combine/FV3GFS_SRC/atmos_cubed_sphere/tools/fv_io.F90
    ! /work/05012/tg843115/stampede2/FV3/Harvey/FV3_combine/FV3GFS_SRC/atmos_cubed_sphere/model/fv_arrays.F90

    !-----------------------------------------------------------------------
    ! Five prognostic state variables for the f-v dynamics
    !-----------------------------------------------------------------------
    ! dyn_state:
    ! D-grid prognostatic variables: u, v, and delp (and other scalars)
    !
    !     o--------u(i,j+1)----------o
    !     |           |              |
    !     |           |              |
    !  v(i,j)------scalar(i,j)----v(i+1,j)
    !     |           |              |
    !     |           |              |
    !     o--------u(i,j)------------o
    !
    ! The C grid component is "diagnostic" in that it is predicted every time step
    ! from the D grid variables.
    !    real, _ALLOCATABLE :: u(:,:,:)    _NULL  ! D grid zonal wind (m/s)
    !    real, _ALLOCATABLE :: v(:,:,:)    _NULL  ! D grid meridional wind (m/s)
    !    real, _ALLOCATABLE :: pt(:,:,:)   _NULL  ! temperature (K)
    !    real, _ALLOCATABLE :: delp(:,:,:) _NULL  ! pressure thickness (pascal)
    !    real, _ALLOCATABLE :: q(:,:,:,:)  _NULL  ! specific humidity and prognostic constituents
    !    real, _ALLOCATABLE :: qdiag(:,:,:,:)  _NULL  ! diagnostic tracers

    !----------------------
    ! non-hydrostatic state:
    !----------------------------------------------------------------------
    !    real, _ALLOCATABLE ::     w(:,:,:)  _NULL  ! cell center vertical wind (m/s)
    !    real, _ALLOCATABLE ::  delz(:,:,:)  _NULL  ! layer thickness (meters)
    !    real, _ALLOCATABLE ::   ze0(:,:,:)  _NULL  ! height at layer edges for remapping
    !    real, _ALLOCATABLE ::  q_con(:,:,:) _NULL  ! total condensates

    !-----------------------------------------------------------------------
    ! Others:
    !-----------------------------------------------------------------------
    !    real, _ALLOCATABLE :: phis(:,:)     _NULL  ! Surface geopotential (g*Z_surf)
    call open_file( filename_f3r_core, nf_nowrite, fid)

    call get_variable3d_local(fid,'T',    xbeg1, xend1, ybeg1, yend1, 1, zmax, time_idx, state%temperature) ! temperature (K)
    call get_variable3d_local(fid,'DZ',   xbeg1, xend1, ybeg1, yend1, 1, zmax, time_idx, delz) ! actually it is thickness in meters
    call get_variable3d_local(fid,'delp', xbeg1, xend1, ybeg1, yend1, 1, zmax, time_idx, delp) ! actually it is thickness in pressure (pascal)
    call get_variable2d_local(fid,'phis', xbeg1, xend1, ybeg1, yend1,          time_idx, state%hgt)
    
    call get_variable3d_local(fid,'u',     xbeg1, xend1,   ybeg1, yend1+1,    zmax, zmax,  time_idx, u1D)
    call get_variable3d_local(fid,'v',     xbeg1, xend1+1, ybeg1, yend1,      zmax, zmax,  time_idx, v1D)
    u1(xbeg1:xend1, ybeg1:yend1) = ( u1D(xbeg1:xend1, ybeg1:yend1,1) + u1D(xbeg1:xend1, (ybeg1+1):(yend1+1),1) )/2
    v1(xbeg1:xend1, ybeg1:yend1) = ( v1D(xbeg1:xend1, ybeg1:yend1,1) + v1D((xbeg1+1):(xend1+1), ybeg1:yend1,1) )/2

    call close_file(fid)
    
    ! calculate pressure levels using ptop
    call open_file( filename_f3r_ak, nf_nowrite, fid)
    call get_variable1d_local(fid,'ak', 1, zmax+1,  time_idx, ak)
    ptop = ak(1)
    call close_file(fid)
    ! surface properties. Check GFS_sfcprop_type defination
    call open_file( filename_f3r_sfc, nf_nowrite, fid)
    call get_variable2d_local(fid,'slmsk',xbeg1, xend1, ybeg1, yend1,           time_idx, slmsk)
    where (slmsk == 1) 
      state%landmask(xbeg1:xend1,ybeg1:yend1) = 1.
    elsewhere
      state%landmask(xbeg1:xend1,ybeg1:yend1) = 0.
    endwhere
    state%lakemask(xbeg1:xend1,ybeg1:yend1) = 0.
    ! 'tsea' in sfc_ file is 'tsfc'
    call get_variable2d_local(fid,'tsea',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%tsk)
    call get_variable2d_local(fid,'f10m',     xbeg1, xend1, ybeg1, yend1,           time_idx, f10m)
    state%u10(xbeg1:xend1, ybeg1:yend1) = u1(xbeg1:xend1, ybeg1:yend1) * f10m(xbeg1:xend1, ybeg1:yend1) ! gsmphys/sfc_diag.f 
    state%v10(xbeg1:xend1, ybeg1:yend1) = v1(xbeg1:xend1, ybeg1:yend1) * f10m(xbeg1:xend1, ybeg1:yend1)
    call close_file(fid)
    
    ! Qvapor and hydrometers 
    call open_file( filename_f3r_tracer, nf_nowrite, fid)
    call get_variable3d_local(fid,'sphum',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qvapor)
    if (nml_l_include_cloud) then
      call get_variable3d_local(fid,'liq_wat',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qcloud)
      call get_variable3d_local(fid,'rainwat',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qrain)
      call get_variable3d_local(fid,'ice_wat',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qice)
      call get_variable3d_local(fid,'snowwat',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qsnow)
      call get_variable3d_local(fid,'graupel',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qgraup)
      ! Not implemented yet
      !call get_variable3d_local(fid,'o3mr',     xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%ozone)
    end if
    call close_file(fid)

    ! latitude and longitude
    call open_file( filename_f3r_grid, nf_nowrite, fid)
    call get_variable2d_local(fid,'grid_latt',   xbeg1, xend1, ybeg1, yend1,           time_idx, state%lat)
    call get_variable2d_local(fid,'grid_lont',   xbeg1, xend1, ybeg1, yend1,           time_idx, state%lon)
    call close_file(fid)
  

    ! revert z-coordinate
    state%temperature(:,:,1:zmax) = state%temperature(:,:,zmax:1:-1)
    state%qvapor     (:,:,1:zmax) = state%qvapor     (:,:,zmax:1:-1)
    state%qcloud     (:,:,1:zmax) = state%qcloud     (:,:,zmax:1:-1)
    state%qrain      (:,:,1:zmax) = state%qrain      (:,:,zmax:1:-1)
    state%qice       (:,:,1:zmax) = state%qice       (:,:,zmax:1:-1)
    state%qsnow      (:,:,1:zmax) = state%qsnow      (:,:,zmax:1:-1)
    state%qgraup     (:,:,1:zmax) = state%qgraup     (:,:,zmax:1:-1)
    delz     (:,:,1:zmax) = delz    (:,:,zmax:1:-1)
    delp     (:,:,1:zmax) = delp    (:,:,zmax:1:-1)

!    error_status = nf_get_att_int(fid, nf_global, 'MP_PHYSICS', mp_physics)


    !where(state%lon > 180.) state%lon = state%lon - 360.
    ! Remove negative hydro
    where(state%qvapor .lt. 1.0e-10) state%qvapor =1.0e-10
    if (nml_l_include_cloud) then
      where(state%qcloud .lt. 0.0) state%qcloud =0.0
      where(state%qice   .lt. 0.0) state%qice   =0.0
      where(state%qrain  .lt. 0.0) state%qrain  =0.0
      where(state%qsnow  .lt. 0.0) state%qsnow  =0.0
      where(state%qgraup .lt. 0.0) state%qgraup =0.0
    end if
    ! calculate other variables not directly saved in WRF output
    
    state%delz = - delz
    
    ! top layer
    state%level_Pressure(:,:,zmax) = ptop
    do z = zmax, 1, -1
      state%level_pressure(:,:,z-1) = state%level_pressure(:,:,z) + delp(:,:,z)
      state%pressure(:,:,z) = ( state%level_pressure(:,:,z-1) + state%level_pressure(:,:,z) ) /2
    end do ! z
   
!  if (mp_physics == 6) then
!    state%mp_physics_name = 'WSM6'
!  else if (mp_physics == 7) then
!    state%mp_physics_name = 'Goddard'
!  else if (mp_physics == 8) then
    state%mp_physics_name = 'Thompson08'
!  end if
  
  print *, 'qsnow', state%qsnow(xbeg1,ybeg1,:)


  ! set projection
!  call set_domain_proj(filename_core, state%proj)
!  state%dx = state%proj%dx / 1000.

  do ix = state%xbeg1, state%xend1
    do iy = state%ybeg1, state%yend1
      state%x(ix,iy) = ix
      state%y(ix,iy) = iy
    end do
  end do

  end subroutine model_state_read_fv3restart
  !---------------


  !------------------------------------------------------
  ! Read HRRR model output into memory (can also read RAP output) 
  !------------------------------------------------------ 
  subroutine model_state_read_hrrr(state, filename_hrrr, time_idx, modeltype, error_status, buffer)
    type(model_state), intent(inout) :: state
    character(len=*), intent(in) :: filename_hrrr
    integer, intent(in) :: time_idx
    character(len=*), intent(in) :: modeltype 
    integer, intent(out) :: error_status
    integer, intent(in), optional :: buffer
    
    integer :: xmax, ymax, zmax
    integer :: xbeg0, xend0, ybeg0, yend0
    integer :: xbeg1, xend1, ybeg1, yend1
    integer :: xbeg, xend, ybeg, yend
    integer :: fid

    REAL, PARAMETER :: P1000MB=100000.D0
    REAL, PARAMETER :: R_D=287.D0
    REAL, PARAMETER :: CP=7.D0*R_D/2.D0
    ! temporary variables
    real, allocatable :: ph  (:,:,:) ! at full model layer
    real, allocatable :: t   (:,:,:)
    real, allocatable :: psfc(:,:)

    integer :: z
    integer :: mp_physics
    integer :: ix, iy, iz
    character(len=4) :: numstr
    character(len=40) :: fcst_hour_str

    integer :: iret
    real, allocatable :: grid(:,:), lat(:,:), lon(:,:)
    character (len=300) :: metadata, grid_info
    real :: lat1, lon1
    
!    real, parameter :: g=9.80665
    
    ! set projection
    !call set_domain_proj(filename_hrrr, state%proj)
    if (modeltype == 'hrrr') then
      call map_set(PROJ_LC, state%proj, lat1=lat1, lon1=lon1, dx=3000., stdlon=-97.5, &
                  truelat1=38.5, truelat2=38.5)
      state%dx = 3
    else if (modeltype == 'rap') then
!      call map_set(PROJ_ROTLL, state%proj, lat1=47.49999, lon1=-104., &
!                  phi=34.478738791416568, lambda=46.113783422176112, ixdim=758, jydim=567,&
!                stagger=HH)
      call map_set(PROJ_ROTLL, state%proj, lat1=54., lon1=-106., &
                  phi=0.121813*(834.-1)/2, lambda=0.121813*(953.-1)/2, ixdim=953, jydim=834,&
                stagger=HH)
      state%dx = 13
    else
      print *, 'modeltype must be hrrr or rap, current is ', modeltype
      stop
    end if
    
!    call ijll_rotlatlon(1.,1.,state%proj,lat1,lon1)
!    print *, lat1, lon1
!    stop



    fcst_hour_str = ':6 hour fcst:'

    if (my_proc_id == 0) then
      iret = grb2_mk_inv(filename_hrrr,trim(filename_hrrr) // '.inv')
      if (iret .ne. 0) stop 1
    end if
    call MPI_Barrier(comm,ierr)
    
    if (modeltype == 'hrrr') then
      iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':PRES:',':1 hybrid level:',  data2=grid,nx=xmax, ny=ymax,lat=lat, lon=lon, desc=metadata, grid_desc=grid_info)
      where (lon > 180.)
        lon = lon -360
      end where
    else if (modeltype == 'rap') then
      iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':PRES:',':1 hybrid level:',  data2=grid,nx=xmax, ny=ymax, desc=metadata, grid_desc=grid_info)
    end if
    if (iret .ne. 1) stop 1
    if (my_proc_id == 0) then
      print *, metadata
      print *, grid_info
      print *, shape(grid)
      print *, 'xmax, ymax', xmax, ymax
    end if
    zmax=50
    !call getxyzmax_wrfout(filename_hrrr, xmax, ymax, zmax)

    if (should_print(1)) then
      write(*,*) 'xmax,ymax,zmax: ', xmax, ymax, zmax
    end if
    ! calculate x,y ranges according to namelist
    nml_i_x_calc_beg = max(1,    nml_i_x_calc_beg)
    nml_i_x_calc_end = min(xmax, nml_i_x_calc_end)
    nml_i_y_calc_beg = max(1,    nml_i_y_calc_beg)
    nml_i_y_calc_end = min(ymax, nml_i_y_calc_end)

    xbeg0 = nml_i_x_calc_beg
    xend0 = nml_i_x_calc_end
    ybeg0 = nml_i_y_calc_beg
    yend0 = nml_i_y_calc_end

    if (should_print(1)) then
      write(*,*) 'xbeg0,xend0,ybeg0,yend0: ', xbeg0,xend0,ybeg0,yend0
    end if
    
    call calc_proc_xyrange_simple(xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend )
    if (should_print(10)) then
      write(*,*) 'xbeg,xend,ybeg,yend: ', xbeg,xend,ybeg,yend
    end if
    ! initialize state
    if (present(buffer)) then
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status, buffer)
    else
      call model_state_alloc( state, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, zmax, error_status)
    end if
    ! allocate temporary arrays
    xbeg1 = state%xbeg1
    xend1 = state%xend1
    ybeg1 = state%ybeg1
    yend1 = state%yend1
    
    
    allocate( ph  (xbeg1:xend1, ybeg1:yend1, zmax), &
              t   (xbeg1:xend1, ybeg1:yend1, zmax), &
              psfc(xbeg1:xend1, ybeg1:yend1), &
              stat=error_status)


    ! 2-D field
    if (modeltype == 'hrrr') then
      state%lat(xbeg1:xend1, ybeg1:yend1) =  lat(xbeg1:xend1, ybeg1:yend1) 
      state%lon(xbeg1:xend1, ybeg1:yend1) =  lon(xbeg1:xend1, ybeg1:yend1) 
      lat1 = lat(1,1)
      lon1 = lon(1,1)
    else if (modeltype == 'rap') then
      
      do ix = state%xbeg1, state%xend1
        do iy = state%ybeg1, state%yend1
          call ijll_rotlatlon( real(ix), real(iy), state%proj, state%lat(ix, iy), state%lon(ix, iy) )
          state%x(ix,iy) = ix
          state%y(ix,iy) = iy
        end do
      end do
      call ijll_rotlatlon(1.,1.,state%proj, lat1, lon1)
    end if

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':HGT:',':surface:',data2=grid)
    if (iret .ne. 1) stop 1
    state%hgt(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':TMP:',':surface:',data2=grid)
    if (iret .ne. 1) stop 2
    state%tsk(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':PRES:',':surface:',data2=grid)
    if (iret .ne. 1) stop 3
    psfc(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':LAND:',':surface:',data2=grid)
    if (iret .ne. 1) stop 4
    state%landmask(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)
    state%lakemask(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)
    
    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':ICEC:',':surface:',data2=grid)
    if (iret .ne. 1) stop 4
    state%icecover(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':UGRD:',':10 m above ground:',data2=grid)
    if (iret .ne. 1) stop 4
    state%u10(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)

    iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':VGRD:',':10 m above ground:',data2=grid)
    if (iret .ne. 1) stop 4
    state%v10(xbeg1:xend1, ybeg1:yend1) = grid(xbeg1:xend1, ybeg1:yend1)
    ! read variables
    ! 2-D fields
    ! 3-D fields


    do iz = 1, zmax
      write(numstr,'(i4)') iz
      ! pressure
      iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':PRES:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
      if (iret .ne. 1) stop 11
      state%pressure(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)
      ! altitude
      iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':HGT:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
      if (iret .ne. 1) stop 12
      ph(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)
      ! Temperature
      iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':TMP:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
      if (iret .ne. 1) stop 13
      state%temperature(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

      if (nml_l_include_cloud) then
        ! Q
        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':SPFH:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 14
        state%qvapor(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':CLMR:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 15
        state%qcloud(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':CIMIXR:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 16
        state%qice(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':RWMR:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 17
        state%qrain(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':SNMR:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 18
        state%qsnow(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)

        iret = grb2_inq(filename_hrrr, trim(filename_hrrr) // '.inv', ':GRLE:',':'// trim(adjustl(numstr)) // ' hybrid level:',data2=grid)
        if (iret .ne. 1) stop 19
        state%qgraup(xbeg1:xend1, ybeg1:yend1, iz) = grid(xbeg1:xend1, ybeg1:yend1)
      end if
!
    end do

!    if (my_proc_id == 0) print *, state%hgt(510,510)
!    if (my_proc_id == 0) print *, ph(510,510,:)
   
    do iz = 1, zmax
      
    end do
!    call get_variable3d_local(fid,'P',       xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, p)
!    call get_variable3d_local(fid,'PB',      xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, pb)
!    call get_variable3d_local(fid,'PH',      xbeg1, xend1, ybeg1, yend1, 1, zmax+1,time_idx, ph)
!    call get_variable3d_local(fid,'PHB',     xbeg1, xend1, ybeg1, yend1, 1, zmax+1,time_idx, phb)
!    call get_variable3d_local(fid,'T',       xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, t)
!    call get_variable3d_local(fid,'QVAPOR',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qvapor)
!    call get_variable2d_local(fid,'PSFC',    xbeg1, xend1, ybeg1, yend1,           time_idx, psfc)
!    call get_variable2d_local(fid,'TSK',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%tsk)
!    call get_variable2d_local(fid,'HGT',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%hgt)
!    call get_variable2d_local(fid,'LANDMASK',xbeg1, xend1, ybeg1, yend1,           time_idx, state%landmask)
!    call get_variable2d_local(fid,'LAKEMASK',xbeg1, xend1, ybeg1, yend1,           time_idx, state%lakemask)
!    call get_variable2d_local(fid,'U10',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%u10)
!    call get_variable2d_local(fid,'V10',     xbeg1, xend1, ybeg1, yend1,           time_idx, state%v10)
!    
!    if (nml_l_include_cloud) then
!      call get_variable3d_local(fid,'QCLOUD',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qcloud)
!      call get_variable3d_local(fid,'QRAIN',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qrain)
!      call get_variable3d_local(fid,'QICE',    xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qice)
!      call get_variable3d_local(fid,'QSNOW',   xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qsnow)
!      call get_variable3d_local(fid,'QGRAUP',  xbeg1, xend1, ybeg1, yend1, 1, zmax,  time_idx, state%qgraup)
!    end if
!
!    error_status = nf_get_att_int(fid, nf_global, 'MP_PHYSICS', mp_physics)
!    call close_file(fid)
    

    !where(state%lon > 180.) state%lon = state%lon - 360.
    ! Remove negative hydro
    where(state%qvapor .lt. 1.0e-10) state%qvapor =1.0e-10
    if (nml_l_include_cloud) then
      where(state%qcloud .lt. 0.0) state%qcloud =0.0
      where(state%qice   .lt. 0.0) state%qice   =0.0
      where(state%qrain  .lt. 0.0) state%qrain  =0.0
      where(state%qsnow  .lt. 0.0) state%qsnow  =0.0
      where(state%qgraup .lt. 0.0) state%qgraup =0.0
    end if
    ! calculate other variables not directly saved in WRF output
    
!    state%pressure = P + PB
    
!    state%temperature = (T + 300.0) * ( (state%pressure / P1000MB) ** (R_D/CP) )
    
    ! top layer
    state%level_Pressure(:,:,zmax) = state%pressure(:,:,zmax)*1.5 - state%pressure(:,:,zmax-1)*0.5 
    ! lowest layer
    state%level_pressure(:,:,0) =  max(psfc, state%pressure(:,:,1)*1.5 - state%pressure(:,:,2)*0.5)
    do z = 2, zmax
      state%level_pressure(:,:,z-1) = (state%pressure(:,:,z) + state%pressure(:,:,z-1) ) * 0.5
    end do ! z

    state%delz          (:,:,1) = (PH(:,:,1) + PH(:,:,2)) / 2.0 - state%hgt(:,:)
    do z = 2, zmax-1
      state%delz(:,:,z) = (PH(:,:,z+1) - PH(:,:,z-1)) / 2.0
    end do ! z
    state%delz(:,:,zmax) = (PH(:,:,zmax) - PH(:,:,zmax-1) ) / 2.0

   
  if (mp_physics == 6) then
    state%mp_physics_name = 'WSM6'
  else if (mp_physics == 7) then
    state%mp_physics_name = 'Goddard'
  else if (mp_physics == 8) then
    state%mp_physics_name = 'Thompson08'
  else if (mp_physics == 28) then
    state%mp_physics_name = 'Thompson08'
  end if


  do ix = state%xbeg1, state%xend1
    do iy = state%ybeg1, state%yend1
      state%x(ix,iy) = ix
      state%y(ix,iy) = iy
    end do
  end do

  end subroutine model_state_read_hrrr
  !---------------


  !------------------------------------------------------
  ! 
  !------------------------------------------------------
  subroutine model_state_add_weight(state1, x1, y1, z1, state2, x2, y2, z2, weight)
    type(model_state), intent(in) :: state1
    type(model_state), intent(inout) :: state2
    integer, intent(in) :: x1, y1, z1, x2, y2, z2
    real, intent(in) :: weight
    
!    state2%lat(x2,y2) = state2%lat(x2,y2) + state1%lat(x1,y1) * weight
!    state2%lon(x2,y2) = state2%lon(x2,y2) + state1%lon(x1,y1) * weight
!    state2%x  (x2,y2) = state2%x  (x2,y2) + state1%x  (x1,y1) * weight
!    state2%y  (x2,y2) = state2%y  (x2,y2) + state1%y  (x1,y1) * weight
!    state2%tsk(x2,y2) = state2%tsk(x2,y2) + state1%tsk(x1,y1) * weight
!    state2%u10(x2,y2) = state2%u10(x2,y2) + state1%u10(x1,y1) * weight
!    state2%v10(x2,y2) = state2%v10(x2,y2) + state1%v10(x1,y1) * weight
!    state2%landmask(x2,y2)   = state2%landmask(x2,y2)   + state1%landmask(x1,y1) * weight
!    state2%lakemask(x2,y2)   = state2%lakemask(x2,y2)   + state1%lakemask(x1,y1) * weight
!    state2%soiltype(x2,y2)   = state2%soiltype(x2,y2)   + state1%soiltype(x1,y1) * weight
!    state2%hgt     (x2,y2)   = state2%hgt     (x2,y2)   + state1%hgt     (x1,y1) * weight
    state2%delz    (x2,y2,z2) = state2%delz    (x2,y2,z2) + state1%delz    (x1,y1,z1) * weight
    state2%level_pressure(x2,y2,z2) = state2%level_pressure(x2,y2,z2) + state1%level_pressure(x1,y1,z1) * weight
    state2%pressure      (x2,y2,z2) = state2%pressure      (x2,y2,z2) + state1%pressure      (x1,y1,z1) * weight
    state2%temperature   (x2,y2,z2) = state2%temperature   (x2,y2,z2) + state1%temperature   (x1,y1,z1) * weight
    state2%qvapor        (x2,y2,z2) = state2%qvapor        (x2,y2,z2) + state1%qvapor        (x1,y1,z1) * weight

    if (nml_l_include_cloud) then
      state2%qcloud(x2,y2,z2) =  state2%qcloud(x2,y2,z2) + state1%qcloud(x1,y1,z1) * weight
      state2%qrain (x2,y2,z2) =  state2%qrain (x2,y2,z2) + state1%qrain (x1,y1,z1) * weight
      state2%qice  (x2,y2,z2) =  state2%qice  (x2,y2,z2) + state1%qice  (x1,y1,z1) * weight
      state2%qsnow (x2,y2,z2) =  state2%qsnow (x2,y2,z2) + state1%qsnow (x1,y1,z1) * weight
      state2%qgraup(x2,y2,z2) =  state2%qgraup(x2,y2,z2) + state1%qgraup(x1,y1,z1) * weight
      state2%nrain (x2,y2,z2) =  state2%nrain (x2,y2,z2) + state1%nrain (x1,y1,z1) * weight
    end if

  end subroutine model_state_add_weight
  !---------------

end module rt_state_module
