module rt_namelist_module
  use rt_constant_module 
  use netcdf
  use mpi_module
  implicit none
  save

  !-------------------------------------
  ! settings
  !-------------------------------------
  character(len=10)  :: nml_s_rt_program = 'crtm'
  integer :: nml_i_debug_level = 1
  logical :: nml_l_include_land = .true.
  logical :: nml_l_include_icecover = .true.
  logical :: nml_l_include_cloud = .true.
  character(len=LENGTH_SENSOR_NAME) :: nml_s_sensor_id(MAX_SENSORS)
  integer :: nml_i_nchannels(MAX_SENSORS)
  integer :: nml_a_channels(MAX_CHANNELS, MAX_SENSORS)
  character(len=20) :: nml_s_reff_method = 'default'
  integer :: nml_i_x_calc_beg  = 1
  integer :: nml_i_x_calc_end  = 999999
  integer :: nml_i_x_calc_step = 1
  integer :: nml_i_y_calc_beg  = 1
  integer :: nml_i_y_calc_end  = 999999
  integer :: nml_i_y_calc_step = 1
  integer :: nml_i_nicpu = -1
  integer :: nml_i_njcpu = -1
  ! used for RTTOV
  character(len=256) :: nml_s_rttov_mietable = ''
  ! used for CRTM
  character(len=256) :: nml_s_crtm_rainLUT = ''
  character(len=256) :: nml_s_crtm_snowLUT = ''
  character(len=256) :: nml_s_crtm_graupelLUT = ''
  integer :: nml_i_nstream = -1

  namelist /rt_settings/ nml_s_rt_program, nml_i_debug_level, &
                         nml_l_include_cloud, nml_l_include_land, nml_l_include_icecover, &
                         nml_s_sensor_id, nml_i_nchannels, nml_a_channels, &
                         nml_s_reff_method, &
                         nml_i_x_calc_beg, nml_i_x_calc_end, nml_i_x_calc_step, &
                         nml_i_y_calc_beg, nml_i_y_calc_end, nml_i_y_calc_step, &
                         nml_i_nicpu,      nml_i_njcpu, &
                         nml_s_rttov_mietable, &
                         nml_s_crtm_rainLUT, nml_s_crtm_snowLUT, nml_s_crtm_graupelLUT, nml_i_nstream



  !-------------------------------------
  ! input
  !-------------------------------------
  character(len=256) :: nml_s_filename_input  = ''
  character(len=256) :: nml_s_filename_obs    = ''
  ! used by FV3_restart
  character(len=256) :: nml_s_filename_f3r_core   = ''
  character(len=256) :: nml_s_filename_f3r_ak     = ''
  character(len=256) :: nml_s_filename_f3r_tracer = ''
  character(len=256) :: nml_s_filename_f3r_phy    = ''
  character(len=256) :: nml_s_filename_f3r_sfc    = ''
  character(len=256) :: nml_s_filename_f3r_grid   = ''


  namelist /rt_input/ nml_s_filename_input,      nml_s_filename_obs,     &
                      nml_s_filename_f3r_core,   nml_s_filename_f3r_ak,  &
                      nml_s_filename_f3r_tracer, nml_s_filename_f3r_phy, &
                      nml_s_filename_f3r_sfc,    nml_s_filename_f3r_grid


  !-------------------------------------
  ! output
  !-------------------------------------
  character(len=256) :: nml_s_filename_output = ''
  logical :: nml_l_output_reff = .false.

  namelist /rt_output/ nml_s_filename_output, &
                       nml_l_output_reff

  ! non-namelist variables but should be saved
  integer :: global_my_iid = -1
  integer :: global_my_jid = -1
contains

  subroutine namelist_read(filename_nml)
    character(len=*), intent(in) :: filename_nml
    
    nml_s_sensor_id(:) = ''
    nml_i_nchannels(:) = 0
    nml_a_channels(:, :) = 0
    open(10, file=filename_nml)
    read(10, nml=rt_settings)
    read(10, nml=rt_input)
    read(10, nml=rt_output)

  end subroutine namelist_read
  ! --------------
  subroutine print_namelist()
    print rt_settings
    print rt_input
    print rt_output
  end subroutine print_namelist
  ! --------------
  subroutine namelist_check( error_status )
    integer, intent(out) :: error_status
    error_status = SUCCESS

  end subroutine namelist_check
  ! --------------
  subroutine namelist_handle_args( error_status )
    integer, intent(out) :: error_status

    integer :: n, i, l
    character(len=256) :: str

    !character(len=10) :: test_str
    !integer :: test_int
    !real :: test_real
    !logical :: test_bool
    n = command_argument_count()
    do i = 1, n
      call get_command_argument(i, str, l, error_status)
      !print *, i, trim(str)
      !call namelist_check_arg(str, 'abc', value_str=test_str )
      !call namelist_check_arg(str, 'int', value_int=test_int )
      !call namelist_check_arg(str, 'real', value_real=test_real )
      !call namelist_check_arg(str, 'bool', value_bool=test_bool )
      !print *, test_str, test_int, test_real, test_bool
      call namelist_check_arg(str, 'nml_s_rt_program',    value_str =nml_s_rt_program   )
      call namelist_check_arg(str, 'nml_i_debug_level',   value_int =nml_i_debug_level  )
      call namelist_check_arg(str, 'nml_l_include_land',  value_bool=nml_l_include_land )
      call namelist_check_arg(str, 'nml_l_include_icecover',  value_bool=nml_l_include_icecover )
      call namelist_check_arg(str, 'nml_l_include_cloud', value_bool=nml_l_include_cloud)
!      call namelist_check_arg(str, 'nml_s_sensor_id',     value_str =nml_s_sensor_id    )
!      call namelist_check_arg(str, 'nml_i_nchannels',     value_int =nml_i_nchannels    )
      call namelist_check_arg(str, 'nml_s_reff_method',   value_str =nml_s_reff_method  )
      call namelist_check_arg(str, 'nml_i_x_calc_beg',    value_int =nml_i_x_calc_beg   )
      call namelist_check_arg(str, 'nml_i_x_calc_end',    value_int =nml_i_x_calc_end   )
      call namelist_check_arg(str, 'nml_i_x_calc_step',   value_int =nml_i_x_calc_step  )
      call namelist_check_arg(str, 'nml_i_y_calc_beg',    value_int =nml_i_y_calc_beg   )
      call namelist_check_arg(str, 'nml_i_y_calc_end',    value_int =nml_i_y_calc_end   )
      call namelist_check_arg(str, 'nml_i_y_calc_step',   value_int =nml_i_y_calc_step  )
      call namelist_check_arg(str, 'nml_i_nicpu',         value_int =nml_i_nicpu)
      call namelist_check_arg(str, 'nml_i_njcpu',         value_int =nml_i_njcpu)
      call namelist_check_arg(str, 'nml_s_rttov_mietable',value_str =nml_s_rttov_mietable)
      call namelist_check_arg(str, 'nml_s_crtm_rainLUT', value_str =nml_s_crtm_rainLUT)
      call namelist_check_arg(str, 'nml_s_crtm_snowLUT', value_str =nml_s_crtm_snowLUT)
      call namelist_check_arg(str, 'nml_s_crtm_graupelLUT', value_str =nml_s_crtm_graupelLUT)
      call namelist_check_arg(str, 'nml_s_filename_input',value_str =nml_s_filename_input)
      call namelist_check_arg(str, 'nml_s_filename_obs',  value_str =nml_s_filename_obs)
      
      call namelist_check_arg(str, 'nml_s_filename_f3r_core',value_str   =nml_s_filename_f3r_core)
      call namelist_check_arg(str, 'nml_s_filename_f3r_ak',value_str     =nml_s_filename_f3r_ak)
      call namelist_check_arg(str, 'nml_s_filename_f3r_tracer',value_str =nml_s_filename_f3r_tracer)
      call namelist_check_arg(str, 'nml_s_filename_f3r_phy',value_str    =nml_s_filename_f3r_phy)
      call namelist_check_arg(str, 'nml_s_filename_f3r_sfc',value_str    =nml_s_filename_f3r_sfc)
      call namelist_check_arg(str, 'nml_s_filename_f3r_grid',value_str   =nml_s_filename_f3r_grid)
      call namelist_check_arg(str, 'nml_i_nstream',            value_int =nml_i_nstream)
      
      call namelist_check_arg(str, 'nml_s_filename_output',value_str =nml_s_filename_output)
      
    end do
  end subroutine namelist_handle_args
  ! --------------
  subroutine namelist_check_arg(arg, myname, value_int, value_real, value_str, value_bool )
    character(len=*), intent(in) :: arg
    character(len=*), intent(in) :: myname
    integer,          optional, intent(inout) :: value_int
    real,             optional, intent(inout) :: value_real
    character(len=*), optional, intent(inout) :: value_str
    logical,          optional, intent(inout) :: value_bool

    integer :: idx
    idx = index(arg, '=')
    if (idx == 0) return ! return if not find '='
!    print *, arg(:idx-1), trim(arg(idx+1:))
    if (arg(:idx-1) == myname) then
      if ( present( value_str)) read(arg(idx+1:), '(A)') value_str
      if ( present( value_int)) read(arg(idx+1:), '(I)') value_int
      if ( present( value_real)) read(arg(idx+1:),'(F)') value_real
      if ( present( value_bool)) read(arg(idx+1:),'(L)') value_bool
    end if
  end subroutine namelist_check_arg

  ! --------------

  subroutine namelist_handle_dependency()
    if (nml_l_include_cloud == .false.) then
      nml_l_output_reff = .false.
    end if

    ! CPU partition
    if (nml_i_nicpu <= 0 .or. nml_i_njcpu <= 0 .or. nml_i_nicpu * nml_i_njcpu > nprocs) then
      nml_i_nicpu = int(sqrt(real(nprocs)))
      nml_i_njcpu = nprocs/nml_i_nicpu
    end if
    if (should_print(1)) print *, 'nicpu, njcpu=', nml_i_nicpu, nml_i_njcpu
    if (nml_i_nicpu * nml_i_njcpu .ne. nprocs) then
      print *, "nicpu*njcpu != nprocs, exit"
      stop
    end if
    global_my_iid = mod( my_proc_id,nml_i_nicpu)
    global_my_jid = int( my_proc_id/nml_i_nicpu)
    
    if (should_print(10))  print *, 'my_proc_id, my_iid, my_jid=', my_proc_id, global_my_iid, global_my_jid

  end subroutine namelist_handle_dependency
  ! --------------
  subroutine namelist_output_nc(fid, error_status)
    integer, intent(in) :: fid
    integer, intent(out) :: error_status
!  logical :: nml_l_output_reff = .false.
    ! settings

    call namelist_output_nc_one(fid, 'nml_s_rt_program',      error_status, value_str  = nml_s_rt_program )
    call namelist_output_nc_one(fid, 'nml_l_include_land',    error_status, value_bool = nml_l_include_land )
    call namelist_output_nc_one(fid, 'nml_l_include_icecover',    error_status, value_bool = nml_l_include_icecover )
    call namelist_output_nc_one(fid, 'nml_l_include_cloud',   error_status, value_bool = nml_l_include_cloud )
!    call namelist_output_nc_one(fid, 'nml_s_sensor_id',       error_status, value_str  = nml_s_sensor_id)
    call namelist_output_nc_one(fid, 'nml_s_reff_method',     error_status, value_str  = nml_s_reff_method)
    call namelist_output_nc_one(fid, 'nml_i_x_calc_beg',      error_status, value_int  = nml_i_x_calc_beg)
    call namelist_output_nc_one(fid, 'nml_i_x_calc_end',      error_status, value_int  = nml_i_x_calc_end)
    call namelist_output_nc_one(fid, 'nml_i_x_calc_step',     error_status, value_int  = nml_i_x_calc_step)
    call namelist_output_nc_one(fid, 'nml_i_y_calc_beg',      error_status, value_int  = nml_i_y_calc_beg)
    call namelist_output_nc_one(fid, 'nml_i_y_calc_end',      error_status, value_int  = nml_i_y_calc_end)
    call namelist_output_nc_one(fid, 'nml_i_y_calc_step',     error_status, value_int  = nml_i_y_calc_step)
    call namelist_output_nc_one(fid, 'nml_s_filename_input',  error_status, value_str  = nml_s_filename_input)
    call namelist_output_nc_one(fid, 'nml_s_filename_obs',    error_status, value_str  = nml_s_filename_obs)
    
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_core',  error_status, value_str   = nml_s_filename_f3r_core)
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_ak',  error_status, value_str     = nml_s_filename_f3r_ak)
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_tracer',  error_status, value_str = nml_s_filename_f3r_tracer)
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_phy',  error_status, value_str    = nml_s_filename_f3r_phy)
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_sfc',  error_status, value_str    = nml_s_filename_f3r_sfc)
    call namelist_output_nc_one(fid, 'nml_s_filename_f3r_grid',  error_status, value_str   = nml_s_filename_f3r_grid)
    
    call namelist_output_nc_one(fid, 'nml_s_filename_output', error_status, value_str  = nml_s_filename_output)
    call namelist_output_nc_one(fid, 'nml_l_output_reff',     error_status, value_bool = nml_l_output_reff )
    call namelist_output_nc_one(fid, 'nml_i_nicpu',           error_status, value_int  = nml_i_nicpu)
    call namelist_output_nc_one(fid, 'nml_i_njcpu',           error_status, value_int  = nml_i_njcpu)
    call namelist_output_nc_one(fid, 'nml_s_rttov_mietable',  error_status, value_str  = nml_s_rttov_mietable)
    call namelist_output_nc_one(fid, 'nml_s_crtm_rainLUT',    error_status, value_str  = nml_s_crtm_rainLUT)
    call namelist_output_nc_one(fid, 'nml_s_crtm_snowLUT',    error_status, value_str  = nml_s_crtm_snowLUT)
    call namelist_output_nc_one(fid, 'nml_s_crtm_graupelLUT',    error_status, value_str  = nml_s_crtm_graupelLUT)
    call namelist_output_nc_one(fid, 'nml_i_nstream',           error_status, value_int  = nml_i_nstream)
  end subroutine namelist_output_nc

  ! --------------
  function nf_put_att_logical(fid, varid, sname, n, var)
    integer, intent(in) :: fid, varid, n
    character(len=*), intent(in) :: sname
    logical, intent(in) :: var

    logical :: nf_put_att_logical
    if (var) then
      nf_put_att_logical = nf_put_att_text(fid, varid, sname, n, 'T')
    else
      nf_put_att_logical = nf_put_att_text(fid, varid, sname, n, 'F')
    end if
  end function nf_put_att_logical

  ! --------------
  subroutine namelist_output_nc_one(fid, myname, error_status, value_int, value_real, value_str, value_bool )
    integer,          intent(in) :: fid
    character(len=*), intent(in) :: myname
    integer,          intent(out) :: error_status
    integer,          optional, intent(inout) :: value_int
    real,             optional, intent(inout) :: value_real
    character(len=*), optional, intent(inout) :: value_str
    logical,          optional, intent(inout) :: value_bool

    if ( present( value_str)) then
      error_status = nf_put_att_text   (fid, nf_global, myname,           len(trim(value_str)), value_str)
    else if ( present( value_int)) then
      error_status = nf_put_att_int    (fid, nf_global, myname, NF_INT,   1,                    value_int)
    else if ( present( value_real)) then
      error_status = nf_put_att_real   (fid, nf_global, myname, NF_FLOAT, 1,                    value_real)
    else if ( present( value_bool)) then
      error_status = nf_put_att_logical(fid, nf_global, myname,           1,                    value_bool)
    end if
  end subroutine namelist_output_nc_one 
  ! --------------
  function should_print(level)
    integer, intent(in) :: level
    logical :: should_print
    
    should_print = nml_i_debug_level >= level .and. (level >=10 .or. my_proc_id == 0)
  end function should_print
end module rt_namelist_module
