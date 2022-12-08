program test

  use rt_namelist_module
  use rt_state_module
  use rt_result_module
  use rt_crtm_module
  use rt_rttov_module
  use rt_zenith_module
  use rt_util_module
  use rt_obs_module
  implicit none

  character(len=256) :: filename_nml!='/home1/05012/tg843115/src/forward_rt/test/test_crtm.nml'
  character(len=256) :: filename_out
  character(len=256) :: filename_in
  character(len=256) :: filename_obs
  type(model_state) :: state, slant1, slant2
  type(rt_result) :: res(2)
  integer :: error_status
  integer :: n_channels, ich, isensor

  logical :: use_crtm 
  integer :: l
  type(microwave_data_type) :: obs_mw 
  
  real, allocatable :: tb_conv(:)
  
  
  ! test Gaussian function
  real :: w(100,100)
  real :: x0, y0
  real :: x(100,100), y(100,100)
  integer :: i, j
  real :: efov_a, efov_c, dx, azi


  type(obs_info_type) :: obs_info

  call obs_info_init(obs_info)
!  call obs_info_print(obs_info)

!  call obs_info_add(obs_info,'ccc',3)
!  call obs_info_print(obs_info)
!  call obs_info_add(obs_info,'aaa',1)
!  call obs_info_print(obs_info)
!  call obs_info_add(obs_info,'ddd',4)
!  call obs_info_print(obs_info)
!  call obs_info_add(obs_info,'bbb',2)
!  call obs_info_print(obs_info)
!  call obs_info_add(obs_info,'aaa',5)
!  call obs_info_print(obs_info)
 
!  stop

  efov_a = 25 
  efov_c = 25
  dx = 1
  azi = 30
  
  do i = 1, 100
    do j = 1, 100
      x(i,j) = i
      y(i,j) = j
    end do
  end do
  x0 = 50
  y0 = 50
  call gaussian_ellipse_weight(x0, y0, 100,100, x, y, efov_a, efov_c, dx, azi, w)
  
  
  call parallel_start()
  call get_command_argument(1, filename_nml, l, error_status)

  call namelist_read(filename_nml)
  call namelist_handle_args( error_status )
  call namelist_handle_dependency()
  
  if (should_print(1)) call print_namelist()
  filename_in = nml_s_filename_input
  filename_out = nml_s_filename_output
  call model_state_read_fv3restart(state,  nml_s_filename_f3r_core, nml_s_filename_f3r_ak, nml_s_filename_f3r_tracer, nml_s_filename_f3r_phy, nml_s_filename_f3r_sfc, nml_s_filename_f3r_grid, 1, error_status)
  filename_obs = nml_s_filename_obs
  
!  call get_microwave(filename_obs, state, obs_mw, error_status)
  
!  call get_obs_info(obs_info, obs_mw)
!  call get_obs_info_fromfile(obs_info, filename_obs, error_status)
  do isensor = 1, MAX_SENSORS
    do ich = 1, MAX_CHANNELS
      call obs_info_add(obs_info, nml_s_sensor_id(isensor), nml_a_channels(ich, isensor))
    end do
  end do
  
  
  if (my_proc_id==0) call obs_info_print(obs_info)
  if (nml_s_rt_program == 'crtm') then
    !call rt_crtm_init( nml_s_sensor_id, state%nz, error_status, nml_a_channels)
    call rt_crtm_init( nml_s_sensor_id, state%nz, error_status, obs_info%channels)
    if (should_print(1)) call print_namelist()
    n_channels = rt_crtm_n_channels(1)
    print *, 'nch1 ', n_channels
    call rt_result_alloc_from_state( res(1), state, n_channels, error_status)

    call calc_zenith_GEO(state, res(1), nml_s_sensor_id(1), error_status) 
    
!    call calc_zenith_obs_microwave( state, res(1), obs_mw, 'gmi_gpm_hf', error_status, 'boxsearch')
!    call rt_crtm_main(1, slant1, res(1), error_status )
    
    call rt_crtm_main(1, state, res(1), error_status )
    
    call rt_crtm_destroy( error_status )
    call rt_result_output_nc(res(1), state, filename_out, error_status)
  else if (nml_s_rt_program == 'rttov') then
    
    do isensor=1, 2    
      if (nml_s_sensor_id(isensor) == 'abi_gr') then
        call rt_rttov_init(state, nml_s_sensor_id(1),  (/7,8,9,10,11,12,13,14,15,16/) )
      else if (nml_s_sensor_id(isensor) == 'gmi_gpm_lf') then
        call rt_rttov_init(state, 'gmi_gpm_lf',  (/1,2,3,4,5,6,7,8,9/) )
      else if (nml_s_sensor_id(isensor) == 'gmi_gpm_hf') then
        call rt_rttov_init(state, 'gmi_gpm_hf',  (/10,11,12,13/) )
      else if (nml_s_sensor_id(isensor) == 'ssmis_f17') then
        call rt_rttov_init(state, nml_s_sensor_id(isensor),  (/9,12,13,15,16,17,18/) )
      end if
      n_channels = rt_rttov_n_channels()
      call rt_result_alloc_from_state( res(isensor), state, n_channels, error_status)
      
      call calc_zenith_obs_microwave( state, res(isensor), obs_mw, nml_s_sensor_id(isensor), error_status, 'boxsearch')
  !    call calc_zenith_obs_microwave( state, res(2), obs_mw, 'gmi_gpm_lf', error_status, 'boxsearch')
     
      call rt_rttov_main( state, res(isensor), error_status )
      call rt_rttov_destroy( error_status )
      call rt_result_output_nc(res(isensor), state, filename_out, error_status)
    end do 
  end if
  
  if (should_print(1)) call print_namelist()
    
!  call rt_result_output_nc(res(2), state, filename_out, error_status)

  call rt_result_dealloc(res(1), error_status)

  
  if (my_proc_id == 0) print *, "SUCCESS"
  call parallel_finish()

end program test


