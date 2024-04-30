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


  ! CAUTION: running MW simulations with slant-path rad transfer need to 
  ! set use_microwave_obs to .true. and provide ascii MW observation file
  ! with nml_s_filename_obs. Otherwise, slant-path rad transfer will be 
  ! turned off.
  logical, parameter :: use_microwave_obs = .true.

  call obs_info_init(obs_info)
  
  call parallel_start()
  call get_command_argument(1, filename_nml, l, error_status)

  call namelist_read(filename_nml)
  call namelist_handle_args( error_status )
  call namelist_handle_dependency()
  
  if (should_print(1)) call print_namelist()
  filename_in = nml_s_filename_input
  filename_out = nml_s_filename_output
  call model_state_read_wrf(state, filename_in, 1, error_status, buffer=15)
  !call calc_zenith_GEO( state, nml_s_sensor_id, error_status)
  
  if (use_microwave_obs) then
  filename_obs = nml_s_filename_obs
    call get_microwave(filename_obs, state, obs_mw, error_status)
    call get_obs_info_fromfile(obs_info, filename_obs, error_status)
  end if
  
  do isensor = 1, MAX_SENSORS
    do ich = 1, MAX_CHANNELS
      call obs_info_add(obs_info, nml_s_sensor_id(isensor), nml_a_channels(ich, isensor))
    end do
  end do
  
  
  if (my_proc_id==0) call obs_info_print(obs_info)
  call rt_crtm_init( nml_s_sensor_id, state%nz, error_status, obs_info%channels)
  if (should_print(1)) call print_namelist()
  n_channels = rt_crtm_n_channels(1)
  print *, 'nch1 ', n_channels
  call rt_result_alloc_from_state( res(1), state, n_channels, error_status)
  n_channels = rt_crtm_n_channels(2)
  print *, 'nch2 ', n_channels
  call rt_result_alloc_from_state( res(2), state, n_channels, error_status)
  
  if (use_microwave_obs) then
    call calc_zenith_obs_microwave( state, res(1), obs_mw, 'ssmis_f17', error_status, 'boxsearch')
    call calc_zenith_obs_microwave( state, res(2), obs_mw, 'ssmis_f17', error_status, 'boxsearch')
    
    call calc_slant_path(state, res(1), slant1, error_status)
    call calc_slant_path(state, res(2), slant2, error_status)
   
    call rt_crtm_main(1, slant1, res(1), error_status )
    call rt_crtm_main(2, slant2, res(2), error_status )
  else
    res(1)%zenith_angle = 49.1
    res(1)%azimuth_angle = 0.0
    res(2)%zenith_angle = 52.8
    res(2)%azimuth_angle = 0.0
    call rt_crtm_main(1, state, res(1), error_status )
    call rt_crtm_main(2, state, res(2), error_status )
  end if
  
  call rt_crtm_destroy( error_status )
  
  if (should_print(1)) call print_namelist()
    
  call rt_result_output_nc(res(1), state, filename_out, error_status)
  call rt_result_output_nc(res(2), state, filename_out, error_status)

  if (use_microwave_obs) then
    allocate( tb_conv(obs_mw%num) )
    call rt_result_convolution( res(1), state,obs_mw%num, obs_mw%ch, obs_mw%ii, obs_mw%jj, obs_mw%efov_aScan, obs_mw%efov_cscan, obs_mw%azimuth_angle, tb_conv, error_status)
    call rt_result_convolution( res(2), state,obs_mw%num, obs_mw%ch, obs_mw%ii, obs_mw%jj, obs_mw%efov_aScan, obs_mw%efov_cscan, obs_mw%azimuth_angle, tb_conv, error_status)
    call rt_conv_output_txt( obs_info, tb_conv, filename_out, error_status, obs_mw)
    deallocate( tb_conv )
  end if
  
  
  call rt_result_dealloc(res(1), error_status)
  call rt_result_dealloc(res(2), error_status)
  

  call parallel_finish()
end program test


