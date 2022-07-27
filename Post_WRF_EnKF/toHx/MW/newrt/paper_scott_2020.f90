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
  type(model_state) :: state
  type(model_state) :: slant(2)
  type(rt_result) :: res(2)
  integer :: error_status
  integer :: n_channels, ich, isensor

  logical :: use_crtm 
  integer :: l
  type(microwave_data_type) :: obs_mw 
  
  real, allocatable :: tb_conv(:)
  
  


  type(obs_info_type) :: obs_info
  
  call parallel_start()

  call obs_info_init(obs_info)
  
  call get_command_argument(1, filename_nml, l, error_status)

  call namelist_read(filename_nml)
  call namelist_handle_args( error_status )
  call namelist_handle_dependency()
  
  if (should_print(1)) call print_namelist()
  filename_in = nml_s_filename_input
  filename_out = nml_s_filename_output
  
  
  
  call model_state_read_wrf(state, filename_in, 1, error_status, buffer=20)
  filename_obs = nml_s_filename_obs
  
  call get_microwave(filename_obs, state, obs_mw, error_status)
  
!  call get_obs_info(obs_info, obs_mw)
  call get_obs_info_fromfile(obs_info, filename_obs, error_status)
  do isensor = 1, MAX_SENSORS
    do ich = 1, MAX_CHANNELS
      call obs_info_add(obs_info, nml_s_sensor_id(isensor), nml_a_channels(ich, isensor))
    end do
  end do
  
  
  if (my_proc_id==0) call obs_info_print(obs_info)

  
  call rt_crtm_init( obs_info%sensors(1:obs_info%nsensors), state%nz, error_status, obs_info%channels)
  if (should_print(1)) call print_namelist()
  
  
  
  do isensor = 1, obs_info%nsensors
    n_channels = rt_crtm_n_channels(isensor)
    if (my_proc_id == 0) print *, 'nch', isensor, n_channels
    call rt_result_alloc_from_state( res(isensor), state, n_channels, error_status)
    
    call calc_zenith_obs_microwave( state, res(isensor), obs_mw, obs_info%sensors(isensor), error_status, 'boxsearch')
    call calc_slant_path(state, res(isensor), slant(isensor), error_status)
    
    !slant(isensor)%qcloud = 0.0
    !slant(isensor)%qice   = 0.0
    !slant(isensor)%qrain  = 0.0
    !slant(isensor)%qsnow  = 0.0
    !slant(isensor)%qgraup = 0.0
    
    call rt_crtm_main(isensor, slant(isensor), res(isensor), error_status )
  end do
  !call calc_zenith_obs_microwave( state, res(1), obs_mw, 'gmi_gpm_hf', error_status, 'boxsearch')
  !call calc_zenith_obs_microwave( state, res(2), obs_mw, 'gmi_gpm_lf', error_status, 'boxsearch')
  !call calc_slant_path(state, res(1), slant(1), error_status)
  !call calc_slant_path(state, res(2), slant(2), error_status)
  !call rt_crtm_main(1, slant(1), res(1), error_status )
  !call rt_crtm_main(2, slant(2), res(2), error_status )
  
!    call rt_crtm_main(1, state, res(1), error_status )
!    call rt_crtm_main(2, state, res(2), error_status )
  
  call rt_crtm_destroy( error_status )
  

  allocate( tb_conv(obs_mw%num) )
  
  do isensor = 1, obs_info%nsensors
    call rt_result_convolution( res(isensor), state,obs_mw%num, obs_mw%ch, obs_mw%ii, obs_mw%jj, obs_mw%efov_aScan, obs_mw%efov_cscan, obs_mw%azimuth_angle, tb_conv, error_status)
  end do

  call rt_conv_output_txt( obs_info, tb_conv, filename_out, error_status, obs_mw)
  deallocate( tb_conv )
  
  do isensor = 1, obs_info%nsensors
    call rt_result_output_nc(res(isensor), state, filename_out, error_status)
    call rt_result_dealloc(res(isensor), error_status)
  end do
  
  if (my_proc_id == 0) print *, "SUCCESS"
  call parallel_finish()

end program test


