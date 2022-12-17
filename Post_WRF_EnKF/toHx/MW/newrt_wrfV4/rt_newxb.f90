module rt_newxb
  
  use rt_namelist_module
  use rt_state_module
  use rt_result_module
  use rt_crtm_module
  use rt_rttov_module
  use rt_zenith_module
  use rt_util_module
  use rt_obs_module
  implicit none

contains

  subroutine newrt_init(filename_nml, nz,  iob_radmin, iob_radmax )
    character(len=*), intent(in) :: filename_nml
    integer, intent(in) :: nz
    integer, intent(out) :: iob_radmin
    integer, intent(out) :: iob_radmax

    type(obs_info_type) :: obs_info
    integer :: iob, isensor, ich
    character(len=20) :: obstype

    integer :: error_status

    call obs_info_init(obs_info)
    call namelist_read(filename_nml)
    call namelist_handle_dependency()
    
    if (my_proc_id == 0) print *, 'newrt init'
    
    if (should_print(1)) call print_namelist()
    iob_radmin = -1
    do iob = 1, obs%num
      obstype = obs%type(iob)
      if (obstype(1:9) == 'Microwave' ) then
        if (iob_radmin < 0) iob_radmin = iob
        iob_radmax = iob
        call obs_info_add( obs_info, obs%sat(iob), obs%ch(iob) )
      end if
    end do

!    do isensor = 1, MAX_SENSORS
!      do ich = 1, MAX_CHANNELS
!        call obs_info_add(obs_info, nml_s_sensor_id(isensor), nml_a_channels(ich, isensor))
!      end do
!    end do
    
    if (my_proc_id==0) call obs_info_print(obs_info)

    if (nml_s_rt_program == 'crtm') then
      call rt_crtm_init( nml_s_sensor_id, nz, error_status, obs_info%channels)
      if (should_print(1)) call print_namelist()
    end if
    if (my_proc_id==0) print *, 'newrt_init finished'

  end subroutine newrt_init

  subroutine xb_to_newrt(filename_wrf, tb_conv)
    character(len=*), intent(in) :: filename_wrf
    real, intent(inout) :: tb_conv(:)

    character(len=256) :: filename_out
    type(model_state) :: state, slant1, slant2
    type(rt_result) :: res
    integer :: error_status
    integer :: n_channels, n_sensors, ich, isensor

    logical :: use_crtm 
    integer :: l
    type(microwave_data_type) :: obs_mw 

    if (my_proc_id == 0) print *, 'newrt used'

    n_sensors = rt_crtm_n_sensors()
    
    filename_out = './tb_out/' // trim(adjustl(filename_wrf))
    call model_state_read_wrf(state, filename_wrf, 1, error_status, buffer=15)
    
    if (nml_s_rt_program == 'crtm') then
      do isensor = 1, n_sensors
          n_channels = rt_crtm_n_channels(isensor)
          print *, 'nch1 ', n_channels
          call rt_result_alloc_from_state( res, state, n_channels, error_status)
          
          call calc_zenith_obs( state, res, obs, rt_crtm_sensor_id(isensor), error_status, 'boxsearch')

          !call calc_slant_path(state, res, slant1, error_status)
         
          call rt_crtm_main(isensor, state, res, error_status )
          
          
          if (should_print(1)) call print_namelist()
            
          call rt_result_output_nc(res, state, filename_out, error_status)


          call rt_result_convolution( res, state, obs%num, obs%ch, obs%position(:,1), obs%position(:,2), obs%efov_aScan, obs%efov_cscan, obs%azimuth_angle,tb_conv, error_status)

          !call rt_conv_output_txt( obs_info, obs_mw, tb_conv, filename_out, error_status)
          call rt_result_dealloc(res, error_status)
        end do
      end if
      
  end subroutine xb_to_newrt

  subroutine newrt_destroy()
    
    integer :: error_status
    call rt_crtm_destroy( error_status )

  end subroutine newrt_destroy

end module rt_newxb
