module rt_result_module

  use rt_namelist_module
  use rt_constant_module
  use rt_state_module
  use rt_util_module
  use obs_define
  use rt_obs_module
  implicit none
  type rt_result
    logical :: l_allocated = .false.

    character(len=LENGTH_SENSOR_NAME) :: sensor_id = ''
    integer :: channels(MAX_CHANNELS)
    ! dimensions
    integer :: xbeg0 = -1
    integer :: xend0 = -1
    integer :: ybeg0 = -1
    integer :: yend0 = -1
    integer :: xbeg = -1
    integer :: xend = -1
    integer :: ybeg = -1
    integer :: yend = -1

    integer :: nch = -1
    integer :: nz = -1

    ! channel-dependent var
    real, allocatable :: TB(:,:,:) ! nx*ny*nch
    real, allocatable :: emis(:,:,:) ! nx*ny*nch

    ! z-dependent var
    real, allocatable :: reff_cloud(:,:,:) ! nx*ny*nz
    real, allocatable :: reff_rain (:,:,:) ! nx*ny*nz
    real, allocatable :: reff_ice  (:,:,:) ! nx*ny*nz
    real, allocatable :: reff_snow (:,:,:) ! nx*ny*nz
    real, allocatable :: reff_graup(:,:,:) ! nx*ny*nz

    ! sensor-dependent var
    real, allocatable :: zenith_angle (:,:)
    real, allocatable :: azimuth_angle(:,:)

    ! Geo var
    real, allocatable :: lat(:,:)
    real, allocatable :: lon(:,:)


  end type rt_result

contains

  subroutine rt_result_alloc( res, xbeg0, xend0, ybeg0, yend0, xbeg, xend, ybeg, yend, nz, nch, error_status)   
    type(rt_result), intent(inout) :: res
    integer, intent(in) :: xbeg0, xend0, ybeg0, yend0, nz, nch
    integer, intent(in) :: xbeg, xend, ybeg, yend
    integer, intent(out) :: error_status


    if (should_print(10)) then
      write(*,*) '--- Initializing rt_result ---'
      write(*,*) 'xbeg, xend: ', xbeg, xend
      write(*,*) 'ybeg, yend: ', ybeg, yend
      write(*,*) 'nz: ', nz
      write(*,*) '--- Initializing rt_result finish ---'
    end if

    if (res%l_allocated .eq. .true.) then
      write(*,*) 'Error in rt_result_module.f90 rt_result_alloc():'
      write(*,*) 'structure already allocated'
      error_status = FAILURE
      return
    end if

    res%xbeg0 = xbeg0
    res%xend0 = xend0
    res%ybeg0 = ybeg0
    res%yend0 = yend0
    res%xbeg = xbeg
    res%xend = xend
    res%ybeg = ybeg
    res%yend = yend
    res%nz   = nz
    res%nch  = nch

    allocate( res%TB           (xbeg:xend, ybeg:yend, nch), &
              res%emis         (xbeg:xend, ybeg:yend, nch), &
              res%zenith_angle (xbeg:xend, ybeg:yend), &
              res%azimuth_angle(xbeg:xend, ybeg:yend), &
              res%lat(xbeg:xend, ybeg:yend), &
              res%lon(xbeg:xend, ybeg:yend), &
              stat=error_status )
    res%TB(:,:,:)            = NF_FILL_FLOAT
    res%emis(:,:,:)          = NF_FILL_FLOAT
    res%zenith_angle(:,:)  = NF_FILL_FLOAT
    res%lat(:,:) = NF_FILL_FLOAT
    res%lon(:,:) = NF_FILL_FLOAT
    res%azimuth_angle = 999.9

    if (nml_l_output_reff) then
      allocate( &
        res%reff_cloud   (xbeg:xend, ybeg:yend, nz), &
        res%reff_rain    (xbeg:xend, ybeg:yend, nz), &
        res%reff_ice     (xbeg:xend, ybeg:yend, nz), &
        res%reff_snow    (xbeg:xend, ybeg:yend, nz), &
        res%reff_graup   (xbeg:xend, ybeg:yend, nz), &
        stat=error_status )
        
      res%reff_cloud(:,:,:)  = NF_FILL_FLOAT
      res%reff_rain (:,:,:)  = NF_FILL_FLOAT
      res%reff_ice  (:,:,:)  = NF_FILL_FLOAT
      res%reff_snow (:,:,:)  = NF_FILL_FLOAT
      res%reff_graup(:,:,:)  = NF_FILL_FLOAT
    end if
    res%l_allocated = .true.
  end subroutine rt_result_alloc
  !------------------------
  subroutine rt_result_alloc_from_state( res, state, nch, error_status)
    type(rt_result), intent(inout) :: res
    type(model_state), intent(in) :: state
    integer, intent(in) :: nch
    integer, intent(out) :: error_status
    call rt_result_alloc( res, state%xbeg0, state%xend0, state%ybeg0, state%yend0, &
                          state%xbeg, state%xend, state%ybeg, state%yend, &
                          state%nz, nch, error_status) 
    res%lat(res%xbeg:res%xend,res%ybeg:res%yend) = state%lat(res%xbeg:res%xend,res%ybeg:res%yend)
    res%lon(res%xbeg:res%xend,res%ybeg:res%yend) = state%lon(res%xbeg:res%xend,res%ybeg:res%yend)

  end subroutine rt_result_alloc_from_state
  !------------------------
  subroutine rt_result_dealloc(res, error_status)
    type(rt_result), intent(inout) :: res 
    integer, intent(out) :: error_status

    res%xbeg0 = -1
    res%xend0 = -1
    res%ybeg0 = -1
    res%yend0 = -1
    res%xbeg  = -1
    res%xend  = -1
    res%ybeg  = -1
    res%yend  = -1
    res%nz    = -1
    res%nch   = -1

    deallocate( res%TB, &
                res%emis, &
                res%zenith_angle, &
                res%azimuth_angle, &
                res%lat, &
                res%lon, &
                stat=error_status )
    if (nml_l_output_reff) then
      deallocate( &
        res%reff_cloud, &
        res%reff_rain , &
        res%reff_ice  , &
        res%reff_snow , &
        res%reff_graup, &
        stat=error_status )
    end if
    res%l_allocated = .false.
  end subroutine rt_result_dealloc

  ! -------------

  subroutine rt_result_output_nc(res, state, filename_out, error_status)
    type(rt_result),   intent(in) :: res
    type(model_state), intent(in) :: state
    character(len=*),  intent(in) :: filename_out
    integer,           intent(out) :: error_status
    
    integer :: fid, x_dimid, y_dimid, ch_dimid, t_dimid, z_dimid
    integer :: xlat_varid, xlon_varid, tb_varid, emis_varid, varid
    integer :: dimids(4), mystarts(4), mycounts(4)
    integer :: rcode
    integer, parameter :: ntime = 1
    character(len=256) :: filename_out_full
    integer :: i

    filename_out_full = trim(adjustl(filename_out)) // '.' // trim(adjustl(res%sensor_id)) // '.' // trim(adjustl(nml_s_rt_program)) // '.nc'
     
    call handle_err( nf_create_par( filename_out_full, IOR(IOR(nf_clobber, NF_NETCDF4), nf_mpiio), &
      MPI_COMM_WORLD,  MPI_INFO_NULL, fid) )

    ! define dimensions
    call handle_err( nf_def_dim(fid, 'x',    res%xend0-res%xbeg0+1, x_dimid), 'nf_def_dim(x)' )
    call handle_err( nf_def_dim(fid, 'y',    res%yend0-res%ybeg0+1, y_dimid), 'nf_def_dim(y)' )
    call handle_err( nf_def_dim(fid, 'z',    res%nz,                z_dimid), 'nf_def_dim(z)' )
    call handle_err( nf_def_dim(fid, 'ch',   res%nch,              ch_dimid), 'nf_def_dim(ch)' )
    call handle_err( nf_def_dim(fid, 'time', ntime,                 t_dimid), 'nf_def_dim(time)' )

    ! define output variables
    ! 1-D variables
    call nf_def_var_with_fill(fid, 'ch', NF_INT,   1,    (/ch_dimid/),     varid)
    ! channel-dependent variables
    dimids = (/x_dimid, y_dimid, ch_dimid, t_dimid/)
    call nf_def_var_with_fill(fid, 'Brightness_Temperature', NF_FLOAT, 4, dimids, tb_varid)
    call nf_def_var_with_fill(fid, 'Surface_Emissivity', NF_FLOAT, 4, dimids, emis_varid)
    ! 2-D var
    dimids(1:3) = (/x_dimid, y_dimid, t_dimid/)
    call nf_def_var_with_fill(fid, 'XLAT',          NF_FLOAT, 3, dimids(1:3), varid)
    call nf_def_var_with_fill(fid, 'XLONG',         NF_FLOAT, 3, dimids(1:3), varid)
    call nf_def_var_with_fill(fid, 'Zenith_Angle',  NF_FLOAT, 3, dimids(1:3), varid)
    call nf_def_var_with_fill(fid, 'Azimuth_Angle', NF_FLOAT, 3, dimids(1:3), varid)
    ! optional outputs
    ! reff 
    if (nml_l_output_reff) then
      dimids = (/x_dimid, y_dimid, z_dimid, t_dimid/)
      call nf_def_var_with_fill(fid, 'reff_cloud',  NF_FLOAT, 4, dimids, varid)
      call nf_def_var_with_fill(fid, 'reff_rain',   NF_FLOAT, 4, dimids, varid)
      call nf_def_var_with_fill(fid, 'reff_ice',    NF_FLOAT, 4, dimids, varid)
      call nf_def_var_with_fill(fid, 'reff_snow',   NF_FLOAT, 4, dimids, varid)
      call nf_def_var_with_fill(fid, 'reff_graup',  NF_FLOAT, 4, dimids, varid)
    end if
   
    call namelist_output_nc(fid, error_status)
    call handle_err( nf_enddef(fid))

    ! write variables to file
    ! channel-dependent variables
    mystarts(1) = res%xbeg - res%xbeg0 + 1
    mystarts(2) = res%ybeg - res%ybeg0 + 1
    mystarts(3) = 1
    mystarts(4) = 1
    mycounts(1) = res%xend - res%xbeg + 1
    mycounts(2) = res%yend - res%ybeg + 1
    mycounts(3) = res%nch
    mycounts(4) = 1
    if (should_print(10)) then
      write(*,*) 'mystarts: ', mystarts
      write(*,*) 'mycounts: ', mycounts
    end if
    call handle_err( nf_put_vara_real(fid, tb_varid, mystarts, mycounts, res%TB(res%xbeg:res%xend, res%ybeg:res%yend,:) ) )
    call handle_err( nf_put_vara_real(fid, emis_varid, mystarts, mycounts, res%emis(res%xbeg:res%xend, res%ybeg:res%yend,:) ) )
    ! ch
    rcode = nf_inq_varid(fid, 'ch', varid)
    call handle_err( nf_put_vara_int(fid, varid, mystarts(3:3), mycounts(3:3), res%channels(1:res%nch)) )
    ! 2-D var
    mystarts(1) = res%xbeg - res%xbeg0 + 1
    mystarts(2) = res%ybeg - res%ybeg0 + 1
    mystarts(3) = 1
    mycounts(1) = res%xend - res%xbeg + 1
    mycounts(2) = res%yend - res%ybeg + 1
    mycounts(3) = 1
    rcode = nf_inq_varid(fid, 'XLAT', varid)
    call handle_err( nf_put_vara_real(fid, varid, mystarts(1:3), mycounts(1:3), state%lat(res%xbeg:res%xend, res%ybeg:res%yend) ) )
    rcode = nf_inq_varid(fid, 'XLONG', varid)
    call handle_err( nf_put_vara_real(fid, varid, mystarts(1:3), mycounts(1:3), state%lon(res%xbeg:res%xend, res%ybeg:res%yend)) )
    rcode = nf_inq_varid(fid, 'Zenith_Angle', varid)
    call handle_err( nf_put_vara_real(fid, varid, mystarts(1:3), mycounts(1:3), res%zenith_angle(res%xbeg:res%xend, res%ybeg:res%yend) ) )
    rcode = nf_inq_varid(fid, 'Azimuth_Angle', varid)
    call handle_err( nf_put_vara_real(fid, varid, mystarts(1:3), mycounts(1:3), res%azimuth_angle(res%xbeg:res%xend, res%ybeg:res%yend) ) )
    if (nml_l_output_reff) then
      mystarts(1) = res%xbeg - res%xbeg0 + 1
      mystarts(2) = res%ybeg - res%ybeg0 + 1
      mystarts(3) = 1
      mystarts(4) = 1
      mycounts(1) = res%xend - res%xbeg + 1
      mycounts(2) = res%yend - res%ybeg + 1
      mycounts(3) = res%nz
      mycounts(4) = 1
      rcode = nf_inq_varid(fid, 'reff_cloud', varid)
      call handle_err( nf_put_vara_real(fid, varid, mystarts, mycounts, res%reff_cloud(res%xbeg:res%xend, res%ybeg:res%yend, :) ) )
      rcode = nf_inq_varid(fid, 'reff_rain',  varid)
      call handle_err( nf_put_vara_real(fid, varid, mystarts, mycounts, res%reff_rain(res%xbeg:res%xend, res%ybeg:res%yend, :)) )
      rcode = nf_inq_varid(fid, 'reff_ice', varid)
      call handle_err( nf_put_vara_real(fid, varid, mystarts, mycounts, res%reff_ice(res%xbeg:res%xend, res%ybeg:res%yend, :) ) )
      rcode = nf_inq_varid(fid, 'reff_snow', varid)
      call handle_err( nf_put_vara_real(fid, varid, mystarts, mycounts, res%reff_snow(res%xbeg:res%xend, res%ybeg:res%yend, :) ) )
      rcode = nf_inq_varid(fid, 'reff_graup', varid)
      call handle_err( nf_put_vara_real(fid, varid, mystarts, mycounts, res%reff_graup(res%xbeg:res%xend, res%ybeg:res%yend, :) ) )
    end if    
    
    call handle_err( nf_close(fid) )
  end subroutine rt_result_output_nc
  ! -------------------
  subroutine rt_result_convolution(res, state, obs_num, obs_ch, obs_x, obs_y, obs_efov_aScan, obs_efov_cscan, obs_azimuth_angle, tb_conv, error_status, reduce)
    type(rt_result),           intent(inout) :: res 
    type(model_state),         intent(inout) :: state
    integer, intent(in) :: obs_num
    integer, intent(in) :: obs_ch(:)
    real,    intent(in) :: obs_x(:)
    real,    intent(in) :: obs_y(:)
    real,    intent(in) :: obs_efov_aScan(:)
    real,    intent(in) :: obs_efov_cScan(:)
    real,    intent(in) :: obs_azimuth_angle(:)
    real,                      intent(inout) :: tb_conv(:)
    integer,                   intent(out)   :: error_status
    logical, intent(in), optional :: reduce

    real, parameter :: NaN = 0.0 / 0.0
    ! To make it simple, grab TB simulation for whole domain before convolution
    ! We calculate channel by channel so that it should not take too much memory
    integer :: ich
    real, allocatable :: tb(:,:) 
    real, allocatable :: xx(:,:), yy(:,:)
    real, allocatable :: w(:,:) 
    integer :: obs_beg, obs_end
    integer :: iobs
    integer :: xmin, xmax, ymin, ymax ! region that needs to calculate weights
    real :: obs_ii, obs_jj 
    real :: max_distance

    real, allocatable :: tb_conv_local(:)

    allocate( tb (res%xbeg0:res%xend0, res%ybeg0:res%yend0), &
              xx (res%xbeg0:res%xend0, res%ybeg0:res%yend0), &
              yy (res%xbeg0:res%xend0, res%ybeg0:res%yend0), &
              w  (res%xbeg0:res%xend0, res%ybeg0:res%yend0), &
              stat = error_status)

    allocate( tb_conv_local(obs_num), stat=error_status )
    tb_conv_local = 0.
    xx = 0.
    yy = 0.
    xx(res%xbeg:res%xend, res%ybeg:res%yend) = state%x(res%xbeg:res%xend, res%ybeg:res%yend)
    yy(res%xbeg:res%xend, res%ybeg:res%yend) = state%y(res%xbeg:res%xend, res%ybeg:res%yend)
    call MPI_Allreduce(MPI_IN_PLACE,xx, (res%xend0-res%xbeg0+1)*(res%yend0-res%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
    call MPI_Allreduce(MPI_IN_PLACE,yy, (res%xend0-res%xbeg0+1)*(res%yend0-res%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
    
!    if (my_proc_id==0) print *, 'channels', res%channels
    
    do ich = 1, res%nch
      tb = 0.
      tb(res%xbeg:res%xend, res%ybeg:res%yend) = res%tb (res%xbeg:res%xend, res%ybeg:res%yend, ich)
      call MPI_Allreduce(MPI_IN_PLACE,tb, (res%xend0-res%xbeg0+1)*(res%yend0-res%ybeg0+1),MPI_REAL,MPI_SUM,comm,error_status)
      ! each processor only process some of the observations
      call calc_proc_helper( nprocs, my_proc_id, 1, obs_num, obs_beg, obs_end )
      do iobs = obs_beg, obs_end
        if ( obs_ch(iobs) == res%channels(ich) ) then
          obs_ii = obs_x(iobs)
          obs_jj = obs_y(iobs)
          print *, 'obs_ii, obs_jj', obs_ii, obs_jj
          ! max_distance = (max(obs_efov_aScan(iobs), obs_efov_cScan(iobs))/state%dx * 2)  # Yinghui's definition of max_distance
          max_distance = ( ((obs_efov_aScan(iobs) + obs_efov_cScan(iobs)) / 2) / 1.8 ) / (state%dx * 2) ! # Scott's definition of max_distance
    
          xmin = int(max(res%xbeg0,   floor(obs_ii) - ceiling(max_distance)))
          xmax = int(min(res%xend0, ceiling(obs_ii) + ceiling(max_distance)))
          ymin = int(max(res%ybeg0,   floor(obs_jj) - ceiling(max_distance)))
          ymax = int(min(res%yend0, ceiling(obs_jj) + ceiling(max_distance)))

          if (xmin .le. xmax .and. ymin .le. ymax) then

            if ( any( tb(xmin:xmax,ymin:ymax) == NF_FILL_FLOAT ) ) then
              tb_conv_local(iobs) = NaN
            else
              w(:,:) = 0.
              print *, 'xmin, xmax, ymin, ymax', xmin, xmax, ymin, ymax
              call beam_conv_gaussian_simple(obs_ii, obs_jj, xmax-xmin+1, ymax-ymin+1, &
                      xx(xmin:xmax,ymin:ymax), yy(xmin:xmax,ymin:ymax), &
                      obs_efov_aScan(iobs), obs_efov_cScan(iobs), state%dx, &
                      w(xmin:xmax,ymin:ymax) )
              write(*,*) 'weights: ', w(xmin:xmax,ymin:ymax)
              write(*,*) 'sum of weights: ', sum( w(xmin:xmax,ymin:ymax) ) 
              !call gaussian_ellipse_weight(obs_ii, obs_jj, xmax-xmin+1, ymax-ymin+1, &
              !          xx(xmin:xmax,ymin:ymax), yy(xmin:xmax,ymin:ymax), &
              !          obs_efov_aScan(iobs), obs_efov_cScan(iobs), state%dx, obs_azimuth_angle(iobs), &
              !          w(xmin:xmax,ymin:ymax) )
              tb_conv_local(iobs) = sum(tb(xmin:xmax,ymin:ymax) * w(xmin:xmax,ymin:ymax) ) / sum( w(xmin:xmax,ymin:ymax) )
              write(*,*) 'tb_conv_local: ', tb_conv_local(iobs)

  !            print *,obs_ch(iobs) ,tb_conv_local(iobs) , sum(tb(xmin:xmax,ymin:ymax) * w(xmin:xmax,ymin:ymax) ) , sum( w(xmin:xmax,ymin:ymax) )
            end if
          end if
        end if
      end do ! observation loop
      ! send tb_conv cross processors
    
      
    end do ! channel loop
    if (.not. present(reduce) .or. reduce == .true.) then
      call MPI_Allreduce(MPI_IN_PLACE,tb_conv_local, obs_num,MPI_REAL,MPI_SUM,comm,error_status)
    end if
    where(tb_conv_local .ne. 0) tb_conv = tb_conv_local
    deallocate( tb, xx, yy, w, tb_conv_local )
  end subroutine rt_result_convolution
  
  ! -------------------
  subroutine rt_conv_output_txt(obs_info, tb_conv, filename_out, error_status, obs_mw)
    type(obs_info_type),   intent(in) :: obs_info
    real,              intent(in)  :: tb_conv(:)
    character(len=*),  intent(in)  :: filename_out
    integer,           intent(out) :: error_status
    type(microwave_data_type), intent(in), optional :: obs_mw
    
    character(len=256) :: filename_out_full

    integer :: i_sensor, i
    integer :: o_unit 
    integer :: iost
    
    o_unit = 11
    if (my_proc_id == 0) then
      do i_sensor = 1, obs_info%nsensors
        filename_out_full = trim(adjustl(filename_out)) // '.' // trim(adjustl(obs_info%sensors(i_sensor))) // '.' // trim(adjustl(nml_s_rt_program)) // '.conv.txt'
        print *, filename_out_full
        open(o_unit, file=filename_out_full, status='replace', form='formatted', iostat=iost)
        do i = 1, obs_mw%num
          if ( trim(adjustl(obs_mw%platform(i))) == trim(adjustl(obs_info%sensors(i_sensor))) ) then
!            write(*, *)trim(adjustl(obs_mw%platform(i))), trim(adjustl(obs_info%sensors(i_sensor))),  obs_mw%lat(i), obs_mw%lon(i)
            write(o_unit, '(2F10.3,I5,5F10.3)') obs_mw%lat(i), obs_mw%lon(i), obs_mw%ch(i), obs_mw%tb(i), tb_conv(i), obs_mw%efov_aScan(i), obs_mw%efov_cScan(i), obs_mw%azimuth_angle(i)
          end if
        end do
        close(o_unit)
      end do
    end if
  end subroutine rt_conv_output_txt




end module rt_result_module

