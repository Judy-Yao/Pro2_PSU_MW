module rt_util_module
  use netcdf
  use rt_constant_module
  implicit none
  

contains
  subroutine handle_err(errcode, fname)
    !include 'netcdf.inc'
    integer, intent(in) ::  errcode
    character(len=*), intent(in), optional :: fname
    if (errcode .ne. 0) then
      if (present(fname)) then
        print *, 'NETCDF Error: ', fname, nf_strerror(errcode)
      else
        print *, 'NETCDF Error: ', nf_strerror(errcode)
      end if
      stop 2
    end if
  end subroutine handle_err


!---------
  subroutine getxyzmax_wrfout(filename, xmax, ymax, zmax)
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: xmax
    integer,          intent(out) :: ymax
    integer,          intent(out) :: zmax

    integer :: ncid
    integer :: dimid
    
    call handle_err ( NF_OPEN(filename, NF_NOWRITE, ncid) )
    call handle_err ( NF_INQ_DIMID (NCID, 'west_east',   dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, xmax) )
    call handle_err ( NF_INQ_DIMID (NCID, 'south_north', dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, ymax) )
    call handle_err ( NF_INQ_DIMID (NCID, 'bottom_top',  dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, zmax) )
  end subroutine getxyzmax_wrfout
!---------
  subroutine getxyzmax_fv3restart(filename, xmax, ymax, zmax)
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: xmax
    integer,          intent(out) :: ymax
    integer,          intent(out) :: zmax

    integer :: ncid
    integer :: dimid
    
    call handle_err ( NF_OPEN(filename, NF_NOWRITE, ncid) )
    call handle_err ( NF_INQ_DIMID (NCID, 'xaxis_1',   dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, xmax) )
    call handle_err ( NF_INQ_DIMID (NCID, 'yaxis_1', dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, ymax) )
    call handle_err ( NF_INQ_DIMID (NCID, 'zaxis_1',  dimid) )
    call handle_err ( NF_INQ_DIMLEN(NCID, dimid, zmax) )
  end subroutine getxyzmax_fv3restart
!---------
  subroutine nf_def_var_with_fill(fid, fieldname,  datatype, ndim, dimids, varid) 
    integer, intent(in) :: fid
    character(len=*), intent(in) :: fieldname
    integer, intent(in) :: datatype
    integer, intent(in) :: ndim
    integer, intent(in) :: dimids(:)
    integer, intent(out) :: varid

    call handle_err( nf_def_var(fid, fieldname,  datatype, ndim, dimids, varid) )
    if (datatype == NF_FLOAT) then
      call handle_err( nf_put_att_real(fid, varid, '_FillValue', NF_FLOAT, 1, NF_FILL_FLOAT) )
    end if

    !call handle_err( nf_def_var_deflate(fid, varid, 0, 1, 1) )

!     STATUS = NF_PUT_ATT_DOUBLE (NCID, RHID, 'valid_range', NF_DOUBLE, &
!                                        2, RHRNGE)
  end subroutine nf_def_var_with_fill


!---------
  subroutine beam_conv_gaussian_simple(x0, y0, nx, ny, x, y, efov_a, efov_c, dx, weight)
    real, intent(in) :: x0, y0 ! (x0,y0) for the center of ellipse
    integer, intent(in) :: nx, ny
    real, intent(in) :: x(nx,ny),  y(nx,ny)  ! (x, y)  for the location in consideration
    real, intent(in) :: efov_a, efov_c ! effective FOV along/cross track in km
    real, intent(in) :: dx ! model grid spacing in km
    real, intent(inout) :: weight(nx,ny) ! weighting at the location

    !real :: la, lc
    real :: sigma, distance(nx,ny), distance_x(nx,ny), distance_y(nx,ny)

    sigma = 0.5 * ((((efov_a+efov_c)/2) / 1.18) / dx)

    distance_x = x-x0
    distance_y = y-y0
    distance = sqrt(distance_x**2 + distance_y**2)    

    weight = exp(-1*(distance**2)) / (2*sigma**2)
    
  end subroutine beam_conv_gaussian_simple

!---------
  subroutine gaussian_ellipse_weight(x0, y0, nx, ny, x, y, efov_a, efov_c, dx, azi, w)
    real, intent(in) :: x0, y0 ! (x0,y0) for the center of ellipse
    integer, intent(in) :: nx, ny
    real, intent(in) :: x(nx,ny),  y(nx,ny)  ! (x, y)  for the location in consideration
    real, intent(in) :: efov_a, efov_c ! effective FOV along/cross track in km
    real, intent(in) :: dx ! model grid spacing in km
    real, intent(in) :: azi ! observation azimuth angle in degrees
    real, intent(inout) :: w(nx,ny) ! weighting at the location

    !real :: la, lc
    real, parameter :: cc = (2*1.18)**2 /2.
    
    !la = sind(azi)*(x-x0) + cosd(azi)*(y-y0)
    !lc = cosd(azi)*(x-x0) - sind(azi)*(y-y0)
    !w = exp( - ( la**2/efov_a**2 + lc**2/efov_c**2) * cc * dx**2) / (PI/cc) / (efov_a/dx * efov_c/dx)
    
    w = exp( - ( (sind(azi)*(x-x0) + cosd(azi)*(y-y0))**2/efov_a**2 + (cosd(azi)*(x-x0) - sind(azi)*(y-y0))**2/efov_c**2) * cc * dx**2) &
      / (PI/cc) / (efov_a/dx * efov_c/dx)
  end subroutine gaussian_ellipse_weight

!---------
  function get_default_sensor_zenith(sensor_id)
    character(len=*), intent(In) :: sensor_id
    real :: get_default_sensor_zenith

    if (sensor_id == 'gmi_gpm_hf') then
      get_default_sensor_zenith = 49.1
    else if (sensor_id == 'gmi_gpm_lf') then
      get_default_sensor_zenith = 52.8
    else if (sensor_id(1:5) == 'amsr2') then
      get_default_sensor_zenith = 55
    else if (sensor_id(1:5) == 'ssmis') then
      get_default_sensor_zenith = 53.1
    else if (sensor_id(1:4) == 'ssmi') then
      get_default_sensor_zenith = 53.1
    else
      get_default_sensor_zenith = 0.0
    end if

  end function

end module rt_util_module
