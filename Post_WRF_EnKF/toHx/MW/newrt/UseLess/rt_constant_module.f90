module rt_constant_module
  implicit none
  save

  integer, parameter :: SUCCESS = 0
  integer, parameter :: FAILURE = 1
  
  !real, parameter :: PI=4.0 * atan(1.0)
  real, parameter :: PI=3.14159
  real, parameter :: deg2rad = PI/180.
  real, parameter :: rad2deg = 180./PI

  integer, parameter :: LENGTH_SENSOR_NAME = 20
  integer, parameter :: MAX_SENSORS  = 10
  integer, parameter :: MAX_CHANNELS = 50

  real, parameter :: DEFAULT_AZIMUTH_CRTM = 999.9
  real, parameter :: DEFAULT_AZIMUTH_RTTOV = 0.0

end module rt_constant_module
