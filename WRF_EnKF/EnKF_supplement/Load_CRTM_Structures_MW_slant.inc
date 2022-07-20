  SUBROUTINE Load_CRTM_Structures_MW_slant(atm,sfc, x1,y1,zmax, zenith,azimuth,dx, debug_atm_data)

   USE CRTM_Module
   use qn2re_wsm6
   use qn2re_gfdlfv3
   use qn2re_thompson08
   
   type(CRTM_Atmosphere_type), intent(in out), dimension(1) :: atm
   type(CRTM_Surface_type), intent(in out), dimension(1) :: sfc

   integer, intent(in) :: x1, y1, zmax

   ! number of layers of ice to modify
   !integer, intent(in) :: N_CLOUDS
   real(fp), intent(in) :: zenith, azimuth
   real, intent(in) :: dx

   logical, intent(in) :: debug_atm_data

   integer :: cloud_num, z_wrf
   integer :: xx,yy,n,x_min,x_max,y_min,y_max

   real, dimension(zmax) :: qcloud_avg, qrain_avg, qice_avg, qsnow_avg, qgraup_avg, nrain_avg
   real, dimension(zmax-1) :: pressure_diff
   real, dimension(zmax)   :: level_pressure_diff

   real :: alt(0:zmax)
   real :: alt_wrf(0:zmax)
   real :: altitude
   real :: xoff, yoff
   real :: weight
   integer :: zz

   real :: T_air, rho_air

!*******************************************************************************
! satellite information
!*******************************************************************************

!  sat_dis = sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!  SCAN_ANGLE = 180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))**2))
!  ZENITH_ANGLE = SCAN_ANGLE+180.0/3.14159*acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!write(*,*) xlong(x,y),lon(x,y)
!write(*,*) lon(x,y)-sat_lon
!write(*,*) cos(lon(x,y)-sat_lon)
!write(*,*) acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
!write(*,*) 'SCAN_ANGLE',SCAN_ANGLE
!write(*,*) 'ZENITH_ANGLE',ZENITH_ANGLE

!*******************************************************************************
! load WRF data into CRTM structures
!*******************************************************************************


!--- calculating alt and delz at x1,y1

   n = 1
   alt(0) = (PH(x1,y1,zmax+1) + PHB(x1,y1,zmax+1)) /G
   do z = 1, zmax-1, 1
     z_wrf = zmax+1-z
     alt(z) = (PH(x1,y1,z_wrf) + PHB(x1,y1,z_wrf)) / G
     delz(z) = alt(z-1) - alt(z)
   end do
   alt(zmax) = hgt(x1,y1)
   delz(z) = alt(zmax-1) - alt(zmax)
   if (delz(zmax) <= 0.) then
     alt(zmax) = alt(zmax-1)-1
     delz(zmax) = 1
   end if

!---Writing Atmospheric Profile
   atm(n)%Climatology         = TROPICAL
   atm(n)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
   atm(n)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
   atm(n)%Level_Pressure(:) = 0.
   atm(n)%Pressure(:) = 0.
   atm(n)%Temperature(:) = 0.
   atm(n)%Absorber(:,:) = 0.
   qcloud_avg(:) = 0.
   qrain_avg(:) = 0.
   qice_avg(:) = 0.
   qsnow_avg(:) = 0.
   qgraup_avg(:) = 0.
   nrain_avg(:) = 0.

   do z=1,zmax,1

     yoff = (alt(z)+alt(z-1))/2 * TAN( zenith / 180. * PI ) * COS( azimuth / 180. * PI) / dx / 1000
     xoff = (alt(z)+alt(z-1))/2 * TAN( zenith / 180. * PI ) * SIN( azimuth / 180. * PI) / dx / 1000

     x_min = x1 + int(floor( xoff ))
     x_max = x1 + int(floor( xoff )) + 1
     if (x_min<1 ) x_min = 1
     if (x_min>ix) x_min = ix
     if (x_max<1 ) x_max = 1
     if (x_max>ix) x_max = ix

     y_min = y1 + int(floor( yoff ))
     y_max = y1 + int(floor( yoff )) + 1
     if (y_min<1 ) y_min = 1
     if (y_min>jx) y_min = jx
     if (y_max<1 ) y_max = 1
     if (y_max>jx) y_max = jx

     !print *, 'proc_id', my_proc_id, 'x_min', x_min, ' y_min', y_min, 'x_max', x_max, 'y_max', y_max
     do xx = x_min, x_max
     do yy = y_min, y_max

       weight = 1
       if (x_min .ne. x_max) weight = weight * abs( x1+xoff -xx )
       if (y_min .ne. y_max) weight = weight * abs( y1+yoff -yy )

       ! altitude of wrf grid at (x,y)
       do zz = 0, zmax, 1
         z_wrf = zmax+1-zz
         alt_wrf(zz) = ((PH(xx,yy,z_wrf) + PHB(xx,yy,z_wrf) + PH(xx,yy,z_wrf-1) + PHB(xx,yy,z_wrf-1))) / G / 2.
       end do
       alt_wrf(zmax) = hgt(xx,yy)
       
       ! get altitude at z (unstaggered grid)
       altitude = ( alt(z-1) + alt(z) )/2.

       ! if altitude is above (below) maximum (minimum) altitude at (x,y), then
       ! use the atmosphere at (x1,y1)
       if (altitude > alt_wrf(0) .or. altitude < alt_wrf(zmax) ) then
         z_wrf = zmax-z+1
         atm(n)%Pressure(z)       = atm(n)%Pressure(z)    + weight * pres(x1,y1,z_wrf) / 100.    ! convert from Pa to hPA
         atm(n)%Temperature(z)    = atm(n)%Temperature(z) + weight * tk(x1,y1,z_wrf)
         atm(n)%Absorber(z,1)     = atm(n)%Absorber(z,1)  + weight * qvapor(x1,y1,z_wrf)*1000.0  ! convert from "kg kg-1" to "g kg-1"; remove negatives

         qcloud_avg(z) = qcloud_avg(z) + weight * qcloud(x1,y1,z_wrf)
         qrain_avg(z)  = qrain_avg(z)  + weight *  qrain(x1,y1,z_wrf)
         qice_avg(z)   = qice_avg(z)   + weight *   qice(x1,y1,z_wrf)
         qsnow_avg(z)  = qsnow_avg(z)  + weight *  qsnow(x1,y1,z_wrf)
         qgraup_avg(z) = qgraup_avg(z) + weight * qgraup(x1,y1,z_wrf)
         nrain_avg(z)  = nrain_avg(z)  + weight *  nrain(x1,y1,z_wrf)
       else
         ! find the layer at (x,y) which has an altitude closest to that of the
         ! layer z
         z_wrf = -10
         do zz = 1, zmax
           if ( altitude > alt_wrf(zz) .and. altitude <= alt_wrf(zz-1) ) then
             z_wrf = zmax - zz + 1
             exit
           end if
         end do

         atm(n)%Pressure(z)       = atm(n)%Pressure(z)    + weight * pres(xx,yy,z_wrf) / 100.    ! convert from Pa to hPA
         atm(n)%Temperature(z)    = atm(n)%Temperature(z) + weight * tk(xx,yy,z_wrf)
         atm(n)%Absorber(z,1)     = atm(n)%Absorber(z,1)  + weight * qvapor(xx,yy,z_wrf)*1000.0  ! convert from "kg kg-1" to "g kg-1"; remove negatives

         qcloud_avg(z) = qcloud_avg(z) + weight * qcloud(xx,yy,z_wrf)
         qrain_avg(z)  = qrain_avg(z)  + weight *  qrain(xx,yy,z_wrf)
         qice_avg(z)   = qice_avg(z)   + weight *   qice(xx,yy,z_wrf)
         qsnow_avg(z)  = qsnow_avg(z)  + weight *  qsnow(xx,yy,z_wrf)
         qgraup_avg(z) = qgraup_avg(z) + weight * qgraup(xx,yy,z_wrf)
         nrain_avg(z)  = nrain_avg(z)  + weight *  nrain(xx,yy,z_wrf)
       end if
     enddo
     enddo ! loop across the zone in x and y which is in the slant path 
   enddo ! loop vertically 

   ! Check if pressure do not increase as it goes down
   ! If not, use column at (x1, y1) instead
   pressure_diff = atm(n)%Pressure(zmax:2:-1) - atm(n)%Pressure(zmax-1:1:-1)

   ! Bad TOA pressure diffs are the most common problem it seems. This is a
   ! bodge.
   if (pressure_diff(zmax-1) .LT. 1) then
     !write(*,*) x1, y1, 'TOA pressures: ', atm(n)%Pressure(2), atm(n)%Pressure(1) 
     atm(n)%Pressure(1) = atm(n)%Pressure(2) - 1
     pressure_diff(zmax-1) = 1
   endif

   ! If bad pressure diffs elsewhere, this is an even more drastic bodge.
   if (ANY(pressure_diff .LE. 0.99)) then
     !write(*,*) x1, y1, ' had their slant path removed due to negative pressure diffs'
     !write(*,*) pressure_diff
!     sfc(n)%land_coverage = 1.5
!     return
     do z=1,zmax,1
       z_wrf = zmax-z+1
       atm(n)%Pressure(z)       = pres(x1,y1,z_wrf) / 100.
! convert from Pa to hPA
       atm(n)%Temperature(z)    = tk(x1,y1,z_wrf)
       atm(n)%Absorber(z,1)     = qvapor(x1,y1,z_wrf)*1000.0
! convert from "kg kg-1" to "g kg-1"; remove negatives

       qcloud_avg(z) = qcloud(x1,y1,z_wrf)
       qrain_avg(z)  =  qrain(x1,y1,z_wrf)
       qice_avg(z)   =   qice(x1,y1,z_wrf)
       qsnow_avg(z)  =  qsnow(x1,y1,z_wrf)
       qgraup_avg(z) = qgraup(x1,y1,z_wrf)
       nrain_avg(z)  =  nrain(x1,y1,z_wrf)
     end do

     ! If the vertical column has a bad TOA pressure diff, then correct it
     pressure_diff = atm(n)%Pressure(zmax:2:-1) - atm(n)%Pressure(zmax-1:1:-1)
     if (pressure_diff(zmax-1) .LT. 1) then
       write(*,*) x1, y1, 'post-adjustment TOA pressures: ', atm(n)%Pressure(2), atm(n)%Pressure(1) 
       atm(n)%Pressure(1) = atm(n)%Pressure(2) - 1
       pressure_diff(zmax-1) = 1
     endif
   endif

   ! put level pressures as the halfway mark between layer pressures (except at
   ! the vertical boundaries)
   atm(n)%Level_Pressure(0) = atm(n)%Pressure(1) * 3.0/2.0 - atm(n)%Pressure(2)/2.0
   do z = 1, zmax-1
     atm(n)%Level_Pressure(z) = (atm(n)%Pressure(z) + atm(n)%Pressure(z+1) )/2.0
   end do
   atm(n)%Level_Pressure(zmax) = max(psfc(x1,y1), pres(x1,y1,1)*3.0/2.0-pres(x1,y1,2)/2.0)/100.0

   ! see if 
   if (ANY(atm(n)%Level_Pressure(:) .LT. 0.01)) then
     write(*,*) x1, y1, 'post_adjustment level pressures negative: '
     write(*,*) atm(n)%Level_Pressure(0:2)
     atm(n)%Level_Pressure(0) = max(0.01,atm(n)%Level_Pressure(2) - atm(n)%Level_Pressure(1))
     if (atm(n)%Level_Pressure(0) < 0.011) then
       write(*,*) '!!!!!',x1,y1,'Level_Pressure(0) had to be set to 0.01!!!!!'
     endif
     if (ANY(atm(n)%Level_Pressure(:) .LT. 0.01)) then
       write(*,*) '!!!!!',x1,y1,'Level_Pressure(0) not the only negative value!!!!!'
     endif
   endif

   level_pressure_diff = atm(n)%Level_Pressure(zmax:1:-1) - atm(n)%Level_Pressure(zmax-1:0:-1)
   if (ANY(level_pressure_diff .LE. 0.99)) then
     write(*,*) x1, y1, 'post_adjustment level pressures diffs: '
     write(*,*) level_pressure_diff
   endif

   !pressure_diff = atm(n)%Pressure(zmax:2:-1) - atm(n)%Pressure(zmax-1:1:-1)
   !if (ANY(pressure_diff .LE. 0.99)) then

!--- writing Cloud data
!
 
   do z=1,N_CLOUDS
     atm(n)%Cloud(z)%Type = 0
     atm(n)%Cloud(z)%Water_Content = 0.0
     atm(n)%Cloud(z)%Effective_Radius = 0.0
   enddo

   icl = 0
!---calculating # of clouds (cloud and rain)
   IF ( atm(n)%n_Clouds > 0 ) THEN
!---Data for cloud

   icl = 1
   atm(n)%Cloud(icl)%Type = WATER_CLOUD
   atm(n)%Cloud(icl)%Effective_Radius(:) = 0
   atm(n)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qcloud_avg(z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = atm(n)%Pressure(z)*100/287.2/(atm(n)%Temperature(z)+0.61*((atm(n)%Absorber(z,1)/1000.)/(1+(atm(n)%Absorber(z,1)/1000.))))
!       if (l_use_default_re) then
         atm(n)%Cloud(icl)%Effective_Radius(k1) =  16.8_fp   ! water cloud effective radius does not matter to the CRTM
!       else
!         atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_cloud(qcloud_avg(z), rho_air) *1e6 !  16.8_fp * radius_ratio
!       end if
       atm(n)%Cloud(icl)%Water_Content(k1)    = &
           qcloud_avg(z)* rho_air *delz(z)
     endif
   enddo

   icl = 2
   atm(n)%Cloud(icl)%Type = RAIN_CLOUD
   atm(n)%Cloud(icl)%Effective_Radius(:) = 0
   atm(n)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qrain_avg(z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = atm(n)%Pressure(z)*100/287.2/(atm(n)%Temperature(z)+0.61*((atm(n)%Absorber(z,1)/1000.)/(1+(atm(n)%Absorber(z,1)/1000.))))
       if (l_use_default_re) then
         atm(n)%Cloud(icl)%Effective_Radius(k1) =  1000.0_fp 
       else
         if (mp_scheme == 0) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_rain(qrain_avg(z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_rain(qrain_avg(z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_rain(qrain_avg(z), rho_air, nrain_avg(z)) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(n)%Cloud(icl)%Water_Content(k1)    = &
           qrain_avg(z)* rho_air *delz(z)

       !write (*,*) x, y, z, atm(n)%Cloud(icl)%Effective_Radius(k1)
     endif
   enddo

   icl = 3
   atm(n)%Cloud(icl)%Type = ICE_CLOUD
   atm(n)%Cloud(icl)%Effective_Radius(:) = 0
   atm(n)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qice_avg(z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = atm(n)%Pressure(z)*100/287.2/(atm(n)%Temperature(z)+0.61*((atm(n)%Absorber(z,1)/1000.)/(1+(atm(n)%Absorber(z,1)/1000.))))
!       if (l_use_default_re) then
         atm(n)%Cloud(icl)%Effective_Radius(k1) =  25.0_fp   ! cloud ice effective radius does not matter to the CRTM
!       else
!         atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_ice(qice_avg(z), rho_air) *1e6 ! 125.0_fp * radius_ratio
!       end if
       
       atm(n)%Cloud(icl)%Water_Content(k1)    = & 
           qice_avg(z)* rho_air *delz(z)
     endif
   enddo

   icl = 4
   atm(n)%Cloud(icl)%Type = SNOW_CLOUD
   atm(n)%Cloud(icl)%Effective_Radius(:) = 0
   atm(n)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qsnow_avg(z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = atm(n)%Pressure(z)*100/287.2/(atm(n)%Temperature(z)+0.61*((atm(n)%Absorber(z,1)/1000.)/(1+(atm(n)%Absorber(z,1)/1000.))))
       if (l_use_default_re) then
         atm(n)%Cloud(icl)%Effective_Radius(k1) = 500.0_fp 
       else
         T_air = atm(n)%Temperature(k1)
         if (mp_scheme == 0) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_snow(qsnow_avg(z), rho_air, T_air)  *1e6!500.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_snow(qsnow_avg(z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_snow(qsnow_avg(z), rho_air, T_air) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(n)%Cloud(icl)%Water_Content(k1)    = & 
           qsnow_avg(z)* rho_air *delz(z)
     endif
   enddo

   icl = 5
   atm(n)%Cloud(icl)%Type = GRAUPEL_CLOUD
   atm(n)%Cloud(icl)%Effective_Radius(:) = 0
   atm(n)%Cloud(icl)%Water_Content(:)    = 0
   do z=zmax,1,-1
     if(qgraup_avg(z).gt.0.0) then
       k1 = zmax-z+1
       rho_air = atm(n)%Pressure(z)*100/287.2/(atm(n)%Temperature(z)+0.61*((atm(n)%Absorber(z,1)/1000.)/(1+(atm(n)%Absorber(z,1)/1000.))))
       if (l_use_default_re) then
         atm(n)%Cloud(icl)%Effective_Radius(k1) = 1500.0_fp 
       else
         if (mp_scheme == 0) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_WSM6_graup(qgraup_avg(z), rho_air)  *1e6!500.0_fp * radius_ratio
         else if (mp_scheme == 1) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_GFDLFV3_graup(qgraup_avg(z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         else if (mp_scheme == 2) then
           atm(n)%Cloud(icl)%Effective_Radius(k1) = qn2re_Thompson08_graup(qgraup_avg(z), rho_air) *1e6 ! 1000.0_fp * radius_ratio
         end if
       end if
       atm(n)%Cloud(icl)%Water_Content(k1)    = & 
           qgraup_avg(z)* rho_air *delz(z)
     endif
   enddo
   ENDIF

   ! remove negative water content
   do icl = 1, 5
     where (atm(n)%Cloud(icl)%Water_Content < 0) atm(n)%Cloud(icl)%Water_Content = 0.
   end do

!   print *, maxval(atm(n)%Cloud(1)%Effective_Radius(:)) &
!          , maxval(atm(n)%Cloud(2)%Effective_Radius(:)) &
!          , maxval(atm(n)%Cloud(3)%Effective_Radius(:)) &
!          , maxval(atm(n)%Cloud(4)%Effective_Radius(:)) &
!          , maxval(atm(n)%Cloud(5)%Effective_Radius(:)) 
!*******************************************************************************
! Surface data
!*************
   sfc(n)%Wind_speed     = windspeed(x1,y1)
   sfc(n)%Wind_direction = winddir(x1,y1)
   if(landmask(x1,y1).eq.1.0) then
    sfc(n)%Water_Coverage = 0.0_fp
    sfc(n)%Land_Coverage = 1.0_fp
    sfc(n)%Land_Temperature = tsk(x1,y1)
    sfc(n)%Soil_Temperature = tsk(x1,y1)
   else
    sfc(n)%Water_Coverage = 1.0_fp
    sfc(n)%Land_Coverage = 0.0_fp
    sfc(n)%Water_Type = 1  ! Sea water
    sfc(n)%Water_Temperature = tsk(x1,y1)
   endif


  END SUBROUTINE load_crtm_structures_mw_slant

