! =========================================================================
!              PARALLELIZED OFFLINE GEOSTATIONARY IR CRTM CODE
! =========================================================================
! Written by M.-Y. Chan
! Based on xb_to_radiance subroutine in EnSRF/xb.f
! =========================================================================
! IMPORTANT NOTES:
! 1) Always check if the desired satellite and channels are set up properly.
! 2) To run after compilation (crtm.exe):
!    >>> PARALLEL_RUN crtm.exe wrf_file crtm_output_name.bin
!    where PARALLEL_RUN can be mpirun, srun, ibrun (whichever suitable)
! =========================================================================



PROGRAM XbtoIR

  ! -----------------------------------------------------------------------
  ! 1. Modules used in this program
  ! -----------------------------------------------------------------------
  USE netcdf
  USE mpi_module
  USE CRTM_Module

  ! -----------------------------------------------------------------------
  ! 2. Predefining variables (no implicit variables)
  ! -----------------------------------------------------------------------

  implicit none

  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'XbtoIR'

  ! CRTM coefficient directory path
  CHARACTER(*), PARAMETER :: CRTM_COEFF_DIR = 'coefficients/'

  ! Satellite-related constants (Tuned to GOES-16 for 2017 before Nov)
  REAL, PARAMETER :: Re      = 6378000.0
  REAL, PARAMETER :: sat_h   = 35780000.0 !35786000.0 
  REAL, PARAMETER :: sat_lon = -89.5/180.0*3.14159 !0./180*3.14159 
  INTEGER, PARAMETER :: n_ch = 1
  INTEGER, DIMENSION(n_ch) :: chlist = (/8/)
  CHARACTER(256) :: Sensor_Id = ADJUSTL('abi_gr')   ! what is the difference between abi_gr and abi_g16


  ! Input and output file names
  CHARACTER(256) :: wrf_fname
  CHARACTER(256) :: out_fname

  ! WRF domain scalars
  INTEGER        :: xmax, ymax, zmax

  ! WRF n-dimensional variables
  REAL, allocatable, dimension(:) :: delz
  REAL, allocatable, dimension(:,:) :: xlat, xlong, lat, lon, psfc, hgt, tsk, landmask
  REAL, allocatable, dimension(:,:,:) :: p, pb, pres, ph, phb, t, tk
  REAL, allocatable, dimension(:,:,:) :: qvapor, qcloud, qrain, qice, qsnow, qgraup

  ! CRTM constants 
  REAL, PARAMETER :: P1000MB=100000.D0
  REAL, PARAMETER :: R_D=287.D0
  REAL, PARAMETER :: Cpd=7.D0*R_D/2.D0

  ! CRTM profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1 
  INTEGER, PARAMETER :: N_ABSORBERS = 2 
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  INTEGER :: N_LAYERS, N_CLOUDS

  ! Other CRTM variables
  REAL(fp) :: ZENITH_ANGLE, SCAN_ANGLE, sat_dis
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, m, irec, pt_end, pt_start, n_pts
  integer :: ncid,ncrcode
  integer :: x,y,tt,v,z,n,reci,ens,n_ec
  integer :: pt_id, tot_pts
  INTEGER :: ncl,icl,k1,k2
  real, allocatable, dimension(:,:,:) :: BTsend, BT

  ! CRTM structures
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)

  ! -----------------------------------------------------------------------
  ! 3. Fill in domain variables and allocate arrays
  ! -----------------------------------------------------------------------
  ! Read in file names
  call getarg( 1, wrf_fname )
  call getarg( 2, out_fname )
  
  ! Determine WRF dimensions
  call get_ij( wrf_fname, xmax, ymax, zmax)
  N_LAYERS = zmax
  N_CLOUDS = zmax*5

  ! Allocate arrays
  allocate(  xlat(xmax,ymax)  )  ! latitude
  allocate(  xlong(xmax,ymax) )  ! longitude
  allocate(  lat(xmax,ymax)   )  ! in radian
  allocate(  lon(xmax,ymax)   )  ! in radian
  allocate(  p(xmax,ymax,zmax)  )
  allocate(  pb(xmax,ymax,zmax)  )
  allocate(  pres(xmax,ymax,zmax)  )
  allocate(  ph(xmax,ymax,zmax+1)  )
  allocate(  phb(xmax,ymax,zmax+1)  )
  allocate(  delz(zmax)  )
  allocate(  t(xmax,ymax,zmax)  )
  allocate(  tk(xmax,ymax,zmax)  )
  allocate(  qvapor(xmax,ymax,zmax)  )
  allocate(  qcloud(xmax,ymax,zmax)  )
  allocate(  qrain(xmax,ymax,zmax)  )
  allocate(  qice(xmax,ymax,zmax)  )
  allocate(  qsnow(xmax,ymax,zmax)  )
  allocate(  qgraup(xmax,ymax,zmax)  )
  allocate(  psfc(xmax,ymax)  )
  allocate(  hgt(xmax,ymax)  )
  allocate(  tsk(xmax,ymax)  )
  allocate(  landmask(xmax,ymax)  )

  ! Initialize parallization
  call parallel_start()

 
  ! -----------------------------------------------------------------------
  ! 4. Initialize the CRTM for the sensor
  ! -----------------------------------------------------------------------
  Error_Status = CRTM_Init( (/Sensor_Id/), &  ! Input... must be an array, hencethe (/../)
                            ChannelInfo  , &  ! Output
                            IRwaterCoeff_File='WuSmith.IRwater.EmisCoeff.bin',&
                            IRlandCoeff_File='IGBP.IRland.EmisCoeff.bin',&
                            File_Path= CRTM_COEFF_DIR, &
                            Quiet=.true.)
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  
  ! -----------------------------------------------------------------------
  ! 5. Initialize CRTM sensor channels
  ! -----------------------------------------------------------------------
  Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset = chlist )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM sensor channels'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  allocate(BT(xmax,ymax,n_Channels))
  allocate(BTsend(xmax,ymax,n_Channels))

  ! Zero things ahead of time
  BT = 0.
  BTsend = 0.


  ! -----------------------------------------------------------------------
  ! 6. Allocate CRTM structure arrays
  ! -----------------------------------------------------------------------

  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF

  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF


  ! -----------------------------------------------------------------------
  ! 7. Import and process WRF data
  ! ----------------------------------------------------------------------- 
  ! Importing WRF data
  call get_variable2d(wrf_fname,'XLAT',xmax,ymax,1,xlat)
  call get_variable2d(wrf_fname,'XLONG',xmax,ymax,1,xlong)
  call get_variable3d(wrf_fname,'P',xmax,ymax,zmax,1,p)
  call get_variable3d(wrf_fname,'PB',xmax,ymax,zmax,1,pb)
  call get_variable3d(wrf_fname,'PH',xmax,ymax,zmax+1,1,ph)
  call get_variable3d(wrf_fname,'PHB',xmax,ymax,zmax+1,1,phb)
  call get_variable3d(wrf_fname,'T',xmax,ymax,zmax,1,t)
  call get_variable3d(wrf_fname,'QVAPOR',xmax,ymax,zmax,1,qvapor)
  call get_variable2d(wrf_fname,'PSFC',xmax,ymax,1,psfc)
  call get_variable2d(wrf_fname,'TSK',xmax,ymax,1,tsk)
  call get_variable2d(wrf_fname,'HGT',xmax,ymax,1,hgt)
  call get_variable3d(wrf_fname,'QCLOUD',xmax,ymax,zmax,1,qcloud)
  call get_variable3d(wrf_fname,'QRAIN',xmax,ymax,zmax,1,qrain)
  call get_variable3d(wrf_fname,'QICE',xmax,ymax,zmax,1,qice)
  call get_variable3d(wrf_fname,'QSNOW',xmax,ymax,zmax,1,qsnow)
  call get_variable3d(wrf_fname,'QGRAUP',xmax,ymax,zmax,1,qgraup)
  call get_variable2d(wrf_fname,'XLAND',xmax,ymax,1,landmask)

  !!!!!!!!!! clear sky
  !qcloud = 0.0
  !qrain = 0.0
  !qice = 0.0
  !qsnow = 0.0
  !qgraup = 0.0

  ! Processing WRF data
  lat = xlat/180.0*3.14159
  lon = xlong/180.0*3.14159
  pres = P + PB
  tk = (T + 300.0) * ( (pres / P1000MB) ** (R_D/Cpd) )

  ! Removing spurious values
  where(qvapor.lt.0.0) qvapor=1.0e-8
  where(qcloud.lt.0.0) qcloud=0.0
  where(qice.lt.0.0) qice=0.0
  where(qrain.lt.0.0) qrain=0.0
  where(qsnow.lt.0.0) qsnow=0.0
  where(qgraup.lt.0.0) qgraup=0.0


  ! -----------------------------------------------------------------------
  ! 8. Set up parallelization
  ! ----------------------------------------------------------------------- 

  ! Total number of points on grid
  tot_pts = xmax * ymax

  ! Determine how many points to allocate to each process
  if(mod(tot_pts,nprocs).eq.0) then
     n_pts = tot_pts / nprocs
  else
     n_pts = tot_pts / nprocs+1
  endif

  ! Determine range of points handled by each process
  pt_start=my_proc_id*n_pts+1
  pt_end=min(tot_pts,(my_proc_id+1)*n_pts)


  ! -----------------------------------------------------------------------
  ! 9. Running the CRTM on each domain point
  ! ----------------------------------------------------------------------- 

  ! For each domain point specified to each process
  do pt_id = pt_start, pt_end

    ! Determine coordinate index of the domain point
    x = mod( pt_id-1, xmax )+1
    y = int( (pt_id-1) / xmax )+1
    
  
    ! Satellite viewing geometry
    sat_dis=sqrt(Re**2.0+(Re+sat_h)**2.0-2.0*Re*(Re+sat_h)*cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
    SCAN_ANGLE=180.0/3.14159*asin(Re/sat_dis*sqrt(1-(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))**2))
    ZENITH_ANGLE=SCAN_ANGLE+180.0/3.14159*acos(cos(lon(x,y)-sat_lon)*cos(lat(x,y)))
    if ( abs( ZENITH_ANGLE ) >= 80 ) then
      write(*,*) x, y, lon(x,y), lat(x,y)
      write(*,*) pt_id, sat_dis, SCAN_ANGLE, ZENITH_ANGLE
      
      stop
    endif
  
    ! --------------------------------------------------------------------
    ! 9-1. Load WRF data into CRTM structures
    ! --------------------------------------------------------------------
    ! calculating delz
    do z=1,zmax
      if(z.eq.1) then
        delz(z) = (PH(x,y,z+1) + PHB(x,y,z+1)) / 9.806 - hgt(x,y)
      else
        delz(z) = ((PH(x,y,z+1) + PHB(x,y,z+1))-(PH(x,y,z) + PHB(x,y,z)))/9.806
      endif
    enddo
    if (delz(1) <= 0.) delz(1) = delz(2)
    !---Atmospheric Profile
    !- Level pressure is defined from 0 to zmax from the top of atmosphere to the surface
    atm(1)%Climatology         = TROPICAL
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID, O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS,VOLUME_MIXING_RATIO_UNITS /)
    atm(1)%Level_Pressure(0) = (pres(x,y,zmax)*3.0/2.0 - pres(x,y,zmax-1)/2.0)/100.0  ! convert from Pa to hPA %?? why
    do z=zmax,1,-1
      if(z.eq.1) then
        atm(1)%Level_Pressure(zmax-z+1) = psfc(x,y)/100.0    ! Pa -> hPa
        !max(psfc(x,y), pres(x,y,1)*3.0/2.0-pres(x,y,2)/2.0)/100.0
      else
        atm(1)%Level_Pressure(zmax-z+1) = ((pres(x,y,z-1)+pres(x,y,z))/2.0)/100.0  ! convert from Pa to hPA
      endif
      atm(1)%Pressure(zmax-z+1)       = pres(x,y,z) / 100.0
      atm(1)%Temperature(zmax-z+1)    = tk(x,y,z)
      atm(1)%Absorber(zmax-z+1,1)     = qvapor(x,y,z)*1000.0
    enddo
    atm(1)%Absorber(:,2) = 5.0E-02 
    !---Cloud Profile
    do z=1,N_CLOUDS
     atm(1)%Cloud(z)%Type = 0
     atm(1)%Cloud(z)%Effective_Radius = 0.0
     atm(1)%Cloud(z)%Water_Content = 0.0
    enddo
    ncl = 0
    icl = 0
    !--calculating # of clouds (cloud and rain)
    do z=zmax,1,-1
      if(qcloud(x,y,z).gt.0.0) then
        ncl = ncl + 1
      endif
      if(qrain(x,y,z).gt.0.0) then
        ncl = ncl + 1
      endif
      if(qice(x,y,z).gt.0.0) then
        ncl = ncl + 1
      endif
      if(qsnow(x,y,z).gt.0.0) then
        ncl = ncl + 1
      endif
      if(qgraup(x,y,z).gt.0.0) then
        ncl = ncl + 1
      endif
    enddo
    !--Data for cloud
    atm(1)%n_Clouds         = ncl
    IF ( atm(1)%n_Clouds > 0 ) THEN
    do z=zmax,1,-1
      if(qcloud(x,y,z).gt.0.0) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = WATER_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 16.8_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qcloud(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qrain(x,y,z).gt.0.0) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = RAIN_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1000.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qrain(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qice(x,y,z).gt.0.0) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = ICE_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 25.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qice(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qsnow(x,y,z).gt.0.0) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = SNOW_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 750.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qsnow(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
      endif
    enddo
    do z=zmax,1,-1
      if(qgraup(x,y,z).gt.0.0) then
        icl = icl + 1
        k1 = zmax-z+1
        k2 = zmax-z+1
        atm(1)%Cloud(icl)%Type = GRAUPEL_CLOUD
        atm(1)%Cloud(icl)%Effective_Radius(k1:k2) = 1500.0_fp
        atm(1)%Cloud(icl)%Water_Content(k1:k2)    = &
            qgraup(x,y,z)*pres(x,y,z)/287.2/(tk(x,y,z)+0.61*(qvapor(x,y,z)/(1+qvapor(x,y,z))))*delz(z)
      endif
    enddo
    ENDIF
  
    !---Surface data
    if(landmask(x,y).eq.1.0) then
     sfc(1)%Water_Coverage = 0.0_fp
     sfc(1)%Land_Coverage = 1.0_fp
     sfc(1)%Land_Temperature = tsk(x,y)
     sfc(1)%Soil_Temperature = tsk(x,y)
    else
     sfc(1)%Water_Coverage = 1.0_fp
     sfc(1)%Land_Coverage = 0.0_fp
     sfc(1)%Water_Type = 1  ! Sea water
     sfc(1)%Water_Temperature = tsk(x,y)
    endif
  
    ! --------------------------------------------------------------------
    ! 9-2. Input satellite viewing geometry
    ! --------------------------------------------------------------------
    !  The Sensor_Scan_Angle is optional.
    CALL CRTM_Geometry_SetValue( Geometry, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE )
  
  
    ! --------------------------------------------------------------------
    ! 9-3. Use the SOI radiative transfer algorithm
    ! --------------------------------------------------------------------
    Options%RT_Algorithm_ID = RT_SOI

    ! --------------------------------------------------------------------
    ! 9-4. Run CRTM forward model
    ! --------------------------------------------------------------------
    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution , &
                                 Options = Options )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP
    END IF


    ! --------------------------------------------------------------------
    ! 9-5. Collecting CRTM forward model output
    ! --------------------------------------------------------------------
    ! User should read the user guide or the source code of the routine
    ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
    ! select the needed variables for outputs.  These variables are contained
    ! in the structure RTSolution.

    !---for file output, edited 2014.9.26
    do l = 1, n_Channels
        BTsend(x,y,l) = real(RTSolution(l,1)%Brightness_Temperature)
        if(BTsend(x,y,l) /= BTsend(x,y,l) .or. BTsend(x,y,l)>HUGE(BTsend(x,y,l)) &
           .or. BTsend(x,y,l) < 100 .or. BTsend(x,y,l) > 400 ) then
          BTsend(x,y,l)=-888888.
        endif
    enddo
  
    !--- end of pt_id(x,y)-loop
  enddo

  ! -----------------------------------------------------------------------
  ! 10. Broadcast BT data across all processes via MPI_Allreduce
  ! ----------------------------------------------------------------------- 
  call MPI_Barrier( comm, ierr )
  if(my_proc_id==0)  write(*,*) "Finished running CRTM on all points"
  CALL MPI_Allreduce(BTsend,BT,xmax*ymax*n_Channels,MPI_REAL,MPI_SUM,comm,ierr)
  if(my_proc_id==0)  write(*,*) "CRTM solutions sent to all processes"

  if(my_proc_id==0)  write(*,*) "location 5"

  ! -----------------------------------------------------------------------
  ! 11. Writing output
  ! ----------------------------------------------------------------------- 
  if(my_proc_id==0) then
    open(10,file=out_fname,&
           form='unformatted',access='direct',recl=1)
    irec = 0
    do y = 1, ymax
    do x = 1, xmax
      irec= irec +1
      write( 10, rec=irec) xlong(x,y)
    enddo
    enddo
    if(my_proc_id==0)  write(*,*) "Wrote Longitudes"
    if(my_proc_id==0)  write(*,*) "Length of xlong:", irec

    do y = 1, ymax
    do x = 1, xmax
      irec= irec +1
      write( 10, rec=irec) xlat(x,y)
    enddo
    enddo
    if(my_proc_id==0)  write(*,*) "Wrote Latitudes"
    if(my_proc_id==0)  write(*,*) "Length of xlong and xlat:", irec

    do l = 1, n_ch
      do y = 1, ymax
      do x = 1, xmax
        irec= irec +1
        write( 10, rec=irec) BT(x,y,l)
      enddo
      enddo
    enddo
    if(my_proc_id==0)  write(*,*) "Wrote BT to files"
    if(my_proc_id==0)  write(*,*) "Total length of records:", irec

    close (10)
  endif

  ! -----------------------------------------------------------------------
  ! 12. Destroy CRTM and end program
  ! ----------------------------------------------------------------------- 
  if(my_proc_id==0) WRITE( *, *)"Destroying the CRTM..."
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP
  END IF
  ! ============================================================================

  call parallel_finish()

END PROGRAM XbtoIR

