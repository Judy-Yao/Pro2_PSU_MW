PROGRAM test

    USE Message_Handler
    USE SpcCoeff_Define
    USE SpcCoeff_Binary_IO, ONLY:SpcCoeff_Binary_InquireFile,SpcCoeff_Binary_ReadFile
    USE Binary_File_Utility, ONLY: Open_Binary_File

    ! Disable all implicit typing
    IMPLICIT NONE

    ! Function result
    INTEGER :: err_stat
    INTEGER :: FileID
    INTEGER :: n_Channels
    INTEGER :: ERR  
    INTEGER :: icount
    CHARACTER(len=120) :: spccoeff_file
    CHARACTER(len=30) :: msg
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: PLR_char
    CHARACTER(len=*), PARAMETER :: Vpol = 'Vpol'
    CHARACTER(len=*), PARAMETER :: Hpol = 'Hpol'
    LOGICAL :: Quiet
    LOGICAL :: noisy
    TYPE(SpcCoeff_type) :: SC


    spccoeff_file = "/work2/06191/tg854905/stampede2/opt/CRTM/PSU_EnKF_CRTM/coefficients/amsr2_gcom-w1.SpcCoeff.bin"
    err_stat = SpcCoeff_Binary_InquireFile(spccoeff_file,n_Channels) 
    IF ( err_stat == SUCCESS ) THEN
        print *, 'Number of Channels:'
        print *, n_Channels
    END IF
    err_stat = SpcCoeff_Binary_ReadFile( &
        spccoeff_file      , &
        SC              , &
        Quiet = .NOT. noisy  )


    IF  ( err_stat == SUCCESS ) THEN
        ALLOCATE( PLR_char(n_Channels), STAT= ERR )
        DO icount  = 1, n_Channels
            Polarization_Type: SELECT CASE( SC%Polarization(icount) ) 
                CASE( 5 )
                    PLR_char(icount) = 'V'
                CASE( 6 )
                    PLR_char(icount) = 'H'
                CASE DEFAULT
                    print *, 'Unrecognised polarization type!'
            END SELECT Polarization_Type
        END DO 

    DO icount  = 1, n_Channels
        write (*,'(a, i2, a, f7.3,a)' ) 'Channel index: ',SC%Sensor_Channel(icount),'---',SC%Frequency(icount), PLR_char(icount)
    
    END DO
    !print *, SC%Release, SC%Version
    !print *, 'Channel index (used in CRTM library): '
    !print *, SC%Sensor_Channel
    
    !print *, 'Polarization type:'
    !print *, PLR_char
    !print *, 'Red flag for channels:'
    !print *, SC%Channel_Flag
    !print *, 'Channel frequency: '
    !print *, SC%Frequency
   END IF

   DEALLOCATE( PLR_char)
END PROGRAM



! Polarization interpretation: CRTM_SfcOptics.f90
! 5: Vertical linear polarization; 6: Horizontal linear polarization
