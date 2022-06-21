PROGRAM test
    
    ! Disable all implicit typing
    IMPLICIT NONE
    
    !
    TYPE string
        CHARACTER(LEN=:), ALLOCATABLE :: s
    END TYPE string
    
    INTEGER :: ifile
    INTEGER :: num_file=12 
    CHARACTER(len=100) :: path = '/work2/06191/tg854905/stampede2/opt/CRTM/PSU_EnKF_CRTM/coefficients/'
    CHARACTER(len=150) :: spccoeff_file

    TYPE(string) :: spcfile(11)
    spcfile(1)%s = 'amsr2_gcom-w1.SpcCoeff.bin'
    spcfile(2)%s = 'atms_npp.SpcCoeff.bin'
    spcfile(3)%s = 'gmi_gpm.SpcCoeff.bin'
    spcfile(4)%s = 'mhs_metop-a.SpcCoeff.bin'
    spcfile(5)%s = 'mhs_metop-b.SpcCoeff.bin'
    spcfile(6)%s = 'mhs_n18.SpcCoeff.bin'
    spcfile(7)%s = 'mhs_n19.SpcCoeff.bin'
    spcfile(8)%s = 'saphir_meghat.SpcCoeff.bin'
    spcfile(9)%s = 'ssmi_f15.SpcCoeff.bin' 
    spcfile(10)%s = 'ssmis_f16.SpcCoeff.bin'
    spcfile(11)%s = 'ssmis_f17.SpcCoeff.bin'
    spcfile(12)%s = 'ssmis_f18.SpcCoeff.bin'


    DO ifile  = 1, num_file
        print *, 'Sensor: ', trim( spcfile(ifile)%s )
        spccoeff_file = trim( path ) // trim( spcfile(ifile)%s ) 
        CALL output_info(spccoeff_file, spcfile(ifile)%s)
    END DO
    !DEALLOCATE( spcfile)

END PROGRAM test

SUBROUTINE output_info(spccoeff_file, filename)

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
    CHARACTER(*), INTENT(IN) :: filename
    CHARACTER(len=2), DIMENSION(:), ALLOCATABLE :: PLR_char
    LOGICAL :: Quiet
    LOGICAL :: noisy
    TYPE(SpcCoeff_type) :: SC


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
                    PLR_char = '/'
            END SELECT Polarization_Type
        END DO 

    DO icount  = 1, n_Channels
        write (*,'(a, i2, a, f10.6,a)' ) 'Channel index: ',SC%Sensor_Channel(icount),'---',SC%Frequency(icount), PLR_char(icount)
    END DO
   
    END IF

   DEALLOCATE( PLR_char)

END SUBROUTINE




! Polarization interpretation: CRTM_SfcOptics.f90
! 5: Vertical linear polarization; 6: Horizontal linear polarization
