PROGRAM MAIN_RADAR_QC
!=======================================================================
!
! [PURPOSE:] To perform radar quality control.
!
! [HISTORY:]
!   02/09/2014 Juan Ruiz created
! TODO
! -Put all the options into a namelist file.
! -Consider other radar sources and other variables (i.e. polarimetric)! 
! -Include a flexible treatment of the topography fields 
!  direct interpolation from tiff files.
! -Detection of changes in radar strategy (time or vertical resolution) 
!=======================================================================

USE common
USE common_radar_tools
USE common_qc_tools

 IMPLICIT NONE

! INTEGER     ,ALLOCATABLE :: radarqc(:,:,:)
 REAL(r_size)             :: rtimer00 , rtimer
 INTEGER                  :: i

 INTEGER, PARAMETER   ::  NTIMES = 10  !Number of times that will be read.

 INTEGER(8), PARAMETER   ::  INIT_DATE = 20130713150010     !20130713080010 !Initial time yyyymmddhhMMSS
 INTEGER(8), PARAMETER   ::  END_DATE  = 20130713200010 !End time yyyymmddhhMMSS
 INTEGER(8)              ::  C_DATE

 CHARACTER(200),PARAMETER ::  DATA_PATH='./'
 CHARACTER(200),PARAMETER ::  OUT_PATH ='./'
 CHARACTER(21)            ::  RADAR_FILE_REF ='Z_yyyymmdd-HHMMSS.dat' 
 CHARACTER(22)            ::  RADAR_FILE_WIND='VR_yyyymmdd-HHMMSS.dat'
 CHARACTER(38)            ::  RADAR_FILE_OUT ='PAWR_LETKF_INPUTV3_yyyymmdd-HHMMSS.dat'
 CHARACTER(300)           ::  INPUT_FILE
 CHARACTER(300)           ::  OUTPUT_FILE

 INTEGER, PARAMETER   ::  DT        = 30  !Time interval in seconds
 INTEGER              ::  ITIME

 INTEGER              :: CYYYY,CMM,CDD,CHH,CMN,CSS, NX, NY, NZ

 REAL(r_size),ALLOCATABLE :: reflectivity(:,:,:,:) !Azimuth, range, elevation and time.
 REAL(r_size),ALLOCATABLE :: wind(:,:,:,:)         !Azimuth, range, elevation and time.

 REAL(r_size),ALLOCATABLE :: qcedreflectivity(:,:,:) !Azimuth, range, elevation and time.
 REAL(r_size),ALLOCATABLE :: qcedwind(:,:,:)         !Azimuth, range, elevation and time.

 LOGICAL  :: file_exist
 LOGICAL  :: do_qc


 REAL(r_size) :: tmp
 INTEGER      :: ia , ir , ie
! !------------------------------------------------------------------
! !  CONFIGURATION PARAMETERS

  USE_PRIOR=.TRUE.
!
  !--------GENERAL PARAMETERS ------------------------
  DEBUG_OUTPUT= .FALSE.     !Turns on individual parameter output.
  !Generate a sample that will be used for verification.
  STATISTIC_OUTPUT=.FALSE.  !Turns on output for statistical analysis of the parameters.
  SKIP_AZ=1
  SKIP_R=2
  SKIP_ELEV=2
  ELEVMAX=20

  !VERTICAL COLUMN PARAMETERS 
  USE_ECHO_TOP=.FALSE. 
  USE_ECHO_BASE=.FALSE.
  USE_ECHO_DEPTH=.FALSE. 
  USE_MAX_DBZ=.FALSE. 
  USE_MAX_DBZ_Z=.FALSE. 
  USE_REF_VGRAD=.TRUE.
  !-------- LOCAL MINIMUM OF LOCAL VARIANCE OF DBZ  --------------------
  USE_LDBZVAR = .FALSE.
  !-------- LOCAL "ISOLATION" OF THE REFLECTIVITY  --------------------
  USE_SPECKLE = .TRUE.
  !-------- REMOVE BLOCKED BEANS --------------------
  USE_BLOCKING = .TRUE.
  !-------- TEXTURE OF REFLECITIVITY (TDBZ)  --------------------
  USE_TDBZ = .FALSE.
  !-------- SIGN OF REFLECITIVITY (SIGN)  --------------------
  USE_SIGN = .FALSE.
  !-------- SPIN OF REFLECITIVITY (SPIN)  --------------------
  USE_RSPIN = .FALSE.
  !-------- TEMPORAL SIGN OF REFLECITIVITY (TSIGN)  --------------------
  USE_TSIGN = .FALSE.
  !-------- TEMPORAL WIND AVERAGE  --------------------
  USE_TWA = .TRUE.
  !-------- TEMPORAL WIND STD  --------------------
  USE_TWS = .FALSE.
  !-------- TEMPORAL WIND CORRELATION --------------------
  USE_TWC = .FALSE.
  !-------- TEMPORAL REFLECITIVITY CORRELATION --------------------
  USE_TRC = .FALSE.
  !-------- TEMPORAL REFLECITIVITY CORRELATION TEXTURE--------------------
  USE_TRCT = .TRUE.
  !-------- TEMPORAL REF MAX TEMPORAL STD --------------------
  USE_TRMS = .FALSE.
  !-------- TEMPORAL REF ANOMALY MSD --------------------
  USE_TRAM = .FALSE.
  !-------- SPATIAL REF ANOMALY MSD --------------------
  USE_SRAM = .FALSE.
  !-------- TEMPORAL WIND ANOMALY MSD --------------------
  USE_TWAM = .FALSE.
  !-------- SPATIAL WIND ANOMALY MSD --------------------
  USE_SWAM = .FALSE.

  ! 1 means Naive Bayes, 2 means Fuzzy logic.
  WEIGHT_TYPE=1

  !-------- OTHER FILTERS --------------------
  USE_ECHO_TOP_FILTER=.TRUE.
  ECHO_TOP_FILTER_THRESHOLD=3000
  USE_ECHO_DEPTH_FILTER=.TRUE.
  ECHO_DEPTH_FILTER_THRESHOLD=2000

 CALL CPU_TIME(rtimer00)

 ITIME=1

 !Main time loop 
 C_DATE = INIT_DATE

 DO WHILE ( C_DATE <= END_DATE )
 do_qc=.true.

  WRITE(6,*)'CURRENT DATE IS ',C_DATE

  CALL com_time2ymdhms(C_DATE,CYYYY,CMM,CDD,CHH,CMN,CSS)

  !Generate file names.
  WRITE(RADAR_FILE_REF(3:10),'(I4.4,I2.2,I2.2)')CYYYY,CMM,CDD
  WRITE(RADAR_FILE_REF(12:17),'(I2.2,I2.2,I2.2)')CHH,CMN,CSS

  WRITE(RADAR_FILE_WIND(4:11),'(I4.4,I2.2,I2.2)')CYYYY,CMM,CDD
  WRITE(RADAR_FILE_WIND(13:18),'(I2.2,I2.2,I2.2)')CHH,CMN,CSS
  
  WRITE(RADAR_FILE_OUT(20:27),'(I4.4,I2.2,I2.2)')CYYYY,CMM,CDD
  WRITE(RADAR_FILE_OUT(29:34),'(I2.2,I2.2,I2.2)')CHH,CMN,CSS

  !ADVANCE TIME
  CALL com_timeinc_sec(C_DATE,DT)

  !READ THE DATA
  INPUT_FILE= TRIM(DATA_PATH) // RADAR_FILE_REF

  

   WRITE(6,*)'WILL READ A FILE ',RADAR_FILE_REF
   INQUIRE( FILE=INPUT_FILE, EXIST=file_exist )
   IF( file_exist )THEN
    CALL radar_read_data_old( RADAR1 , INPUT_FILE )
    radar1%year  =REAL(CYYYY,r_size)
    radar1%month =REAL(CMM,r_size)
    radar1%day   =REAL(CDD,r_size)
    radar1%hour  =REAL(CHH,r_size)
    radar1%minute=REAL(CMN,r_size)
    radar1%second=REAL(CSS,r_size)
   ELSE
    WRITE(*,*)"[Warining]: File ",RADAR_FILE_REF, "does not exist"
    do_qc=.false.
    IF( ITIME == 1 )THEN
     WRITE(*,*)"[Error]: First file needs to be present"
     STOP
    ENDIF
   ENDIF

  IF( ITIME == 1 )THEN
    ALLOCATE(reflectivity(RADAR1%na,RADAR1%nr,RADAR1%ne,NTIMES))
    reflectivity=UNDEF
    ALLOCATE(qcedreflectivity(RADAR1%na,RADAR1%nr,RADAR1%ne))
    qcedreflectivity=UNDEF
  ENDIF

  !Change low reflectivity values to undef (qc will ignore those)
  WHERE( RADAR1%radarv3d(:,:,:,1) <= minz )
    RADAR1%radarv3d(:,:,:,1)=undef
  ENDWHERE

  IF( do_qc )THEN
  reflectivity(:,:,:,NTIMES)=RADAR1%radarv3d(:,:,:,1)
  ENDIF

  WRITE(6,*)'WILL READ A FILE ',RADAR_FILE_WIND
  INPUT_FILE= TRIM(DATA_PATH) // RADAR_FILE_WIND
  INQUIRE( FILE=INPUT_FILE, EXIST=file_exist )
   IF( file_exist )THEN
    CALL radar_read_data_old( RADAR1 , INPUT_FILE )
   ELSE
    WRITE(*,*)"[Warining]: File ",RADAR_FILE_REF, "does not exist"
    do_qc=.false.
    IF( ITIME == 1 )THEN
     WRITE(*,*)"[Error]: First file needs to be present"
     STOP
    ENDIF
   ENDIF

  !Change undef wind values to our undef 
  WHERE( RADAR1%radarv3d(:,:,:,1) <= -299.0d0 )
    RADAR1%radarv3d(:,:,:,1)=undef
  ENDWHERE


  IF( ITIME == 1 )THEN
    ALLOCATE(wind(RADAR1%na,RADAR1%nr,RADAR1%ne,NTIMES))
    wind=UNDEF
    ALLOCATE(qcedwind(RADAR1%na,RADAR1%nr,RADAR1%ne))
    qcedwind=UNDEF
  ENDIF

  IF( do_qc )THEN
  wind(:,:,:,NTIMES)        =RADAR1%radarv3d(:,:,:,1)
  ENDIF


  IF( ITIME == 1 )THEN !Georeference radar data.
    WRITE(6,*)'GEOREFERENCING RADAR DATA'

    CALL radar_georeference( RADAR1 )
     NX= RADAR1%na
     NY= RADAR1%nr
     NZ= RADAR1%ne

  ENDIF

  !Perform QC
  IF( ITIME >= NTIMES .AND. do_qc )THEN  !We need at least ntimes volumes to perform QC.
  WRITE(6,*)'PERFORMING RADAR QC'
    OUTPUT_FILE= TRIM(OUT_PATH) // RADAR_FILE_OUT

    WRITE(6,*)'WILL WRITE THE RESULT TO: ',OUTPUT_FILE

    TIMEOUT=C_DATE
    CALL RADAR_QC(reflectivity,wind,NTIMES,qcedreflectivity,qcedwind)  

    WRITE(6,*)"SUCCESSFULY PERFORM QC FOR ",CYYYY,CMM,CDD,CHH,CMN,CSS

    !Write QCED DATA
    CALL RADAR_WRITE_FILE( RADAR1 , qcedreflectivity(:,:,:) , qcedwind(:,:,:) , qcflag , attenuation , OUTPUT_FILE )
  ELSE
   WRITE(6,*)"QC NOT PERFORMED FOR ",CYYYY,CMM,CDD,CHH,CMN,CSS

  ENDIF

  !Shift the data in time.
  IF( do_qc )THEN
    reflectivity(:,:,:,1:NTIMES-1)=reflectivity(:,:,:,2:NTIMES) 
    wind(:,:,:,1:NTIMES-1)        =wind(:,:,:,2:NTIMES)
  ENDIF

  ITIME=ITIME+1
 ENDDO



END PROGRAM MAIN_RADAR_QC
