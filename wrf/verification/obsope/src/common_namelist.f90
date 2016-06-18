MODULE common_namelist
!=======================================================================
!
! [PURPOSE:] Namelist input for LETKF analysis
!
! [HISTORY:]
!   25/09/2009 Juan Ruiz  created
!
!=======================================================================
  USE common
  USE common_wrf
  IMPLICIT NONE
  PUBLIC

  !GENERAL
  INTEGER :: nslots=2 ! number of time slots for 4D-LETKF
  INTEGER :: nbslot=2 ! basetime slot
  INTEGER :: nbv=40   ! Number of ensemble members
 
  !OBSERVATIONS
  REAL(r_size) :: threshold_dz=1000.0d0
  REAL(r_size) :: gross_error=15.0d0
  REAL(r_size) :: gross_error_tycll=10.0d0        !degree
  REAL(r_size) :: gross_error_tycmip=60.0d2       !Pa
  REAL(r_size) :: rainratio_threshold = 0.3d0     !Percentaje of ensemble members with reflectivity greather than 0.
  REAL(r_size) :: minrefdbz=0.0d0                 !Reflectivity values below this threshold won't be assimilated.
  REAL(r_size) :: pseudo_rh_error=0.1             !Obserational error for pseudo RH observations.

  !RADAR DA
  INTEGER            :: interpolation_technique=1
  INTEGER            :: nradar = 1                !Number of assimilated radars
  LOGICAL            :: use_wt=.true.             !Use or not terminal velocity in the computation of VR

  !OBSGRID SECTION
  LOGICAL            :: do_obsgrid=.false.        !Wheater obsgrid output will be generated.
  REAL(r_size)       :: regrid_res=1.0d0          !Horizontal resolution of obsgrid grid (degree)
  REAL(r_size)       :: regrid_vert_res=10000.0d0 !Vertical resolution of obsgrid grid (hPa)
  REAL(r_size)       :: regrid_vert_zres=1000.0d0 !Vertical resolution of obsgrid grid (m)
  INTEGER            :: nlevreg=10                !Number of vertical levels.
  INTEGER, PARAMETER :: maxarea=20                !Maximum number of verification areas
  INTEGER            :: narea=1
  REAL(r_size)       :: vlon1(maxarea)=0.0d0
  REAL(r_size)       :: vlon2(maxarea)=360.0d0
  REAL(r_size)       :: vlat1(maxarea)=-90.0d0
  REAL(r_size)       :: vlat2(maxarea)=90.0d0
  LOGICAL            :: regrid_output=.false.
  INTEGER            :: slotstep=1                !Verification times can be gropued using this parameter.
  INTEGER            :: slotoffset=0              !In which time slot will the verification start.


  !Verification variables flags
  !3D variables
  integer , parameter :: maxverifvars=20
  integer :: variable_list(maxverifvars)     !Use observation ID to define variables
  CHARACTER(LEN=20) :: NAMELIST_FILE='./obsope.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / GENERAL / nslots , nbslot , nbv
NAMELIST / OBSERVATIONS / threshold_dz , gross_error , gross_error_tycll , gross_error_tycmip , &
&   rainratio_threshold , minrefdbz   
NAMELIST / RADAR_DA  / interpolation_technique ,  nradar , use_wt
NAMELIST / OBSGRID /  narea , vlon1 , vlon2 , vlat1 , vlat2 ,   regrid_output, regrid_res , regrid_vert_res , &
                     regrid_vert_zres , slotstep , slotoffset , do_obsgrid , variable_list , nlevreg


INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=OBSERVATIONS,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at OBSERVATIONS section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=RADAR_DA,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at RADAR_DA section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=OBSGRID,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at OBSGRID section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)


CLOSE(54)

END SUBROUTINE READ_NAMELIST

END MODULE common_namelist
