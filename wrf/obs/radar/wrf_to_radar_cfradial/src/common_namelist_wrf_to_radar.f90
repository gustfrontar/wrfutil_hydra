MODULE common_namelist
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
!$USE OMP_LIB
  USE common

  IMPLICIT NONE
  PUBLIC

  INTEGER, PARAMETER :: max_vert_levels=100
  INTEGER, PARAMETER :: max_radars=100

  !Namelist variables
  LOGICAL            :: ADD_OBS_ERROR=.true.
  REAL(r_size)       :: REFLECTIVITY_ERROR=1.0d0 !Standard deviation of reflectivity error in DBZ
  REAL(r_size)       :: RADIALWIND_ERROR=0.5d0   !Standard deviation of radialwind error in m/s
  REAL(r_size)       :: min_ref = 1.0d0          !Minimum reflectivity value (in Power units)
  INTEGER            :: METHOD_REF_CALC=2        !Method used for reflectivity computaton.

  !--------------------------------------------FAKE RADAR PARAMETERS
  LOGICAL   ::   COMPUTE_ATTENUATION=.FALSE. !If attenuation will be simulated or not when converting from model to radar.
  INTEGER :: n_radar =1 !Number of radars.
  INTEGER :: n_times =1 !Number of model files that will be interpolated to the radar.
  LOGICAL   :: use_wt                  !If terminal velocity will be used in radial velocity computations.

  LOGICAL :: USE_PSEUDORH = .FALSE.    !This flag is included to make this code compatible with LETKF_WRF observation operator.

  CHARACTER(100) :: w2rnamelist='w2r.namelist'

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE get_namelist_vars
IMPLICIT NONE
LOGICAL :: file_exist
INTEGER :: ierr

NAMELIST / GENERAL / add_obs_error , reflectivity_error , radialwind_error ,  &
                     use_wt , n_radar , n_times , method_ref_calc

INQUIRE(FILE=w2rnamelist,EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=w2rnamelist)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

END SUBROUTINE get_namelist_vars

END MODULE common_namelist
