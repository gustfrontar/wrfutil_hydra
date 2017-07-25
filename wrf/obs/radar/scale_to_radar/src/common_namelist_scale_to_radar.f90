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

  INTEGER, PARAMETER :: max_radars=100

  !Namelist variables
  LOGICAL            :: ADD_OBS_ERROR=.false.
  REAL(r_size)       :: REFLECTIVITY_ERROR=1.0d0 !Standard deviation of reflectivity error in DBZ
  REAL(r_size)       :: RADIALWIND_ERROR=0.5d0   !Standard deviation of radialwind error in m/s
  INTEGER            :: METHOD_REF_CALC=2        !Method used for reflectivity computaton.
  REAL(r_size)       :: min_ref = 1.0d0          !Minimum reflectivity value (in dBz)

  !--------------------------------------------FAKE RADAR PARAMETERS
  LOGICAL   :: compute_attenuation =.false. !If attenuation will be simulated or not when converting from model to radar.
  INTEGER   :: n_radar =1 !Number of radars.
  INTEGER   :: n_times =1 !Number of times that will be processed.
  LOGICAL   :: use_wt                  !If terminal velocity will be used in radial velocity computations.
  INTEGER   :: input_type=1  !1-restart files, 2-history files.

  !Projection parameters
  character(len=4), public :: MPRJ_type = 'NONE' !< map projection type
                                               ! 'NONE'
                                               ! 'LC'
                                               ! 'PS'
                                               ! 'MER'
                                               ! 'EC'
  real(r_size), public :: MPRJ_basepoint_x        ! position of base point in the model [m]
  real(r_size), public :: MPRJ_basepoint_y        ! position of base point in the model [m]
  real(r_size), public :: MPRJ_rotation =  0.0d0 ! rotation factor (only for 'NONE' type)
  real(r_size), public :: MPRJ_LC_lat1  = 30.0d0 ! standard latitude1 for L.C. projection [deg]
  real(r_size), public :: MPRJ_LC_lat2  = 60.0d0 ! standard latitude2 for L.C. projection [deg]
  real(r_size), public :: MPRJ_PS_lat             ! standard latitude1 for P.S. projection [deg]
  real(r_size), public :: MPRJ_M_lat    =  0.0d0 ! standard latitude1 for Mer. projection [deg] 
  real(r_size), public :: MPRJ_EC_lat   =  0.0d0 ! standard latitude1 for E.C. projection [deg]

  real(r_size), public :: MPRJ_basepoint_lon = 135.221d0 ! position of base point
  real(r_size), public :: MPRJ_basepoint_lat =  34.653d0 ! position of base point

  !Index parameters
  INTEGER  , public :: IHALO = 2 , JHALO = 2 , KHALO = 2


  !Define the type of forward operator.
  CHARACTER(100) :: s2rnamelist='s2r.namelist'

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE get_namelist_vars
IMPLICIT NONE
LOGICAL :: file_exist
INTEGER :: ierr


NAMELIST / GENERAL / add_obs_error , reflectivity_error , radialwind_error ,              &
                     use_wt , n_radar , n_times , method_ref_calc , input_type  


!Scale requires projection parameters
NAMELIST / PARAM_MAPPROJ / &
       MPRJ_basepoint_lon, &
       MPRJ_basepoint_lat, &
       MPRJ_basepoint_x,   &
       MPRJ_basepoint_y,   &
       MPRJ_type,          &
       MPRJ_rotation,      &
       MPRJ_LC_lat1,       &
       MPRJ_LC_lat2,       &
       MPRJ_PS_lat,        &
       MPRJ_M_lat,         &
       MPRJ_EC_lat

NAMELIST / PARAM_INDEX /   &
       IHALO ,             &
       JHALO ,             &
       KHALO 


INQUIRE(FILE=s2rnamelist,EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=s2rnamelist)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

READ(54,NML=PARAM_MAPPROJ,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at PARAM_MAPPROJ section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

READ(54,NML=PARAM_INDEX,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at PARAM_MAPPROJ section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)


END SUBROUTINE get_namelist_vars

END MODULE common_namelist
