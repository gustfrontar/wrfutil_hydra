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

  !--------------------------------------------FAKE RADAR PARAMETERS
  LOGICAL   ::   FAKE_RADAR=.TRUE.           !We will create a non existant radar.
  LOGICAL   ::   COMPUTE_ATTENUATION=.FALSE. !If attenuation will be simulated or not when converting from model to radar.
  INTEGER :: n_radar =1 !Number of radars.
  INTEGER :: n_model =1 !Number of model files that will be interpolated to the radar.
  REAL(r_size) :: fradar_lon(max_radars)=undef , fradar_lat(max_radars)=undef !Grid point position of radar.
  REAL(r_size) :: fradar_z(max_radars)=0.0d0
  REAL(r_size) :: radar_az_res(max_radars) =1.0d0 , radar_r_res(max_radars) =500.0d0 , radar_el_res(max_radars)=1.25d0
  REAL(r_size) :: radar_min_az(max_radars) =0.0d0 , radar_max_az(max_radars)=360.0d0
  REAL(r_size) :: radar_min_r(max_radars)  =0.0d0 , radar_max_r(max_radars) =240000
  REAL(r_size) :: radar_min_el(max_radars) =0.1d0 , radar_max_el(max_radars)=18.0d0
  LOGICAL   :: use_level_list=.false.  !If true then we can specifiy a list of vertical lelves. (this will be applied to all radars)
  LOGICAL   :: use_wt                  !If terminal velocity will be used in radial velocity computations.
  INTEGER   :: n_vertical_levels=1     !Number of vertical levels (in case use_level_list == true )
  REAL(r_size) :: level_list(max_vert_levels)=-9.99d0  !List of vertical levels.

  !Define the type of forward operator.
  REAL(r_size) :: radar_lambda(max_radars)=3.0d0  !Wave length of the radar (to be used in the computation of reflectivity.
                                                  !3 Xband , 5 Cband , 10 Sband
  CHARACTER(1) :: endian = 'b'                    !'b' -> big_endian , 'l' -> little_endian
  CHARACTER(100) :: w2rnamelist='w2r.namelist'

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE get_namelist_vars
IMPLICIT NONE
LOGICAL :: file_exist
INTEGER :: ierr


NAMELIST / GENERAL / add_obs_error , reflectivity_error , radialwind_error , &
                     fake_radar , fradar_lon , fradar_lat , fradar_z , &
                     radar_az_res  , radar_r_res , radar_el_res      , &
                     radar_min_az  , radar_max_az   , &
                     radar_min_el  , radar_max_el   , &
                     radar_min_r   , radar_max_r    , &
                     use_wt , radar_lambda          , &
                     use_level_list , n_vertical_levels , level_list  , &
                     n_radar , n_model                                  !Number of radars that will be processed.

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
