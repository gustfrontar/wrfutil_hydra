MODULE common_namelist
!=======================================================================
!
! [PURPOSE:] Observational procedures
!
! [HISTORY:]
!   01/23/2009 Takemasa MIYOSHI  created
!
!=======================================================================
  USE common

  IMPLICIT NONE
  PUBLIC

  !Namelist variables
  LOGICAL            :: ADD_OBS_ERROR=.true. !Flag to include or not an unbiased gaussian perturbation in the observation.

  !Magnitude of observational error for different observations.
  REAL(r_size)       :: U_ERROR=1.0d0 !Standard deviation for observation errors.
  REAL(r_size)       :: V_ERROR=1.0d0   
  REAL(r_size)       :: T_ERROR=1.0d0
  REAL(r_size)       :: Q_ERROR=1.0d-3
  REAL(r_size)       :: PS_ERROR=1.0d2
  !NOTE: Even if ADD_OBS_ERROR is false, the error information will be included in the observation file and the values.
  !listed above will be indicated as the observation error.

  !Flags to include or not different observation types.
  LOGICAL    :: UOBS  = .true.
  LOGICAL    :: VOBS  = .true.
  LOGICAL    :: TOBS  = .true.
  LOGICAL    :: QOBS  = .true.
  LOGICAL    :: PSOBS = .true.

  INTEGER      :: OBSTYPE=1        !Observation spatial distribution type. 1 is random, 2 is regular and 3 is read from input file.
  REAL(r_size) :: OBSDENSITY=1.0d0 !For random obs only. Observation density 1 means one observation per grid point.
  INTEGER      :: SKIP=1           !For regular obs only. How many horizontal grid points to skip between observations.
  INTEGER      :: SKIPZ=1          !For regular obs only. How many vertical grid points to skip between observations.

  LOGICAL      :: use_wt=.true.    !This variable has to be here to keep compatibility with letkf-wrf modules.

  CHARACTER(1) :: endian = 'b'                    !'b' -> big_endian , 'l' -> little_endian
  CHARACTER(100) :: w2bnamelist='w2b.namelist'

CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE get_namelist_vars
IMPLICIT NONE
LOGICAL :: file_exist
INTEGER :: ierr


NAMELIST / GENERAL / add_obs_error , u_error,v_error,t_error,q_error,ps_error , &
                     uobs,vobs,tobs,qobs,psobs,     &
                     obstype,obsdensity,skip,skipz, &
                     endian

INQUIRE(FILE=w2bnamelist,EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=w2bnamelist)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

END SUBROUTINE get_namelist_vars

END MODULE common_namelist
