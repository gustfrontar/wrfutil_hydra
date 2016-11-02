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
  IMPLICIT NONE
  PUBLIC

!-----------------------------------------------------------------------------
! RANDOM PERTURBATION PARAMETERS 
!-----------------------------------------------------------------------------

 !Wheter variables will be randomly perturbed
 LOGICAL :: PERTURB_T = .TRUE.
 LOGICAL :: PERTURB_RH = .TRUE.
 LOGICAL :: PERTURB_WIND = .TRUE.

 !Random perturbation amplitude.
 REAL(r_size) :: PERTURB_T_AMP = 0.5d0
 REAL(r_size) :: PERTURB_RH_AMP = 5.0d0
 REAL(r_size) :: PERTURB_WIND_AMP = 0.5d0

 REAL(r_size) :: RANDOM_AMP_FACTOR=0.5d0

 !Random perturbation lenght scale (horizontal)
 REAL(r_size) :: PERTURB_T_SCLH = 40000d0
 REAL(r_size) :: PERTURB_RH_SCLH = 40000d0
 REAL(r_size) :: PERTURB_WIND_SCLH = 40000d0

 !Random perturbation length scale (vertical)
 REAL(r_size) :: PERTURB_T_SCLV = 5000d0
 REAL(r_size) :: PERTURB_RH_SCLV = 5000d0
 REAL(r_size) :: PERTURB_WIND_SCLV = 5000d0

!-----------------------------------------------------------------------------
! BALANCED PERTURBATION PARAMETERS 
!-----------------------------------------------------------------------------

 LOGICAL :: PERTURB_ATMOSPHERE= .TRUE.  !If atmospheric variables will be perturbed
 LOGICAL :: PERTURB_SST       = .TRUE.  !If sst will be perturbed
 LOGICAL :: PERTURB_SOIL      = .TRUE.  !If soiltemp and soilmoisture will be perturbed

 REAL(r_size) :: AMP_FACTOR !Balanced perturbation amplitude

!-----------------------------------------------------------------------------
!  GENERAL PARAMETERS
!-----------------------------------------------------------------------------

 INTEGER       :: current_time=1 , max_time=1

 CHARACTER(LEN=20) :: NAMELIST_FILE='./pertmetem.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / BALANCED    / current_time , max_time , amp_factor,                  &
                         perturb_atmosphere , perturb_sst , perturb_soil 

NAMELIST / RANDOM      / perturb_t , perturb_rh , perturb_wind,                 &
                         perturb_t_amp  , perturb_rh_amp  , perturb_wind_amp,   &
                         perturb_t_sclh , perturb_rh_sclh , perturb_wind_sclh,  &
                         perturb_t_sclv , perturb_rh_sclv , perturb_wind_sclv,  &
                         random_amp_factor


INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=BALANCED,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at BALANCED section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=RANDOM,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at RANDOM section"
WRITE(*,*)"Using default values"
ENDIF

CLOSE(54)


END SUBROUTINE READ_NAMELIST

END MODULE common_namelist
