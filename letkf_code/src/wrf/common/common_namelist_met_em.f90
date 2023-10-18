MODULE common_namelist_met_em
!=======================================================================
!
! [PURPOSE:] Namelist input for LETKF analysis
!
! [HISTORY:]
!   25/09/2009 Juan Ruiz  created
!
!=======================================================================
  USE common
  USE common_met_em
  IMPLICIT NONE
  PUBLIC

  !GENERAL
  INTEGER :: nbv_ori=20   ! Original number of ensemble members
  INTEGER :: nbv_tar=80   ! Target number of perturbed ensemble members
  INTEGER :: ntimes=1     ! Number of times to perturb
  INTEGER :: method=1     ! RejuGauss 1 , RejuUniform 2 , Specular 3 
  REAL(r_singl) :: sigma = 1.0e-3  !Perturbation amplitude (methods 1 and 2 only)

  CHARACTER(LEN=50) :: NAMELIST_FILE='./pertmetem.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / GENERAL / nbv_ori , nbv_tar , ntimes , method , sigma

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
  STOP 1
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
CLOSE(54)

END SUBROUTINE READ_NAMELIST

END MODULE common_namelist_met_em
