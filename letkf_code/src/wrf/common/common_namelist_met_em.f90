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
  INTEGER :: nbv_ori=20    ! Original number of ensemble members
  INTEGER :: nbv_tar=80    ! Target number of perturbed ensemble members
  INTEGER :: ntimes=1      ! Number of times to perturb
  INTEGER :: method=1      ! RejuGauss 1 , RejuUniform 2 , Specular 3 
  INTEGER :: niter=10      ! Number of sections in which we will divide the output ensemble
                           ! to reduce memory use. 
  LOGICAL :: recenter_mean=.FALSE. !If true the original ensemble mean will be recenter around 
                                   !a different state (e.g. a higher resolution deterministic run) 
  REAL(r_sngl) :: sigma = 1.0e-3   !Perturbation amplitude (methods 1 and 2 only)

  INTEGER :: iseed=10              !Random seed for the perturbation weigth generator.
                                   !When the pert_met_em program is called multiple times in the same experiment, 
                                   !fixing the random seed ensures that the combination of perturbations is always the same.

  CHARACTER(500) :: file_ini=' ' , file_end=' ' , file_tar=' '          !File names to interpolate.
                                                                        !File ini corresponds to the initial time, file_end to the
                                                                        !final time and file_tar to the target time.
  CHARACTER(19)  :: date_tar=' '                                           
  REAL(r_sngl)   :: time_ini=0.0e0 , time_end=0.0e0 , time_tar=0.0e0    !Times corresponding to the initial, final and target times
                                                                        !in any time continous unit (eg seconds)


  CHARACTER(LEN=50) :: NAMELIST_FILE='./pertmetem.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / GENERAL / nbv_ori , nbv_tar , ntimes , method , sigma , niter , recenter_mean , iseed 
NAMELIST / INTERP / time_ini , time_end , time_tar , date_tar , file_ini , file_end , file_tar

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
REWIND(54)
READ(54,NML=INTERP,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at INTERP section"
WRITE(*,*)"Using default values"
ENDIF
CLOSE(54)

END SUBROUTINE READ_NAMELIST

END MODULE common_namelist_met_em
