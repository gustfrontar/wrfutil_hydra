MODULE common_namelist
!=======================================================================
!
! [PURPOSE:] Namelist input for breeding rescaling module
!
! [HISTORY:]
!   25/09/2009 Juan Ruiz  created
!
!=======================================================================
  USE common
  IMPLICIT NONE
  PUBLIC

  !BREEDING SECTION

  character(10) :: norm_type = 'UVT'    !UV , T , UVT
  character(10) :: breeding_type = 'double'       !single (only pp and ctrl are used), double (uses pp, pn and ctrl to center the perturbation)
  character(10) :: input_file_prefix_ctrl = 'ctrl'
  character(10) :: input_file_prefix_pp   = 'pp'
  character(10) :: input_file_prefix_pn   = 'pn'

  character(10) :: input_file_prefix_pp_old   = 'pp_old' !Input file prefix for positive perturbation with breeding_in_place
  character(10) :: input_file_prefix_pn_old   = 'pn_old' !Input file prefix for negative perturbation with breeding_in_place

  logical :: breeding_in_place = .false.  !If true the current perturbation will be also stored in the perturbed restart files
                                          !corresponding to the beginning of the window.


  real(r_size) :: wind_scale = 1.0d0 , t_scale = 1.0d0 , p_scale = 1.0d0   !Scalling factor to combine multiple variables.
  real(r_size) :: norm_threshold = 1.0d0 !Norm threshold for rescalling.
  logical :: ave_norm_xy = .false. , ave_norm_z = .false.  !Averaging flags for horizontal and vertical averaging.

  !Average over a sub domain. Default equal 0 means average over the entire domain.
  integer :: is = 0  , ie = 0 , js = 0 , je = 0 , ks = 0 , ke = 0

  !Smoothing flags. Aplies a Lanczos smoothing. Smoothing in the vertical asumes dz = constant and equal to the mean
  !vertical resolution.
  
  logical :: smooth_norm_xy = .false. , smooth_norm_z = .false.

  !Horizontal and vertical smoothing scales. Variability at scales shorter than these will be removed by the filter.
  real( r_size) :: smooth_sx = 3.0d3 , smooth_sz = 3.0d3

  logical :: init_bv = .false.    !Initialize the BV as a random smoothed perturbation.
  logical :: perturb_t = .true. , perturb_wind = .true. !Enable or disable initial perturbations in temperature and wind.
  
  real(r_size) :: perturb_t_amp = 0.5 , perturb_wind_amp = 0.25 !Initial temperature and wind perturbation amplitude.
  real(r_size) :: perturb_t_sclh = 0.5d3 , perturb_t_sclv = 0.5d3  !Horizontal scales for initial perturbations.
  real(r_size) :: perturb_wind_sclh = 0.5d3 , perturb_wind_sclv = 0.5d3  !Vertical scale for initial perturbations.

  integer :: bdy_width = 10 !Region in which perturbation amplitude will be reduced near the boundaries.

 CHARACTER(LEN=20) :: NAMELIST_FILE='./nml.breeding'


CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / BREEDING /                          &
 & norm_type , breeding_type , input_file_prefix_ctrl , input_file_prefix_pp , input_file_prefix_pn , &
 & input_file_prefix_pp_old , input_file_prefix_pn_old , breeding_in_place ,                          &
 & wind_scale , t_scale , p_scale ,                                                                   &
 & norm_threshold ,                                                                                   &
 & ave_norm_xy , ave_norm_z , is , ie , js , je , ks , ke ,                                           &
 & smooth_norm_xy , smooth_norm_z ,  smooth_sx , smooth_sz ,                                          &
 & init_bv , perturb_t , perturb_wind ,                                                               &
 & perturb_t_amp , perturb_wind_amp , perturb_t_sclh , perturb_t_sclv ,                               &
 & perturb_wind_sclh , perturb_wind_sclv  , bdy_width                                                                     

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=BREEDING,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at BREEDING section"
WRITE(*,*)"Using default values"
ENDIF
CLOSE(54)




END SUBROUTINE READ_NAMELIST

END MODULE common_namelist
