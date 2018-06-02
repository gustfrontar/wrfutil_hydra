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

  REAL(r_size) ::  var_local(nv3d+nv2d,nid_obs) = 1.0d0
  REAL(r_size) ::  var_local_par(np2d,nid_obs) = 1.0d0

  !GENERAL
  INTEGER :: nslots=2 ! number of time slots for 4D-LETKF
  INTEGER :: nbslot=2 ! basetime slot
  INTEGER :: nbv=40   ! Number of ensemble members

  !LOCALIZATION
  REAL(r_size) :: sigma_obs=2000.0d0
  REAL(r_size) :: sigma_obsv=0.2d0
  REAL(r_size) :: sigma_obsz=2000.0d0 !To perform vertical localization in Z.
  REAL(r_size) :: sigma_obst=3.0d0

  INTEGER      :: lev_update_q = 40  !Condensate won't be updated for levels over this level. 

  REAL(r_size) :: var_local_uv(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_t(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_tv(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_moist(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_ps(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_ref(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_dop(nv3d+nv2d)=1.0d0
  REAL(r_size) :: var_local_co(nv3d+nv2d)=1.0d0

  REAL(r_size) :: var_localp_uv(np2d)=1.0d0
  REAL(r_size) :: var_localp_t(np2d)=1.0d0
  REAL(r_size) :: var_localp_tv(np2d)=1.0d0
  REAL(r_size) :: var_localp_moist(np2d)=1.0d0
  REAL(r_size) :: var_localp_ps(np2d)=1.0d0
  REAL(r_size) :: var_localp_ref(np2d)=1.0d0
  REAL(r_size) :: var_localp_dop(np2d)=1.0d0
  REAL(r_size) :: var_localp_co(np2d)=1.0d0
 
  !OBSERVATIONS
  REAL(r_size) :: threshold_dz=1000.0d0
  REAL(r_size) :: gross_error=15.0d0
  REAL(r_size) :: gross_error_tycll=10.0d0    ! degree
  REAL(r_size) :: gross_error_tycmip=60.0d2   ! Pa
  REAL(r_size) :: undef_obs=9.99d9            ! Missing obs code
  REAL(r_size) :: rainratio_threshold = 0.3d0 !Percentaje of ensemble members with reflectivity greather than 0.
  LOGICAL      :: force_norain_assimilation = .false. 

  !PARAMETER ESTIMATION
  LOGICAL :: estpar = .false. 
  LOGICAL :: smooth_par_update_flag = .false.
  INTEGER :: update_parameter_2d(np2d) = 0     !WHETER PARAMETERS WILL BE ESTIMATED OR NOT (HFX_FACTOR, QFX_FACTOR, UST_FACTOR) 
  INTEGER :: update_parameter_0d(np2d) = 0     !WHETER PARAMETERS WILL BE ESTIMATED AS OD(1) PARAMETERS.
  INTEGER :: parameter_localization_type_0d(np2d)=0  !1 means localize as a surface variable
                                                                !2 means no vertical localization.
  INTEGER :: parameter_localization_type_2d(np2d)=1  !1 mean localize as a surface variable
                                                                !2 means no vertical localization.
  REAL(r_size)  :: param_sprd_init(np2d)= 0.05       !This is the initial parameter spread.
  LOGICAL       :: TRANSPAR=.false.                          !Wheter parameters will be transformed using the hiperbolic tangent function.
  LOGICAL       :: ADDINFPAR=.false.                         !Wheter to use additive inflation for the parameters.
  INTEGER       :: parameter_inflation_type=1                !1=Constant parameter sprea
                                                                        !2=Use state 2d variables estimated inflation, 
                                                                        !3= Use Parameter estimated inflation.
  REAL(r_size)  :: parameter_fixinflation(np2d)=1.1
  REAL(r_size)  :: additive_inflation_factor=0.1d0           !Additive noise as a percentage of the total parameter spread.
  REAL(r_size)  :: param_default_value(np2d)=1.0 !If parameter is not estimated then this value will be assigned to parameter analysis.
  REAL(r_size)  :: param_max_value(np2d)= 2.0  !This is the maximum allowed value for the parameters.
  REAL(r_size)  :: param_min_value(np2d)= 0.5  !This is the minimum allowed value for the parameters.

  !INFLATION
  REAL(r_size) :: cov_infl_mul = 1.1d0 !multiplicative inflation
                                                 ! > 1: globally constant covariance inflation
                                                 ! < 0: 3D inflation values input from a GPV file "infl_mul.grd"
  REAL(r_size) :: sp_infl_add = 0.d0   !additive inflation
  REAL(r_size) :: relax_alpha_spread = 0.0d0 !Relaxation to prior (Hamill and Whitaker 2012)
  REAL(r_size) :: relax_alpha        = 0.0d0 !Relaxation to prior (Zhang et al 2005)
  REAL(r_size) :: min_infl_mul       = 1.0d0 !Minimum multiplicative covariance inflation.

  !RADAR DA
  INTEGER            :: interpolation_technique=1
  INTEGER            :: nradar = 1               !Number of assimilated radars
  LOGICAL            :: use_wt=.true.            !Use or not terminal velocity in the computation of VR
  LOGICAL            :: use_pseudorh=.false.     !Enable pseudo RH observations.
  REAL(r_size)       :: minrefdbz=0.0d0          !Reflectivity values below this threshold won't be assimilated.
  REAL(r_size)       :: pseudorh_error=0.1       !Obserational error for pseudo RH observations.

  !OBSGRID SECTION
  REAL(r_size)       :: regrid_res=1.0d0         !Horizontal resolution of obsgrid grid (degree)
  REAL(r_size)       :: regrid_vert_res=100.0d0  !Vertical resolution of obsgrid grid (hPa)
  INTEGER, PARAMETER :: maxarea=20               !Maximum number of verification areas
  INTEGER            :: narea=1
  REAL(r_size)       :: vlon1(maxarea)=0.0d0
  REAL(r_size)       :: vlon2(maxarea)=360.0d0
  REAL(r_size)       :: vlat1(maxarea)=-90.0d0
  REAL(r_size)       :: vlat2(maxarea)=90.0d0
  LOGICAL            :: regrid_output=.false.
  LOGICAL            :: compute_stat=.false.
  LOGICAL            :: filter_input=.false.

  CHARACTER(LEN=20) :: NAMELIST_FILE='./letkf.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / GENERAL / nslots , nbslot , nbv
NAMELIST / LOCALIZATION / sigma_obs , sigma_obsv  , sigma_obsz , sigma_obst ,  &
 & var_local_uv , var_local_t , var_local_tv , var_local_moist , var_local_ps , var_local_ref , var_local_dop , var_local_co , &
 & var_localp_uv , var_localp_t , var_localp_tv , var_localp_moist , var_localp_ps , var_localp_ref , var_localp_dop ,         &
   var_localp_co  
NAMELIST / OBSERVATIONS / threshold_dz , gross_error , gross_error_tycll , gross_error_tycmip , &
&  undef_obs , force_norain_assimilation  
NAMELIST / PARAMETER_ESTIMATION / estpar , smooth_par_update_flag , update_parameter_2d , update_parameter_0d , &
 & parameter_localization_type_0d , parameter_localization_type_2d , param_sprd_init , transpar ,  &
 & addinfpar , parameter_inflation_type , parameter_fixinflation , additive_inflation_factor , &
 & param_default_value , param_min_value , param_max_value 
NAMELIST / INFLATION / cov_infl_mul , sp_infl_add , relax_alpha_spread , relax_alpha ,min_infl_mul
NAMELIST / RADAR_DA  / interpolation_technique ,  nradar , use_wt , use_pseudorh , pseudorh_error , minrefdbz , rainratio_threshold
NAMELIST / OBSGRID / nbv , narea , vlon1 , vlon2 , vlat1 , vlat2 , compute_stat, filter_input , regrid_output, regrid_res ,  &
                     regrid_vert_res


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
READ(54,NML=LOCALIZATION,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at LOCALIZATION section"
WRITE(*,*)"Using default values"
ENDIF

REWIND(54)
READ(54,NML=OBSERVATIONS,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at OBSERVATIONS section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=PARAMETER_ESTIMATION,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at PARAMETER_ESTIMATION section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=INFLATION,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at INFLATION section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=RADAR_DA,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at RADAR_DA section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

CLOSE(54)


  var_local(:,1)=var_local_uv(1:nv3d+nv2d)
  var_local(:,2)=var_local_t(1:nv3d+nv2d)
  var_local(:,3)=var_local_tv(1:nv3d+nv2d)
  var_local(:,4)=var_local_moist(1:nv3d+nv2d)
  var_local(:,5)=var_local_ps(1:nv3d+nv2d)
  var_local(:,6)=var_local_ref(1:nv3d+nv2d)
  var_local(:,7)=var_local_dop(1:nv3d+nv2d)
  var_local(:,8)=var_local_co(1:nv3d+nv2d)

  var_local_par(:,1)=var_localp_uv(1:np2d)
  var_local_par(:,2)=var_localp_t(1:np2d)
  var_local_par(:,3)=var_localp_tv(1:np2d)
  var_local_par(:,4)=var_localp_moist(1:np2d)
  var_local_par(:,5)=var_localp_ps(1:np2d)
  var_local_par(:,6)=var_localp_ref(1:np2d)
  var_local_par(:,7)=var_localp_dop(1:np2d)
  var_local_par(:,8)=var_localp_co(1:np2d)

END SUBROUTINE READ_NAMELIST

END MODULE common_namelist
