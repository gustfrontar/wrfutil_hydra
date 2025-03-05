#include "scalelib.h"
!OCL SERIAL
MODULE class_Cloud_Micro
  use scale_io
!  use tic_toc
  use maxdims
  use class_Group, only: &
     Group, &
     col_lut_aux, &
     vap_igp_aux
  use class_AirGroup, only: &
     AirGroup
  use class_Mass_Bin, only: &
     Mass_Bin, &
     data1d_lut, &
     data1d_lut_big
  use mod_amps_utility, only: &
     random_genvar
  use acc_amps
  use par_amps
  use mod_amps_const
  use com_amps, only: &
     fid_alog, &
     debug
  implicit none
  private

  public :: make_cloud_micro_cnfg
  public :: set_cloud_micro_cur
  public :: ini_cloud_micro
  public :: reality_check
  public :: check_water_apmass
  public :: cal_micro_tendency
  public :: print_cloud_tendency
  public :: return_output
  public :: update_terminal_vel
  public :: cal_dmtend
  public :: cal_dmtend_scale
  public :: cal_dmtend_kid
  public :: cal_dmtend_pm
  public :: cal_dmtend_z

  public :: Cloud_Micro

  TYPE Cloud_Micro
     ! memory is aligned as declared below for use in common block
     sequence
     ! *********************************************************************
     !                         Aerosol  Groups
     type (Group),dimension(4)  :: aerosol
!     type (Group),dimension(2)  :: aerosol
     !  1: Cloud Condensation Nuclei, CCN
     !  2: Ice Nucei, IN
     ! *********************************************************************

     ! *********************************************************************
     !                         Hydrometeor  Groups
     type (Group)       :: rain, solid_hydro
!     type (Group)       :: rain
     ! *********************************************************************

     !
     ! *********************************************************************
     !                     Thermodynamic property of air
     type (AirGroup)    :: air
     ! *********************************************************************

     ! message from reality_check subroutine
     ! 0. vapor and hydrometeors do not exist.
     ! 1. vapor only. no hydrometeors.
     ! 2. vapor and liquid hydrometeors exist.
     ! 3. vapor and solid hydrometeoors exist.
     ! 4. vapor, liquid and solid hydrometeors exist.
     !
     !integer,dimension(LMAX)                  :: mes_rc
     integer, allocatable                  :: mes_rc(:)
     !
     ! *********************************************************************
     !                              options
     ! level of complexity
     !   level : 5 = 12 prognostic variables
     !   level : 4 = 10 prognostic variables
     !   level : 3 = 10 prognostic variables
     !   level : 2 =  8 prognostic variables
     !   level : 1 =  4 prognostic variables
     integer     :: level_comp
     !
     ! debugging level
     integer     :: debug_level
     !
     ! collision level
     integer     :: coll_level
     !
     ! output from this microphysics scheme to dynamic model
     !   1 : only tendency is updated.
     !   2 : tendency and variable are updated.
     integer     :: out_type

     ! token (type of hydrometeor) for rain and solid hydro.
     integer     :: token_r, token_s, token_a
     !
     ! number of bins for rain, and solid hydrometeor groups
     integer     :: nbin_r, nbin_s, nbin_a
     !
     ! number of categories for aerosols
     integer     :: ncat_a
     !
     ! number of grids that have microphyscis in it.
     integer     :: L
     !
     ! distribution type inside of a bin
     integer     :: dtype_c, dtype_r, dtype_s, dtype_a(4)
     !
     ! hydrodynamic breakup method
     !   o number of method that you are using for hydrodynamic breakup
     !   1: The particle that became larger than maximum size shed the extra
     !       mass into the same bin (under construction!).
     !   2: Use distribution of fragments and probability of breakup.
     integer     :: hbreak_c, hbreak_r, hbreak_s
     !
     ! flag for prediction
     !  -1: concentration and mass for the group are zero.
     !   0: not predicted (initial setting).
     !   1: mixing ratio is predicted, and concentration is diagnosed.
     !   2: mixing ratio and concentration are predicted.
     integer     :: flagp_c, flagp_r, flagp_s, flagp_a
     !
     ! model type of CCN activation
     integer     :: act_type
     !
     ! time period for printing (seconds)
     integer     :: T_print_period

     ! random generaion: 1, max frequency: 0
     integer           :: ihabit_gm_random

     ! output format
     character(len=16) :: output_format
     !
     ! fixed concentration for cloud droplets
     real(PS)    :: fcon_c
     !
     !
     ! current time
     real(PS)        :: cur_time

     ! *********************************************************************

     ! coefficients for calculating total mass from total concentration
     real(PS),dimension(4)      :: coef_ap

     ! fraction of soluble mass in aerosols
     real(PS),dimension(4)      :: eps_ap0

     ! bulk dry density of aerosol particles at initialization
     real(PS),dimension(4)      :: den_apt0
     real(PS),dimension(4)      :: den_aps0, den_api0


     ! total mass for each grid box
!tmp     real(PS),pointer,dimension(:) :: mapt_r,mapt_s,mt_s,mt_r
!tmp     real(PS),pointer,dimension(:,:) :: mapt_ac
     !real(PS),dimension(LMAX) :: mapt_r,mapt_s,mt_s,mt_r
     !real(PS),dimension(ncamx,LMAX) :: mapt_ac
     real(PS), allocatable :: mapt_r(:), mapt_s(:), mt_s(:), mt_r(:)
     real(PS), allocatable :: mapt_ac(:,:)

     ! execution flag
     integer,dimension(20)     :: micexfg

     ! collisional-breakup variables
     integer :: imin_bk,imax_bk,jmin_bk,jmax_bk
! this is for 30 bins
!     real(PS) :: bu_fd(2,17835),bu_tmass(435)
! this is for 80 bins
     real(PS) :: bu_fd(2,62400),bu_tmass(780)

     real(PS) :: CCNMAX, CRIC_RN_IMM, frac_dust

     ! aerosol variables
     real(PS),dimension(4) :: nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean &
              ,ap_sig_cp,ap_mean_cp,cdf_cp_0,cdf_cp_180m0

     ! collision LUT
     type (col_lut_aux) :: adrpdrp
     real(PS),dimension(201,201) :: drpdrp
     type (col_lut_aux) :: ahexdrp
     real(PS),dimension(64,71) :: hexdrp
     type (col_lut_aux) :: abbcdrp
     real(PS),dimension(64,71) :: bbcdrp
     type (col_lut_aux) :: acoldrp
     real(PS),dimension(62,71) :: coldrp
     type (col_lut_aux) :: agp1drp
     real(PS),dimension(37,125) :: gp1drp
     type (col_lut_aux) :: agp4drp
     real(PS),dimension(27,125) :: gp4drp
     type (col_lut_aux) :: agp8drp
     real(PS),dimension(21,125) :: gp8drp

     ! chemical notation for aerosols
     character (len=16) :: APSNAME(4)

     ! time step for dynamics
     real(PS) :: dt_dy

     ! time step for collection and vapor dep processes
     real(PS) :: dt_cl,dt_vp
     ! number of steps over the dynamical time step
     integer :: n_step_cl, n_step_vp

     ! Inherent Growth parameterization
     type (vap_igp_aux) :: vigp

     ! Osmotic coefficients
     type (data1d_lut) :: osm_nhs4,osm_sdch

     ! random generator vars
     type (random_genvar) :: rdsd

     ! lookup table for standard normal distribution
     type (data1d_lut_big) :: snrml

     ! lookup table for invserse standard normal distribution
     type (data1d_lut_big) :: isnrml
  end type Cloud_Micro

CONTAINS

  subroutine make_Cloud_Micro_cnfg( &
       CM &
      ,dt_step, n_step_cl, n_step_vp, c_time  &
      ,NRTYPE,NRBIN, NRCAT  &
      ,NSTYPE,NSBIN, NSCAT  &
      ,NATYPE,NABIN, NACAT  &
      ,estbar,esitbar  &
      ,level_comp,debug_level,coll_level_IN,out_type,T_print_period,output_format &
      ,token_r,token_s,token_a,dtype_c,dtype_r,dtype_s,dtype_a &
      ,hbreak_c,hbreak_r,hbreak_s,flagp_c,flagp_r,flagp_s,flagp_a &
      ,act_type,micexfg &
      ,fcon_c,coef_ap,eps_ap &
      ,srat_c,srat_r,srat_s,srat_a &
      ,sadd_c,sadd_r,sadd_s,sadd_a &
      ,minmass_r,minmass_s,minmass_a &
      ,den_apt,den_aps,den_api &
      ,binbr,binbi &
      ,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_fd,bu_tmass &
      ,APSNAME,nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean &
      ,ap_sig_cp,ap_mean_cp,cdf_cp_0,cdf_cp_180m0 &
      ,nr_drpdrp,nc_drpdrp,xs_drpdrp,dx_drpdrp,ys_drpdrp,dy_drpdrp,drpdrp &
      ,nr_hexdrp,nc_hexdrp,xs_hexdrp,dx_hexdrp,ys_hexdrp,dy_hexdrp,hexdrp &
      ,nr_bbcdrp,nc_bbcdrp,xs_bbcdrp,dx_bbcdrp,ys_bbcdrp,dy_bbcdrp,bbcdrp &
      ,nr_coldrp,nc_coldrp,xs_coldrp,dx_coldrp,ys_coldrp,dy_coldrp,coldrp &
      ,nr_gp1drp,nc_gp1drp,xs_gp1drp,dx_gp1drp,ys_gp1drp,dy_gp1drp,gp1drp &
      ,nr_gp4drp,nc_gp4drp,xs_gp4drp,dx_gp4drp,ys_gp4drp,dy_gp4drp,gp4drp &
      ,nr_gp8drp,nc_gp8drp,xs_gp8drp,dx_gp8drp,ys_gp8drp,dy_gp8drp,gp8drp &
      ,nok_igp,x_igp,a_igp,b_igp &
      ,n_osm_nh42so4,xs_osm_nh42so4,dx_osm_nh42so4,y_osm_nh42so4 &
      ,n_osm_sodchl,xs_osm_sodchl,dx_osm_sodchl,y_osm_sodchl &
      ,n_snrml,xs_snrml,dx_snrml,y_snrml &
      ,n_isnrml,xs_isnrml,dx_isnrml,y_isnrml &
      ,ihabit_gm_random &
      ,CCNMAX,CRIC_RN_IMM,frac_dust &
      ,lbin)
!!c       lbin) result (CM)
!!c       binbr,binbi) result (CM)
    use class_AirGroup, only: &
       make_AirGroup_alloc
    use class_Group, only: &
       make_Group, &
       make_col_lut

    type (Cloud_Micro),intent(inout)    :: CM
!tmp    type (Cloud_Micro)    :: CM
!tmp    COMMON /ICMPRIV/CM
!tmp!$omp threadprivate(/ICMPRIV/)

    ! model time step in second
    real(DS), intent(in) :: dt_step
    real(PS) :: time_step
    ! number of time steps relative to one dynamic time step
    integer, intent(in) :: n_step_cl, n_step_vp
    ! current model time in second
    real(DS), intent(in) :: c_time
    !
    integer, intent(in) :: NRBIN, NSBIN, NABIN
    integer, intent(in) :: NRTYPE, NSTYPE, NATYPE
    integer, intent(in) :: NRCAT, NSCAT, NACAT
    integer, intent(in) :: lbin

    ! random generaion: 1, max frequency: 0
    integer, intent(in) :: ihabit_gm_random

    real(DS), intent(in) :: estbar(150), esitbar(111)

    integer, intent(in) :: coll_level_IN

    integer :: level_comp,debug_level,out_type,T_print_period,&
       token_r,token_s,token_a,dtype_c,dtype_r,dtype_s,dtype_a(*),&
       hbreak_c,hbreak_r,hbreak_s,flagp_c,flagp_r,flagp_s,flagp_a,&
       act_type,micexfg(*)
    real(PS) :: srat_c,srat_r,srat_s,srat_a,sadd_c,sadd_r,sadd_s,sadd_a,&
       fcon_c,minmass_r,minmass_s,minmass_a,&
       den_apt(*),den_aps(*),den_api(*),&
       binbr(*),binbi(*),&
       coef_ap(*),eps_ap(*)
    real(PS), intent(in) :: CCNMAX, CRIC_RN_IMM, frac_dust

    ! output format
    character(len=16), intent(in) :: output_format

    ! bin boundaries for aerosols have to be introduced here.
    real(PS) :: binba(82)

    ! collisional-breakup variables
    integer, intent(in) :: imin_bk,imax_bk,jmin_bk,jmax_bk
    real(PS) :: bu_fd(2,*),bu_tmass(*)
    ! chemical notation for aerosols
    character (len=16) :: APSNAME(*)
    ! aerosol variables
    real(PS),dimension(*) :: nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean &
           ,ap_sig_cp,ap_mean_cp,cdf_cp_0,cdf_cp_180m0

    ! collision variables
    integer,intent(in) :: nr_drpdrp,nc_drpdrp
    real(PS),intent(in) :: xs_drpdrp,dx_drpdrp,ys_drpdrp,dy_drpdrp
    real(PS),dimension(nr_drpdrp,*) :: drpdrp
    integer,intent(in) :: nr_hexdrp,nc_hexdrp
    real(PS),intent(in) :: xs_hexdrp,dx_hexdrp,ys_hexdrp,dy_hexdrp
    real(PS),dimension(nr_hexdrp,*) :: hexdrp
    integer,intent(in) :: nr_bbcdrp,nc_bbcdrp
    real(PS),intent(in) :: xs_bbcdrp,dx_bbcdrp,ys_bbcdrp,dy_bbcdrp
    real(PS),dimension(nr_bbcdrp,*) :: bbcdrp
    integer,intent(in) :: nr_coldrp,nc_coldrp
    real(PS),intent(in) :: xs_coldrp,dx_coldrp,ys_coldrp,dy_coldrp
    real(PS),dimension(nr_coldrp,*) :: coldrp
    integer,intent(in) :: nr_gp1drp,nc_gp1drp
    real(PS),intent(in) :: xs_gp1drp,dx_gp1drp,ys_gp1drp,dy_gp1drp
    real(PS),dimension(nr_gp1drp,*) :: gp1drp
    integer,intent(in) :: nr_gp4drp,nc_gp4drp
    real(PS),intent(in) :: xs_gp4drp,dx_gp4drp,ys_gp4drp,dy_gp4drp
    real(PS),dimension(nr_gp4drp,*) :: gp4drp
    integer,intent(in) :: nr_gp8drp,nc_gp8drp
    real(PS),intent(in) :: xs_gp8drp,dx_gp8drp,ys_gp8drp,dy_gp8drp
    real(PS),dimension(nr_gp8drp,*) :: gp8drp

    integer,intent(in) :: nok_igp
    real(PS),dimension(*),intent(in) :: x_igp
    real(PS),dimension(nok_igp,*),intent(in) :: a_igp,b_igp

    integer,intent(in) :: n_osm_nh42so4
    real(PS),intent(in) :: xs_osm_nh42so4,dx_osm_nh42so4
    real(PS),dimension(*),intent(in) :: y_osm_nh42so4
    integer,intent(in) :: n_osm_sodchl
    real(PS),intent(in) :: xs_osm_sodchl,dx_osm_sodchl
    real(PS),dimension(*),intent(in) :: y_osm_sodchl

    integer,intent(in) :: n_snrml
    real(PS),intent(in) :: xs_snrml,dx_snrml
    real(PS),dimension(*),intent(in) :: y_snrml

    integer,intent(in) :: n_isnrml
    real(PS),intent(in) :: xs_isnrml,dx_isnrml
    real(PS),dimension(*),intent(in) :: y_isnrml

    integer :: i!, var_Status,j
    ! max element for breakup arrays
    integer :: i1d_pair_max,kk_max
!!$    integer :: omp_in_parallel,omp_get_num_threads,omp_get_thread_num
!!$    integer(8) :: loc

!!c    write(fid_alog,*) "in make_cloud, initialization"

    allocate( CM%mes_rc(LMAX) )
    allocate( CM%mapt_r(LMAX) )
    allocate( CM%mapt_s(LMAX) )
    allocate( CM%mt_s(LMAX) )
    allocate( CM%mt_r(LMAX) )
    allocate( CM%mapt_ac(ncamx,LMAX) )

    ! set time step
    time_step = dt_step

    CM%dt_dy = dt_step
    !  time step for collection process
    CM%dt_cl = dt_step/real(n_step_cl,PS_KIND)
    !  time step for vapor deposition process
    !CM%dt_vp = dt_step/real(n_step_vp)
    CM%dt_vp = CM%dt_cl/real(n_step_vp,PS_KIND)
    CM%n_step_cl=n_step_cl
    CM%n_step_vp=n_step_vp


    ! +++ input options +++
    CM%level_comp=level_comp
    CM%debug_level=debug_level
    CM%coll_level=coll_level_IN
    CM%out_type=out_type
    CM%T_print_period=T_print_period
    CM%output_format=output_format
    CM%token_r=token_r; CM%token_s=token_s; CM%token_a=token_a
    CM%nbin_r=NRBIN; CM%nbin_s=NSBIN; CM%nbin_a=NABIN
    CM%ncat_a=NACAT
    CM%dtype_c=dtype_c; CM%dtype_r=dtype_r; CM%dtype_s=dtype_s
    CM%dtype_a(1:nacat)=dtype_a(1:nacat)
    CM%hbreak_c=hbreak_c;CM%hbreak_r=hbreak_r;CM%hbreak_s=hbreak_s
    CM%flagp_c=flagp_c;CM%flagp_r=flagp_r;CM%flagp_s=flagp_s;CM%flagp_a=flagp_a
    CM%act_type=act_type
!!c    CM%srat_c=srat_c;CM%srat_r=srat_r;CM%srat_s=srat_s;CM%srat_a=srat_a
!!c    CM%sadd_c=sadd_c;CM%sadd_r=sadd_r;CM%sadd_s=sadd_s;CM%sadd_a=sadd_a
    CM%fcon_c=fcon_c
!!c    CM%minmass_c=minmass_c; CM%minmass_r=minmass_r; CM%minmass_s=minmass_s
!!c    CM%minmass_a=minmass_a
    CM%coef_ap(1:nacat)=coef_ap(1:nacat)

    ! ++++++ current time  +++++
    CM%cur_time = c_time
!!c    write(fid_alog,*) "time is ", CM%cur_time

    CM%eps_ap0(1:nacat)=eps_ap(1:nacat)
    CM%den_apt0(1:nacat)=den_apt(1:nacat)
    CM%den_aps0(1:nacat)=den_aps(1:nacat)
    CM%den_api0(1:nacat)=den_api(1:nacat)

!write(fid_alog,*) "densdebug1", CM%eps_ap0(1:nacat)
!write(fid_alog,*) "densdebug2", CM%den_apt0(1:nacat)
!write(fid_alog,*) "densdebug3", CM%den_aps0(1:nacat)
!write(fid_alog,*) "densdebug4", CM%den_api0(1:nacat)
!write(fid_alog,*) "densdebug5", CM%dt_dy, CM%dt_cl, CM%n_step_cl
!write(fid_alog,*) "densdebug6", CM%dt_vp, CM%n_step_vp

!!c    write(fid_alog,*) "0"

    CM%micexfg(1:20)=micexfg(1:20)

!!c    write(fid_alog,*) "1"
    CM%imin_bk=imin_bk
    CM%imax_bk=imax_bk
    CM%jmin_bk=jmin_bk
    CM%jmax_bk=jmax_bk
!!c    write(fid_alog,*) "2"
    i1d_pair_max=(imax_bk-1)-jmin_bk+1+(imax_bk-imin_bk)*(1+imax_bk-imin_bk)/2
    kk_max=i1d_pair_max*NRBIN
    CM%bu_fd(1:2,1:kk_max)=bu_fd(1:2,1:kk_max)
!!c    write(fid_alog,*) "3"
    CM%bu_tmass(1:i1d_pair_max)=bu_tmass(1:i1d_pair_max)
!!c    write(fid_alog,*) "4"
    CM%APSNAME(1:nacat)=APSNAME(1:nacat)
!!c    write(fid_alog,*) "5"
    CM%nu_aps(1:nacat)=nu_aps(1:nacat)
!!c    write(fid_alog,*) "6"
    CM%phi_aps(1:nacat)=phi_aps(1:nacat)
!!c    write(fid_alog,*) "7"
    CM%m_aps(1:nacat)=m_aps(1:nacat)
!!c    write(fid_alog,*) "8"
    CM%ap_lnsig(1:nacat)=ap_lnsig(1:nacat)
!!c    write(fid_alog,*) "9"
    CM%ap_mean(1:nacat)=ap_mean(1:nacat)

    CM%ap_sig_cp(1:nacat)=ap_sig_cp(1:nacat)
    CM%ap_mean_cp(1:nacat)=ap_mean_cp(1:nacat)
    CM%cdf_cp_0(1:nacat)=cdf_cp_0(1:nacat)
    CM%cdf_cp_180m0(1:nacat)=cdf_cp_180m0(1:nacat)

    CM%CCNMAX = CCNMAX
    CM%CRIC_RN_IMM = CRIC_RN_IMM
    CM%frac_dust = frac_dust

    CM%adrpdrp=make_col_lut(nr_drpdrp,nc_drpdrp,xs_drpdrp,ys_drpdrp &
                           ,dx_drpdrp,dy_drpdrp)
    CM%drpdrp(1:nr_drpdrp,1:nc_drpdrp)=drpdrp(1:nr_drpdrp,1:nc_drpdrp)

    CM%ahexdrp=make_col_lut(nr_hexdrp,nc_hexdrp,xs_hexdrp,ys_hexdrp &
                           ,dx_hexdrp,dy_hexdrp)
    CM%hexdrp(1:nr_hexdrp,1:nc_hexdrp)=hexdrp(1:nr_hexdrp,1:nc_hexdrp)

    CM%abbcdrp=make_col_lut(nr_bbcdrp,nc_bbcdrp,xs_bbcdrp,ys_bbcdrp &
                           ,dx_bbcdrp,dy_bbcdrp)
    CM%bbcdrp(1:nr_bbcdrp,1:nc_bbcdrp)=bbcdrp(1:nr_bbcdrp,1:nc_bbcdrp)

    CM%acoldrp=make_col_lut(nr_coldrp,nc_coldrp,xs_coldrp,ys_coldrp &
                           ,dx_coldrp,dy_coldrp)
    CM%coldrp(1:nr_coldrp,1:nc_coldrp)=coldrp(1:nr_coldrp,1:nc_coldrp)

    CM%agp1drp=make_col_lut(nr_gp1drp,nc_gp1drp,xs_gp1drp,ys_gp1drp &
                           ,dx_gp1drp,dy_gp1drp)
    CM%gp1drp(1:nr_gp1drp,1:nc_gp1drp)=gp1drp(1:nr_gp1drp,1:nc_gp1drp)

    CM%agp4drp=make_col_lut(nr_gp4drp,nc_gp4drp,xs_gp4drp,ys_gp4drp &
                           ,dx_gp4drp,dy_gp4drp)
    CM%gp4drp(1:nr_gp4drp,1:nc_gp4drp)=gp4drp(1:nr_gp4drp,1:nc_gp4drp)

    CM%agp8drp=make_col_lut(nr_gp8drp,nc_gp8drp,xs_gp8drp,ys_gp8drp &
                           ,dx_gp8drp,dy_gp8drp)
    CM%gp8drp(1:nr_gp8drp,1:nc_gp8drp)=gp8drp(1:nr_gp8drp,1:nc_gp8drp)

    CM%vigp%nok=nok_igp
    CM%vigp%x(1:nok_igp)=x_igp(1:nok_igp)
    CM%vigp%a(1:nok_igp,1:4)=a_igp(1:nok_igp,1:4)
    CM%vigp%b(1:nok_igp,1:4)=b_igp(1:nok_igp,1:4)

    CM%osm_nhs4%n=n_osm_nh42so4
    CM%osm_nhs4%dx=dx_osm_nh42so4
    CM%osm_nhs4%xs=xs_osm_nh42so4
    CM%osm_nhs4%y(1:n_osm_nh42so4)=y_osm_nh42so4(1:n_osm_nh42so4)
    CM%osm_sdch%n=n_osm_sodchl
    CM%osm_sdch%dx=dx_osm_sodchl
    CM%osm_sdch%xs=xs_osm_sodchl
    CM%osm_sdch%y(1:n_osm_sodchl)=y_osm_sodchl(1:n_osm_sodchl)

    CM%snrml%n=n_snrml
    CM%snrml%dx=dx_snrml
    CM%snrml%xs=xs_snrml
    CM%snrml%y(1:n_snrml)=y_snrml(1:n_snrml)

    CM%isnrml%n=n_isnrml
    CM%isnrml%dx=dx_isnrml
    CM%isnrml%xs=xs_isnrml
    CM%isnrml%y(1:n_isnrml)=y_isnrml(1:n_isnrml)

    CM%rdsd%jseed=0
    CM%rdsd%ifrst=0
    CM%rdsd%isect_seed=0
    CM%rdsd%nextn=0

    CM%ihabit_gm_random=ihabit_gm_random

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! +++ allocate memory to each group +++
!tmp    nullify(CM%mes_rc)
!tmp    allocate( CM%mes_rc(lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mes_rc &
!tmp         in class_cloud_micro"
!tmp
!tmp    nullify(CM%mapt_r)
!tmp    allocate( CM%mapt_r(lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mapt_r &
!tmp         in class_cloud_micro"
!tmp
!tmp    nullify(CM%mapt_s)
!tmp    allocate( CM%mapt_s(lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mapt_s &
!tmp         in class_cloud_micro"
!tmp
!tmp    nullify(CM%mapt_ac)
!tmp    allocate( CM%mapt_ac(2,lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mapt_ac &
!tmp         in class_cloud_micro"
!tmp
!tmp    nullify(CM%mt_r)
!tmp    allocate( CM%mt_r(lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mt_r &
!tmp         in class_cloud_micro"
!tmp
!tmp    nullify(CM%mt_s)
!tmp    allocate( CM%mt_s(lbin), stat = var_Status)
!tmp    if(var_Status /= 0 ) stop "Memory not available for CM%mt_s &
!tmp         in class_cloud_micro"

    CM%mes_rc=1

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!!c    write(fid_alog,*) "bf make air"
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! construct airgroup object
    CM%air = make_AirGroup_alloc(lbin,estbar,esitbar)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! construct Group object
    ! rain object
!tmp    write(fid_alog,*) "lbin",lbin
!tmp    write(fid_alog,*) "bf make rain ice, token:",CM%token_r,CM%token_s
!tmp    write(fid_alog,*) "bf make rain ice, nbin:",CM%nbin_r,CM%nbin_s

!tmp    CM%rain = make_Group (lbin, time_step, CM%token_r, CM%nbin_r, &
!tmp         srat_r, sadd_r, minmass_r, CM%dtype_r,1.0_PS,den_aps(1),den_api(1),CM%eps_ap0(1),binbr)
    call make_Group (CM%rain, lbin, time_step, CM%token_r, CM%nbin_r, &
         srat_r, sadd_r, minmass_r, CM%dtype_r,1.0_PS,den_aps(1),den_api(1),CM%eps_ap0(1),binbr)

!!c   write(fid_alog,*) "ck0 rain token",CM%rain%token

    ! solid hydrometer object
!!c    write(fid_alog,*) "bf make solid"
!tmp    CM%solid_hydro = make_Group (lbin, time_step, CM%token_s, CM%nbin_s, &
!tmp         srat_s, sadd_s, minmass_s, CM%dtype_s,den_i,den_aps(2),den_api(2),CM%eps_ap0(2),binbi)
     call make_Group (CM%solid_hydro, lbin, time_step, CM%token_s, CM%nbin_s, &
         srat_s, sadd_s, minmass_s, CM%dtype_s,den_i,den_aps(2),den_api(2),CM%eps_ap0(2),binbi)

!!c   write(fid_alog,*) "ck0 ice token",CM%solid_hydro%token

    ! aerosol object
!!c    write(fid_alog,*) "bf make aerosols"
!!c    write(fid_alog,*) "bin aerosol",CM%nbin_a,CM%dtype_a(1:2)
     binba(:) = 0.0_PS ! tentative
    do i=1,CM%ncat_a
!tmp       CM%aerosol(i) = make_Group (lbin, time_step, CM%token_a, CM%nbin_a, &
!tmp            srat_a, sadd_a, minmass_a, CM%dtype_a(i),den_apt(i),den_aps(i),den_api(i),CM%eps_ap0(i))
       call make_Group (CM%aerosol(i),lbin, time_step, CM%token_a, CM%nbin_a, &
            srat_a, sadd_a, minmass_a, CM%dtype_a(i),den_apt(i),den_aps(i),den_api(i),CM%eps_ap0(i),binba)
!!c      write(fid_alog,*) "ck0 ap token",CM%aerosol(i)%token
    end do
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!ito    if(c_time.le.dt_step) then
!ito       write(fid_alog,*) "make_cloud1",lbin, time_step, CM%token_r, CM%nbin_r, &
!ito         srat_r, sadd_r, minmass_r, CM%dtype_r
!ito       write(fid_alog,*) "timesteps",CM%dt_cl,CM%dt_vp,CM%n_step_cl,CM%n_step_vp
!ito       write(fid_alog,*) "make_cloud2",1.0_PS,&
!ito         den_aps(1),den_api(1),CM%eps_ap0(1)
!ito       write(fid_alog,*) "make_cloud3",binbr(1:CM%nbin_r+1)
!ito       write(fid_alog,*) "make_cloud4",CM%APSNAME,CM%m_aps
!ito       write(fid_alog,*) "make_cloud5",CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk
!ito       write(fid_alog,*) "ck rain,ice,aerosol token",CM%rain%token,CM%solid_hydro%token,CM%aerosol(1:CM%ncat_a)%token
!ito!!c       write(fid_alog,*) "bu_fd"
!ito!!c       write(fid_alog,*) CM%bu_fd
!ito!!c       write(fid_alog,*) "bu_tmass"
!ito!!c       write(fid_alog,*) CM%bu_tmass
!ito
!ito    endif
!!c    write(fid_alog,'("omp in make cnfg",3I5)') omp_in_parallel(),omp_get_num_threads(),omp_get_thread_num()
!!c    write(fid_alog,'("memory loc in make cnfg",5I15)') loc(CM%air),loc(CM%rain),loc(CM%solid_hydro),&
!!c                  loc(CM%rain%MS(1,1)),loc(CM%solid_hydro%MS(1,1))

!!c  end function make_Cloud_Micro_cnfg
  end subroutine make_Cloud_Micro_cnfg

  subroutine set_Cloud_Micro_cur(CM,c_time)
    type (Cloud_Micro),intent(inout)    :: CM
!tmp  subroutine set_Cloud_Micro_cur(c_time)
!tmp    type (Cloud_Micro)   :: CM
!tmp    COMMON /ICMPRIV/CM
!tmp!$omp threadprivate(/ICMPRIV/)
    ! current model time in second
    real(DS)               :: c_time

    ! ++++++ current time  +++++
    CM%cur_time = c_time
!!c    write(fid_alog,*) "time is ", CM%cur_time
    CM%mes_rc=1

  end subroutine set_Cloud_Micro_cur

!tmp  subroutine ini_cloud_micro( &
!tmp       XC &
  subroutine ini_cloud_micro( &
       CM &
      ,XC &
      ,XR, NRTYPE,NRBIN, NRCAT  &
      ,XS, NSTYPE,NSBIN, NSCAT  &
      ,XA, NATYPE,NABIN, NACAT  &
      ,RV, DEN, PT, T, W  &
      ,L,qtp,ID,JD,KD)
!tmp      ,c_time)
    use class_AirGroup, only: &
       ini_airgroup
    use class_Group, only: &
       ini_group_all

    type (Cloud_Micro),intent(inout)    :: CM
!tmp    type (Cloud_Micro)   :: CM
!tmp    COMMON /ICMPRIV/CM
!tmp!$omp threadprivate(/ICMPRIV/)

    ! C: cloud_drop
    ! R: rain
    ! S: solid hydrometeor
    ! NXBIN       : number of bins
    ! NXTYPE      : number of variables that each bin contain
    ! X(IBIN,ITYPE)    : quantity for ITYPE of IBIN-th bin
    ! DXDT(ITYPE) : tendency of the quantity
    ! ITYPE : 1. mixing ratio (g/g)
    !         2. concentration (#/g)
    !         3. volume of circumscribing sphere (cm^3/g)
    !         4. mixing ratio of a-axis production (g/g)
    !         5. mixing ratio of c-axis production (g/g)
    !         6. mixing ratio of d-axis production (g/g)
    !         7. mixing ratio of rossetta production (g/g)
    !         8. mixing ratio of riming production (g/g)
    !         9. volume of riming production (cm^3/g)
    !        10. mixing ratio of melt water (g/g)
    !
!tmp    integer, intent(in)    :: NRBIN, NSBIN, NABIN
!tmp    integer, intent(in)    :: NRTYPE, NSTYPE,NATYPE
!tmp    integer, intent(in)    :: NRCAT, NSCAT,NACAT
!tmp    integer, intent(in)    :: L
    integer    :: NRBIN, NSBIN, NABIN
    integer    :: NRTYPE, NSTYPE,NATYPE
    integer    :: NRCAT, NSCAT,NACAT
    integer    :: L
    real(MP_KIND)  :: XC(*)
    real(MP_KIND)  :: XR(NRTYPE,NRBIN,NRCAT,*), XS(NSTYPE,NSBIN,NSCAT,*),XA(NATYPE,NABIN,NACAT,*),qtp(*)

    integer :: ID(*),JD(*),KD(*)
    ! mixing ratio of vapor, density of ambient air, total pressure
    ! temperature, and vertical velocity
    real(MP_KIND)     :: RV(*), DEN(*), PT(*), T(*), W(*)
!tmp    ! current model time in second
!tmp    real,intent(in) :: c_time

    !integer :: i, var_Status,j
!!c    integer :: omp_in_parallel,omp_get_num_threads,omp_get_thread_num
!!c    integer(8) :: loc

!!c    write(fid_alog,*) "in ini_cloud_micro"

    ! +++ input options +++
!tmp    CM%cur_time=c_time
    CM%L=L
    CM%rain%L=L
    CM%solid_hydro%L=L
    CM%aerosol(1:CM%ncat_a)%L=L

    ! ++++++ current time  +++++
!!c    write(fid_alog,*) "omp_threads, thread num, time:", omp_get_num_threads(),omp_get_thread_num(),CM%cur_time

    CM%mes_rc=1

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize thermo_var object
    call ini_AirGroup(CM%air,L,RV,DEN,PT,T,W,CM%rdsd,CM%ihabit_gm_random)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfinigroup')

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! initialize the concentration and mass in each bin of all the groups
    call ini_group_all(CM%level_comp,CM%air,CM%rain,CM%solid_hydro,&
         CM%aerosol,CM%fcon_c,&
         L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
         XC,XR,XS,XA,CM%flagp_c,CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%coef_ap,&
         CM%eps_ap0,CM%den_apt0,CM%den_aps0,CM%den_api0,ID,JD,KD)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'afinigroup')

!ito    if(CM%cur_time.lt.CM%rain%dt) then
!ito       write(fid_alog,*) "ini_cloud1",CM%L, CM%rain%dt, CM%token_r, CM%nbin_r, &
!ito         CM%rain%a_b, CM%rain%b_b, CM%rain%mbinb, CM%rain%org_dtype
!ito       write(fid_alog,*) "ini_cloud2",1.0_PS,&
!ito         CM%rain%MS(1,1)%den_as,CM%rain%MS(1,1)%den_ai,CM%eps_ap0(1)
!ito       write(fid_alog,*) "ini_cloud3",CM%rain%binb
!ito       write(fid_alog,*) "ini_cloud4",CM%APSNAME,CM%m_aps
!ito       write(fid_alog,*) "ini_cloud5",CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk
!ito    endif
!!c       write(fid_alog,'("omp in ini_cloud",3I5)') omp_in_parallel(),omp_get_num_threads(),omp_get_thread_num()
!!c       write(fid_alog,'("memory loc in ini_cloud",5I15)') loc(CM%air),loc(CM%rain),loc(CM%solid_hydro),&
!!c                  loc(CM%rain%MS(1,1)),loc(CM%solid_hydro%MS(1,1))
  end subroutine ini_Cloud_Micro

!tmp  subroutine delete_Cloud_Micro (CM)
!tmp    type (Cloud_Micro), intent(inout)   :: CM
!tmp    integer                     :: ok
!tmp    integer                     :: i
!tmp
!tmp
!tmp    ! +++ deallocate memory used for the Group Object +++
!tmp    call delete_AirGroup( CM%air)
!tmp    call delete_Group( CM%solid_hydro)
!tmp    call delete_Group( CM%rain)
!tmp    do i=1,2
!tmp       call delete_Group( CM%aerosol(i))
!tmp    end do
!tmp
!tmp
!tmp    deallocate(CM%mt_s)
!tmp    nullify(CM%mt_s)
!tmp    deallocate(CM%mt_r)
!tmp    nullify(CM%mt_r)
!tmp    deallocate(CM%mapt_ac)
!tmp    nullify(CM%mapt_ac)
!tmp    deallocate(CM%mapt_s)
!tmp    nullify(CM%mapt_s)
!tmp    deallocate(CM%mapt_r)
!tmp    nullify(CM%mapt_r)
!tmp    deallocate( CM%mes_rc)
!tmp    nullify(CM%mes_rc)
!tmp
!tmp
!tmp    ! +++ deallocate memory used for the lookup_table object +++
!tmp!!c    call delete_Lookup_Table( CM%LT)
!tmp
!tmp  end subroutine delete_Cloud_Micro
!mark2
  subroutine cal_micro_tendency(CM,ID,JD,KD,qtp,thil,iproc,istrt &
                               ,LL,dmtendl,dcontendl,dbintendl &
! <<< 2014/10 T. Hashino added for KiD
!                 ,dM_auto_liq,dM_accr_liq &
!                 ,dM_auto_ice,dM_accr_ice &
!                 ,dM_auto_rim,dM_accr_rim)
!                 ,cptime_mtd)
                               )

    use scale_prc, only: &
       PRC_abort
    use class_AirGroup, only: &
       update_AirGroup
    use class_Group, only: &
       update_group_all, &
       ini_tendency, &
       ini_tendency_ag, &
       update_mesrc, &
       update_modelvars_all
    use mod_amps_core, only: &
       coalescence, &
       melting_shedding, &
       hydrodyn_breakup, &
       ice_nucleation1, &
       ice_nucleation2, &
       cal_aptact_var8_kc04dep, &
       cal_aptact_var8_vec, &
       vapor_deposition, &
       mv_ice2liq, &
       diag_t, &
       diag_pq
    use mod_amps_check, only: &
       repair
    implicit none
! >>> 2014/10 T. Hashino added for KiD
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate tendencies by cloud microphysical processes
    !
    ! This version marches the collision processes and vapor-related processes
    ! with seperate time step, and update the particle properties as it proceeds.
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Cloud_Micro), intent(inout)    :: CM
    integer :: ID(*),JD(*),KD(*)
    real(MP_KIND)    :: qtp(*),thil(*)
    integer,intent(in) :: iproc,istrt,LL
    real(MP_KIND), intent(inout) :: dmtendl(10,2,LL), dcontendl(10,2,LL), dbintendl(3,2,mxnbin,LL)
! <<< 2014/10 T. Hashino added for KiD
    real(PS), dimension(LMAX) :: &
                  dM_auto_liq,dM_accr_liq &
                 ,dM_auto_ice,dM_accr_ice &
                 ,dM_auto_rim,dM_accr_rim
! >>> 2014/10 T. Hashino added for KiD
    !real,dimension(20) :: cptime_mtd
    !
    !real,dimension(20) :: cptime
    !real :: cpsum
    !character(len=20),dimension(20)  :: proname
    integer,dimension(mxntend) :: iupdate_gr,iupdate_gs
    integer,dimension(mxntend,CM%ncat_a) :: iupdate_ga
    integer :: i,it_cl,it_vp!,n

    !
    real(PS),dimension(LMAX) :: T_a_r, RV_r
    ! error message
    integer   :: em
    !real(PS) :: s1,r1

    !integer :: i_gg, j_gg ! CHIARUI
    !type (Group), save :: ggg1, ggg2 ! CHIARUI

    integer :: ierr

    ! initialize
    em=0
    !cptime=0.0


!dbg    write(fid_alog,'("current time",F12.2)') CM%cur_time
    if(debug .and. CM%cur_time.eq.0.0) then
       write(fid_alog,'("micexfg",20I2)') CM%micexfg(1:20)
    endif

    RV_r(1:CM%L)=CM%air%TV(1:CM%L)%rv

!!    n=178
!!    i=7
!!    call denchk1(CM%solid_hydro,i,n,id,jd,kd,'in_caltend')
!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfcol')

    col_loop1: do it_cl=1,CM%n_step_cl
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 1. Collection processes
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! set time step size
      CM%rain%dt=CM%dt_cl
      CM%solid_hydro%dt=CM%dt_cl
      do i=1,CM%ncat_a
        CM%aerosol(i)%dt=CM%dt_cl
      enddo
      ! set the tendency mask for repair and update
      !                           1 2 3 4 5 6 7 8 9 10 11 12
          iupdate_gr(1:mxntend)=(/0,1,0,1,0,1,0,1,1, 1, 1, 1/)
          iupdate_gs(1:mxntend)=(/0,1,0,1,0,1,0,1,1, 1, 1, 1/)
      do i=1,CM%ncat_a
        iupdate_ga(1:mxntend,i)=(/0,1,0,1,0,0,0,1,1, 1, 1, 1/)
      enddo
      ! initialize all tendencies that are being updated.
      call ini_tendency(CM%rain,iupdate_gr,1)
      call ini_tendency(CM%solid_hydro,iupdate_gs,1)
      do i=1,CM%ncat_a
        call ini_tendency(CM%aerosol(i),iupdate_ga(1,i),1)
      enddo
      call ini_tendency_ag(CM%air,3)

      !call tic(s1,r1)
!      call PROF_rapstart('amps_update_diagnose',3)
      if(it_cl>1) then
        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)
        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

      endif
      !cptime(1)=cptime(1)+toc(s1,r1)
    ! diagnose physical quantities of each group object
!!c    write(fid_alog,*) "diag"
!!c      write(fid_alog,*) "ck rdsd01",CM%rdsd
      !call tic(s1,r1)
      if( CM%flagp_s > 0 ) &
        call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps, &
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      !cptime(2)=cptime(2)+toc(s1,r1)


!!c      write(fid_alog,*) "ck rdsd02",CM%rdsd
      !call tic(s1,r1)
      if(CM%level_comp>=6.and.it_cl==1) then
!         ! move 75% melt ice particles
        call mv_ice2liq(CM%solid_hydro,CM%rain,CM%aerosol,CM%air,CM%mes_rc)

        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)

        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

        if( CM%flagp_s > 0 ) &
          call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps, &
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      endif
      !cptime(3)=cptime(3)+toc(s1,r1)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfdiagpq2')

      !call tic(s1,r1)
      if( CM%flagp_r > 0 ) &
        call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
      !cptime(4)=cptime(4)+toc(s1,r1)

      !call tic(s1,r1)
      if( CM%flagp_a /= 0 ) then
        do i=1,CM%ncat_a
          call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                        i,CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
        end do
      end if
      !cptime(5)=cptime(5)+toc(s1,r1)
      if( em /= 0 ) then
        LOG_ERROR("cal_micro_tendency",*) "diag_pg > error!"
        call PRC_abort
      endif
!      call PROF_rapend('amps_update_diagnose',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !    PRINT OUT
      !call tic(s1,r1)
      !if(it_cl==1) then
      !  if(CM%micexfg(1)==1.and.mod(int(CM%cur_time*100), CM%T_print_period*100)==0) then

!!c        write(fid_alog,*) "print rain"
      !    if( CM%flagp_r > 0 ) &
      !      call print_out( CM%rain, CM%air, CM%cur_time, CM%output_format,&
      !                      ID,JD,KD,iproc,istrt)
!!c        write(fid_alog,*) "print aerosols"
      !    if( CM%flagp_a /= 0 ) &
      !      call print_out_ap( CM%aerosol, CM%air, CM%cur_time, CM%output_format,&
      !                         ID,JD,KD,iproc,istrt,CM%ncat_a )
!!c        write(fid_alog,*) "print ices"
      !    if( CM%flagp_s > 0 ) &
      !      call print_out( CM%solid_hydro, CM%air, CM%cur_time, CM%output_format,&
      !                      ID,JD,KD,iproc,istrt)
      !  end if
      !endif
      !cptime(20)=cptime(20)+toc(s1,r1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfcoale1')
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! collision-coalescence process
      ! NOTE:
!!c    write(fid_alog,*) "coalescence"
      ! +++ rain and rain process +++
      !call tic(s1,r1)
!      call PROF_rapstart('amps_collisioncoal',3)
      if(CM%micexfg(2)==1.and.CM%flagp_r > 0 ) &
         call coalescence( CM%rain,CM%rain,CM%air,CM%level_comp,CM%coll_level,CM%mes_rc &
                ,CM%micexfg(18) &  ! 2014/10 T. Hashino changed for KID
                ,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_liq,dM_accr_liq)
! >>> 2014/10 T. Hashino added for KiD

!!c         call coalescence( CM%rain,CM%rain,CM%air,CM%level_comp,CM%mes_rc,CM%cur_time,0,&
!!c                 CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd)
      !cptime(6)=cptime(6)+toc(s1,r1)
!      call PROF_rapend('amps_collisioncoal',3)

      ! +++ auto conversion of cloud droplet bin +++
      !call tic(s1,r1)
!tmp      if(CM%micexfg(12)==1.and.CM%flagp_r > 0 ) &
!tmp           call auto_conversion( CM%rain,CM%air,CM%level_comp,CM%mes_rc)
      !cptime(7)=cptime(7)+toc(s1,r1)


      ! +++ aggregation process +++
!!c    write(fid_alog,*) "in of agg"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_aggregation',3)
      if(CM%micexfg(3)==1.and.CM%flagp_s > 0 ) &
         call coalescence(CM%solid_hydro,CM%solid_hydro,CM%air,CM%level_comp,CM%coll_level,CM%mes_rc &
                ,0,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_ice,dM_accr_ice)
! >>> 2014/10 T. Hashino added for KiD
      !cptime(8)=cptime(8)+toc(s1,r1)
!      call PROF_rapend('amps_aggregation',3)
!!c    write(fid_alog,*) "out of agg"
!!c    call check_tendency_ap(CM%solid_hydro,CM%aerosol,CM%air,-1)
    ! +++ riming process +++
!!c    write(fid_alog,*) "in of rim rain"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_riming',3)
      if(CM%micexfg(4)==1.and.CM%flagp_s > 0 .and. CM%flagp_r > 0 ) then
         call coalescence( CM%solid_hydro, CM%rain, CM%air, CM%level_comp,CM%coll_level,CM%mes_rc &
                ,0,CM%imin_bk,CM%imax_bk,CM%jmin_bk,CM%jmax_bk,CM%bu_tmass,CM%bu_fd &
                ,ID,JD,KD &
! <<< 2014/10 T. Hashino added for KiD
                ,CM%adrpdrp,CM%drpdrp &
                ,CM%ahexdrp,CM%hexdrp &
                ,CM%abbcdrp,CM%bbcdrp &
                ,CM%acoldrp,CM%coldrp &
                ,CM%agp1drp,CM%gp1drp &
                ,CM%agp4drp,CM%gp4drp &
                ,CM%agp8drp,CM%gp8drp &
                ,dM_auto_rim,dM_accr_rim)
! >>> 2014/10 T. Hashino added for KiD
      !cptime(9)=cptime(9)+toc(s1,r1)
      endif
!      call PROF_rapend('amps_riming',3)
!!c    write(fid_alog,*) "out of rim rain"

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! melting and shedding process
!!c    write(fid_alog,*) "ms"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_meltshed',3)
      if(CM%micexfg(8)==1.and. CM%flagp_s > 0 .and. CM%flagp_r > 0) &
         call melting_shedding( CM%solid_hydro, CM%air, CM%level_comp, &
         CM%rain, CM%mes_rc)
      !cptime(18)=cptime(18)+toc(s1,r1)
!      call PROF_rapend('amps_meltshed',3)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,1,KD,ID,JD)
!!c    call check_tendency_ap(CM%solid_hydro,CM%aerosol,CM%air,1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!c    call ck_mlttend(CM%rain,CM%solid_hydro,ID,JD,KD)


      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! hydrodynamic breakup process
      ! NOTE:
      ! upgrade the density for hydrodynamic breakup process
      ! calculate the volume of a spheroid circumscribing the ice crystals, and
      ! identify aggregates and grauples.
!!c    write(fid_alog,*) "hb"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_hydrobreakup',3)
      if(CM%micexfg(11)==1.and. CM%flagp_r > 0 .and. CM%hbreak_r > 0 ) &
         call hydrodyn_breakup( CM%rain, CM%hbreak_r,CM%mes_rc)
      !cptime(10)=cptime(10)+toc(s1,r1)
      !call tic(s1,r1)
!!c    write(fid_alog,*) "hb ice"
      if(CM%micexfg(9)==1.and. CM%flagp_s > 0 .and. CM%hbreak_s > 0 ) &
           call hydrodyn_breakup( CM%solid_hydro, CM%hbreak_s,CM%mes_rc)
      !cptime(11)=cptime(11)+toc(s1,r1)
!      call PROF_rapend('amps_hydrobreakup',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfinuc1')

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Ice nucleation 1
!!c    write(fid_alog,*) "in1"
      !call tic(s1,r1)
!      call PROF_rapstart('amps_icenucleation',3)
      if(CM%micexfg(10)==1.and. CM%flagp_s > 0 ) &
        call Ice_Nucleation1( CM%solid_hydro, CM%rain, CM%aerosol,CM%air, CM%level_comp,CM%mes_rc,&
          CM%APSNAME,CM%nu_aps,CM%M_aps,&!,CM%phi_aps
          CM%ap_sig_cp,CM%ap_mean_cp,CM%CRIC_RN_IMM,CM%frac_dust, &!,CM%cdf_cp_180m0
          CM%micexfg(14),CM%micexfg(15),CM%micexfg(16),CM%micexfg(17),&
          ID,JD,KD, &
          CM%osm_nhs4,CM%osm_sdch)
      !cptime(13)=cptime(13)+toc(s1,r1)
!      call PROF_rapend('amps_icenucleation',3)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfrepair1')

      !call tic(s1,r1)
!      call PROF_rapstart("amps_repair",3)
      call repair(CM%level_comp, CM%air, CM%rain, CM%solid_hydro, CM%ncat_a, CM%aerosol,&
         CM%flagp_r, CM%flagp_s,CM%flagp_a,CM%mes_rc,& !,CM%act_type
         .false.,.true.,&
         iupdate_gr,iupdate_gs,iupdate_ga,&
         qtp,ID,JD,KD,'af_col',it_cl)
      !cptime(12)=cptime(12)+toc(s1,r1)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'bfupdategroup1')

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! update the group vars
      !call tic(s1,r1)
      call update_group_all(CM%rain,CM%solid_hydro,CM%aerosol,CM%air, &
                 CM%level_comp,CM%ncat_a,0,CM%mes_rc, &
                 qtp,RV_r,&
                 iupdate_gr,iupdate_gs,iupdate_ga,&
                 CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%eps_ap0,ID,JD,KD)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !cptime(19)=cptime(19)+toc(s1,r1)
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call cal_dmtend_scale(dmtendl,dcontendl,dbintendl,CM%L,CM)

!dbg    call lenchk1(CM%solid_hydro,id,jd,kd,'afupdategroup1')
!      call PROF_rapend("amps_repair",3)

      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! 2. vapor growth related processes
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! set time step size
      CM%rain%dt=CM%dt_vp
      CM%solid_hydro%dt=CM%dt_vp
      do i=1,CM%ncat_a
        CM%aerosol(i)%dt=CM%dt_vp
      enddo
      ! set the tendency mask for repair and update
      !                             1 2 3 4 5 6 7 8 9 10 11 12
            iupdate_gr(1:mxntend)=(/1,0,1,0,0,0,0,0,0, 0, 0, 0/)
            iupdate_gs(1:mxntend)=(/1,0,1,0,0,0,1,0,0, 0, 0, 0/)
      do i=1,CM%ncat_a
        if(i/=2) then
          iupdate_ga(1:mxntend,i)=(/1,0,1,0,0,1,0,0,0, 0, 0, 0/)
        elseif(i==2) then
          iupdate_ga(1:mxntend,i)=(/1,0,0,0,0,0,1,0,0, 0, 0, 0/)
        endif
      enddo

      vap_loop: do it_vp=1,CM%n_step_vp
        ! initialize all tendencies that are being updated.
        call ini_tendency(CM%rain,iupdate_gr,0)
        call ini_tendency(CM%solid_hydro,iupdate_gs,0)
        do i=1,CM%ncat_a
          call ini_tendency(CM%aerosol(i),iupdate_ga(1,i),0)
        enddo
        call ini_tendency_ag(CM%air,1)
        call ini_tendency_ag(CM%air,2)

        !call tic(s1,r1)
!        call PROF_rapstart('amps_update_diagnose',3)
        ! update hydrometeor flags
        call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                          CM%mes_rc,CM%ncat_a)
        ! re-analyze temperature field because T is changed by qr and qi.
        call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
            ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
            )

        ! update air group variables
        call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

        !cptime(1)=cptime(1)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_s > 0 ) &
          call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
        !cptime(2)=cptime(2)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_r > 0 ) &
          call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                        CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                        CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                        CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
        !cptime(4)=cptime(4)+toc(s1,r1)

        !call tic(s1,r1)
        if( CM%flagp_a /= 0 ) then
          do i=1,CM%ncat_a
            call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                          CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                          CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                          CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                          i, CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
          end do
        end if
        !cptime(5)=cptime(5)+toc(s1,r1)
        if( em /= 0 ) then
          LOG_ERROR("cal_micro_tendency",*) "diag_pq_vap > error!"
          call PRC_abort
        endif
!        call PROF_rapend('amps_update_diagnose',3)

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Vapor variable advancement and Activation of CCN
        !call tic(s1,r1)
!        call PROF_rapstart('amps_ccn',3)
        if(CM%act_type==2) then
          !
          ! cloud droplets activation from dry aerosols
          !   stable version
          call cal_aptact_var8_vec(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
              ,CM%flagp_a,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps,CM%CCNMAX&
              ,CM%snrml,ID,JD,KD)

        elseif(CM%act_type==3) then
          ! cloud droplets activation from dry aerosols with saturation adjustment
!org          call cal_aptact_var9(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
!org              ,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
!org              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
!org              ,ID,JD,KD)

        elseif(CM%act_type==4) then
          ! haze activation from dry aerosols with saturation adjustment
!          call cal_aptact_var11(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
!              ,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%mes_rc&
!              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
!              ,ID,JD,KD)

        else
          !
          ! haze activation from dry aerosols and
          ! ice nucleation process based on classical nucleation process approach.
          !
          ! requires use of small time steps (less than 1 sec).
          call cal_aptact_var8_kc04dep(CM%level_comp,CM%air,CM%ncat_a,CM%aerosol,CM%rain,CM%solid_hydro&
              ,CM%flagp_a,CM%flagp_r,CM%flagp_s,CM%micexfg(10),CM%micexfg(13),CM%micexfg(19)&
              ,CM%mes_rc&
              ,qtp,thil,CM%nu_aps,CM%phi_aps,CM%M_aps&
              ,CM%ap_sig_cp,CM%ap_mean_cp,CM%CCNMAX,CM%CRIC_RN_IMM,CM%frac_dust& !,CM%cdf_cp_180m0
              ,CM%snrml,ID,JD,KD)
        endif
        !cptime(14)=cptime(14)+toc(s1,r1)
!        call PROF_rapend('amps_ccn',3)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! vapor deposition process
        ! 1. vapor deposition process on drops first
!!c    write(fid_alog,*) "vr"
        !call tic(s1,r1)
!        call PROF_rapstart('amps_vapordep',3)
        if(CM%micexfg(6)==1.and. CM%flagp_r > 0 ) &
          call vapor_deposition( CM%rain,CM%aerosol, CM%air, CM%level_comp,CM%mes_rc &
                                ,ID,JD,KD &
                                ,CM%vigp,CM%rdsd,CM%ihabit_gm_random)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,-1,KD,ID,JD)
        !cptime(15)=cptime(15)+toc(s1,r1)
!!c    write(fid_alog,*) "vs"

        !call tic(s1,r1)
        if(CM%micexfg(7)==1.and. CM%flagp_s > 0 ) &
          call vapor_deposition(CM%solid_hydro,CM%aerosol,CM%air,CM%level_comp,CM%mes_rc &
                               ,ID,JD,KD &
                               ,CM%vigp,CM%rdsd,CM%ihabit_gm_random)
        !cptime(16)=cptime(16)+toc(s1,r1)
!!c    call check_tendency_ap(CM%rain,CM%aerosol,CM%air,0,KD,ID,JD)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! Ice nucleation 2
        !   This one include the mode of depositional growth
!!c    write(fid_alog,*) "in2"
        !call tic(s1,r1)
        if(CM%micexfg(10)==1.and. CM%flagp_s > 0 ) &
          call Ice_Nucleation2( CM%solid_hydro, CM%rain, CM%aerosol,CM%air, CM%level_comp,CM%mes_rc &
           ,CM%APSNAME,CM%nu_aps,CM%M_aps &
           ,CM%micexfg(13) &
           ,CM%flagp_a  &
           ,ID,JD,KD &
           ,CM%vigp,CM%rdsd,CM%ihabit_gm_random)
        !write(fid_alog,*) CM%cur_time, "ice2"
        !cptime(17)=cptime(17)+toc(s1,r1)
!        call PROF_rapend('amps_vapordep',3)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! repair routine
!!c    write(fid_alog,*) " front repair"
        !call tic(s1,r1)
!        call PROF_rapstart("amps_repair",3)
        call repair(CM%level_comp, CM%air, CM%rain, CM%solid_hydro, CM%ncat_a, CM%aerosol,&
             CM%flagp_r, CM%flagp_s,CM%flagp_a,CM%mes_rc,& !,CM%act_type
             .true.,.false.,&
             iupdate_gr,iupdate_gs,iupdate_ga,&
             qtp,ID,JD,KD,'af_vap',it_vp)
        !cptime(12)=cptime(12)+toc(s1,r1)
!!c    write(fid_alog,*) "out of repair"
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ! update the group vars
        !call tic(s1,r1)
        call update_group_all(CM%rain,CM%solid_hydro,CM%aerosol,CM%air, &
                 CM%level_comp,CM%ncat_a,1,CM%mes_rc, &
                 qtp,RV_r,&
                 iupdate_gr,iupdate_gs,iupdate_ga,&
                 CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%eps_ap0,ID,JD,KD)

        !cptime(19)=cptime(19)+toc(s1,r1)
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        call cal_dmtend_scale(dmtendl,dcontendl,dbintendl,CM%L,CM)
!        call PROF_rapend("amps_repair",3)

!!c        stop

      enddo vap_loop
    enddo col_loop1

    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Finally, diagnose the new variables
    ! for calculation of weighted terminal velocity.
    !call tic(s1,r1)
    ! update hydrometeor flags
    call update_mesrc(CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
                      CM%mes_rc,CM%ncat_a)
    ! re-analyze temperature field because T is changed by qr and qi.
    call diag_t(T_a_r,CM%air,CM%rain,CM%solid_hydro&
         ,thil,CM%mes_rc,CM%flagp_r,CM%flagp_s&
         )
    ! update air group variables
    call update_airgroup(CM%air,RV_r,T_a_r,CM%rdsd,CM%ihabit_gm_random)

    !cptime(1)=cptime(1)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_s > 0 ) &
      call diag_pq( CM%solid_hydro, CM%air, CM%level_comp, em, &
                    CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                    CM%eps_ap0(2),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                    CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
    !cptime(2)=cptime(2)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_r > 0 ) &
      call diag_pq( CM%rain, CM%air, CM%level_comp, em, &
                    CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                    CM%eps_ap0(1),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                    CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml)
    !cptime(4)=cptime(4)+toc(s1,r1)

    !call tic(s1,r1)
    if( CM%flagp_a /= 0 ) then
      do i=1,CM%ncat_a
        call diag_pq( CM%aerosol(i), CM%air, CM%level_comp, em, &
                      CM%mes_rc,ID,JD,KD,CM%rdsd,CM%ihabit_gm_random, &
                      CM%eps_ap0(i),CM%nu_aps,CM%phi_aps,CM%m_aps,&
                      CM%ap_sig_cp,CM%ap_mean_cp,CM%cdf_cp_0,CM%isnrml,&
                      i, CM%flagp_a,CM%ap_lnsig,CM%ap_mean)
      end do
    end if
    !cptime(5)=cptime(5)+toc(s1,r1)
    if( em /= 0 ) then
      LOG_ERROR("cal_micro_tendency",*) "diag_pq_vap > error!"
      call PRC_abort
    endif
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




    !cptime_mtd(1:20)=cptime_mtd(1:20)+cptime(1:20)

  end subroutine cal_micro_tendency



  subroutine return_output (CM,L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       NATYPE,NABIN,NACAT,RV,XC,XR,XS,XA,qtp)
    use class_Group, only: &
       update_modelvars_all, &
       cal_needgive
    use class_Thermo_Var, only: &
       renew_rv_var
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! calculate tendency of each variable and update variables (option).
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Cloud_Micro), intent(in) :: CM
    integer,intent(in) :: L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT
    real(MP_KIND)    :: RV(*),XC(*),XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),&
         XA(NATYPE,NABIN,NACAT,*),qtp(*)
    !integer :: ID(*),JD(*),KD(*)
    integer                     :: i

    call update_modelvars_all(XC,XR,XS,XA,RV, &
             L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
             CM%level_comp,CM%air,CM%rain,CM%solid_hydro,CM%aerosol,&
             CM%flagp_r,CM%flagp_s,CM%flagp_a)

!old    call cal_model_tendency_all(CM%level_comp,CM%out_type,CM%air,CM%mes_rc,&
!old         CM%rain,CM%solid_hydro,CM%aerosol,&
!old         L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,NATYPE,NABIN,NACAT,&
!old         qtp,RV,XC,XR,XS,XA,CM%flagp_r,CM%flagp_s,CM%flagp_a,CM%eps_ap0,ID,JD,KD)


  end subroutine return_output

  subroutine reality_check( CM,qtp,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use class_Group, only: &
       get_totalmass, &
       cal_needgive
    use class_Thermo_Var, only: &
       renew_rv_var
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! reality-check routine
    !
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Cloud_Micro), intent(inout)         :: CM

    real(PS),dimension(CM%L)         :: need_M, give_M
    real(PS),dimension(CM%L)         :: M_tot, M_v, M_tr, M_ts, M_tc
    ! necessary mass from vapor and mass given to vapor
    real(PS)                  :: ccloud
    real(PS), parameter       :: min_mcloud = (4.0_PS*PI/3.0_PS)*1.0e-9
    real(PS), parameter       :: min_rvapor = 1.0e-9 !  0.000001 g/kg
!!c    real(PS), parameter       :: min_rvapor = 0.0
    integer :: ID(*),JD(*),KD(*)
    real(MP_KIND) :: qtp(*)
    !
    integer :: i,j
    integer,dimension(CM%L) ::em,ierror

    ! +++ initialization +++
    CM%mes_rc = 1

    do j=1,CM%L
      ! first, check the total
!!c       write(fid_alog,*) " "
!!c       write(fid_alog,'(3I5)') ID(j),JD(j),KD(j)
!!c       write(fid_alog,'("mixing ratio of vapor",I5,ES15.6,&
!!c            /,"cloud drop",ES15.6,&
!!c            /,"rain",10ES15.6,/,"solid",5ES15.6)') &
!!c            j,CM%air(j)%rv,&
!!c            CM%cloud_drop%MS(1,j)%mass(1),&
!!c            (CM%rain%MS(i,j)%mass(1),i=1,CM%rain%N_BIN),&
!!c            (CM%solid_hydro%MS(i,j)%mass(1),i=1,CM%solid_hydro%N_BIN)
!!c       write(fid_alog,'("real",I5,10ES15.6)') j,(CM%solid_hydro%MS(i,j)%a_len,i=1,5),&
!!c            (CM%solid_hydro%MS(i,j)%c_len,i=1,5)


      M_v(j) = CM%air%TV(j)%rv*CM%air%TV(j)%den
      M_tr(j) = get_totalmass(CM%rain,j)
      M_ts(j) = get_totalmass(CM%solid_hydro,j)
    enddo
    do j=1,CM%L
      em(j)=0
      M_tc(j) = 0.0
      M_tot(j) = M_v(j) + M_tr(j) + M_ts(j)
    enddo

    do j=1,CM%L
      if( M_tot(j) <= 0.0 ) then
        CM%mes_rc(j) = 0
      else
        if(M_tr(j)>0.0_PS) then
          if(M_ts(j)>0.0_PS) then
            CM%mes_rc(j)=4
          else
            CM%solid_hydro%mark_cm(j)=3
            CM%mes_rc(j)=2
          end if
        else
          CM%rain%mark_cm(j)=3
          if(M_ts(j)>0.0_PS) then
            CM%mes_rc(j)=3
          else
            CM%solid_hydro%mark_cm(j)=3
            CM%mes_rc(j)=1
          end if
        end if
      endif
    enddo

    need_M = 0.0_PS
    give_M = 0.0_PS
    if( CM%flagp_s > 0 ) then
      do j=1,CM%L
        if(CM%mes_rc(j)>=2) then
          !
          ! in cases where hydrometeors exist in a grid box.
          !
          call cal_needgive( CM%solid_hydro, CM%rain, j,M_ts(j), M_tr(j), &
                   need_M(j), give_M(j))
        endif
      enddo
    endif
    if( CM%flagp_r > 0 ) then
      do j=1,CM%L
        if(CM%mes_rc(j)>=2) then
          !
          ! in cases where hydrometeors exist in a grid box.
          !
          call cal_needgive(CM%rain, CM%rain, j,M_tr(j), M_tc(j), &
                   need_M(j), give_M(j))
        endif
      enddo
    endif

    ierror(1:CM%L)=0
    do j=1,CM%L
      if(CM%mes_rc(j)>=2) then
        !
        ! in cases where hydrometeors exist in a grid box.
        !
        if( need_M(j) > 0.0 .or. give_M(j) > 0.0 ) then
          ierror(j)=1
        else
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! NOTE: if the vapor is not enough to maintain the minimum mixing ratio,
          !       the mass is produced, therefore, mass conservation breaks down.
!!c             M_v = max( M_v + give_M - need_M, min_rvapor*CM%air%TV(j)%den)
          M_v(j)=M_v(j)+give_M(j)-need_M(j)
        endif
      endif
    enddo
    if(any(ierror>0)) then
      do j=1,CM%L
        if(ierror(j)==1) then
           LOG_ERROR("reality_check",*) "need and give_M are not 0", j,need_M(j),give_M(j)
        endif
      enddo
      call PRC_abort
    endif

    do j=1,CM%L
      if(CM%mes_rc(j)>=2) then
        !
        ! in cases where hydrometeors exist in a grid box.
        !
        if(M_v(j)<min_rvapor*CM%air%TV(j)%den) then
!!c                write(fid_alog,'("vapor mass is less than realistic value.",3I5,10ES15.6)') &
!!c                            KD(j),ID(j),JD(j),&
!!c                            M_v(j),&
!!c                            min_rvapor*CM%air%TV(j)%den,M_tr(j),M_ts(j)

          M_v(j)=min_rvapor*CM%air%TV(j)%den
          em(j)=1
        endif
      endif
    enddo
    do j=1,CM%L
      if(CM%mes_rc(j)>=2) then
        !
        ! in cases where hydrometeors exist in a grid box.
        !
        if( give_M(j) /= 0.0_PS .or. need_M(j) /= 0.0_PS .or. em(j)/=0) then
          qtp(j)=qtp(j)+max(0.0_PS,M_v(j)/CM%air%TV(j)%den-CM%air%TV(j)%rv)
          call renew_rv_var(CM%air%TV(j), M_v(j)/CM%air%TV(j)%den)
        end if
        ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      end if
      CM%mt_r(j)=M_tr(j)/CM%air%TV(j)%den
      CM%mt_s(j)=M_ts(j)/CM%air%TV(j)%den
    end do

  end subroutine reality_check

  subroutine reality_check_scl( CM,qtp,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    use class_Group, only: &
       get_totalmass, &
       cal_needgive
    use class_Thermo_Var, only: &
       renew_rv_var
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! reality-check routine
    !
    !
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    type (Cloud_Micro), intent(inout)         :: CM

    real(PS)                       :: M_tot, M_v, M_tr, M_ts, M_tc
    ! necessary mass from vapor and mass given to vapor
    real(PS)                       :: need_M, give_M
    real(PS)                  :: ccloud
    real(PS), parameter       :: min_mcloud = (4.0_PS*PI/3.0_PS)*1.0e-9
    real(PS), parameter       :: min_rvapor = 1.0e-9 !  0.000001 g/kg
!!c    real(PS), parameter       :: min_rvapor = 0.0
    integer :: ID(*),JD(*),KD(*)
    real(PS) :: qtp(*)
    !
    integer                        :: i,j,em

    ! +++ initialization +++
    CM%mes_rc = 1

    do j=1,CM%L
       ! first, check the total
!!c       write(fid_alog,*) " "
!!c       write(fid_alog,'(3I5)') ID(j),JD(j),KD(j)
!!c       write(fid_alog,'("mixing ratio of vapor",I5,ES15.6,&
!!c            /,"cloud drop",ES15.6,&
!!c            /,"rain",10ES15.6,/,"solid",5ES15.6)') &
!!c            j,CM%air(j)%rv,&
!!c            CM%cloud_drop%MS(1,j)%mass(1),&
!!c            (CM%rain%MS(i,j)%mass(1),i=1,CM%rain%N_BIN),&
!!c            (CM%solid_hydro%MS(i,j)%mass(1),i=1,CM%solid_hydro%N_BIN)
!!c       write(fid_alog,'("real",I5,10ES15.6)') j,(CM%solid_hydro%MS(i,j)%a_len,i=1,5),&
!!c            (CM%solid_hydro%MS(i,j)%c_len,i=1,5)

       em=0

       M_v = CM%air%TV(j)%rv*CM%air%TV(j)%den
       M_tot = 0.0_PS
       M_tc = 0.0_PS

!!c       M_tr = CM%mt_r(j)*CM%air%TV(j)%den
!!c       M_ts = CM%mt_s(j)*CM%air%TV(j)%den
       M_tr = get_totalmass(CM%rain,j)
       M_ts = get_totalmass(CM%solid_hydro,j)

       M_tot = M_v + M_tr + M_ts
       if( M_tot <= 0.0_PS ) then
          CM%mes_rc(j) = 0
!!c       else if( M_v>0.0_PS ) then
       else
          if(M_tr>0.0_PS) then
             if(M_ts>0.0_PS) then
                CM%mes_rc(j)=4
             else
                CM%solid_hydro%mark_cm(j)=3
                CM%mes_rc(j)=2
             end if
          else
             CM%rain%mark_cm(j)=3
             if(M_ts>0.0_PS) then
                CM%mes_rc(j)=3
             else
                CM%solid_hydro%mark_cm(j)=3
                CM%mes_rc(j)=1
             end if
          end if
          if(CM%mes_rc(j)>=2) then
             ! in cases where hydrometeors exist in a grid box.
             need_M = 0.0_PS
             give_M = 0.0_PS
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if( CM%flagp_s > 0 ) then
                call cal_needgive( CM%solid_hydro, CM%rain, j,M_ts, M_tr, &
                     need_M, give_M)
!!!c          call check_length(CM%solid_hydro)
             end if
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if( CM%flagp_r > 0 ) then
                call cal_needgive(CM%rain, CM%rain, j,M_tr, M_tc, &
                     need_M, give_M)
             end if
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             if( need_M > 0.0 .or. give_M > 0.0 ) then
                LOG_ERROR("reality_check_scl",*) "need and give_M are not 0", j,need_M,give_M
                call PRC_abort
             end if

             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             ! NOTE: if the vapor is not enough to maintain the minimum mixing ratio,
             !       the mass is produced, therefore, mass conservation breaks down.
!!c             M_v = max( M_v + give_M - need_M, min_rvapor*CM%air%TV(j)%den)
             M_v=M_v+give_M-need_M
             if(M_v<min_rvapor*CM%air%TV(j)%den) then
!!c                write(fid_alog,'("vapor mass is less than realistic value.",3I5,10ES15.6)') &
!!c                            KD(j),ID(j),JD(j),&
!!c                            M_v,&
!!c                            min_rvapor*CM%air%TV(j)%den,M_tr,M_ts

                M_v=min_rvapor*CM%air%TV(j)%den
                em=1
             endif
             if( give_M /= 0.0_PS .or. need_M /= 0.0_PS .or. em/=0) then
                qtp(j)=qtp(j)+max(0.0_PS,M_v/CM%air%TV(j)%den-CM%air%TV(j)%rv)
                call renew_rv_var(CM%air%TV(j), M_v/CM%air%TV(j)%den)
             end if
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          end if

       end if
       CM%mt_r(j)=M_tr/CM%air%TV(j)%den
       CM%mt_s(j)=M_ts/CM%air%TV(j)%den

    end do

  end subroutine reality_check_scl


  subroutine print_cloud_tendency(CM,ID,JD,KD,iproc,istrt)
    use class_Group, only: &
       print_tendency, &
       print_tendency_ap
    type (Cloud_Micro), intent(in)         :: CM
    integer :: ID(*),JD(*),KD(*)
    integer,intent(in) :: iproc,istrt
!!c    return

    if( mod(int(CM%cur_time*100), CM%T_print_period*100) == 0 ) then

!!c       if( CM%flagp_c > 0 ) &
!!c            call print_tendency( CM%cloud_drop, CM%cur_time,ID,JD,KD,iproc,istrt)

       if( CM%flagp_r > 0 ) &
            call print_tendency( CM%rain, CM%cur_time, CM%output_format,ID,JD,KD,iproc,istrt)

       if( CM%flagp_a /= 0 ) &
            call print_tendency_ap( CM%aerosol, CM%cur_time, CM%output_format,ID,JD,KD,iproc,istrt,CM%ncat_a )

       if( CM%flagp_s > 0 ) &
            call print_tendency( CM%solid_hydro, CM%cur_time, CM%output_format,ID,JD,KD,iproc,istrt)
    end if
  end subroutine print_cloud_tendency

  subroutine update_terminal_vel(CM,L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       vtcloud,XR,XS)
    use mod_amps_core, only: &
       w_terminal_vel
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in) :: L
    integer,intent(in) :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(MP_KIND) :: vtcloud(*)
    real(MP_KIND) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*)
    ! mass-weighted terminal velocity
    ! 1st argument
    !    1: mass weighted terminal velocity
    !    2: con weighted terminal velocity
    ! 2nd argument
    !    bin number
    ! 3rd argument
    !    grid number
!tmp    real, pointer, dimension(:,:,:)  :: VTR,VTS
    real(PS), dimension(2,mxnbin,LMAX)  :: VTR,VTS
    real(PS) :: artificial_tv
    integer,parameter :: nwv=2
    integer :: i,k,n,in,item(2)

    if( CM%flagp_r > 0 ) &
         call w_terminal_vel(CM%rain, CM%air, VTR)

    if( CM%flagp_s > 0 ) &
         call w_terminal_vel(CM%solid_hydro, CM%air, VTS)


    ! The terminal veloctiy is assigned to the last argument of parameter.
    ! The unit is also fixed to MKS.
    !!!! NOTE: terminal velocity has to be negative in the dynamic model!!!
    !!!! now the velocity calculated in here is positive is downward.

    do n=1,CM%L
      vtcloud(n)=-VTR(1,1,n)*1.0e-2
    enddo

    artificial_tv = 0.95

!CDIR NODEP
    do in=1,NRBIN*CM%L
      n=(in-1)/NRBIN+1
      i=in-(n-1)*NRBIN

      k=1
      if(XR(rmt_q,i,1,n)>0.0) then
        XR(NRTYPE-nwv+k,i,1,n)=-VTR(k,i,n)*1.0e-2
      else
        XR(NRTYPE-nwv+k,i,1,n)=0.0
      end if

      k=2
      if(XR(rcon_q,i,1,n)>0.0) then
        XR(NRTYPE-nwv+k,i,1,n)=-VTR(k,i,n)*1.0e-2
      else
        XR(NRTYPE-nwv+k,i,1,n)=0.0
      end if

    end do

!CDIR NODEP
    do in=1,NSBIN*CM%L
      n=(in-1)/NSBIN+1
      i=in-(n-1)*NSBIN

      k=1
      if(XS(imt_q,i,1,n)>0.0) then
        XS(NSTYPE-nwv+k,i,1,n)=-VTS(k,i,n)*1.0e-2
      else
        XS(NSTYPE-nwv+k,i,1,n)=0.0
      end if

      k=2
      if(XS(icon_q,i,1,n)>0.0) then
        XS(NSTYPE-nwv+k,i,1,n)=-VTS(k,i,n)*1.0e-2
      else
        XS(NSTYPE-nwv+k,i,1,n)=0.0
      end if
    end do

  end subroutine update_terminal_vel

  subroutine chk_ice(CM, L,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       vtcloud, XR, XS, loc)
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in) :: L
    integer,intent(in) :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: vtcloud(*)
    real(PS) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*)
    integer :: i,n
    character *(*) loc

    if ( debug ) then
       write(fid_alog,*) loc
       do n=1,L
          do i=1,NSBIN
             if(XS(1,i,NSCAT,n)-XS(11,i,NSCAT,n)-XS(10,i,NSCAT,n)>XS(1,i,NSCAT,n)/100.0) then
                write(fid_alog,'("ch_qi 1",2I5,3ES15.6)') i,n,XS(1,i,NSCAT,n),XS(10,i,NSCAT,n),XS(11,i,NSCAT,n)
             end if
             if(XS(11,i,NSCAT,n)<0.0) then
                write(fid_alog,'("ch_qi 2",2I5,3ES15.6)') i,n,XS(1,i,NSCAT,n),XS(10,i,NSCAT,n),XS(11,i,NSCAT,n)
             end if
          end do
       end do
    end if
  end subroutine chk_ice

  subroutine check_water_apmass_scl( &
       CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XA,XR,XS,qtp,RV,from,iswitch,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    type (Cloud_Micro), intent(inout)         :: CM
    integer,intent(in)  :: NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: qtp(*),&
         RV(*),XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),XA(NATYPE,NABIN,NACAT,*)

    integer,intent(in) :: iswitch
    integer :: ID(*),JD(*),KD(*)

    integer :: n,i,j, k,iter
    character (len=*) :: from

    real(PS) :: maps_r,maps_s,mapt_r2,maps_r2,mapt_s2,maps_s2,mapt_a2,sum_r,sum_s,diff
    integer :: mark_r(NRBIN),mark_s(NSBIN),ick_r(NRBIN),ick_s(NSBIN)
    integer :: inor,inos
    real(PS) :: frac,frac_s,ov
    real(PS),parameter :: mlmt=1.0e-30,ok_frac=1.0e-2

    !     The minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    real(PS),parameter :: n_lmt_ap=1.0e-5
    !     The minimum radius possible for large particles
    real(PS),parameter :: r3_lmt=1.0e-18
    real(PS),parameter :: mx_aprat=1.0
!    real(PS),parameter :: mx_aprat=0.99
    !     4.0*pi/3.0
    real(PS),parameter :: coef3=4.18879020478639

    real(PS), parameter              :: Rdvchiarui = 287.04_RP/461.50_RP ! origin 0.622_PS

    real(PS) :: m_lmt_ap,T0,T1,mtr0,mts0,massagg

    if(iswitch==0) then
!!c       write(fid_alog,*) "here1?"
       do n=1,CM%L
          ! +++ loop over grids +++

          ! +++ for aerosols +++
          do i = 1, nacat
             CM%mapt_ac(i,n)=CM%aerosol(i)%MS(1,n)%mass(amt)/CM%air%TV(n)%den
          enddo

          CM%mt_r(n)=0.0_PS
          CM%mapt_r(n)=0.0_PS
          maps_r=0.0_PS

          CM%mt_s(n)=0.0_PS
          CM%mapt_s(n)=0.0_PS
          maps_s=0.0_PS

          ! +++ for rain +++
          if(CM%rain%mark_cm(n)/=3) then
             do i=1,NRBIN
                CM%mt_r(n)=CM%mt_r(n)&
                     +max(0.0_PS,CM%rain%MS(i,n)%mass(rmt)-CM%rain%MS(i,n)%mass(rmat))
                CM%mapt_r(n)=CM%mapt_r(n)+CM%rain%MS(i,n)%mass(rmat)
                maps_r=maps_r+CM%rain%MS(i,n)%mass(rmas)
             end do
             CM%mt_r(n)=CM%mt_r(n)/CM%air%TV(n)%den
             CM%mapt_r(n)=CM%mapt_r(n)/CM%air%TV(n)%den
             maps_r=maps_r/CM%air%TV(n)%den
          end if

          ! +++ for ice +++
          if(CM%solid_hydro%mark_cm(n)/=3) then
             do i=1,NSBIN
                CM%mt_s(n)=CM%mt_s(n)&
                     +max(0.0_PS,CM%solid_hydro%MS(i,n)%mass(imt)-CM%solid_hydro%MS(i,n)%mass(imat))
                CM%mapt_s(n)=CM%mapt_s(n)+CM%solid_hydro%MS(i,n)%mass(imat)
                maps_s=maps_s+CM%solid_hydro%MS(i,n)%mass(imas)
             end do
             CM%mt_s(n)=CM%mt_s(n)/CM%air%TV(n)%den
             CM%mapt_s(n)=CM%mapt_s(n)/CM%air%TV(n)%den
             maps_s=maps_s/CM%air%TV(n)%den
          end if

!!c          write(fid_alog,*) "From ",from
!!c          write(fid_alog,*) "mtr, msr, mti, msi"
!!c          write(fid_alog,'(4ES15.6)') CM%mapt_r(n),maps_r,CM%mapt_s(n),maps_s
!!c          write(fid_alog,*) "mtap1, msap1, mtap2, msap2"
!!c          write(fid_alog,'(4ES15.6)') CM%aerosol(1)%MS(1,n)%mass(amt)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(1)%MS(1,n)%mass(ams)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(2)%MS(1,n)%mass(amt)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(2)%MS(1,n)%mass(ams)/CM%air%TV(n)%den
!!c          write(fid_alog,*) "total aerosol mass"
!!c          write(fid_alog,'(ES15.6)') CM%mapt_r(n)+CM%mapt_s(n)+CM%mapt_a(n)

!!c          write(fid_alog,'("SHIPS 1:k,es,qs,qv,t,p,den,qex,qtp",I5,8ES15.6)') &
!!c                  n,CM%air%TV(n)%e_sat(1),CM%air%TV(n)%rv_sat(1),CM%air%TV(n)%rv,CM%air%TV(n)%T,CM%air%TV(n)%P,CM%air%TV(n)%den,qex(n),qtp(n)

!!c          write(fid_alog,'("SHIPS 1:k,i,j,sv1,sv2,qv,t,p,den,qtp,mtr,mts",3I5,20ES15.6)') &
!!c                  KD(n),ID(n),JD(n),CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%rv,CM%air%TV(n)%T,&
!!c                  CM%air%TV(n)%P,CM%air%TV(n)%den,qtp(n),CM%mt_r(n),CM%mt_s(n)

          ! +++ check total water +++
          CM%air%TV(n)%rv=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
!!c          CM%air%TV(n)%rv=max(0.0_PS,qtp(n)-qex(n)-CM%mt_r(n)-CM%mt_s(n))
          CM%air%TV(n)%e=CM%air%TV(n)%P*CM%air%TV(n)%rv/(Rdvchiarui+CM%air%TV(n)%rv)
          do j=1,2
             CM%air%TV(n)%s_v(j)=CM%air%TV(n)%e/CM%air%TV(n)%e_sat(j)-1.0_PS
             CM%air%TV(n)%s_v_n(j)=CM%air%TV(n)%s_v(j)
          end do
!!c         if(CM%mt_r(n)+CM%mt_s(n)+CM%air%TV(n)%rv<qtp(n)) then
!!c             write(fid_alog,*) "check_water_apmass > total is less than qtp at time 0"
!!c             write(fid_alog,*) "n,mtr,mts,rv,qtp",n,CM%mt_r(n),CM%mt_s(n),CM%air%TV(n)%rv,qtp(n)
!!c          end if

!!c          write(fid_alog,'("SHIPS 2:k,i,j,sv1,sv2,qv,t,p,den,qtp,mtr,mts",3I5,20ES15.6)') &
!!c                  KD(n),ID(n),JD(n),CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%rv,CM%air%TV(n)%T,&
!!c                  CM%air%TV(n)%P,CM%air%TV(n)%den,qtp(n),CM%mt_r(n),CM%mt_s(n)

       end do
    elseif(iswitch==1) then
!!c       write(fid_alog,*) "flagp_a is",CM%flagp_a
       do n=1,CM%L

          mtr0=CM%mt_r(n)
          mts0=CM%mt_s(n)

          ! +++ for aerosols +++
          mapt_a2=0.0_PS
          do i=1,nacat
             mapt_a2=mapt_a2+XA(1,1,i,n)
          end do



          CM%mt_r(n)=0.0_PS
          mapt_r2=0.0_PS
          maps_r2=0.0_PS
          mark_r=0
          sum_r=0.0_PS
          ick_r=0
          inor=0

          CM%mt_s(n)=0.0_PS
          mapt_s2=0.0_PS
          maps_s2=0.0_PS
          mark_s=0
          sum_s=0.0_PS
          inos=0
          ick_s=0

          if( CM%flagp_a/=0 ) then
             ! +++ for rain +++
             do i=1,NRBIN
                CM%mt_r(n)=CM%mt_r(n)+XR(rmt_q,i,NRCAT,n)
                mapt_r2=mapt_r2+XR(rmat_q,i,NRCAT,n)
                maps_r2=maps_r2+XR(rmas_q,i,NRCAT,n)
                if(XR(rmt_q,i,NRCAT,n)>0.0_PS.and.XR(rmat_q,i,NRCAT,n)<=0.0_PS) then
                   inor=inor+1
                   ick_r(i)=1
                   mark_r(inor)=i
                   sum_r=sum_r+XR(rmt_q,i,NRCAT,n)
                   !write(fid_alog,*) "xr1,xr3",XR(rmt_q,i,NRCAT,n),XR(rmat_q,i,NRCAT,n)
                end if
             end do

             do i=1,NSBIN
                CM%mt_s(n)=CM%mt_s(n)+XS(imt_q,i,NSCAT,n)
                mapt_s2=mapt_s2+XS(imat_q,i,NSCAT,n)
                maps_s2=maps_s2+XS(imas_q,i,NSCAT,n)
                if(XS(imt_q,i,NSCAT,n)>0.0_PS.and.XS(imat_q,i,NSCAT,n)<=0.0_PS) then
                   inos=inos+1
                   ick_s(i)=1
                   mark_s(inos)=i
                   sum_s=sum_s+XS(imt_q,i,NSCAT,n)
                   if (debug) write(fid_alog,*) "xs1,xs13",XS(imt_q,i,NSCAT,n),XS(imat_q,i,NSCAT,n)
                end if
             end do

             if(inos>0.or.inor>0) then
                ! non positive total aerosol masses should be fixed in
                ! subroutine "cal_model_tendency_all"
                LOG_ERROR("check_water_apmass_scl",*) "check_total_apmass > something is wrong, inor or inos is not 0"
                LOG_ERROR_CONT(*) inor,inos
                call PRC_abort
             end if

             if(mapt_a2+mapt_r2+mapt_s2<mlmt) then
                LOG_ERROR("check_water_apmass_scl",*) "mass of AP is too low"
                call PRC_abort
             else
                T0=0.0_PS
                do i=1,NACAT
                   T0=T0+CM%mapt_ac(i,n)
                end do
                T0=T0+CM%mapt_r(n)+CM%mapt_s(n)
                T1=mapt_a2+mapt_r2+mapt_s2
                frac=T0/T1
             end if
             if(debug .and. frac>=5.0) then
                write(fid_alog,'("water_apmass>frac is large",3I5,ES15.6)') KD(n),ID(n),JD(n),frac
                !write(fid_alog,'("0:T0,mta(cat),mtar,mtai",10ES15.6)') T0,CM%mapt_ac(1:NACAT,n)&
                !   ,CM%mapt_r(n),CM%mapt_s(n)
                !write(fid_alog,'("1:T1,mta,mtar,mtai",10ES15.6)') T1,mapt_a2,mapt_r2,mapt_s2
                !write(fid_alog,'("sv1,sv2,W",3ES15.6)') CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%W
                !write(fid_alog,'("dsv1,dsv2",3ES15.6)') CM%air%TV(n)%s_v_n(1),CM%air%TV(n)%s_v_n(2)
                !write(fid_alog,'("dmdt1 r",22ES15.6)') (CM%rain%MS(i,n)%dmassdt(1,1),i=1,CM%rain%N_BIN)
                !write(fid_alog,'("dmdt1 s",22ES15.6)') (CM%solid_hydro%MS(i,n)%dmassdt(1,1),i=1,CM%solid_hydro%N_BIN)
                !write(fid_alog,'("sum dmdt at r",22ES15.6)') sum(CM%rain%MS(:,n)%dmassdt(rmat,1))
                !write(fid_alog,'("sum dmdt at s",22ES15.6)') sum(CM%solid_hydro%MS(:,n)%dmassdt(imat,1))
                !write(fid_alog,'("sum dmdt at a",22ES15.6)') sum(CM%aerosol(1)%MS(:,n)%dmassdt(amt,1)),&
                !       sum(CM%aerosol(2)%MS(:,n)%dmassdt(amt,1)),sum(CM%aerosol(3)%MS(:,n)%dmassdt(amt,1))
                !Component 'CM' to right of part references with nonzero rank must not have ALLOCATABLE attribute.
             end if

             if(abs(abs(frac)-1.0_PS)>=ok_frac) then
                if(frac>1.0_PS) then
!!c                   write(fid_alog,'("water_apmass>apmass lost:grid,frac,T0,T1,ap1m,ap2m",I5,9ES15.6)')&
!!c                          n,frac,T0,T1,XA(1,1,1,n),XA(1,1,2,n)
                else
                   if(CM%flagp_a==-1.or.CM%flagp_a==-4) then
                     ! predict only CCN
                      i=1
                      XA(amt_q,1,i,n)=frac*XA(amt_q,1,i,n)
                      XA(acon_q,1,i,n)=XA(acon_q,1,i,n)*frac
                      XA(ams_q,1,i,n)=XA(ams_q,1,i,n)*frac
                      do i=3,nacat
                         XA(amt_q,1,i,n)=frac*XA(amt_q,1,i,n)
                         XA(acon_q,1,i,n)=XA(acon_q,1,i,n)*frac
                         XA(ams_q,1,i,n)=XA(ams_q,1,i,n)*frac
                      end do

                   elseif(CM%flagp_a==-2.or.CM%flagp_a==-5) then
                     ! predict only IN
                      i=2
                      XA(amt_q,1,i,n)=frac*XA(amt_q,1,i,n)
                      XA(acon_q,1,i,n)=XA(acon_q,1,i,n)*frac
                      XA(ams_q,1,i,n)=XA(ams_q,1,i,n)*frac
                   elseif(CM%flagp_a==-3.or.CM%flagp_a==-6) then
                     ! do not pass anything to dynamics

                   else
                      do i=1,nacat
                         XA(amt_q,1,i,n)=frac*XA(amt_q,1,i,n)
                         XA(acon_q,1,i,n)=XA(acon_q,1,i,n)*frac
                         XA(ams_q,1,i,n)=XA(ams_q,1,i,n)*frac
                      end do
                   endif

                   mapt_a2=0.0_PS
                   do i=1,nacat
                      mapt_a2=mapt_a2+XA(amt_q,1,i,n)
                   end do

                   mapt_r2=0.0_PS
                   mapt_s2=0.0_PS
                   do i=1,NRBIN
                      XR(rmat_q,i,NRCAT,n)=min(frac*XR(rmat_q,i,NRCAT,n),mx_aprat*XR(rmt_q,i,NRCAT,n))
                      XR(rmas_q,i,NRCAT,n)=min(frac*XR(rmas_q,i,NRCAT,n),XR(rmat_q,i,NRCAT,n))
                      mapt_r2=mapt_r2+XR(rmat_q,i,NRCAT,n)
                   end do
                   do i=1,NSBIN
                      XS(imat_q,i,NSCAT,n)=min(frac*XS(imat_q,i,NSCAT,n),mx_aprat*XS(imt_q,i,NSCAT,n))
                      XS(imas_q,i,NSCAT,n)=min(frac*XS(imas_q,i,NSCAT,n),XS(imat_q,i,NSCAT,n))
                      mapt_s2=mapt_s2+XS(imat_q,i,NSCAT,n)
                   end do

                end if

             end if
          end if

!!c          return

          CM%mt_r(n)=0.0_PS
          CM%mt_s(n)=0.0_PS
          ! +++ for rain +++
          do i=1,NRBIN
             CM%mt_r(n)=CM%mt_r(n)+max(0.0_PS,XR(rmt_q,i,NRCAT,n)-XR(rmat_q,i,NRCAT,n))
          end do
          do i=1,NSBIN
             CM%mt_s(n)=CM%mt_s(n)+max(0.0_PS,XS(imt_q,i,NSCAT,n)-XS(imat_q,i,NSCAT,n))
          end do

!!c          RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
!!c          cycle

          ! check total water
          if(CM%mt_r(n)+CM%mt_s(n)>mlmt) then
             if(qtp(n)-RV(n)<0.0_PS) then
                RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
             end if
             frac=(qtp(n)-RV(n))/(CM%mt_r(n)+CM%mt_s(n))

             ! **************************************************************************************
             !                ************* NOTE **************
             ! frac can be far from 1 in cases where mass of liquid and solid hydrometeors
             ! are much smaller than vapor mixing ratio (less than one million-th).
             ! So that is a numerical limitation, instead of errors in microphysical precitions.
             ! **************************************************************************************

!!c             write(fid_alog,'("total water> n frac,qtp,rv,mtr,mts,W:",I4,7ES15.6)')&
!!c                        n,frac,qtp(n),RV(n),CM%mt_r(n), CM%mt_s(n),CM%air%TV(n)%W
             if(abs(abs(frac)-1.0_PS)>=ok_frac) then
!!c                if(frac>2.0) then
!!c                   write(fid_alog,*) "ERROR"
!!c                   write(fid_alog,'("check_total water> k,i,j,frac,qtp,rv,mtr,mts,W:",3I5,7ES15.6)')&
!!c                        KD(n),ID(n),JD(n),&
!!c                        frac,qtp(n),RV(n),CM%mt_r(n), CM%mt_s(n),CM%air%TV(n)%W
!!c                   write(fid_alog,*) "W,sv1",CM%air%TV(n)%W,CM%air%TV(n)%s_v(1)
!!c                   write(fid_alog,*) "rv,rv1",CM%air%TV(n)%rv,CM%air%TV(n)%rv_sat(1)
!!c                   write(fid_alog,*) "1t",(CM%rain%MS(i,n)%dmassdt(1,1),i=1,CM%rain%N_BIN)
!!c                   write(fid_alog,*) "3t",(CM%rain%MS(i,n)%dmassdt(1,3),i=1,CM%rain%N_BIN)
!!c                   write(fid_alog,*) "1ap",(CM%rain%MS(i,n)%dmassdt(2,1),i=1,CM%rain%N_BIN)
!!c                   write(fid_alog,*) "3ap",(CM%rain%MS(i,n)%dmassdt(2,3),i=1,CM%rain%N_BIN)
!!c                   write(fid_alog,*) "binmass_r",(CM%rain%MS(i,n)%mass(rmt),i=1,CM%rain%N_BIN)
!!c                   write(fid_alog,*) "binmass_a",(CM%aerosol(1)%MS(i,n)%mass(amt),i=1,CM%aerosol(1)%N_BIN)
!!c                end if

                if(frac>1.0_PS.or.frac==0.0_PS) then
                   ! this can create high supersaturation over water if mt_r and mt_s are not well predicted.
                   RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
                else
                   ! +++ for liquid +++
                   do i=1,NRBIN
                      do j=1,NRTYPE
                         if(XR(j,i,NRCAT,n)>0.0_PS) then
                            XR(j,i,NRCAT,n)=XR(j,i,NRCAT,n)*frac
                         end if
                      end do
                   end do
                   ! +++ for ice +++
                   do i=1,NSBIN
                      do j=1,NSTYPE
                         if(XS(j,i,NSCAT,n)>0.0_PS) then
                            XS(j,i,NSCAT,n)=XS(j,i,NSCAT,n)*frac
                         end if
                      end do
                   end do
                end if
             end if
          else
             RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
          end if

!!c          do i=1,NSBIN
!!c             massagg=XS(1,i,NSCAT,n)-XS(10,i,NSCAT,n)-XS(11,i,NSCAT,n)
!!c             if(XS(11,i,NSCAT,n)<massagg) then
!!c                write(fid_alog,*) "water_apmass>agg forming",massagg,XS(1,i,NSCAT,n),XS(10,i,NSCAT,n),XS(11,i,NSCAT,n)
!!c                write(fid_alog,*) "i,n",i,n
!!c             end if
!!c          end do
       end do
    end if

  end subroutine check_water_apmass_scl

  subroutine check_water_apmass( &
       CM,NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XA,XR,XS,qtp,RV,from,iswitch,ID,JD,KD)
    use scale_prc, only: &
       PRC_abort
    type (Cloud_Micro), intent(inout)         :: CM
    integer,intent(in)  :: NATYPE,NABIN,NACAT,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(MP_KIND) :: qtp(*),&
         RV(*),XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),XA(NATYPE,NABIN,NACAT,*)

    integer,intent(in) :: iswitch
    integer :: ID(*),JD(*),KD(*)

    integer :: n,i,j,k,ijn,in,kn,iter
    character (len=*) :: from

    real(PS) :: maps_r,maps_s,maps_r2,maps_s2,sum_r,sum_s,diff
    real(PS),dimension(CM%L) :: mapt_r2,mapt_s2,mapt_a2
    integer,dimension(CM%L) :: inor,inos
    real(PS) :: frac_s,ov,rv_n,tmass1
    real(PS),parameter :: mlmt=1.0e-30,ok_frac=1.0e-2,trunc=1.0e-5

    !     The minimum possible background concentration for large particle
    !     is assumed to be 1.0e-5 cm^-3
    real(PS),parameter :: n_lmt_ap=1.0e-5
    !     The minimum radius possible for large particles
    real(PS),parameter :: r3_lmt=1.0e-18
    real(PS),parameter :: mx_aprat=1.0
!    real(PS),parameter :: mx_aprat=0.99
    !     4.0*pi/3.0
    real(PS),parameter :: coef3=4.18879020478639

    real(PS), parameter              :: Rdvchiarui = 287.04_RP/461.50_RP ! origin 0.622_PS

    real(PS) :: m_lmt_ap,mtr0,mts0,massagg
    real(PS),dimension(CM%L) :: T0,T1,frac
    integer,dimension(CM%L) :: ierror,icond1

    if(iswitch==0) then
!!c       write(fid_alog,*) "here1?"

      ! +++ for aerosols +++
      do kn=1,nacat*CM%L
        k=(kn-1)/CM%L+1
        n=kn-(k-1)*CM%L

        CM%mapt_ac(k,n)=CM%aerosol(k)%MS(1,n)%mass(amt)/CM%air%TV(n)%den
      enddo

      do n=1,CM%L
        CM%mt_r(n)=0.0_PS
        CM%mapt_r(n)=0.0_PS
        CM%mt_s(n)=0.0_PS
        CM%mapt_s(n)=0.0_PS
      enddo

      do i=1,NRBIN
        do n=1,CM%L
          ! +++ for rain +++
          if(CM%rain%mark_cm(n)/=3) then
            CM%mt_r(n)=CM%mt_r(n)+ &
                max(0.0_PS,CM%rain%MS(i,n)%mass(rmt)-CM%rain%MS(i,n)%mass(rmat))
            CM%mapt_r(n)=CM%mapt_r(n)+ &
                CM%rain%MS(i,n)%mass(rmat)
          endif
        enddo
      enddo

      do i=1,NSBIN
        do n=1,CM%L
          ! +++ for ice +++
          if(CM%solid_hydro%mark_cm(n)/=3) then
            CM%mt_s(n)=CM%mt_s(n)+ &
              max(0.0_PS,CM%solid_hydro%MS(i,n)%mass(imt)-CM%solid_hydro%MS(i,n)%mass(imat))
            CM%mapt_s(n)=CM%mapt_s(n)+ &
              CM%solid_hydro%MS(i,n)%mass(imat)
          endif
        enddo
      enddo

      do n=1,CM%L
        CM%mt_r(n)=CM%mt_r(n)/CM%air%TV(n)%den
        CM%mapt_r(n)=CM%mapt_r(n)/CM%air%TV(n)%den

        CM%mt_s(n)=CM%mt_s(n)/CM%air%TV(n)%den
        CM%mapt_s(n)=CM%mapt_s(n)/CM%air%TV(n)%den
      end do

      do n=1,CM%L
!!c          write(fid_alog,*) "From ",from
!!c          write(fid_alog,*) "mtr, msr, mti, msi"
!!c          write(fid_alog,'(4ES15.6)') CM%mapt_r(n),maps_r,CM%mapt_s(n),maps_s
!!c          write(fid_alog,*) "mtap1, msap1, mtap2, msap2"
!!c          write(fid_alog,'(4ES15.6)') CM%aerosol(1)%MS(1,n)%mass(amt)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(1)%MS(1,n)%mass(ams)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(2)%MS(1,n)%mass(amt)/CM%air%TV(n)%den&
!!c               ,CM%aerosol(2)%MS(1,n)%mass(ams)/CM%air%TV(n)%den
!!c          write(fid_alog,*) "total aerosol mass"
!!c          write(fid_alog,'(ES15.6)') CM%mapt_r(n)+CM%mapt_s(n)+CM%mapt_a(n)

!!c          write(fid_alog,'("SHIPS 1:k,es,qs,qv,t,p,den,qex,qtp",I5,8ES15.6)') &
!!c                  n,CM%air%TV(n)%e_sat(1),CM%air%TV(n)%rv_sat(1),CM%air%TV(n)%rv,CM%air%TV(n)%T,CM%air%TV(n)%P,CM%air%TV(n)%den,qex(n),qtp(n)

!!c          write(fid_alog,*) "SHIPS 1:k,i,j,sv1,sv2,qv,t,p,den,qtp,mtr,mts",&
!!c                  KD(n),ID(n),JD(n),CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%rv,CM%air%TV(n)%T,&
!!c                  CM%air%TV(n)%P,CM%air%TV(n)%den,qtp(n),CM%mt_r(n),CM%mt_s(n)

        ! +++ check total water +++
        rv_n=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
        CM%air%TV(n)%rv=rv_n
        CM%air%TV(n)%e=CM%air%TV(n)%P*rv_n/(Rdvchiarui+rv_n)

        CM%air%TV(n)%s_v(1)=CM%air%TV(n)%e/CM%air%TV(n)%e_sat(1)-1.0_PS
        CM%air%TV(n)%s_v_n(1)=CM%air%TV(n)%s_v(1)
        CM%air%TV(n)%s_v(2)=CM%air%TV(n)%e/CM%air%TV(n)%e_sat(2)-1.0_PS
        CM%air%TV(n)%s_v_n(2)=CM%air%TV(n)%s_v(2)
!!c         if(CM%mt_r(n)+CM%mt_s(n)+CM%air%TV(n)%rv<qtp(n)) then
!!c             write(fid_alog,*) "check_water_apmass > total is less than qtp at time 0"
!!c             write(fid_alog,*) "n,mtr,mts,rv,qtp",n,CM%mt_r(n),CM%mt_s(n),CM%air%TV(n)%rv,qtp(n)
!!c          end if

!!c          write(fid_alog,*) "SHIPS 2:k,i,j,sv1,sv2,qv,t,p,den,qtp,mtr,mts",&
!!c                  KD(n),ID(n),JD(n),CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%rv,CM%air%TV(n)%T,&
!!c                  CM%air%TV(n)%P,CM%air%TV(n)%den,qtp(n),CM%mt_r(n),CM%mt_s(n)

      end do
    elseif(iswitch==1) then
!!c       write(fid_alog,*) "flagp_a is",CM%flagp_a

      ! +++ for aerosols +++
      mapt_a2=0.0_PS
      do i=1,nacat
        do n=1,CM%L
          mapt_a2(n)=mapt_a2(n)+CM%aerosol(i)%MS(1,n)%mass(amt)/CM%air%tv(n)%den
        end do
      enddo

      if( CM%flagp_a>0 ) then
        ! +++ for rain +++
        do n=1,CM%L
          CM%mt_r(n)=0.0_PS
          mapt_r2(n)=0.0_PS
!!!          maps_r2(n)=0.0_PS
          CM%mt_s(n)=0.0_PS
          mapt_s2(n)=0.0_PS
!!!          maps_s2(n)=0.0_PS
          inor(n)=0
          inos(n)=0
          T0(n)=0.0_PS
          frac(n)=1.0_PS
        enddo
        do i=1,NRBIN
          do n=1,CM%L
            CM%mt_r(n)=CM%mt_r(n)+XR(rmt_q,i,NRCAT,n)
            mapt_r2(n)=mapt_r2(n)+XR(rmat_q,i,NRCAT,n)
!!!            maps_r2(n)=maps_r2(n)+XR(rmas_q,i,NRCAT,n)
          end do
        enddo
        do i=1,NSBIN
          do n=1,CM%L
            CM%mt_s(n)=CM%mt_s(n)+XS(imt_q,i,NSCAT,n)
            mapt_s2(n)=mapt_s2(n)+XS(imat_q,i,NSCAT,n)
!!!            maps_s2(n)=maps_s2(n)+XS(imas_q,i,NSCAT,n)
          enddo
        enddo

        do i=1,NACAT
          do n=1,CM%L
            T0(n)=T0(n)+CM%mapt_ac(i,n)
          end do
        end do
        do n=1,CM%L
          T0(n)=T0(n)+CM%mapt_r(n)+CM%mapt_s(n)
        end do

        ! +++ for rain +++
        do i=1,NRBIN
          do n=1,CM%L
            if(XR(rmt_q,i,NRCAT,n)>0.0_PS.and.XR(rmat_q,i,NRCAT,n)<=0.0_PS) then
              inor(n)=inor(n)+1
            end if
            !write(fid_alog,*) "xr1,xr3",n, i, XR(rmt_q,i,NRCAT,n),XR(rmat_q,i,NRCAT,n)
          end do
        end do
        !do n=1,CM%L
        !  write(fid_alog,'(a,I3,4ES15.6)') "xa1,xa3",n, XA(amt_q,1,1,n),XA(acon_q,1,1,n),XA(amt_q,1,2,n),XA(acon_q,1,2,n)
        !enddo
        do i=1,NSBIN
          do n=1,CM%L
            if(XS(imt_q,i,NSCAT,n)>0.0_PS.and.XS(imat_q,i,NSCAT,n)<=0.0_PS) then
              inos(n)=inos(n)+1
!!!              write(fid_alog,*) "xs1,xs13",XS(imt_q,i,NSCAT,n),XS(imat_q,i,NSCAT,n)
            end if
          end do
        end do
        if(any(inos>0).or.any(inor>0)) then
          do n=1,CM%L
            if(inos(n)>0.or.inor(n)>0) then
              ! non positive total aerosol masses should be fixed in
              ! subroutine "cal_model_tendency_all"
               LOG_ERROR("check_water_apmass",*) "check_total_apmass > something is wrong, inor or inos is not 0"
               LOG_ERROR_CONT(*) n,inor(:),inos(:)
            end if
          enddo
          call PRC_abort
        endif

        ierror=0
        do n=1,CM%L
          if(mapt_a2(n)+mapt_r2(n)+mapt_s2(n)<mlmt) then
            ierror(n)=1
          else
            T1(n)=mapt_a2(n)+mapt_r2(n)+mapt_s2(n)
            frac(n)=T0(n)/T1(n)
          end if
        enddo
        if(any(ierror>0)) then
          do n=1,CM%L
            if(ierror(n)>0) then
               LOG_ERROR("check_water_apmass",*) "mass of AP is too low",n
            endif
          enddo
          call PRC_abort
        endif
        if(any(frac>5.0)) then
          do n=1,CM%L
            if(debug .and. frac(n)>5.0) then
              write(fid_alog,'("water_apmass>frac is large",4I5,ES15.6)') n,KD(n),ID(n),JD(n),frac(n)
              !write(fid_alog,'("0:T0,mta(cat),mtar,mtai",10ES15.6)') T0(n)*CM%air%TV(n)%den,CM%mapt_ac(1:NACAT,n)*CM%air%TV(n)%den&
              !     ,CM%mapt_r(n)*CM%air%TV(n)%den,CM%mapt_s(n)*CM%air%TV(n)%den
!!!              write(fid_alog,'("0:T0,mta(cat),mtar,mtai",10ES15.6)') T0(n),CM%mapt_ac(1:NACAT,n)&
!!!                   ,CM%mapt_r(n),CM%mapt_s(n)
!org              write(fid_alog,'("1:T1,mta,mtar,mtai",10ES15.6)') T1(n),mapt_a2(n),mapt_r2(n),mapt_s2(n)
              !write(fid_alog,'("1:T1,mta(cat),mtar,mtai",10ES15.6)') T1(n)*CM%air%TV(n)%den &
              !            ,CM%aerosol(1:nacat)%MS(1,n)%mass(amt) &
              !            ,mapt_r2(n)*CM%air%TV(n)%den,mapt_s2(n)*CM%air%TV(n)%den
              !write(fid_alog,'("sv1,sv2,W",3ES15.6)') CM%air%TV(n)%s_v(1),CM%air%TV(n)%s_v(2),CM%air%TV(n)%W
              !write(fid_alog,'("dsv1,dsv2",3ES15.6)') CM%air%TV(n)%s_v_n(1),CM%air%TV(n)%s_v_n(2)
              !write(fid_alog,*) "dmdt1 r", (CM%rain%MS(i,n)%dmassdt(1,1),i=1,CM%rain%N_BIN)
              !write(fid_alog,*) "dmdt1 s", (CM%solid_hydro%MS(i,n)%dmassdt(1,1),i=1,CM%solid_hydro%N_BIN)
              !write(fid_alog,*) "sum dmdt at r 1", sum(CM%rain%MS(:,n)%dmassdt(rmat,1))
              !write(fid_alog,*) "sum dmdt at s 1", sum(CM%solid_hydro%MS(:,n)%dmassdt(imat,1))
              !write(fid_alog,*) "sum dmdt at a 1", sum(CM%aerosol(1)%MS(:,n)%dmassdt(amt,1)),&
              !         sum(CM%aerosol(2)%MS(:,n)%dmassdt(amt,1)),sum(CM%aerosol(3)%MS(:,n)%dmassdt(amt,1))
              !write(fid_alog,*) "dmdt r rmt", sum(CM%rain%MS(:,n)%dmassdt(rmt,1))
              !write(fid_alog,*) "dmdt s imt", sum(CM%solid_hydro%MS(:,n)%dmassdt(imt,1))
              !Component 'CM' to right of part references with nonzero rank must not have ALLOCATABLE attribute.
            endif
          enddo
          !stop
        endif

        do i=1,nacat
          do n=1,CM%L
            if(abs(abs(frac(n))-1.0_PS)>=ok_frac.and.frac(n)<=1.0_PS) then
              XA(amt_q,1,i,n)=frac(n)*XA(amt_q,1,i,n)
              XA(acon_q,1,i,n)=XA(acon_q,1,i,n)*frac(n)
              XA(ams_q,1,i,n)=XA(ams_q,1,i,n)*frac(n)
            endif
          enddo
        enddo

!            mapt_a2=0.0_PS
!            do i=1,nacat
!              mapt_a2=mapt_a2+XA(amt_q,1,i,n)
!            end do
!
!            mapt_r2(n)=0.0_PS
!            mapt_s2(n)=0.0_PS

!CDIR NODEP
        do in=1,NRBIN*CM%L
          n=(in-1)/NRBIN+1
          i=in-(n-1)*NRBIN

          if(abs(abs(frac(n))-1.0_PS)>=ok_frac.and.frac(n)<=1.0_PS) then
            tmass1=min(frac(n)*XR(rmat_q,i,NRCAT,n),mx_aprat*XR(rmt_q,i,NRCAT,n))
            XR(rmat_q,i,NRCAT,n)=tmass1
            XR(rmas_q,i,NRCAT,n)=min(frac(n)*XR(rmas_q,i,NRCAT,n),tmass1)
          endif
        end do

!CDIR NODEP
        do in=1,NSBIN*CM%L
          n=(in-1)/NSBIN+1
          i=in-(n-1)*NSBIN
          if(abs(abs(frac(n))-1.0_PS)>=ok_frac.and.frac(n)<=1.0_PS) then
            tmass1=min(frac(n)*XS(imat_q,i,NSCAT,n),mx_aprat*XS(imt_q,i,NSCAT,n))
            XS(imat_q,i,NSCAT,n)=tmass1
            XS(imas_q,i,NSCAT,n)=min(frac(n)*XS(imas_q,i,NSCAT,n),tmass1)
          endif
        end do
      endif

      do n=1,CM%L
        CM%mt_r(n)=0.0_PS
        CM%mt_s(n)=0.0_PS
        ierror(n)=0
      enddo
      ! +++ for rain +++
      if( CM%flagp_r > 0 ) then
        do i=1,NRBIN
          do n=1,CM%L
            CM%mt_r(n)=CM%mt_r(n)+max(0.0_PS,XR(rmt_q,i,NRCAT,n)-XR(rmat_q,i,NRCAT,n))
          end do
        end do
      endif
      if( CM%flagp_s > 0 ) then
        do i=1,NSBIN
          do n=1,CM%L
            CM%mt_s(n)=CM%mt_s(n)+max(0.0_PS,XS(imt_q,i,NSCAT,n)-XS(imat_q,i,NSCAT,n))
          end do
        end do
      endif

      ! **************************************************************************************
      !                ************* NOTE **************
      ! frac can be far from 1 in cases where mass of liquid and solid hydrometeors
      ! are much smaller than vapor mixing ratio (less than one million-th).
      ! So that is a numerical limitation, instead of errors in microphysical precitions.
      ! **************************************************************************************

      ! check total water
      do n=1,CM%L
        if(CM%mt_r(n)+CM%mt_s(n)>trunc*RV(n)) then
          if(qtp(n)-RV(n)<0.0_PS) then
            RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
          end if
          frac(n)=(qtp(n)-RV(n))/(CM%mt_r(n)+CM%mt_s(n))
        else
          RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
          frac(n)=1.0
        endif
      enddo

      if(debug) then
        do n=1,CM%L
          write(fid_alog,'("frac",I5,10ES15.6)') n,frac(n),qtp(n),rv(n),CM%mt_r(n),CM%mt_s(n) &
                 ,CM%mt_r(n)*CM%air%TV(n)%den,CM%mt_s(n)*CM%air%TV(n)%den
        enddo
      endif

      icond1=0
      do n=1,CM%L
        if(CM%mt_r(n)+CM%mt_s(n)>trunc*RV(n)) then
          if(abs(abs(frac(n))-1.0_PS)>=ok_frac) then
            if(frac(n)>1.0_PS.or.frac(n)==0.0_PS) then
              ! this can create high supersaturation over water if mt_r and mt_s
              !   are not well predicted.
              RV(n)=max(0.0_PS,qtp(n)-CM%mt_r(n)-CM%mt_s(n))
            else
              icond1(n)=1
            endif
          endif
        endif
      enddo

      ! +++ for liquid +++
!CDIR NODEP
      do ijn=1,NRTYPE*NRBIN*CM%L
        in=(ijn-1)/NRTYPE+1
        j=ijn-(in-1)*NRTYPE
        n=(in-1)/NRBIN+1
        i=in-(n-1)*NRBIN
        if(icond1(n)==1) then
          if(XR(j,i,NRCAT,n)>0.0_PS) then
            XR(j,i,NRCAT,n)=XR(j,i,NRCAT,n)*frac(n)
          end if
        end if
      end do
      ! +++ for ice +++
!CDIR NODEP
      do ijn=1,NSTYPE*NSBIN*CM%L
        in=(ijn-1)/NSTYPE+1
        j=ijn-(in-1)*NSTYPE
        n=(in-1)/NSBIN+1
        i=in-(n-1)*NSBIN
        if(icond1(n)==1) then
          if(XS(j,i,NSCAT,n)>0.0_PS) then
            XS(j,i,NSCAT,n)=XS(j,i,NSCAT,n)*frac(n)
          end if
        end if
      end do

    end if

  end subroutine check_water_apmass

  subroutine cal_dmtend(dmtend,tcon,tmass,tvar,dctend_in,iproc,&
       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth,&
       mpirank_w,mpirank_e,mpirank_s,mpirank_n)
! <<< 2014/10 T. Hashino added for KiD
  use par_amps, only : isplit_bin_liq
! >>> 2014/10 T. Hashino added for KiD
    real(PS),intent(inout) :: dmtend(11,9,*),tcon(9,*),tmass(9,*)
    real(PS),intent(inout) :: tvar(4,*),dctend_in(7,*)
    integer,intent(in) :: iproc
    integer,intent(in) :: ID(*),JD(*),KD(*),nzp,nxp,nyp,nfpt,nfpth
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in)  :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),&
            zstv(*),zz(*),xx(*),yy(*)
    integer,intent(in) :: mpirank_w,mpirank_e,mpirank_s,mpirank_n

    integer :: jproc
    integer :: n,i,j,ish,i1,i2,j1,j2
    real(PS) :: dvol,dcon(9),tconl(9),tmassl(9)
    real(PS),parameter :: cloudsz=50.0e-4

!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    if(nxp==3) then
      i1=1
      i2=3
    else
      if(mpirank_w==-9) then
        i1=3
      else
        i1=1
      endif
      if(mpirank_e==-9) then
        i2=nxp-3
      else
        i2=nxp
      endif
    endif
    if(mpirank_s==-9) then
      j1=3
    else
      j1=1
    endif
    if(mpirank_n==-9) then
      j2=nyp-3
    else
      j2=nyp
    endif

    jproc=iproc+1

    do n=1,CM%L
       ! +++ loop over grids +++
       if( (ID(n)<=i1.or.ID(n)>=i2).or.(JD(n)<=j1.or.JD(n)>=j2).or.&
           (KD(n)>=nzp-nfpt+1)) cycle

       !  1.0e+6 make the unit to #/m^3, g/m^3, g/m^3/s
       dvol=(xx(ID(n))-xx(ID(n)-1))*(yy(JD(n))-yy(JD(n)-1))&
               *(zz(KD(n))-max(zstv(n),zz(KD(n)-1)))*1.0e+6


       ! +++ for rain +++

          dcon=0.0_PS
          tconl=0.0_PS
          tmassl=0.0_PS

!tmp          dvol=1.0

          do i=1,CM%rain%N_BIN
             if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XR(1,i,NRCAT,n)<=1.0e-25.or.XR(2,i,NRCAT,n)<=1.0e-25)) cycle
!tmp             if(CM%rain%MS(i,n)%len<=cloudsz) then
!tmp                ish=1
!tmp             else
!tmp                ish=2
!tmp             endif
             if(i<=isplit_bin_liq) then
                ish=1
             else
                ish=2
             endif
             ! vapor deposition
             dmtend(1,ish,jproc)=dmtend(1,ish,jproc)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,jproc)=dmtend(2,ish,jproc)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dmtend(3,ish,jproc)=dmtend(3,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! CCN activation
             dmtend(4,ish,jproc)=dmtend(4,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)&
                   *dvol
             !
             ! riming
             dmtend(5,ish,jproc)=dmtend(5,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,jproc)=dmtend(6,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,jproc)=dmtend(7,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,8)&
                   *dvol
             !
             ! ice nucleation (immersion)
             dmtend(8,ish,jproc)=dmtend(8,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,11)&
                   *dvol
             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,jproc)=dmtend(9,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,12)&
                   *dvol
             !
             ! ice nucleation (secondary)
             dmtend(10,ish,jproc)=dmtend(10,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,9)&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dcon(ish)=dcon(ish)&
                   +CM%rain%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%rain%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%rain%MS(i,n)%mass(1)*dvol

             ! mean capacitance
             tvar(1,jproc)=tvar(1,jproc)+CM%rain%MS(i,n)%con*CM%rain%MS(i,n)%CAP*dvol
          enddo

       ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25.or.CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XS(1,i,NSCAT,n)<=1.0e-25.or.XS(2,i,NSCAT,n)<=1.0e-25)) cycle

             ish=CM%solid_hydro%IS(i,n)%sh_type
             ish=ish+2
!             if(ish<=2) then
!                ish=3
!             elseif(ish<=4) then
!                ish=4
!             else
!                ish=5
!             endif

             ! vapor deposition
             dmtend(1,ish,jproc)=dmtend(1,ish,jproc)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,jproc)=dmtend(2,ish,jproc)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! aggregation
             dmtend(3,ish,jproc)=dmtend(3,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! ice nucleation (deposition)
             dmtend(4,ish,jproc)=dmtend(4,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)&
                   *dvol

             dctend_in(2,jproc)=dctend_in(2,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(7)&
                   *dvol

             !
             ! riming
             dmtend(5,ish,jproc)=dmtend(5,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,jproc)=dmtend(6,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,jproc)=dmtend(7,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)&
                   *dvol
             dctend_in(3,jproc)=dctend_in(3,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(8)&
                   *dvol

             !
             ! ice nucleation (immersion)
             dmtend(8,ish,jproc)=dmtend(8,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)&
                   *dvol
             dctend_in(4,jproc)=dctend_in(4,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(11)&
                   *dvol

             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,jproc)=dmtend(9,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)&
                   *dvol
             dctend_in(5,jproc)=dctend_in(5,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(12)&
                   *dvol

             !
             ! ice nucleation (secondary)
             dmtend(10,ish,jproc)=dmtend(10,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)&
                   *dvol
             dctend_in(6,jproc)=dctend_in(6,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(9)&
                   *dvol


             ! aggregation
             dcon(ish)=dcon(ish)&
                   +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol

             ! DHF freezing
             dctend_in(7,jproc)=dctend_in(7,jproc)&
                   +CM%solid_hydro%MS(i,n)%dcondt(3)&
                   *dvol


             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%solid_hydro%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%solid_hydro%MS(i,n)%mass(1)*dvol

             ! mean capacitance
             tvar(2,jproc)=tvar(2,jproc)+&
                       CM%solid_hydro%MS(i,n)%con*CM%solid_hydro%MS(i,n)%CAP*dvol
             tvar(3,jproc)=tvar(3,jproc)+&
                       CM%solid_hydro%MS(i,n)%con*CM%solid_hydro%MS(i,n)%den*dvol
             tvar(4,jproc)=tvar(4,jproc)+&
                       CM%solid_hydro%MS(i,n)%con* &
                       CM%solid_hydro%IS(i,n)%semi_cip/max(1.0e-4_RP,CM%solid_hydro%IS(i,n)%semi_aip)*&
!org                       CM%solid_hydro%MS(i,n)%c_len/max(1.0e-4_RP,CM%solid_hydro%MS(i,n)%a_len)*&
                       dvol

          enddo

          do ish=1,2
             tconl(8)=tconl(8)+tconl(ish)
             tmassl(8)=tmassl(8)+tmassl(ish)
             dcon(8)=dcon(8)+dcon(ish)
          enddo
          do ish=3,7
             tconl(9)=tconl(9)+tconl(ish)
             tmassl(9)=tmassl(9)+tmassl(ish)
             dcon(9)=dcon(9)+dcon(ish)
          enddo
! commentted out for SHEBA
!          do ish=1,9
!             if(tconl(ish)>1.0e-20_PS) then
!               dmtend(11,ish,jproc)=dmtend(11,ish,jproc)+&
!                   (-dcon(ish)*tmassl(ish)/&
!                   (tconl(ish)+dcon(ish)*CM%rain%dt))
!             endif
!          enddo
! end commentted out for SHEBA
          do ish=1,7
             tcon(ish,jproc)=tcon(ish,jproc)+tconl(ish)
             tmass(ish,jproc)=tmass(ish,jproc)+tmassl(ish)
          enddo

          dctend_in(1,jproc)=sum(dctend_in(2:7,jproc))


    enddo

  end subroutine cal_dmtend

  subroutine cal_dmtend_pm(dmtend,tcon,tmass,tvar,dctend_in,iproc,&
       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XR,XS,KD,ID,JD)
! <<< 2014/10 T. Hashino added for KiD
  use par_amps, only : isplit_bin_liq
! >>> 2014/10 T. Hashino added for KiD
    real(PS) :: dmtend(11,9,*),tcon(9,*),tmass(9,*),tvar(*),dctend_in(*)
    integer,intent(in) :: iproc
    integer,intent(in) :: ID(*),JD(*),KD(*)
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in)  :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*)

    integer :: jproc
    integer :: n,i,j,ish
    real(PS) :: dvol,dcon(9),tconl(9),tmassl(9)
    real(PS),parameter :: cloudsz=50.0e-4

!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    jproc=iproc+1

    do n=1,CM%L

       dvol=1.0


       ! +++ for rain +++

          dcon=0.0_PS
          tconl=0.0_PS
          tmassl=0.0_PS


          do i=1,CM%rain%N_BIN
             if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XR(1,i,NRCAT,n)<=1.0e-25.or.XR(2,i,NRCAT,n)<=1.0e-25)) cycle
             if(CM%rain%MS(i,n)%len<=cloudsz) then
               ish=1
             else
               ish=2
             endif
             if(i<=isplit_bin_liq) then
                ish=1
             else
                ish=2
             endif
             ! vapor deposition
             dmtend(1,ish,jproc)=dmtend(1,ish,jproc)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,jproc)=dmtend(2,ish,jproc)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dmtend(3,ish,jproc)=dmtend(3,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! CCN activation
             dmtend(4,ish,jproc)=dmtend(4,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)&
                   *dvol
             !
             ! riming
             dmtend(5,ish,jproc)=dmtend(5,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,jproc)=dmtend(6,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,jproc)=dmtend(7,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,8)&
                   *dvol
             !
             ! ice nucleation (immersion)
             dmtend(8,ish,jproc)=dmtend(8,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,11)&
                   *dvol
             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,jproc)=dmtend(9,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,12)&
                   *dvol
             !
             ! ice nucleation (secondary)
             dmtend(10,ish,jproc)=dmtend(10,ish,jproc)&
                   +CM%rain%MS(i,n)%dmassdt(1,9)&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dcon(ish)=dcon(ish)&
                   +CM%rain%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%rain%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%rain%MS(i,n)%mass(1)*dvol

             ! mean capacitance
             tvar(1)=tvar(1)+CM%rain%MS(i,n)%con*CM%rain%MS(i,n)%CAP
          enddo

       ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25.or.CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XS(1,i,NSCAT,n)<=1.0e-25.or.XS(2,i,NSCAT,n)<=1.0e-25)) cycle

             ish=CM%solid_hydro%IS(i,n)%sh_type
             ish=ish+2
!             if(ish<=2) then
!                ish=3
!             elseif(ish<=4) then
!                ish=4
!             else
!                ish=5
!             endif

             ! vapor deposition
             dmtend(1,ish,jproc)=dmtend(1,ish,jproc)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,jproc)=dmtend(2,ish,jproc)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! aggregation
             dmtend(3,ish,jproc)=dmtend(3,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! ice nucleation (deposition)
             dmtend(4,ish,jproc)=dmtend(4,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)&
                   *dvol
             dctend_in(2)=dctend_in(2)&
                   +CM%solid_hydro%MS(i,n)%dcondt(7)&
                   *dvol
             !
             ! riming
             dmtend(5,ish,jproc)=dmtend(5,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,jproc)=dmtend(6,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,jproc)=dmtend(7,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)&
                   *dvol
             dctend_in(3)=dctend_in(3)&
                   +CM%solid_hydro%MS(i,n)%dcondt(8)&
                   *dvol
             !
             ! ice nucleation (immersion)
             dmtend(8,ish,jproc)=dmtend(8,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)&
                   *dvol
             dctend_in(4)=dctend_in(4)&
                   +CM%solid_hydro%MS(i,n)%dcondt(11)&
                   *dvol
             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,jproc)=dmtend(9,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)&
                   *dvol
             dctend_in(5)=dctend_in(5)&
                   +CM%solid_hydro%MS(i,n)%dcondt(12)&
                   *dvol
             !
             ! ice nucleation (secondary)
             dmtend(10,ish,jproc)=dmtend(10,ish,jproc)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)&
                   *dvol
             dctend_in(6)=dctend_in(6)&
                   +CM%solid_hydro%MS(i,n)%dcondt(9)&
                   *dvol

             ! aggregation
             dcon(ish)=dcon(ish)&
                   +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol

             ! DHF freezing
             dctend_in(7)=dctend_in(7)&
                   +CM%solid_hydro%MS(i,n)%dcondt(3)&
                   *dvol

             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%solid_hydro%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%solid_hydro%MS(i,n)%mass(1)*dvol

             ! mean capacitance
             tvar(2)=tvar(2)+CM%solid_hydro%MS(i,n)%con*CM%solid_hydro%MS(i,n)%CAP
             tvar(3)=tvar(3)+CM%solid_hydro%MS(i,n)%con*CM%solid_hydro%MS(i,n)%den
             tvar(4)=tvar(4)+CM%solid_hydro%MS(i,n)%con* &
                 CM%solid_hydro%MS(i,n)%c_len/max(1.0e-4_RP,CM%solid_hydro%MS(i,n)%a_len)
          enddo

          do ish=1,2
             tconl(8)=tconl(8)+tconl(ish)
             tmassl(8)=tmassl(8)+tmassl(ish)
             dcon(8)=dcon(8)+dcon(ish)
          enddo
          do ish=3,7
             tconl(9)=tconl(9)+tconl(ish)
             tmassl(9)=tmassl(9)+tmassl(ish)
             dcon(9)=dcon(9)+dcon(ish)
          enddo
! commentted out for SHEBA
!          do ish=1,9
!             if(tconl(ish)>1.0e-20_PS) then
!               dmtend(11,ish,jproc)=dmtend(11,ish,jproc)+&
!                   (-dcon(ish)*tmassl(ish)/&
!                   (tconl(ish)+dcon(ish)*CM%rain%dt))
!             endif
!          enddo
! end commentted out for SHEBA
          do ish=1,9
             tcon(ish,jproc)=tcon(ish,jproc)+tconl(ish)
             tmass(ish,jproc)=tmass(ish,jproc)+tmassl(ish)
          enddo

          if(sum(tcon(1:2,jproc))>1.0e-20) then
            tvar(1)=tvar(1)/(tcon(1,jproc)+tcon(2,jproc))
          endif
          if(sum(tcon(3:7,jproc))>1.0e-20) then
            tvar(2)=tvar(2)/(sum(tcon(3:7,jproc)))
            tvar(3)=tvar(3)/(sum(tcon(3:7,jproc)))
            tvar(4)=tvar(4)/(sum(tcon(3:7,jproc)))
          endif

          dctend_in(1)=sum(dctend_in(2:7))

    enddo

  end subroutine cal_dmtend_pm


  subroutine cal_dmtend_scale(dmtendl,dcontendl,dbintendl,L,CM)

    use par_amps, only : isplit_bin_liq
    implicit none
    type (Cloud_Micro), intent(in) :: CM

    !integer,intent(in)  :: NSCAT
    !real(PS) :: XS(NSTYPE,NSBIN,NSCAT,*)
    integer, intent(in)    :: L

    real(MP_KIND),    intent(inout) :: dmtendl(10,2,L)
    real(MP_KIND),    intent(inout) :: dcontendl(10,2,L)
    real(MP_KIND),    intent(inout) :: dbintendl(3,2,mxnbin,L)

    integer :: n,i,j,ish
    real(PS) :: dvol,dcon(9),tconl(9),tmassl(9)
    real(PS),parameter :: cloudsz=50.0e-4

!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    ! flowchart of mass tendency calculation
    ! get Np=newly formed number concentration and Mp=newly formed total mass
    ! call assign_Qp_v3p to obtain non-mass PPVs as Qp
    ! call cal_xxx_p_v5 to check whether Qp make sense, and get espect ratio, habits, and etc.
    ! call cal_linprms_vec_s calculate the coefficients of cubic or linear distribution in each bin
    ! call cal_transbin_vec to transfer bin mass
    ! dmassdt, dcondt, dvoldt ( non-mass ice PPV calculated from dvoldt)
    ! ---
    ! dmassdt(2)  = aggregation (liq-liq or ice-ice)         = dmtendl(3) [ coalescence ]
    ! dmassdt(4)  = riming (liq-ice)                         = dmtendl(5) [ coalescence ]
    ! dmassdt(10) = melting-shedding                         = dmtendl(6) [ coalescence ]
    ! dmassdt(6)  = hydrodynamics breakup                    = dmtendl(none) [ hydrodyn_breakup ]
    ! dmassdt(8)  = contact freezing                         = dmtendl(7) [ ice_nucleation1 ]
    ! dmassdt(9)  = secondary (Hallett-Mossop) freezing      = dmtendl(10) [ ice_nucleation1 ]
    ! dmassdt(10) = melting                                  = dmtendl(11)[ coalescence ]
    ! dmassdt(11) = immersion freezing (Meyers or KC theory) = dmtendl(8) [ ice_nucleation1 ]
    ! dmassdt(12) = homogeneous freezing (Heymsfield)        = dmtendl(9) [ ice_nucleation1 ]
    ! dmassdt(3)  = CCN activation (KC theory)               = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(3)  = Deliquescence freezing (KC theory)       = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(7)  = Deposition freezing (KC theory)          = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(1)  = drop condensation                        = dmtendl(1) [ vapor_deposition ]
    ! dmassdt(1)  = ice deposition                           = dmtendl(1) [ vapor_deposition ]
    ! dmassdt(1)  = drop evaporation                         = dmtendl(2) [ vapor_deposition ]
    ! dmassdt(1)  = ice evaporation                          = dmtendl(2) [ vapor_deposition ]
    ! ---

    dvol = 0.0
    dcon = 0.0
    tconl = 0.0
    tmassl = 0.0
    do n=1,L
       ! +++ loop over grids +++
       ! +++ for rain +++


          dvol=1.0

          do i=1,CM%rain%N_BIN
             !if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
             !   (XR(rmt_q,i,NRCAT,n)<=1.0e-25.or.XR(rcon_q,i,NRCAT,n)<=1.0e-25)) cycle
!tmp             if(CM%rain%MS(i,n)%len<=cloudsz) then
!tmp                ish=1
!tmp             else
!tmp                ish=2
!tmp             endif
             if(i<=isplit_bin_liq) then
!             if(i<=19) then
                ish=1
             else
                ish=2
             endif
             ! vapor deposition
             dbintendl(2,1,i,n)=dbintendl(2,1,i,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dbintendl(3,1,i,n)=dbintendl(3,1,i,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! riming
             dbintendl(1,1,i,n)=dbintendl(1,1,i,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol

             ! vapor deposition
             dmtendl(1,1,n) = dmtendl(1,1,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(1,1,n) = dcontendl(1,1,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dcondt(1))&
                   *dvol

             ! evaporation
             dmtendl(2,1,n) = dmtendl(2,1,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(2,1,n) = dcontendl(2,1,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dcondt(1))&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dmtendl(3,1,n) = dmtendl(3,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)&
                   *dvol

             dcontendl(3,1,n) = dcontendl(3,1,n)&
                   +CM%rain%MS(i,n)%dcondt(2)&
                   *dvol

             ! CCN activation
             dmtendl(4,1,n) = dmtendl(4,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)&
                   *dvol

             dcontendl(4,1,n) = dcontendl(4,1,n)&
                   +CM%rain%MS(i,n)%dcondt(3)&
                   *dvol

             ! riming
             dmtendl(5,1,n) = dmtendl(5,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol

             dcontendl(5,1,n) = dcontendl(5,1,n)&
                   +CM%rain%MS(i,n)%dcondt(4)&
                   *dvol

             ! melting-shedding
             dmtendl(6,1,n) = dmtendl(6,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)&
                   *dvol

             dcontendl(6,1,n) = dcontendl(6,1,n)&
                   +CM%rain%MS(i,n)%dcondt(10)&
                   *dvol

             ! ice nucleation (contact)
             dmtendl(7,1,n) = dmtendl(7,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,8)&
                   *dvol

             dcontendl(7,1,n) = dcontendl(7,1,n)&
                   +CM%rain%MS(i,n)%dcondt(8)&
                   *dvol

             ! ice nucleation (immersion)
             dmtendl(8,1,n) = dmtendl(8,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,11)&
                   *dvol

             dcontendl(8,1,n) = dcontendl(8,1,n)&
                   +CM%rain%MS(i,n)%dcondt(11)&
                   *dvol

             ! ice nucleation (homo freezing)
             dmtendl(9,1,n) = dmtendl(9,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,12)&
                   *dvol

             dcontendl(9,1,n) = dcontendl(9,1,n)&
                   +CM%rain%MS(i,n)%dcondt(12)&
                   *dvol

             ! ice nucleation (secondary)
             dmtendl(10,1,n) = dmtendl(10,1,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,9)&
                   *dvol

             dcontendl(10,1,n) = dcontendl(10,1,n)&
                   +CM%rain%MS(i,n)%dcondt(9)&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             !dcon(ish)=dcon(ish)&
             !      +CM%rain%MS(i,n)%dcondt(2)*dvol
             !dmtendl(11,n)=dmtendl(11,n)&
             !      +CM%rain%MS(i,n)%dmassdt(1,2)*dvol

             !dcontendl(11,n)=dcontendl(11,n)&
             !      +CM%rain%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             !tconl(ish)=tconl(ish)+CM%rain%MS(i,n)%con*dvol
             !tmassl(ish)=tmassl(ish)+CM%rain%MS(i,n)%mass(1)*dvol
          enddo

          ! total liquid
          !do ish=1,2
          !   tconl(8)=tconl(8)+tconl(ish)
          !   tmassl(8)=tmassl(8)+tmassl(ish)
          !   dcon(8)=dcon(8)+dcon(ish)
          !enddo
          !do ish=1,2
          !   tcon(ish,n)=tcon(ish,n)+tconl(ish)
          !   tmass(ish,n)=tmass(ish,n)+tmassl(ish)
          !enddo

          ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25 .or. &
                 CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25) & ! .and. &
                !(XS(imt_q,i,NSCAT,n)<=1.0e-25.or.XS(icon_q,i,NSCAT,n)<=1.0e-25) .and. &
               ) cycle

             ! vapor deposition
             dbintendl(2,2,i,n)=dbintendl(2,2,i,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dbintendl(3,2,i,n)=dbintendl(3,2,i,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! riming
             dbintendl(1,2,i,n)=dbintendl(1,2,i,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol

             ish=CM%solid_hydro%IS(i,n)%sh_type
             ish=ish+2
             ! type of solid hydrometeor
             ! 1+2 : pristine crystal
             ! 2+2 : rimed crystal
             ! 3+2 : aggregate
             ! 4+2 : rimed aggregate
             ! 5+2 : graupel
             ! 6+2 : hail

             ! vapor deposition
             dmtendl(1,2,n)=dmtendl(1,2,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(1,2,n)=dcontendl(1,2,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dcondt(1))&
                   *dvol

             ! evaporation
             dmtendl(2,2,n)=dmtendl(2,2,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(2,2,n)=dcontendl(2,2,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dcondt(1))&
                   *dvol

             ! aggregation
             dmtendl(3,2,n)=dmtendl(3,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,2)&
                   *dvol

             dcontendl(3,2,n)=dcontendl(3,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(2)&
                   *dvol

             ! ice nucleation (deposition)
             dmtendl(4,2,n)=dmtendl(4,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)&
                   *dvol&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,3)&
                   *dvol

             dcontendl(4,2,n)=dcontendl(4,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(7)&
                   *dvol&
                   +CM%solid_hydro%MS(i,n)%dcondt(3)&
                   *dvol

             ! riming
             dmtendl(5,2,n)=dmtendl(5,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol

             dcontendl(5,2,n)=dcontendl(5,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(4)&
                   *dvol

             ! melting-shedding
             dmtendl(6,2,n)=dmtendl(6,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,10)&
                   *dvol

             dcontendl(6,2,n)=dcontendl(6,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(10)&
                   *dvol

             ! ice nucleation (contact)
             dmtendl(7,2,n)=dmtendl(7,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)&
                   *dvol

             dcontendl(7,2,n)=dcontendl(7,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(8)&
                   *dvol

             ! ice nucleation (immersion)
             dmtendl(8,2,n)=dmtendl(8,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)&
                   *dvol

             dcontendl(8,2,n)=dcontendl(8,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(11)&
                   *dvol

             ! ice nucleation (homo freezing)
             dmtendl(9,2,n)=dmtendl(9,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)&
                   *dvol

             dcontendl(9,2,n)=dcontendl(9,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(12)&
                   *dvol

             ! ice nucleation (secondary)
             dmtendl(10,2,n)=dmtendl(10,2,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)&
                   *dvol

             dcontendl(10,2,n)=dcontendl(10,2,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(9)&
                   *dvol

             ! aggregation -- this is temporarily replaced by melting process
             !dcon(ish)=dcon(ish)&
             !      +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol
             !dmtendl(11,ish,n)=dmtendl(11,ish,n)&
             !      +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol
             !dmtendl(11,ish,n)=dmtendl(11,ish,n)+(&
             !       CM%solid_hydro%MS(i,n)%dmassdt(2,10)&
             !      +CM%solid_hydro%MS(i,n)%dmassdt(3,10)&
             !      +CM%solid_hydro%MS(i,n)%dmassdt(4,10))&
             !      *dvol

             ! total concentration and mass
             !tconl(ish)=tconl(ish)+CM%solid_hydro%MS(i,n)%con*dvol
             !tmassl(ish)=tmassl(ish)+CM%solid_hydro%MS(i,n)%mass(1)*dvol
          enddo

          !do ish=3,7
          !   tconl(9)=tconl(9)+tconl(ish)
          !   tmassl(9)=tmassl(9)+tmassl(ish)
          !   dcon(9)=dcon(9)+dcon(ish)
          !enddo
          !do ish=1,9
          !   if(tconl(ish)>1.0e-20_PS) then
          !     dmtend(11,ish,n)=dmtend(11,ish,n)+&
          !         (-dcon(ish)*tmassl(ish)/&
          !         (tconl(ish)+dcon(ish)*CM%rain%dt))
          !   endif
          !enddo
          !do ish=1,7
          !   tcon(ish,n)=tcon(ish,n)+tconl(ish)
          !   tmass(ish,n)=tmass(ish,n)+tmassl(ish)
          !enddo

    enddo

  end subroutine cal_dmtend_scale


  subroutine cal_dmtend_scale_original(dmtendl,dcontendl,L,CM)

    use par_amps, only : isplit_bin_liq

    type (Cloud_Micro), intent(in) :: CM

    !integer,intent(in)  :: NSCAT
    !real(PS) :: XS(NSTYPE,NSBIN,NSCAT,*)

    integer, intent(in)    :: L
    real(PS),    intent(inout) :: dmtendl(11,8,L)
    real(PS),    intent(inout) :: dcontendl(11,8,L)

    integer :: n,i,j,ish
    real(PS) :: dvol,dcon(9),tconl(9),tmassl(9)
    real(PS),parameter :: cloudsz=50.0e-4

!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    ! flowchart of mass tendency calculation
    ! get Np=newly formed number concentration and Mp=newly formed total mass
    ! call assign_Qp_v3p to obtain non-mass PPVs as Qp
    ! call cal_xxx_p_v5 to check whether Qp make sense, and get espect ratio, habits, and etc.
    ! call cal_linprms_vec_s calculate the coefficients of cubic or linear distribution in each bin
    ! call cal_transbin_vec to transfer bin mass
    ! dmassdt, dcondt, dvoldt ( non-mass ice PPV calculated from dvoldt)
    ! ---
    ! dmassdt(2)  = aggregation (liq-liq or ice-ice)         = dmtendl(3) [ coalescence ]
    ! dmassdt(4)  = riming (liq-ice)                         = dmtendl(5) [ coalescence ]
    ! dmassdt(10) = melting-shedding                         = dmtendl(6) [ coalescence ]
    ! dmassdt(6)  = hydrodynamics breakup                    = dmtendl(none) [ hydrodyn_breakup ]
    ! dmassdt(8)  = contact freezing                         = dmtendl(7) [ ice_nucleation1 ]
    ! dmassdt(9)  = secondary (Hallett-Mossop) freezing      = dmtendl(10) [ ice_nucleation1 ]
    ! dmassdt(10) = melting                                  = dmtendl(11)[ coalescence ]
    ! dmassdt(11) = immersion freezing (Meyers or KC theory) = dmtendl(8) [ ice_nucleation1 ]
    ! dmassdt(12) = homogeneous freezing (Heymsfield)        = dmtendl(9) [ ice_nucleation1 ]
    ! dmassdt(3)  = CCN activation (KC theory)               = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(3)  = Deliquescence freezing (KC theory)       = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(7)  = Deposition freezing (KC theory)          = dmtendl(4) [ cal_aptact_var8_kc04dep ]
    ! dmassdt(1)  = drop condensation                        = dmtendl(1) [ vapor_deposition ]
    ! dmassdt(1)  = ice deposition                           = dmtendl(1) [ vapor_deposition ]
    ! dmassdt(1)  = drop evaporation                         = dmtendl(2) [ vapor_deposition ]
    ! dmassdt(1)  = ice evaporation                          = dmtendl(2) [ vapor_deposition ]
    ! ---

    dvol = 0.0
    dcon = 0.0
    tconl = 0.0
    tmassl = 0.0
    do n=1,L
       ! +++ loop over grids +++
       ! +++ for rain +++


          dvol=1.0

          do i=1,CM%rain%N_BIN
             !if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
             !   (XR(rmt_q,i,NRCAT,n)<=1.0e-25.or.XR(rcon_q,i,NRCAT,n)<=1.0e-25)) cycle
!tmp             if(CM%rain%MS(i,n)%len<=cloudsz) then
!tmp                ish=1
!tmp             else
!tmp                ish=2
!tmp             endif
             if(i<=isplit_bin_liq) then
                ish=1
             else
                ish=2
             endif
             ! vapor deposition
             dmtendl(1,ish,n) = dmtendl(1,ish,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(1,ish,n) = dcontendl(1,ish,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dcondt(1))&
                   *dvol

             ! evaporation
             dmtendl(2,ish,n) = dmtendl(2,ish,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(2,ish,n) = dcontendl(2,ish,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dcondt(1))&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dmtendl(3,ish,n) = dmtendl(3,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)&
                   *dvol

             dcontendl(3,ish,n) = dcontendl(3,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(2)&
                   *dvol

             ! CCN activation
             dmtendl(4,ish,n) = dmtendl(4,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)&
                   *dvol

             dcontendl(4,ish,n) = dcontendl(4,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(3)&
                   *dvol

             ! riming
             dmtendl(5,ish,n) = dmtendl(5,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol

             dcontendl(5,ish,n) = dcontendl(5,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(4)&
                   *dvol

             ! melting-shedding
             dmtendl(6,ish,n) = dmtendl(6,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)&
                   *dvol

             dcontendl(6,ish,n) = dcontendl(6,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(10)&
                   *dvol

             ! ice nucleation (contact)
             dmtendl(7,ish,n) = dmtendl(7,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,8)&
                   *dvol

             dcontendl(7,ish,n) = dcontendl(7,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(8)&
                   *dvol

             ! ice nucleation (immersion)
             dmtendl(8,ish,n) = dmtendl(8,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,11)&
                   *dvol

             dcontendl(8,ish,n) = dcontendl(8,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(11)&
                   *dvol

             ! ice nucleation (homo freezing)
             dmtendl(9,ish,n) = dmtendl(9,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,12)&
                   *dvol

             dcontendl(9,ish,n) = dcontendl(9,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(12)&
                   *dvol

             ! ice nucleation (secondary)
             dmtendl(10,ish,n) = dmtendl(10,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,9)&
                   *dvol

             dcontendl(10,ish,n) = dcontendl(10,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(9)&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             !dcon(ish)=dcon(ish)&
             !      +CM%rain%MS(i,n)%dcondt(2)*dvol
             dmtendl(11,ish,n)=dmtendl(11,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)*dvol

             dcontendl(11,ish,n)=dcontendl(11,ish,n)&
                   +CM%rain%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             !tconl(ish)=tconl(ish)+CM%rain%MS(i,n)%con*dvol
             !tmassl(ish)=tmassl(ish)+CM%rain%MS(i,n)%mass(1)*dvol
          enddo

          ! total liquid
          !do ish=1,2
          !   tconl(8)=tconl(8)+tconl(ish)
          !   tmassl(8)=tmassl(8)+tmassl(ish)
          !   dcon(8)=dcon(8)+dcon(ish)
          !enddo
          !do ish=1,2
          !   tcon(ish,n)=tcon(ish,n)+tconl(ish)
          !   tmass(ish,n)=tmass(ish,n)+tmassl(ish)
          !enddo

          ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25 .or. &
                 CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25) & ! .and. &
                !(XS(imt_q,i,NSCAT,n)<=1.0e-25.or.XS(icon_q,i,NSCAT,n)<=1.0e-25) .and. &
               ) cycle

             ish=CM%solid_hydro%IS(i,n)%sh_type
             ish=ish+2
             ! type of solid hydrometeor
             ! 1+2 : pristine crystal
             ! 2+2 : rimed crystal
             ! 3+2 : aggregate
             ! 4+2 : rimed aggregate
             ! 5+2 : graupel
             ! 6+2 : hail

             ! vapor deposition
             dmtendl(1,ish,n)=dmtendl(1,ish,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(1,ish,n)=dcontendl(1,ish,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dcondt(1))&
                   *dvol

             ! evaporation
             dmtendl(2,ish,n)=dmtendl(2,ish,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             dcontendl(2,ish,n)=dcontendl(2,ish,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dcondt(1))&
                   *dvol

             ! aggregation
             dmtendl(3,ish,n)=dmtendl(3,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,2)&
                   *dvol

             dcontendl(3,ish,n)=dcontendl(3,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(2)&
                   *dvol

             ! ice nucleation (deposition)
             dmtendl(4,ish,n)=dmtendl(4,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)&
                   *dvol&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,3)&
                   *dvol

             dcontendl(4,ish,n)=dcontendl(4,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(7)&
                   *dvol&
                   +CM%solid_hydro%MS(i,n)%dcondt(3)&
                   *dvol

             ! riming
             dmtendl(5,ish,n)=dmtendl(5,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol

             dcontendl(5,ish,n)=dcontendl(5,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(4)&
                   *dvol

             ! melting-shedding
             dmtendl(6,ish,n)=dmtendl(6,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,10)&
                   *dvol

             dcontendl(6,ish,n)=dcontendl(6,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(10)&
                   *dvol

             ! ice nucleation (contact)
             dmtendl(7,ish,n)=dmtendl(7,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)&
                   *dvol

             dcontendl(7,ish,n)=dcontendl(7,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(8)&
                   *dvol

             ! ice nucleation (immersion)
             dmtendl(8,ish,n)=dmtendl(8,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)&
                   *dvol

             dcontendl(8,ish,n)=dcontendl(8,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(11)&
                   *dvol

             ! ice nucleation (homo freezing)
             dmtendl(9,ish,n)=dmtendl(9,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)&
                   *dvol

             dcontendl(9,ish,n)=dcontendl(9,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(12)&
                   *dvol

             ! ice nucleation (secondary)
             dmtendl(10,ish,n)=dmtendl(10,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)&
                   *dvol

             dcontendl(10,ish,n)=dcontendl(10,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dcondt(9)&
                   *dvol

             ! aggregation -- this is temporarily replaced by melting process
             !dcon(ish)=dcon(ish)&
             !      +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol
             !dmtendl(11,ish,n)=dmtendl(11,ish,n)&
             !      +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol
             dmtendl(11,ish,n)=dmtendl(11,ish,n)+(&
                    CM%solid_hydro%MS(i,n)%dmassdt(2,10)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(3,10)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(4,10))&
                   *dvol

             ! total concentration and mass
             !tconl(ish)=tconl(ish)+CM%solid_hydro%MS(i,n)%con*dvol
             !tmassl(ish)=tmassl(ish)+CM%solid_hydro%MS(i,n)%mass(1)*dvol
          enddo

          !do ish=3,7
          !   tconl(9)=tconl(9)+tconl(ish)
          !   tmassl(9)=tmassl(9)+tmassl(ish)
          !   dcon(9)=dcon(9)+dcon(ish)
          !enddo
          !do ish=1,9
          !   if(tconl(ish)>1.0e-20_PS) then
          !     dmtend(11,ish,n)=dmtend(11,ish,n)+&
          !         (-dcon(ish)*tmassl(ish)/&
          !         (tconl(ish)+dcon(ish)*CM%rain%dt))
          !   endif
          !enddo
          !do ish=1,7
          !   tcon(ish,n)=tcon(ish,n)+tconl(ish)
          !   tmass(ish,n)=tmass(ish,n)+tmassl(ish)
          !enddo

    enddo

  end subroutine cal_dmtend_scale_original


  subroutine cal_dmtend_kid(dmtend,tcon,tmass,iproc,&
       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth)
  use par_amps, only : isplit_bin_liq
    real(PS) :: dmtend(11,9,*),tcon(9,*),tmass(9,*)
    integer :: iproc
    integer :: ID(*),JD(*),KD(*),nzp,nxp,nyp,nfpt,nfpth
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in)  :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),&
            zstv(*),zz(*),xx(*),yy(*)

    integer :: n,i,j,ish
    real(PS) :: dvol,dcon(9),tconl(9),tmassl(9)
    real(PS),parameter :: cloudsz=50.0e-4

!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    do n=1,CM%L
       ! +++ loop over grids +++
       ! +++ for rain +++

          dcon=0.0_PS
          tconl=0.0_PS
          tmassl=0.0_PS

          dvol=1.0

          do i=1,CM%rain%N_BIN
             if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XR(rmt_q,i,NRCAT,n)<=1.0e-25.or.XR(rcon_q,i,NRCAT,n)<=1.0e-25)) cycle
!tmp             if(CM%rain%MS(i,n)%len<=cloudsz) then
!tmp                ish=1
!tmp             else
!tmp                ish=2
!tmp             endif
             if(i<=isplit_bin_liq) then
                ish=1
             else
                ish=2
             endif
             ! vapor deposition
             dmtend(1,ish,n)=dmtend(1,ish,n)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,n)=dmtend(2,ish,n)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dmtend(3,ish,n)=dmtend(3,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! CCN activation
             dmtend(4,ish,n)=dmtend(4,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)&
                   *dvol
             !
             ! riming
             dmtend(5,ish,n)=dmtend(5,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,n)=dmtend(6,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,n)=dmtend(7,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,8)&
                   *dvol
             !
             ! ice nucleation (immersion)
             dmtend(8,ish,n)=dmtend(8,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,11)&
                   *dvol
             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,n)=dmtend(9,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,12)&
                   *dvol
             !
             ! ice nucleation (secondary)
             dmtend(10,ish,n)=dmtend(10,ish,n)&
                   +CM%rain%MS(i,n)%dmassdt(1,9)&
                   *dvol

             ! collision-coalescence of liquid hydrometeors
             dcon(ish)=dcon(ish)&
                   +CM%rain%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%rain%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%rain%MS(i,n)%mass(1)*dvol
          enddo

       ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25.or.CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XS(imt_q,i,NSCAT,n)<=1.0e-25.or.XS(icon_q,i,NSCAT,n)<=1.0e-25)) cycle

             ish=CM%solid_hydro%IS(i,n)%sh_type
             ish=ish+2
!             if(ish<=2) then
!                ish=3
!             elseif(ish<=4) then
!                ish=4
!             else
!                ish=5
!             endif

             ! vapor deposition
             dmtend(1,ish,n)=dmtend(1,ish,n)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! evaporation
             dmtend(2,ish,n)=dmtend(2,ish,n)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))&
                   *dvol

             ! aggregation
             dmtend(3,ish,n)=dmtend(3,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,2)&
                   *dvol
             !
             ! ice nucleation (deposition)
             dmtend(4,ish,n)=dmtend(4,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)&
                   *dvol
             !
             ! riming
             dmtend(5,ish,n)=dmtend(5,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)&
                   *dvol
             !
             ! melting-shedding
             dmtend(6,ish,n)=dmtend(6,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,10)&
                   *dvol
             !
             ! ice nucleation (contact)
             dmtend(7,ish,n)=dmtend(7,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)&
                   *dvol
             !
             ! ice nucleation (immersion)
             dmtend(8,ish,n)=dmtend(8,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)&
                   *dvol
             !
             ! ice nucleation (homo freezing)
             dmtend(9,ish,n)=dmtend(9,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)&
                   *dvol
             !
             ! ice nucleation (secondary)
             dmtend(10,ish,n)=dmtend(10,ish,n)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)&
                   *dvol

             ! aggregation
             dcon(ish)=dcon(ish)&
                   +CM%solid_hydro%MS(i,n)%dcondt(2)*dvol

             ! total concentraion and mass
             tconl(ish)=tconl(ish)+CM%solid_hydro%MS(i,n)%con*dvol
             tmassl(ish)=tmassl(ish)+CM%solid_hydro%MS(i,n)%mass(1)*dvol
          enddo

          do ish=1,2
             tconl(8)=tconl(8)+tconl(ish)
             tmassl(8)=tmassl(8)+tmassl(ish)
             dcon(8)=dcon(8)+dcon(ish)
          enddo
          do ish=3,7
             tconl(9)=tconl(9)+tconl(ish)
             tmassl(9)=tmassl(9)+tmassl(ish)
             dcon(9)=dcon(9)+dcon(ish)
          enddo
          do ish=1,9
             if(tconl(ish)>1.0e-20_PS) then
               dmtend(11,ish,n)=dmtend(11,ish,n)+&
                   (-dcon(ish)*tmassl(ish)/&
                   (tconl(ish)+dcon(ish)*CM%rain%dt))
             endif
          enddo
          do ish=1,7
             tcon(ish,n)=tcon(ish,n)+tconl(ish)
             tmass(ish,n)=tmass(ish,n)+tmassl(ish)
          enddo

    enddo

  end subroutine cal_dmtend_kid



  subroutine cal_dmtend_z(dmtend_z,tcon_z,tmass_z,iproc,&
       CM,NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT,&
       dM_auto_liq,dM_accr_liq, &
       XR,XS,KD,ID,JD,zstv,zz,xx,yy,nzp,nxp,nyp,nfpt,nfpth,&
       n1mx,&
       mpirank_w,mpirank_e,mpirank_s,mpirank_n)
    integer :: n1mx
    real(PS) :: dmtend_z(15,n1mx,*),tcon_z(2,n1mx,*),tmass_z(2,n1mx,*)
    integer :: iproc
    integer :: ID(*),JD(*),KD(*),nzp,nxp,nyp,nfpt,nfpth
    type (Cloud_Micro), intent(in)         :: CM
    integer,intent(in)  :: NRTYPE,NRBIN,NRCAT,NSTYPE,NSBIN,NSCAT
    real(PS) :: XR(NRTYPE,NRBIN,NRCAT,*),XS(NSTYPE,NSBIN,NSCAT,*),&
            zstv(*),zz(*),xx(*),yy(*)
    real(PS), dimension(*),intent(in) :: &
                  dM_auto_liq,dM_accr_liq
    integer,intent(in) :: mpirank_w,mpirank_e,mpirank_s,mpirank_n

    integer :: n,i,j,ip,i1,i2,j1,j2

    if(nxp==3) then
      i1=1
      i2=3
    else
      if(mpirank_w==-9) then
        i1=3
      else
        i1=1
      endif
      if(mpirank_e==-9) then
        i2=nxp-3
      else
        i2=nxp
      endif
    endif
    if(mpirank_s==-9) then
      j1=3
    else
      j1=1
    endif
    if(mpirank_n==-9) then
      j2=nyp-3
    else
      j2=nyp
    endif
!    write(fid_alog,*) "nfpt,nfpth",nzp,nxp,nyp,nfpt,nfpth

    j=iproc+1
    do n=1,CM%L
       ! +++ loop over grids +++
       if( (ID(n)<=i1.or.ID(n)>=i2).or.(JD(n)<=j1.or.JD(n)>=j2).or.&
           (KD(n)>=nzp-nfpt+1)) cycle


       ! +++ for rain +++
          do i=1,CM%rain%N_BIN
             if((CM%rain%MS(i,n)%con<=1.0e-25.or.CM%rain%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XR(1,i,NRCAT,n)<=1.0e-25.or.XR(2,i,NRCAT,n)<=1.0e-25)) cycle

             ! vapor deposition
             ip=1
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +max(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))/CM%air%TV(n)%den

             ! evaporation
             ip=2
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +min(0.0_DS,CM%rain%MS(i,n)%dmassdt(1,1))/CM%air%TV(n)%den

             !
             ! CCN activation
             ip=3
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%rain%MS(i,n)%dmassdt(1,3)/CM%air%TV(n)%den

             !
             ! melting-shedding
             ip=4
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%rain%MS(i,n)%dmassdt(1,10)/CM%air%TV(n)%den

             !
             ! autoconversion rate
             ip=5
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +dM_auto_liq(n)/CM%air%TV(n)%den

             !
             ! accretion rate
             ip=6
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +dM_accr_liq(n)/CM%air%TV(n)%den

             tcon_z(1,KD(n),j)=tcon_z(1,KD(n),j)&
                   +CM%rain%MS(i,n)%con/CM%air%TV(n)%den
             tmass_z(1,KD(n),j)=tmass_z(1,KD(n),j)&
                   +CM%rain%MS(i,n)%mass(1)/CM%air%TV(n)%den
          enddo

       ! +++ for ice +++
          do i=1,CM%solid_hydro%N_BIN
             if((CM%solid_hydro%MS(i,n)%con<=1.0e-25.or.CM%solid_hydro%MS(i,n)%mass(1)<=1.0e-25).and.&
                (XS(1,i,NSCAT,n)<=1.0e-25.or.XS(2,i,NSCAT,n)<=1.0e-25)) cycle

             ! vapor deposition
             ip=7
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +max(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))/CM%air%TV(n)%den

             ! evaporation
             ip=8
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +min(0.0_DS,CM%solid_hydro%MS(i,n)%dmassdt(1,1))/CM%air%TV(n)%den

             !
             ! ice nucleation (deposition)
             ip=9
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,7)/CM%air%TV(n)%den
             !
             ! riming
             ip=10
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,4)/CM%air%TV(n)%den
             !
             !
             ! ice nucleation (contact)
             ip=11
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,8)/CM%air%TV(n)%den

             !
             ! ice nucleation (immersion)
             ip=12
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,11)/CM%air%TV(n)%den
             !
             ! ice nucleation (homo freezing)
             ip=13
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,12)/CM%air%TV(n)%den
             !
             ! ice nucleation (secondary)
             ip=14
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,9)/CM%air%TV(n)%den

             ! DHF freezing
             ip=15
             dmtend_z(ip,KD(n),j)=dmtend_z(ip,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%dmassdt(1,3)/CM%air%TV(n)%den

             tcon_z(2,KD(n),j)=tcon_z(2,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%con/CM%air%TV(n)%den
             tmass_z(2,KD(n),j)=tmass_z(2,KD(n),j)&
                   +CM%solid_hydro%MS(i,n)%mass(1)/CM%air%TV(n)%den

          enddo

    enddo

  end subroutine cal_dmtend_z

end MODULE class_Cloud_Micro
