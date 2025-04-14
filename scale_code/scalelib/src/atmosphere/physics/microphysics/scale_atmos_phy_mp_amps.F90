!-------------------------------------------------------------------------------
!> module Advanced Microphysics Prediction System (AMPS)
!!
!! @par Description:
!!      This module contains subroutines for the AMPS model
!!
!! - Reference
!!  - Hashino and Tripoli (2007;JAS) doi:10.1175/JAS39631
!!  - Ong et al. (2022;JAMES) doi:10.1029/2021MS002887
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_mp_amps
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prof

#ifndef DISABLE_AMPS
  use maxdims, only: mxnbin,npamx,nbamx,ncamx
  use scale_const, only: &
       EPS    => CONST_EPS, &
       PI     => CONST_PI,  &
       CPdry  => CONST_CPdry,  &
       Rdry   => CONST_Rdry,   &
       GRAV   => CONST_GRAV,   &
       LHV0   => CONST_LHV0,   &
       LHS0   => CONST_LHS0,   &
       PRE00  => CONST_PRE00,  &
       EPSvap => CONST_EPSvap

  use class_Cloud_Micro, only: Cloud_Micro
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_MP_amps_tracer_setup
  public :: ATMOS_PHY_MP_amps_setup
  public :: ATMOS_PHY_MP_amps_finalize
  public :: ATMOS_PHY_MP_amps_tendency
  public :: ATMOS_PHY_MP_amps_cloud_fraction
  public :: ATMOS_PHY_MP_amps_effective_radius
  public :: ATMOS_PHY_MP_amps_qtrc2qhyd
  public :: ATMOS_PHY_MP_amps_qtrc2nhyd
  public :: ATMOS_PHY_MP_amps_qhyd2qtrc
  public :: ATMOS_PHY_MP_amps_init_qtrc_1DMODEL
  public :: ATMOS_PHY_MP_amps_init_qtrc_BOX

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  integer, public :: ATMOS_PHY_MP_amps_ntracers
  integer, public :: ATMOS_PHY_MP_amps_nmasses
  integer, public :: ATMOS_PHY_MP_amps_nwaters
  integer, public :: ATMOS_PHY_MP_amps_nices
  integer, public :: ATMOS_PHY_MP_amps_nfrz
  integer, public :: ATMOS_PHY_MP_amps_nccn
  integer, public :: ATMOS_PHY_MP_amps_nppv

  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_MP_amps_tracer_names(:)
  character(len=H_MID)  , public, allocatable :: ATMOS_PHY_MP_amps_tracer_descriptions(:)
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_MP_amps_tracer_units(:)
  integer,                public, allocatable :: ATMOS_PHY_MP_amps_tracer_ref(:)

  character(len=H_SHORT), public :: AMPS_mt_NAME(20), & !< name  of the variables
                                    AMPS_ct_NAME(20), & !< name  of the variables
                                    AMPS_bt_NAME(3)     !< name  of the variables

  data AMPS_mt_NAME / 'ampsl_vapor', 'ampsl_eva', 'ampsl_col', 'ampsl_act', &
                      'ampsl_rim', 'ampsl_melt', 'ampsl_cont', 'ampsl_imm', &
                      'ampsl_homo', 'ampsl_hmproc', &
                      'ampsi_vapor', 'ampsi_eva', 'ampsi_col', 'ampsi_ddf', &
                      'ampsi_rim', 'ampsi_mlsh', 'ampsi_cont', 'ampsi_imm', &
                      'ampsi_homo', 'ampsi_hmproc' /

  data AMPS_ct_NAME / 'ampscl_vapor', 'ampscl_eva', 'ampscl_col', 'ampscl_act', &
                      'ampscl_rim', 'ampscl_melt', 'ampscl_cont', 'ampscl_imm', &
                      'ampscl_homo', 'ampscl_hmproc', &
                      'ampsci_vapor', 'ampsci_eva', 'ampsci_col', 'ampsci_ddf', &
                      'ampsci_rim', 'ampsci_mlsh', 'ampsci_cont', 'ampsci_imm', &
                      'ampsci_homo', 'ampsci_hmproc' /

  data AMPS_bt_NAME / 'riming_bin', 'vapor_bin', 'eva_bin' /

  character(len=H_MID),   public :: AMPS_t_DESC(20) !< desc. of the variables

  data AMPS_t_DESC / 'liq vapor deposition tendency', 'liq evaporation tendency', &
                     'liq collection tendency', 'liq CCN activation tendency', &
                     'liq riming tendency', 'liq melting tendency', &
                     'liq contact tendency', 'liq immersion tendency', &
                     'liq homogeneous f tendency', 'liq Hallet Mossop tendency', &
                     'ice vapor deposition tendency', 'ice evaporation tendency', &
                     'ice collection tendency', 'ice deliques depo tendency', &
                     'ice riming tendency', 'ice melting-shedding tendency', &
                     'ice contact tendency', 'ice immersion tendency', &
                     'ice homogeneous f tendency', 'ice Hallet Mossop tendency' /

  character(len=H_SHORT), public :: AMPS_t_UNIT(2) !< unit  of the variables

  data AMPS_t_UNIT / '/g cm3/s', '/cm3/s' /



#ifdef DISABLE_AMPS
  integer, parameter, private :: QA = 2
#else
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters
  !
  integer, parameter :: k1eta = 2

  integer :: QA, I_QV, I_QL, I_QW, I_QI, I_QPPVL, I_QPPVI, I_QPPVA
  integer :: ivis = 1
  !integer, dimension(max_nmoments_liq)  :: I_scl2ship_l
  !integer, dimension(max_nmoments_ice)  :: I_scl2ship_i
  !integer, dimension(max_nmoments_aero) :: I_scl2ship_a

  ! internal time step
  integer :: TIME_AMPS

  logical  :: amps_debug         = .false.
  logical  :: amps_ignore        = .false.
  logical  :: l_fix_aerosols     = .true.
  logical  :: l_sediment         = .true.
  logical  :: l_fill_aerosols    = .false.
  logical  :: l_bin_shift        = .false.
  logical  :: l_axis_limit       = .true.
  integer  :: ini_aerosol_prf    = 3
  integer  :: l_gaxis_version    = 1
  integer  :: l_aadv_version     = 2
  integer  :: l_reff_version     = 2
  integer  :: iadvv              = 1

  integer  :: nx, ny, nz, nzh
  integer, parameter :: max_nmoments_liq=4, max_nmoments_ice=16, max_nmoments_aero=3
  integer, parameter :: nspecies=2, max_char_len=200
  integer, parameter :: naerosol=4 ! number of aerosol species
  integer, parameter  :: num_aero_bins(naerosol) = 1
  integer, parameter  :: num_aero_moments(naerosol) = 3
  integer, parameter  :: num_h_moments(nspecies) = (/  &
              4 & ! liquid
              ,16 & ! ice
              /)
  integer  :: num_h_bins(nspecies)= (/ &
  ! for fine
  !         80 & ! liquid
  !        ,83 &  ! ice
  ! for midium
  !         60 & ! liquid
  !        ,40 &  ! ice
  ! for coarse (default)
            40 & ! liquid
           ,20 &  ! ice
           /)

  character(len=9) :: mom_names_ice(max_nmoments_ice)= &
       (/ 'ice_imass'   &                        ! ice total mass
       ,  'rim_imass'   &                        ! rime mass
       ,  'agg_imass'   &                        ! aggragate mass
       ,  'cry_imass'   &                        ! crystal mass
       ,  'mtw_imass'   &                        ! melt water mass
       ,  'frw_imass'   &                        ! frozen water mass (nucleation)
       ,  'tae_imass'   &                        ! total aerosol mass
       ,  'sat_imass'   &                        ! soluble aerosol mass
       ,  'con_imass'   &                        ! number concentration
       ,  'vol_imass'   &                        ! circumscribing volume
       ,  'a_axis   '   &                        ! a-axis length
       ,  'c_axis   '   &                        ! c-axis length
       ,  'd_axis   '   &                        ! d-axis length
       ,  'ag_axis  '   &                        ! center of gravity a (polycrystals)
       ,  'cg_axis  '   &                        ! center of gravity c (polycrystals)
       ,  'ex_cry   '   &                        ! extra crystalline structure
       /)

  character(len=34) :: mom_descriptions_ice(max_nmoments_ice)= &
       (/ 'ice total mass                    '    &
       ,  'rime mass                         '    &
       ,  'aggragate mass                    '    &
       ,  'crystal mass                      '    &
       ,  'melt water mass                   '    &
       ,  'frozen water mass (nucleation)    '    &
       ,  'total aerosol mass                '    &
       ,  'soluble aerosol mass              '    &
       ,  'number concentration              '    &
       ,  'circumscribing volume             '    &
       ,  'a-axis length                     '    &
       ,  'c-axis length                     '    &
       ,  'd-axis length                     '    &
       ,  'center of gravity a (polycrystals)'    &
       ,  'center of gravity c (polycrystals)'    &
       ,  'extra crystalline structure       '    &
       /)

  character(len=6) :: mom_units_ice(max_nmoments_ice)= &
       (/ 'kg/kg '   &                        ! ice total mass
       ,  'kg/kg '   &                        ! rime mass
       ,  'kg/kg '   &                        ! aggregate mass
       ,  'kg/kg '   &                        ! crystal mass
       ,  'kg/kg '   &                        ! melt water mass
       ,  'kg/kg '   &                        ! frozen water mass (nucleation)
       ,  'kg/kg '   &                        ! total aerosol mass
       ,  'kg/kg '   &                        ! soluble aerosol mass
       ,  '1/kg  '   &                        ! number concentration
       ,  'cm3/kg'   &                        ! circumscribing volume
       ,  '1/kg  '   &                        ! a-axis length
       ,  '1/kg  '   &                        ! c-axis length
       ,  '1/kg  '   &                        ! d-axis length
       ,  '1/kg  '   &                        ! center of gravity a (polycrystals)
       ,  '1/kg  '   &                        ! center of gravity c (polycrystals)
       ,  '1/kg  '   &                        ! extra crystalline structure
       /)

  character(len=9) :: mom_names_liq(max_nmoments_liq)= &
       (/ 'liq_lmass'   &                        ! liquid total mass
       ,  'tae_lmass'   &                        ! total aerosol mass
       ,  'sae_lmass'   &                        ! soluble aerosol mass
       ,  'con_lmass'   &                        ! number concentration
       /)

  character(len=20) :: mom_descriptions_liq(max_nmoments_liq)= &
       (/ 'liquid total mass   '   &
       ,  'total aerosol mass  '   &
       ,  'soluble aerosol mass'   &
       ,  'number concentration'   &
       /)

  character(len=5) :: mom_units_liq(max_nmoments_liq)= &
       (/ 'kg/kg'   &                        ! liquid total mass, rho_r/rho_total
       ,  'kg/kg'   &                        ! total aerosol mass
       ,  'kg/kg'   &                        ! soluble aerosol mass
       ,  '1/kg '   &                        ! number concentration
       /)

  character(len=4) :: aero_names(naerosol)= &
       (/ 'ccn1'   &
       ,  'in1 '   &
       ,  'ccn2'   &
       ,  'ccn3'   &
       /)

  character(len=9) :: mom_names_aero(max_nmoments_aero)= &
       (/ 'ae_amass '   &
       ,  'con_amass'   &
       ,  'sae_amass'   &
       /)

  character(len=20) :: mom_descriptions_aero(max_nmoments_aero)= &
       (/ 'total aerosol mass  '  &
       ,  'number concentration'  &
       ,  'soluble aerosol mass'  &
       /)

  character(len=5) :: mom_units_aero(max_nmoments_aero)= &
      (/  'kg/kg'   &
       ,  '1/kg '   &
       ,  'kg/kg'   &
       /)

  character(len=32),dimension(20) :: proname_ifc,proname_mtd

  ! for liquid category
  integer :: split_bins = 21
  ! for ice category
  integer :: ibin_is
  integer :: ibin_sg
  integer :: ibin_gh


  ! table for saturation vapor pressure over super cooled liquid and ice
  real(DP), dimension(150) ::  estbar
  real(DP), dimension(111) ::  esitbar

  ! lagrangian parabolic advection method
  real(RP), allocatable :: ADVPPMZ(:,:,:,:), ADVPPMZE(:,:,:)

  ! indices about the moment, bin, and categories
  integer :: npr       ! number of moments of a liquid bin
  integer :: nbr       ! number of liquid bins
  integer :: ncr       ! number of categories for liquid, always = 1
  integer :: npi       ! number of moments of a ice bin
  integer :: nbi       ! number of ice bins
  integer :: nci       ! number of categories for ice, always = 1
  integer :: npa       ! number of moments of a aerosols bin
  integer :: nba       ! number of aerosols bins
  integer :: nca       ! number of categories for aerosols

  ! number of liquid and ice PPV stored in SCALE
  integer :: numberPPVL
  integer :: numberPPVI

  ! local temperory variables for level of complexity
  integer :: level

  ! for PPM sedimentation scheme
  integer :: ngrid

  ! local variable, number of bins of haze, set in binmicrosetup_scale (after config_varindex)
  integer :: nbhzcl

  ! switches
  integer, parameter :: max_cpu=16
  logical, save :: micro_io_strt(max_cpu)=.True.

  ! intial aerosol profiles, 1D -> 3D
  real(RP), dimension(:,:,:,:,:,:), allocatable :: qapv_ini

  ! intial density profiles, 1D -> 3D
  real(RP), dimension(:,:,:), allocatable :: den_ini

  ! mean mass, used in bin shift calculation after dynamics advection and after sedimentation
  real(RP),dimension(:,:,:,:,:), allocatable :: mmassrv_global
  real(RP),dimension(:,:,:,:,:), allocatable :: mmassiv_global

  !type (Cloud_Micro)  :: CM
  type (Cloud_Micro), dimension(:), allocatable :: CM

#endif
  !----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Config
  subroutine ATMOS_PHY_MP_amps_tracer_setup
#ifndef DISABLE_AMPS
    use mod_amps_utility, only: &
       config_varindex
    use com_amps
    use par_amps
#endif
    use scale_prc, only: &
       PRC_abort, &
       PRC_IsMaster
    implicit none


#ifdef DISABLE_AMPS
    LOG_ERROR("ATMOS_PHY_MP_amps_tracer_setup",*) "Scalelib was compiles without AMPS"
    LOG_ERROR_CONT(*) "Recompile scalelib with SCALE_DISABLE_AMPS=F"
    call PRC_abort
#else
    namelist / PARAM_ATMOS_PHY_MP_AMPS_bin / &
       num_h_bins,         & ! 40 or 80 bins for liquid, ONLY 20 bins for ice
       nbin_h,             & ! number of bins for haze particles, e.x. 20 for 40 liq. bins
       iadvv,              & ! 1 for Euler, 2 for PPM sedimentation scheme
       ini_aerosol_prf,    & ! 1 for SHEBA, 2 for MPACE, 3 for general use (1 and 2 OVERRIDE aerosol settings in AMPSTASK.F)
       l_fix_aerosols,     & ! invariant aerosols throughout integration
       l_sediment,         & ! sediment on or off
       l_fill_aerosols,    & ! fill aerosols in cloud-free region or not
       l_bin_shift,        & ! whether bin shift is performed after advection, default is false for testing phase, it is still under checking
       l_gaxis_version,    & ! axis definition, 1 (a, c, d, ag, cg), 2 (a, c/a, d/a, ag, cg/ag), 3 (a, c, d, ag/a, cg/a), default is 1
       l_aadv_version,     & ! 1 for all PPVs advected in the same fasion, 2 for modified advection for non-mass PPVs
       l_axis_limit,       & ! whether center of gravity axis (ag and cg) limit=1 is applied, default is true
       l_reff_version,     & ! using (1) maximum dimension or (2) equivalent spherical radius for calculating radiation
       amps_debug,         & ! debugging on or off
       amps_ignore           ! ignore amps microphysics or not

    integer :: i, m, n, o, ierr

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*) 'Tracers setup for SHIPS (2007) Spectral BIN model'

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*) 'READ BIN NUMBER'

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_AMPS_bin,iostat=ierr)

    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*)  'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_amps_tracer_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_AMPS_bin, Check!'
       call PRC_abort
    end if

    LOG_NML(PARAM_ATMOS_PHY_MP_AMPS_bin)

    level = 5 ! default is 5, for future broader use

    ! liquid
    ! parameters
    npr = num_h_moments(1) + 2   ! two terminal velocity variables
    ! bin numbers
    nbr = num_h_bins(1)
    ! category numbers
    ncr = 1

    ! ice
    ! parameters
    npi = num_h_moments(2)+2   ! two terminal velocity variables
    ! bin numbers
    nbi = num_h_bins(2)
    ! category numbers
    nci = 1

    ! aerosol
    ! parameters
    npa = num_aero_moments(1)+2   ! two terminal velocity variables
    ! bin numbers
    nba = num_aero_bins(1)
    ! category numbers
    !   by default, category 1, 3, 4 are CCN
    !   by default, category 2 IN

    ! set up the constants for the distribution of aerosols
    select case(ini_aerosol_prf)
    case(1) ! SHEBA
      nca = 3
    case(2) ! MPACE
      nca = 3
    case(3) ! general three-type aerosol case (fine CCN, IN, coarse CCN)
      nca = 3
    case default
      nca = 4
    end select

    LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*) "flagp_c,flagp_r,flagp_s,flagp_a:",flagp_c,flagp_r,flagp_s,flagp_a
    LOG_INFO_CONT(*) "npr,nbr,ncr,npi,nbi,nci,npa,nba,nca:",npr,nbr,ncr,npi,nbi,nci,npa,nba,nca

    !
    !  Configure the variable indexes
    !
    call config_varindex

    ! m is number of liquid bin, n is number of ice bin
    !
    !        total rain mass            melt water mass         rimed mass    ... agg, crystals, ...
    !  | v | l1 | l2 | ...... | lm | w1 | w2 | ...... | wn | r1 | r2 | ...... | rn |   ......
    !      frozen water mass       other liquid PPVs        other ice PPVs
    !  | f1 | f2 | ...... | fn | p1 | p2 | ...... | px | p1 | p2 | ...... | px |
    !        aerosol PPVs
    !  | p1 | p2 | ...... | px |
    ! vapor|   no. of liquid bins  |    no. of ice bins    | no. of PPV X no of bins, liquid -> ice
    !
    !  ccn are then put after the last PPV is defined above

    ! no. of total liquid and ice mass tracers
    ! 1 = vapor mixing ratio
    ATMOS_PHY_MP_amps_nmasses  = 1 + num_h_bins(1) + 2*num_h_bins(2)
    ! number of all tracers
    ATMOS_PHY_MP_amps_ntracers = 1 + num_h_bins(1)*num_h_moments(1) &
                                   + num_h_bins(2)*num_h_moments(2) &
                                   + num_aero_bins(1)*num_aero_moments(1)*nca
    ! number of liquid water = rain water + melt water
    ATMOS_PHY_MP_amps_nwaters  = num_h_bins(1) + num_h_bins(2)
    ! number of ice water = rimed + aggregate + crystals [+ frozen water (no more)]
    ATMOS_PHY_MP_amps_nices    = num_h_bins(2)
    ! number of frozen ice
    ATMOS_PHY_MP_amps_nfrz     = num_h_bins(2)
    ! number of nonmass PPVs and aerosol PPVs = total and soluble aerosol mass, no. con. (liquid) +
    ! total mass, total and soluble aerosol mass, no. con., volume and length related PPVs (ice)  +
    ! CCN PPVs
    ATMOS_PHY_MP_amps_nppv     = num_h_bins(1)*(num_h_moments(1) - 1) + &
                                 num_h_bins(2)*(num_h_moments(2) - 2) + &
                                 num_aero_bins(1)*num_aero_moments(1)*nca

    I_QV    = 1                                               ! starting index of vapor
    I_QL    = 2                                               ! starting index of liquid water
    I_QW    = I_QL    + num_h_bins(1)                         ! starting index of melt water
    I_QI    = I_QW    + num_h_bins(2)                         ! starting index of ice water
    I_QPPVL = I_QI    + num_h_bins(2)                         ! starting index of liquid PPVs
    I_QPPVI = I_QPPVL + num_h_bins(1)*(num_h_moments(1) - 1)  ! starting index of ice PPVs
    I_QPPVA = I_QPPVI + num_h_bins(2)*(num_h_moments(2) - 2)  ! starting index of aerosol PPVs

    numberPPVL = num_h_moments(1) - 1
    numberPPVI = num_h_moments(2) - 2

    allocate( ATMOS_PHY_MP_amps_tracer_names       (ATMOS_PHY_MP_amps_ntracers) )
    allocate( ATMOS_PHY_MP_amps_tracer_descriptions(ATMOS_PHY_MP_amps_ntracers) )
    allocate( ATMOS_PHY_MP_amps_tracer_units       (ATMOS_PHY_MP_amps_ntracers) )
    allocate( ATMOS_PHY_MP_amps_tracer_ref         (ATMOS_PHY_MP_amps_ntracers) )

    LOG_INFO("ATMOS_PHY_MP_amps_tracer_setup",*) "----------AMPS SCALE variable indices-----------"
    LOG_INFO_CONT(*) "I_QV                       : ", I_QV
    LOG_INFO_CONT(*) "I_QL                       : ", I_QL
    LOG_INFO_CONT(*) "I_QW                       : ", I_QW
    LOG_INFO_CONT(*) "I_QI                       : ", I_QI
    LOG_INFO_CONT(*) "I_QPPVL                    : ", I_QPPVL
    LOG_INFO_CONT(*) "I_QPPVI                    : ", I_QPPVI
    LOG_INFO_CONT(*) "I_QPPVA                    : ", I_QPPVA
    LOG_INFO_CONT(*) "ATMOS_PHY_MP_amps_nmasses  : ", ATMOS_PHY_MP_amps_nmasses
    LOG_INFO_CONT(*) "ATMOS_PHY_MP_amps_ntracers : ", ATMOS_PHY_MP_amps_ntracers
    LOG_INFO_CONT(*) "ATMOS_PHY_MP_amps_nwaters  : ", ATMOS_PHY_MP_amps_nwaters
    LOG_INFO_CONT(*) "ATMOS_PHY_MP_amps_nices    : ", ATMOS_PHY_MP_amps_nices
    LOG_INFO_CONT(*) "ATMOS_PHY_MP_amps_nppv     : ", ATMOS_PHY_MP_amps_nppv
    LOG_INFO_CONT(*) "precision                  : ", DP, RP
    LOG_INFO_CONT(*) "-------end of AMPS SCALE variable indices-------"

    LOG_INFO_CONT(*) ""

    LOG_INFO_CONT(*) "----------AMPS water variable indices-----------"
    LOG_INFO_CONT(*) "rmt_q                      : ", rmt_q
    LOG_INFO_CONT(*) "rmat_q                     : ", rmat_q
    LOG_INFO_CONT(*) "rmas_q                     : ", rmas_q
    LOG_INFO_CONT(*) "rcon_q                     : ", rcon_q
    LOG_INFO_CONT(*) "-------end of AMPS water variable indices-------"

    LOG_INFO_CONT(*) ""

    LOG_INFO_CONT(*) "-----------AMPS ice variable indices------------"
    LOG_INFO_CONT(*) "imt_q                      : ", imt_q
    LOG_INFO_CONT(*) "imr_q                      : ", imr_q
    LOG_INFO_CONT(*) "ima_q                      : ", ima_q
    LOG_INFO_CONT(*) "imc_q                      : ", imc_q
    LOG_INFO_CONT(*) "imw_q                      : ", imw_q
    LOG_INFO_CONT(*) "imf_q                      : ", imf_q
    LOG_INFO_CONT(*) "imat_q                     : ", imat_q
    LOG_INFO_CONT(*) "imas_q                     : ", imas_q
    LOG_INFO_CONT(*) "icon_q                     : ", icon_q
    LOG_INFO_CONT(*) "ivcs_q                     : ", ivcs_q
    LOG_INFO_CONT(*) "iacr_q                     : ", iacr_q
    LOG_INFO_CONT(*) "iccr_q                     : ", iccr_q
    LOG_INFO_CONT(*) "idcr_q                     : ", idcr_q
    LOG_INFO_CONT(*) "iag_q                      : ", iag_q
    LOG_INFO_CONT(*) "icg_q                      : ", icg_q
    LOG_INFO_CONT(*) "inex_q                     : ", inex_q
    LOG_INFO_CONT(*) "--------end of AMPS ice variable indices--------"


    !
    !  Synchorize the variable indices in AMPS and indices in SCALE
    !
    !I_scl2ship_l(1) = I_QL       ! rmt_q
    !I_scl2ship_l(2) = I_QPPVL    ! rmat_q
    !I_scl2ship_l(3) = I_QPPVL+2  ! rmas_q
    !I_scl2ship_l(4) = I_QPPVL+3  ! rcon_q

    !I_scl2ship_i(1) = I_QPPVI    ! imt_q
    !I_scl2ship_i(2) = I_QI       ! imr_q
    !I_scl2ship_i(3) = I_QPPVL+2  ! ima_q
    !I_scl2ship_i(4) = I_QPPVL+3  ! imc_q
    !I_scl2ship_i(5) = I_QW       ! imw_q
    !I_scl2ship_i(6) = I_QPPVL+3  ! imf_q

    !---------------------------------------------------------------------------
    !
    !++ calculate each category and aerosol
    !
    !---------------------------------------------------------------------------

    ATMOS_PHY_MP_amps_tracer_names(1) = 'QV'
    ATMOS_PHY_MP_amps_tracer_units(1) = 'kg/kg'
    ATMOS_PHY_MP_amps_tracer_descriptions(1) = 'Water Vapor mixing ratio'
    ATMOS_PHY_MP_amps_tracer_ref(1) = 0

    ! liquid
    i = I_QL
    do n = 1, num_h_bins(1)
       ! total rain mass (without aerosol)
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_liq(rmt_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_liq(rmt_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_liq(rmt_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1
    enddo
    i = I_QW
    do n = 1, num_h_bins(2)
       ! melt water mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imw_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imw_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imw_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1
    enddo

    ! ice
    i = I_QI
    do n = 1, num_h_bins(2)
       ! total ice mass (without melt water and aerosol)
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imt_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imt_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imt_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1
    enddo

    ! other liquid PPVs
    i = I_QPPVL
    do n = 1, num_h_bins(1)
       ! total aerosol mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_liq(rmat_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_liq(rmat_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_liq(rmat_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! soluble aerosol mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_liq(rmas_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_liq(rmas_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_liq(rmas_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! number concentration
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_liq(rcon_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_liq(rcon_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_liq(rcon_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1
    enddo


    ! other ice PPVs
    i = I_QPPVI
    do n = 1, num_h_bins(2)
       ! 1. rimed mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imr_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imr_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imr_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 2. aggregate mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(ima_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(ima_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(ima_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 3. crystal mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imc_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imc_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imc_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 4. frozen water mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imf_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imf_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imf_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 5. total aerosol mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imat_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imat_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imat_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 6. soluble aerosol mass
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(imas_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(imas_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(imas_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 7. number concentration
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(icon_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(icon_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(icon_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = 0
       i = i + 1

       ! 8. circumscribing volume
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(ivcs_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(ivcs_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(ivcs_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -1
       i = i + 1

       ! 9. a axis length
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(iacr_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(iacr_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(iacr_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -2
       i = i + 1

       ! 10. c axis length
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(iccr_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(iccr_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(iccr_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -3
       i = i + 1

       ! 11. d axis length
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(idcr_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(idcr_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(idcr_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -4
       i = i + 1

       ! 12. a center of gravity
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(iag_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(iag_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(iag_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -5
       i = i + 1

       ! 13. c center of gravity
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(icg_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(icg_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(icg_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -6
       i = i + 1

       ! 14. crystalline structure
       ATMOS_PHY_MP_amps_tracer_units(i) = mom_units_ice(inex_q)
       write(ATMOS_PHY_MP_amps_tracer_names(i),       '(a,i0)')  trim(mom_names_ice(inex_q)), n
       write(ATMOS_PHY_MP_amps_tracer_descriptions(i),'(a,i0)')  trim(mom_descriptions_ice(inex_q)), n
       ATMOS_PHY_MP_amps_tracer_ref(i) = -7
       i = i + 1
    enddo

    ! aerosol
    o = I_QPPVA
    do i = 1, nca
      do n = 1, num_aero_bins(1)
        do m = 1, num_aero_moments(1)
           ATMOS_PHY_MP_amps_tracer_units(o) = mom_units_aero(m)
           write(ATMOS_PHY_MP_amps_tracer_names(o),       '(a,2i0)')  trim(mom_names_aero(m)), n, i
           write(ATMOS_PHY_MP_amps_tracer_descriptions(o),'(a,2i0)')  trim(mom_descriptions_aero(m)), n, i
           ATMOS_PHY_MP_amps_tracer_ref(o) = 0
           o = o + 1
        enddo
      enddo
    enddo

    if (l_aadv_version == 1) then
       ATMOS_PHY_MP_amps_tracer_ref(:) = 0
    endif

    QA = ATMOS_PHY_MP_amps_ntracers

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_MP_amps_setup(KA,KS,KE,IA,IS,IE,JA,JS,JE)
#ifndef DISABLE_AMPS
    use mod_amps_utility, only: &
       qsparm2
    use mod_amps_lib, only:  &
       ini_aerosols_mpace, &
       ini_aerosols_const, &
       ini_aerosols_sheba
    use mod_scladv_slppm, only: &
       cal_ppmcoef, &
       cal_ppmcoef_eta
    use scale_const, only: &
       CONST_THERMODYN_TYPE
    ! real_z is the coordinates transformed from no-topograpy (flat) to topography
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
       REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    ! this is not real lat_lon cooedinates, this is the x-y coordinates of a rectangular grid
    use scale_atmos_grid_cartesC, only: &
       REAL_FX => ATMOS_GRID_CARTESC_FX, &
       REAL_FY => ATMOS_GRID_CARTESC_FY
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_myrank, &
       PRC_abort
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE, IA, IS, IE, JA, JS, JE

#ifndef DISABLE_AMPS
    integer  :: i, j, k


    ! check thermodyn type
    if ( CONST_THERMODYN_TYPE .ne. "SIMPLE" ) then
       LOG_ERROR("ATMOS_PHY_MP_amps_setup",*) "CONST_THERMODYN_TYPE must be 'SIMPLE'!", trim(CONST_THERMODYN_TYPE)
       call PRC_abort
    end if

    ! compute the tables of supersaturated pressure over ice and supercooled water
    call qsparm2(estbar,esitbar)

    ! initialize bin structure, variable indices, coefficients, and lookup tables
    call init_AMPS(KA,          &    ! [IN]
                   KS,          &    ! [IN]
                   KE,          &    ! [IN]
                   IA,          &    ! [IN]
                   IS,          &    ! [IN]
                   IE,          &    ! [IN]
                   JA,          &    ! [IN]
                   JS,          &    ! [IN]
                   JE,          &    ! [IN]
                   REAL_CZ,     &    ! [IN]
                   REAL_FZ,     &    ! [IN]
                   REAL_FX,     &    ! [IN]
                   REAL_FY)          ! [IN]


    ! set up the constants for the distribution of aerosols
    select case(ini_aerosol_prf)
    case(1) ! SHEBA
      call ini_aerosols_sheba(nca)
    case(2) ! MPACE
      call ini_aerosols_mpace(nca)
    case(3) ! general three-type aerosol case (fine CCN, IN, coarse CCN)
      call ini_aerosols_const(nca)
    case default ! under constrcution, do not use
      call ini_aerosols_const(nca)
    end select

    ! initialize parabolic advection method
    !$omp parallel do
    do j = JS, JE
    do i = IS, IE
       call cal_ppmcoef(2,nzh,nzh,REAL_FZ(KS-1:KE,i,j),ADVPPMZ(:,:,i,j))
       call cal_ppmcoef_eta(nzh,REAL_FZ(KS-1:KE,i,j),k1eta,ADVPPMZE(:,i,j))
    end do
    end do

    TIME_AMPS = 0

    ! initialize global microphysical variables, for bin shift after advection
    allocate(mmassrv_global(nbr,ncr,KS-1:KE,IS:IE,JS:JE) &
            ,mmassiv_global(nbi,nci,KS-1:KE,IS:IE,JS:JE))

    ! initialize initial profiles of aerosol distribution
    allocate(qapv_ini(npamx,nbamx,ncamx,KS-1:KE,IS:IE,JS:JE))
    ! initialize 1D initial density profile
    allocate(den_ini(KS-1:KE,IS:IE,JS:JE))

    !$omp workshare
    mmassrv_global(:,:,:,:,:) = 0.0_RP
    mmassiv_global(:,:,:,:,:) = 0.0_RP
    qapv_ini(:,:,:,:,:,:) = 0.0_RP
    den_ini(:,:,:) = 0.0_RP
    !$omp end workshare

    if(level.eq.5) then
       proname_ifc(1)="ini_cloud_micro"
       proname_ifc(2)="reality_check"
       proname_ifc(3)="check_water_apmass 1"
       proname_ifc(4)="cal_micro_tendency"
       proname_ifc(5)="print_cloud_tendency"
       proname_ifc(6)="return_output"
       proname_ifc(7)="check_water_apmass 2"
       proname_ifc(8)="update_terminal_vel"
       proname_ifc(9)="cal_dmtend"
       proname_mtd(1)="diag-thermo"
       proname_mtd(2)="diag-solid"
       proname_mtd(3)="mv ice2liq"
       proname_mtd(4)="diag-rain"
       proname_mtd(5)="diag-aerosols"
       proname_mtd(6)="rain-rain"
       proname_mtd(7)="auto-conv"
       proname_mtd(8)="ice_ice"
       proname_mtd(9)="ice_rain"
       proname_mtd(10)="hyd-bre rain"
       proname_mtd(11)="hyd-bre ice"
       proname_mtd(12)="repair"
       proname_mtd(13)="ice nuc1"
       proname_mtd(14)="CCN activation"
       proname_mtd(15)="vap-dep rain"
       proname_mtd(16)="vap-dep ice"
       proname_mtd(17)="ice nuc2"
       proname_mtd(18)="melting"
       proname_mtd(19)="update_group"
       proname_mtd(20)="print out"
    endif

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_setup


  subroutine init_AMPS(KA,KS,KE,IA,IS,IE,JA,JS,JE,REAL_CZ,REAL_FZ,REAL_FX,REAL_FY)
#ifndef DISABLE_AMPS
    use mod_amps_utility, only: &
       read_seed
    use mod_amps_lib, only: &
       binmicrosetup_scale
    use par_amps, only: &
       isplit_bin_liq, &
       isplit_bin_ice
    use scale_prc, only: &
        PRC_abort, &
        PRC_IsMaster
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE, IA, IS, IE, JA, JS, JE
    real(RP), dimension(KA,IA,JA),   intent(in) :: REAL_CZ
    real(RP), dimension(0:KA,IA,JA), intent(in) :: REAL_FZ
    real(RP), dimension(0:IA), intent(in) :: REAL_FX
    real(RP), dimension(0:JA), intent(in) :: REAL_FY

#ifndef DISABLE_AMPS
    ! local variables
    integer :: ier, i, j, k


    ! grid info
    nz = KE - KS + 1
    nzh = nz + 1
    nx = IE - IS + 1
    ny = JE - JS + 1
    ngrid = 1

    if ( PRC_IsMaster ) then
       LOG_INFO("init_AMPS",*) "nx,ny,nz",nx,ny,nz, IS, IE, JS, JE, KS, KE
    end if


    call binmicrosetup_scale(npr, &        ! [IN]
                             nbr, &        ! [IN]
                             ncr, &        ! [IN]
                             npi, &        ! [IN]
                             nbi, &        ! [IN]
                             nci, &        ! [IN]
                             npa, &        ! [IN]
                             nba, &        ! [IN]
                             nca, &        ! [IN]
                             nz,  &        ! [IN]
                             nx,  &        ! [IN]
                             ny,  &        ! [IN]
                             nbhzcl, &     ! [OUT]
                             estbar, &     ! [IN]
                             esitbar)      ! [IN]

    ! set split bin boundary, com_par
    split_bins = min(isplit_bin_liq, nbr)

    ! category bound (tentative)
    ibin_is = min(isplit_bin_ice, nbi)
    ibin_sg = ibin_is + (nbi-ibin_is) / 3
    ibin_gh = ibin_is + (nbi-ibin_is) * 2 / 3


    ! read seeds for random generator
    call read_seed(0.0_RP)


    allocate(ADVPPMZ(11,KS-1:KE,IS:IE,JS:JE), ADVPPMZE(11,IS:IE,JS:JE), &
            stat=ier)
    if(ier/=0) then
      LOG_ERROR("init_AMPS",*) "Allocation failed for Lagrangian PPM"
      call PRC_abort
    endif

#endif
    return
  end subroutine init_AMPS


!merk
!OCL SERIAL
  subroutine moistthermo2_scale( &
       KS, KE, &
       npr,nbr,ncr,npi,nbi,nci,npa,nba,nca, &
       level,nbhzcl,  &
       i, j,   &
       thp,    &    ! [INOUT]  theta il
       th,     &    ! [INOUT]  potential temperature (from SCALE)
       den,    &    ! [IN]     density for moist air XX(from SCALE)XX
       !                                den0,   &    ! [IN]     density for moist air XX(from SCALE)XX
       t,      &    ! [INOUT]  temperature           (from SCALE)
       !                                t0,     &    ! [IN]     original AMPS temperature
       qtp,    &    ! [INOUT]  total hydrometeors
       qv,     &    ! [INOUT]  vapor mixing ratio
       qc,     &    ! [INOUT]  cloud droplets mixing ratio
       qr,     &    ! [INOUT]  rain droplets mixing ratio
       qi,     &    ! [INOUT]  ice mixing ratio
       micptr, &    ! [INOUT]  indicator of microphysical processes
       qrp,    &    ! [INOUT]  rain PPVs
       qip,    &    ! [INOUT]  ice PPVs
       qap,    &    ! [INOUT]  aerosol PPVs
       pi,     &    ! [IN]     exner function
       ptot,   &    ! [IN]     pressure              (from SCALE)
       ivis,isect,from)
#ifndef DISABLE_AMPS
    use par_amps
    !use common_physics, only: qsaturation, qisaturation
    use mod_amps_lib, only: get_sat_vapor_pres_lk
    use scale_prc, only: PRC_abort
#endif
    implicit none

    ! arguments
    integer,  intent(in) :: KS, KE
    integer,  intent(in) :: npr,nbr,ncr,npi,nbi,nci,npa,nba,nca
    integer,  intent(in) :: level
    integer,  intent(in) :: nbhzcl
    integer,  intent(in) :: i, j

    integer,  intent(in) :: ivis,isect

    integer,  intent(inout), dimension(KS-1:KE) :: micptr

    real(RP), intent(inout), dimension(KS-1:KE) :: t, thp, th, qtp, qv, qc, qr, qi!, t0
    real(RP), intent(in),    dimension(KS-1:KE) :: den, pi, ptot!, den0
    real(RP), intent(inout), dimension(npr,nbr,ncr,KS-1:KE) :: qrp
    real(RP), intent(inout), dimension(npi,nbi,nci,KS-1:KE) :: qip
    real(RP), intent(inout), dimension(npa,nba,nca,KS-1:KE) :: qap

    character(len=*), intent(in) :: from

#ifndef DISABLE_AMPS
! newspace
    integer :: k,iter,  &
              ipr,ibr,icr,ipi,ibi,ici,ipa,iba,ica,icit,niter1, &
              ibr_st, nwv
    real(RP) :: fct, t_tmp, qtotal
    real(RP) :: theta,til,th1,th2
    real(RP) :: qvc
    real(RP) :: qs, qvis, es, sw, cnr
    real(RP) :: ov,error
    real(RP) :: tol = 0.001
    real(RP) :: cpr
    real(RP) :: aklv
    real(RP) :: akiv
    real(RP) :: svw_lmt=0.0e+0

!
! NOTE: these limit parameters are very important. This makes the model
!       slower or faster, and also affect predicted fields!
!
    ! mass mixing ratio limit for cloud droplets
    !  rad=0.5 micron and 1 1/m^3
    real(DP), parameter :: RCLMT=0.5d-15
    ! mass mixing ratio limit for rain
    !  rad=100 micron and 0.025 1/m^3
    real(DP), parameter :: RRLMT=1.0d-10
    ! mass mxing ratio limit for ice
    !  rad=0.5 micron and 1 1/m^3
    real(DP), parameter :: RILMT=1.0d-15
    ! for bin
    real(DP), parameter :: RRLMTB=1.0d-22
    real(DP), parameter :: RILMTB=1.0d-22
    !real(DP), parameter :: RILMTB=1.0e-18


    cpr  = CPdry / Rdry
    aklv = LHV0 / CPdry
    akiv = LHS0 / CPdry


    if ( level.ge.5 ) then
!     qr starts from nbhzcl+1 bin of qrp for bin microphysics (SHIPS)
!     liquid bin spectrum that uses nbhzcl bins as haze and cloud droplet.

      ibr_st=nbhzcl+1
      nwv=2
      do k = KS-1, KE

        qr(k)=0.0_RP
        qi(k)=0.0_RP
        qc(k)=0.0_RP
        micptr(k)=0
        cnr=0.0_RP
        qtotal = 0.0_RP
        do icr=1,ncr
          do ibr=1,ibr_st-1
            qtotal = qtotal + qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k)
            if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
              qc(k)=qc(k)   &
                    +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
              micptr(k)=1
              cnr=cnr+qrp(rcon_q,ibr,icr,k)
            else
              do ipr=1,npr
                qrp(ipr,ibr,icr,k)=0.0_RP
              end do
            end if
          end do

          do ibr=ibr_st,nbr
            qtotal = qtotal + qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k)
            if(qrp(rmt_q,ibr,icr,k).ge.RRLMTB) then
              qr(k)=qr(k) &
                   +max(0.0_RP,qrp(rmt_q,ibr,icr,k)-qrp(rmat_q,ibr,icr,k))
              micptr(k)=1
              cnr=cnr+qrp(rcon_q,ibr,icr,k)
            else
              do ipr=1,npr
                qrp(ipr,ibr,icr,k)=0.0_RP
              end do
            endif
          end do
        enddo

        do ibi=1,nbi
          do ici=1,nci
            qtotal = qtotal + qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)
            if(qip(imt_q,ibi,ici,k).ge.RILMTB) then
              qi(k)=qi(k) &
                   +max(0.0_RP,qip(imt_q,ibi,ici,k)-qip(imat_q,ibi,ici,k)-qip(imw_q,ibi,ici,k))

              qr(k)=qr(k) &
                   +max(0.0_RP,qip(imw_q,ibi,ici,k))
              micptr(k)=1
            else
              do ipi=1,npi
                qip(ipi,ibi,ici,k)=0.0_RP
              end do
            endif
          enddo
        enddo

        ! calculate the total water
        qtotal = qtotal + qv(k)
        qtp(k) = qtotal

        ! calculate the theta il
        thp(k)=th(k)/(1.0_RP  +(aklv*(qr(k)+qc(k))   &
                     +akiv*qi(k))/max(t(k),253.0_RP))

!------ This is initial guess ------
        th1=th(k)
        ov=qv(k)
        !deb_t = t(k)
        if(ivis<=1) then
          qv(k)=max(0.0_RP,qtp(k)-qr(k)-qi(k)-qc(k))
          !qv(k)=max(0.D0,qtotal-qr(k)-qi(k)-qc(k))

          til=thp(k)*pi(k)/CPdry
          t(k)=til*(1.0_RP+  &
               (aklv*(qr(k)+qc(k))+akiv*qi(k))/253.0_RP)
          if(t(k)>253.0_RP) then
            t(k)=0.5_RP*(til+sqrt(til**2+4.0_RP*til*  &
                 (aklv*(qr(k)+qc(k))+akiv*qi(k))))
            if(t(k)<253.0_RP-1.0e-3_RP) then
               LOG_ERROR("moistthermo2_scale",*) "mo2 t(k) is not right",t(k)
              call PRC_abort
            endif
          endif
          th1=thp(k)*(1.0_RP+(aklv*(qr(k)+qc(k)) &
              +akiv*qi(k))/max(t(k),253.0_RP))

          th(k)=th1
          theta=th1*(EPSvap+qv(k))/EPSvap/(1.0_RP+qv(k))

          es = get_sat_vapor_pres_lk(1, t(k), estbar, esitbar )/10.0_RP
          qs = es/(ptot(k) - es)*0.622_RP
          sw=ptot(k)*qv(k) / ( 0.622_RP * ( 1.0_RP+0.61_RP*qv(k) ) ) / es - 1.0_RP
          if(sw.gt.0.50_RP) then
             LOG_INFO("moistthermo2_scale",'("mo2 sup:",a,4I5,20ES15.6)')  &
                  from,isect,k,i,j, &
                  ptot(k),qv(k),qr(k)+qc(k),qi(k),t(k), &
                  thp(k),pi(k), &
                  sw
          endif

          if(ivis==1) then
            if(qv(k).gt.qs) micptr(k)=1
          end if
        end if

        if(ivis>=1.and.t(k).lt.273.16_RP) then
          es = get_sat_vapor_pres_lk(2, t(k), estbar, esitbar )/10.0_RP
          qvis = es/(ptot(k) - es)*0.622_RP
          if(qv(k).gt.qvis) micptr(k)=1
        end if

      enddo

    end if

#endif
    return
  end subroutine moistthermo2_scale

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_MP_amps_finalize

#ifndef DISABLE_AMPS
    deallocate( ADVPPMZ, ADVPPMZE )

    deallocate( mmassrv_global, mmassiv_global )
    deallocate( qapv_ini )
    deallocate( den_ini )

    deallocate( ATMOS_PHY_MP_amps_tracer_names )
    deallocate( ATMOS_PHY_MP_amps_tracer_descriptions )
    deallocate( ATMOS_PHY_MP_amps_tracer_units )
    deallocate( ATMOS_PHY_MP_amps_tracer_ref )
#endif

    return
  end subroutine ATMOS_PHY_MP_amps_finalize

  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  subroutine ATMOS_PHY_MP_amps_tendency( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dt, U, V, &
       W, MOMZ, &
       DENS, PRES, TEMP, &
       POTT, EXNER,      &
       QTRC, QDRY,       &
       CPtot, CVtot,     &
       CCN,              &
       CZ, FZ,           &
       RHOQ_t,           & ! hydrometeor tendency output
       RHOE_t,           & ! diabatic heating tendency output
       CPtot_t,          & ! constant pressure specific tendency output
       CVtot_t,          & ! constant volume specific tendency output
       DENS_t,           & ! density tendency output, this includes dry air, vapor, liquid, and ice
       MOMZ_t,           & ! vertical momentum tendency output
       RHOU_t,           & ! u momentum tendency output
       RHOV_t,           & ! v momentum tendency output
       SFLX_rain,        & ! surface rain flux output
       SFLX_snow         ) ! surface snow flux output
#ifndef DISABLE_AMPS
    use scale_prc, only: &
       PRC_IsMaster, &
       PRC_myrank, &
       PRC_abort
    use scale_const, only: &
       LHV0 => CONST_LHV0, &
       LHF0 => CONST_LHF0
    use scale_atmos_hydrometeor, only: &
       CP_VAPOR, &
       CP_WATER, &
       CP_ICE, &
       CV_VAPOR, &
       CV_WATER, &
       CV_ICE
    use scale_file_history, only: &
       FILE_HISTORY_in

    use com_amps
    use par_amps
    use mod_amps_lib, only:  &
       ini_aerosols_profile_mpace, &
       ini_aerosols_profile_const, &
       ini_aerosols_profile_sheba, &
       get_sat_vapor_pres_lk

    use class_Cloud_Micro, only: &
       make_cloud_micro_cnfg, &
       set_cloud_micro_cur

    use mod_amps_lib, only: ifc_cloud_micro

!    use mod_scladv_slppm
    use mod_amps_utility, only: &
       assign_seeds,moistthermo2,findcond1,findcond2,diabatic,negadj,negadj_bin2, &
       negadj_bin_ap,cal_mmass,cal_mmass_scale,shift_bin_vec, &
       fill_bulk_ap,fill_bin_ap,integ_mictend1,integ_mictend2, &
       srand2, &
       sclsedprz_original,sclsedaer
#endif
    implicit none

!tst    use module_mp_clr_two_moment, only: clr_two_moment_nms

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(DP), intent(in) :: dt
    real(RP), intent(in) :: U    (KA,IA,JA) ! cell center velocity
    real(RP), intent(in) :: V    (KA,IA,JA) ! cell center velocity
    real(RP), intent(in) :: W    (KA,IA,JA) ! cell center velocity
    real(RP), intent(in) :: MOMZ (KA,IA,JA) ! cell center velocity
    real(RP), intent(in) :: DENS (KA,IA,JA)
    real(RP), intent(in) :: PRES (KA,IA,JA)
    real(RP), intent(in) :: TEMP (KA,IA,JA)
    real(RP), intent(in) :: POTT (KA,IA,JA)
    real(RP), intent(in) :: EXNER(KA,IA,JA)
    real(RP), intent(in) :: QTRC (KA,IA,JA,QA)
    real(RP), intent(in) :: QDRY (KA,IA,JA)
    real(RP), intent(in) :: CPtot(KA,IA,JA)
    real(RP), intent(in) :: CVtot(KA,IA,JA)
    real(RP), intent(in) :: CCN  (KA,IA,JA)
    real(RP), intent(in) :: CZ   (KA,IA,JA)
    real(RP), intent(in) :: FZ   (0:KA,IA,JA)

    real(RP), intent(out) :: RHOQ_t (KA,IA,JA,QA) ! tendency of tracers
    real(RP), intent(out) :: RHOE_t (KA,IA,JA)    ! diabatic heating
    real(RP), intent(out) :: CPtot_t(KA,IA,JA)    ! specific heat tendency
    real(RP), intent(out) :: CVtot_t(KA,IA,JA)    ! specific heat tendency
    real(RP), intent(out) :: DENS_t(KA,IA,JA)     ! density tendency
    real(RP), intent(out) :: MOMZ_t(KA,IA,JA)     ! z momentum tendency
    real(RP), intent(out) :: RHOU_t(KA,IA,JA)     ! rho u tendency
    real(RP), intent(out) :: RHOV_t(KA,IA,JA)     ! rho v tendency

    real(RP), intent(out) :: SFLX_rain(IA,JA)     ! surface rain flux
    real(RP), intent(out) :: SFLX_snow(IA,JA)     ! surface snow flux

#ifndef DISABLE_AMPS
    real(RP),dimension(KS-1:KE) :: thv,thetav,moist_denv,tv,piv,wbv,momv,ptotv,pbv, &
                                   qvv,qcv,qrv,qiv,qtpv, &
                                   waccv
    real(RP),dimension(KS-1:KE,2) :: trpv, ftrpv

    real(RP) :: spdsfcv,thskinv

    real(RP) :: factor_mxr1, factor_mxr2, qtotal(KS-1:KE), qtotal2(KS-1:KE), qtotal3(KS-1:KE)
    real(RP) :: factor_e1, factor_e2, factor_e1_, factor_e2_

    integer,dimension(KS-1:KE) :: micptrv,imicv,jmicv,kmicv
    integer :: kmic(1:KE+1)

    real(RP) :: dzzmv(KS-1:KE), dzvmv(KS-1:KE)

    real(RP),dimension(npr,nbr,ncr,KS-1:KE) :: qrpv
    real(RP),dimension(npi,nbi,nci,KS-1:KE) :: qipv
    real(RP),dimension(npa,nba,nca,KS-1:KE) :: qapv
    real(RP),dimension(nbr,ncr,KS-1:KE) :: qrov
    real(RP),dimension(nbi,nci,KS-1:KE) :: qiov

    real(RP), dimension(KS-1:KE) :: den_t                         ! for scale sedimentation

    real(RP), dimension(KS-1:KE) :: den_diff                      ! density error in AMPS
    real(RP), dimension(KS-1:KE) :: DENS_NEW                      ! density change in SCALE

    integer,dimension(nbr,ncr) :: k1br,k2br
    integer,dimension(nbi,nci) :: k1bi,k2bi
    integer,dimension(nba,nca) :: k1ba,k2ba

    real(RP),dimension(nbr,ncr,KS-1:KE) :: mmassrv
    real(RP),dimension(nbi,nci,KS-1:KE) :: mmassiv

    real(RP),dimension(nzh) :: thetavm,moist_denvm,tvm,wbvm,ptotvm,pivm, &
                               pbvm,qvvm,qcvm,qrvm,qivm,v3v,zstv,dzzmvm
    real(RP),dimension(nzh,2) :: trpvm
    integer, dimension(nzh) :: micptrvm,imicvm,jmicvm,kmicvm

    real(RP),dimension(npr,nbr,ncr,nzh) :: qrpvm
    real(RP),dimension(npi,nbi,nci,nzh) :: qipvm
    real(RP),dimension(npa,nba,nca,nzh) :: qapvm

    real(RP),dimension(10,2,nzh)        :: dmtendlm, dcontendlm
    real(RP),dimension(10,2,KS-1:KE)    :: dmtendl, dcontendl
    real(RP),dimension(3,2,mxnbin,nzh)     :: dbintendlm
    real(RP),dimension(3,2,mxnbin,KS-1:KE) :: dbintendl

    real(RP) :: dz1v
    real(RP) :: pgnd

    integer,dimension(1) :: k1etatv

    integer :: i, j, k, m
    integer :: k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i,k1a,k2a
    integer :: ipr_qpr,ipr,ibr,icr,ipi_qi,ipi,ibi,ici,icit,ipa_qpa,ipa,iba,ica
    integer :: idir, isect, nsect, istrt
    integer :: isprayr,isprayi,isn,iph

    integer  :: nmic

    !character(max_char_len) :: name, units, name_dmtend(7)
    !integer,dimension(7) :: iswt_dmtend=0

    !$ integer :: omp_get_thread_num,omp_in_parallel,omp_get_num_threads,omp_get_max_threads
    integer :: var_status

!!!    integer,dimension(max_cpu) :: seed_sec

    integer,allocatable :: seed_sec(:)
    integer,allocatable,dimension(:) :: JSEED,IFRST,isect_seed,NEXTN

    integer,parameter :: iproc_t=0

    real(RP) :: Emoist(KS:KE,2)
    real(RP) :: dq

    ! history
    real(RP) :: AMPS_mt(KA,IA,JA,20) ! mass tendency
    real(RP) :: AMPS_ct(KA,IA,JA,20) ! concentration tendency
    real(RP) :: AMPS_bt(KA,IA,JA,nbr+nbi,3) ! mass tendency
    real(RP) :: AMPS_tv(KA,IA,JA,nbr+nbi,2) ! terminal velocity
    character(7) :: AMPS_tv_NAME, AMPS_bin_NAME
    real(RP) :: QLIQ(KA,IA,JA)
    real(RP) :: QICE(KA,IA,JA,6) ! total, rime, aggregate, crystal, melt water, frozen water

    LOG_PROGRESS(*) 'atmosphere / physics / microphysics / AMPS'

    ! initialize output array
    !$omp workshare
    RHOQ_t(:,:,:,:) = 0.0_RP
    DENS_t(:,:,:) = 0.0_RP
    MOMZ_t(:,:,:) = 0.0_RP
    RHOU_t(:,:,:) = 0.0_RP
    RHOV_t(:,:,:) = 0.0_RP
    CPtot_t(:,:,:) = 0.0_RP
    CVtot_t(:,:,:) = 0.0_RP
    RHOE_t(:,:,:) = 0.0_RP
    SFLX_rain(:,:) = 0.0_RP
    SFLX_snow(:,:) = 0.0_RP

    ! history
    AMPS_mt(:,:,:,:) = 0.0_RP
    AMPS_ct(:,:,:,:) = 0.0_RP
    AMPS_bt(:,:,:,:,:) = 0.0_RP
    AMPS_tv(:,:,:,:,:) = 0.0_RP
    !$omp end workshare

    if (amps_ignore) return

    ! initialize
    call PROF_rapstart("amps_make",3)
    call PROF_rapend("amps_make",3)
    call PROF_rapstart("amps_micro",3)
    call PROF_rapend("amps_micro",3)
    call PROF_rapstart("amps_tendency",3)
    call PROF_rapend("amps_tendency",3)
    call PROF_rapstart("amps_update_diagnose",3)
    call PROF_rapend("amps_update_diagnose",3)
    call PROF_rapstart("amps_repair",3)
    call PROF_rapend("amps_repair",3)
    call PROF_rapstart("amps_collisioncoal",3)
    call PROF_rapend("amps_collisioncoal",3)
    call PROF_rapstart("amps_aggregation",3)
    call PROF_rapend("amps_aggregation",3)
    call PROF_rapstart("amps_riming",3)
    call PROF_rapend("amps_riming",3)
    call PROF_rapstart("amps_meltshed",3)
    call PROF_rapend("amps_meltshed",3)
    call PROF_rapstart("amps_hydrobreakup",3)
    call PROF_rapend("amps_hydrobreakup",3)
    call PROF_rapstart("amps_icenucleation",3)
    call PROF_rapend("amps_icenucleation",3)
    call PROF_rapstart("amps_ccn",3)
    call PROF_rapend("amps_ccn",3)
    call PROF_rapstart("amps_vapordep",3)
    call PROF_rapend("amps_vapordep",3)
    call PROF_rapstart("amps_meltshed",3)
    call PROF_rapend("amps_meltshed",3)
    call PROF_rapstart("amps_sedimentation",3)
    call PROF_rapend("amps_sedimentation",3)


    idir = 3
    waccv = -1.0e+0 ! negative one, so that the velocity is negative downwards

    isprayr = 0
    isprayi = 0

    ! redundant, to be removed
    micro_io_strt = .true.

    nsect = 1
    !$ nsect = omp_get_max_threads()

    if (.not. allocated(CM)) then
       allocate(CM(nsect))
    endif

    if(level.eq.5) then
       ! assign seeds for each section
       if(.not.allocated(seed_sec)) then
          allocate(seed_sec(nsect),jseed(nsect),ifrst(nsect) &
                  ,isect_seed(nsect),nextn(nsect),stat=var_status)
          if(var_status.ne.0) then
             LOG_ERROR("ATMOS_PHY_MP_amps_tendency",*) "seed_sec not allocated",var_status, nsect
             call PRC_abort
          endif
          nextn(:) = 0
       endif
       call assign_seeds(nsect,seed_sec)
    endif

    !call PROF_rapstart("amps_micro",3)


    if(level==5) then
!   copy seed variables for the given section
       do i = 1, nsect
          call srand2(jseed(i),ifrst(i) &
                     ,isect_seed(i) &
                     ,seed_sec(i),i)
       enddo

!#######################################################################################
!-------- construct AMPS variables --------
! A cloud object CM is constructed. Physical coefficients are transferred into global
! variables in CM. Four objects under CM, i.e. aerosols (group), liquid (group), ice
! (group), and air (air) are each constructed. The shape are mass bin are stored in the
! group object. Token is used to differentiate which types of hydrometeors in group.
!#######################################################################################
       call PROF_rapstart("amps_make",3)
       if(TIME_AMPS == 0) then
          do i = 1, nsect
             call make_Cloud_Micro_cnfg(CM(i) &
                      ,dt, n_step_cl,n_step_vp, TIME_AMPS*dt  &
                      ,npr,nbr,ncr  &
                      ,npi,nbi,nci  &
                      ,npa,nba,nca  &
                      ,estbar,esitbar  &
                      ,level_comp,debug_level,coll_level,out_type,T_print_period,output_format  &
                      ,token_r,token_s,token_a,dtype_c,dtype_r,dtype_s  &
                      ,dtype_a  &
                      ,hbreak_c,hbreak_r,hbreak_s,flagp_c,flagp_r,flagp_s  &
                      ,flagp_a,act_type  &
                      ,micexfg  &
                      ,fcon_c,coef_ap,eps_ap  &
                      ,srat_c,srat_r,srat_s,srat_a  &
                      ,sadd_c,sadd_r,sadd_s,sadd_a  &
                      ,minmass_r,minmass_s,minmass_a  &
                      ,den_apt,den_aps,den_api  &
                      ,binbr,binbi  &
                      ,imin_bk,imax_bk,jmin_bk,jmax_bk,bu_fd,bu_tmass  &
                      ,APSNAME,nu_aps,phi_aps,m_aps,ap_lnsig,ap_mean  &
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
                      ,nzh)
          enddo
       else
          do i = 1, nsect
             call set_Cloud_Micro_cur(CM(i),TIME_AMPS*dt)
          enddo
       endif

    endif
    call PROF_rapend("amps_make",3)

    !$omp workshare
    QLIQ(:,:,:) = 0.0_RP
    QICE(:,:,:,:) = 0.0_RP
    !$omp end workshare


    !$omp parallel do OMP_SCHEDULE_ collapse(2) default(none) &
    !$omp private(i,j,k,m, &
    !$omp         qrpv,qipv,qapv,qcv,qrv,qiv,qrpvm,qipvm,qapvm,qcvm,qrvm,qivm,qrov,qiov, &
    !$omp         den_t,den_diff,DENS_NEW,moist_denv, &
    !$omp         moist_denvm, &
    !$omp         dz1v,k1etatv, &
    !$omp         k1r,k2r,k1br,k2br,k1i,k2i,k1bi,k2bi,k1m,k2m,k1c,k2c, &
    !$omp         dmtendl,dmtendlm,dcontendl,dcontendlm,dbintendl,dbintendlm, &
    !$omp         factor_mxr1,factor_mxr2, &
    !$omp         ptotv,tv,thv,piv,pbv,qvv,thetav,wbv,momv,trpv, &
    !$omp         ptotvm,tvm,v3v,pivm,pbvm,qvvm,thetavm,wbvm,trpvm,zstv,dzzmvm, &
    !$omp         qtotal,qtotal2,qtotal3, &
    !$omp         mmassrv,mmassiv, &
    !$omp         imicv,kmicv,jmicv,imicvm,kmicvm,jmicvm,ftrpv, &
    !$omp         icr,ipr,ipr_qpr,ibr,ici,ipi,ipi_qi,ibi,ica,ipa,ipa_qpa,iba, &
    !$omp         iph,ivis,istrt,isn, &
    !$omp         micptrv,micptrvm,nmic,kmic, &
    !$omp         isect, &
    !$omp         pgnd,thskinv,spdsfcv, &
    !$omp         dzzmv,dzvmv, &
    !$omp         Emoist,dq) &
    !$omp shared(CM, &
    !$omp        nz,nzh, &
    !$omp        IS,JS,KS,IE,JE,KE,IA,JA,KA, &
    !$omp        level,l_gaxis_version,l_bin_shift,l_axis_limit,l_fix_aerosols,l_sediment,l_fill_aerosols,ini_aerosol_prf,amps_debug, &
    !$omp        jseed,isect_seed,nextn,ifrst,seed_sec, &
    !$omp        TIME_AMPS,dt, &
    !$omp        QDRY,QTRC,DENS,W,MOMZ,PRES,TEMP,U,V,CVtot,SFLX_rain,SFLX_snow, &
    !$omp        GRAV,PRE00,Rdry,CPdry,CP_VAPOR,CP_WATER,CV_VAPOR,CV_WATER,CP_ICE,CV_ICE,EPS,LHF0, &
    !$omp        AMPS_mt,AMPS_bt,AMPS_ct,AMPS_tv,QLIQ,QICE, &
    !$omp        CPtot_t,CVtot_t,MOMZ_t,RHOU_t,RHOV_t,RHOE_t,RHOQ_t,DENS_t, &
    !$omp        CZ,FZ, &
    !$omp        estbar,esitbar, &
    !$omp        rmt_q,rmat_q,rmas_q,rcon_q,imt_q,imr_q,ima_q,imc_q,imw_q,imf_q,imat_q,imas_q, &
    !$omp        icon_q,iacr_q,iccr_q,idcr_q,iag_q,icg_q,ivcs_q,inex_q, &
    !$omp        amt_q,acon_q,ams_q, &
    !$omp        ncr,nbr,npr,nci,nbi,npi,nca,nba,npa, &
    !$omp        binbr,binbi,nbin_h,nbhzcl, &
    !$omp        I_QV,I_QL,I_QPPVL,I_QI,I_QW,I_QPPVI,I_QPPVA, &
    !$omp        qapv_ini,den_ini, &
    !$omp        N_ap_ini,coef_ap,eps_ap,den_apt, &
    !$omp        advppmz,advppmze, &
    !$omp        mmassrv_global,mmassiv_global, &
    !$omp        level_comp,micro_io_strt,idir,isprayr,isprayi,iadvv,waccv, &
    !$omp        nsect)
    do j = JS, JE
    do i = IS, IE

       isect = 1
       !$ isect = omp_get_thread_num()+1

       qrpv = 0.0_RP
       qipv = 0.0_RP
       qapv = 0.0_RP
       qrpvm = 0.0_RP
       qipvm = 0.0_RP
       qapvm = 0.0_RP

       den_t = 0.0_RP
       den_diff = 0.0_RP
       DENS_NEW = 0.0_RP
       moist_denv = 0.0_RP

       dz1v = FZ(KS,i,j) - FZ(KS-1,i,j)

       k1etatv(:) = k1eta

       k1r = 10000
       k2r = -10000
       k1i = 10000
       k2i = -10000

       dmtendl = 0.0_RP
       dmtendlm = 0.0_RP
       dcontendl = 0.0_RP
       dcontendlm = 0.0_RP
       dbintendl = 0.0_RP
       dbintendlm = 0.0_RP

       ! total water content
       qtotal = 0.0_RP

!mark1
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     1.1. Initialize over a vertical column
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
       Z_LOOP_01: do k = KS, KE
          factor_mxr1 = (QDRY(k,i,j) + &
                         QTRC(k,i,j,I_QV)) ! moist air mixing ratio
          factor_mxr2 = (QDRY(k,i,j)*DENS(k,i,j) + &
                         QTRC(k,i,j,I_QV)*DENS(k,i,j))
          ! rho[kg m-3] = 0.001 rho[g cm-3]

          !
          ! set the thermodynamics variables to the updated state.
          !
          !k=k_kid+1

          ! pressure
          ptotv(k) = PRES(k,i,j)

          ! temperature
          tv(k) = TEMP(k,i,j)

          ! potential temperature
          thv(k) = tv(k)*(PRE00/ptotv(k))**(Rdry/CPdry)

          ! exner function
          piv(k) = tv(k)/thv(k)*CPdry

          ! unknown
          pbv(k) = 0.0_RP

          ! density
          moist_denv(k) = DENS(k,i,j)*factor_mxr1

          ! vapor mixing ratio
          qvv(k) = QTRC(k,i,j,I_QV)/factor_mxr1
          qtotal(k) = qtotal(k) + qvv(k)*factor_mxr2
          Emoist(k,1) = - LHV0 * QTRC(k,i,j,I_QV) * DENS(k,i,j)

          ! virtual potential temperature
          thetav(k) = thv(k)*(1.0e+0+0.61*qvv(k))

          ! w at full grid
          wbv(k) = W(k,i,j)
          ! w at half grid
          momv(k) = MOMZ(k,i,j)

          kmicv(k) = k-KS+1
          imicv(k) = i
          jmicv(k) = j

          ftrpv(k,:) = 0.0_RP

! mark2
          !
          ! Initialize liquid spectrum
          !
          do icr = 1, ncr ! assume to be 1
             ipr_qpr = 0
             do ibr = 1, nbr
                qrpv(rmt_q,ibr,icr,k) = ( QTRC(k,i,j,I_QL+ibr-1) + QTRC(k,i,j,I_QPPVL+ipr_qpr) ) / factor_mxr1
                do ipr = rmat_q, rmas_q
                   qrpv(ipr,ibr,icr,k) = QTRC(k,i,j,I_QPPVL+ipr_qpr)/factor_mxr1
                   ipr_qpr = ipr_qpr + 1
                enddo
                ! number concentration, non-mass variable
                qrpv(rcon_q,ibr,icr,k) = QTRC(k,i,j,I_QPPVL+ipr_qpr)/factor_mxr1/0.001_RP
                ipr_qpr = ipr_qpr + 1
                do ipr = npr-1, npr
                   qrpv(ipr,ibr,icr,k) = 0.0
                enddo
                mmassrv(ibr,icr,k) = 0.0
                if(qrpv(rmt_q,ibr,icr,k)>1.0e-25) then
                   k1r = min(k1r,k-KS+1)
                   k2r = max(k2r,k-KS+1)
                endif
                qrov(ibr,icr,k) = qrpv(rmt_q,ibr,icr,k)
                ! total water content
                qtotal(k) = qtotal(k) + (qrpv(rmt_q,ibr,icr,k) - qrpv(rmat_q,ibr,icr,k))*factor_mxr2
             enddo

          enddo

          !
          ! Initialize ice spectrum
          !
          do ici = 1, nci ! assumed to be 1
             ipi_qi = 0
             do ibi = 1, nbi
                ! total ice mass
                qipv(imt_q,ibi,ici,k) = ( QTRC(k,i,j,I_QI+ibi-1) &
                                        + QTRC(k,i,j,I_QW+ibi-1) &
                                        + QTRC(k,i,j,I_QPPVI+ipi_qi+4) ) / factor_mxr1
                ! rimed, aggregate, and crystal
                do ipi = imr_q, imc_q
                   qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                   ipi_qi = ipi_qi + 1
                enddo
                ! melt water
                qipv(imw_q,ibi,ici,k) = QTRC(k,i,j,I_QW+ibi-1)/factor_mxr1
                ! frozen water
                qipv(imf_q,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                ipi_qi = ipi_qi + 1

                ! aerosol mass
                do ipi = imat_q, imas_q
                   qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1
                   ipi_qi = ipi_qi + 1
                enddo

                ! number concentration and other non-mass variables
                if (l_gaxis_version == 1) then
                   ! -- ORIGINAL

                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo

                else if (l_gaxis_version == 2) then
                   ! -- METHOD 1
                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo

                   if (qipv(icon_q,ibi,ici,k) < 1.e-22_RP) then

                      qipv(iccr_q,ibi,ici,k) = 0.0_RP
                      qipv(idcr_q,ibi,ici,k) = 0.0_RP
                      qipv(icg_q,ibi,ici,k) = 0.0_RP

                   else
                      qipv(iccr_q,ibi,ici,k) = qipv(iccr_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(idcr_q,ibi,ici,k) = qipv(idcr_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(icg_q,ibi,ici,k) = qipv(icg_q,ibi,ici,k)*qipv(iag_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                   endif

                else if (l_gaxis_version == 3) then
                   ! -- METHOD 2
                   do ipi = imas_q+1, npi-2
                      qipv(ipi,ibi,ici,k) = QTRC(k,i,j,I_QPPVI+ipi_qi)/factor_mxr1/0.001_RP
                      ipi_qi = ipi_qi + 1
                   enddo
                   if (qipv(icon_q,ibi,ici,k) < 1.e-22_RP) then
                      qipv(iag_q,ibi,ici,k) = 0.0_RP
                      qipv(icg_q,ibi,ici,k) = 0.0_RP
                   else
                      qipv(iag_q,ibi,ici,k) = qipv(iag_q,ibi,ici,k)*qipv(iacr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                      qipv(icg_q,ibi,ici,k) = qipv(icg_q,ibi,ici,k)*qipv(iccr_q,ibi,ici,k)/qipv(icon_q,ibi,ici,k)
                   endif

                endif

             enddo

             do ibi = 1, nbi
                do ipi = npi-1, npi
                   qipv(ipi,ibi,ici,k) = 0.0_RP
                enddo
                mmassiv(ibi,ici,k)=0.0_RP
                if(qipv(imt_q,ibi,ici,k)>1.0e-25) then
                   k1i=min(k1i,k)
                   k2i=max(k2i,k)
                endif
                qiov(ibi,ici,k)=qipv(imt_q,ibi,ici,k)
                ! total water content
                qtotal(k) = qtotal(k) + (qipv(imt_q,ibi,ici,k) - qipv(imat_q,ibi,ici,k))*factor_mxr2
             enddo

          enddo

          !
          ! Initialize aerosol spectrum
          !
          ipa_qpa = 0
          do ica = 1, nca
             do iba = 1, nba
                !
                !do ipa=1,npa-2
                qapv(amt_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa)/factor_mxr1
                qapv(acon_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa+1)/factor_mxr1/0.001_RP
                qapv(ams_q,iba,ica,k) = QTRC(k,i,j,I_QPPVA+ipa_qpa+2)/factor_mxr1
                ipa_qpa = ipa_qpa + 3
                !enddo
                do ipa = npa-1, npa
                   qapv(ipa,iba,ica,k)=0.0_RP
                enddo
             enddo
          enddo


          do ibi = 1, nbi
             Emoist(k,1) = Emoist(k,1) + LHF0 * QTRC(k,i,j,I_QI+ibi-1) * DENS(k,i,j)
          end do

       enddo Z_LOOP_01
       ! set underground, this is used for surface flux
       tv(KS-1) = tv(KS)
       pbv(KS-1) = 0.0_RP
       qvv(KS-1) = qvv(KS)
       ptotv(KS-1) = ptotv(KS)*exp( GRAV/Rdry/(tv(KS-1)*(1.0_RP + 0.61_RP*qvv(KS-1))) * ( CZ(KS,i,j) - CZ(KS-1,i,j) ) )
       moist_denv(KS-1) = ptotv(KS-1)/(Rdry*tv(KS-1)*(1.0_RP+0.61_RP*qvv(KS-1)))
       piv(KS-1) = (ptotv(KS-1)/PRE00)**(Rdry/CPdry)*CPdry
       thv(KS-1) = tv(KS-1)*CPdry/piv(KS-1)
       thetav(KS-1) = thv(KS-1)*(1.0_RP+0.61_RP*qvv(KS-1)/QDRY(KS,i,j))

       pgnd=0.5_RP*(ptotv(KS-1)+ptotv(KS))

       thskinv = thv(KS-1)
       spdsfcv = 0.0_RP
       kmicv(KS-1) = 1
       imicv(KS-1) = i
       jmicv(KS-1) = j


       ! grid vertical spacing
       do k = KS, KE
          dzzmv(k) = 1.0_RP / ( FZ(k,i,j) - FZ(k-1,i,j) )
       enddo
       dzzmv(KS-1) = dzzmv(KS)

       ! half grid vertical spacing (for MOMZ)
       do k = KS-1, KE-1
          dzvmv(k) = 1.0_RP / ( CZ(k+1,i,j) - CZ(k,i,j) )
       enddo

!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     1.2   retrieve mean mass calculated at the end of previous call of microphysics scheme
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
       if (level.ge.5) then
          do k = KS-1, KE
             do icr = 1, ncr
             do ibr = 1, nbr
                mmassrv(ibr,icr,k) = 0.0_RP
             enddo
             enddo
             do ici = 1, nci
             do ibi = 1, nbi
                mmassiv(ibi,ici,k) = 0.0_RP
             enddo
             enddo
          enddo
          if (TIME_AMPS == 0) then
             if (k2r.gt.0) call cal_mmass(mmassrv &
                                         ,k1r,k2r,npr,nbr,ncr,qrpv,nzh &
                                         ,0)
             if (k2i.gt.0) call cal_mmass(mmassiv &
                                         ,k1i,k2i,npi,nbi,nci,qipv,nzh &
                                         ,1)
          else
             if (k2r.gt.0) then
                do k = KS, KE
                   mmassrv(:,:,k) = mmassrv_global(:,:,k,i,j)
                enddo
             endif
             if (k2i.gt.0) then
                do k = KS, KE
                   mmassiv(:,:,k) = mmassiv_global(:,:,k,i,j)
                enddo
             endif
          endif
       endif

!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     1.3   shift bins in rain and ice spectra for AMPS
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
       if(level.eq.5 .and. TIME_AMPS > 0 .and. l_bin_shift) then
          if (k2r.gt.0) then
             iph=1
             call shift_bin_vec(nzh,npr,nbr,ncr,k1r,k2r &
                               ,binbr,qrpv,mmassrv &
                               ,moist_denv,iph,kmicv,imicv,jmicv)
          endif
          if (k2i.gt.0) then
             iph=2
             call shift_bin_vec(nzh,npi,nbi,nci,k1i,k2i &
                               ,binbi,qipv,mmassiv &
                               ,moist_denv,iph,kmicv,imicv,jmicv)
          endif
       endif

!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     1.4   diagnose thermodynamical variables
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
!     diagnose the vapor content from the (total hydrometeors - liquid+ice)
!     diagnose the total hydrometeors, density, and ice-liquid potential temperature
       ivis = 1
       call moistthermo2_scale( KS, KE, &
                                npr,nbr,ncr,npi,nbi,nci,npa,nba,nca,   &
                                level,nbhzcl,   &
                                i,j, &
                                trpv(:,1),thv,moist_denv,tv,  &
                                trpv(:,2),qvv,qcv,qrv,qiv,micptrv,   &
                                qrpv,qipv,qapv,piv,ptotv,  &
                                ivis, isect,'first')

       if (TIME_AMPS == 0) then
          ! aerosol initial profiles (KS:KE,IS:IE,JS:JE)
          ! the input is kgm-3/kgm-3 and m-3/kgm-3
          ! the output's unit is gcm-3/gcm-3 and cm-3/gcm-3
! mark3
          select case(ini_aerosol_prf)
          case(1)
             call ini_aerosols_profile_sheba(qapv_ini(:,:,:,:,i,j),nzh,nca,CZ(KS-1:KE,i,j),moist_denv(:))
          case(2)
             call ini_aerosols_profile_mpace(qapv_ini(:,:,:,:,i,j),nzh,nca,CZ(KS-1:KE,i,j),moist_denv(:))
          case default
             call ini_aerosols_profile_const(qapv_ini(:,:,:,:,i,j),nzh,nca,CZ(KS-1:KE,i,j),moist_denv(:))
          end select

          den_ini(:,i,j) = moist_denv(:)
       endif

       ! refill aerosols to initial condition
       if (l_fix_aerosols) then
          do k = KS, KE
             do ica = 1, nca
                do iba = 1, nba
                   do ipa = 1, npa-2
                      qapv(ipa,iba,ica,k) = qapv_ini(ipa,iba,ica,k,i,j)*den_ini(k,i,j)/moist_denv(k)
                   enddo
                enddo
             enddo
          enddo
       endif

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.  Cloud Microphysics
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
       nmic=0
       do k = KS-1, KE
          if(micptrv(k)>0) then
             nmic=nmic+1
             kmic(nmic)=k
          endif
       enddo
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.1.  Gather cloudy grids
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------

       MIC_IF: if(nmic>0) then
          if(micro_io_strt(isect)) then
             istrt=1
          else
             istrt=0
          endif

          do m = 1, nmic
             k=kmic(m)
             kmicvm(m)=k-KS+2
             imicvm(m)=i
             jmicvm(m)=j

             trpvm(m,1)=trpv(k,1)
             trpvm(m,2)=trpv(k,2)
             thetavm(m)=thetav(k)
             moist_denvm(m)=moist_denv(k)
             !if(moist_denvm(m).lt.0.) then
             !  write(*,*) "part4 :density is negative:",m,moist_denvm(m)
             !end if

             qcvm(m)=qcv(k)
             v3v(m)=0.0_RP
             qvvm(m)=qvv(k)
             pivm(m)=piv(k)
             pbvm(m)=pbv(k)
             ptotvm(m)=ptotv(k)
             tvm(m)=tv(k)
             wbvm(m)=wbv(k)

             qrvm(m)=qrv(k)
             qivm(m)=qiv(k)

             zstv(m)=0.0_RP
             dzzmvm(m)=dzzmv(k)

             micptrvm(m)=micptrv(k)
             do icr = 1, ncr
                do ibr = 1, nbr
                   do ipr = 1, npr
                      qrpvm(ipr,ibr,icr,m)=qrpv(ipr,ibr,icr,k)
                   enddo
                enddo
             enddo
             do ici = 1, nci
                do ibi = 1, nbi
                   do ipi = 1, npi
                      qipvm(ipi,ibi,ici,m)=qipv(ipi,ibi,ici,k)
                   enddo
                enddo
             enddo
             do ica = 1, nca
                do iba = 1, nba
                   do ipa = 1, npa
                      qapvm(ipa,iba,ica,m)=qapv(ipa,iba,ica,k)
                   enddo
                enddo
             enddo

          enddo

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.2.  Call aerosol-cloud microphysics schemes
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
          if(level.ge.5) then
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.2.1.  AMPS
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------

             call PROF_rapstart("amps_micro",3)

             call ifc_cloud_micro( &
                                  CM(isect) &
                                  ,qcvm,v3v &
                                  ,qrpvm,npr,nbr,ncr  &
                                  ,qipvm,npi,nbi,nci &
                                  ,qapvm,npa,nba,nca &
                                  ,qvvm, moist_denvm, ptotvm, TVm, wbvm &
                                  ,nmic,imicvm,jmicvm,kmicvm,iproc_t,isect-1,istrt &
                                  ,trpvm(:,2),trpvm(:,1) &
                                  ,jseed(isect),ifrst(isect),isect_seed(isect) &
                                  ,nextn(isect)  &
                                  ,dmtendlm,dcontendlm,dbintendlm &
                                  )

             call PROF_rapend("amps_micro",3)

          end if

          micro_io_strt(isect)=.false.
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     2.3   diagnose thermodynamical variables
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
          ivis=0
          ! microphysics variables here are mixing ratio over moist air density
          call moistthermo2(thetavm,moist_denvm,tvm,ptotvm  &
                           ,trpvm(:,2),qvvm,qcvm,qrvm,qivm,micptrvm   &
                           ,npr,nbr,ncr,npi,nbi,nci   &
                           ,level,nbhzcl   &
                           ,trpvm(:,1),qrpvm,qipvm, pbvm  &
                           ,estbar,esitbar,pivm &
                           ,nmic,imicvm,jmicvm,kmicvm &
                           ,ivis,'af_cm', isect)

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.4.  Scatter back cloudy grids
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
          k1br = 10000
          k2br = -10000
          k1bi = 10000
          k2bi = -10000
          do m = 1, nmic
             k = kmic(m)

             trpv(k,1) = trpvm(m,1)
             trpv(k,2) = trpvm(m,2)
             moist_denv(k) = moist_denvm(m)
             thetav(k) = thetavm(m)

             qcv(k) = qcvm(m)
             qvv(k) = qvvm(m)
             ptotv(k) = ptotvm(m)
             tv(k) = tvm(m)
             qrv(k) = qrvm(m)
             qiv(k) = qivm(m)

             do icr = 1, ncr
                do ibr = 1, nbr
                   do ipr = 1, npr
                      qrpv(ipr,ibr,icr,k) = qrpvm(ipr,ibr,icr,m)
                   enddo
                   if(qrpv(1,ibr,icr,k).ne.0.0) then
                      k1br(ibr,icr) = min(k1br(ibr,icr),k)
                      k2br(ibr,icr) = max(k2br(ibr,icr),k)
                   end if
                enddo
             enddo
             do ici = 1, nci
                do ibi = 1, nbi
                   do ipi = 1, npi
                      qipv(ipi,ibi,ici,k) = qipvm(ipi,ibi,ici,m)
                   enddo
                   if(qipv(1,ibi,ici,k).ne.0.0) then
                      k1bi(ibi,ici) = min(k1bi(ibi,ici),k)
                      k2bi(ibi,ici) = max(k2bi(ibi,ici),k)
                   end if
                enddo
             enddo
             do ica = 1, nca
                do iba = 1, nba
                   do ipa = 1, npa
                      qapv(ipa,iba,ica,k) = qapvm(ipa,iba,ica,m)
                   enddo
                enddo
             enddo

             dmtendl(:,:,k) = dmtendlm(:,:,m)
             dcontendl(:,:,k) = dcontendlm(:,:,m)
             dbintendl(:,:,:,k) = dbintendlm(:,:,:,m)

          enddo

          ! mass tendencies output
          do m = 1, 10 ! liquid
             do k = KS, KE
                AMPS_mt(k,i,j,m) = dmtendl(m,1,k)
             enddo
          enddo
          do m = 11, 20 ! ice
             do k = KS, KE
                AMPS_mt(k,i,j,m) = dmtendl(m-10,2,k)
             enddo
          enddo
          do m = 1, 10 ! liquid
             do k = KS, KE
                AMPS_ct(k,i,j,m) = dcontendl(m,1,k)
             enddo
          enddo
          do m = 11, 20 ! ice
             do k = KS, KE
                AMPS_ct(k,i,j,m) = dcontendl(m-10,2,k)
             enddo
          enddo
          do m = 1, 3 ! bin-wise ice
             do k = KS, KE
                do ibr = 1, nbr
                   AMPS_bt(k,i,j,ibr,m) = dbintendl(m,1,ibr,k)
                enddo
                do ibi = 1, nbi
                   AMPS_bt(k,i,j,nbr+ibi,m) = dbintendl(m,2,ibi,k)
                enddo
             enddo
          enddo

          ! terminal velocity output
          do ibr = 1, nbr ! liquid
             do k = KS, KE
                do m = 2, 1, -1 ! mass and conc. weighted terminal velocities
                   AMPS_tv(k,i,j,ibr,m) = qrpv(npr-m+1,ibr,1,k) ! category 1
                enddo
             enddo
          enddo
          do ibi = 1, nbi ! ice
             do k = KS, KE
                do m = 2, 1, -1 ! mass and conc. weighted terminal velocities
                   AMPS_tv(k,i,j,nbr+ibi,m) = qipv(npi-m+1,ibi,1,k) ! category 1
                enddo
             enddo
          enddo

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.6.  Find limits
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
          call findcond1(k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i,micptrv,qrv,qiv &
                        ,nzh,qcv,k1etatv(:),ncr,nci,i,j,level,0,idir,'bf_sed')

! k1m = lower limit of grid index containing any hydrometeors
! k2m = upper limit of grid index containing any hydrometeors
! k1c = lower limit of grid index containing cloud hydrometeors (if level=5, cloud = rain)
! k2c = upper limit of grid index containing cloud hydrometeors (if level=5, cloud = rain)
! k1r = lower limit of grid index containing rain hydrometeors
! k2r = upper limit of grid index containing rain hydrometeors
! k1i = lower limit of grid index containing ice hydrometeors
! k2i = upper limit of grid index containing ice hydrometeors
! j = 1
! idir = 3
! kbnd = 0
! k1eta = 2

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.7.  Check total mass conservation
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
          if (amps_debug) then
             ! check total water content
             qtotal2 = 0.0_RP
             qtotal3 = 0.0_RP
             do k = KS, KE
                qtotal2(k) = qtotal2(k) + qvv(k)*moist_denv(k)
                do ibr=1,nbr
                   qtotal2(k) = qtotal2(k) + (qrpv(rmt_q,ibr,1,k) - qrpv(rmat_q,ibr,1,k))*moist_denv(k)
                enddo
                do ibi = 1, nbi
                   qtotal2(k) = qtotal2(k) + (qipv(imt_q,ibi,1,k) - qipv(imat_q,ibi,1,k))*moist_denv(k)
                   qtotal3(k) = qtotal3(k) + (qipv(imt_q,ibi,1,k) - qipv(imat_q,ibi,1,k) - qipv(imw_q,ibi,1,k))*moist_denv(k)
                enddo
             enddo

             do k = KS, KE
                den_diff(k) = den_diff(k) + qtotal2(k) - qtotal(k)
                if (abs(qtotal(k) - qtotal2(k)) > 0.001_RP*qtotal(k)) then
                   LOG_WARN("ATMOS_PHY_MP_amps_tendency",'(a,3I3,3ES15.6)') "CHECKTOTAL ERROR:", k, i, j, qtotal(k), qtotal2(k), qtotal3(k)
                endif
             enddo
          endif

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     2.8.  Energy (C_V t) tendency of AMPS microphysical processes before sedimentation
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------

          do k = KS, KE
             ! vapor difference
             Emoist(k,2) = - LHV0 * qvv(k) * moist_denv(k)

             ! ice difference
             do ibi = 1, nbi
                Emoist(k,2) = Emoist(k,2) &
                     + LHF0 * ( qipv(imt_q,ibi,1,k) - qipv(imw_q,ibi,1,k) - qipv(imat_q,ibi,1,k) ) * moist_denv(k)
             enddo
          enddo

          ! diabatic heating tendency, potential energy rho g h to be calculated in sedimentation later
          do k = KS, KE
             RHOE_t(k,i,j) = ( Emoist(k,2) - Emoist(k,1) ) / dt
!!$             RHOE_t(k,i,j) = RHOE_t(k,i,j) + (tv(k)* &
!!$                  (CVtot(k,i,j) + &
!!$                  CVtot_t(k,i,j)*dt) - &
!!$                  TEMP(k,i,j)*CVtot(k,i,j))* &
!!$                  DENS(k,i,j)/dt
          enddo

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     3.  Sedimentation
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
          Sedimention_IF: if(l_sediment) then
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     3.1.  Calculate mean masses for all the bins for AMPS
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
             if(level.ge.5) then
                do k = KS-1, KE
                   do icr=1,ncr
                   do ibr=1,nbr
                      mmassrv(ibr,icr,k) = 0.0_RP
                   enddo
                   enddo
                   do ici=1,nci
                   do ibi=1,nbi
                      mmassiv(ibi,ici,k) = 0.0_RP
                   enddo
                   enddo
                enddo
                if(k2r.gt.0) call cal_mmass(mmassrv &
                                           ,k1r,k2r,npr,nbr,ncr,qrpv,nzh &
                                           ,0)
                if(k2i.gt.0) call cal_mmass(mmassiv &
                                           ,k1i,k2i,npi,nbi,nci,qipv,nzh &
                                           ,1)
             end if
!
!--------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!     3.2 Implement sedimentation
!--------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!
             call PROF_rapstart("amps_sedimentation",3)
!
!-----------------------------------------------------------------------------------------------
!     3.2.1.     Now the liq fields
!-----------------------------------------------------------------------------------------------
!

             if(k2r.gt.0) then
                isn=0
                ! original
                call sclsedprz_original(qrv,qiv,qcv,trpv(:,2) &
                                       ,isprayr &
                                       ,nzh,1,1,1,npr,nbr,ncr,k1r,k2r,k1m,k2m &
                                       ,k1br,k2br &
                                       ,qrpv,moist_denv,momv,waccv &
                                       ,spdsfcv &
                                       ,thetav,qvv,thskinv,pgnd &
                                       ,k1etatv(1),CZ(KS-1:KE,i,j),FZ(KS-1:KE,i,j),dzzmv(:),dzvmv(:),dz1v &
                                       ,dt,i,j,advppmz(:,:,i,j),advppmze(:,i,j),isn,iadvv &
                                       ,level,nbhzcl,kmicv,imicv,jmicv &
                                       ,U(KS-1:KE,i,j) &          ! [IN]     u velocity,  cell face
                                       ,V(KS-1:KE,i,j) &          ! [IN]     v velocity,  cell face
                                       ,tv(:)                       &       ! [IN]     temperature
                                       ,DENS(KS-1:KE,i,j) &       ! [IN]     dry+moist density
                                       ,MOMZ(KS-1:KE,i,j) &       ! [IN]     z momentum
                                       ,CV_WATER &                          ! [IN]     specific heat water
                                       ,CV_ICE &                            ! [IN]     specific heat ice
                                       ,den_t(:)                      &     ! [INOUT]  density tendency
                                       ,MOMZ_t(KS-1:KE,i,j) &     ! [INOUT]  z momentum tendency
                                       ,RHOU_t(KS-1:KE,i,j) &     ! [INOUT]  rho u velocity tendency
                                       ,RHOV_t(KS-1:KE,i,j) &     ! [INOUT]  rho v velocity tendency
                                       ,RHOE_t(KS-1:KE,i,j) &     ! [INOUT]  rho e velocity tendency
                                       ,SFLX_rain(i,j))           ! [INOUT]  surface flux

                SFLX_rain(i,j) = abs(SFLX_rain(i,j))
             end if
!
!-----------------------------------------------------------------------------------------------
!     3.2.2.     Now the ice fields
!-----------------------------------------------------------------------------------------------
!

             if(k2i.gt.0) then
                isn=1
                call sclsedprz_original(qiv,qrv,qcv,trpv(:,2) &
                                       ,isprayi &
                                       ,nzh,1,1,1,npi,nbi,nci,k1i,k2i,k1m,k2m &
                                       ,k1bi,k2bi &
                                       ,qipv,moist_denv,momv,waccv &
                                       ,spdsfcv &
                                       ,thetav,qvv,thskinv,pgnd &
                                       ,k1etatv(1),CZ(KS-1:KE,i,j),FZ(KS-1:KE,i,j),dzzmv(:),dzvmv(:),dz1v &
                                       ,dt,i,j,advppmz(:,:,i,j),advppmze(:,i,j),isn,iadvv &
                                       ,level,nbhzcl,kmicv,imicv,jmicv &
                                       ,U(KS-1:KE,i,j) &          ! [IN]     u velocity,  cell face
                                       ,V(KS-1:KE,i,j) &          ! [IN]     v velocity,  cell face
                                       ,tv(:)                       &       ! [IN]     temperature
                                       ,DENS(KS-1:KE,i,j) &       ! [IN]     dry+moist density
                                       ,MOMZ(KS-1:KE,i,j) &       ! [IN]     z momentum
                                       ,CV_WATER &                          ! [IN]     specific heat water
                                       ,CV_ICE &                            ! [IN]     specific heat ice
                                       ,den_t(:)                      &     ! [INOUT]  density tendency
                                       ,MOMZ_t(KS-1:KE,i,j) &     ! [INOUT]  z momentum tendency
                                       ,RHOU_t(KS-1:KE,i,j) &     ! [INOUT]  rho u velocity tendency
                                       ,RHOV_t(KS-1:KE,i,j) &     ! [INOUT]  rho v velocity tendency
                                       ,RHOE_t(KS-1:KE,i,j) &     ! [INOUT]  rho e velocity tendency
                                       ,SFLX_snow(i,j))           ! [INOUT]  surface flux


                SFLX_snow(i,j) = abs(SFLX_snow(i,j))
             endif

!-----------------------------------------------------------------------------------------------
!     3.2.3.     Now the aerosol fields
!-----------------------------------------------------------------------------------------------
!
          ! NOT IN USE, NOT TESTED YET
          !if(level.eq.4.and.(k2r.gt.0.or.k2i.gt.0)) then
            ! for level 4 CLR scheme
            !   aerosol particles inside of hydrometeors are
            !   predicted with qapv.
          !  k1a=min(k1r,k1i)
          !  k2a=max(k2r,k2i)
          !  call sclsedaer(ica_sed1,ica_sed2 &
          !                ,nzh,1,1,npa,nba,nca,k1a,k2a &
          !                ,qapv,moist_denv,waccv &
          !                ,k1etatv(1),REAL_CZ(KS-1:KE,i,j),FZ(KS-1:KE,i,j),dzzmv(:),dz1v &
          !                ,dt,i,j,advppmz(:,:,i,j),advppmze(:,i,j) &
          !                ,iadvv &
          !                ,level,kmicv,imicv,jmicv)
          !end if

             call PROF_rapend("amps_sedimentation",3)

!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     3.3   shift bins
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
! bin no.     1        2        3
!         |       |        |        |
!     binbi(1) binbi(2) binbi(3) binbi(4)
!           +dM_1    +dM_2    +dM_3
! The shift of boundaries (left and right) in a particular bin(k) due to the
! tendency term are equal to max(binbi(k) + dM_k) and max(binbi(k) + dM_k),
! respectively.

             if(level.eq.5) then
                iph=1
                call shift_bin_vec(nzh,npr,nbr,ncr,k1r,k2r &
                                  ,binbr,qrpv,mmassrv &
                                  ,moist_denv,iph,kmicv,imicv,jmicv)
                iph=2
                call shift_bin_vec(nzh,npi,nbi,nci,k1i,k2i &
                                  ,binbi,qipv,mmassiv &
                                  ,moist_denv,iph,kmicv,imicv,jmicv)
             endif
          endif Sedimention_IF

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     3.4.     Negative adjustment
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
          if(level.eq.5) then
!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     3.4.1   AMPS
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!

!    ircnfl = 1   ipcnfl = 1    iscnfl = 0    iacnfl = 0    igcnfl = 0
             call negadj_bin2(nzh,nzh,npr,nbr,ncr,npi,nbi,nci,i,j &
                             ,trpv(:,2),trpv(:,1),pbv(:),thetav(:),qvv(:) &
                             ,qcv,qrpv,qipv &
                             ,piv &
                             ,1,1 &
                             ,imicv,jmicv,kmicv)

             call negadj_bin_ap(nzh,npa,nba,nca &
                               ,npr,nbr,ncr,npi,nbi,nci &
                               ,j,k &
                               ,qapv,qrpv,qipv,moist_denv &
                               ,1,1 &
                               ,imicv,jmicv,kmicv)

          end if

       else

          do k = KS, KE
             RHOE_t(k,i,j) = 0.0_RP
          end do

       endif MIC_IF

!
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!     3.5.   Re diagnose the thermo variables
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!
      ! diagnose thermodynamical variables
       ivis=1
       call moistthermo2(thetav,moist_denv,tv,ptotv  &
                        ,trpv(:,2),qvv,qcv,qrv,qiv,micptrv   &
                        ,npr,nbr,ncr,npi,nbi,nci   &
                        ,level,nbhzcl   &
                        ,trpv(:,1),qrpv,qipv, pbv  &
                        ,estbar,esitbar, piv &
                        ,nmic,imicv,jmicv,kmicv &
                        ,ivis,'final',isect)

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     3.6.  Find limits
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
!

       call findcond2(k1m,k2m,k1c,k2c,k1r,k2r,k1i,k2i,micptrv   &
                     ,k1br,k2br,k1bi,k2bi   &
                     ,qrv,qiv   &
                     ,nzh,npr,nbr,ncr,npi,nbi,nci,level,nbhzcl   &
                     ,qcv,qrpv,qipv,trpv(:,2),qrov,qiov   &
                     ,k1etatv,i,j,'find2',idir)


!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     3.7.  Fill aerosol variables
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
!

       if(l_fill_aerosols) then
          ! fill the cloud-free environment with reference values
          if(level.eq.5) then
             call fill_bin_ap(nzh,npa,nba,nca &
                             ,npr,nbr,ncr,npi,nbi,nci  &
                             ,qapv,qrpv,qipv,qapv_ini(:,:,:,:,i,j),moist_denv,den_ini(:,i,j) &
                             ,1,nzh,imicv,jmicv,kmicv)
          endif
       endif

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     3.8.  Calculate mean masses for all the bins before sending for dynamics
!-------------------------------------------------------------------------
!--------------------------------------------------------------------------
!
       if(level.ge.5) then
          do k = KS, KE
             do icr = 1, ncr
             do ibr = 1, nbr
                mmassrv_global(ibr,icr,k,i,j) = 0.0_RP
             enddo
             enddo
             do ici = 1, nci
             do ibi = 1, nbi
                mmassiv_global(ibi,ici,k,i,j) = 0.0_RP
             enddo
             enddo
          enddo
          call cal_mmass_scale(mmassrv_global(:,:,:,i,j),nzh, &
                              2,nzh,npr,nbr,ncr,qrpv, &
                              0)
          call cal_mmass_scale(mmassiv_global(:,:,:,i,j),nzh, &
                              2,nzh,npi,nbi,nci,qipv, &
                              1)
       end if
!mark4
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!     4. Calculate tendencies
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
       ! tendency calculation notebook:
       ! X(n+1) = [ X(n) DEN(n) + RHOQ dt ] / DEN(n+1)
       ! RHOQ = [ X(n+1) DEN(n+1) - X(n) DEN(n) ] / dt
       ! X is SCALE variable, x is AMPS variable
       ! -----  mass  -----
       ! x is mixing ratio, its unit in AMPS is g cm-3/g cm-3
       ! X is mixing ratio, its unit in SCALE is kg m-3/kg m-3
       ! X is SCALE variable, x is AMPS variable
       ! no conversion needs to be done, as they are unitless, i.e. x = X
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1)  - X(n) DEN(n) ] / dt
       ! -----PPV mass-----
       ! x is mixing ratio, its unit in AMPS is g cm-3/g cm-3
       ! X is density, its unit in SCALE is kg m-3
       ! X is SCALE variable, x is AMPS variable
       ! conversion from x -> X must be performed by multiplying DEN (kg m-3)
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1) DEN(n+1) - X(n) DEN(n) ] / dt
       ! -----non-mass-----
       ! x is PPV per unit volume (cm-3) per density of moist air DEN (g cm-3)
       ! X is PPV per unit volume of cm-3, it is treated as a tracer
       ! conversion from x -> X must be performed by multiplying DEN (g cm-3)
       ! 1 kg m-3 = 0.001 g cm-3, thus X(n+1) = x(n+1)*DEN(n+1)*0.001
       ! thus tendency is RHOQ = [ x(n+1) DEN(n+1) 0.001 DEN(n+1)  - X(n) DEN(n) ] / dt

       ! density tendency, this includes all ice and liquid mass tracers,
       ! excluding PPVs and aerosols
       ! this only considers the density change due to sedimentation,
       ! since the total density is conserved in a cell in microphysics processes

       do k = KS, KE
          DENS_t(k,i,j) = DENS_t(k,i,j) + den_t(k)
          DENS_NEW(k) = DENS(k,i,j) + DENS_t(k,i,j)*dt
       enddo

       ! vapor tendency, mixing ratio -> density
       do k = KS, KE
          RHOQ_t(k,i,j,I_QV) = (qvv(k)*moist_denv(k) - &
               QTRC(k,i,j,I_QV)*DENS(k,i,j))/dt
       enddo

       ! liquid tendency
       do k = KS, KE
          do icr = 1, ncr
             ipr_qpr = 0
             do ibr = 1, nbr
                RHOQ_t(k,i,j,I_QL+ibr-1) = ((qrpv(rmt_q,ibr,icr,k) - &
                     qrpv(rmat_q,ibr,icr,k))*moist_denv(k) -   &
                     QTRC(k,i,j,I_QL+ibr-1)*DENS(k,i,j))/dt

                do ipr = rmat_q, rmas_q
                   RHOQ_t(k,i,j,I_QPPVL+ipr_qpr) = &
                        (qrpv(ipr,ibr,icr,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVL+ipr_qpr)*DENS(k,i,j))/dt
                   ipr_qpr = ipr_qpr + 1
                enddo
                ! number concentration, non-mass variable
                RHOQ_t(k,i,j,I_QPPVL+ipr_qpr) = &
                     (qrpv(rcon_q,ibr,icr,k)*moist_denv(k)*0.001_RP - &
                     QTRC(k,i,j,I_QPPVL+ipr_qpr)*DENS(k,i,j))/dt
                ipr_qpr = ipr_qpr + 1
             enddo
          enddo
       enddo

       ! ice tendency
       do k = KS, KE
          do ici = 1, nci
             ipi_qi = 0
             do ibi = 1, nbi
                ! total ice mass
                RHOQ_t(k,i,j,I_QI+ibi-1) = &
                     ((qipv(imt_q,ibi,ici,k) - qipv(imat_q,ibi,ici,k) - &
                     qipv(imw_q,ibi,ici,k))*moist_denv(k) - &
                     QTRC(k,i,j,I_QI+ibi-1)*DENS(k,i,j))/dt

                ! melt water mass
                RHOQ_t(k,i,j,I_QW+ibi-1) = &
                     (qipv(imw_q,ibi,ici,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QW+ibi-1)*DENS(k,i,j))/dt

                ! rimed ,aggregate, and crystal
                do ipi = imr_q, imc_q
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ipi,ibi,ici,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1
                enddo
                ! frozen water
                RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                     (qipv(imf_q,ibi,ici,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                ipi_qi = ipi_qi + 1

                ! aerosol mass
                do ipi = imat_q, imas_q
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ipi,ibi,ici,k)*moist_denv(k) - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1
                enddo

                ! number concentration and other non-mass variables
                if (l_gaxis_version == 1) then
                   ! -- ORIGINAL
                   do ipi = imas_q+1, npi-2
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(ipi,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt

                      if ( ipi == iag_q .and. l_axis_limit .and. qipv(iacr_q,ibi,ici,k) > EPS ) then
                         if ( qipv(iag_q,ibi,ici,k) / qipv(iacr_q,ibi,ici,k) > 1.0_RP) then
                            RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                                 (qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP &
                                 - QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         endif
                      endif

                      if ( ipi == icg_q .and. l_axis_limit .and. qipv(iccr_q,ibi,ici,k) > EPS ) then
                         if ( qipv(icg_q,ibi,ici,k) / qipv(iccr_q,ibi,ici,k) > 1.0_RP) then
                            RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                                 (qipv(iccr_q,ibi,ici,k)*moist_denv(k)*0.001_RP &
                                 - QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         endif
                      endif

                      ipi_qi = ipi_qi + 1
                   enddo

                else if (l_gaxis_version == 2) then
                   ! -- METHOD 1
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(icon_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ivcs_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   if (qipv(iacr_q,ibi,ici,k) < 1.e-22_RP) then
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (0.0_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1

                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (0.0_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   else
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(iccr_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                           qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1

                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(idcr_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                           qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   endif

                   if ( l_axis_limit .and. (qipv(iag_q,ibi,ici,k)/qipv(iacr_q,ibi,ici,k))**(1.0_RP/3.0_RP) > 1.0_RP) then
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   else
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (qipv(iag_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   endif

                   if (qipv(iag_q,ibi,ici,k) < 1.e-22_RP) then
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (0.0_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   else
                      if ( l_axis_limit .and. (qipv(icg_q,ibi,ici,k)/qipv(iccr_q,ibi,ici,k))**(1.0_RP/3.0_RP) > 1.0_RP) then
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(iccr_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                              qipv(iag_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         ipi_qi = ipi_qi + 1
                      else
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(icg_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                              qipv(iag_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                         ipi_qi = ipi_qi + 1
                      endif
                   endif

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(inex_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                else if (l_gaxis_version == 3) then
                   ! -- METHOD 2
                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(icon_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(ivcs_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(iccr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(idcr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                   if (qipv(iacr_q,ibi,ici,k) < 1.e-22_RP) then
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (0.0_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   else
                      if ( l_axis_limit .and. (qipv(iag_q,ibi,ici,k)/qipv(iacr_q,ibi,ici,k))**(1.0_RP/3.0_RP) > 1.0_RP) then
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(icon_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      else
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(iag_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                              qipv(iacr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      endif
                      ipi_qi = ipi_qi + 1
                   endif

                   if (qipv(iccr_q,ibi,ici,k) < 1.e-22_RP) then
                      RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                           (0.0_RP - &
                           QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      ipi_qi = ipi_qi + 1
                   else
                      if ( l_axis_limit .and. (qipv(icg_q,ibi,ici,k)/qipv(iccr_q,ibi,ici,k))**(1.0_RP/3.0_RP) > 1.0_RP) then
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(icon_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      else
                         RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                              (qipv(icg_q,ibi,ici,k)*qipv(icon_q,ibi,ici,k)/ &
                              qipv(iccr_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                              QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                      endif
                      ipi_qi = ipi_qi + 1
                   endif

                   RHOQ_t(k,i,j,I_QPPVI+ipi_qi) = &
                        (qipv(inex_q,ibi,ici,k)*moist_denv(k)*0.001_RP - &
                        QTRC(k,i,j,I_QPPVI+ipi_qi)*DENS(k,i,j))/dt
                   ipi_qi = ipi_qi + 1

                endif

             enddo
          enddo
       enddo

       ! aerosol tendency
       do k = KS, KE
          ipa_qpa = 0
          do ica = 1, nca
             do iba = 1, nba
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa) = &
                     (qapv(amt_q,iba,ica,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QPPVA+ipa_qpa)*DENS(k,i,j))/dt
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa+1) = &
                     (qapv(acon_q,iba,ica,k)*moist_denv(k)*0.001_RP - &
                     QTRC(k,i,j,I_QPPVA+ipa_qpa+1)*DENS(k,i,j))/dt
                RHOQ_t(k,i,j,I_QPPVA+ipa_qpa+2) = &
                     (qapv(ams_q,iba,ica,k)*moist_denv(k) - &
                     QTRC(k,i,j,I_QPPVA+ipa_qpa+2)*DENS(k,i,j))/dt
                ipa_qpa = ipa_qpa + 3
             enddo
          enddo
       enddo

       ! specific heat tendency, loop over liquid and ice mass tendency, (mixing ratio)
       do k = KS, KE
          ! vapor difference
          dq = qvv(k) * moist_denv(k) / DENS_NEW(k) - QTRC(k,i,j,I_QV)
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_VAPOR * dq / dt
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_VAPOR * dq / dt

          ! liquid difference
          dq = 0.0_RP
          do ibr = 1, nbr
             ! liquid drop mass
             dq = dq + ( qrpv(rmt_q,ibr,1,k) - qrpv(rmat_q,ibr,1,k) ) * moist_denv(k) / DENS_NEW(k) &
                       - QTRC(k,i,j,I_QL+ibr-1)
          enddo
          do ibi = 1, nbi
             ! melt water mass
             dq = dq  + qipv(imw_q,ibi,1,k) * moist_denv(k) / DENS_NEW(k) - QTRC(k,i,j,I_QW+ibi-1)
          enddo
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_WATER * dq / dt
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_WATER * dq / dt

          ! ice difference
          dq = 0.0_RP
          do ibi = 1, nbi
             dq = dq + ( qipv(imt_q,ibi,1,k) - qipv(imw_q,ibi,1,k) - qipv(imat_q,ibi,1,k) ) * moist_denv(k) / DENS_NEW(k) &
                       - QTRC(k,i,j,I_QI+ibi-1)
          enddo
          CPtot_t(k,i,j) = CPtot_t(k,i,j) + CP_ICE * dq / dt
          CVtot_t(k,i,j) = CVtot_t(k,i,j) + CV_ICE * dq / dt
       enddo

       ! momentum flux, rhou_t, rhov_t, and surface flux were calculated in the sedimentation process

       do icr = 1, ncr
       do ibr = 1, nbr
       do k = KS, IE
          QLIQ(k,i,j) = QLIQ(k,i,j) + qrpv(rmt_q,ibr,icr,k) * moist_denv(k) / DENS_NEW(k)
       end do
       end do
       end do

       do ici = 1, nci
       do ibi = 1, nbi
       do k = KS, IE
          ! total
          QICE(k,i,j,1) = QICE(k,i,j,1) + qipv(imt_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
          ! rime
          QICE(k,i,j,2) = QICE(k,i,j,2) + qipv(imr_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
          ! aggregate
          QICE(k,i,j,3) = QICE(k,i,j,3) + qipv(ima_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
          ! crystal
          QICE(k,i,j,4) = QICE(k,i,j,4) + qipv(imc_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
          ! melt water
          QICE(k,i,j,5) = QICE(k,i,j,5) + qipv(imw_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
          ! frozen water
          QICE(k,i,j,6) = QICE(k,i,j,6) + qipv(imf_q,ibi,ici,k) * moist_denv(k) / DENS_NEW(k)
       end do
       end do
       end do


    enddo
    enddo

    TIME_AMPS = TIME_AMPS + 1

    do i = 1, 20
       call FILE_HISTORY_in( AMPS_mt(:,:,:,i), AMPS_mt_NAME(i), &
                             AMPS_t_DESC(i), AMPS_t_UNIT(1), fill_halo=.true. )
    enddo
    do i = 1, 20
       call FILE_HISTORY_in( AMPS_ct(:,:,:,i), AMPS_ct_NAME(i), &
                             AMPS_t_DESC(i), AMPS_t_UNIT(2), fill_halo=.true. )
    enddo
    do i = 1, ATMOS_PHY_MP_amps_nwaters
       if (i < 10) then
          write (AMPS_bin_NAME, "(I1)") i
       else
          write (AMPS_bin_NAME, "(I2)") i
       endif
       call FILE_HISTORY_in( AMPS_bt(:,:,:,i,1), trim(AMPS_bt_NAME(1))//trim(AMPS_bin_NAME), &
                             'Bin-wise riming mass rate', 'g /cm3/s ' , fill_halo=.true. )
       call FILE_HISTORY_in( AMPS_bt(:,:,:,i,2), trim(AMPS_bt_NAME(2))//trim(AMPS_bin_NAME), &
                             'Bin-wise vapor mass rate', 'g /cm3/s ' , fill_halo=.true. )
       call FILE_HISTORY_in( AMPS_bt(:,:,:,i,3), trim(AMPS_bt_NAME(3))//trim(AMPS_bin_NAME), &
                             'Bin-wise evaporation mass rate', 'g /cm3/s ' , fill_halo=.true. )
    enddo

    do i = 1, ATMOS_PHY_MP_amps_nwaters
       if (i < 10) then
          write (AMPS_tv_NAME, "(I1)") i
       else
          write (AMPS_tv_NAME, "(I2)") i
       endif
       if (i <= ATMOS_PHY_MP_amps_nwaters - ATMOS_PHY_MP_amps_nices) then
          AMPS_tv_NAME = 'liq_'//trim(AMPS_tv_NAME)
       else
          AMPS_tv_NAME = 'ice_'//trim(AMPS_tv_NAME)
       endif

       call FILE_HISTORY_in( AMPS_tv(:,:,:,i,1), 'mass_terminalV_'//trim(AMPS_tv_NAME), &
                             'Mass weighted terminal velocity', 'm/s', fill_halo=.true. )
       call FILE_HISTORY_in( AMPS_tv(:,:,:,i,2), 'conc_terminalV_'//trim(AMPS_tv_NAME), &
                             'Concentration weighted terminal velocity', 'm/s', fill_halo=.true. )
    enddo

    call FILE_HISTORY_in( QLIQ(:,:,:), 'TOTAL_LIQ', &
                         'Mass of total liquid water', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,1), 'TOTAL_ICE', &
                         'Mass of total ice', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,2), 'RIME', &
                         'Mass of rime ice', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,3), 'AGGREGATE', &
                         'Mass of aggregate ice', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,4), 'CRYSTAL', &
                         'Mass of crystal ice', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,5), 'MELT_WATER', &
                         'Mass of melt water', 'kg/kg', fill_halo=.true. )
    call FILE_HISTORY_in( QICE(:,:,:,6), 'FROZEN_WATER', &
                         'Mass of frozon water', 'kg/kg', fill_halo=.true. )

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_tendency
  !____________________________________________________________________________________


  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_amps_cloud_fraction( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       QTRC,           &
       mask_criterion, &
       cldfrac         )
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: QTRC(KA,IA,JA,QA-1)
    real(RP), intent(in)  :: mask_criterion

    real(RP), intent(out) :: cldfrac(KA,IA,JA)

#ifndef DISABLE_AMPS
    real(RP) :: qhydro
    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    !$omp parallel do collapse(2) default(none) &
    !$omp private(i,j,k,iq,qhydro) &
    !$omp shared(JS,JE,IS,IE,KS,KE, &
    !$omp        nbr,nbi,I_QI,I_QW, &
    !$omp        mask_criterion, &
    !$omp        cldfrac,QTRC)
    do j  = JS, JE
    do i  = IS, IE
    do k  = KS, KE
       qhydro = 0.0_RP
       do iq = 1, nbr
          qhydro = qhydro + QTRC(k,i,j,iq)
       enddo
       do iq = 1, nbi
          qhydro = qhydro + QTRC(k,i,j,I_QI+iq-2) + QTRC(k,i,j,I_QW+iq-2)
       enddo

       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP,qhydro-mask_criterion)
    enddo
    enddo
    enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_cloud_fraction
  !____________________________________________________________________________________


  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_amps_effective_radius( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS0,      &
       TEMP0,      &
       QTRC0,      &
       Re          )
    use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in)  :: DENS0(KA,IA,JA)      ! density      [IN]
    real(RP), intent(in)  :: TEMP0(KA,IA,JA)      ! temperature  [IN]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA-1) ! hydrometeor mixing ratio (excluding vapor) [IN]
    real(RP), intent(out) :: Re(KA,IA,JA,N_HYD)   ! hydrometeor effective radius (cm) [OUT]

#ifndef DISABLE_AMPS
    ! following AMPS, the density of water is default 1000 kg/m3
    real(RP), parameter :: den_water=0.001_RP      ! kg/cm3
    ! following AMPS, the density of ice is default 916.68 kg/m3
    real(RP), parameter :: den_ice=0.00091668_RP   ! kg/cm3

    ! lowest limitation, a tuning value
    real(DP), parameter :: RRLMTB=1.0d-22 ! 18
    real(DP), parameter :: RILMTB=1.0d-22 ! 18
    real(DP), parameter :: RLMTB=1.0d-15 ! 10

    real(RP), parameter :: multiplicativeFactor = 4.0_RP*PI/3.0_RP
    real(RP), parameter :: cylinderFactor = 2.0_RP*PI
    real(RP), parameter :: oneThirdFactor = 1.0_RP/3.0_RP
    real(RP), parameter :: twoThirdFactor = 2.0_RP/3.0_RP
    real(RP), parameter :: M32CM3 = 1000000.0_RP

    real(RP) :: rho_ice, rho_rim, rho_agg, rho_cry, rho_con, rho_a, rho_c, rho_vol, rho_nex, alpha, kp
    real(RP) :: temp1, temp2
    real(RP) :: kcsw, vcs, vcsw, vim, vspace, semi_a, semi_c, micore
    real(RP) :: Rmax, sum3, sum2, conc, sum_q
    real(RP) :: r3, factor
    integer  :: i, j, k, iq, iq1, iq2
    integer  :: ic
    integer  :: liqConc_index, iceConc_index, rim_index, agg_index, cry_index, snow_type
    integer  :: a_index, c_index, vol_index, nexice_index


    liqConc_index = I_QPPVL + 2 - 1
    iceConc_index = I_QPPVI + 6 - 1

    !$omp parallel do collapse(2) default(none) &
    !$omp private(i,j,k,iq, &
    !$omp         sum2,sum3,conc,sum_q,rho_con,r3) &
    !$omp shared(JS,JE,IS,IE,KS,KE,split_bins,nbr, &
    !$omp        liqConc_index,numberPPVL, &
    !$omp        Re,QTRC0,DENS0)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       ! effective radius of cloud water
       sum3 = 0.0_RP
       sum2 = 0.0_RP
       conc = 0.0_RP
       sum_q = 0.0_RP
       do iq = 1, split_bins
          rho_con = QTRC0(k,i,j,liqConc_index+(iq-1)*numberPPVL) * DENS0(k,i,j)
          if ( rho_con         <= RRLMTB .or. &
               QTRC0(k,i,j,iq) <= RRLMTB) cycle
          ! r = ( mass / den_w / (4 pi / 3) )^1/3
          ! need to convert from m-3 to cm-3
          r3 = QTRC0(k,i,j,iq) / QTRC0(k,i,j,liqConc_index+(iq-1)*numberPPVL) &
               / M32CM3 / den_water / multiplicativeFactor
          ! r^3
          sum3 = sum3 + r3 * rho_con
          ! r^2
          sum2 = sum2 + r3**twoThirdFactor * rho_con

          conc = conc + rho_con
          sum_q = sum_q + QTRC0(k,i,j,iq)
       end do

       if ( conc > 0.0_RP .and. sum_q > RLMTB ) then
          Re(k,i,j,I_HC) = sum3/sum2
       else
          Re(k,i,j,I_HC) = 0.0_RP
       endif

       ! effective radius of rain water
       sum3 = 0.0_RP
       sum2 = 0.0_RP
       conc = 0.0_RP
       sum_q = 0.0_RP
       do iq = split_bins+1, nbr
          rho_con = QTRC0(k,i,j,liqConc_index+(iq-1)*numberPPVL) * DENS0(k,i,j)
          if ( rho_con         <= RRLMTB .or. &
               QTRC0(k,i,j,iq) <= RRLMTB) cycle
          ! r = ( mass / den_w / (4 pi / 3) )^1/3
          ! need to convert from m-3 to cm-3
          r3 = QTRC0(k,i,j,iq) / QTRC0(k,i,j,liqConc_index+(iq-1)*numberPPVL) &
               / M32CM3 / den_water / multiplicativeFactor
          ! r^3
          sum3 = sum3 + r3 * rho_con
          ! r^2
          sum2 = sum2 + r3**twoThirdFactor

          conc = conc + rho_con
          sum_q = sum_q + QTRC0(k,i,j,iq)
       end do

       if ( conc > 0.0_RP .and. sum_q > RLMTB ) then
          Re(k,i,j,I_HR) = sum3/sum2                                                ! cm
       else
          Re(k,i,j,I_HR) = 0.0_RP
       endif

    enddo
    enddo
    enddo


    rim_index    = I_QPPVI     - 1 ! rimed mass index
    agg_index    = I_QPPVI + 1 - 1 ! aggregate mass index
    cry_index    = I_QPPVI + 2 - 1 ! crystal mass index
    a_index      = I_QPPVI + 8 - 1 ! a-axis index
    c_index      = I_QPPVI + 9 - 1 ! c-axis index
    vol_index    = I_QPPVI + 7 - 1 ! circumscribing volume index
    nexice_index = I_QPPVI + 13 - 1 ! polycrystal index

    ! effective radius of snow
    if (l_reff_version == 1) then
       ! under construction, calculate maximum dimension according to their shape

       !$omp parallel do collapse(2) default(none) &
       !$omp private(i,j,k,ic,iq,iq1,iq2, &
       !$omp         sum2,sum3,conc,sum_q,rho_rim,rho_agg,rho_cry,rho_ice,rho_con, &
       !$omp         rho_a,rho_c,rho_vol,rho_nex,snow_type,alpha,kp,temp1,temp2, &
       !$omp         vcs,micore,semi_a,semi_c,vim,vspace,kcsw,vcsw,Rmax,factor) &
       !$omp shared(JS,JE,IS,IE,KS,KE,nbi,ibin_is,ibin_sg,ibin_gh, &
       !$omp        I_QI,I_QW,rim_index,agg_index,cry_index,iceConc_index,a_index,c_index, &
       !$omp        vol_index,nexice_index,numberPPVI, &
       !$omp        l_gaxis_version, &
       !$omp        Re,QTRC0,DENS0)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          do ic = 1, 4

             if ( ic == 1 ) then
                iq1 = 1
                iq2 = ibin_is
             else if ( ic == 2 ) then
                iq1 = ibin_is + 1
                iq2 = ibin_sg
             else if ( ic == 3 ) then
                iq1 = ibin_sg + 1
                iq2 = ibin_gh
             else
                iq1 = ibin_gh + 1
                iq2 = nbi
             end if

             sum3 = 0.0_RP
             sum2 = 0.0_RP
             conc = 0.0_RP
             sum_q = 0.0_RP

             do iq = iq1, iq2
                rho_rim = QTRC0(k,i,j,rim_index+(iq-1)*numberPPVI)     * DENS0(k,i,j)  ! kg/m3
                rho_agg = QTRC0(k,i,j,agg_index+(iq-1)*numberPPVI)     * DENS0(k,i,j)  ! kg/m3
                rho_cry = QTRC0(k,i,j,cry_index+(iq-1)*numberPPVI)     * DENS0(k,i,j)  ! kg/m3
                rho_ice = QTRC0(k,i,j,I_QI+iq-2)                       * DENS0(k,i,j)  ! kg/m3
                rho_con = QTRC0(k,i,j,iceConc_index+(iq-1)*numberPPVI) * DENS0(k,i,j)  ! /cm3
                rho_a   = QTRC0(k,i,j,a_index+(iq-1)*numberPPVI)       * DENS0(k,i,j)  ! cm3/cm3
                if (l_gaxis_version == 1 .or.l_gaxis_version == 3) then
                   rho_c   = QTRC0(k,i,j,c_index+(iq-1)*numberPPVI)    * DENS0(k,i,j)  ! cm3/cm3
                else
                   if (QTRC0(k,i,j,iceConc_index+(iq-1)*numberPPVI) < 1.e-22_RP) then
                      rho_c   = 0.0_RP                                                 ! cm3/cm3
                   else
                      rho_c   = QTRC0(k,i,j,c_index+(iq-1)*numberPPVI) * &
                                QTRC0(k,i,j,a_index+(iq-1)*numberPPVI) / &
                                QTRC0(k,i,j,iceConc_index+(iq-1)*numberPPVI) * DENS0(k,i,j)  ! cm3/cm3
                   endif
                endif
                rho_vol = QTRC0(k,i,j,vol_index+(iq-1)*numberPPVI)     * DENS0(k,i,j)  ! cm3/cm3
                rho_nex = QTRC0(k,i,j,nexice_index+(iq-1)*numberPPVI)  * DENS0(k,i,j)  ! /cm3
                if ( rho_ice <= RILMTB .or. &
                     rho_con <= RILMTB .or. &
                     rho_a   <= RILMTB .or. &
                     rho_c   <= RILMTB .or. &
                     rho_vol <= RILMTB) then
                   cycle
                endif
                ! first, we decide the type of snow
                ! if rimed mass is greater than 10% of total mass
                if ( rho_rim > 0.1_RP*rho_ice ) then
                   ! if 10% of total mass is greater than aggregate mass
                   if ( 0.1_RP*rho_ice > rho_agg ) then
                      ! if crystal mass is greater than rimed mass
                      if ( rho_cry > rho_rim ) then
                         ! rimed crystal
                         snow_type = 3
                      else
                         ! graupel
                         snow_type = 5
                      endif
                   else
                      ! if rimed mass is greater than aggregate mass
                      if ( rho_rim > rho_agg ) then
                         ! graupel
                         snow_type = 5
                      else
                         ! rimed aggregate
                         snow_type = 4
                      endif
                   endif
                else
                   ! if 10% of total mass is greater than aggregate mass
                   if ( 0.1_RP*rho_ice > rho_agg ) then
                      ! pristine crystal
                      snow_type = 1
                   else
                      ! aggregate
                      snow_type = 2
                   endif
                endif

                ! calculate the aspect ratio
                if ( snow_type == 1 ) then ! pristine crystal
                   if ( rho_nex / rho_con < 0.5_RP ) then
                      alpha = ( rho_a / rho_c )**oneThirdFactor
                   else
                      alpha = 0.25_RP
                   endif
                else if ( snow_type == 2 .or. snow_type == 4 ) then ! aggregate or rimed aggregate
                   temp1 = rho_agg + rho_cry
                   temp2 = 10.0_RP * rho_cry
                   if ( temp1 >= temp2 ) then
                      alpha = oneThirdFactor
                   else
                      alpha = ( rho_a / rho_c )**oneThirdFactor
                      kp = log10(oneThirdFactor/alpha)
                      alpha = ((temp1/temp2)**kp)/3.0_RP
                   endif
                else if ( snow_type == 3 ) then ! rimed crystal
                   ! SUBJECT TO MODIFICATION IN THE FUTURE,
                   ! FOR NOW ALPHA FOR RIMED CRYSTALS IS THE SAME AS FOR CRYSTALS
                   if ( rho_nex / rho_con < 0.5_RP  ) then
                      alpha = ( rho_a / rho_c )**oneThirdFactor
                   else
                      alpha = 0.25_RP
                   endif
                else if ( snow_type == 5 ) then ! graupel
                   temp1 = rho_agg + rho_cry + rho_rim
                   temp2 = 10.0_RP * (rho_Cry + rho_agg)
                   if ( temp1 >= temp2 ) then
                      alpha = 0.8_RP
                   else
                      alpha = ( rho_a / rho_c )**oneThirdFactor
                      kp = log10(0.8_RP/alpha)
                      alpha = ((temp1/temp2)**kp)*0.8_RP
                   endif
                endif

                !                                                            unit:
                vcs = rho_vol / rho_con                                    ! cm3
                temp1 = QTRC0(k,i,j,I_QW+iq-2) * DENS0(k,i,j) / &          ! kg/m3
                     rho_con / M32CM3                                   ! kg
                if ( snow_type == 5 ) then ! graupel
                   factor = cylinderFactor
                   semi_a = max(1.0_RP, alpha**3)
                else
                   factor = multiplicativeFactor
                   semi_a = sqrt(1.0_RP + alpha**2)**3
                endif
                semi_a = (vcs / factor / semi_a)**oneThirdFactor           ! cm
                semi_c = semi_a/alpha                                      ! cm
                micore = rho_ice &                                         ! kg/m3
                     / rho_con / M32CM3 &                                  ! kg
                     - temp1
                vim = factor*semi_c*semi_a**2                              ! cm3
                vspace = max(0.0_RP,vim - micore/den_ice)                  ! cm3
                kcsw = 1.0_RP + max(0.0_RP, temp1/den_water - vspace)/vim
                vcsw = kcsw*vcs                                            ! cm3

                Rmax = (vcsw/multiplicativeFactor)**oneThirdFactor         ! cm

                sum3 = sum3 + Rmax**3*rho_con    ! cm3/cm3
                sum2 = sum2 + Rmax**2*rho_con    ! cm2/cm3
                conc = conc + rho_con            ! /cm3
                sum_q = sum_q + QTRC0(k,i,j,I_QI+iq-2)

             end do ! iq

             if ( conc > 0.0_RP .and. sum_q > RLMTB ) then
                Re(k,i,j,I_HI+ic-1) = sum3/sum2
             else
                Re(k,i,j,I_HI+ic-1) = 0.0_RP
             endif

          end do ! ic

       enddo
       enddo
       enddo

    else if (l_reff_version == 2) then
       ! assume all ice crystals have constant density 0.91668 g cm-3 and calculate their maximum dimension

       !$omp parallel do collapse(2) default(none) &
       !$omp private(i,j,k,ic,iq,iq1,iq2, &
       !$omp         sum2,sum3,conc,sum_q,rho_ice,rho_con, &
       !$omp         vcsw,Rmax) &
       !$omp shared(JS,JE,IS,IE,KS,KE,nbi,ibin_is,ibin_sg,ibin_gh, &
       !$omp        I_QI,rim_index,agg_index,cry_index,iceConc_index, &
       !$omp        numberPPVI, &
       !$omp        Re,QTRC0,DENS0)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          do ic = 1, 4

             if ( ic == 1 ) then
                iq1 = 1
                iq2 = ibin_is
             else if ( ic == 2 ) then
                iq1 = ibin_is + 1
                iq2 = ibin_sg
             else if ( ic == 3 ) then
                iq1 = ibin_sg + 1
                iq2 = ibin_gh
             else
                iq1 = ibin_gh + 1
                iq2 = nbi
             end if

             sum3 = 0.0_RP
             sum2 = 0.0_RP
             conc = 0.0_RP
             sum_q = 0.0_RP

             do iq = iq1, iq2

                rho_ice = QTRC0(k,i,j,I_QI+iq-2)                       * DENS0(k,i,j)  ! kg/m3
                rho_con = QTRC0(k,i,j,iceConc_index+(iq-1)*numberPPVI) * DENS0(k,i,j)  ! /cm3

                if ( rho_ice <= RILMTB .or. &
                     rho_con <= RILMTB) then
                   cycle
                endif

                vcsw = rho_ice / M32CM3 / rho_con / den_ice

                Rmax = (vcsw/multiplicativeFactor)**oneThirdFactor                    ! cm

                sum3 = sum3 + Rmax**3*rho_con    ! cm3/cm3
                sum2 = sum2 + Rmax**2*rho_con    ! cm2/cm3
                conc = conc + rho_con            ! /cm3

                sum_q = sum_q + QTRC0(k,i,j,I_QI+iq-2)
             enddo

             if ( conc > 0.0_RP .and. sum_q > RLMTB ) then
                Re(k,i,j,I_HI+ic-1) = sum3/sum2                                             ! cm
             else
                Re(k,i,j,I_HI+ic-1) = 0.0_RP
             endif

          end do

       enddo
       enddo
       enddo

    else if (l_reff_version == 3) then
       ! just use circumscribing volume to calculate the maximum dimension without taking shape into account

       !$omp parallel do collapse(2) default(none) &
       !$omp private(i,j,k,ic,iq,iq1,iq2,ibin_is,ibin_sg,ibin_gh, &
       !$omp         sum2,sum3,conc,sum_q,rho_ice,rho_con,vcsw,Rmax) &
       !$omp shared(JS,JE,IS,IE,KS,KE,nbi, &
       !$omp        I_QI,iceConc_index,vol_index, &
       !$omp        numberPPVI, &
       !$omp        Re,QTRC0,DENS0)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE

          do ic = 1, 4

             if ( ic == 1 ) then
                iq1 = 1
                iq2 = ibin_is
             else if ( ic == 2 ) then
                iq1 = ibin_is + 1
                iq2 = ibin_sg
             else if ( ic == 3 ) then
                iq1 = ibin_sg + 1
                iq2 = ibin_gh
             else
                iq1 = ibin_gh + 1
                iq2 = nbi
             end if

             sum3 = 0.0_RP
             sum2 = 0.0_RP
             conc = 0.0_RP
             sum_q = 0.0_RP

             do iq = iq1, iq2

                rho_ice = QTRC0(k,i,j,I_QI+iq-2)                       * DENS0(k,i,j)  ! kg/m3
                rho_con = QTRC0(k,i,j,iceConc_index+(iq-1)*numberPPVI) * DENS0(k,i,j)  ! /cm3
                vcsw    = QTRC0(k,i,j,vol_index+(iq-1)*numberPPVI)     * DENS0(k,i,j)  ! cm3

                if ( rho_ice <= RILMTB .or. &
                     rho_con <= RILMTB) then
                   cycle
                endif

                Rmax = (vcsw/multiplicativeFactor)**oneThirdFactor                    ! cm

                sum3 = sum3 + Rmax**3*rho_con    ! cm3/cm3
                sum2 = sum2 + Rmax**2*rho_con    ! cm2/cm3
                conc = conc + rho_con            ! /cm3

                sum_q = sum_q + QTRC0(k,i,j,I_QI+iq-2)
             enddo

             if ( conc > 0.0_RP .and. sum_q > RLMTB ) then
                Re(k,i,j,I_HI+ic-1) = sum3/sum2                                             ! cm
             else
                Re(k,i,j,I_HI+ic-1) = 0.0_RP
             endif

          end do ! ic

       enddo
       enddo
       enddo

    endif

#endif
    return
  end subroutine
  !____________________________________________________________________________________



  !-----------------------------------------------------------------------------
  !> get mass ratio of each category
  subroutine ATMOS_PHY_MP_amps_qhyd2qtrc( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, &
       Qe,   &
       QTRC, &
       QNUM  )
#ifndef DISABLE_AMPS
    use com_amps, only: &
       nbin_h
#endif
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG, &
       I_HH
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)
    real(RP), intent(in) :: Qe(KA,IA,JA,N_HYD) ! mass ratio of each cateory [kg/kg]

    real(RP), intent(out)  :: QTRC(KA,IA,JA,QA-1)

    real(RP), intent(in), optional :: QNUM(KA,IA,JA,N_HYD) ! number concentration

#ifndef DISABLE_AMPS
    ! liquid bins
    real(RP),parameter :: c_mnm=4.188790205e-15, c_max=6.54498e-8, maxmass_r=5.2359870E-01
    real(RP) :: dsrat, dsrat_h
    real(RP) :: local_binbr(nbr+1)

    ! ice bins
    real(RP),parameter :: minmass_s = 4.18879020478639E-12
    real(RP),parameter :: mg1=4.18879020478639E-12, mg2=1.0E-6, mg3=1.0E-2, mg4=1.0E+1
    real(RP) :: srat_s
    real(RP) :: local_binbi(nbi+1)

    integer :: ibin_check=1
    real(RP), parameter :: den_aps=1.790_RP, den_api=2.60_RP, eps_ap=0.7_RP, r_t=0.3e-4_RP
    real(RP) :: molality, r_i, apt_ini, aps_ini, r_t2

    real(RP) :: dropletMass
    logical :: hdl_apt
    integer :: k, i, j, ibin, bin_check, lcon_index, lamt_index, lams_index

    real(RP) :: CM32M3

    real(RP) :: q, qn
    integer  :: ic


    CM32M3 = 1000000.0_RP

    ! liquid bins
    local_binbr(1) = c_mnm
    dsrat_h = (c_max/c_mnm)**(1.0/(nbin_h))
    do i = 2, nbin_h
       local_binbr(i) = local_binbr(i-1)*dsrat_h
    enddo
    local_binbr(nbin_h+1) = c_max
    dsrat=(maxmass_r/c_max)**(1.0/(nbr - nbin_h))

    local_binbr(nbr+1) = maxmass_r

    do i = nbin_h+2,nbr+1
       local_binbr(i) = local_binbr(i-1)*dsrat
    end do

    ! ice bins
    if (nbi == 40) then
       srat_s = (mg2/mg1)**(1.0_RP/8.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 8+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/20.0_RP)
       do i = 8+2, 8+20+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/12.0_RP)
       do i = 8+20+2, 8+20+12+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    elseif (nbi == 20) then
       srat_s = (mg2/mg1)**(1.0_RP/4.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 4+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/10.0_RP)
       do i = 4+2, 4+10+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/6.0_RP)
       do i = 4+10+2,4+10+6+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    elseif (nbi == 10) then
       srat_s = (mg2/mg1)**(1.0_RP/2.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 2+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/5.0_RP)
       do i = 2+2,2+5+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/3.0_RP)
       do i = 2+5+2, 2+5+3+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    else
       LOG_ERROR("ATMOS_PHY_MP_amps_qhyd2qtrc",*) 'no. of bins for ice is not 10, 20, or 40, please check.', nbi
       call PRC_abort
    endif

    lcon_index = I_QPPVL + 1
    lamt_index = I_QPPVL - 1
    lams_index = I_QPPVL

    !$omp workshare
    QTRC(:,:,:,:) = 0.0_RP
    !$omp end workshare

    if (present(QNUM)) then

       !$omp parallel do collapse(2) default(none) &
       !$omp private(i,j,k,ibin,bin_check,r_t2,dropletMass,hdl_apt,r_i,apt_ini, &
       !$omp         ic,q,qn) &
       !$omp shared(JS,JE,IS,IE,KS,KE,nbr, &
       !$omp        lcon_index,lamt_index,lams_index,numberPPVL, &
       !$omp        local_binbr,CM32M3, &
       !$omp        Qe,QNUM,QTRC,DENS)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          bin_check = 0
          do ic = 1, 2
             if ( ic==1 ) then
                q = Qe(k,i,j,I_HC)
             else
                q = Qe(k,i,j,I_HR)
             end if
             if ( Q == 0.0_RP .or. QN == 0.0_RP ) cycle
             bin_check = 1
             r_t2 = r_t
             dropletMass = Q * DENS(k,i,j) * 0.001_RP / QN * CM32M3
             do ibin = 1, nbr
                if ( dropletMass >= local_binbr(ibin) .and. dropletMass < local_binbr(ibin+1) ) then
                   QTRC(k,i,j,ibin) = QTRC(k,i,j,ibin) + Q
                   QTRC(k,i,j,lcon_index+numberPPVL*(ibin-1)) = QTRC(k,i,j,lcon_index+numberPPVL*(ibin-1)) + QN / CM32M3 / DENS(k,i,j)

                   hdl_apt = .true.
                   do while(hdl_apt)
                      r_i = r_t2*((den_aps * (1.0_RP - eps_ap) / eps_ap) &
                           / (den_api + (1.0_RP - eps_ap) / eps_ap * den_aps))**(1.0_RP/3.0_RP)

                      apt_ini = 4.0_RP * 3.14_RP / 3.0_RP &
                           * (r_i**3 * den_api + (r_t2**3 - r_i**3) * den_aps)

                      ! total ap mass
                      QTRC(k,i,j,lamt_index+numberPPVL*(ibin-1)) = &
                           apt_ini*QTRC(k,i,j,lcon_index+numberPPVL*(ibin-1))

                      ! soluble ap mass
                      QTRC(k,i,j,lams_index+numberPPVL*(ibin-1)) = &
                           QTRC(k,i,j,lamt_index+numberPPVL*(ibin-1))*eps_ap

                      if ((dropletMass + &
                           QTRC(k,i,j,lamt_index+numberPPVL*(ibin-1))) / &
                           QTRC(k,i,j,lcon_index+numberPPVL*(ibin-1)) > local_binbr(ibin+1)) then
                         QTRC(k,i,j,lamt_index+numberPPVL*(ibin-1)) = 0.0_RP
                         QTRC(k,i,j,lams_index+numberPPVL*(ibin-1)) = 0.0_RP
                         r_t2 = r_t2*0.9_RP
                      else
                         hdl_apt = .false.
                         bin_check = 1
                      endif
                   enddo

                   exit
                endif
             enddo ! ibin
          enddo ! iq

          !if ( bin_check == 0 ) then
          !   write(*,*) "ATMOS_PHY_MP_amps_qhyd2qtrc: The droplet mass is not within the bin region", Qe(k,i,j,I_HC), QNUM(k,i,j,I_HC)/CM32M3, local_binbr(1), Qe(k,i,j,I_HC)*DENS(k,i,j)*0.001_RP/QNUM(k,i,j,I_HC)*CM32M3, local_binbr(nbr+1), hdl_apt
          !   call PRC_abort
          !endif
       end do
       end do
       end do
    else
       !$omp parallel do collapse(2) default(none) &
       !$omp private(i,j,k) &
       !$omp shared(JS,JE,IS,IE,KS,KE, &
       !$omp        I_QPPVL, &
       !$omp        local_binbr, &
       !$omp        Qe,QTRC,DENS)
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          ! only first bin is filled
          QTRC(k,i,j,1) = Qe(k,i,j,I_HC) + Qe(k,i,j,I_HR)
          QTRC(k,i,j,I_QPPVL+1) = QTRC(k,i,j,1) / (1.05_RP*local_binbr(1)) / DENS(k,i,j)
       end do
       end do
       end do
    endif

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_qhyd2qtrc
  !____________________________________________________________________________________


  !-----------------------------------------------------------------------------
  !> Calculate mass ratio of each category
  subroutine ATMOS_PHY_MP_amps_qtrc2qhyd( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS0,      &
       QTRC0,      &
       Qe          )
    use scale_prc, only: &
       PRC_myrank
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: DENS0(KA,IA,JA)           ! moist air density [kg/m3]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA-1)      ! tracer mass concentration [kg/kg]
    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD)     ! mixing ratio of each cateory [kg/kg]

#ifndef DISABLE_AMPS
    integer :: ibin
    integer :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    !$omp workshare
    Qe(:,:,:,:) = 0.0_RP
    !$omp end workshare

    !$omp parallel do collapse(2) default(none) &
    !$omp private(i,j,k,ibin) &
    !$omp shared(JS,JE,IS,IE,KS,KE, &
    !$omp        split_bins,nbr,nbi,ibin_is,ibin_sg,ibin_gh,I_QI,I_QW, &
    !$omp        Qe,QTRC0)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
      ! --- liquid spectrum
      ! cloud droplets
      do ibin = 1, split_bins
         Qe(k,i,j,I_HC) = Qe(k,i,j,I_HC) + QTRC0(k,i,j,ibin)
      enddo
      do ibin = 1, ibin_is
         Qe(k,i,j,I_HC) = Qe(k,i,j,I_HC) + QTRC0(k,i,j,I_QW+ibin-2)
      end do
      ! rain
      do ibin = split_bins+1, nbr
         Qe(k,i,j,I_HR) = Qe(k,i,j,I_HR) + QTRC0(k,i,j,ibin)
      enddo
      do ibin = ibin_is+1, nbi
         Qe(k,i,j,I_HR) = Qe(k,i,j,I_HR) + QTRC0(k,i,j,I_QW+ibin-2)
      end do

      ! --- ice spectrum
      do ibin = 1, ibin_is
         Qe(k,i,j,I_HI) = Qe(k,i,j,I_HI) + QTRC0(k,i,j,I_QI+ibin-2)
      enddo
      do ibin = ibin_is+1, ibin_sg
         Qe(k,i,j,I_HS) = Qe(k,i,j,I_HS) + QTRC0(k,i,j,I_QI+ibin-2)
      enddo
      do ibin = ibin_sg+1, ibin_gh
         Qe(k,i,j,I_HI) = Qe(k,i,j,I_HI) + QTRC0(k,i,j,I_QI+ibin-2)
      enddo
      do ibin = ibin_gh+1, nbi
         Qe(k,i,j,I_HI) = Qe(k,i,j,I_HI) + QTRC0(k,i,j,I_QI+ibin-2)
      enddo

    enddo
    enddo
    enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_qtrc2qhyd
  !____________________________________________________________________________________


  !-----------------------------------------------------------------------------
  !> Calculate number concentration of each category
  subroutine ATMOS_PHY_MP_amps_qtrc2nhyd( &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE, &
       DENS0,      &
       QTRC0,      &
       Ne          )
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC,  &
       I_HR,  &
       I_HI,  &
       I_HS,  &
       I_HG,  &
       I_HH
    implicit none

    integer,  intent(in)  :: KA, KS, KE
    integer,  intent(in)  :: IA, IS, IE
    integer,  intent(in)  :: JA, JS, JE
    real(RP), intent(in)  :: DENS0(KA,IA,JA)           ! moist air density [kg/m3]
    real(RP), intent(in)  :: QTRC0(KA,IA,JA,QA-1)      ! tracer mass concentration [kg/kg]
    real(RP), intent(out) :: Ne   (KA,IA,JA,N_HYD)     ! number concentration of each cateory [1/m3]

#ifndef DISABLE_AMPS
    real(RP) :: CM32M3
    integer :: ibin, liqConc_index, iceConc_index
    integer :: k, i, j
    !---------------------------------------------------------------------------

!OCL XFILL
    Ne(:,:,:,:) = 0.0_RP

    liqConc_index = I_QPPVL + 2 - 1
    iceConc_index = I_QPPVI + 6 - 1

    CM32M3 = 1000000.0_RP

    !$omp parallel do collapse(2) default(none) &
    !$omp private(i,j,k,ibin) &
    !$omp shared(JS,JE,IS,IE,KS,KE, &
    !$omp        split_bins,nbr,nbi,liqConc_index,iceConc_index,numberPPVL,numberPPVI, &
    !$omp        CM32M3, &
    !$omp        Ne,QTRC0,DENS0)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       ! --- liquid spectrum
       ! cloud droplets
       do ibin = 1, split_bins
          Ne(k,i,j,I_HC) = Ne(k,i,j,I_HC) + QTRC0(k,i,j,liqConc_index+(ibin-1)*numberPPVL)*CM32M3 * DENS0(k,i,j)
       enddo
       ! rain
       do ibin = split_bins+1, nbr
          Ne(k,i,j,I_HR) = Ne(k,i,j,I_HR) + QTRC0(k,i,j,liqConc_index+(ibin-1)*numberPPVL)*CM32M3 * DENS0(k,i,j)
       enddo

       ! --- ice spectrum
       ! pristine ice, aggregate, rimed
       do ibin = 1, nbi
          Ne(k,i,j,I_HS) = Ne(k,i,j,I_HS) + QTRC0(k,i,j,iceConc_index+(ibin-1)*numberPPVI)*CM32M3 * DENS0(k,i,j)
       enddo
    enddo
    enddo
    enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_qtrc2nhyd
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !> get mass ratio of each category for 1D advection test
  subroutine ATMOS_PHY_MP_amps_init_qtrc_1DMODEL( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,    &
       QTRC,    &
       y_range, &
       case_number)
#ifndef DISABLE_AMPS
    use com_amps, only: &
       nbin_h
    use scale_atmos_grid_cartesC, only: &
       ATMOS_GRID_CARTESC_CYG
    use scale_prc_cartesC, only: &
       PRC_NUM_Y
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
#endif
    implicit none

    integer,  intent(in) :: KA, KS, KE
    integer,  intent(in) :: IA, IS, IE
    integer,  intent(in) :: JA, JS, JE
    real(RP), intent(in) :: y_range
    integer,  intent(in) :: case_number

    real(RP), intent(in) :: DENS(KA,IA,JA)

    real(RP), intent(out)  :: QTRC(KA,IA,JA,QA-1)

#ifndef DISABLE_AMPS
    real(RP) :: rj
    integer :: JMAX
    integer  :: i, j, k, ibin
    real(RP) :: ice_meanmass, ice_meancon, vol_axis, a_axis, c_axis, ag_axis, cg_axis, ex_axis
    ! ice bins
    real(RP),parameter :: minmass_s = 4.18879020478639E-12
    real(RP),parameter :: mg1=4.18879020478639E-12, mg2=1.0E-6, mg3=1.0E-2, mg4=1.0E+1
    real(RP) :: srat_s
    real(RP) :: local_binbi(nbi+1)


    ! define ice bins
    srat_s=(mg2/mg1)**(1.0_RP/4.0_RP)
    local_binbi(1) = minmass_s
    do i=2,4+1
       local_binbi(i)=srat_s*local_binbi(i-1)
    end do
    srat_s=(mg3/mg2)**(1.0_RP/10.0_RP)
    do i=4+2,4+10+1
       local_binbi(i)=srat_s*local_binbi(i-1)
    end do
    srat_s=(mg4/mg3)**(1.0_RP/6.0_RP)
    do i=4+10+2,4+10+6+1
       local_binbi(i)=srat_s*local_binbi(i-1)
    end do

    JMAX = JE - JS + 1

    !$omp parallel do &
    !$omp private (rj,ice_meanmass,ice_meancon,vol_axis,a_axis,c_axis,ag_axis,cg_axis,ex_axis)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       rj = ATMOS_GRID_CARTESC_CYG(JMAX*PRC_myrank+j) / y_range
       if ( rj <= 1.0_RP ) then

          do ibin = 1, 20
             ! a typical meanmass encountered in SHEBA case, can be changed to other values
             ice_meanmass = 2.0_RP*3.0_RP*sqrt(3.0_RP)*0.09126_RP*0.009_RP*0.0125_RP**2
             ice_meancon = 0.001_RP ! cm-3

             if (ice_meanmass < local_binbi(ibin) .or. ice_meanmass > local_binbi(ibin+1)) then
                cycle
             endif

             if (case_number == 1) then
                vol_axis = ( 2.0_RP*sqrt(0.0125_RP**2 + 0.002_RP**2) + rj*0.035_RP )
                a_axis = 0.0125_RP
                c_axis = 0.009_RP - rj*0.007_RP
                ag_axis  = rj * 0.00625_RP
                cg_axis  = rj * 0.001_RP
                ex_axis  = rj
             else if (case_number == 2) then
                vol_axis = ( 2.0_RP*sqrt(0.028125_RP**2 + 0.0045_RP**2) + rj * 0.02_RP )
                a_axis = 0.028125_RP - rj * 0.021875_RP
                c_axis = 0.0045_RP
                ag_axis  = rj * 0.003125_RP
                cg_axis  = rj * 0.00225_RP
                ex_axis  = rj
             endif
             !if (k == KS .and. i == IS) then
             !   write(*,'(I3, 8E15.6)') j, ATMOS_GRID_CARTESC_CYG(JMAX*PRC_myrank+j), rj, vol_axis, a_axis, c_axis, ag_axis, cg_axis, ex_axis
             !endif

             ! total mass
             QTRC(k,i,j,I_QI+ibin-2) = 0.0_RP

             ! rimed mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)-1) = 0.0_RP

             ! aggregate mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1) ) = 0.0_RP

             ! crystals  mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+1) = ice_meanmass*0.001_RP*ice_meancon*1000000.0_RP / DENS(k,i,j)

             ! frozen mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+2) = 0.0_RP

             ! total aerosol mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+3) = 0.0_RP

             ! soluble aerosol mass
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+4) = 0.0_RP

             ! number conc.
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+5) = ice_meancon / DENS(k,i,j)

             ! circumscribing volume, just a random number
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+6) = ice_meancon*PI/6.0_RP*(vol_axis)**3 / DENS(k,i,j)

             ! a axis
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+7) = ice_meancon*(a_axis)**3 / DENS(k,i,j)

             ! c axis
             if (l_gaxis_version == 2) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+8) = ice_meancon*(c_axis/a_axis)**3 / DENS(k,i,j)
             else
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+8) = ice_meancon*(c_axis)**3 / DENS(k,i,j)
             endif

             ! d axis
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+9) = 0.0_RP

             ! ag axis
             if (l_gaxis_version == 1 .or. l_gaxis_version == 2) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+10) = ice_meancon*(ag_axis)**3 / DENS(k,i,j)
             else if (l_gaxis_version == 3) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+10) = ice_meancon*(ag_axis/a_axis)**3 / DENS(k,i,j)
             endif

             ! cg axis
             if (l_gaxis_version == 1) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+11) = ice_meancon*(cg_axis)**3 / DENS(k,i,j)
             else if (l_gaxis_version == 2) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+11) = ice_meancon*(cg_axis/ag_axis)**3 / DENS(k,i,j)
             else if (l_gaxis_version == 3) then
                QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+11) = ice_meancon*(cg_axis/c_axis)**3 / DENS(k,i,j)
             endif

             ! ex
             QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+12) = ice_meancon*ex_axis / DENS(k,i,j)

             !if (k == KS .and. i == IS) then
             !   write(*,'(I3, 5E15.6)') j, ATMOS_GRID_CARTESC_CYG(JMAX*PRC_myrank+j), rj, QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+5), QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+11), QTRC(k,i,j,I_QPPVI+numberPPVI*(ibin-1)+12)
             !endif
          enddo
       endif
    enddo
    enddo
    enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_init_qtrc_1DMODEL
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !> get mass ratio of each category for box model
  subroutine ATMOS_PHY_MP_amps_init_qtrc_BOX( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, &
       QTRC, &
       box_initial_conc, &
       box_m_knot )
#ifndef DISABLE_AMPS
    use com_amps, only: &
       nbin_h
    use scale_prc, only: &
       PRC_abort, &
       PRC_myrank
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HR, &
       I_HI, &
       I_HS, &
       I_HG, &
       I_HH
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)

    real(RP), intent(out)  :: QTRC(KA,IA,JA,QA-1)

    real(RP), intent(in) :: box_initial_conc, box_m_knot

#ifndef DISABLE_AMPS
    ! liquid bins
    real(RP),parameter :: c_mnm=4.188790205e-15, c_max=6.54498e-8, maxmass_r=5.2359870E-01
    real(RP) :: dsrat, dsrat_h
    real(RP) :: local_binbr(nbr+1)

    ! ice bins
    real(RP),parameter :: minmass_s = 4.18879020478639E-12
    real(RP),parameter :: mg1=4.18879020478639E-12, mg2=1.0E-6, mg3=1.0E-2, mg4=1.0E+1
    real(RP) :: srat_s
    real(RP) :: local_binbi(nbi+1)

    integer :: ibin_check=1
    real(RP),parameter :: den_aps=1.790_RP, den_api=2.60_RP, eps_ap=0.7_RP, r_t=0.3e-4
    real(RP) :: molality, r_i, apt_ini, aps_ini, r_t2

    real(RP) :: dropletMass
    logical :: hdl_apt
    integer :: k, i, j, ibin, bin_check, lcon_index, lamt_index, lams_index


    ! liquid bins
    local_binbr(1) = c_mnm
    dsrat_h = (c_max/c_mnm)**(1.0/(nbin_h))
    do i = 2, nbin_h
       local_binbr(i) = local_binbr(i-1)*dsrat_h
    enddo
    local_binbr(nbin_h+1) = c_max
    dsrat=(maxmass_r/c_max)**(1.0/(nbr - nbin_h))

    local_binbr(nbr+1) = maxmass_r

    do i = nbin_h+2,nbr+1
       local_binbr(i) = local_binbr(i-1)*dsrat
    end do

    ! ice bins
    if (nbi == 40) then
       srat_s = (mg2/mg1)**(1.0_RP/8.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 8+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/20.0_RP)
       do i = 8+2, 8+20+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/12.0_RP)
       do i = 8+20+2, 8+20+12+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    elseif (nbi == 20) then
       srat_s = (mg2/mg1)**(1.0_RP/4.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 4+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/10.0_RP)
       do i = 4+2, 4+10+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/6.0_RP)
       do i = 4+10+2,4+10+6+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    elseif (nbi == 10) then
       srat_s = (mg2/mg1)**(1.0_RP/2.0_RP)
       local_binbi(1) = minmass_s
       do i = 2, 2+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg3/mg2)**(1.0_RP/5.0_RP)
       do i = 2+2,2+5+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
       srat_s = (mg4/mg3)**(1.0_RP/3.0_RP)
       do i = 2+5+2, 2+5+3+1
          local_binbi(i) = srat_s*local_binbi(i-1)
       enddo
    else
       LOG_ERROR("ATMOS_PHY_MP_amps_qhyd2qtrc",*) 'no. of bins for ice is not 10, 20, or 40, please check.', nbi
       call PRC_abort
    endif

    lcon_index = I_QPPVL + 1
    lamt_index = I_QPPVL - 1
    lams_index = I_QPPVL

    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       do ibin = 1, nbr
          ! m0 = 4.18 x 10-9 g
          ! f(X) = N / m0 exp(-2X), X = m/m0
          !m_knot = 8.3682D-9

          dropletMass = box_initial_conc * (exp(-local_binbr(ibin)/box_m_knot) - exp(-local_binbr(ibin+1)/box_m_knot))
          QTRC(k,i,j,lcon_index+numberPPVL*(ibin-1)) = dropletMass / DENS(k,i,j)

          dropletMass = box_initial_conc * (local_binbr(ibin)*exp(-local_binbr(ibin)/box_m_knot) - local_binbr(ibin+1)*exp(-local_binbr(ibin+1)/box_m_knot)) + &
                        box_initial_conc * box_m_knot * (exp(-local_binbr(ibin)/box_m_knot) - exp(-local_binbr(ibin+1)/box_m_knot))
          QTRC(k,i,j,ibin) = dropletMass * 1000.0_RP / DENS(k,i,j)
       enddo

       !write(*,*) k, i, j, sum(QTRC(k,i,j,1:nbr)), box_initial_conc
    enddo
    enddo
    enddo
    !do ibin = 1, nbr
    !   write(*,'(I3, 8ES15.6)') ibin, &
    !        QTRC(KS,IS,JS,lcon_index+numberPPVL*(ibin-1)) * DENS(KS,IS,JS), &
    !        QTRC(KS+1,IS,JS,lcon_index+numberPPVL*(ibin-1)) * DENS(KS+1,IS,JS), &
    !        QTRC(KS,IS+1,JS,lcon_index+numberPPVL*(ibin-1)) * DENS(KS,IS+1,JS), &
    !        QTRC(KS,IS,JS+1,lcon_index+numberPPVL*(ibin-1)) * DENS(KS,IS,JS+1), &
    !        QTRC(KS+1,IS+1,JS,lcon_index+numberPPVL*(ibin-1)) * DENS(KS+1,IS+1,JS), &
    !        QTRC(KS+1,IS,JS+1,lcon_index+numberPPVL*(ibin-1)) * DENS(KS+1,IS,JS+1), &
    !        QTRC(KS,IS+1,JS+1,lcon_index+numberPPVL*(ibin-1)) * DENS(KS,IS+1,JS+1), &
    !        QTRC(KS+1,IS+1,JS+1,lcon_index+numberPPVL*(ibin-1)) * DENS(KS+1,IS+1,JS+1)
    !enddo
    !write(*,*)
    !do ibin = 1, nbr
    !   write(*,'(I3, 8ES15.6)') ibin, &
    !        QTRC(KS,IS,JS,ibin) * DENS(KS,IS,JS), &
    !        QTRC(KS+1,IS,JS,ibin) * DENS(KS+1,IS,JS), &
    !        QTRC(KS,IS+1,JS,ibin) * DENS(KS,IS+1,JS), &
    !        QTRC(KS,IS,JS+1,ibin) * DENS(KS,IS,JS+1), &
    !        QTRC(KS+1,IS+1,JS,ibin) * DENS(KS+1,IS+1,JS), &
    !        QTRC(KS+1,IS,JS+1,ibin) * DENS(KS+1,IS,JS+1), &
    !        QTRC(KS,IS+1,JS+1,ibin) * DENS(KS,IS+1,JS+1), &
    !        QTRC(KS+1,IS+1,JS+1,ibin) * DENS(KS+1,IS+1,JS+1)
    !enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_amps_init_qtrc_BOX

end module scale_atmos_phy_mp_amps
