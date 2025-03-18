!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer

  use scale_prc, only: &
     PRC_abort
  use scale_const, only: &
     PI     => CONST_PI,     &
     GRAV   => CONST_GRAV,   &
     Pstd   => CONST_Pstd,   &
     Rdry   => CONST_Rdry,   &
     Rvap   => CONST_Rvap,   &
     CPdry  => CONST_CPdry,  &
     P00    => CONST_PRE00
  use scale_random, only: &
     RANDOM_uniform
  use scale_comm_cartesC, only: &
     COMM_vars8, &
     COMM_wait
  use scale_atmos_grid_cartesC, only: &
     CZ  => ATMOS_GRID_CARTESC_CZ,  &
     CX  => ATMOS_GRID_CARTESC_CX,  &
     CY  => ATMOS_GRID_CARTESC_CY,  &
     FZ  => ATMOS_GRID_CARTESC_FZ,  &
     FX  => ATMOS_GRID_CARTESC_FX,  &
     FY  => ATMOS_GRID_CARTESC_FY,  &
     CXG => ATMOS_GRID_CARTESC_CXG, &
     FXG => ATMOS_GRID_CARTESC_FXG, &
     FYG => ATMOS_GRID_CARTESC_FYG
  use scale_atmos_grid_cartesC_real, only: &
     REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
     REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
     AREA    => ATMOS_GRID_CARTESC_REAL_AREA
  use scale_atmos_profile, only: &
     PROFILE_isa => ATMOS_PROFILE_isa
  use scale_atmos_hydrometeor, only: &
     HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV
  use scale_atmos_hydrostatic, only: &
     HYDROSTATIC_buildrho        => ATMOS_HYDROSTATIC_buildrho,       &
     HYDROSTATIC_buildrho_atmos  => ATMOS_HYDROSTATIC_buildrho_atmos, &
     HYDROSTATIC_buildrho_bytemp => ATMOS_HYDROSTATIC_buildrho_bytemp
  use scale_atmos_saturation, only: &
     SATURATION_pres2qsat_all => ATMOS_SATURATION_pres2qsat_all, &
     SATURATION_psat_all => ATMOS_SATURATION_psat_all
  use mod_atmos_vars, only: &
     DENS, &
     MOMX, &
     MOMY, &
     MOMZ, &
     RHOT, &
     QTRC
  use mod_atmos_phy_ae_vars, only: &
     CCN => ATMOS_PHY_AE_CCN

  use mod_atmos_admin, only: &
     ATMOS_sw_dyn_DGM
  use mod_atmos_dyn_vars, only: &
     dyn_dg_vars
  use scale_atmos_dyn_dgm_operator, only: &
     fem_mesh3D => mesh3D, &
     fem_elem3D => elem3D
  use scale_fem_localmesh_3d, only: &
     FEM_LocalMesh3D => LocalMesh3D
  use scale_fem_meshfield_base, only: &
     FEM_MeshField3D => MeshField3D

  use mod_mkinit_common, only: &
     MKINIT_common_BUBBLE_setup, &
     MKINIT_common_BUBBLE_setup_DG, &
     MKINIT_common_RECT_setup,  &
     MKINIT_common_flux_setup,  &
     MKINIT_common_land_setup,  &
     MKINIT_common_ocean_setup, &
     MKINIT_common_urban_setup, &
     MKINIT_common_read_sounding,       &
     MKINIT_common_read_sounding_mpace, &
     THETAstd
   
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKINIT_setup
  public :: MKINIT_finalize
  public :: MKINIT

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public            :: MKINIT_TYPE        = -1
  integer, public, parameter :: I_IGNORE           =  0

  integer, public, parameter :: I_PLANESTATE       =  1
  integer, public, parameter :: I_TRACERBUBBLE     =  2
  integer, public, parameter :: I_COLDBUBBLE       =  3

  integer, public, parameter :: I_LAMBWAVE         =  4
  integer, public, parameter :: I_GRAVITYWAVE      =  5
  integer, public, parameter :: I_KHWAVE           =  6
  integer, public, parameter :: I_TURBULENCE       =  7
  integer, public, parameter :: I_MOUNTAINWAVE     =  8

  integer, public, parameter :: I_WARMBUBBLE       =  9
  integer, public, parameter :: I_SUPERCELL        = 10
  integer, public, parameter :: I_SQUALLLINE       = 11
  integer, public, parameter :: I_WK1982           = 12
  integer, public, parameter :: I_DYCOMS2_RF01     = 13
  integer, public, parameter :: I_DYCOMS2_RF02     = 14
  integer, public, parameter :: I_RICO             = 15

  integer, public, parameter :: I_INTERPORATION    = 16

  integer, public, parameter :: I_LANDCOUPLE       = 17
  integer, public, parameter :: I_OCEANCOUPLE      = 18
  integer, public, parameter :: I_URBANCOUPLE      = 19
  integer, public, parameter :: I_TRIPLECOUPLE     = 20
  integer, public, parameter :: I_BUBBLECOUPLE     = 21

  integer, public, parameter :: I_SEABREEZE        = 22
  integer, public, parameter :: I_HEATISLAND       = 23

  integer, public, parameter :: I_DYCOMS2_RF02_DNS = 24

  integer, public, parameter :: I_REAL             = 25

  integer, public, parameter :: I_GRAYZONE         = 26
  integer, public, parameter :: I_BOXAERO          = 27
  integer, public, parameter :: I_WARMBUBBLEAERO   = 28

  integer, public, parameter :: I_CAVITYFLOW       = 29
  integer, public, parameter :: I_BAROCWAVE        = 30
  integer, public, parameter :: I_BOMEX            = 31

  integer, public, parameter :: I_MPACE            = 32
!  integer, public, parameter :: I_BOXAMPS          = 33

  integer, public, parameter :: I_SONDE_PERTURB    = 34

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: SBMAERO_setup
  private :: AEROSOL_setup

  private :: MKINIT_planestate
  private :: MKINIT_tracerbubble
  private :: MKINIT_coldbubble
  private :: MKINIT_lambwave
  private :: MKINIT_gravitywave
  private :: MKINIT_khwave
  private :: MKINIT_turbulence
  private :: MKINIT_cavityflow
  private :: MKINIT_mountainwave
  private :: MKINIT_barocwave

  private :: MKINIT_warmbubble
  private :: MKINIT_supercell
  private :: MKINIT_squallline
  private :: MKINIT_wk1982
  private :: MKINIT_DYCOMS2_RF01
  private :: MKINIT_DYCOMS2_RF02
  private :: MKINIT_RICO
  private :: MKINIT_BOMEX

  private :: MKINIT_landcouple
  private :: MKINIT_oceancouple
  private :: MKINIT_urbancouple
  private :: MKINIT_seabreeze
  private :: MKINIT_heatisland

  private :: MKINIT_DYCOMS2_RF02_DNS

  private :: MKINIT_real

  private :: MKINIT_grayzone

  private :: MKINIT_boxaero
  private :: MKINIT_warmbubbleaero

  private :: MKINIT_MPACE
!  private :: MKINIT_boxamps

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter           :: NITER_RH = 4

  real(RP), private, allocatable         :: pres    (:,:,:) ! pressure [Pa]
  real(RP), private, allocatable         :: temp    (:,:,:) ! temperature [K]
  real(RP), private, allocatable         :: pott    (:,:,:) ! potential temperature [K]
  real(RP), private, allocatable         :: psat    (:,:,:) ! satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv      (:,:,:) ! water vapor [kg/kg]
  real(RP), private, allocatable         :: qc      (:,:,:) ! cloud water [kg/kg]
  real(RP), private, allocatable         :: nc      (:,:,:) ! cloud water number density [1/kg]
  real(RP), private, allocatable         :: velx    (:,:,:) ! velocity u [m/s]
  real(RP), private, allocatable         :: vely    (:,:,:) ! velocity v [m/s]
  real(RP), private, allocatable         :: ptrc    (:,:,:) ! passive tracer

  real(RP), private, allocatable         :: pres_sfc(:,:) ! surface pressure [Pa]
  real(RP), private, allocatable         :: temp_sfc(:,:) ! surface temperature [K]
  real(RP), private, allocatable         :: pott_sfc(:,:) ! surface potential temperature [K]
  real(RP), private, allocatable         :: qsat_sfc(:,:) ! surface satulated water vapor [kg/kg]
  real(RP), private, allocatable         :: qv_sfc  (:,:) ! surface water vapor [kg/kg]
  real(RP), private, allocatable         :: qc_sfc  (:,:) ! surface cloud water [kg/kg]

  real(RP), private, allocatable         :: rndm    (:,:,:) ! random    number (0-1)
  real(RP), private, allocatable, target :: bubble  (:,:,:) ! bubble    factor (0-1)
  real(RP), private, allocatable, target :: rect    (:,:,:) ! rectangle factor (0-1)
  real(RP), private, allocatable         :: gan     (:)     ! gamma     factor (0-1)

  type(FEM_MeshField3D), pointer :: DG_MOMX, DG_MOMY, DG_MOMZ, DG_DDENS, DG_DRHOT
  type(FEM_MeshField3D), pointer :: DG_DENS_hyd, DG_PRES_hyd
  type(FEM_MeshField3D) :: DG_QV, DG_QC, DG_NC
  type(FEM_MeshField3D) :: DG_ptrc
  type(FEM_MeshField3D), target :: DG_bubble

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine MKINIT_setup
    implicit none

    character(len=H_SHORT) :: MKINIT_initname = 'NONE'

    namelist / PARAM_MKINIT / &
       MKINIT_initname

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_setup",*) 'Not appropriate names in namelist PARAM_MKINIT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT)

    allocate( pres(KA,IA,JA) )
    allocate( temp(KA,IA,JA) )
    allocate( pott(KA,IA,JA) )
    allocate( psat(KA,IA,JA) )
    allocate( qv  (KA,IA,JA) )
    allocate( qc  (KA,IA,JA) )
    allocate( nc  (KA,IA,JA) )
    allocate( velx(KA,IA,JA) )
    allocate( vely(KA,IA,JA) )
    allocate( ptrc(KA,IA,JA) )

    allocate( pres_sfc(IA,JA) )
    allocate( temp_sfc(IA,JA) )
    allocate( pott_sfc(IA,JA) )
    allocate( qsat_sfc(IA,JA) )
    allocate( qv_sfc  (IA,JA) )
    allocate( qc_sfc  (IA,JA) )

    allocate( rndm  (KA,IA,JA) )
    allocate( bubble(KA,IA,JA) )
    allocate( rect  (KA,IA,JA) )

    if ( ATMOS_sw_dyn_DGM ) then
      call DG_ptrc%Init( "ptrc", "", fem_mesh3d )
      call DG_bubble%Init( "bubble", "", fem_mesh3D )
      call DG_QV%Init( "QV", "", fem_mesh3d, fill_value=0.0_RP )
      call DG_QC%Init( "QC", "", fem_mesh3d, fill_value=0.0_RP )
      call DG_NC%Init( "NC", "", fem_mesh3d, fill_value=0.0_RP )
    end if

    !$acc enter data create(pres,temp,pott,psat,qv,qc,nc,velx,vely,ptrc,pres_sfc,temp_sfc,pott_sfc,qsat_sfc,qv_sfc,qc_sfc,rndm,bubble,rect)

    select case(trim(MKINIT_initname))
    case('NONE')
       MKINIT_TYPE = I_IGNORE
    case('PLANESTATE')
       MKINIT_TYPE = I_PLANESTATE
    case('TRACERBUBBLE')
       MKINIT_TYPE = I_TRACERBUBBLE
    case('COLDBUBBLE')
       MKINIT_TYPE = I_COLDBUBBLE
       call MKINIT_common_BUBBLE_setup( bubble )
    case('LAMBWAVE')
       MKINIT_TYPE = I_LAMBWAVE
       call MKINIT_common_BUBBLE_setup( bubble )
    case('GRAVITYWAVE')
       MKINIT_TYPE = I_GRAVITYWAVE
       call MKINIT_common_BUBBLE_setup( bubble )
    case('KHWAVE')
       MKINIT_TYPE = I_KHWAVE
    case('TURBULENCE')
       MKINIT_TYPE = I_TURBULENCE
    case('MOUNTAINWAVE')
       MKINIT_TYPE = I_MOUNTAINWAVE
       call MKINIT_common_BUBBLE_setup( bubble )
       if ( ATMOS_sw_dyn_DGM ) call MKINIT_common_BUBBLE_setup_DG( DG_bubble )
    case('WARMBUBBLE')
       MKINIT_TYPE = I_WARMBUBBLE
       call MKINIT_common_BUBBLE_setup( bubble )
       if ( ATMOS_sw_dyn_DGM ) call MKINIT_common_BUBBLE_setup_DG( DG_bubble )
    case('SUPERCELL')
       MKINIT_TYPE = I_SUPERCELL
       call MKINIT_common_BUBBLE_setup( bubble )
    case('SQUALLLINE')
       MKINIT_TYPE = I_SQUALLLINE
    case('WK1982')
       MKINIT_TYPE = I_WK1982
       call MKINIT_common_BUBBLE_setup( bubble )
    case('DYCOMS2_RF01')
       MKINIT_TYPE = I_DYCOMS2_RF01
    case('DYCOMS2_RF02')
       MKINIT_TYPE = I_DYCOMS2_RF02
    case('RICO')
       MKINIT_TYPE = I_RICO
    case('BOMEX')
       MKINIT_TYPE = I_BOMEX
    case('INTERPORATION')
       MKINIT_TYPE = I_INTERPORATION
    case('LANDCOUPLE')
       MKINIT_TYPE = I_LANDCOUPLE
    case('OCEANCOUPLE')
       MKINIT_TYPE = I_OCEANCOUPLE
    case('URBANCOUPLE')
       MKINIT_TYPE = I_URBANCOUPLE
    case('TRIPLECOUPLE')
       MKINIT_TYPE = I_TRIPLECOUPLE
    case('BUBBLECOUPLE')
       MKINIT_TYPE = I_BUBBLECOUPLE
       call MKINIT_common_BUBBLE_setup( bubble )
    case('SEABREEZE')
       MKINIT_TYPE = I_SEABREEZE
    case('HEATISLAND')
       MKINIT_TYPE = I_HEATISLAND
    case('DYCOMS2_RF02_DNS')
       MKINIT_TYPE = I_DYCOMS2_RF02_DNS
    case('REAL')
       MKINIT_TYPE = I_REAL
    case('GRAYZONE')
       MKINIT_TYPE = I_GRAYZONE
    case('BOXAERO')
       MKINIT_TYPE = I_BOXAERO
    case('WARMBUBBLEAERO')
       MKINIT_TYPE = I_WARMBUBBLEAERO
       call MKINIT_common_BUBBLE_setup( bubble )
    case('CAVITYFLOW')
       MKINIT_TYPE = I_CAVITYFLOW
    case('BAROCWAVE')
       MKINIT_TYPE = I_BAROCWAVE
    case('MPACE')
       MKINIT_TYPE = I_MPACE
!    case('BOXAMPS')
!       MKINIT_TYPE = I_BOXAMPS
    case('SONDE_PERTURB')
       MKINIT_TYPE = I_SONDE_PERTURB
       call MKINIT_common_BUBBLE_setup( bubble )
    case default
       LOG_ERROR("MKINIT_setup",*) 'Unsupported TYPE:', trim(MKINIT_initname)
       call PRC_abort
    endselect

    return
  end subroutine MKINIT_setup

  !-----------------------------------------------------------------------------
  !> Finalize
  subroutine MKINIT_finalize
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_finalize",*) 'Finalize'

    !$acc exit data delete(pres,temp,pott,psat,qv,qc,nc,velx,vely,ptrc,pres_sfc,temp_sfc,pott_sfc,qsat_sfc,qv_sfc,qc_sfc,rndm,bubble,rect)

    deallocate( pres )
    deallocate( temp )
    deallocate( pott )
    deallocate( psat )
    deallocate( qv   )
    deallocate( qc   )
    deallocate( nc   )
    deallocate( velx )
    deallocate( vely )
    deallocate( ptrc )

    deallocate( pres_sfc )
    deallocate( temp_sfc )
    deallocate( pott_sfc )
    deallocate( qsat_sfc )
    deallocate( qv_sfc   )
    deallocate( qc_sfc   )

    deallocate( rndm   )
    deallocate( bubble )
    deallocate( rect   )

    if ( ATMOS_sw_dyn_DGM ) then
      call DG_ptrc%Final()
      call DG_bubble%Final()
      call DG_QV%Final(); call DG_QC%Final(); call DG_NC%Final()
   endif

    return
  end subroutine MKINIT_finalize

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine MKINIT( output )
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       N_HYD, &
       I_HC
    use mod_atmos_phy_mp_vars, only: &
       QS_MP, &
       QE_MP
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use mod_atmos_phy_mp_driver, only: &
       ATMOS_PHY_MP_driver_qhyd2qtrc
!    use scale_atmos_phy_mp_sdm, only: &
!       ATMOS_PHY_MP_sdm_random_setup
    implicit none
    logical, intent(out) :: output

    real(RP) :: QHYD(KA,IA,JA,N_HYD)
    real(RP) :: QNUM(KA,IA,JA,N_HYD)

    logical :: convert_qtrc
    integer :: k, i, j, iq

    integer :: ldomID, kelem
    type(FEM_LocalMesh3D), pointer :: lmesh
    real(RP), allocatable :: DG_DENS(:,:)
    real(RP), allocatable :: DG_QHYD(:,:,:)
    real(RP), allocatable :: DG_QNUM(:,:,:)
    real(RP), allocatable :: DG_QTRC_out(:,:,:)
    !---------------------------------------------------------------------------

    if ( MKINIT_TYPE == I_IGNORE ) then
      LOG_NEWLINE
      LOG_PROGRESS(*) 'skip  making initial data'
      output = .false.
    else
      LOG_NEWLINE
      LOG_PROGRESS(*) 'start making initial data'

      !--- Initialize variables
      !$omp workshare
      !$acc kernels
      pres(:,:,:) = UNDEF
      temp(:,:,:) = UNDEF
      pott(:,:,:) = UNDEF
      psat(:,:,:) = UNDEF
      velx(:,:,:) = UNDEF
      vely(:,:,:) = UNDEF

      rndm  (:,:,:) = UNDEF
      !$acc end kernels

      !$acc kernels
      pres_sfc(:,:) = UNDEF
      temp_sfc(:,:) = UNDEF
      pott_sfc(:,:) = UNDEF
      qsat_sfc(:,:) = UNDEF
      !$acc end kernels

      !$acc kernels
      qv    (:,:,:) = 0.0_RP
      qc    (:,:,:) = 0.0_RP
      nc    (:,:,:) = 0.0_RP
      !$acc end kernels
      !$acc kernels
      qv_sfc(:,:) = 0.0_RP
      qc_sfc(:,:) = 0.0_RP
      !$acc end kernels

      !$acc kernels
      ptrc(:,:,:) = UNDEF
      !$acc end kernels

      !$acc kernels
!OCL XFILL
      QTRC(:,:,:,:) = UNDEF
!OCL XFILL
      QHYD(:,:,:,:) = 0.0_RP
!OCL XFILL
      QNUM(:,:,:,:) = 0.0_RP
      !$acc end kernels
      !$omp end workshare

      if ( ATMOS_sw_dyn_DGM ) then
         call dyn_dg_vars%GetVarsPtr( DG_DDENS, DG_MOMX, DG_MOMY, DG_MOMZ, DG_DRHOT, &
            DG_DENS_hyd, DG_PRES_hyd )
      end if

      call PROF_rapstart('_MkInit_main',3)

      convert_qtrc = .true.

!      if( ATMOS_PHY_MP_TYPE .eq. 'SDM' ) then
!         call ATMOS_PHY_MP_sdm_random_setup
!      end if

      select case(MKINIT_TYPE)
      case(I_PLANESTATE)
         call MKINIT_planestate
      case(I_TRACERBUBBLE)
         call MKINIT_tracerbubble
      case(I_COLDBUBBLE)
         call MKINIT_coldbubble
      case(I_LAMBWAVE)
         call MKINIT_lambwave
      case(I_GRAVITYWAVE)
         call MKINIT_gravitywave
      case(I_KHWAVE)
         call MKINIT_khwave
      case(I_TURBULENCE)
         call MKINIT_turbulence
      case(I_MOUNTAINWAVE)
         call MKINIT_mountainwave
      case(I_WARMBUBBLE)
         call MKINIT_warmbubble
      case(I_SUPERCELL)
         call MKINIT_supercell
      case(I_SQUALLLINE)
         call MKINIT_squallline
      case(I_WK1982)
         call MKINIT_wk1982
      case(I_DYCOMS2_RF01)
         call MKINIT_DYCOMS2_RF01
      case(I_DYCOMS2_RF02)
         call MKINIT_DYCOMS2_RF02
      case(I_RICO)
         call MKINIT_RICO
      case(I_BOMEX)
         call MKINIT_BOMEX
      case(I_OCEANCOUPLE)
         call MKINIT_planestate
         call MKINIT_oceancouple
      case(I_LANDCOUPLE)
         call MKINIT_planestate
         call MKINIT_landcouple
      case(I_URBANCOUPLE)
         call MKINIT_planestate
         call MKINIT_urbancouple
      case(I_TRIPLECOUPLE)
         call MKINIT_planestate
         call MKINIT_oceancouple
         call MKINIT_landcouple
         call MKINIT_urbancouple
      case(I_BUBBLECOUPLE)
         call MKINIT_planestate
         call MKINIT_warmbubble
         call MKINIT_oceancouple
         call MKINIT_landcouple
         call MKINIT_urbancouple
      case(I_SEABREEZE)
         call MKINIT_planestate
         call MKINIT_seabreeze
      case(I_HEATISLAND)
         call MKINIT_planestate
         call MKINIT_heatisland
      case(I_DYCOMS2_RF02_DNS)
         call MKINIT_DYCOMS2_RF02_DNS
      case(I_REAL)
         call MKINIT_real
         convert_qtrc = .false.
      case(I_GRAYZONE)
         call MKINIT_grayzone
      case(I_BOXAERO)
         call MKINIT_boxaero
      case(I_WARMBUBBLEAERO)
         call MKINIT_warmbubbleaero
      case(I_CAVITYFLOW)
         call MKINIT_cavityflow
      case(I_BAROCWAVE)
         call MKINIT_barocwave
      case(I_MPACE)
         call MKINIT_MPACE
!      case(I_BOXAMPS)
!         call MKINIT_boxamps
         convert_qtrc = .false.
      case(I_SONDE_PERTURB)
         call MKINIT_sonde_perturb
      case default
         LOG_ERROR("MKINIT",*) 'Unsupported TYPE:', MKINIT_TYPE
         call PRC_abort
      endselect

      ! water content
      if ( ( .not. ATMOS_HYDROMETEOR_dry ) .AND. convert_qtrc ) then
         !$acc kernels
!OCL XFILL
         QHYD(:,:,:,I_HC) = qc(:,:,:)
!OCL XFILL
         QNUM(:,:,:,I_HC) = nc(:,:,:)
         !$acc end kernels
         call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                             qv(:,:,:), QHYD(:,:,:,:), & ! [IN]
                                             DENS(:,:,:),              & ! [IN]
                                             QTRC(:,:,:,QS_MP:QE_MP),  & ! [OUT]
                                             QNUM=QNUM(:,:,:,:)        ) ! [IN]

         if ( ATMOS_sw_dyn_DGM ) then
           do ldomID=1, fem_mesh3d%LOCAL_MESH_NUM
             lmesh => fem_mesh3d%lcmesh_list(ldomID)
             allocate( DG_DENS(fem_elem3D%Np,lmesh%NeA), DG_QHYD(fem_elem3d%Np,lmesh%NeA,N_HYD), DG_QNUM(fem_elem3d%Np,lmesh%NeA,N_HYD) )
             allocate( DG_QTRC_out(fem_elem3d%Np,lmesh%NeA,QS_MP:QE_MP) )

             !$omp parallel 
             !$omp workshare
             DG_QHYD(:,:,:) = 0.0_RP
             DG_QNUM(:,:,:) = 0.0_RP
             !$omp end workshare
             !$omp do
             do kelem=lmesh%NeS, lmesh%NeE
               DG_DENS(:,kelem) = DG_DENS_hyd%local(ldomID)%val(:,kelem) + DG_DDENS%local(ldomID)%val(:,kelem)
               DG_QHYD(:,kelem,I_HC) = DG_QC%local(ldomID)%val(:,kelem)
               DG_QNUM(:,kelem,I_HC) = DG_NC%local(ldomID)%val(:,kelem)
             enddo
             !$omp end do
             !$omp end parallel
             call ATMOS_PHY_MP_driver_qhyd2qtrc( &
               fem_elem3D%Np, 1, fem_elem3D%Np, lmesh%NeA, lmesh%NeS, lmesh%NeE, 1, 1, 1, &
               DG_qv%local(ldomID)%val(:,:), DG_QHYD(:,:,:), DG_DENS(:,:),                & ! [IN]
               DG_QTRC_out(:,:,QS_MP:QE_MP),  & ! [OUT]
               QNUM=DG_QNUM(:,:,:)            ) ! [IN]
               
             do iq=QS_MP, QE_MP
               !$omp parallel do
               do kelem=lmesh%NeS, lmesh%NeE
                  dyn_dg_vars%qtrc_vars(iq)%local(ldomID)%val(:,kelem) = DG_QTRC_out(:,kelem,iq)
               enddo
             enddo

             deallocate( DG_DENS, DG_QHYD, DG_QNUM )
             deallocate( DG_QTRC_out )
           enddo
         end if
      end if

      call tke_setup

      call AEROSOL_setup

      call SBMAERO_setup( convert_qtrc ) ! [INOUT]

      ! passive tracer
      call TRACER_inq_id( "PTracer", iq )
      if ( iq > 0 ) then
         QTRC(:,:,:,iq) = ptrc(:,:,:)
         if ( ATMOS_sw_dyn_DGM ) then
            do ldomID=1, fem_mesh3d%LOCAL_MESH_NUM
               dyn_dg_vars%qtrc_vars(iq)%local(ldomID)%val(:,:) = DG_ptrc%local(ldomID)%val(:,:)
            enddo
         end if
      endif

      !$omp parallel do collapse(3)
      !$acc kernels
      do iq = 1, QA
      do j = 1, JA
      do i = 1, iA
      do k = 1, KA
         if ( QTRC(k,i,j,iq) == UNDEF ) then
            QTRC(k,i,j,iq) = 0.0_RP
         end if
      end do
      end do
      end do
      end do
      !$acc end kernels

      call PROF_rapend  ('_MkInit_main',3)

      LOG_PROGRESS(*) 'end   making initial data'

      output = .true.

    endif

    return
  end subroutine MKINIT

  !-----------------------------------------------------------------------------
  !> Setup aerosol condition for Kajino 2013 or SPRINTARS scheme
  subroutine AEROSOL_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_mkinit
    use scale_atmos_phy_ae_offline, only: &
       ATMOS_PHY_AE_offline_mkinit
    use scale_atmos_phy_ae_sprintars, only: &
       ATMOS_PHY_AE_SPRINTARS_mkinit, &
       ATMOS_PHY_AE_SPRINTARS_mkemission
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    implicit none

    real(RP), parameter :: d_min_def = 1.e-9_RP ! default lower bound of 1st size bin
    real(RP), parameter :: d_max_def = 1.e-5_RP ! upper bound of last size bin
    integer,  parameter :: n_kap_def = 1        ! number of kappa bins
    real(RP), parameter :: k_min_def = 0.e0_RP  ! lower bound of 1st kappa bin
    real(RP), parameter :: k_max_def = 1.e0_RP  ! upper bound of last kappa bin

    real(RP) :: ccn_init = 50.E+6_RP ! initial cloud condensation nucrei [#/m3]

    real(RP) :: m0_init = 0.0_RP    ! initial total num. conc. of modes (Atk,Acm,Cor) [#/m3]
    real(RP) :: dg_init = 80.e-9_RP ! initial number equivalen diameters of modes     [m]
    real(RP) :: sg_init = 1.6_RP    ! initial standard deviation                      [-]

    real(RP) :: d_min_inp(3) = d_min_def
    real(RP) :: d_max_inp(3) = d_max_def
    real(RP) :: k_min_inp(3) = k_min_def
    real(RP) :: k_max_inp(3) = k_max_def
    integer  :: n_kap_inp(3) = n_kap_def

    real(RP) :: qdry(KA,IA,JA)

    namelist / PARAM_AERO / &
       ccn_init,  &
       m0_init,   &
       dg_init,   &
       sg_init,   &
       d_min_inp, &
       d_max_inp, &
       k_min_inp, &
       k_max_inp, &
       n_kap_inp

    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_AE_TYPE /= 'OFF' .AND. ATMOS_PHY_AE_TYPE /= 'NONE' ) then

       LOG_NEWLINE
       LOG_INFO("AEROSOL_setup",*) 'Setup'

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_AERO,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          LOG_INFO("AEROSOL_setup",*) 'Not found namelist. Default used!'
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("AEROSOL_setup",*) 'Not appropriate names in namelist PARAM_AERO. Check!'
          call PRC_abort
       endif
       LOG_NML(PARAM_AERO)

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          !$acc data create(qdry)
          !$acc kernels
          qdry(:,:,:) = 1.0_RP - qv(:,:,:) - qc(:,:,:)
          !$acc end kernels
          !$acc update host(dens,temp,pres,qdry,qv)
          call ATMOS_PHY_AE_kajino13_mkinit( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
                                             QA_AE,                   & ! (in)
                                             DENS(:,:,:),             & ! (in)
                                             TEMP(:,:,:),             & ! (in)
                                             PRES(:,:,:),             & ! (in)
                                             QDRY(:,:,:),             & ! (in)
                                             QV  (:,:,:),             & ! (in)
                                             m0_init,                 & ! (in)
                                             dg_init,                 & ! (in)
                                             sg_init,                 & ! (in)
                                             d_min_inp(:),            & ! (in)
                                             d_max_inp(:),            & ! (in)
                                             k_min_inp(:),            & ! (in)
                                             k_max_inp(:),            & ! (in)
                                             n_kap_inp(:),            & ! (in)
                                             QTRC(:,:,:,QS_AE:QE_AE), & ! (out)
                                             CCN(:,:,:)               ) ! (out)
          !$acc update device(QTRC(:,:,:,QS_AE:QE_AE),CCN)
          !$acc end data
       case ( 'OFFLINE' )
          call ATMOS_PHY_AE_offline_mkinit ( KA, KS, KE, IA, IS, IE, JA, JS, JE, & ! (in)
                                             ccn_init,                & ! (in)
                                             CCN(:,:,:)               ) ! (out)
       case ( 'SPRINTARS' )
          call ATMOS_PHY_AE_SPRINTARS_mkinit( QA_AE,                   & ! (in)
                                              ccn_init,                & ! (in)
                                              QTRC(:,:,:,QS_AE:QE_AE), & ! (inout)
                                              CCN (:,:,:)              ) ! (inout)
          !$acc update device(QTRC(:,:,:,QS_AE:QE_AE),CCN)
          call ATMOS_PHY_AE_SPRINTARS_mkemission
       case default
          !$acc kernels
          CCN(:,:,:) = ccn_init
          !$acc end kernels
       end select

    endif

    return
  end subroutine AEROSOL_setup

  !-----------------------------------------------------------------------------
  !> Setup aerosol condition for Spectral Bin Microphysics (SBM) model
  subroutine SBMAERO_setup( convert_qtrc )
    use scale_atmos_hydrometeor, only: &
       I_QV, &
       QHE
    use mod_atmos_admin, only: &
       ATMOS_PHY_MP_TYPE
    use scale_atmos_phy_mp_suzuki10, only: &
       nccn
    implicit none

    logical, intent(inout) :: convert_qtrc

    integer :: iq, i, j, k
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_MP_TYPE /= 'SUZUKI10' ) return

    if ( .not. convert_qtrc ) return

    !--- Super saturated air at initial
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k  = KS, KE
       QTRC(k,i,j,I_QV) = qv(k,i,j) + qc(k,i,j)
    end do
    end do
    end do
    !$acc end kernels

    !-- Aerosol distribution
    if ( nccn /= 0 ) then
       !$acc kernels
       !$acc loop collapse(4) independent
       do iq = 1, nccn
       do j = JSB, JEB
       do i = ISB, IEB
       do k  = KS, KE
          QTRC(k,i,j,QHE+iq) = gan(iq) / DENS(k,i,j) ! [note] gan is never set.
       enddo
       enddo
       enddo
       enddo
       !$acc end kernels
    endif

    convert_qtrc = .false.

    return
  end subroutine SBMAERO_setup

  !-----------------------------------------------------------------------------
  function faero( f0,r0,x,alpha,rhoa )
    use scale_const, only: &
       pi => CONST_PI
    implicit none

    real(RP), intent(in) ::  x, f0, r0, alpha, rhoa
    real(RP) :: faero
    real(RP) :: rad
    !---------------------------------------------------------------------------

    rad = ( exp(x) * 3.0_RP / 4.0_RP / pi / rhoa )**(1.0_RP/3.0_RP)

    faero = f0 * (rad/r0)**(-alpha)

    return
  end function faero

  !-----------------------------------------------------------------------------
  !> TKE setup
  subroutine tke_setup
    use scale_const, only: &
       EPS   => CONST_EPS, &
       UNDEF => CONST_UNDEF
    use mod_atmos_phy_tb_vars, only: &
       I_TKE
    use mod_atmos_phy_bl_vars, only: &
       Zi => ATMOS_PHY_BL_Zi, &
       QS_BL => QS, &
       QE_BL => QE
    use mod_atmos_phy_bl_driver, only: &
       atmos_phy_bl_driver_mkinit
    use scale_fem_meshfield_util, only: &
       MeshFieldUtil_interp_FVtoDG_3D
    implicit none

    real(RP) :: TKE_CONST
    real(RP) :: Zi_CONST

    namelist / PARAM_MKINIT_TKE / &
       TKE_CONST, &
       Zi_CONST

    integer :: k, i, j, iq
    integer :: ierr

    integer :: ldomID
    class(FEM_LocalMesh3D), pointer :: lcmesh3D
    !---------------------------------------------------------------------------

    TKE_CONST = EPS
    Zi_CONST = 100.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TKE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("tke_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("tke_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_TKE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TKE)

    if ( I_TKE > 0 ) then ! TB
       !$omp parallel do collapse(2)
       !$acc kernels
       do j = JSB, JEB
       do i = ISB, IEB
       do k = 1, KA
          if ( QTRC(k,i,j,I_TKE) == UNDEF ) then
             QTRC(k,i,j,I_TKE) = TKE_CONST
          end if
       enddo
       enddo
       enddo
       !$acc end kernels
    end if

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       Zi(i,j) = Zi_CONST
    end do
    end do
    !$acc end kernels

    ! BL
    call atmos_phy_bl_driver_mkinit( TKE_CONST )

    if ( ATMOS_sw_dyn_DGM ) then
      ldomID = 1
      lcmesh3D => fem_mesh3d%lcmesh_list(ldomID)
      do iq=QS_BL, QE_BL
         call MeshFieldUtil_interp_FVtoDG_3D( dyn_dg_vars%qtrc_vars(iq)%local(ldomID)%val, &
            QTRC(:,:,:,iq), lcmesh3D, fem_elem3d, IS, IE, IA, IHALO, JE, JE, JA, JHALO, KS, KE, KA, KHALO, &
            0 )
      end do
    end if

    return
  end subroutine tke_setup

  !-----------------------------------------------------------------------------
  !> Make initial state for MPACE experiment
  subroutine MKINIT_MPACE
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC
    use mod_atmos_phy_mp_vars, only: &
       QA_MP, &
       QS_MP, &
       QE_MP
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)
    real(RP) :: QCI1D(KA)
    real(RP) :: QNUM1D(KA)

    !real(RP) :: QHYD(KA,IA,JA,N_HYD)
    !real(RP) :: QNUM(KA,IA,JA,N_HYD)
    
    real(RP) :: bubbles(KA,IA,JA), temp

    integer  :: MPACE_fluctuationNumberLayers = 10
    real(RP) :: MPACE_fluctuation = 0.5D0
    namelist / PARAM_MKINIT_MPACE / &
       MPACE_fluctuation, &
       MPACE_fluctuationNumberLayers

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_MPACE",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_MPACE",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_MPACE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_MPACE",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_MPACE",*) 'Not appropriate names in namelist PARAM_MKINIT_MPACE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_MPACE)

    call MKINIT_common_read_sounding_mpace( RHO, VELX, VELY, POTT, QV1D, QCI1D, QNUM1D ) ! (out)

    ! initiate small fluctuation in potential temperature
    do i = ISB, IEB
    do k = KS, KS+MPACE_fluctuationNumberLayers-1
       call random_number(temp)
       bubbles(k,i,JSB) = MPACE_fluctuation * temp
    enddo
    enddo
    do j = JSB+1, JEB
      bubbles(:,:,j) = bubbles(:,:,JSB)
    enddo

    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = RHO(k) * bubbles(k,i,j)
       MOMX(k,i,j) = RHO(k) * VELX(k)
       MOMY(k,i,j) = RHO(k) * VELY(k)

       RHOT(k,i,j) = RHO(k) * ( POTT(k) )

       qv  (k,i,j) = QV1D(k)

       qc(k,i,j) = QCI1D(k) ! only the smallest (bin) cloud droplets are initialized
       if (QNUM1D(k) /= 0.0_RP) then
          nc(k,i,j) = QNUM1D(k)*1000000.0_RP
       else
          ! number concentration (cm-3->m-3) of smallest droplets
          ! assuming that each droplet mass is 1 x 10^-14 g (smallest bin)
          nc(k,i,j) = QCI1D(k)*RHO(k)/1.E-17_RP*1.E-6_RP *1000000.0_RP
       endif
    enddo
    enddo
    enddo
    
    !call ATMOS_PHY_MP_driver_qhyd2qtrc( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
    !                                    qv(:,:,:), QHYD(:,:,:,:), & ! [IN]
    !                                    QTRC(:,:,:,QS_MP:QE_MP),  & ! [OUT]
    !                                    QNUM=QNUM(:,:,:,:)        ) ! [IN]

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_MPACE

!  !-----------------------------------------------------------------------------
!  !> Make initial state of Box model experiment for collision-coalescence in AMPS
!  subroutine MKINIT_boxamps
!    use scale_const, only: &
!       Rdry  => CONST_Rdry, &
!       Rvap  => CONST_Rvap, &
!       CVdry => CONST_CVdry, &
!       CVvap => CONST_CVvap, &
!       CPdry => CONST_CPdry, &
!       CPvap => CONST_CPvap, &
!       CPliq => CONST_CL, &
!       EPSvap => CONST_EPSvap
!    use scale_atmos_hydrometeor, only: &
!       ATMOS_HYDROMETEOR_dry
!    use scale_atmos_thermodyn, only: &
!       ATMOS_THERMODYN_rhot2temp_pres, &
!       ATMOS_THERMODYN_temp_pres2pott
!    use mod_atmos_admin, only: &
!       ATMOS_PHY_MP_TYPE
!    use mod_atmos_phy_mp_vars, only: &
!       QS_MP, &
!       QE_MP
!    use scale_atmos_phy_mp_amps, only: &
!       ATMOS_PHY_MP_amps_init_qtrc_BOX
!    implicit none
!
!    real(RP) :: init_dens  = 1.00_RP   ![kg/m3]
!    real(RP) :: init_temp  = 300.00_RP ![K]
!    real(RP) :: init_pres  = 1.E+5_RP  ![Pa]
!    real(RP) :: init_ssliq = 0.0_RP   ![%]
!    real(RP) :: box_initial_conc = 100.0_RP
!    real(RP) :: box_m_knot = 4.18e-9_RP
!    integer  :: switch_rh_or_sq = 0
!
!    namelist / PARAM_MKINIT_BOXAMPS / &
!       init_dens, &
!       init_temp, &
!       init_pres, &
!       init_ssliq, &
!       switch_rh_or_sq, &
!       box_initial_conc, &
!       box_m_knot
!
!    real(RP) :: rtot (KA,IA,JA)
!    real(RP) :: cvtot(KA,IA,JA)
!    real(RP) :: cptot(KA,IA,JA)
!    real(RP) :: qdry
!    real(RP) :: psat, qsat
!    integer  :: i, j, k, ierr
!
!    real(RP) :: estbar(150), esitbar(111), psat_AMPS, qs_AMPS, qv_AMPS
!    !---------------------------------------------------------------------------
!
!    if ( ATMOS_PHY_MP_TYPE /= 'AMPS' ) then
!       LOG_INFO("MKINIT_boxamps",*) 'For [Box model of AMPS],'
!       LOG_INFO("MKINIT_boxamps",*) 'ATMOS_PHY_MP_TYPE should be AMPS. Stop! ', trim(ATMOS_PHY_MP_TYPE)
!       call PRC_abort
!    endif
!
!    if ( ATMOS_HYDROMETEOR_dry ) then
!       LOG_ERROR("MKINIT_boxamps",*) 'QV is not registered'
!       call PRC_abort
!    end if
!
!    LOG_NEWLINE
!    LOG_INFO("MKINIT_boxamps",*) 'Setup initial state'
!
!    !--- read namelist
!    rewind(IO_FID_CONF)
!    read(IO_FID_CONF,nml=PARAM_MKINIT_BOXAMPS,iostat=ierr)
!    if( ierr < 0 ) then !--- missing
!       LOG_INFO("MKINIT_boxamps",*) 'Not found namelist. Default used.'
!    elseif( ierr > 0 ) then !--- fatal error
!       LOG_ERROR("MKINIT_boxamps",*) 'Not appropriate names in namelist PARAM_MKINIT_BOXAMPS. Check!'
!       call PRC_abort
!    endif
!    LOG_NML(PARAM_MKINIT_BOXAMPS)
!
!    ! SCALE way
!    call SATURATION_psat_all( init_temp, psat )
!    qsat = EPSvap * psat / ( init_pres - ( 1.0_RP-EPSvap ) * psat )
!    ! AMPS way
!    call QSPARM2(estbar,esitbar)
!    psat_AMPS = get_sat_vapor_pres_lk2(1, init_temp, estbar, esitbar )/10.0_RP
!    qs_AMPS = psat_AMPS/(init_pres - psat_AMPS)*0.622D0
!    qv_AMPS = psat_AMPS*0.622_RP*init_ssliq/100.0_RP / (init_pres - psat_AMPS*init_ssliq/100.0_RP)
! 
!    do j = 1, JA
!    do i = 1, IA
!    do k = 1, KA
!       MOMX(k,i,j) = 0.0_RP
!       MOMY(k,i,j) = 0.0_RP
!       MOMZ(k,i,j) = 0.0_RP
!
!       if (switch_rh_or_sq == 1) then
!          ! from specific humidity
!          qv(k,i,j) = ( init_ssliq + 1.0_RP ) * qsat
!       else if (switch_rh_or_sq == 0) then
!          ! from relative humidity
!          !rh = 100*pq/(e+q)/s
!          !rh/100*(e+q)s = pq
!          !q (p - rh*s/100) = s*e*rh/100
!          !q = s*e*rh/100 / (p - rh*s/100)
!          qv(k,i,j) = psat_AMPS*0.622_RP*init_ssliq/100.0_RP / (init_pres - psat_AMPS*init_ssliq/100.0_RP)
!       endif
!       QTRC(k,i,j,QS_MP) = qv(k,i,j)
!       
!       qdry = 1.0 - qv(k,i,j) - box_initial_conc*box_m_knot*1000.0/init_dens
!       rtot (k,i,j) = Rdry  * qdry + Rvap  * qv(k,i,j)
!       cvtot(k,i,j) = CVdry * qdry + CVvap * qv(k,i,j) + CPliq*box_initial_conc*box_m_knot*1000.0/init_dens
!       cptot(k,i,j) = CPdry * qdry + CPvap * qv(k,i,j) + CPliq*box_initial_conc*box_m_knot*1000.0/init_dens
!       call ATMOS_THERMODYN_temp_pres2pott(init_temp, init_pres, cptot(k,i,j), rtot(k,i,j), pott(k,i,j))
!       DENS(k,i,j) = init_pres / (init_temp * rtot(k,i,j))
!       RHOT(k,i,j) = DENS(k,i,j) * pott(k,i,j)
!
!       !nc(k,i,j) = -1.0_RP ! temperory indicator for box model calculations
!    enddo
!    enddo
!    enddo
!
!    !call ATMOS_THERMODYN_rhot2temp_pres( KA, 1, KA, IA, 1, IA, JA, 1, JA, &
!    !                                     DENS(:,:,:), RHOT(:,:,:),                & ! (in)
!    !                                     rtot(:,:,:), cvtot(:,:,:), cptot(:,:,:), & ! (in)
!    !                                     temp(:,:,:), pres(:,:,:)                 ) ! (out)
!
!    call ATMOS_PHY_MP_amps_init_qtrc_BOX( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
!                                          DENS(:,:,:), QTRC(:,:,:,QS_MP+1:QE_MP), &
!                                          box_initial_conc, &
!                                          box_m_knot )
!
!    return
!  end subroutine MKINIT_boxamps

  SUBROUTINE QSPARM2(ESTBAR,ESITBAR)
!        based on Murphy and Koop (2005)
    implicit none
    real(RP), DIMENSION(150), INTENT(INOUT) ::  ESTBAR
    real(RP), DIMENSION(111), INTENT(INOUT) ::  ESITBAR
    real(RP) :: T
    integer  :: K, JD

    T=163.0_RP
    DO K = 1 , 111
       T = T + 1.0_RP
       ESITBAR(K) = exp(9.550426_RP - 5723.265_RP/T + 3.53068_RP*log(T) - 0.00728332_RP*T)
    ENDDO

    T=163.0_RP
    DO JD = 1, 150
       T = T + 1.0_RP
       ESTBAR(JD) = exp(54.842763_RP - 6763.22_RP/T - 4.210_RP*log(T) + 0.000367_RP*T &
                    + tanh(0.0415_RP*(T - 218.8_RP))*(53.878_RP - 1331.22_RP/T &
                    - 9.44523_RP*log(T) + 0.014025_RP*T))
    END DO
    RETURN
  END subroutine QSPARM2

  function get_sat_vapor_pres_lk2(phase, T, estbar, esitbar ) result(e_sat)
  !_______________________________________________________________________
  !     THIS ROUTINE COMPUTES SATURATION MIXING RATIO OVER ICE
  !     WATER BASED ON A TABLE LOOK UP PROCEEDURE DEVELOPED BY
  !     DERICKSON AND COTTON (1977) FOR THE LIQUID PHASE
  !_______________________________________________________________________
    implicit none
    integer,  intent(in) :: phase
    real(RP), intent(in) :: T
    real(RP), intent(in) :: estbar(150), esitbar(111)
    integer  :: I, J
    real(RP) :: wt,e_sat

    ! g/s^2/cm
    if( phase == 1 ) then
      I = MAX(1,MIN(int(T) - 163,149))
      wt = MAX(MIN(T-real(I + 163,RP),1.0_RP),0.0_RP)
      e_sat = (estbar(I)*(1.0-wt) + estbar(I+1)*wt)*10.0_RP
    else if( phase == 2 ) then
      J = MAX(1,MIN(int(T) - 163,110))
      wt = MAX(MIN(T-real(J + 163,RP),1.0_RP),0.0_RP)
      e_sat = (esitbar(J)*(1.0-wt) + esitbar(J+1)*wt)*10.0_RP
    end if
  end function get_sat_vapor_pres_lk2

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform + random disturbance )
  subroutine MKINIT_planestate
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =   0.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    real(RP) :: ENV_TLAPS    =   0.0_RP ! Lapse rate of THETA [K/m]
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =   0.0_RP ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =   0.0_RP ! amplitude of random disturbance theta
    real(RP) :: RANDOM_U     =   0.0_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =   0.0_RP ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =   0.0_RP ! amplitude of random disturbance RH

    namelist / PARAM_MKINIT_PLANESTATE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_THETA,    &
       ENV_TLAPS,    &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       RANDOM_THETA, &
       RANDOM_U,     &
       RANDOM_V,     &
       RANDOM_RH

    integer :: ierr
    integer :: k, i, j, itr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_planestate",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_PLANESTATE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_planestate",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR_CONT(*) 'Not appropriate names in namelist PARAM_MKINIT_PLANESTATE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_PLANESTATE)

    ! calc in dry condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       pott_sfc(i,j) = SFC_THETA
       pres_sfc(i,j) = SFC_PRES
    enddo
    enddo
    !$acc end kernels

    if ( ENV_THETA < 0.0_RP ) then ! use isa profile

       call PROFILE_isa( KA, KS, KE,      & ! [IN]
                         IA, ISB, IEB,    & ! [IN]
                         JA, JSB, JEB,    & ! [IN]
                         pott_sfc(:,:), & ! [IN]
                         pres_sfc(:,:), & ! [IN]
                         REAL_CZ (:,:,:), & ! [IN]
                         pott    (:,:,:)  ) ! [OUT]

    else

       !$acc kernels
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * REAL_CZ(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels

    endif

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       pott_sfc(i,j) = pott_sfc(i,j) + ( rndm(KS-1,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_THETA
       enddo
    enddo
    enddo
    !$acc end kernels

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       
       ! Calculate QV from RH. 
       ! Note that the RH consequently obtained by following calculations is not precisely identical with the RH set by namelist, 
       ! because the iteration is not performed in the calculation of qv and density is re-built after including moisture. 
       
       call SATURATION_pres2qsat_all( IA, ISB, IEB, JA, JSB, JEB, &
                                      temp_sfc(:,:), pres_sfc(:,:), & ! [IN]
                                      qsat_sfc(:,:)                 ) ! [OUT]

       call RANDOM_uniform(rndm) ! make random
       !$acc kernels
       do j = JSB, JEB
       do i = ISB, IEB
          qv_sfc(i,j) = max( 0.0_RP, SFC_RH + ( rndm(KS-1,i,j) * 2.0_RP  - 1.0_RP ) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(i,j)
       enddo
       enddo
       !$acc end kernels
    end if


    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       do itr = 1, NITER_RH
          call SATURATION_psat_all( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                    temp(:,:,:), & ! [IN]
                                    psat(:,:,:)  ) ! [OUT]

          !$acc kernels
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             qv(k,i,j) = max( 0.0_RP, ENV_RH + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_RH ) * 1.E-2_RP * psat(k,i,j) / ( dens(k,i,j) * Rvap * temp(k,i,j) )
          enddo
          enddo
          enddo
          !$acc end kernels

          ! make density & pressure profile in moist condition
          call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                     pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                                     pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                                     REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                                     DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]
       end do
    end if

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, min(IEB,IA-1)
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, min(JEB,JA-1)
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_planestate

  !-----------------------------------------------------------------------------
  !> Make initial state for tracer bubble experiment
  subroutine MKINIT_tracerbubble
   use scale_fem_element_base, only: ElementBase2D
   use scale_atmos_dyn_dgm_hydrostatic, only: &
      dg_hydrostatic_calc_basicstate_constPT => hydrostatic_calc_basicstate_constPT
    implicit none

#ifndef DRY
    ! Surface state
    real(RP)               :: SFC_THETA            ! surface potential temperature [K]
    real(RP)               :: SFC_PRES             ! surface pressure [Pa]
    ! Environment state
    real(RP)               :: ENV_THETA            ! potential temperature of environment [K]
    real(RP)               :: ENV_U     =   0.0_RP ! velocity u of environment [m/s]
    real(RP)               :: ENV_V     =   0.0_RP ! velocity v of environment [m/s]
    ! Bubble
    character(len=H_SHORT) :: SHAPE_PTracer = 'BUBBLE' ! BUBBLE or RECT
    real(RP)               :: BBL_PTracer   = 1.0_RP   ! extremum of passive tracer in bubble [kg/kg]
    namelist / PARAM_MKINIT_TRACERBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       ENV_U,     &
       ENV_V,     &
       SHAPE_PTracer, &
       BBL_PTracer

    real(RP), pointer :: shapeFac(:,:,:) => null()

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: k, i, j
    integer :: ierr

    type(FEM_LocalMesh3D), pointer :: lmesh
    integer :: ldomid, kelem, ke_x, ke_y, ke_z, ke2D, p
    class(ElementBase2D), pointer :: fem_elem2D
    type(FEM_MeshField3D), pointer :: shapeFac_DG

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_tracerbubble",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TRACERBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_tracerbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_tracerbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_TRACERBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TRACERBUBBLE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    do k = KS, KE
       pott(k,1,1) = ENV_THETA
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V       * DENS(k,1,1)
       RHOT(k,i,j) = pott(k,1,1) * DENS(k,1,1)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make tracer bubble
    select case(SHAPE_PTracer)
    case('BUBBLE')
       call MKINIT_common_BUBBLE_setup( bubble )
       shapeFac => bubble
       if ( ATMOS_sw_dyn_DGM ) then
         call MKINIT_common_BUBBLE_setup_DG( DG_bubble )
         shapeFac_DG => DG_bubble
       endif
    case('RECT')
       call MKINIT_common_RECT_setup( rect )
       shapeFac => rect
    case default
       LOG_ERROR("MKINIT_tracerbubble",*) 'SHAPE_PTracer=', trim(SHAPE_PTracer), ' cannot be used on advect. Check!'
       call PRC_abort
    end select

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       ptrc(k,i,j) = BBL_PTracer * shapeFac(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ATMOS_sw_dyn_DGM ) then
      do ldomid=1, fem_mesh3D%LOCAL_MESH_NUM
         lmesh => fem_mesh3D%lcmesh_list(ldomid)
         fem_elem2D => lmesh%lcmesh2D%refElem2D

         call dg_hydrostatic_calc_basicstate_constPT( DG_DENS_hyd%local(ldomid)%val, DG_PRES_hyd%local(ldomid)%val, &
            SFC_THETA, SFC_PRES, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%zlev(:,:),                         &
            lmesh, lmesh%refElem3D )
   
         !$omp parallel do collapse(3) private( kelem, ke2D, ke_x, ke_y, ke_z )
         do ke_y=1, lmesh%NeY
         do ke_x=1, lmesh%NeX
         do ke_z=1, lmesh%NeZ
            ke2D = ke_x + (ke_y-1)*lmesh%NeX
            kelem = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY
            
            DG_DDENS%local(ldomid)%val(:,kelem) = 0.0_RP            
            DG_DRHOT%local(ldomid)%val(:,kelem) = 0.0_RP
            DG_MOMX%local(ldomid)%val(:,kelem) = DG_DENS_hyd%local(ldomid)%val(:,kelem) * ENV_U  
            DG_MOMY%local(ldomid)%val(:,kelem) = DG_DENS_hyd%local(ldomid)%val(:,kelem) * ENV_V
            DG_MOMZ%local(ldomid)%val(:,kelem) = 0.0_RP

            DG_ptrc%local(ldomid)%val(:,kelem) = BBL_PTracer * shapeFac_DG%local(ldomid)%val(:,kelem)
         enddo
         enddo
         enddo            
      enddo
    endif

#endif

    return
  end subroutine MKINIT_tracerbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  !! Default values are following by Straka et al. (1993)
  !!  BBL_TEMP = -15.0_RP   ! in temperature  [K]
  !!  BBL_CZ   =   3.0E3_RP ! center location [m]: z
  !!  BBL_CX   =  19.2E3_RP ! center location [m]: x
  !!  BBL_CY   =   1.0E2_RP ! center location [m]: y
  !!  BBL_RZ   =   2.0E3_RP ! bubble radius   [m]: z
  !!  BBL_RX   =   4.0E3_RP ! bubble radius   [m]: x
  !!  BBL_RY   =   1.0E3_RP ! bubble radius   [m]: y
  subroutine MKINIT_coldbubble
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    ! Bubble
    real(RP) :: BBL_TEMP     = -15.0_RP ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_COLDBUBBLE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_THETA, &
       BBL_TEMP

    real(RP) :: RovCP

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_coldbubble",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_COLDBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_coldbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_coldbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_COLDBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_COLDBUBBLE)

    RovCP = Rdry / CPdry

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    do k = KS, KE
       pott(k,1,1) = ENV_THETA
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP

       ! make cold bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1)                                           &
                                   + BBL_TEMP * ( P00/pres(k,1,1) )**RovCP * bubble(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_coldbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for cold bubble experiment
  subroutine MKINIT_lambwave
    implicit none

    ! Surface state
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_TEMP     = 300.0_RP ! temperature of environment [K]
    ! Bubble
    real(RP) :: BBL_PRES     =  100._RP ! extremum of pressure in bubble [Pa]

    namelist / PARAM_MKINIT_LAMBWAVE / &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_TEMP,  &
       BBL_PRES

    real(RP) :: RovCP

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_lambwave",*) 'Setup initial state'

    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAMBWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_lambwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_lambwave",*) 'Not appropriate names in namelist PARAM_MKINIT_LAMBWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_LAMBWAVE)

    RovCP = Rdry / CPdry

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = SFC_PRES/(Rdry*ENV_TEMP) * exp( - GRAV/(Rdry*ENV_TEMP) * CZ(k) )
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make pressure bubble
       pres(k,i,j) = DENS(k,i,j) * ENV_TEMP * Rdry + BBL_PRES * bubble(k,i,j)

       RHOT(k,i,j) = DENS(k,i,j) * ENV_TEMP * ( P00/pres(k,i,j) )**RovCP
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_lambwave

  !-----------------------------------------------------------------------------
  !> Make initial state for gravity wave experiment
  !! Default values are following by Skamarock and Klemp (1994)
  subroutine MKINIT_gravitywave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U        =  20.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_BVF      =  0.01_RP ! Brunt Vaisala frequencies of environment [1/s]
    ! Bubble
    real(RP) :: BBL_THETA    =  0.01_RP ! extremum of potential temperature in bubble [K]

    namelist / PARAM_MKINIT_GRAVITYWAVE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       ENV_BVF,   &
       BBL_THETA

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_gravitywave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAVITYWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_gravitywave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_gravitywave",*) 'Not appropriate names in namelist PARAM_MKINIT_GRAVITYWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_GRAVITYWAVE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    do k = KS, KE
       pott(k,1,1) = SFC_THETA * exp( ENV_BVF*ENV_BVF / GRAV * CZ(k) )
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,1,1)
       MOMY(k,i,j) = ENV_V * DENS(k,1,1)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_gravitywave

  !-----------------------------------------------------------------------------
  !> Make initial state for Kelvin-Helmholtz wave experiment
  subroutine MKINIT_khwave
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA                  ! surface potential temperature [K]
    real(RP) :: SFC_PRES                   ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_L1_ZTOP    = 1900.0_RP ! top    height of the layer1 (low  THETA) [m]
    real(RP) :: ENV_L3_ZBOTTOM = 2100.0_RP ! bottom height of the layer3 (high THETA) [m]
    real(RP) :: ENV_L1_THETA   =  300.0_RP ! THETA in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_THETA   =  301.0_RP ! THETA in the layer3 (high THETA) [K]
    real(RP) :: ENV_L1_U       =    0.0_RP ! velocity u in the layer1 (low  THETA) [K]
    real(RP) :: ENV_L3_U       =   20.0_RP ! velocity u in the layer3 (high THETA) [K]
    ! Disturbance
    real(RP) :: RANDOM_U       =    0.0_RP ! amplitude of random disturbance u

    namelist / PARAM_MKINIT_KHWAVE / &
       SFC_THETA,      &
       SFC_PRES,       &
       ENV_L1_ZTOP,    &
       ENV_L3_ZBOTTOM, &
       ENV_L1_THETA,   &
       ENV_L3_THETA,   &
       ENV_L1_U,       &
       ENV_L3_U,       &
       RANDOM_U

    real(RP) :: fact

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_khwave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_KHWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_khwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_khwave",*) 'Not appropriate names in namelist PARAM_MKINIT_KHWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_KHWAVE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    do k = KS, KE
       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       pott(k,1,1) = ENV_L1_THETA * ( 1.0_RP - fact ) &
                   + ENV_L3_THETA * (          fact )
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       RHOT(k,i,j) = DENS(k,1,1) * pott(k,1,1)
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       fact = ( CZ(k)-ENV_L1_ZTOP ) / ( ENV_L3_ZBOTTOM-ENV_L1_ZTOP )
       fact = max( min( fact, 1.0_RP ), 0.0_RP )

       MOMX(k,i,j) = ( ENV_L1_U * ( 1.0_RP - fact )                 &
                     + ENV_L3_U * (          fact )                 &
                     + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_U &
                     ) * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_khwave

  !-----------------------------------------------------------------------------
  !> Make initial state for turbulence experiment
  subroutine MKINIT_turbulence
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =   0.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_THETA               ! potential temperature of environment [K]
    real(RP) :: ENV_TLAPS    = 4.E-3_RP ! Lapse rate of THETA [K/m]
    real(RP) :: ENV_U        =   5.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =   0.0_RP ! relative humidity of environment [%]
    ! Disturbance
    real(RP) :: RANDOM_THETA =   1.0_RP ! amplitude of random disturbance theta
    real(RP) :: RANDOM_U     =   0.0_RP ! amplitude of random disturbance u
    real(RP) :: RANDOM_V     =   0.0_RP ! amplitude of random disturbance v
    real(RP) :: RANDOM_RH    =   0.0_RP ! amplitude of random disturbance RH

    namelist / PARAM_MKINIT_TURBULENCE / &
       SFC_THETA,    &
       SFC_PRES,     &
       SFC_RH,       &
       ENV_THETA,    &
       ENV_TLAPS,    &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       RANDOM_THETA, &
       RANDOM_U,     &
       RANDOM_V,     &
       RANDOM_RH

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j, itr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_turbulence",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd
    ENV_THETA = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_TURBULENCE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_turbulence",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_turbulence",*) 'Not appropriate names in namelist PARAM_MKINIT_TURBULENCE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_TURBULENCE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    do k = KS, KE
       pott(k,1,1) = ENV_THETA + ENV_TLAPS * CZ(k)
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA + ( rndm(KS-1,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_THETA

       do k = KS, KE
          pott(k,i,j) = ENV_THETA + ENV_TLAPS * CZ(k) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_THETA
       enddo
    enddo
    enddo
    !$acc end kernels

    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       ! calc QV from RH
       call SATURATION_pres2qsat_all( IA, ISB, IEB, JA, JSB, JEB, &
                                      temp_sfc(:,:), pres_sfc(:,:), & ! [IN]
                                      qsat_sfc(:,:)                 ) ! [OUT]

       call RANDOM_uniform(rndm) ! make random
       !$acc kernels
       !$acc loop collapse(2) independent
       do j = JSB, JEB
       do i = ISB, IEB
          qv_sfc(i,j) = min( 0.0_RP, SFC_RH + ( rndm(KS-1,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_RH ) * 1.E-2_RP * qsat_sfc(1,1)
       enddo
       enddo
       !$acc end kernels

       do itr = 1, NITER_RH
          call SATURATION_psat_all( KA, KS, KE, IA, IS, JE, JA, JS, JE, &
                                    temp(:,:,:), & ! [IN]
                                    psat(:,:,:)  ) ! [OUT]

          !$acc kernels
          !$acc loop collapse(2) independent
          do j = JSB, JEB
          do i = ISB, IEB
          do k = KS, KE
             qv(k,i,j) = min( 0.0_RP, ENV_RH + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_RH ) * 1.E-2_RP * psat(k,i,j) / ( dens(k,i,j) * Rvap * temp(k,i,j) )
          enddo
          enddo
          enddo
          !$acc end kernels

          ! make density & pressure profile in moist condition
          call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                     pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                                     pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                                     REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                                     DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]
       end do
    end if

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = ( ENV_U + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_U ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = ( ENV_V + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_V ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_turbulence

  !-----------------------------------------------------------------------------
  !> Make initial state for cavity flow
  subroutine MKINIT_cavityflow
    implicit none

    ! Nondimenstional numbers for a cavity flow problem
    real(RP) :: REYNOLDS_NUM = 1.D03
    real(RP) :: MACH_NUM     = 3.D-2
    real(RP) :: Ulid         = 1.D01
    real(RP) :: PRES0        = 1.D05

    namelist / PARAM_MKINIT_CAVITYFLOW / &
         Ulid        , &
         PRES0       , &
         REYNOLDS_NUM, &
         MACH_NUM

    real(RP) :: DENS0
    real(RP) :: TEMP
    real(RP) :: Gam
    real(RP) :: Cs2

    integer :: k, i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_cavityflow",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_CAVITYFLOW,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_cavityflow",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_cavityflow",*) 'Not appropriate names in namelist PARAM_MKINIT_CAVITYFLOW. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_CAVITYFLOW)

    Gam   = CPdry / ( CPdry - Rdry )
    Cs2   = ( Ulid / MACH_NUM )**2
    TEMP  = Cs2   / ( Gam * Rdry )
    DENS0 = PRES0 / ( Rdry * TEMP )

    LOG_INFO("MKINIT_cavityflow",*) "DENS = ", DENS0
    LOG_INFO("MKINIT_cavityflow",*) "PRES = ", PRES0
    LOG_INFO("MKINIT_cavityflow",*) "TEMP = ", RHOT(10,10,4)/DENS0, TEMP
    LOG_INFO("MKINIT_cavityflow",*) "Ulid = ", Ulid
    LOG_INFO("MKINIT_cavityflow",*) "Cs   = ", sqrt(Cs2)

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = DENS0
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       PRES(k,i,j) = PRES0
       RHOT(k,i,j) = P00/Rdry * (P00/PRES0)**((Rdry - CPdry)/CPdry)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    MOMX(KE+1:KA,:,:) = DENS0 * Ulid
    !$acc end kernels

    return
  end subroutine MKINIT_cavityflow

  !-----------------------------------------------------------------------------
  !> Make initial state ( horizontally uniform )
  subroutine MKINIT_mountainwave
   use scale_atmos_dyn_dgm_hydrostatic, only: &
      dg_hydrostatic_calc_basicstate_constBVFreq => hydrostatic_calc_basicstate_constBVFreq
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA      ! surface potential temperature [K]
    real(RP) :: SFC_PRES       ! surface pressure [Pa]
    ! Environment state
    real(RP) :: ENV_U = 0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V = 0.0_RP ! velocity v of environment [m/s]

    real(RP) :: SCORER      = 2.E-3_RP ! Scorer parameter (~=N/U) [1/m]
    real(RP) :: BBL_PTracer =   0.0_RP ! extremum of passive tracer in bubble [kg/kg]

    namelist / PARAM_MKINIT_MOUNTAINWAVE / &
       SFC_THETA, &
       SFC_PRES,  &
       ENV_U,     &
       ENV_V,     &
       SCORER,    &
       BBL_PTracer

    real(RP) :: Ustar2, N2

    integer :: ierr
    integer :: k, i, j

    type(FEM_LocalMesh3D), pointer :: lmesh
    integer :: ldomid
    integer :: kelem, ke_x, ke_y, ke_z
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_mountainwave",*) 'Setup initial state'

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_MOUNTAINWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_mountainwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_mountainwave",*) 'Not appropriate names in namelist PARAM_MKINIT_MOUNTAINWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_MOUNTAINWAVE)

    ! calc in dry condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       Ustar2 = ENV_U * ENV_U + ENV_V * ENV_V
       N2     = Ustar2 * (SCORER*SCORER)

       pott(k,i,j) = SFC_THETA * exp( N2 / GRAV * REAL_CZ(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,i,j)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U       * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V       * DENS(k,i,j)
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ATMOS_sw_dyn_DGM ) then
      do ldomid=1, fem_mesh3D%LOCAL_MESH_NUM
         lmesh => fem_mesh3D%lcmesh_list(ldomid)

         call dg_hydrostatic_calc_basicstate_constBVFreq( DG_DENS_hyd%local(ldomid)%val, DG_PRES_hyd%local(ldomid)%val, &
            sqrt(N2), SFC_THETA, SFC_PRES, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%zlev(:,:),                   &
            lmesh, lmesh%refElem3D )
         
         !$omp parallel do collapse(3) private( kelem, ke_x, ke_y, ke_z )
         do ke_y=1, lmesh%NeY
         do ke_x=1, lmesh%NeX
         do ke_z=1, lmesh%NeZ
            kelem = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
            DG_DDENS%local(ldomid)%val(:,kelem) = 0.0_RP
            DG_DRHOT%local(ldomid)%val(:,kelem) = 0.0_RP
            DG_MOMX%local(ldomid)%val(:,kelem) = ENV_U * DG_DENS_hyd%local(ldomid)%val(:,kelem)
            DG_MOMY%local(ldomid)%val(:,kelem) = ENV_V * DG_DENS_hyd%local(ldomid)%val(:,kelem)
            DG_MOMZ%local(ldomid)%val(:,kelem) = 0.0_RP
         enddo
         enddo
         enddo
      enddo
    endif

    ! optional : add tracer bubble
    if ( BBL_PTracer > 0.0_RP ) then
       !$acc kernels
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          ptrc(k,i,j) = BBL_PTracer * bubble(k,i,j)
       enddo
       enddo
       enddo
       !$acc end kernels
    endif

    return
  end subroutine MKINIT_mountainwave

  !-----------------------------------------------------------------------------
  !> Make initial state
  !!
  !! The initial state for a baroclinic wave test described by Ullrich and Jablonowski(2012)
  !! is generated.
  subroutine MKINIT_barocwave
    use scale_const, only: &
      OHM => CONST_OHM,        &
      RPlanet => CONST_RADIUS, &
      GRAV    => CONST_GRAV
    use scale_prc
    use scale_atmos_grid_cartesC, only: &
         y0  => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y, &
         FYG => ATMOS_GRID_CARTESC_FYG
    implicit none

    ! Parameters for global domain size
    real(RP) :: Ly                      ! The domain size in y-direction [m]

    ! Parameters for inital stratification
    real(RP) :: REF_TEMP    = 288.E0_RP ! The reference temperature [K]
    real(RP) :: REF_PRES    = 1.E5_RP   ! The reference pressure [Pa]
    real(RP) :: LAPSE_RATE  = 5.E-3_RP  ! The lapse rate [K/m]

    ! Parameters associated with coriolis parameter on a beta-plane
    real(RP) :: Phi0Deg     = 45.E0_RP  ! The central latitude [degree_north]

    ! Parameters for background zonal jet
    real(RP) :: U0 = 35.E0_RP          ! The parameter associated with zonal jet maximum amplitude  [m/s]
    real(RP) :: b  = 2.E0_RP           ! The vertical half-width [1]

    ! Parameters for inital perturbation of zonal wind with a Gaussian profile
    !
    real(RP) :: Up  = 1.E0_RP         ! The maximum amplitude of zonal wind perturbation [m/s]
    real(RP) :: Lp  = 600.E3_RP       ! The width of Gaussian profile
    real(RP) :: Xc  = 2000.E3_RP      ! The center point (x) of inital perturbation
    real(RP) :: Yc  = 2500.E3_RP      ! The center point (y) of inital perturbation

    namelist / PARAM_MKINIT_BAROCWAVE / &
       REF_TEMP, REF_PRES, LAPSE_RATE, &
       phi0Deg,                        &
       U0, b,                          &
       Up, Lp, Xc, Yc

    real(RP) :: f0, beta0

    real(RP) :: geopot(KA,IA,JA)
    real(RP) :: eta(KA,IA,JA)
    real(RP) :: temp(KA,IA,JA)

    real(RP) :: y
    real(RP) :: ln_eta
    real(RP) :: del_eta
    real(RP) :: yphase
    real(RP) :: temp_vfunc
    real(RP) :: geopot_hvari

    logical :: converged

    integer :: ierr
    integer :: k, i, j

    integer :: itr

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    logical :: error

    integer,  parameter :: ITRMAX = 1000
    real(RP), parameter :: CONV_EPS = 1E-15_RP
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_barocwave",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BAROCWAVE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_barocwave",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_barocwave",*) 'Not appropriate names in namelist PARAM_MKINIT_BAROCWAVE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BAROCWAVE)

    Ly = FYG(JAG-JHALO) - FYG(JHALO)

    ! Set coriolis parameters
    f0 = 2.0_RP*OHM*sin(phi0Deg*PI/180.0_RP)
    beta0 = (2.0_RP*OHM/RPlanet)*cos(phi0Deg*PI/180.0_RP)

    ! Calculate eta(=p/p_s) level corresponding to z level of each (y,z) grid point
    ! using Newton's iteration method

    !$acc data create(eta)

    !$omp workshare
    !$acc kernels
    eta(:,:,:) = 1.0E-8_RP   ! Set first guess of eta
    !$acc end kernels
    !$omp end workshare

    error = .false.

    !$omp parallel do private(y,yphase,geopot_hvari,del_eta,itr,ln_eta,temp_vfunc,converged) reduction(.or.:error)
    !$acc kernels
    !$acc loop independent collapse(2) &
    !$acc private(work1,work2,work3) reduction(.or.:error)
    do j = JSB, JEB
    do i = ISB, IEB            ! Note that initial fields are zonaly symmetric

       y = CY(j)
       yphase  = 2.0_RP*PI*y/Ly

       ! Calc horizontal variation of geopotential height
       geopot_hvari = 0.5_RP*U0*(                                                                          &
            (f0 - beta0*y0)*(y - 0.5_RP*Ly*(1.0_RP + sin(yphase)/PI))                                      &
            + 0.5_RP*beta0*( y**2 - Ly*y/PI*sin(yphase) - 0.5_RP*(Ly/PI)**2*(cos(yphase) + 1.0_RP)         &
                             - Ly**2/3.0_RP                                                        )       &
            )

       ! Set surface pressure and temperature
       pres_sfc(i,j) = REF_PRES
       pott_sfc(i,j) = REF_TEMP - geopot_hvari/Rdry

       do k = KS, KE
          del_eta = 1.0_RP

          !-- The loop for iteration
          itr = 0
          do while( abs(del_eta) > CONV_EPS )
             ln_eta = log(eta(k,i,j))

             temp_vfunc = eta(k,i,j)**(Rdry*LAPSE_RATE/Grav)
             temp(k,i,j) = &
                  REF_TEMP*temp_vfunc &
                  + geopot_hvari/Rdry*(2.0_RP*(ln_eta/b)**2 - 1.0_RP)*exp(-(ln_eta/b)**2)
             geopot(k,i,j) = &
                  REF_TEMP*GRAV/LAPSE_RATE*(1.0_RP - temp_vfunc)  &
                  + geopot_hvari*ln_eta*exp(-(ln_eta/b)**2)

             del_eta = -  ( - Grav*CZ(k) + geopot(k,i,j) )    & ! <- F
                  &      *( - eta(k,i,j)/(Rdry*temp(k,i,j))   ) ! <- (dF/deta)^-1

             eta(k,i,j) = eta(k,i,j) + del_eta
             itr = itr + 1

             if ( itr > ITRMAX ) then
                LOG_ERROR("MKINIT_barocwave",*) "Fail the convergence of iteration. Check!"
                LOG_ERROR_CONT(*) "* (X,Y,Z)=", CX(i), CY(j), CZ(k)
                LOG_ERROR_CONT(*) "itr=", itr, "del_eta=", del_eta, "eta=", eta(k,i,j), "temp=", temp(k,i,j)
#ifdef _OPENACC
                error = .true.
#else
                call PRC_abort
#endif
             end if
          enddo !- End of loop for iteration ----------------------------

          PRES(k,i,j) = eta(k,i,j)*REF_PRES
          DENS(k,i,j) = PRES(k,i,j)/(Rdry*temp(k,i,j))
          pott(k,i,j) = temp(k,i,j)*eta(k,i,j)**(-Rdry/CPdry)

       enddo

       ! Make density & pressure profile in dry condition using the profile of
       ! potential temperature calculated above.
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  pott(:,i,j), qv(:,i,j), qc(:,i,j),                      & ! [IN]
                                  pres_sfc(i,j), pott_sfc(i,j), qv_sfc(i,j), qc_sfc(i,j), & ! [IN]
                                  REAL_CZ(:,i,j), REAL_FZ(:,i,j),                         & ! [IN]
#ifdef _OPENACC
                                  work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                                  DENS(:,i,j), temp(:,i,j), pres(:,i,j), temp_sfc(i,j),   & ! [OUT]
                                  converged                                               ) ! [OUT]
       if ( .not. converged ) then
          LOG_ERROR("MKINIT_barocwave",*) "failed to obtain a state in hydrostatic balance", i, j
#ifdef _OPENACC
          error = .true.
#else
          call PRC_abort
#endif
       end if

    enddo
    enddo
    !$acc end kernels

    if ( error ) then
       call PRC_abort
    end if

    !-----------------------------------------------------------------------------------

    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do k = KS, KE

       eta(k,IS,j) = pres(k,IS,j)/REF_PRES
       ln_eta = log(eta(k,IS,j))
       yphase = 2.0_RP*PI*CY(j)/Ly
!!$       PRES(k,IS:IE,j) = eta(k,IS,j)*REF_PRES
!!$       DENS(k,IS:IE,j) = PRES(k,IS,j)/(Rdry*temp(k,IS,j))
       !$acc loop independent
       do i = ISB, IEB
          if ( i .ne. IS ) then
             DENS(k,i,j) = DENS(k,IS,j)
             PRES(k,i,j) = PRES(k,IS,j)
          end if
          RHOT(k,i,j) = DENS(k,IS,j)*pott(k,IS,j) !temp(k,IS,j)*eta(k,IS,j)**(-Rdry/CPdry)
       end do
       !$acc loop independent
       do i = ISB, IEB
          MOMX(k,i,j) = DENS(k,IS,j)*(-U0*sin(0.5_RP*yphase)**2*ln_eta*exp(-(ln_eta/b)**2))
       end do
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    MOMY(:,:,:) = 0.0_RP
    MOMZ(:,:,:) = 0.0_RP
    !$acc end kernels

    !---------------------------------------------------------------------------------------

    ! Add the inital perturbation for zonal velocity
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = MOMX(k,i,j) &
           +  DENS(k,i,j)* Up*exp( - ((FX(i) - Xc)**2 + (CY(j) - Yc)**2)/Lp**2 )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

    return
  end subroutine MKINIT_barocwave

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubble
    use scale_const, only: &
       CVdry  => CONST_CVdry   
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry, &
       CV_VAPOR, CP_VAPOR       
    use scale_atmos_saturation, only: &
       ATMOS_SATURATION_pres2qsat_all       
    use scale_fem_element_base, only: ElementBase2D
    use scale_atmos_dyn_dgm_hydrostatic, only: &
       dg_hydrostatic_build_rho => hydrostaic_build_rho_XYZ, &
       dg_hydrostatic_calc_basicstate_constPT => hydrostatic_calc_basicstate_constPT
  
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 14.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =   1.0_RP ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j, itr

    type(FEM_LocalMesh3D), pointer :: lmesh
    integer :: ldomid
    integer :: kelem, ke_x, ke_y, ke_z
    integer :: ke2D, p, p2D, p3
    class(ElementBase2D), pointer :: fem_elem2D
    real(RP), allocatable :: DG_DENS(:), DG_POT(:,:,:,:), DG_QV_tmp(:)
    real(RP), allocatable :: DG_sfc_rhot(:), DG_bnd_SFC_PRES(:,:)
    integer :: iq_QV
    real(RP), allocatable :: DG_psat_z(:,:), DG_pres_z(:,:), DG_temp_z(:,:), DG_qsat_z(:,:)
    real(RP), allocatable :: DG_Rtot(:,:,:,:), DG_CPtot(:,:,:,:), DG_CPtot_ov_CVtot(:,:,:,:)
    real(RP), allocatable :: int_weight(:)
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_warmbubble",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_warmbubble",*) 'QV is not registered'
       call PRC_abort
    end if

    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_warmbubble",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_warmbubble",*) 'Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    !$acc loop seq
    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
!          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
          pott(k,1,1) = SFC_THETA + ENV_L2_TLAPS * ( CZ(k) - ENV_L1_ZTOP )
       else                                ! Layer 3
!          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( CZ(k) -CZ(k-1) )
          pott(k,1,1) = SFC_THETA + ENV_L2_TLAPS * ( ENV_L2_ZTOP - ENV_L1_ZTOP ) &
            + ENV_L3_TLAPS * ( CZ(k) - ENV_L2_ZTOP )
       endif
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    !$acc kernels
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]
    !$acc end kernels

    ! calc QV from RH
    !$acc kernels
    call SATURATION_pres2qsat_all( temp_sfc(1,1), pres_sfc(1,1), & ! [IN]
                                   qsat_sfc(1,1)                 ) ! [OUT]
    !$acc end kernels

    !$acc kernels
    qv_sfc(1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1)
    !$acc end kernels

    do itr = 1, NITER_RH
       !$acc kernels
       call SATURATION_psat_all( KA, KS, KE, &
                                 temp(:,1,1), & ! [IN]
                                 psat(:,1,1)  ) ! [OUT]
       !$acc end kernels

       !$acc kernels
       do k = KS, KE
          if( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 1 and 2
             qv(k,1,1) = ENV_RH * 1.E-2_RP * psat(k,1,1) / ( dens(k,1,1) * Rvap * temp(k,1,1) )
          else                            ! Layer 3
             qv(k,1,1) = 0.0_RP
          endif
       enddo
       !$acc end kernels

       ! make density & pressure profile in moist condition
       !$acc kernels
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                                  pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                                  CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                                  work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                                  DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                                  converged                                               ) ! [OUT]
       !$acc end kernels

    end do

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = qv(k,1,1)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( ATMOS_sw_dyn_DGM ) then
      do ldomid=1, fem_mesh3D%LOCAL_MESH_NUM
         lmesh => fem_mesh3D%lcmesh_list(ldomid)
         fem_elem2D => lmesh%lcmesh2D%refElem2D
         allocate( DG_POT(fem_elem3D%Np,lmesh%NeZ,lmesh%NeX,lmesh%NeY) )
         allocate( DG_DENS(fem_elem3D%Np), DG_sfc_rhot(fem_elem2D%Np) )
         allocate( DG_bnd_SFC_PRES(fem_elem2D%Np,lmesh%lcmesh2D%NeA) )
         allocate( int_weight(fem_elem3d%Np) )

         call dg_hydrostatic_calc_basicstate_constPT( DG_DENS_hyd%local(ldomid)%val, DG_PRES_hyd%local(ldomid)%val, &
            SFC_THETA, SFC_PRES, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%zlev(:,:),                         &
            lmesh, lmesh%refElem3D )
   
         !$omp parallel do collapse(3) private( kelem, ke2D, ke_x, ke_y, ke_z, DG_sfc_rhot )
         do ke_y=1, lmesh%NeY
         do ke_x=1, lmesh%NeX
         do ke_z=1, lmesh%NeZ
            ke2D = ke_x + (ke_y-1)*lmesh%NeX
            kelem = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY

            where ( lmesh%zlev(:,kelem) <= ENV_L1_ZTOP )
               DG_POT(:,ke_z,ke_x,ke_y) = SFC_THETA
            elsewhere ( lmesh%zlev(:,kelem) < ENV_L2_ZTOP  )
               DG_POT(:,ke_z,ke_x,ke_y) = SFC_THETA + ( lmesh%zlev(:,kelem) - ENV_L1_ZTOP ) * ENV_L2_TLAPS
            elsewhere
               DG_POT(:,ke_z,ke_x,ke_y) = SFC_THETA + ( ENV_L2_ZTOP - ENV_L1_ZTOP ) * ENV_L2_TLAPS       &
                                                    + ( lmesh%zlev(:,kelem) - ENV_L2_ZTOP ) * ENV_L3_TLAPS
            end where

            if ( ke_z == 1 ) then
               DG_sfc_rhot(:) = DG_DENS_hyd%local(ldomid)%val(fem_elem3D%Hslice(:,1),ke2D) * DG_POT(fem_elem3D%Hslice(:,1),ke_z,ke_x,ke_y)        
               DG_bnd_SFC_PRES(:,ke2D) = pres_sfc(1,1) ! P00 * ( Rdry * DG_sfc_rhot(:) / P00 )**( CpDry / CVdry )
            endif      
         enddo
         enddo
         enddo

         call dg_hydrostatic_build_rho( DG_DDENS%local(ldomid)%val, &
            DG_DENS_hyd%local(ldomid)%val, DG_PRES_hyd%local(ldomid)%val,                         &
            DG_POT, lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%zlev, lmesh, fem_elem3d, &
            bnd_SFC_PRES=DG_bnd_SFC_PRES )

         !$omp parallel do collapse(3) private(DG_DENS, kelem, ke_x,ke_y,ke_z)
         do ke_y=1, lmesh%NeY
         do ke_x=1, lmesh%NeX
         do ke_z=1, lmesh%NeZ
            kelem = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1) * lmesh%Ne2D

            DG_DENS(:) = DG_DENS_hyd%local(ldomid)%val(:,kelem) + DG_DDENS%local(ldomid)%val(:,kelem)
            DG_PRES_hyd%local(ldomid)%val(:,kelem) = P00 * ( Rdry / P00 * DG_DENS(:) * DG_POT(:,ke_z,ke_x,ke_y) )**(CpDry/CvDry)
            DG_DENS_hyd%local(ldomid)%val(:,kelem) = DG_DENS(:)
            DG_DDENS%local(ldomid)%val(:,kelem) = 0.0_RP
         enddo
         enddo
         enddo            
         
         call TRACER_inq_id( "QV", iq_QV )
     
         ! calc QV from RH
         allocate( DG_QV_tmp(fem_elem3d%Np) )
         allocate( DG_pres_z(fem_elem3d%Nnode_v,lmesh%NeZ), DG_temp_z(fem_elem3d%Nnode_v,lmesh%NeZ) )
         allocate( DG_psat_z(fem_elem3d%Nnode_v,lmesh%NeZ), DG_qsat_z(fem_elem3d%Nnode_v,lmesh%NeZ) )
         allocate( DG_Rtot(fem_elem3D%Np,lmesh%NeZ,lmesh%NeX,lmesh%NeY) )
         allocate( DG_CPtot(fem_elem3D%Np,lmesh%NeZ,lmesh%NeX,lmesh%NeY) )
         allocate( DG_CPtot_ov_CVtot(fem_elem3D%Np,lmesh%NeZ,lmesh%NeX,lmesh%NeY) )

         DG_Rtot(:,:,1,1) = Rdry
         DG_CPtot_ov_CVtot(:,:,1,1) = CpDry / CvDry
         do itr = 1, NITER_RH
            do ke_z = 1, lmesh%NeZ
            do p3=1, fem_elem3D%Nnode_v
               kelem = 1 + (ke_z - 1)*lmesh%Ne2D
               p  = 1 + (p3 - 1)*fem_elem3D%Nnode_h1D**2
               DG_pres_z(p3,ke_z) = P00 * ( DG_Rtot(p,ke_z,1,1) * ( DG_DENS_hyd%local(ldomid)%val(p,kelem) + DG_DDENS%local(ldomid)%val(p,kelem) ) * DG_POT(p,ke_z,1,1) / P00 )**DG_CPtot_ov_CVtot(p,ke_z,1,1)
               DG_temp_z(p3,ke_z) = DG_pres_z(p3,ke_z) &
                  / ( DG_Rtot(p,ke_z,1,1) * ( DG_DENS_hyd%local(ldomid)%val(p,kelem) + DG_DDENS%local(ldomid)%val(p,kelem) ) )
            end do
            end do
            call SATURATION_psat_all( &
               fem_elem3d%Nnode_v, 1, fem_elem3d%Nnode_v, lmesh%NeZ, 1, lmesh%NeZ, &
               DG_temp_z,                                                          & ! [IN]
               DG_psat_z                                                           ) ! [OUT]

            !$omp parallel do collapse(3) private(ke_z,ke_x,ke_y,kelem,ke2D,p3,p2D,p, DG_QV_tmp, DG_sfc_rhot, DG_DENS )
            do ke_y=1, lmesh%NeY
            do ke_x=1, lmesh%NeX
            do ke_z=1, lmesh%NeZ
               ke2D = ke_x + (ke_y-1)*lmesh%NeX
               kelem = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY

               DG_DENS(:) = DG_DENS_hyd%local(ldomid)%val(:,kelem) + DG_DDENS%local(ldomid)%val(:,kelem)

               do p3=1, fem_elem3d%Nnode_v
               do p2D=1, fem_elem3d%Nnode_h1D**2
                  p = p2D + (p3 - 1)*fem_elem3d%Nnode_h1D**2
                  if ( lmesh%zlev(p,kelem) > ENV_L2_ZTOP ) then
                     DG_QV_tmp(p) = 0.0_RP
                  else
                     DG_QV_tmp(p) = ENV_RH * 1.0E-2_RP * DG_psat_z(p3,ke_z) / ( DG_DENS(p) * Rvap * DG_temp_z(p3,ke_z) )
                  end if
               end do
               end do
               DG_QV%local(ldomid)%val(:,kelem)= DG_QV_tmp(:)

               DG_Rtot          (:,ke_z,ke_x,ke_y) = Rdry  * ( 1.0_RP - DG_QV_tmp(:) ) + Rvap     * DG_QV_tmp(:)
               DG_CPtot         (:,ke_z,ke_x,ke_y) = CPdry * ( 1.0_RP - DG_QV_tmp(:) ) + CP_VAPOR * DG_QV_tmp(:)
               DG_CPtot_ov_CVtot(:,ke_z,ke_x,ke_y) = DG_CPtot(:,ke_z,ke_x,ke_y)                         &
                                                / ( CVdry * ( 1.0_RP - DG_QV_tmp(:) ) + CV_VAPOR * DG_QV_tmp(:) ) 
               if ( ke_z == 1 ) then
                  DG_sfc_rhot(:) = DG_DENS(fem_elem3D%Hslice(:,1)) * DG_POT(fem_elem3D%Hslice(:,1),ke_z,ke_x,ke_y)        
                  DG_bnd_SFC_PRES(:,ke2D) = pres_sfc(1,1) ! P00 * ( DG_Rtot(fem_elem3D%Hslice(:,1),ke_z,ke_x,ke_y) * DG_sfc_rhot(:) / P00 )**DG_CPtot_ov_CVtot(fem_elem3D%Hslice(:,1),ke_z,ke_x,ke_y)
               endif      
            end do
            end do
            end do

            call dg_hydrostatic_build_rho( DG_DDENS%local(ldomid)%val, & ! (out)
               DG_DENS_hyd%local(ldomid)%val, DG_PRES_hyd%local(ldomid)%val, & ! (in)
               DG_POT, DG_Rtot, DG_CPtot_ov_CVtot,                           & ! (in)
               lmesh%pos_en(:,:,1), lmesh%pos_en(:,:,2), lmesh%zlev(:,:),    & ! (in)
               lmesh, fem_elem3d, bnd_SFC_PRES=DG_bnd_SFC_PRES ) ! (in)
         enddo

         !$omp parallel do collapse(3) private( kelem, ke2D, ke_x,ke_y,ke_z, k,i,j, DG_DENS, int_weight )
         do ke_y=1, lmesh%NeY
         do ke_x=1, lmesh%NeX
         do ke_z=1, lmesh%NeZ
            ke2D = ke_x + (ke_y-1)*lmesh%NeX
            kelem = ke2D + (ke_z-1)*lmesh%NeX*lmesh%NeY

            DG_DENS(:) =  DG_DENS_hyd%local(ldomid)%val(:,kelem) + DG_DDENS%local(ldomid)%val(:,kelem)
            DG_POT(:,ke_z,ke_x,ke_y) = DG_POT(:,ke_z,ke_x,ke_y) + BBL_THETA * DG_bubble%local(ldomid)%val(:,kelem)

            DG_DRHOT%local(ldomid)%val(:,kelem) = DG_DENS(:) * DG_POT(:,ke_z,ke_x,ke_y)               &
                        - P00 / Rdry * ( DG_PRES_hyd%local(ldomid)%val(:,kelem) / P00 )**(CVdry/CPdry)    
            DG_MOMX%local(ldomid)%val(:,kelem) = DG_DENS(:) * ENV_U  
            DG_MOMY%local(ldomid)%val(:,kelem) = DG_DENS(:) * ENV_V
            DG_MOMZ%local(ldomid)%val(:,kelem) = 0.0_RP

            int_weight(:) = fem_elem3d%IntWeight_lgl(:) * lmesh%J(:,kelem) * lmesh%Gsqrt(:,kelem)
            int_weight(:) = int_weight(:) / sum(int_weight(:) )

            k = KHALO + ke_z; i = IHALO + ke_x; j = JHALO + ke_y
            DENS(k,i,j) = sum( int_weight(:) * DG_DENS(:) )
            MOMX(k,i,j) = sum( int_weight(:) * DG_MOMX%local(ldomid)%val(:,kelem) )
            MOMY(k,i,j) = sum( int_weight(:) * DG_MOMY%local(ldomid)%val(:,kelem) )
            MOMZ(k,i,j) = sum( int_weight(:) * DG_MOMZ%local(ldomid)%val(:,kelem) )
            RHOT(k,i,j) = sum( int_weight(:) * DG_DENS(:) * DG_POT(:,ke_z,ke_x,ke_y)  )
            qv(k,i,j) = sum( int_weight(:) * DG_DENS(:) * DG_QV%local(ldomid)%val(:,kelem) ) / DENS(k,i,j)
         enddo
         enddo
         enddo

         deallocate( DG_POT, DG_DENS, DG_QV_tmp, DG_sfc_rhot, DG_bnd_SFC_PRES )
         deallocate( DG_pres_z, DG_temp_z, DG_psat_z, DG_qsat_z )
         deallocate( DG_Rtot, DG_CPtot, DG_CPtot_ov_CVtot )
         deallocate( int_weight )
      enddo
    endif

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_warmbubble

  !-----------------------------------------------------------------------------
  !> Make initial state for supercell experiment
  subroutine MKINIT_supercell
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_SUPERCELL / &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_supercell",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_supercell",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SUPERCELL,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_supercell",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_supercell",*) 'Not appropriate names in namelist PARAM_MKINIT_SUPERCELL. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SUPERCELL)

    call MKINIT_common_read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = RHO(k) * VELX(k)
       MOMY(k,i,j) = RHO(k) * VELY(k)

       ! make warm bubble
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo
    !$acc end kernels

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_supercell

  !-----------------------------------------------------------------------------
  !> Make initial state for squallline experiment
  subroutine MKINIT_squallline
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    real(RP) :: RANDOM_THETA =  0.01_RP
    real(RP) :: OFFSET_velx  = 12.0_RP
    real(RP) :: OFFSET_vely  = -2.0_RP

    namelist / PARAM_MKINIT_SQUALLLINE / &
       RANDOM_THETA,         &
       OFFSET_velx,          &
       OFFSET_vely

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_squallline",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_squallline",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SQUALLLINE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_squallline",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_squallline",*) 'Not appropriate names in namelist PARAM_MKINIT_SQUALLLINE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SQUALLLINE)

    call MKINIT_common_read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ( VELX(k) - OFFSET_velx ) * RHO(k)
       MOMY(k,i,j) = ( VELY(k) - OFFSET_vely ) * RHO(k)
       RHOT(k,i,j) = RHO(k) * ( POTT(k) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * RANDOM_THETA )
       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo
    !$acc end kernels

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_squallline

  !-----------------------------------------------------------------------------
  !> Make initial state by Weisman and Klemp (1982)
  subroutine MKINIT_wk1982
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA = 300.0_RP   ! surface pot. temperature [K]
    real(RP) :: SFC_PRES               ! surface pressure         [Pa]
    ! Parameter in Weisman and Klemp (1982)
    real(RP) :: TR_Z      = 12000.0_RP ! height           of tropopause  [m]
    real(RP) :: TR_THETA  =   343.0_RP ! pot. temperature at tropopause  [K]
    real(RP) :: TR_TEMP   =   213.0_RP ! temperature      at tropopause  [K]
    real(RP) :: SHEAR_Z   =  3000.0_RP ! center height of shear layer    [m]
    real(RP) :: SHEAR_U   =    15.0_RP ! velocity u over the shear layer [m/s]
    real(RP) :: QV0       =    14.0_RP ! maximum vapor mixing ration     [g/kg]
    ! Bubble
    real(RP) :: BBL_THETA = 3.D0 ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_WK1982 / &
       SFC_THETA, &
       SFC_PRES,  &
       TR_Z,      &
       TR_THETA,  &
       TR_TEMP,   &
       SHEAR_Z,   &
       SHEAR_U,   &
       QV0,       &
       BBL_THETA

    real(RP) :: rh    (KA,IA,JA)
    real(RP) :: rh_sfc(   IA,JA)

    integer :: ierr
    integer :: k, i, j, itr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_wk1982",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_wk1982",*) 'QV is not registered'
       call PRC_abort
    end if

    SFC_PRES  = Pstd

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WK1982,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_wk1982",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_wk1982",*) 'Not appropriate names in namelist PARAM_MKINIT_WK1982. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WK1982)

    ! calc in dry condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       pres_sfc(i,j) = SFC_PRES
       pott_sfc(i,j) = SFC_THETA

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             pott(k,i,j) = pott_sfc(i,j) &
                         + ( TR_THETA - pott_sfc(i,j) ) * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             pott(k,i,j) = TR_THETA * exp( GRAV * ( REAL_CZ(k,i,j) - TR_Z ) / CPdry / TR_TEMP )
          endif
       enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    !$acc data create(rh,rh_sfc)

    ! calc QV from RH
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       rh_sfc(i,j) = 1.0_RP - 0.75_RP * ( REAL_FZ(KS-1,i,j) / TR_Z )**1.25_RP

       do k = KS, KE
          if ( REAL_CZ(k,i,j) <= TR_Z ) then ! below initial cloud top
             rh(k,i,j) = 1.0_RP - 0.75_RP * ( REAL_CZ(k,i,j) / TR_Z )**1.25_RP
          else
             rh(k,i,j) = 0.25_RP
          endif
       enddo
    enddo
    enddo
    !$acc end kernels

    QV0 = QV0 * 1e-3_RP ! g/kg to kg/kg
    QV0 = QV0 / ( 1.0_RP + QV0 ) ! mixing ratio to specicic humidity

    call SATURATION_pres2qsat_all( IA, ISB, IEB, JA, JSB, JEB, &
                                   temp_sfc(:,:), pres_sfc(:,:), & ! [IN]
                                   qsat_sfc(:,:)                 ) ! [OUT]
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc(i,j) = min( rh_sfc(i,j) * qsat_sfc(i,j), QV0 )
    enddo
    enddo
    !$acc end kernels

    do itr = 1, NITER_RH
       call SATURATION_psat_all( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                 temp(:,:,:), & ! [IN]
                                 psat(:,:,:)  ) ! [OUT]
       !$acc kernels
       do j = JSB, JEB
       do i = ISB, IEB
       do k = KS, KE
          qv(k,i,j) = min( rh(k,i,j) * psat(k,i,j) / ( dens(k,i,j) * Rdry * temp(k,i,j) ), QV0 )
       enddo
       enddo
       enddo
       !$acc end kernels

       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                                  pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                                  pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                                  REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                                  DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]
    end do

    !$acc update host(pres(:,IS,JS),pott(:,IS,JS),rh(:,IS,JS),qv(:,IS,JS))
    do k = KS, KE
       LOG_INFO("MKINIT_wk1982",*) k, REAL_CZ(k,IS,JS), pres(k,IS,JS), pott(k,IS,JS), rh(k,IS,JS), qv(k,IS,JS)*1000
    enddo

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = SHEAR_U * tanh( REAL_CZ(k,i,j) / SHEAR_Z ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = 0.0_RP
       MOMZ(k,i,j) = 0.0_RP
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,i,j) * ( pott(k,i,j) + BBL_THETA * bubble(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call MKINIT_common_flux_setup

    !$acc end data

    return
  end subroutine MKINIT_wk1982

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF01
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w
    logical  :: USE_LWSET    = .false. ! use liq. water. static energy temp.?

    namelist / PARAM_MKINIT_RF01 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG,     &
       USE_LWSET

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint
    real(RP) :: GEOP_sw ! switch for geopotential energy correction

    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP ! pi/2

    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF01",*) 'Setup initial state'

    rewind(IO_FID_CONF)

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF01",*) 'QV is not registered'
       call PRC_abort
    end if

    read(IO_FID_CONF,nml=PARAM_MKINIT_RF01,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF01",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF01",*) 'Not appropriate names in namelist PARAM_MKINIT_RF01. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF01)

    if ( USE_LWSET ) then
       GEOP_sw = 1.0_RP
    else
       GEOP_sw = 0.0_RP
    endif

    !$acc data create(potl,lhv)

    ! calc in dry condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1017.8E2_RP ! [Pa]
       pott_sfc(i,j) = 289.0_RP    ! [K]

       do k = KS, KE
          velx(k,i,j) =   7.0_RP
          vely(k,i,j) =  -5.5_RP
          if ( CZ(k) < 820.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 289.0_RP - GRAV / CPdry * CZ(k) * GEOP_sw
          elseif( CZ(k) <= 860.0_RP ) then
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             potl(k,i,j) = ( 289.0_RP - GRAV / CPdry * CZ(k) * GEOP_sw                          ) * (0.5_RP-sint) &
                         + ( 297.5_RP+sign(abs(CZ(k)-840.0_RP)**(1.0_RP/3.0_RP),CZ(k)-840.0_RP) &
                           - GRAV / CPdry * CZ(k) * GEOP_sw                                     ) * (0.5_RP+sint)
          else
             potl(k,i,j) = 297.5_RP + ( CZ(k)-840.0_RP )**(1.0_RP/3.0_RP) &
                         - GRAV / CPdry * CZ(k) * GEOP_sw
          endif
       enddo

    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    ! calc in moist condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc  (i,j) = 9.0E-3_RP   ! [kg/kg]

       do k = KS, KE
          if    ( CZ(k) <   820.0_RP ) then ! below initial cloud top
             qall = 9.0E-3_RP
          elseif( CZ(k) <=  860.0_RP ) then ! boundary
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             qall = 9.0E-3_RP * (0.5_RP-sint) &
                  + 1.5E-3_RP * (0.5_RP+sint)
          elseif( CZ(k) <= 5000.0_RP ) then
             qall = 1.5E-3_RP
          else
             qall = 0.0_RP
          endif

          if    ( CZ(k) <=  600.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( CZ(k) < 820.0_RP ) then ! in the cloud
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact
          elseif( CZ(k) <= 860.0_RP ) then ! boundary
             sint = sin( pi2 * ( CZ(k)-840.0_RP ) / 20.0_RP ) * 0.5_RP
             fact = ( CZ(k)-600.0_RP ) / ( 840.0_RP-600.0_RP )
             qc(k,i,j) = 0.45E-3_RP * fact * (0.5_RP-sint)
          else
             qc(k,i,j) = 0.0_RP
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo
    !$acc end kernels

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMZ(k,i,j) = ( ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
       else
          MOMZ(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMX(k,i,j) = ( velx(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       else
          MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMY(k,i,j) = ( vely(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       else
          MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          RHOT(k,i,j) = ( pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * DENS(k,i,j)
       else
          RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 120.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       end if
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF01

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP) :: PERTURB_AMP  = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0 ! 0 -> no perturbation
                                 ! 1 -> perturbation for PT
                                 ! 2 -> perturbation for u,v,w

    namelist / PARAM_MKINIT_RF02 / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2
    real(RP) :: sint
    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2
    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF02",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF02",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF02",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF02",*) 'Not appropriate names in namelist PARAM_MKINIT_RF02. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF02)

    !$acc data create(potl,LHV)

    ! calc in dry condition
    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1017.8E2_RP   ! [Pa]
       pott_sfc(i,j) = 288.3_RP      ! [K]

       do k = KS, KE
          velx(k,i,j) =  3.0_RP + 4.3 * CZ(k)*1.E-3_RP
          vely(k,i,j) = -9.0_RP + 5.6 * CZ(k)*1.E-3_RP

          if ( CZ(k) < 775.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
          else if ( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
             potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP &
                         + ( 295.0_RP+sign(abs(CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),CZ(k)-795.0_RP) ) &
                         * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
          endif
       enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    ! calc in moist condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc(i,j) = 9.45E-3_RP

       do k = KS, KE
          if ( CZ(k) < 775.0_RP ) then ! below initial cloud top
             qall = 9.45E-3_RP ! [kg/kg]
          else if ( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
             qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
                  ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ! [kg/kg]
          endif

          if( CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( CZ(k) < 775.0_RP ) then
             fact = ( CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact
          elseif( CZ(k) <= 815.0_RP ) then
             sint = sin( pi2 * ( CZ(k)-795.0_RP )/20.0_RP )
             fact = ( CZ(k)-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.65E-3_RP * fact * (1.0_RP-sint) * 0.5_RP
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo
    !$acc end kernels

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMZ(k,i,j) = ( 0.0_RP + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
     else
       MOMZ(k,i,j) = 0.0_RP
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMX(k,i,j) = ( velx(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     else
       MOMX(k,i,j) = ( velx(k,i,j) ) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMY(k,i,j) = ( vely(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     else
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF02

  !-----------------------------------------------------------------------------
  !> Make initial state for stratocumulus
  subroutine MKINIT_DYCOMS2_RF02_DNS
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP) :: ZB  = 750.0_RP ! domain bottom
!   real(RP) :: ZT  = 900.0_RP ! domain top
    real(RP) :: CONST_U = 0.0_RP
    real(RP) :: CONST_V = 0.0_RP
    real(RP) :: PRES_ZB = 93060.0_RP
    real(RP) :: PERTURB_AMP  = 0.0_RP
    integer  :: RANDOM_LIMIT = 5
    integer  :: RANDOM_FLAG  = 0 ! 0 -> no perturbation
                                 ! 1 -> perturbation for PT
                                 ! 2 -> perturbation for u,v,w

    namelist / PARAM_MKINIT_RF02_DNS / &
       ZB, CONST_U, CONST_V,PRES_ZB,&
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]

    real(RP) :: qall ! QV+QC
    real(RP) :: fact
    real(RP) :: pi2

    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    pi2 = atan(1.0_RP) * 2.0_RP  ! pi/2

    LOG_NEWLINE
    LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_DYCOMS2_RF02_DNS",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RF02_DNS,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_DYCOMS2_RF02_DNS",*) 'Not appropriate names in namelist PARAM_MKINIT_RF02_DNS. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RF02_DNS)

    !$acc data create(potl,LHV)

    ! calc in dry condition
    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = PRES_ZB
!      pott_sfc(i,j) = 288.3_RP      ! [K]
!      qv_sfc  (i,j) = 9.45E-3_RP

       do k = KS, KE

          velx(k,i,j) = CONST_U
          vely(k,i,j) = CONST_V

!         if ( ZB+CZ(k) < 775.0_RP ) then ! below initial cloud top
          if ( ZB+CZ(k) <= 795.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 288.3_RP ! [K]
             qall = 9.45E-3_RP ! [kg/kg]
! necessary?
!         else if ( CZ(k) <= 815.0_RP ) then
!            sint = sin( pi2 * (CZ(k) - 795.0_RP)/20.0_RP )
!            potl(k,i,j) = 288.3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 295.0_RP+sign(abs(CZ(k)-795.0_RP)**(1.0_RP/3.0_RP),CZ(k)-795.0_RP) ) * (1.0_RP+sint)*0.5_RP
!            qall = 9.45E-3_RP * (1.0_RP-sint)*0.5_RP + &
!                  ( 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-CZ(k))/500.0_RP ) ) ) * (1.0_RP+sint)*0.5_RP
          else
             potl(k,i,j) = 295.0_RP + ( zb+CZ(k)-795.0_RP )**(1.0_RP/3.0_RP)
             qall = 5.E-3_RP - 3.E-3_RP * ( 1.0_RP - exp( (795.0_RP-(zb+CZ(k)))/500.0_RP ) ) ! [kg/kg]
          endif

          if( ZB+CZ(k) < 400.0_RP ) then
             qc(k,i,j) = 0.0_RP
          elseif( ZB+CZ(k) <= 795.0_RP ) then
             fact = ( (zb+CZ(k))-400.0_RP ) / ( 795.0_RP-400.0_RP )
             qc(k,i,j) = 0.8E-3_RP * fact
          else
             qc(k,i,j) = 0.0_RP
          endif
          qv(k,i,j) = qall - qc(k,i,j)

          !if(i==is.and.j==js)LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chkk',k,cz(k)+zb,qc(k,i,j),qv(k,i,j)
       enddo
    enddo
    enddo
    !$acc end kernels

    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk3',ks,ke
    ! extrapolation (temtative)
    !$acc kernels
    pott_sfc(:,:) = potl(ks,:,:)-0.5*(potl(ks+1,:,:)-potl(ks,:,:))
    qv_sfc  (:,:) = qv  (ks,:,:)-0.5*(qv  (ks+1,:,:)-qv  (ks,:,:))
    qc_sfc  (:,:) = qc  (ks,:,:)-0.5*(qc  (ks+1,:,:)-qc  (ks,:,:))
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]

    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMZ(k,i,j) = ( 0.0_RP + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
     else
       MOMZ(k,i,j) = 0.0_RP
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk8'
    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMX(k,i,j) = ( velx(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     else
       MOMX(k,i,j) = ( velx(k,i,j) ) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
    !$acc end kernels
    !LOG_INFO("MKINIT_DYCOMS2_RF02_DNS",*)'chk9'

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then
       MOMY(k,i,j) = ( vely(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     else
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
     if( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then
       RHOT(k,i,j) = ( pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                   * DENS(k,i,j)
     else
       RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
     endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 55.0E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

#endif
    return
  end subroutine MKINIT_DYCOMS2_RF02_DNS

  !-----------------------------------------------------------------------------
  !> Make initial state for RICO inter comparison
  subroutine MKINIT_RICO
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP):: PERTURB_AMP_PT = 0.1_RP
    real(RP):: PERTURB_AMP_QV = 2.5E-5_RP

    namelist / PARAM_MKINIT_RICO / &
       PERTURB_AMP_PT, &
       PERTURB_AMP_QV

    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall ! QV+QC
    real(RP) :: fact

    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_RICO",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_RICO",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_RICO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_RICO",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_RICO",*) 'Not appropriate names in namelist PARAM_MKINIT_RICO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_RICO)

    !$acc data create(potl,LHV)

    ! calc in moist condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1015.4E2_RP ! [Pa]
       pott_sfc(i,j) = 297.9_RP

       do k = KS, KE
          !--- potential temperature
          if ( CZ(k) < 740.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 297.9_RP
          else
             fact = ( CZ(k)-740.0_RP ) * ( 317.0_RP-297.9_RP ) / ( 4000.0_RP-740.0_RP )
             potl(k,i,j) = 297.9_RP + fact
          endif

          !--- horizontal wind velocity
          if ( CZ(k) <= 4000.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( -1.9_RP+9.9_RP ) / ( 4000.0_RP-0.0_RP )
             velx(k,i,j) =  -9.9_RP + fact
             vely(k,i,j) =  -3.8_RP
          else
             velx(k,i,j) =  -1.9_RP
             vely(k,i,j) =  -3.8_RP
          endif
       enddo

    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc  (i,j) = 16.0E-3_RP   ! [kg/kg]

       do k = KS, KE
          !--- mixing ratio of vapor
          if ( CZ(k) <= 740.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( 13.8E-3_RP-16.0E-3_RP ) / ( 740.0_RP-0.0_RP )
             qall = 16.0E-3_RP + fact
          elseif ( CZ(k) <= 3260.0_RP ) then ! boundary
             fact = ( CZ(k)-740.0_RP ) * ( 2.4E-3_RP-13.8E-3_RP ) / ( 3260.0_RP-740.0_RP )
             qall = 13.8E-3_RP + fact
          elseif( CZ(k) <= 4000.0_RP ) then
             fact = ( CZ(k)-3260.0_RP ) * ( 1.8E-3_RP-2.4E-3_RP ) / ( 4000.0_RP-3260.0_RP )
             qall = 2.4E-3_RP + fact
          else
             qall = 0.0_RP
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo
    !$acc end kernels

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       temp(k,i,j) = temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j)
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA  ,i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       RHOT(k,i,j) = ( pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP )*PERTURB_AMP_PT ) * DENS(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       qv(k,i,j) = qv(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP_QV
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 70.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

#endif
    return
  end subroutine MKINIT_RICO

  !-----------------------------------------------------------------------------
  !> Make initial state for BOMEX inter comparison
  subroutine MKINIT_BOMEX
    use scale_const, only: &
       Rdry  => CONST_Rdry,  &
       Rvap  => CONST_Rvap,  &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap, &
       CL    => CONST_CL
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

#ifndef DRY
    real(RP):: PERTURB_AMP_PT = 0.1_RP
    real(RP):: PERTURB_AMP_QV = 2.5E-5_RP

    namelist / PARAM_MKINIT_BOMEX / &
       PERTURB_AMP_PT, &
       PERTURB_AMP_QV

    real(RP) :: LHV (KA,IA,JA) ! latent heat of vaporization [J/kg]
    real(RP) :: potl(KA,IA,JA) ! liquid potential temperature
    real(RP) :: qall ! QV+QC
    real(RP) :: fact

    real(RP) :: qdry, Rtot, CPtot

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_BOMEX",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_BOMEX",*) 'QV is not registered'
       call PRC_abort
    end if

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BOMEX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_BOMEX",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_BOMEX",*) 'Not appropriate names in namelist PARAM_MKINIT_BOMEX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BOMEX)

    !$acc data create(potl,LHV)

    ! calc in moist condition
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB

       pres_sfc(i,j) = 1015.E2_RP ! [Pa]
       pott_sfc(i,j) = 299.1_RP

       do k = KS, KE
          !--- potential temperature
          if ( CZ(k) < 520.0_RP ) then ! below initial cloud top
             potl(k,i,j) = 298.7_RP
          elseif( CZ(k) < 1480.0_RP ) then
             fact = ( CZ(k)-520.0_RP ) * ( 302.4_RP-298.7_RP ) / ( 1480.0_RP-520.0_RP )
             potl(k,i,j) = 298.7_RP + fact
          elseif( CZ(k) < 2000.0_RP ) then
             fact = ( CZ(k)-1480.0_RP ) * ( 308.2_RP-302.4_RP ) / ( 2000.0_RP-1480.0_RP )
             potl(k,i,j) = 302.4_RP + fact
          else
             fact = ( CZ(k)-2000.0_RP ) * 3.65E-3_RP
             potl(k,i,j) = 308.2_RP + fact
          endif

          !--- horizontal wind velocity
          if ( CZ(k) <= 700.0_RP ) then ! below initial cloud top
             velx(k,i,j) =  -8.75_RP
             vely(k,i,j) =   0.0_RP
          else
             fact = 1.8E-3_RP * ( CZ(k)-700.0_RP )
             velx(k,i,j) =  -8.75_RP + fact
             vely(k,i,j) =  0.0_RP
          endif
       enddo

    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               potl(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       qv_sfc(i,j) = 22.45E-3_RP   ! [kg/kg]

       do k = KS, KE
          !--- mixing ratio of vapor
          if ( CZ(k) <= 520.0_RP ) then ! below initial cloud top
             fact = ( CZ(k)-0.0_RP ) * ( 16.3E-3_RP-17.0E-3_RP ) / ( 520.0_RP-0.0_RP )
             qall = 17.0E-3_RP + fact
          elseif ( CZ(k) <= 1480.0_RP ) then ! boundary
             fact = ( CZ(k)-520.0_RP ) * ( 10.7E-3_RP-16.3E-3_RP ) / ( 1480.0_RP-520.0_RP )
             qall = 16.3E-3_RP + fact
          elseif( CZ(k) <= 2000.0_RP ) then
             fact = ( CZ(k)-1480.0_RP ) * ( 4.2E-3_RP-10.7E-3_RP ) / ( 2000.0_RP-1480.0_RP )
             qall = 10.7E-3_RP + fact
          else
             fact = ( CZ(k)-2000.0_RP ) * ( -1.2E-6_RP )
             qall = 4.2E-3_RP + fact
          endif

          qv(k,i,j) = qall - qc(k,i,j)
       enddo

    enddo
    enddo
    !$acc end kernels

    call HYDROMETEOR_LHV( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                          temp(:,:,:), LHV(:,:,:) )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       qdry = 1.0_RP - qv(k,i,j) - qc(k,i,j)
       Rtot = Rdry * qdry + Rvap * qv(k,i,j)
       CPtot = CPdry * qdry + CPvap * qv(k,i,j) + CL * qc(k,i,j)
       pott(k,i,j) = ( temp(k,i,j) + LHV(k,i,j) / CPdry * qc(k,i,j) ) * ( P00 / pres(k,i,j) )**(Rtot/CPtot)
    enddo
    enddo
    enddo
    !$acc end kernels

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, IA, ISB, IEB, JA, JSB, JEB, &
                               pott(:,:,:), qv(:,:,:), qc(:,:,:),                      & ! [IN]
                               pres_sfc(:,:), pott_sfc(:,:), qv_sfc(:,:), qc_sfc(:,:), & ! [IN]
                               REAL_CZ(:,:,:), REAL_FZ(:,:,:), AREA(:,:),              & ! [IN]
                               DENS(:,:,:), temp(:,:,:), pres(:,:,:), temp_sfc(:,:)    ) ! [OUT]


    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA  ,i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call COMM_vars8( DENS(:,:,:), 1 )
    call COMM_wait ( DENS(:,:,:), 1 )

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMZ(k,i,j) = 0.0_RP
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMX(k,i,j) = velx(k,i,j) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       MOMY(k,i,j) = vely(k,i,j) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if( CZ(k) <= 1600.0_RP ) then !--- lowest 40 model layer when dz=40m
         RHOT(k,i,j) = ( pott(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP_PT ) * DENS(k,i,j)
       else
         RHOT(k,i,j) = pott(k,i,j) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if( CZ(k) <= 1600.0_RP ) then !--- lowest 40 model layer when dz=40m
          qv(k,i,j) = qv(k,i,j) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP_QV
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( qc(k,i,j) > 0.0_RP ) then
          nc(k,i,j) = 70.E6_RP / DENS(k,i,j) ! [number/m3] / [kg/m3]
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc end data

#endif
    return
  end subroutine MKINIT_BOMEX

  !-----------------------------------------------------------------------------
  !> Make initial state ( ocean variables )
  subroutine MKINIT_oceancouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_oceancouple",*) 'Setup initial state'

    call MKINIT_common_flux_setup

    call MKINIT_common_ocean_setup

    return
  end subroutine MKINIT_oceancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( land variables )
  subroutine MKINIT_landcouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_landcouple",*) 'Setup initial state'

    call MKINIT_common_flux_setup

    call MKINIT_common_land_setup

    return
  end subroutine MKINIT_landcouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( urban variables )
  subroutine MKINIT_urbancouple
    implicit none

    LOG_NEWLINE
    LOG_INFO("MKINIT_urbancouple",*) 'Setup initial state'

    call MKINIT_common_flux_setup

    call MKINIT_common_urban_setup

    return
  end subroutine MKINIT_urbancouple

  !-----------------------------------------------------------------------------
  !> Make initial state ( sea breeze )
  subroutine MKINIT_seabreeze
    use scale_landuse, only: &
       LANDUSE_frac_land, &
       LANDUSE_calc_fact, &
       LANDUSE_fillhalo
    use scale_atmos_grid_cartesC, only: &
       DOMAIN_CENTER_X => ATMOS_GRID_CARTESC_DOMAIN_CENTER_X
    use scale_land_grid_cartesC_real, only: &
       LAND_GRID_CARTESC_REAL_set_areavol
    use scale_ocean_grid_cartesC_real, only: &
       OCEAN_GRID_CARTESC_REAL_set_areavol
    implicit none

    real(RP) :: LAND_SIZE

    namelist / PARAM_MKINIT_SEABREEZE / &
       LAND_SIZE

    integer :: ierr
    integer :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_seabreeze",*) 'Setup initial state'

    LAND_SIZE = 0.0_RP

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SEABREEZE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_seabreeze",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_seabreeze",*) 'Not appropriate names in namelist PARAM_MKINIT_SEABREEZE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SEABREEZE)

    call MKINIT_common_flux_setup

    call MKINIT_common_land_setup

    call MKINIT_common_ocean_setup

    ! make landuse conditions
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       if ( abs( CX(i) - DOMAIN_CENTER_X ) < LAND_SIZE ) then
          LANDUSE_frac_land(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j) = 0.0_RP
       endif
    enddo
    enddo
    !$acc end kernels

    ! calculate landuse factors
    call LANDUSE_fillhalo( FILL_BND=.true. )
    call LANDUSE_calc_fact

    call LAND_GRID_CARTESC_REAL_set_areavol
    call OCEAN_GRID_CARTESC_REAL_set_areavol

    return
  end subroutine MKINIT_seabreeze

  !-----------------------------------------------------------------------------
  !> Make initial state ( heat island )
  subroutine MKINIT_heatisland
    use scale_prc_cartesC, only: &
       PRC_NUM_X
    use scale_landuse, only: &
       LANDUSE_frac_land,  &
       LANDUSE_frac_urban, &
       LANDUSE_calc_fact,  &
       LANDUSE_fillhalo
    implicit none

    real(RP) :: dist

    integer  :: i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_heatisland",*) 'Setup initial state'

    call MKINIT_common_flux_setup

    call MKINIT_common_land_setup

    call MKINIT_common_urban_setup

    ! 1/9 size of domain
    dist = ( CXG(IMAX*PRC_NUM_X) - CXG(1) ) / 9.0_RP

    ! make landuse conditions
    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       if (       CX(i) >= dist * 4.0_RP &
            .AND. CX(i) <  dist * 5.0_RP ) then
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 1.0_RP
       else
          LANDUSE_frac_land(i,j)  = 1.0_RP
          LANDUSE_frac_urban(i,j) = 0.0_RP
       endif
    enddo
    enddo
    !$acc end kernels

    ! calculate landuse factors
    call LANDUSE_fillhalo( FILL_BND=.true. )
    call LANDUSE_calc_fact

    return
  end subroutine MKINIT_heatisland

  !-----------------------------------------------------------------------------
  !> Make initial state for grayzone experiment
  subroutine MKINIT_grayzone
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP) :: RHO(KA)
    real(RP) :: VELX(KA)
    real(RP) :: VELY(KA)
    real(RP) :: POTT(KA)
    real(RP) :: QV1D(KA)

    real(RP) :: PERTURB_AMP = 0.0_RP
    integer  :: RANDOM_LIMIT = 0
    integer  :: RANDOM_FLAG  = 0       ! 0 -> no perturbation
                                       ! 1 -> petrurbation for pt
                                       ! 2 -> perturbation for u, v, w

    namelist / PARAM_MKINIT_GRAYZONE / &
       PERTURB_AMP,     &
       RANDOM_LIMIT,    &
       RANDOM_FLAG

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_grayzone",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_grayzone",*) 'QV is not registered'
       call PRC_abort
    end if

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_GRAYZONE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_grayzone",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_grayzone",*) 'Not appropriate names in namelist PARAM_MKINIT_GRAYZONE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_GRAYZONE)

    call MKINIT_common_read_sounding( RHO, VELX, VELY, POTT, QV1D ) ! (out)

    !$acc kernels
!   do j = JS, JE
!   do i = IS, IE
    do j = 1, ja
    do i = 1, ia
    do k = KS, KE
       DENS(k,i,j) = RHO(k)
!      MOMZ(k,i,j) = 0.0_RP
!      MOMX(k,i,j) = RHO(k) * VELX(k)
!      MOMY(k,i,j) = RHO(k) * VELY(k)

!      RHOT(k,i,j) = RHO(k) * POTT(k)
       qv  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo
    !$acc end kernels

    !$acc kernels
    !$acc loop collapse(2) independent
    do j = JSB, JEB
    do i = ISB, IEB
       DENS(   1:KS-1,i,j) = DENS(KS,i,j)
       DENS(KE+1:KA,  i,j) = DENS(KE,i,j)
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMZ(k,i,j) = ( ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k+1,i,j) + DENS(k,i,j) )
       else
          MOMZ(k,i,j) = 0.0_RP
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMX(k,i,j) = ( velx(k) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       else
          MOMX(k,i,j) = velx(k) * 0.5_RP * ( DENS(k,i+1,j) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 2 .AND. k <= RANDOM_LIMIT ) then ! below initial cloud top
          MOMY(k,i,j) = ( vely(k) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       else
          MOMY(k,i,j) = vely(k) * 0.5_RP * ( DENS(k,i,j+1) + DENS(k,i,j) )
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    call RANDOM_uniform(rndm) ! make random
    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       if ( RANDOM_FLAG == 1 .and. k <= RANDOM_LIMIT ) then ! below initial cloud top
          RHOT(k,i,j) = ( pott(k) + ( rndm(k,i,j) * 2.0_RP - 1.0_RP ) * PERTURB_AMP ) &
                      * DENS(k,i,j)
       else
          RHOT(k,i,j) = pott(k) * DENS(k,i,j)
       endif
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_grayzone

  !-----------------------------------------------------------------------------
  !> Make initial state of Box model experiment for zerochemical module
  subroutine MKINIT_boxaero
    use scale_const, only: &
       Rdry  => CONST_Rdry, &
       Rvap  => CONST_Rvap, &
       CVdry => CONST_CVdry, &
       CVvap => CONST_CVvap, &
       CPdry => CONST_CPdry, &
       CPvap => CONST_CPvap
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_rhot2temp_pres
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE
    implicit none

    real(RP) :: init_dens  = 1.12_RP   ![kg/m3]
    real(RP) :: init_temp  = 298.18_RP ![K]
    real(RP) :: init_pres  = 1.E+5_RP  ![Pa]
    real(RP) :: init_ssliq = 0.01_RP   ![%]

    namelist / PARAM_MKINIT_BOXAERO / &
       init_dens, &
       init_temp, &
       init_pres, &
       init_ssliq

    real(RP) :: rtot (KA,IA,JA)
    real(RP) :: cvtot(KA,IA,JA)
    real(RP) :: cptot(KA,IA,JA)
    real(RP) :: qdry
    real(RP) :: qsat
    integer  :: i, j, k, ierr
    !---------------------------------------------------------------------------

    if ( ATMOS_PHY_AE_TYPE /= 'KAJINO13' ) then
       LOG_INFO("MKINIT_boxaero",*) 'For [Box model of aerosol],'
       LOG_INFO("MKINIT_boxaero",*) 'ATMOS_PHY_AE_TYPE should be KAJINO13. Stop! ', trim(ATMOS_PHY_AE_TYPE)
       call PRC_abort
    endif

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_boxaero",*) 'QV is not registered'
       call PRC_abort
    end if

    LOG_NEWLINE
    LOG_INFO("MKINIT_boxaero",*) 'Setup initial state'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_BOXAERO,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_boxaero",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_boxaero",*) 'Not appropriate names in namelist PARAM_MKINIT_BOXAERO. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_BOXAERO)

    call SATURATION_pres2qsat_all( init_temp, init_pres, qsat )

    !$acc data create(rtot,cvtot,cptot)

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       DENS(k,i,j) = init_dens
       MOMX(k,i,j) = 0.0_RP
       MOMY(k,i,j) = 0.0_RP
       MOMZ(k,i,j) = 0.0_RP
       pott(k,i,j) = init_temp * ( P00/init_pres )**(Rdry/CPdry)
       RHOT(k,i,j) = init_dens * pott(k,i,j)

       qv(k,i,j) = ( init_ssliq + 1.0_RP ) * qsat

       qdry = 1.0 - qv(k,i,j)
       rtot (k,i,j) = Rdry  * qdry + Rvap  * qv(i,i,j)
       cvtot(k,i,j) = CVdry * qdry + CVvap * qv(i,i,j)
       cptot(k,i,j) = CPdry * qdry + CPvap * qv(i,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    call ATMOS_THERMODYN_rhot2temp_pres( KA, 1, KA, IA, 1, IA, JA, 1, JA, &
                                         dens(:,:,:), RHOT(:,:,:),                & ! (in)
                                         rtot(:,:,:), cvtot(:,:,:), cptot(:,:,:), & ! (in)
                                         temp(:,:,:), pres(:,:,:)                 ) ! (out)

    !$acc end data

    return
  end subroutine MKINIT_boxaero

  !-----------------------------------------------------------------------------
  !> Make initial state for warm bubble experiment
  subroutine MKINIT_warmbubbleaero
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    ! Surface state
    real(RP) :: SFC_THETA               ! surface potential temperature [K]
    real(RP) :: SFC_PRES                ! surface pressure [Pa]
    real(RP) :: SFC_RH       =  80.0_RP ! surface relative humidity [%]
    ! Environment state
    real(RP) :: ENV_U        =   0.0_RP ! velocity u of environment [m/s]
    real(RP) :: ENV_V        =   0.0_RP ! velocity v of environment [m/s]
    real(RP) :: ENV_RH       =  80.0_RP ! Relative Humidity of environment [%]
    real(RP) :: ENV_L1_ZTOP  =  1.E3_RP ! top height of the layer1 (constant THETA)       [m]
    real(RP) :: ENV_L2_ZTOP  = 14.E3_RP ! top height of the layer2 (small THETA gradient) [m]
    real(RP) :: ENV_L2_TLAPS = 4.E-3_RP ! Lapse rate of THETA in the layer2 (small THETA gradient) [K/m]
    real(RP) :: ENV_L3_TLAPS = 3.E-2_RP ! Lapse rate of THETA in the layer3 (large THETA gradient) [K/m]
    ! Bubble
    real(RP) :: BBL_THETA    =   1.0_RP ! extremum of temperature in bubble [K]

    namelist / PARAM_MKINIT_WARMBUBBLE / &
       SFC_THETA,    &
       SFC_PRES,     &
       ENV_U,        &
       ENV_V,        &
       ENV_RH,       &
       ENV_L1_ZTOP,  &
       ENV_L2_ZTOP,  &
       ENV_L2_TLAPS, &
       ENV_L3_TLAPS, &
       BBL_THETA

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr
    integer :: k, i, j, itr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("MKINIT_warmbubbleaero",*) 'Setup initial state'

    if ( ATMOS_HYDROMETEOR_dry ) then
       LOG_ERROR("MKINIT_warmbubbleaero",*) 'QV is not registerd'
       call PRC_abort
    end if


    SFC_THETA = THETAstd
    SFC_PRES  = Pstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_WARMBUBBLE,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("MKINIT_warmbubbleaero",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("MKINIT_warmbubbleaero",*) 'Not appropriate names in namelist PARAM_MKINIT_WARMBUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_WARMBUBBLE)

    ! calc in dry condition
    !$acc kernels
    pres_sfc(1,1) = SFC_PRES
    pott_sfc(1,1) = SFC_THETA
    !$acc end kernels

    !$acc kernels
    !$acc loop seq
    do k = KS, KE
       if    ( CZ(k) <= ENV_L1_ZTOP ) then ! Layer 1
          pott(k,1,1) = SFC_THETA
       elseif( CZ(k) <  ENV_L2_ZTOP ) then ! Layer 2
          pott(k,1,1) = pott(k-1,1,1) + ENV_L2_TLAPS * ( CZ(k)-CZ(k-1) )
       else                                ! Layer 3
          pott(k,1,1) = pott(k-1,1,1) + ENV_L3_TLAPS * ( CZ(k)-CZ(k-1) )
       endif
    enddo
    !$acc end kernels

    ! make density & pressure profile in dry condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                               pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                               CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                               DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                               converged                                               ) ! [OUT]

    ! calc QV from RH
    call SATURATION_pres2qsat_all( temp_sfc(1,1), pres_sfc(1,1), & ! [IN]
                                   qsat_sfc(1,1)                 ) ! [OUT]
    !$acc kernels
    qv_sfc(1,1) = SFC_RH * 1.E-2_RP * qsat_sfc(1,1)
    !$acc end kernels

    do itr = 1, NITER_RH
       call SATURATION_psat_all( KA, KS, KE, &
                                 temp(:,1,1), & ! [IN]
                                 psat(:,1,1)  ) ! [OUT]
       !$acc kernels
       do k = KS, KE
          if( CZ(k) <= ENV_L2_ZTOP ) then ! Layer 1 and 2
             qv(k,1,1) = ENV_RH * 1.E-2_RP * psat(k,1,1) / ( dens(k,1,1) * Rvap * temp(k,1,1) )
          else                            ! Layer 3
             qv(k,1,1) = 0.0_RP
          endif
       enddo
       !$acc end kernels

       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  pott(:,1,1), qv(:,1,1), qc(:,1,1),                      & ! [IN]
                                  pres_sfc(1,1), pott_sfc(1,1), qv_sfc(1,1), qc_sfc(1,1), & ! [IN]
                                  CZ(:), FZ(:),                                           & ! [IN]
#ifdef _OPENACC
                                  work1(:), work2(:), work3(:),                           & ! [WORK]
#endif
                                  DENS(:,1,1), temp(:,1,1), pres(:,1,1), temp_sfc(1,1),   & ! [OUT]
                                  converged                                               ) ! [OUT]
    end do

    !$acc kernels
    !$acc loop collapse(3) independent
    do j = JSB, JEB
    do i = ISB, IEB
    do k = KS, KE
       DENS(k,i,j) = DENS(k,1,1)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = ENV_U * DENS(k,i,j)
       MOMY(k,i,j) = ENV_V * DENS(k,i,j)

       ! make warm bubble
       RHOT(k,i,j) = DENS(k,1,1) * ( pott(k,1,1) + BBL_THETA * bubble(k,i,j) )

       qv  (k,i,j) = qv(k,1,1)
    enddo
    enddo
    enddo
    !$acc end kernels

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_warmbubbleaero

  !-----------------------------------------------------------------------------
  !> Make initial state from a sounding file and add perturbation (noise and/or bubble)
  subroutine MKINIT_sonde_perturb
    use scale_atmos_grid_cartesC, only: &
         GRID_CZ => ATMOS_GRID_CARTESC_CZ

    implicit none

    ! 1D working variables to read sonde data
    real(RP) :: RHO1D(KA)
    real(RP) :: VELX1D(KA)
    real(RP) :: VELY1D(KA)
    real(RP) :: POTT1D(KA)
    real(RP) :: QV1D(KA)

    ! Noise
    real(RP) :: PERTURB_AMP_PT  = 0.0_RP
    real(RP) :: PERTURB_AMP_QV =  0.0_RP
    real(RP) :: RANDOM_LIMIT_Z        ! noise is added below this height. Default is the domain height [m]

    ! Bubble
    real(RP) :: BBL_THETA = 0.D0 ! extremum of temperature in bubble [K]

    NAMELIST / PARAM_MKINIT_SONDE_PERTURB / &
       PERTURB_AMP_PT,     &
       PERTURB_AMP_QV,     &
       RANDOM_LIMIT_Z,  &
       BBL_THETA

    integer :: ierr
    integer :: k, i, j
    !---------------------------------------------------------------------------

    if( IO_L ) write(IO_FID_LOG,*)
    if( IO_L ) write(IO_FID_LOG,*) '+++ Module[SONDE_PERTURB]/Categ[INIT]'

    !--- read namelist
    RANDOM_LIMIT_Z = GRID_CZ(KE) ! default value of RANDOM_LIMIT

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SONDE_PERTURB,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       if( IO_L ) write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       write(*,*) 'xxx Not appropriate names in namelist PARAM_MKINIT_SONDE_PERTURB. Check!'
       call PRC_abort
    endif
    if( IO_NML ) write(IO_FID_LOG,nml=PARAM_MKINIT_SONDE_PERTURB)

    !!!!!! read sonde data
    call MKINIT_common_read_sounding( RHO1D, VELX1D, VELY1D, POTT1D, QV1D ) ! (out)

    do j = 1, JA
    do i = 1, IA
    do k = KS, KE
       DENS(k,i,j) = RHO1D(k)
       MOMZ(k,i,j) = 0.0_RP
       MOMX(k,i,j) = RHO1D(k) * VELX1D(k)
       MOMY(k,i,j) = RHO1D(k) * VELY1D(k)
       RHOT(k,i,j) = RHO1D(k) * POTT1D(k)
       QV  (k,i,j) = QV1D(k)
    enddo
    enddo
    enddo

    !!!!!! add perturbation to PT
    call RANDOM_uniform(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if( GRID_CZ(k) <= RANDOM_LIMIT_Z ) then
          RHOT(k,i,j) = RHOT(k,i,j) + (2.0_RP*(rndm(k,i,j)-0.50_RP))*PERTURB_AMP_PT*DENS(k,i,j)
       endif
    enddo
    enddo
    enddo

    !!!!!! add perturbation to QV
    call RANDOM_uniform(rndm) ! make random
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       if( GRID_CZ(k) <= RANDOM_LIMIT_Z ) then
          QV(k,i,j) = QV(k,i,j) + (2.0_RP*(rndm(k,i,j)-0.50_RP))*PERTURB_AMP_QV
       endif
    enddo
    enddo
    enddo

    !!!!!! flux setup (not sure how this works)
    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_sonde_perturb
  !-----------------------------------------------------------------------------
  !> Make initial state ( real case )
  subroutine MKINIT_real
    use mod_realinput, only: &
         REALINPUT_atmos, &
         REALINPUT_surface
    implicit none

    call PROF_rapstart('__Real_Atmos',2)

    call REALINPUT_atmos

    call PROF_rapend  ('__Real_Atmos',2)
    call PROF_rapstart('__Real_Surface',2)

    call REALINPUT_surface

    call PROF_rapend  ('__Real_Surface',2)

    call MKINIT_common_flux_setup

    return
  end subroutine MKINIT_real

end module mod_mkinit
