!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics
!!
!! @par Description
!!          HEVE DGM scheme for large time step in Atmospheric dynamical process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_large_dgm_heve
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_index
  use scale_tracer
  use scale_const, only: &
    Rdry => CONST_Rdry,   &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
    
  use scale_fem_element_base, only: &
    ElementBase3D
  use scale_fem_localmesh_3d, only: &
    LocalMesh3D
  use scale_fem_localmeshfield_base, only: &
    LocalMeshFieldBaseList
  use scale_fem_meshfield_base, only: &
    MeshField3D, MeshField2D

  use scale_atmos_dyn_dgm_operator, only: &
    mesh3D, &
    dyn_bnd
  use scale_atmos_dyn_dgm_vars, only: &
    AtmosDynDGM_vars
#ifdef DEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_large_dgm_heve_setup
  public :: ATMOS_DYN_Tstep_large_dgm_heve_regist_dynvar
  public :: ATMOS_DYN_Tstep_large_dgm_heve_finalize
  public :: ATMOS_DYN_Tstep_large_dgm_heve

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: set_FVvars
  private :: interp_FVtoDG_DAMP_vars
  private :: add_nudge_forcing

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

  ! for communication
  integer :: I_COMM_DENS
  integer :: I_COMM_MOMX
  integer :: I_COMM_MOMY
  integer :: I_COMM_MOMZ
  integer :: I_COMM_RHOT
  integer :: I_COMM_MOMX_c
  integer :: I_COMM_MOMY_c
  real(RP), allocatable :: MOMX_c(:,:,:)
  real(RP), allocatable :: MOMY_c(:,:,:)

  integer, allocatable :: I_COMM_RHOQ_t(:)
  integer, allocatable :: I_COMM_QTRC(:)

  type(AtmosDynDGM_vars), pointer :: dyn_vars

  ! For boundary grid nudging
  integer, public, parameter :: BGN_TEND_DENS_ID     = 1
  integer, public, parameter :: BGN_TEND_MOMX_ID     = 2
  integer, public, parameter :: BGN_TEND_MOMY_ID     = 3
  integer, public, parameter :: BGN_TEND_MOMZ_ID     = 4
  integer, public, parameter :: BGN_TEND_RHOT_ID     = 5
  integer, public, parameter :: BGN_TEND_NUM         = 5

  type(MeshField3D) :: bgn_tend_dg(BGN_TEND_NUM)
  integer :: HIST_bgn_tend_dg(BGN_TEND_NUM)
  integer, allocatable :: HIST_damp_QTRC(:)

  integer :: GRID_NUDGING_FV2DG_INTERP_ORDER = 1
  type(MeshField3D) :: DAMP_DENS_dg
  type(MeshField3D) :: DAMP_RHOT_dg

  !-----------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------------
  !> Setup
!OCL SERIAL
  subroutine ATMOS_DYN_Tstep_large_dgm_heve_setup( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG )
    use scale_prc, only: &
       PRC_abort
    use scale_prc_cartesC, only: &
       PRC_TwoD,  &
       PRC_HAS_E, &
       PRC_HAS_W, &
       PRC_HAS_N, &
       PRC_HAS_S
    use scale_const, only: &
       OHM => CONST_OHM, &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    use scale_file_history, only: &
       FILE_HISTORY_reg, &
       FILE_HISTORY_put
    use scale_file_fem_history, only: &
       FILE_HISTORY_meshfield_put
    use scale_monitor, only: &
       MONITOR_reg, &
       MONITOR_put
    use scale_fem_variableinfo, only: &
      VariableInfo

    implicit none

    ! MPI_RECV_INIT requires intent(inout)
    real(RP),               intent(inout) :: DENS(KA,IA,JA)
    real(RP),               intent(inout) :: MOMZ(KA,IA,JA)
    real(RP),               intent(inout) :: MOMX(KA,IA,JA)
    real(RP),               intent(inout) :: MOMY(KA,IA,JA)
    real(RP),               intent(inout) :: RHOT(KA,IA,JA)
    real(RP),               intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP),               intent(inout) :: PROG(KA,IA,JA,VA)

    integer :: iv, iq

    type(VariableInfo) :: BGN_TEND_VARINFO(BGN_TEND_NUM)
    DATA BGN_TEND_VARINFO / &
      VariableInfo( BGN_TEND_DENS_ID, 'DG_DENS_BGN_t', 'tendency of boundary grid nudging for DENS', &
                    'kg/m3/s',  3, 'FEM_XYZ',  '' ), &
      VariableInfo( BGN_TEND_MOMX_ID, 'DG_MOMX_BGN_t', 'tendency of boundary grid nudging for MOMX', &
                    'kg/m3.m/s.s-1',  3, 'FEM_XYZ',  '' ), &
      VariableInfo( BGN_TEND_MOMY_ID, 'DG_MOMY_BGN_t', 'tendency of boundary grid nudging for MOMY', &
                    'kg/m3.m/s.s-1',  3, 'FEM_XYZ',  '' ), &
      VariableInfo( BGN_TEND_MOMZ_ID, 'DG_MOMZ_BGN_t', 'tendency of boundary grid nudging for MOMZ', &
                    'kg/m3.m/s.s-1',  3, 'FEM_XYZ',  '' ), &
      VariableInfo( BGN_TEND_RHOT_ID, 'DG_RHOT_BGN_t', 'tendency of boundary grid nudging for RHOT', &
                    'K.s-1',  3, 'FEM_XYZ',  '' ) /


    namelist /ATMOS_DYN_TSTEP_LARGE_DGM_HEVE/ &
      GRID_NUDGING_FV2DG_INTERP_ORDER

    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=ATMOS_DYN_TSTEP_LARGE_DGM_HEVE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_DYN_Tstep_large_dgm_heve_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_DYN_Tstep_large_dgm_heve_setup",*) 'Not appropriate names in namelist ATMOS_DYN_TSTEP_LARGE_DGM_HEVE. Check!'
       call PRC_abort
    endif
    LOG_NML(ATMOS_DYN_TSTEP_LARGE_DGM_HEVE)

    I_COMM_DENS = 1
    I_COMM_MOMZ = 2
    I_COMM_MOMX = 3
    I_COMM_MOMY = 4
    I_COMM_RHOT = 5
    I_COMM_MOMX_c = 6
    I_COMM_MOMY_c = 7

    allocate( MOMX_c(KA,IA,JA), MOMY_c(KA,IA,JA) )
    call COMM_vars8_init( 'DENS', DENS, I_COMM_DENS )
    call COMM_vars8_init( 'MOMZ', MOMY, I_COMM_MOMZ )
    call COMM_vars8_init( 'MOMX', MOMX, I_COMM_MOMX )
    call COMM_vars8_init( 'MOMY', MOMY, I_COMM_MOMY )
    call COMM_vars8_init( 'RHOT', RHOT, I_COMM_RHOT )
    call COMM_vars8_init( 'MOMX_c', MOMX_c, I_COMM_MOMX_c )
    call COMM_vars8_init( 'MOMY_c', MOMY_c, I_COMM_MOMY_c )

    allocate( I_COMM_QTRC(QA) )
    do iq = 1, QA
      I_COMM_QTRC(iq) = 5 + VA + iq
      call COMM_vars8_init( 'QTRC',   QTRC  (:,:,:,iq), I_COMM_QTRC(iq) )
    end do

    ! Initialize variables for boundary grid nudging

    do iv = 1, BGN_TEND_NUM
      call bgn_tend_dg(iv)%Init( BGN_TEND_VARINFO(iv)%NAME, BGN_TEND_VARINFO(iv)%UNIT, mesh3D, fill_value=0.0_RP )

      call FILE_HISTORY_reg( BGN_TEND_VARINFO(iv)%NAME, BGN_TEND_VARINFO(iv)%DESC, BGN_TEND_VARINFO(iv)%unit, HIST_bgn_tend_dg(iv), &
        dim_type=BGN_TEND_VARINFO(iv)%dim_type )    
      call FILE_HISTORY_meshfield_put( HIST_bgn_tend_dg(iv), bgn_tend_dg(iv) )  
    end do

    call DAMP_DENS_dg%Init( "DAMP_DENS", "kg/m3", mesh3D )
    call DAMP_RHOT_dg%Init( "DAMP_DENS", "kg/m3", mesh3D )

    ! Put zero data as physics tendencies on DG nodes at t=0
    call dyn_vars%Put_hist_phytend()

    return
  end subroutine ATMOS_DYN_Tstep_large_dgm_heve_setup

  !-----------------------------------------------------------------------------
  !> finalize
!OCL SERIAL
  subroutine ATMOS_DYN_Tstep_large_dgm_heve_finalize
    implicit none
    integer :: iv
    !-----------------------------

    deallocate( MOMX_c, MOMY_c )
    do iv=1, BGN_TEND_NUM
      call bgn_tend_dg(iv)%Final()
    end do
    call DAMP_DENS_dg%Final()
    call DAMP_RHOT_dg%Final()

    deallocate( I_COMM_QTRC )

    return
  end subroutine ATMOS_DYN_Tstep_large_dgm_heve_finalize

  !-----------------------------------------------------------------------------
  !> Regist object for managing DG variables
!OCL SERIAL
  subroutine ATMOS_DYN_Tstep_large_dgm_heve_regist_dynvar( dyn_dg_vars )
    implicit none
    type(AtmosDynDGM_vars), intent(inout), target :: dyn_dg_vars
    !-----------------------------

    dyn_vars => dyn_dg_vars
    return
  end subroutine ATMOS_DYN_Tstep_large_dgm_heve_regist_dynvar

  !-----------------------------------------------------------------------------
  !> Dynamical Process
!OCL SERIAL
  subroutine ATMOS_DYN_Tstep_large_dgm_heve( &
       DENS, MOMZ, MOMX, MOMY, RHOT, QTRC, PROG,             &
       DENS_av, MOMZ_av, MOMX_av, MOMY_av, RHOT_av, QTRC_av, &
       num_diff, num_diff_q,                                 &
       QTRC0,                                                &
       DENS_tp, MOMZ_tp, MOMX_tp, MOMY_tp, RHOT_tp, RHOQ_tp, &
       CORIOLI,                                              &
       CDZ, CDX, CDY, FDZ, FDX, FDY,                         &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,                   &
       PHI, GSQRT,                                           &
       J13G, J23G, J33G, MAPF,                               &
       AQ_R, AQ_CV, AQ_CP, AQ_MASS,                          &
       REF_dens, REF_pott, REF_qv, REF_pres,                 &
       BND_W, BND_E, BND_S, BND_N, TwoD,                     &
       ND_COEF, ND_COEF_Q, ND_LAPLACIAN_NUM,                 &
       ND_SFC_FACT, ND_USE_RS,                               &
       BND_QA, BND_IQ, BND_SMOOTHER_FACT,                    &
       DAMP_DENS,       DAMP_VELZ,       DAMP_VELX,          &
       DAMP_VELY,       DAMP_POTT,       DAMP_QTRC,          &
       DAMP_alpha_DENS, DAMP_alpha_VELZ, DAMP_alpha_VELX,    &
       DAMP_alpha_VELY, DAMP_alpha_POTT, DAMP_alpha_QTRC,    &
       MFLUX_OFFSET_X, MFLUX_OFFSET_Y,                       &
       wdamp_coef, divdmp_coef,                              &
       FLAG_TRACER_SPLIT_TEND,                               &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T, FLAG_FCT_TRACER,       &
       FLAG_FCT_ALONG_STREAM,                                &
       USE_AVERAGE,                                          &
       I_QV,                                                 &
       DTLS, DTSS, Llast                                     )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       EPS    => CONST_EPS
    use scale_file_history, only: &
      FILE_HISTORY_query, &
      FILE_HISTORY_put,   &
      FILE_HISTORY_set_disable
    use scale_file_fem_history, only: &
      FILE_HISTORY_meshfield_put      
    use scale_atmos_dyn_tinteg_short_rkdg, only: &
      ATMOS_DYN_Tinteg_short_rkdg
    use scale_atmos_dyn_tinteg_tracer_rkdg, only: &
      ATMOS_DYN_Tinteg_tracer_short_rkdg
    use scale_atmos_dyn_tstep_tracer_dgm_heve, only: &
      ATMOS_DYN_Tstep_tracer_dgm_heve_cal_alphdens_advtest
    use scale_atmos_dyn_dgm_operator, only: &
      set_fv_phytend => ATMOS_DYN_DGM_operator_set_fv_phytend, &
      mesh3D, elem3D,   &
      ONLY_TRCADV_FLAG, TRCADV_disable_limiter
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

    real(RP), intent(inout) :: DENS_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_av(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_av(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_av(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_av(KA,IA,JA,QA)

    real(RP), intent(out)   :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(out)   :: num_diff_q(KA,IA,JA,3)

    real(RP), intent(in)    :: QTRC0(KA,IA,JA,QA)

    real(RP), intent(in)    :: DENS_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_tp(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_tp(KA,IA,JA)
    real(RP), intent(in)    :: RHOQ_tp(KA,IA,JA,QA)

    real(RP), intent(in)    :: CORIOLI(IA,JA)

    real(RP), intent(in)    :: CDZ (KA)
    real(RP), intent(in)    :: CDX (IA)
    real(RP), intent(in)    :: CDY (JA)
    real(RP), intent(in)    :: FDZ (KA-1)
    real(RP), intent(in)    :: FDX (IA-1)
    real(RP), intent(in)    :: FDY (JA-1)
    real(RP), intent(in)    :: RCDZ(KA)
    real(RP), intent(in)    :: RCDX(IA)
    real(RP), intent(in)    :: RCDY(JA)
    real(RP), intent(in)    :: RFDZ(KA-1)
    real(RP), intent(in)    :: RFDX(IA-1)
    real(RP), intent(in)    :: RFDY(JA-1)

    real(RP), intent(in)    :: PHI  (KA,IA,JA)   !< geopotential
    real(RP), intent(in)    :: GSQRT(KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)    :: J13G (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)    :: J23G (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)    :: J33G              !< (3,3) element of Jacobian matrix
    real(RP), intent(in)    :: MAPF (IA,JA,2,4)  !< map factor

    real(RP), intent(in)    :: AQ_R   (QA)
    real(RP), intent(in)    :: AQ_CV  (QA)
    real(RP), intent(in)    :: AQ_CP  (QA)
    real(RP), intent(in)    :: AQ_MASS(QA)

    real(RP), intent(in)    :: REF_dens(KA,IA,JA)
    real(RP), intent(in)    :: REF_pott(KA,IA,JA)
    real(RP), intent(in)    :: REF_qv  (KA,IA,JA)
    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure

    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical,  intent(in)    :: TwoD

    real(RP), intent(in)    :: ND_COEF
    real(RP), intent(in)    :: ND_COEF_Q
    integer,  intent(in)    :: ND_LAPLACIAN_NUM
    real(RP), intent(in)    :: ND_SFC_FACT
    logical,  intent(in)    :: ND_USE_RS

    integer,  intent(in)    :: BND_QA
    integer,  intent(in)    :: BND_IQ(QA)
    real(RP), intent(in)    :: BND_SMOOTHER_FACT

    real(RP), intent(in)    :: DAMP_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_QTRC(KA,IA,JA,BND_QA)

    real(RP), intent(in)    :: DAMP_alpha_DENS(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELZ(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELX(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_VELY(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_POTT(KA,IA,JA)
    real(RP), intent(in)    :: DAMP_alpha_QTRC(KA,IA,JA,BND_QA)
    real(RP), intent(in)    :: MFLUX_OFFSET_X(KA,JA,2)
    real(RP), intent(in)    :: MFLUX_OFFSET_Y(KA,IA,2)

    real(RP), intent(in)    :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef

    logical,  intent(in)    :: FLAG_TRACER_SPLIT_TEND
    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_TRACER
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    logical,  intent(in)    :: USE_AVERAGE

    integer,  intent(in)    :: I_QV

    real(DP), intent(in)    :: DTLS
    real(DP), intent(in)    :: DTSS
    logical , intent(in)    :: Llast

    real(RP) :: dtl
    real(RP) :: dts
    integer  :: nstep

    integer :: iq, step
    integer :: iv

    integer :: lcdomid
    class(LocalMesh3D), pointer :: lcmesh3D
    type(MeshField3D), pointer :: DDENS_dg, MOMX_dg, MOMY_dg, MOMZ_dg, DRHOT_dg 
    type(MeshField3D), pointer :: DENS_hyd_dg, PRES_hyd_dg
    type(MeshField3D), pointer :: DENS_tp_dg, MOMX_tp_dg, MOMY_tp_dg, MOMZ_tp_dg, RHOT_tp_dg, RHOH_tp_dg

    type(MeshField3D), pointer :: QTRC_dg, RHOQ_tp_dg
    type(MeshField3D), pointer :: QTRC_tmp_dg, FCT_coef_dg, DDENS_TRC_dg, DDENS0_TRC_dg
    type(MeshField3D), pointer :: MFLX_x_tavg_dg, MFLX_y_tavg_dg, MFLX_z_tavg_dg
    type(MeshField3D), pointer :: alphDensM_dg, alphDensP_dg

    type(LocalMeshFieldBaseList) :: lc_QTRC_dg(QA)
    type(LocalMeshFieldBaseList) :: lc_RHOQ_tp_dg(QA)

    integer :: kelem
    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_Large_Preparation", 2)

    dts   = real(DTSS, kind=RP)            ! short time step
    dtl   = real(DTLS, kind=RP)            ! large time step
    nstep = ceiling( ( dtl - eps ) / dts )
    dts   = dtl / nstep                    ! dts is divisor of dtl and smaller or equal to dtss
    if ( ONLY_TRCADV_FLAG ) nstep = -1

    !
    call dyn_vars%GetVarsPtr( DDENS_dg, MOMX_dg, MOMY_dg, MOMZ_dg, DRHOT_dg, &
       DENS_hyd_dg, PRES_hyd_dg )
  
    call dyn_vars%GetPhyTendPtr( DENS_tp_dg, MOMX_tp_dg, MOMY_tp_dg, MOMZ_tp_dg, RHOT_tp_dg, RHOH_tp_dg )
    
    lcdomID=1
    lcmesh3D => mesh3D%lcmesh_list(lcdomID)

    do iq=1, QA
      call dyn_vars%GetTRCVarsPtr( iq, QTRC_dg, RHOQ_tp_dg )
      lc_RHOQ_tp_dg(iq)%ptr => RHOQ_tp_dg%local(lcdomid)
      lc_QTRC_dg(iq)%ptr => QTRC_dg%local(lcdomid)
    end do

    ! Update specific heat
    if ( QA > 0 ) call dyn_vars%Calculate_specific_heat()
    ! Exchange halo data of auxiliary variables
    call dyn_vars%Exchange_auxvars() 

    call PROF_rapend  ("DYN_Large_Preparation", 2)

    !###########################################################################
    ! Update DENS,MONZ,MOMX,MOMY,MOMZ,RHOT
    !###########################################################################

    call PROF_rapstart("DYN_Large_Tendency", 2)
    call PROF_rapend  ("DYN_Large_Tendency", 2)

    call PROF_rapstart("DYN_Large_Boundary", 2)

    lcdomID=1
    lcmesh3D => mesh3D%lcmesh_list(lcdomID)  

    ! Store FV boundary data for boundary grid nuding
    call interp_FVtoDG_DAMP_vars( DAMP_DENS_dg%local(lcdomid)%val, DAMP_RHOT_dg%local(lcdomid)%val, & ! (out)
      DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_VELZ, DAMP_POTT, DAMP_QTRC,  & ! (in)
      BND_QA, lcmesh3D, elem3D )

    ! Set boundary data of computational domain
    call dyn_bnd%Store_nuding_FV_bnddata( &
      DAMP_DENS_dg%local(lcdomid)%val, DAMP_VELX, DAMP_VELY, DAMP_VELZ, DAMP_RHOT_dg%local(lcdomid)%val, DAMP_QTRC, &
      MAPF, DAMP_alpha_VELZ, BND_QA, BND_IQ, lcmesh3D, lcmesh3D%refElem3D )

    ! Calculate the tendencies with boundary grid nudging and add it to DG forcing
    call add_nudge_forcing( &
      DENS_tp_dg%local(lcdomid)%val, MOMX_tp_dg%local(lcdomid)%val, MOMY_tp_dg%local(lcdomid)%val,        & ! (inout)
      MOMZ_tp_dg%local(lcdomid)%val, RHOT_tp_dg%local(lcdomid)%val, lc_RHOQ_tp_dg,                        & ! (inout)
      bgn_tend_dg(BGN_TEND_DENS_ID)%local(lcdomid)%val, bgn_tend_dg(BGN_TEND_MOMX_ID)%local(lcdomid)%val, & ! (inout)
      bgn_tend_dg(BGN_TEND_MOMY_ID)%local(lcdomid)%val, bgn_tend_dg(BGN_TEND_MOMZ_ID)%local(lcdomid)%val, & ! (inout)
      bgn_tend_dg(BGN_TEND_RHOT_ID)%local(lcdomid)%val,                                            & ! (in)
      DENS, MOMX, MOMY, MOMZ, RHOT, QTRC,                                                          & ! (in)
      DAMP_DENS, DAMP_VELX, DAMP_VELY, DAMP_VELZ, DAMP_POTT, DAMP_QTRC,                            & ! (in)
      DAMP_DENS_dg%local(lcdomid)%val, DAMP_RHOT_dg%local(lcdomid)%val,                            & ! (in)
      DAMP_alpha_DENS, DAMP_alpha_VELX, DAMP_alpha_VELY, DAMP_alpha_VELZ, DAMP_alpha_POTT, DAMP_alpha_QTRC, & ! (in)
      BND_QA, BND_IQ, BND_SMOOTHER_FACT, DENS,                                                     & ! (in)
      DDENS_dg%local(lcdomid)%val, MOMX_dg%local(lcdomid)%val, MOMY_dg%local(lcdomid)%val,         & ! (in)
      MOMZ_dg%local(lcdomid)%val, DRHOT_dg%local(lcdomid)%val,                                     & ! (in)
      lc_QTRC_dg, DENS_hyd_dg%local(lcdomid)%val, PRES_hyd_dg%local(lcdomid)%val, lcmesh3D,        & ! (in) 
      lcmesh3D%refElem3D )                     ! (in)

    call PROF_rapend  ("DYN_Large_Boundary", 2)    
    
    do step = 1, nstep
      !-----< prepare tendency >-----

      ! call PROF_rapstart("DYN_Large_Tendency", 2)
      ! call PROF_rapend  ("DYN_Large_Tendency", 2)
      
      call FILE_HISTORY_set_disable( .not. ( Llast .and. ( step==nstep ) ) )

      !- Start short time integration -------------------
      call PROF_rapstart("DYN_Short_Tinteg", 2)
      call ATMOS_DYN_Tinteg_short_rkdg( dyn_vars, &
        step, nstep, Llast )
      call PROF_rapend  ("DYN_Short_Tinteg", 2)

      call FILE_HISTORY_set_disable( .false. )
   enddo ! dynamical steps


   !###########################################################################
   ! Update Tracers
   !###########################################################################
   
    if ( QA > 0 ) then
      call PROF_rapstart( 'ATM_DYN_trc_update_pre', 2)

      !
      call dyn_vars%GetTRCAuxVarsPtr( QTRC_tmp_dg, FCT_coef_dg, DDENS_TRC_dg, DDENS0_TRC_dg, DENS_hyd_dg, &
        MFLX_x_tavg_dg, MFLX_y_tavg_dg, MFLX_z_tavg_dg, alphDensM_dg, alphDensP_dg )
      
      if ( ONLY_TRCADV_FLAG ) then
        call dyn_vars%Exchange_prgvars()
        
        do lcdomid=1, mesh3D%LOCAL_MESH_NUM
          lcmesh3D => mesh3D%lcmesh_list(lcdomID)
          !$omp parallel do
          do kelem=lcmesh3D%NeS, lcmesh3D%NeE
            DDENS0_TRC_dg %local(lcdomid)%val(:,kelem) = DDENS_dg%local(lcdomid)%val(:,kelem) 
            DDENS_TRC_dg  %local(lcdomid)%val(:,kelem) = DDENS_dg%local(lcdomid)%val(:,kelem)
            MFLX_z_tavg_dg%local(lcdomid)%val(:,kelem) = MOMZ_dg%local(lcdomid)%val(:,kelem)
            MFLX_x_tavg_dg%local(lcdomid)%val(:,kelem) = MOMX_dg%local(lcdomid)%val(:,kelem)
            MFLX_y_tavg_dg%local(lcdomid)%val(:,kelem) = MOMY_dg%local(lcdomid)%val(:,kelem)
          enddo
          call ATMOS_DYN_Tstep_tracer_dgm_heve_cal_alphdens_advtest( &
            alphDensM_dg%local(lcdomid)%face_val, alphDensP_dg%local(lcdomid)%face_val,                                      & ! (inout)
            DDENS_dg%local(lcdomid)%val, MOMX_dg%local(lcdomid)%val, MOMY_dg%local(lcdomid)%val, MOMZ_dg%local(lcdomid)%val, & ! (in)
            DENS_hyd_dg%local(lcdomid)%val,                                                           & ! (in)
            lcmesh3D%Gsqrt,                                                                  & ! (in)
            lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3), & ! (in)
            lcmesh3D%VMapM, lcmesh3D%VMapP, lcmesh3D, lcmesh3D%refElem3D                     ) ! (in)
        enddo
      end if

      ! Exchange halo data of mass flux
      call PROF_rapstart( 'ATM_DYN_exchange_mflx', 3)
      call dyn_vars%Exchange_massflux()
      call PROF_rapend( 'ATM_DYN_exchange_mflx', 3)
  
      call PROF_rapstart( 'ATM_DYN_applyBC_mflux', 3)
      do lcdomid=1, mesh3D%LOCAL_MESH_NUM
        lcmesh3D => mesh3D%lcmesh_list(lcdomID)
        call dyn_bnd%ApplyBC_PROGVARS_lc( lcdomid, & ! (in)
          DDENS_TRC_dg%local(lcdomID)%val(:,:),                                                                    & ! (inout)
          MFLX_x_tavg_dg%local(lcdomID)%val, MFLX_y_tavg_dg%local(lcdomID)%val, MFLX_z_tavg_dg%local(lcdomID)%val, & ! (inout)
          DRHOT_dg%local(lcdomID)%val,                                                                             & ! (inout)
          DENS_hyd_dg%local(lcdomID)%val, PRES_hyd_dg%local(lcdomID)%val,                                                 & ! (in)
          lcmesh3D%Gsqrt(:,:), lcmesh3D%GsqrtH(:,:), lcmesh3D%GIJ(:,:,1,1), lcmesh3D%GIJ(:,:,1,2), lcmesh3D%GIJ(:,:,2,2), & ! (in) 
          lcmesh3D%GI3(:,:,1), lcmesh3D%GI3(:,:,2),                                           & ! (in)
          lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3),    & ! (in)
          lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB,                                     & ! (in)
          lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D        ) ! (in)
      end do
      call PROF_rapend( 'ATM_DYN_applyBC_mflux', 3) 

      call PROF_rapend( 'ATM_DYN_trc_update_pre', 2)
    end if
    
    do iq = 1, QA
      !------------------------------------------------------------------------
      ! Update each tracer
      !------------------------------------------------------------------------

      call dyn_vars%GetTRCVarsPtr( iq, QTRC_dg, RHOQ_tp_dg )

      lcdomID=1
      lcmesh3D => mesh3D%lcmesh_list(lcdomID)
      call set_fv_phytend( RHOQ_tp_dg%local(lcdomid)%val, &
        RHOQ_tp(:,:,:,iq), lcmesh3D, elem3D )

      call PROF_rapstart("DYN_Tracer_Tinteg", 2)
      call ATMOS_DYN_Tinteg_tracer_short_rkdg( &
        QTRC_dg, QTRC_tmp_dg, FCT_coef_dg, iq,              &
        MFLX_x_tavg_dg, MFLX_y_tavg_dg, MFLX_z_tavg_dg,     &
        RHOQ_tp_dg, alphDensM_dg, alphDensP_dg,             &
        DENS_hyd_dg, DDENS_dg, DDENS_TRC_dg, DDENS0_TRC_dg, &
        dyn_vars, BND_QA, BND_IQ, Llast )
      call PROF_rapend  ("DYN_Tracer_Tinteg", 2)
    enddo

    call PROF_rapstart( 'ATM_DYN_trc_update_post', 2)

    ! Update specific heat and pressure
    call dyn_vars%Calculate_specific_heat()
    call dyn_vars%Calcuate_pressure()

    if ( QA > 0 ) then
      ! Negative fixer
      if (      .not. TRCADV_disable_limiter &
          .and. (.not. ONLY_TRCADV_FLAG )    ) then
        call dyn_vars%Fix_negative_value_qtrc()
      end if
    end if

    call PROF_rapend( 'ATM_DYN_trc_update_post', 2)

    !---

    if ( Llast ) then
      call PROF_rapstart("ATM_DYN_DG2FVvar", 2)
      lcdomID=1
      lcmesh3D => mesh3D%lcmesh_list(lcdomID)
      call set_FVvars( DENS, MOMX, MOMY, MOMZ, RHOT, QTRC, &
        DDENS_dg%local(lcdomID)%val, MOMX_dg%local(lcdomID)%val, MOMY_dg%local(lcdomID)%val, &
        MOMZ_dg%local(lcdomID)%val, DRHOT_dg%local(lcdomID)%val, lc_QTRC_dg,                 &
        DENS_hyd_dg%local(lcdomID)%val, PRES_hyd_dg%local(lcdomID)%val,                      &
        lcmesh3D, lcmesh3D%refElem3D, BND_W, BND_E, BND_S, BND_N                             )
      call PROF_rapend("ATM_DYN_DG2FVvar", 2)
      
      !- history ---------------------------------------

      call dyn_vars%Put_hist_phytend()

      do iv=1, BGN_TEND_NUM
        call FILE_HISTORY_meshfield_put( HIST_bgn_tend_dg(iv), bgn_tend_dg(iv) )
      end do
    endif

    return
  end subroutine ATMOS_DYN_Tstep_large_dgm_heve

  !-- private subroutines --------------------------------------------

  subroutine set_FVvars( DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT, QTRC, DENS_hyd, PRES_hyd, &
    lcmesh3D, elem3D, BND_W, BND_E, BND_S, BND_N )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: DENS_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_fv(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_fv(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_fv(KA,IA,JA,QA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(QA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh3D%NeA)
    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N

    integer :: ke_x, ke_y, ke_z
    integer :: i, j, k
    integer :: iq
    integer :: kelem, kelem2D
    real(RP) :: IntWeight(elem3D%Np)
    real(RP) :: vol(lcmesh3D%Ne)
    real(RP) :: RHOT_hyd(elem3D%Np)
    real(RP) :: r_dens
    real(RP) :: MOMZ_c(KA,IA,JA)

    integer :: IS_, IE_, JS_, JE_
    !--------------------------------------------------------------------

    !$omp parallel private( &
    !$omp ke_x, ke_y, ke_z, i, j, k, kelem, kelem2D, iq, IntWeight, RHOT_hyd ) 
    
    !$omp do collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      kelem2D = lcmesh3D%EMap3Dto2D(kelem)
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      IntWeight(:) = elem3D%IntWeight_lgl(:) * lcmesh3D%J(:,kelem) * lcmesh3D%Gsqrt(:,kelem)
      vol(kelem) = sum(IntWeight(:))
      IntWeight(:) = IntWeight(:) / vol(kelem)

      DENS_fv(k,i,j) = sum( IntWeight(:) * ( DENS_hyd(:,kelem) + DDENS(:,kelem) ) )
      MOMX_c(k,i,j) = sum( IntWeight(:) * MOMX(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,1,1)) )
      MOMY_c(k,i,j) = sum( IntWeight(:) * MOMY(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,2,2)) )
      MOMZ_c(k,i,j) = sum( IntWeight(:) * MOMZ(:,kelem) )

      RHOT_hyd(:) = PRES00 / Rdry * ( PRES_hyd(:,kelem) / PRES00 )**(CvDry/CpDry)
      RHOT_fv(k,i,j) = sum( IntWeight(:) * ( RHOT_hyd(:) + DRHOT(:,kelem) ) )
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do collapse(3)
    do iq=1, QA    
      do ke_z=1, lcmesh3D%NeZ
      do ke_y=1, lcmesh3D%NeY
      do ke_x=1, lcmesh3D%NeX
        kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
        i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

        IntWeight(:) = elem3D%IntWeight_lgl(:) * lcmesh3D%J(:,kelem) * lcmesh3D%Gsqrt(:,kelem) &
                     / vol(kelem)
        QTRC_fv(k,i,j,iq) = sum( IntWeight(:) * ( DENS_hyd(:,kelem) + DDENS(:,kelem) ) * QTRC(iq)%ptr%val(:,kelem) ) &
                          / DENS_fv(k,i,j)
      enddo
      enddo
      enddo
    enddo
    !$omp end do
    !$omp end parallel

    call COMM_vars8( DENS_fv(:,:,:), I_COMM_DENS )
    call COMM_vars8( MOMX_c(:,:,:), I_COMM_MOMX_c )
    call COMM_vars8( MOMY_c(:,:,:), I_COMM_MOMY_c )
    call COMM_vars8( RHOT_fv(:,:,:), I_COMM_RHOT )

    call COMM_wait ( DENS_fv(:,:,:), I_COMM_DENS, .false. )
    call COMM_wait ( MOMX_c(:,:,:), I_COMM_MOMX_c, .false. )
    call COMM_wait ( MOMY_c(:,:,:), I_COMM_MOMY_c, .false. )
    call COMM_wait ( RHOT_fv(:,:,:), I_COMM_RHOT, .false. )

    do iq=1, QA
      call COMM_vars8( QTRC_fv(:,:,:,iq), I_COMM_QTRC(iq) )
    enddo

    if ( BND_W ) then
      IS_ = IS
    else
      IS_ = IS-1
    end if
    if ( BND_E ) then
      IE_ = IE-1
    else
      IE_ = IE
    end if
    if ( BND_S ) then
      JS_ = JS
    else
      JS_ = JS-1
    end if
    if ( BND_N ) then
      JE_ = JE-1
    else
      JE_ = JE
    end if

    !$omp parallel private(i,j,k)
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS_, IE_
    do k=KS, KE
      MOMX_fv(k,i,j) = 0.5_RP * ( MOMX_c(k,i,j) + MOMX_c(k,i+1,j) )
    enddo
    enddo
    enddo

    !$omp do collapse(2)
    do j=JS_, JE_
    do i=IS, IE
    do k=KS, KE
      MOMY_fv(k,i,j) = 0.5_RP * ( MOMY_c(k,i,j) + MOMY_c(k,i,j+1) )
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS, IE
      do k=KS, KE-1
        MOMZ_fv(k,i,j) = 0.5_RP * ( MOMZ_c(k,i,j) + MOMZ_c(k+1,i,j) )
      enddo
      MOMZ_fv(1:KS-1,i,j) = 0.0_RP; MOMZ_fv(KE:KA,i,j) = 0.0_RP
    enddo
    enddo 
    !$omp end parallel

    do iq = 1, QA
      call COMM_wait ( QTRC_fv(:,:,:,iq), I_COMM_QTRC(iq), .false. )
    enddo
    
    call COMM_vars8( MOMZ_fv(:,:,:), I_COMM_MOMZ )
    call COMM_vars8( MOMX_fv(:,:,:), I_COMM_MOMX )
    call COMM_vars8( MOMY_fv(:,:,:), I_COMM_MOMY )

    call COMM_wait ( MOMZ_fv(:,:,:), I_COMM_MOMZ, .false. )
    call COMM_wait ( MOMX_fv(:,:,:), I_COMM_MOMX, .false. )
    call COMM_wait ( MOMY_fv(:,:,:), I_COMM_MOMY, .false. )

    return
  end subroutine set_FVvars

!OCL SERIAL
  subroutine interp_FVtoDG_DAMP_vars( DAMP_DENS_dg_, DAMP_RHOT_dg_, & ! (out)
    DAMP_DENS_fv, DAMP_VELX_fv, DAMP_VELY_fv, DAMP_VELZ_fv, DAMP_POTT_fv, DAMP_QTRC_fv,  & ! (in)
    BND_QA, lcmesh3D, elem3D )
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO, &
      I_UY, I_XV
    use scale_fem_meshfield_util, only: &
      MeshFieldUtil_interp_FVtoDG_3D
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    integer, intent(in) :: BND_QA
    real(RP), intent(inout) :: DAMP_DENS_dg_(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: DAMP_RHOT_dg_(elem3D%Np,lcmesh3D%NeA)    
    real(RP), intent(in) :: DAMP_DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELX_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELY_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELZ_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_POTT_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_QTRC_fv(KA,IA,JA,BND_QA)   
    
    real(RP) :: DAMP_MOMX_tmp(KA,IA,JA)
    real(RP) :: DAMP_MOMY_tmp(KA,IA,JA)
    real(RP) :: DAMP_MOMZ_tmp(KA,IA,JA)
    real(RP) :: DAMP_RHOT_tmp(KA,IA,JA)
    integer :: i, j, k
    !--------------------------------------------------
    
    !$omp parallel do private(i,j,k) collapse(2)
    do j=JS, JE
    do i=IS, IE
      do k=KS, KE
        DAMP_RHOT_tmp(k,i,j) = DAMP_DENS_fv(k,i,j) * DAMP_POTT_fv(k,i,j)
      end do
    end do
    end do

    call MeshFieldUtil_interp_FVtoDG_3D( DAMP_DENS_dg_, &
      DAMP_DENS_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      GRID_NUDGING_FV2DG_INTERP_ORDER )
    call MeshFieldUtil_interp_FVtoDG_3D( DAMP_RHOT_dg_, &
      DAMP_RHOT_tmp, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      GRID_NUDGING_FV2DG_INTERP_ORDER )

    return
  end subroutine interp_FVtoDG_DAMP_vars

!OCL SERIAL
  subroutine add_nudge_forcing( DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, RHOT_tp, RHOQ_tp,                   & ! (inout)
    DENS_bgn_t, MOMX_bgn_t, MOMY_bgn_t, MOMZ_bgn_t, RHOT_bgn_t,                                         & ! (inout)
    DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv,                                               & ! (in)
    DAMP_DENS_fv, DAMP_VELX_fv, DAMP_VELY_fv, DAMP_VELZ_fv, DAMP_POTT_fv, DAMP_QTRC_fv,                 & ! (in)
    DAMP_DENS_dg_, DAMP_RHOT_dg_,                                                                       & ! (in)
    DAMP_alpha_DENS_fv, DAMP_alpha_VELX_fv, DAMP_alpha_VELY_fv, DAMP_alpha_VELZ_fv, DAMP_alpha_POTT_fv, & ! (in)
    DAMP_alpha_QTRC_fv,                                                                                 & ! (in)
    BND_QA, BND_IQ, BND_SMOOTHER_FACT, DENS00_fv,                                                       & ! (in)
    DG_DDENS, DG_MOMX, DG_MOMY, DG_MOMZ, DG_DRHOT, DG_QTRC, DG_DENS_hyd, DG_PRES_hyd, lcmesh3D, elem3D )  ! (in)
    use scale_const, only: &
       EPS    => CONST_EPS
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    integer,  intent(in) :: BND_QA
    real(RP), intent(inout) :: DENS_tp(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMX_tp(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMY_tp(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMZ_tp(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: RHOT_tp(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: DENS_bgn_t(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMX_bgn_t(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMY_bgn_t(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: MOMZ_bgn_t(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(inout) :: RHOT_bgn_t(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(inout) :: RHOQ_tp(QA)
    real(RP), intent(in) :: DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMX_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMY_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_fv(KA,IA,JA)
    real(RP), intent(in) :: RHOT_fv(KA,IA,JA)
    real(RP), intent(in) :: QTRC_fv(KA,IA,JA,QA)
    real(RP), intent(in) :: DAMP_DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELX_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELY_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_VELZ_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_POTT_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_QTRC_fv(KA,IA,JA,BND_QA)
    real(RP), intent(in) :: DAMP_DENS_dg_(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DAMP_RHOT_dg_(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DAMP_alpha_DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_alpha_VELX_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_alpha_VELY_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_alpha_VELZ_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_alpha_POTT_fv(KA,IA,JA)
    real(RP), intent(in) :: DAMP_alpha_QTRC_fv(KA,IA,JA,BND_QA)
    integer,  intent(in) :: BND_IQ(QA)
    real(RP), intent(in) :: BND_SMOOTHER_FACT
    real(RP), intent(in) :: DENS00_fv(KA,IA,JA)
    real(RP), intent(in) :: DG_DDENS(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DG_MOMX(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DG_MOMY(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DG_MOMZ(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DG_DRHOT(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(inout) :: DG_QTRC(QA)
    real(RP), intent(in) :: DG_DENS_hyd(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DG_PRES_hyd(elem3D%Np,lcmesh3D%NeA)

    integer :: kelem, ke2D
    integer :: ke_x, ke_y, ke_z
    integer :: i, j, k
    integer :: iq, iqb

    real(RP) :: diff_fv(KA,IA,JA)
    real(RP) :: damp_tmp_fv(KA,IA,JA)
    real(RP) :: damp_tmp_qtrc(elem3D%Np)
    real(RP) :: DENS_tq_fv(KA,IA,JA)

    real(RP) :: var_avg, vol
    real(RP) :: int_w(elem3D%Np,lcmesh3D%ne)

    real(RP) :: RHOT_hyd(elem3D%Np)
    real(RP) :: sqrt_G_11(elem3D%Np), sqrt_G_22(elem3D%Np)
    !--------------------------------------------

    !$omp parallel private(kelem,vol)
    !$omp workshare
    DENS_tq_fv(:,:,:) = 0.0_RP
    !$omp end workshare
    !$omp do 
    do kelem=lcmesh3D%NeS, lcmesh3D%NeE
      int_w(:,kelem) = lcmesh3D%Gsqrt(:,kelem)*lcmesh3D%J(:,kelem)*elem3D%IntWeight_lgl(:)
      vol = sum( int_w(:,kelem) )      
      int_w(:,kelem) = int_w(:,kelem) / vol
    end do
    !$omp end parallel

    do iq = 1, QA
      iqb = BND_IQ(iq)
      if ( iqb > 0 ) then
        !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem, var_avg, damp_tmp_qtrc)
        !$omp do OMP_SCHEDULE_ collapse(2)
        do j = JS-1, JE+1
        do i = max(IS-1,1), min(IE+1,IA)
        do k = KS, KE
            diff_fv(k,i,j) = QTRC_fv(k,i,j,iq) - DAMP_QTRC_fv(k,i,j,iqb)
        enddo
        enddo
        enddo
        !$omp do OMP_SCHEDULE_ collapse(2)
        do ke_y=1, lcmesh3D%NeY
        do ke_x=1, lcmesh3D%NeX
        do ke_z=1, lcmesh3D%NeZ
          kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
          i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

          damp_tmp_qtrc(:) = - DAMP_alpha_QTRC_fv(k,i,j,iqb) &
              * ( diff_fv(k,i,j)                                                                                        & ! rayleigh damping
                - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) & ! horizontal smoother
                 * 0.125_RP * BND_SMOOTHER_FACT ) 
          !- Damp higher modes in DG
          ! var_avg = sum( int_w(:,kelem) * ( DG_QTRC(iq)%ptr%val(:,kelem) ) )
          ! damp_tmp_qtrc(:) = damp_tmp_qtrc(:) &
          !   - DAMP_alpha_QTRC_fv(k,i,j,iqb) &
          !     * ( ( DG_QTRC(iq)%ptr%val(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )

          RHOQ_tp(iq)%ptr%val(:,kelem) = RHOQ_tp(iq)%ptr%val(:,kelem) + damp_tmp_qtrc(:)  * DENS00_fv(k,i,j)          
          DENS_tq_fv(k,i,j) = DENS_tq_fv(k,i,j) + sum( int_w(:,kelem) * damp_tmp_qtrc(:) ) * DENS00_fv(k,i,j) * TRACER_MASS(iq) ! only for mass tracer
        enddo
        enddo
        enddo        
        !$omp end parallel
      end if
    end do
    
    !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem,var_avg)
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-1, JE+1
    do i = max(IS-1,1), min(IE+1,IA)
    do k = KS, KE
      diff_fv(k,i,j) = DENS_fv(k,i,j) - DAMP_DENS_fv(k,i,j)
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z
      
      DENS_bgn_t(:,kelem) = - DAMP_alpha_DENS_fv(k,i,j) &
                          * ( diff_fv(k,i,j)                                                                                        & ! rayleigh damping
                            - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) & ! horizontal smoother
                              * 0.125_RP * BND_SMOOTHER_FACT ) &
                          + DENS_tq_fv(k,i,j) * ( 0.5_RP - sign( 0.5_RP, DAMP_alpha_DENS_fv(k,i,j)-EPS ) ) ! density change due to rayleigh damping for tracers

      !- Alternative procedure
      ! DENS_bgn_t(:,kelem) = &
      !     DENS_damp_fv(k,i,j) &
      !   - DAMP_alpha_DENS_fv(k,i,j) * ( &
      !    ( DG_DENS_hyd(:,kelem) + DG_DDENS(:,kelem) - DAMP_DENS_dg_(:,kelem) ) )

      !- Damp higher modes in DG
      ! var_avg = sum( int_w(:,kelem) * ( DG_DENS_hyd(:,kelem) + DG_DDENS(:,kelem) ) )
      ! DENS_bgn_t(:,kelem) = DENS_bgn_t(:,kelem) &
      !   - DAMP_alpha_DENS_fv(k,i,j) * ( DG_DENS_hyd(:,kelem) + DG_DDENS(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )

      DENS_tp(:,kelem) = DENS_tp(:,kelem) + DENS_bgn_t(:,kelem)
    enddo
    enddo
    enddo
    !$omp end parallel
    
    !- MOMZ -------

    !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem,var_avg)
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-1, JE+1
    do i = IS-1, IE+1
    do k = KS, KE-1
        diff_fv(k,i,j) = MOMZ_fv(k,i,j) - DAMP_VELZ_fv(k,i,j) * ( DENS_fv(k,i,j)+DENS_fv(k+1,i,j) ) * 0.5_RP
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS, IE
      do k = KS, KE-1
        damp_tmp_fv(k,i,j) = - DAMP_alpha_VELZ_fv(k,i,j) * ( &
              - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) &
              * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
      enddo
      damp_tmp_fv(KS-1,i,j) = 0.0_RP
      damp_tmp_fv(KE,i,j) = 0.0_RP
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      diff_fv(k,i,j) = 0.5_RP * ( ( MOMZ_fv(k-1,i,j) + MOMZ_fv(k,i,j) ) - DENS_fv(k,i,j) * ( DAMP_VELZ_fv(k-1,i,j) + DAMP_VELZ_fv(k,i,j) ) )
      MOMZ_bgn_t(:,kelem) = &
          - 0.5_RP * ( DAMP_alpha_VELZ_fv(k-1,i,j) + DAMP_alpha_VELZ_fv(k,i,j) ) * ( &
              diff_fv(k,i,j) )                                                       & ! rayleigh damping
          + 0.5_RP * ( damp_tmp_fv(k-1,i,j) + damp_tmp_fv(k,i,j) )                   & ! horizontal smoother
          + DENS_bgn_t(:,kelem) * 0.5_RP * ( MOMZ_fv(k-1,i,j) + MOMZ_fv(k,i,j) ) / DENS_fv(k,i,j)
      
      !- Damp higher modes in DG
      ! var_avg = sum( int_w(:,kelem) * ( DG_MOMZ(:,kelem) ) )
      ! MOMZ_bgn_t(:,kelem) = MOMZ_bgn_t(:,kelem) &
      !   - 0.5_RP * ( DAMP_alpha_VELZ_fv(k-1,i,j) + DAMP_alpha_VELZ_fv(k,i,j) ) &
      !     * ( ( DG_MOMZ(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )

      MOMZ_tp(:,kelem) = MOMZ_tp(:,kelem) + MOMZ_bgn_t(:,kelem)
    enddo
    enddo
    enddo
    !$omp end parallel

    !- MOMX -------
    !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem,ke2D,var_avg,sqrt_G_11)
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-1, JE+1
    do i = IS-2, IE+1
    do k = KS, KE
        diff_fv(k,i,j) = MOMX_fv(k,i,j) - DAMP_VELX_fv(k,i,j) * ( DENS_fv(k,i,j)+DENS_fv(k,i+1,j) ) * 0.5_RP
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS, JE
    do i = IS-1, IE
    do k = KS, KE
      damp_tmp_fv(k,i,j) = - DAMP_alpha_VELX_fv(k,i,j) * ( &
              - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) &
              * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      ke2D = lcmesh3D%EMap3Dto2D(kelem)
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      diff_fv(k,i,j) = 0.5_RP * ( ( MOMX_fv(k,i-1,j) + MOMX_fv(k,i,j) ) - DENS_fv(k,i,j) * ( DAMP_VELX_fv(k,i-1,j) + DAMP_VELX_fv(k,i,j) ) )
      sqrt_G_11(:) = sqrt(lcmesh3D%G_ij(elem3D%IndexH2Dto3D(:),ke2D,1,1))

      MOMX_bgn_t(:,kelem) = &
          - 0.5_RP * ( DAMP_alpha_VELX_fv(k,i-1,j) + DAMP_alpha_VELX_fv(k,i,j) ) * ( &
            diff_fv(k,i,j) * sqrt_G_11(:) )                                          & ! rayleigh damping
          + ( 0.5_RP * ( damp_tmp_fv(k,i-1,j) + damp_tmp_fv(k,i,j) )                 & ! horizontal smoother
          + DENS_bgn_t(:,kelem) * 0.5_RP * ( MOMX_fv(k,i-1,j) + MOMX_fv(k,i,j) ) / DENS_fv(k,i,j) ) * sqrt_G_11(:)

      !- Damp higher modes in DG
      ! var_avg = sum( int_w(:,kelem) * ( DG_MOMX(:,kelem) ) )
      ! MOMX_bgn_t(:,kelem) = MOMX_bgn_t(:,kelem) &
      !   - 0.5_RP * ( DAMP_alpha_VELX_fv(k,i-1,j) + DAMP_alpha_VELX_fv(k,i,j) ) &
      !     * ( ( DG_MOMX(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )
          
      MOMX_tp(:,kelem) = MOMX_tp(:,kelem) + MOMX_bgn_t(:,kelem)
    enddo
    enddo
    enddo
    !$omp end parallel

    !- MOMY -------

    !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem,ke2D,var_avg,sqrt_G_22)
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-2, JE+1
    do i = IS, IE
    do k = KS, KE
        diff_fv(k,i,j) = MOMY_fv(k,i,j) - DAMP_VELY_fv(k,i,j) * ( DENS_fv(k,i,j)+DENS_fv(k,i,j+1) ) * 0.5_RP
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-1, JE
    do i = IS, IE
    do k = KS, KE
        damp_tmp_fv(k,i,j) = - DAMP_alpha_VELY_fv(k,i,j) * ( &
              - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) &
              * 0.125_RP * BND_SMOOTHER_FACT ) ! horizontal smoother
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      ke2D = lcmesh3D%EMap3Dto2D(kelem)
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      diff_fv(k,i,j) = 0.5_RP * ( ( MOMY_fv(k,i,j-1) + MOMY_fv(k,i,j) ) - DENS_fv(k,i,j) * ( DAMP_VELY_fv(k,i,j-1) + DAMP_VELY_fv(k,i,j) ) )
      sqrt_G_22(:) = sqrt(lcmesh3D%G_ij(elem3D%IndexH2Dto3D(:),ke2D,2,2))

      MOMY_bgn_t(:,kelem) = &
          - 0.5_RP * ( DAMP_alpha_VELY_fv(k,i,j-1) + DAMP_alpha_VELY_fv(k,i,j) ) * ( &
            diff_fv(k,i,j) * sqrt_G_22(:) )                            & ! rayleigh damping
          + ( 0.5_RP * ( damp_tmp_fv(k,i,j-1) + damp_tmp_fv(k,i,j) )   & ! horizontal smoother
          + DENS_bgn_t(:,kelem) * 0.5_RP * ( MOMY_fv(k,i,j-1) + MOMY_fv(k,i,j) ) / DENS_fv(k,i,j) ) * sqrt_G_22(:) 

      !- Damp higher modes in DG
      ! var_avg = sum( int_w(:,kelem) * ( DG_MOMY(:,kelem) ) )
      ! MOMY_bgn_t(:,kelem) = MOMY_bgn_t(:,kelem) &
      !   - 0.5_RP * ( DAMP_alpha_VELY_fv(k,i,j-1) + DAMP_alpha_VELY_fv(k,i,j) ) &
      !     * ( ( DG_MOMY(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )

      MOMY_tp(:,kelem) = MOMY_tp(:,kelem) + MOMY_bgn_t(:,kelem)
    enddo
    enddo
    enddo
    !$omp end parallel

    !- RHOT --------

    !$omp parallel private(i,j,k,ke_x,ke_y,ke_z,kelem,var_avg,RHOT_hyd)
    !$omp do OMP_SCHEDULE_ collapse(2)
    do j = JS-1, JE+2
    do i = max(IS-1,1), min(IE+2,IA)
    do k = KS, KE
        diff_fv(k,i,j) = RHOT_fv(k,i,j) - DAMP_POTT_fv(k,i,j) * DENS_fv(k,i,j)
    enddo
    enddo
    enddo
    !$omp do OMP_SCHEDULE_ collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      RHOT_hyd(:) = PRES00 / Rdry * ( DG_PRES_hyd(:,kelem) / PRES00 )**(CvDry/CpDry)

      damp_tmp_fv(k,i,j) = - DAMP_alpha_POTT_fv(k,i,j) &
                      * ( diff_fv(k,i,j)                                                                                        & ! rayleigh damping
                        - ( diff_fv(k,i-1,j) + diff_fv(k,i+1,j) + diff_fv(k,i,j-1) + diff_fv(k,i,j+1) - diff_fv(k,i,j)*4.0_RP ) & ! horizontal smoother
                        * 0.125_RP * BND_SMOOTHER_FACT )
      
      RHOT_bgn_t(:,kelem) = &
          + damp_tmp_fv(k,i,j) &
          + DENS_bgn_t(:,kelem) * RHOT_fv(k,i,j) / DENS_fv(k,i,j)

      !- Alternative procedure
      ! RHOT_bgn_t(:,kelem) = &
      !      + damp_fv(k,i,j) &
      !     - DAMP_alpha_POTT_fv(k,i,j) * ( &
      !       + ( RHOT_hyd(:) + DG_DRHOT(:,kelem) - DAMP_RHOT_dg_(:,kelem) )  ) &
      !     + DENS_damp_dg(:,kelem) * ( RHOT_hyd(:) + DG_DRHOT(:,kelem) ) / ( DG_DENS_hyd(:,kelem) + DG_DDENS(:,kelem) )
          
      !- Damp higher modes in DG
      ! var_avg = sum( int_w(:,kelem) * ( RHOT_hyd(:) + DG_DRHOT(:,kelem) ) )
      ! RHOT_bgn_t(:,kelem) = RHOT_bgn_t(:,kelem) &
      !   - 0.5_RP * DAMP_alpha_POTT_fv(k,i,j) &
      !     * ( ( RHOT_hyd(:) + DG_DRHOT(:,kelem) - var_avg ) * BND_SMOOTHER_FACT )

      RHOT_tp(:,kelem) = RHOT_tp(:,kelem) + RHOT_bgn_t(:,kelem)
    enddo
    enddo
    enddo
    !$omp end parallel

    return
  end subroutine add_nudge_forcing

end module scale_atmos_dyn_tstep_large_dgm_heve
