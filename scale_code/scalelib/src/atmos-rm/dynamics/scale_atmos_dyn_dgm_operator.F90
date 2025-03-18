!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics / DGM operator
!!
!! @par Description
!!          differential opeatros with DGM for Atmospheric dynamical process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_dgm_operator
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_tracer
  use scale_prc

  use scale_sparsemat, only: SparseMat

  use scale_fem_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_fem_element_hexahedral, only:  HexahedralElement
  use scale_fem_localmesh_2d, only: LocalMesh2D
  use scale_fem_localmesh_3d, only: LocalMesh3D
  use scale_fem_mesh_cubedom3d, only: MeshCubeDom3D
  use scale_fem_mesh_topography, only: MeshTopography
  use scale_fem_mesh_mapprojection, only: MeshMapProjection
  use scale_fem_element_modalfilter, only: ModalFilter

  use scale_atmos_dyn_dgm_bnd, only: AtmDynBnd

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
  public :: ATMOS_DYN_DGM_operator_setup
  public :: ATMOS_DYN_DGM_operator_finalize

  public :: ATMOS_DYN_DGM_operator_setup_vi
  public :: ATMOS_DYN_DGM_operator_setup_metric_bc

  public :: ATMOS_DYN_DGM_operator_set_fv_phytend
  public :: ATMOS_DYN_DGM_operator_apply_modalfilter
  public :: ATMOS_DYN_DGM_operator_addtend_spongelayer

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  type(HexahedralElement), public :: elem3D
  type(MeshCubeDom3D), public, target :: mesh3D
  type(AtmDynBnd), public :: dyn_bnd

  type(MeshTopography), public :: topography
  type(MeshMapProjection), public :: mapprojection

  character(len=H_SHORT) :: SpMV_storage_format = 'ELL' ! CSR or ELL  
  type(SparseMat), public :: Dx
  type(SparseMat), public :: Dy
  type(SparseMat), public :: Dz
  type(SparseMat), public :: Lift
  type(SparseMat), public :: FaceIntMat

  real(RP), public, allocatable :: IntrpMat_VPOrdM1(:,:)

  logical, public :: MODALFILTER_FLAG   = .false.
  type(ModalFilter), public :: modal_filter_3d

  logical, public :: ONLY_TRCADV_FLAG          = .false.
  logical, public :: TRCADV_MODALFILTER_FLAG   = .false.
  type(ModalFilter), public :: modal_filter_3d_tracer
  logical, public :: TRCADV_disable_limiter    = .false.

  logical, public :: SPONGELAYER_FLAG = .false.
  real(RP) :: SL_WDAMP_TAU        = -1.0_RP !> the maximum tau for Rayleigh damping of w [s]
  real(RP) :: SL_WDAMP_HEIGHT     = -1.0_RP !> the height to start apply Rayleigh damping [m]
  integer  :: SL_WDAMP_LAYER      = -1      !> the vertical number of finite element to start apply Rayleigh damping [num]
  logical  :: SL_HORIVELDAMP_FLAG = .false. !> Is the horizontal velocity damped? 


  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  logical :: hevi_flag
  integer :: PHYTEND_interp_order

  integer :: vcoord_type_id
  integer :: FV_TOPO_interp_order
  
contains

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_setup()
    use scale_prc, only: &
      PRC_myrank, PRC_nprocs
    use scale_prc_cartesC, only: &
      PRC_NUM_X, PRC_NUM_Y, &
      PRC_PERIODIC_X, PRC_PERIODIC_Y
    use scale_atmos_grid_cartesC_index, only: &
      IMAX, JMAX, KMAX,               &
      IAG, IHALO, JAG, JHALO, KS, KE, &
      IMAXG, JMAXG
    use scale_atmos_grid_cartesC, only: &
      FXG => ATMOS_GRID_CARTESC_FXG, &
      FYG => ATMOS_GRID_CARTESC_FYG, &
      FZ => ATMOS_GRID_CARTESC_FZ
    use scale_fem_mesh_base2d, only: MeshBase2D
    use scale_fem_meshutil_vcoord, only: &
      MeshUtil_get_VCoord_TypeID
    implicit none

    integer  :: PolyOrder_h       = 3
    integer  :: PolyOrder_v       = 3
    logical  :: LumpedMassMatFlag = .false.
    character(len=H_MID) :: VERTICAL_COORD_NAME = "TERRAIN_FOLLOWING"

    namelist / PARAM_ATMOS_DYN_DGM / &
      PolyOrder_h, PolyOrder_v, &
      LumpedMassMatFlag,        &
      VERTICAL_COORD_NAME,      &
      FV_TOPO_interp_order,     &
      MODALFILTER_FLAG,         &
      ONLY_TRCADV_FLAG,         &
      TRCADV_MODALFILTER_FLAG,  &
      TRCADV_disable_limiter,   &
      PHYTEND_interp_order,     &
      SPONGELAYER_FLAG
    
    integer :: ierr

    integer :: p1, p2, p_
    real(RP), allocatable :: invV_VPOrdM1(:,:)

    class(MeshBase2D), pointer :: mesh2D 
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("PARAM_ATMOS_DYN_DGM",*) 'Setup'

    FV_TOPO_interp_order = 0
    PHYTEND_interp_order = 2
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_DGM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("ATMOS_DYN_DGM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("ATMOS_DYN_DGM_setup",*) 'Not appropriate names in namelist PARAM_ATM_DYN_DGM. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_DYN_DGM)
    !--------------

    ! Setup mesh for constructing spatial operators with DGM 
    LOG_INFO("ATMOS_DYN_DGM_operator_setup",*) "Generate DG mesh.."

    call elem3D%Init( PolyOrder_h, PolyOrder_v, LumpedMassMatFlag )
    call mesh3D%Init( IMAXG, JMAXG, KMAX, &
       FXG(IHALO), FXG(IAG-IHALO), FYG(JHALO), FYG(JAG-JHALO), FZ(KS-1), FZ(KE), &
       PRC_PERIODIC_X, PRC_PERIODIC_Y, .false.,                                  &
       elem3D, 1, PRC_NUM_X, PRC_NUM_Y, PRC_nprocs, PRC_myrank, FZ(KS-1:KE) )
    call mesh3D%Generate()

    ! Setup object for managing topography
    call mesh3D%GetMesh2D( mesh2D )
    call topography%Init( "topo", mesh2D )
    vcoord_type_id = MeshUtil_get_VCoord_TypeID(VERTICAL_COORD_NAME)

    ! Setup object for managing map projeciton
    call mapprojection%Init( mesh2D )

    ! Setup spatial operators with DGM 
    call Dx%Init( elem3D%Dx1, storage_format=SpMV_storage_format )
    call Dy%Init( elem3D%Dx2, storage_format=SpMV_storage_format )
    call Dz%Init( elem3D%Dx3, storage_format=SpMV_storage_format )
    call Lift%Init( elem3D%Lift, storage_format=SpMV_storage_format )

    ! Setup modal filter
    call setup_modalfilter( elem3D )
    
    ! Setup modal filter and face integration matrix for tracer advection
    call setup_modalfilter_tracer( elem3D )
    call setup_FaceIntMat()

    ! Setup sponge layer
    call setup_sponge_layer()

    ! Setup vertical filter for buoyancy term
    allocate( IntrpMat_VPOrdM1(elem3D%Np,elem3D%Np) )
    allocate( invV_VPOrdM1(elem3D%Np,elem3D%Np) )
    
    InvV_VPOrdM1(:,:) = elem3D%invV
    do p2=1, elem3D%Nnode_h1D
    do p1=1, elem3D%Nnode_h1D
      p_ = p1 + (p2-1)*elem3D%Nnode_h1D + (elem3D%Nnode_v-1)*elem3D%Nnode_h1D**2
      InvV_VPOrdM1(p_,:) = 0.0_RP
    end do
    end do
    IntrpMat_VPOrdM1(:,:) = matmul(elem3D%V, invV_VPOrdM1)

    ! Initialize a object for boundary conditions for DG dycore
    call dyn_bnd%Init( mesh3D )

    !
    hevi_flag = .false.

    return
  end subroutine ATMOS_DYN_DGM_operator_setup

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_finalize()
    implicit none
    !---------------------------------------------
    
    deallocate( IntrpMat_VPOrdM1 )
    call dyn_bnd%Final()

    call Dx%Final(); call Dy%Final(); call Dz%Final()
    call Lift%Final()

    call topography%Final()
    call mapprojection%Final()

    call mesh3D%Final()
    call elem3D%Final()

    if ( MODALFILTER_FLAG ) then
      call modal_filter_3d%Final()
    end if

    return
  end subroutine ATMOS_DYN_DGM_operator_finalize

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_setup_vi() 
    implicit none
    hevi_flag = .true.
    return
  end subroutine ATMOS_DYN_DGM_operator_setup_vi
  
!> Setup metric & BC information
!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_setup_metric_bc() 
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO
    use scale_comm_fem_meshfield_rectdom2d, only: MeshFieldCommRectDom2D
    use scale_comm_fem_meshfield_cubedom3d, only: MeshFieldCommCubeDom3D

    use scale_const, only: &
      PI => CONST_PI
    implicit none

    ! type(LocalMesh2D), pointer :: lcmesh2D
    type(MeshFieldCommCubeDom3D) :: comm3D
    type(MeshFieldCommRectDom2D) :: comm2D
    !------------------------------------------------------
    
    !--- Setup terrain follwoing coordinate
    
    call comm2D%Init( 1, 0, 0, mesh3D%mesh2D )
    call comm3D%Init( 2, 1, 0, mesh3D )

    call topography%SetVCoordinate( mesh3D, &
      vcoord_type_id, mesh3D%zmax_gl, comm3D, comm2D )
    
    call comm2D%Final()
    call comm3D%Final()    

    !--- Setup map projection

    call comm2D%Init( 2, 0, 0, mesh3D%mesh2D )
    call mapprojection%SetHCoordinate( mesh3D, comm2D )
    call comm2D%Final()

    !--- Setup BC information
        
    call dyn_bnd%SetBCInfo()

    return
  end subroutine ATMOS_DYN_DGM_operator_setup_metric_bc

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_set_fv_phytend( tend_dg, &
      tend_fv, lcmesh3D, elem3D_, phytend_interp_ord )
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO
    use scale_fem_meshfield_util, only: &
      MeshFieldUtil_interp_FVtoDG_3D
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D_
    real(RP), intent(out) :: tend_dg(elem3D_%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: tend_fv(KA,IA,JA)
    integer, intent(in), optional :: phytend_interp_ord

    integer :: phytend_interp_ord_
    !----------------------------------

    if ( present(phytend_interp_ord) ) then
      phytend_interp_ord_ = phytend_interp_ord
    else
      phytend_interp_ord_ = PHYTEND_interp_order
    end if
    call MeshFieldUtil_interp_FVtoDG_3D( tend_dg, &
      tend_fv, lcmesh3D, elem3D_, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      phytend_interp_ord_, out_dg_flush_flag=.false. )

    return
  end subroutine ATMOS_DYN_DGM_operator_set_fv_phytend

!> Apply modal filter
!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_apply_modalfilter(  &
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_,     & ! (inout)
    lmesh, elem, filter, do_weight_Gsqrt     ) ! (in)

    use scale_fem_localmesh_base, only: LocalMeshBase
    use scale_fem_element_base, only: ElementBase
    use scale_fem_element_modalfilter, only: ModalFilter
    implicit none

    class(LocalMeshBase), intent(in) :: lmesh
    class(ElementBase), intent(in) :: elem    
    real(RP), intent(inout)  :: DDENS_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: MOMZ_(elem%Np,lmesh%NeA)
    real(RP), intent(inout)  :: DRHOT_(elem%Np,lmesh%NeA)
    type(ModalFilter), intent(in) :: filter
    logical, intent(in), optional :: do_weight_Gsqrt
    
    integer :: ke
    real(RP) :: tmp(elem%Np,5)
    integer :: ii, kk
    real(RP) :: Mik
    logical :: do_weight_Gsqrt_
    real(RP) :: RGsqrt(elem%Np)
    !------------------------------------

    if ( present( do_weight_Gsqrt ) ) then
      do_weight_Gsqrt_ = do_weight_Gsqrt
    else
      do_weight_Gsqrt_ = .false.
    end if

    if ( do_weight_Gsqrt_ ) then
      !$omp parallel do private( tmp, ii, kk, Mik, RGsqrt )
      do ke=lmesh%NeS, lmesh%NeE

        tmp(:,:) = 0.0_RP
        do ii=1, elem%Np
        do kk=1, elem%Np
          Mik = filter%FilterMat(ii,kk) * lmesh%Gsqrt(kk,ke)

          tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
          tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
          tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
          tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
          tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
        end do
        end do

        RGsqrt(:) = 1.0_RP / lmesh%Gsqrt(:,ke)
        DDENS_(:,ke) = tmp(:,1) * RGsqrt(:)
        MOMX_ (:,ke) = tmp(:,2) * RGsqrt(:)
        MOMY_ (:,ke) = tmp(:,3) * RGsqrt(:)
        MOMZ_ (:,ke) = tmp(:,4) * RGsqrt(:)
        DRHOT_(:,ke) = tmp(:,5) * RGsqrt(:)
      end do    

    else

      !$omp parallel do private( tmp, ii, kk, Mik )
      do ke=lmesh%NeS, lmesh%NeE

        tmp(:,:) = 0.0_RP
        do ii=1, elem%Np
        do kk=1, elem%Np
          Mik = filter%FilterMat(ii,kk)
          tmp(ii,1) = tmp(ii,1) + Mik * DDENS_(kk,ke)
          tmp(ii,2) = tmp(ii,2) + Mik * MOMX_ (kk,ke)
          tmp(ii,3) = tmp(ii,3) + Mik * MOMY_ (kk,ke)
          tmp(ii,4) = tmp(ii,4) + Mik * MOMZ_ (kk,ke)
          tmp(ii,5) = tmp(ii,5) + Mik * DRHOT_(kk,ke)
        end do
        end do      
        DDENS_(:,ke) = tmp(:,1)
        MOMX_ (:,ke) = tmp(:,2)
        MOMY_ (:,ke) = tmp(:,3)
        MOMZ_ (:,ke) = tmp(:,4)
        DRHOT_(:,ke) = tmp(:,5)
      end do

    end if

    return
  end subroutine ATMOS_DYN_DGM_operator_apply_modalfilter

!> Calculate tendencies with sponge layer
!OCL SERIAL
  subroutine ATMOS_DYN_DGM_operator_addtend_spongelayer( MOMX_dt, MOMY_dt, MOMZ_dt, &
    MOMX_, MOMY_, MOMZ_,                                                  &
    lmesh, elem   )

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    real(RP), intent(inout) :: MOMX_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMY_dt(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: MOMZ_dt(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMX_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMY_(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: MOMZ_(elem%Np,lmesh%NeA)

    integer :: ke
    integer :: ke_x, ke_y, ke_z
    integer :: keZtop
    real(RP) :: wdamp_coef(elem%Np)
    real(RP) :: zTop(elem%Nnode_h1D**2)
    real(RP) :: s
    !-----------------------------------------------------------------

    ! if ( this%hveldamp_flag ) then
    !   s = 1.0_RP
    ! else
    !   s = 0.0_RP
    ! end if
    s = 0.0_RP

    !$omp parallel do collapse(3) private(ke,keZtop,zTop,wdamp_coef)
    do ke_z = 1, lmesh%NeZ
    do ke_y = 1, lmesh%NeY
    do ke_x = 1, lmesh%NeX
      ke = ke_x + (ke_y-1)*lmesh%NeX + (ke_z-1)*lmesh%NeX*lmesh%NeY
      keZtop =  ke_x + (ke_y-1)*lmesh%NeX + (lmesh%NeZ-1)*lmesh%NeX*lmesh%NeY
      zTop(:) = lmesh%pos_en(elem%Hslice(:,elem%Nnode_v),keZtop,3)

      call calc_wdampcoef( &
        SL_WDAMP_TAU, SL_WDAMP_HEIGHT, lmesh%pos_en(:,ke,3), zTop(:), &
        elem%Nnode_h1D, elem%Nnode_v,                                 &
        wdamp_coef(:) )

      MOMX_dt(:,ke) = MOMX_dt(:,ke) - s * wdamp_coef(:) * MOMX_(:,ke)
      MOMY_dt(:,ke) = MOMY_dt(:,ke) - s * wdamp_coef(:) * MOMY_(:,ke)
      MOMZ_dt(:,ke) = MOMZ_dt(:,ke) - wdamp_coef(:) * MOMZ_(:,ke)
    end do
    end do
    end do

    return
  end subroutine ATMOS_DYN_DGM_operator_addtend_spongelayer

!- private ----------------------

!> Setup modal filter
!OCL SERIAL
  subroutine setup_modalfilter( refElem3D )
    implicit none
    class(HexahedralElement), target, intent(in) :: refElem3D

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_DGM_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_DGM_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_DGM_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_DGM_setup_modalfilter",*) 'Invalid names in namelist PARAM_ATMOS_DYN_MODALFILTER. Check!'
      call PRC_abort
    endif

    LOG_NML(PARAM_ATMOS_DYN_DGM_MODALFILTER)

    call modal_filter_3d%Init( &
      refElem3D,                           & ! (in)
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)

    return
  end subroutine setup_modalfilter

!> Setup modal filter for tracers
!OCL SERIAL
  subroutine setup_modalfilter_tracer( refElem3D )
    implicit none

    class(HexahedralElement), target, intent(in) :: refElem3D

    real(RP) :: MF_ETAC_h  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_h = 36.0_RP
    integer  :: MF_ORDER_h = 16
    real(RP) :: MF_ETAC_v  = 2.0_RP/3.0_RP
    real(RP) :: MF_ALPHA_v = 36.0_RP
    integer  :: MF_ORDER_v = 16

    namelist /PARAM_ATMOS_DYN_DGM_TRACER_MODALFILTER/ &
      MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   &
      MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    

    integer :: ierr
    !---------------------------------------------------------------

    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_DGM_TRACER_MODALFILTER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_TRCADV3D_setup_modalfilter",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_TRCADV3D_setup_modalfilter",*) 'Invalid names in namelist PARAM_ATMOS_DYN_TRACER_MODALFILTER. Check!'
      call PRC_abort
    endif

    LOG_NML(PARAM_ATMOS_DYN_DGM_TRACER_MODALFILTER)

    call modal_filter_3d_tracer%Init( &
        refElem3D,                           & ! (in)
        MF_ETAC_h, MF_ALPHA_h, MF_ORDER_h,   & ! (in)
        MF_ETAC_v, MF_ALPHA_v, MF_ORDER_v    ) ! (in)

    return
  end subroutine setup_modalfilter_tracer

!OCL SERIAL
  subroutine setup_FaceIntMat()
    use scale_polynominal, only: Polynominal_GenGaussLobattoPtIntWeight
    implicit none

    real(RP) :: intWeight_lgl1DPts_h(elem3D%Nnode_h1D)
    real(RP) :: intWeight_lgl1DPts_v(elem3D%Nnode_v)   
    real(RP) :: intWeight_h(elem3D%Nnode_h1D*elem3D%Nnode_v) 
    real(RP) :: intWeight_v(elem3D%Nnode_h1D**2)
    real(RP) :: IntWeight(elem3D%Nfaces,elem3D%NfpTot)

    integer :: f
    integer :: i, j, k, l
    integer :: is, ie
    !-----------------------------------------------------------

    IntWeight(:,:) = 0.0_RP
    intWeight_lgl1DPts_h(:) = Polynominal_GenGaussLobattoPtIntWeight(elem3D%PolyOrder_h)
    intWeight_lgl1DPts_v(:) = Polynominal_GenGaussLobattoPtIntWeight(elem3D%PolyOrder_v)

    do f=1, elem3D%Nfaces_h
      do k=1, elem3D%Nnode_v
      do i=1, elem3D%Nnode_h1D
        l = i + (k-1)*elem3D%Nnode_h1D
        intWeight_h(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_v(k)
      end do
      end do

      is = (f-1)*elem3D%Nfp_h + 1
      ie = is + elem3D%Nfp_h - 1
      IntWeight(f,is:ie) = intWeight_h(:)
    end do

    do f=1, elem3D%Nfaces_v
      do j=1, elem3D%Nnode_h1D
      do i=1, elem3D%Nnode_h1D
        l = i + (j-1)*elem3D%Nnode_h1D
        intWeight_v(l) = intWeight_lgl1DPts_h(i) * intWeight_lgl1DPts_h(j)
      end do
      end do

      is = elem3D%Nfaces_h*elem3D%Nfp_h + (f-1)*elem3D%Nfp_v + 1
      ie = is + elem3D%Nfp_v - 1
      IntWeight(elem3D%Nfaces_h+f,is:ie) = intWeight_v(:)
    end do

    call FaceIntMat%Init( IntWeight )

    return
  end subroutine setup_FaceIntMat

!> Setup sponge layer
!OCL SERIAL
  subroutine setup_sponge_layer()
    implicit none

    namelist /PARAM_ATMOS_DYN_DGM_SPONGELAYER/ &
      SL_WDAMP_TAU,                            &                
      SL_WDAMP_HEIGHT,                         &
      SL_WDAMP_LAYER,                          &
      SL_HORIVELDAMP_FLAG

    integer :: ierr
    !---------------------------------------------------------------

    ! Sponge layer
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_DYN_DGM_SPONGELAYER,iostat=ierr)
    if( ierr < 0 ) then !--- missing
      LOG_INFO("ATMOS_DYN_DGM_operator_setup_spongelayer",*) 'Not found namelist. Default used.'
    else if( ierr > 0 ) then !--- fatal error
      LOG_ERROR("ATMOS_DYN_DGM_operator_setup_spongelayer",*) 'Not appropriate names in namelist PARAM_ATMOS_DYN_DGM_SPONGELAYER. Check!'
      call PRC_abort
    end if
    LOG_NML(PARAM_ATMOS_DYN_DGM_SPONGELAYER)

    return
  end subroutine setup_sponge_layer 

!OCL SERIAL
  subroutine calc_wdampcoef( &
    wdamp_tau, wdamp_height, z, zTop, Nnode_h1D, Nnode_v, & ! (in)
    wdamp_coef                                            ) ! (out)

    use scale_const, only: &
      PI => CONST_PI
    implicit none

    integer, intent(in) :: Nnode_h1D
    integer, intent(in) :: Nnode_v
    real(RP), intent(out) :: wdamp_coef(Nnode_h1D**2,Nnode_v)
    real(RP), intent(in) :: wdamp_tau
    real(RP), intent(in) :: wdamp_height
    real(RP), intent(in) :: z(Nnode_h1D**2,Nnode_v)
    real(RP), intent(in) :: zTop(Nnode_h1D**2)

    integer :: p_z
    real(RP) :: sw(Nnode_h1D**2)
    real(RP) :: r_wdamp_tau
    !-----------------------------------------------------------------

    r_wdamp_tau = 1.0_RP / wdamp_tau
    do p_z=1, Nnode_v
      wdamp_coef(:,p_z) = 0.25_RP * r_wdamp_tau                                        &
        * ( 1.0_RP + sign( 1.0_RP, z(:,p_z) - wdamp_height )                         ) &
        * ( 1.0_RP - cos( PI * (z(:,p_z) - wdamp_height)/(zTop(:) - wdamp_height) )  )
    end do
    return
  end subroutine calc_wdampcoef
end module scale_atmos_dyn_dgm_operator