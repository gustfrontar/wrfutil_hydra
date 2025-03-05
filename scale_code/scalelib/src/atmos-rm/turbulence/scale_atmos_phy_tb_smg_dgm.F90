!> module ATMOSPHERE / Physics Turbulence (discretized by DGM)
!!
!! @par Description
!!      Sub-grid scale turbulnce process
!!      Smagorinsky-type
!!
!! @author Yuta Kawai, Team SCALE
!!
!! @par Reference
!!  - Brown et al., 1994:
!!    Large-eddy simulaition of stable atmospheric boundary layers with a revised stochastic subgrid model.
!!    Roy. Meteor. Soc., 120, 1485-1512
!!  - Scotti et al., 1993:
!!    Generalized Smagorinsky model for anisotropic grids.
!!    Phys. Fluids A, 5, 2306-2308
!!  - Nishizawa et al., 2015:
!!    Influence of grid aspect ratio on planetary boundary layer turbulence in large-eddy simulations
!!    Geosci. Model Dev., 8, 3393–3419
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_tb_smg_dgm
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use scale_precision
  use scale_io
  use scale_prc
  use scale_prof
  use scale_const, only: &
    EPS  => CONST_EPS,     &
    GRAV => CONST_GRAV,    &
    Rdry => CONST_Rdry,    &
    CPdry => CONST_CPdry,  &
    CVdry => CONST_CVdry,  &
    PRES00 => CONST_PRE00, &
    KARMAN  => CONST_KARMAN, &
    RPlanet => CONST_RADIUS    

  use scale_sparsemat  
  use scale_fem_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_fem_element_hexahedral, only: HexahedralElement
  use scale_fem_localmesh_2d, only: LocalMesh2D  
  use scale_fem_localmesh_3d, only: LocalMesh3D
  use scale_fem_mesh_base3d, only: MeshBase3D
  use scale_fem_localmeshfield_base, only: LocalMeshField3D
  use scale_fem_meshfield_base, only: MeshField3D

  use scale_atmos_phy_tb_common_dgm, only: &
    atm_phy_tb_common_dgm_setup, &
    atm_phy_tb_common_dgm_Final, &
    bnd_info, &
    T11_ID => ATMOS_PHY_TB_AUX_T11_ID, T12_ID => ATMOS_PHY_TB_AUX_T12_ID, T13_ID => ATMOS_PHY_TB_AUX_T13_ID, &
    T21_ID => ATMOS_PHY_TB_AUX_T21_ID, T22_ID => ATMOS_PHY_TB_AUX_T22_ID, T23_ID => ATMOS_PHY_TB_AUX_T23_ID, &
    T31_ID => ATMOS_PHY_TB_AUX_T31_ID, T32_ID => ATMOS_PHY_TB_AUX_T32_ID, T33_ID => ATMOS_PHY_TB_AUX_T33_ID, &
    DF1_ID => ATMOS_PHY_TB_AUX_DIFFFLX1_ID, DF2_ID => ATMOS_PHY_TB_AUX_DIFFFLX2_ID, DF3_ID => ATMOS_PHY_TB_AUX_DIFFFLX3_ID, &
    TKE_ID => ATMOS_PHY_TB_DIAG_TKE_ID, NU_ID => ATMOS_PHY_TB_DIAG_NU_ID, KH_ID => ATMOS_PHY_TB_DIAG_KH_ID, &
    TB_MOMX_t_VID => ATMOS_PHY_TB_MOMX_t_ID,TB_MOMY_t_VID => ATMOS_PHY_TB_MOMY_t_ID,  &
    TB_MOMZ_t_VID => ATMOS_PHY_TB_MOMZ_t_ID, TB_RHOT_t_VID => ATMOS_PHY_TB_RHOT_t_ID, &
    ATMOS_PHY_TB_TENDS_NUM1
  use scale_atmos_phy_tb_dgm_vars, only: &
    DFQ1_ID => ATMOS_PHY_TB_AUXTRC_DFQ1_ID, DFQ2_ID => ATMOS_PHY_TB_AUXTRC_DFQ2_ID, DFQ3_ID => ATMOS_PHY_TB_AUXTRC_DFQ3_ID

    
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: atm_phy_tb_dgm_smg_Init
  public :: atm_phy_tb_dgm_smg_Final
  public :: atm_phy_tb_dgm_smg_cal

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: cal_del_flux_grad

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter   :: OneOverThree  = 1.0_RP / 3.0_RP
  real(RP), private, parameter   :: twoOverThree  = 2.0_RP / 3.0_RP
  real(RP), private, parameter   :: FourOverThree = 4.0_RP / 3.0_RP

  real(RP), private              :: Cs            = 0.13_RP ! Smagorinsky constant (Scotti et al. 1993)
  real(RP), private, parameter   :: PrN           = 0.7_RP  ! Prandtl number in neutral conditions
  real(RP), private, parameter   :: RiC           = 0.25_RP ! critical Richardson number
  real(RP), private, parameter   :: FmC           = 16.0_RP ! fum = sqrt(1 - c*Ri)
  real(RP), private, parameter   :: FhB           = 40.0_RP ! fuh = sqrt(1 - b*Ri)/PrN
  real(RP), private              :: RPrN                    ! 1 / PrN
  real(RP), private              :: RRiC                    ! 1 / RiC
  real(RP), private              :: OnemPrNovRiC            ! PrN / RiC
  
  ! for backscatter
  real(RP), private, parameter   :: CB   = 1.4_RP
  real(RP), private, parameter   :: CBt  = 0.45_RP
  real(RP), private, parameter   :: aN   = 0.47958315233127197_RP ! a_N = sqrt(0.23)
  real(RP), private, parameter   :: atN4 = 0.09_RP                ! a_{\theta N}^4 = 0.3**2
  real(RP), private, parameter   :: C1o  = aN**3
  real(RP), private, parameter   :: D1o  = PrN * atN4 / aN


  real(RP), private              :: filter_fac    = 2.0_RP
  real(RP), private              :: NU_MAX        = 10000.0_RP
  real(RP), private              :: tke_fac       

contains
!OCL SERIAL
  subroutine atm_phy_tb_dgm_smg_Init( mesh )
    implicit none    
    class(MeshBase3D), intent(in) :: mesh

    logical  :: consistent_tke = .true.

    namelist / PARAM_ATMOS_PHY_TB_SMG_DGM / &
      Cs,                                   &
      NU_MAX,                               &
      filter_fac,                           &
      consistent_tke
    
    integer :: ierr
    !--------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Setup'
    LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Smagorinsky-type Eddy Viscocity Model'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_TB_SMG_DGM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_TB_dgm_smg_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_TB_dgm_smg_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_TB_SMG_DGM. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_TB_SMG_DGM)

    RPrN         = 1.0_RP / PrN
    RRiC         = 1.0_RP / RiC
    OnemPrNovRiC = ( 1.0_RP - PrN ) * RRiC

    if ( consistent_tke ) then
      tke_fac = 1.0_RP
    else
      tke_fac = 0.0_RP
    end if

    !-
    call atm_phy_tb_common_dgm_setup( mesh )

    return
  end subroutine atm_phy_tb_dgm_smg_Init

!OCL SERIAL
  subroutine atm_phy_tb_dgm_smg_Final()
    implicit none
    !--------------------------------------------------------------------

    call atm_phy_tb_common_dgm_Final()
    return
  end subroutine atm_phy_tb_dgm_smg_Final
  
!OCL SERIAL
  subroutine atm_phy_tb_dgm_smg_cal( tb_vars, &
    dyn_vars, mesh3d )
    use scale_tracer, only: &
      QA, TRACER_ADVC, TRACER_NAME
    use scale_atmos_phy_tb_dgm_vars, only: &
      AtmosPhyTbDGM_vars
    use scale_atmos_dyn_dgm_vars, only: &
      AtmosDynDGM_vars
    use scale_atmos_dyn_dgm_operator, only: &
      Dx, Dy, Dz, Lift, dyn_bnd
    use scale_atmos_phy_tb_common_dgm, only: &
      atm_phy_tb_common_dgm_cal_tend, &
      atm_phy_tb_common_dgm_cal_grad_qtrc, &
      atm_phy_tb_common_dgm_cal_tend_qtrc
    implicit none
    type(AtmosPhyTbDGM_vars), intent(inout), target :: tb_vars
    type(AtmosDynDGM_vars), intent(inout), target :: dyn_vars
    class(MeshBase3D), intent(in), target :: mesh3D

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: ldomID

    type(MeshField3D), pointer :: DDENS, MOMX, MOMY, MOMZ, DRHOT
    type(MeshField3D), pointer :: DENS_hyd, PRES_hyd
    type(MeshField3D), pointer :: Rtot, CVtot, CPtot, PRES, PT
    type(MeshField3D), pointer :: QTRC, tb_RHOQ_t 

    logical :: cal_grad_flag
    integer :: iq
    !-------------------------------------------------

    call dyn_vars%GetVarsPtr( DDENS, MOMX, MOMY, MOMZ, DRHOT, &
       DENS_hyd, PRES_hyd )
    call dyn_vars%GetAuxVarsPtr( Rtot, CVtot, CPtot, PT, PRES )

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D =>  mesh3D%lcmesh_list(ldomID)

      call dyn_bnd%ApplyBC_Grad_TBVARS_lc( ldomID, &
        DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val,                 & ! (in)
        PT%local(ldomID)%val, PRES%local(ldomID)%val,                                                                    & ! (inout)
        DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val, Rtot%local(ldomID)%val, CPtot%local(ldomID)%val,         & ! (in)
        lcmesh3D%Gsqrt(:,:), lcmesh3D%GsqrtH(:,:), lcmesh3D%GIJ(:,:,1,1), lcmesh3D%GIJ(:,:,1,2), lcmesh3D%GIJ(:,:,2,2),  & ! (in)
        lcmesh3D%GI3(:,:,1), lcmesh3D%GI3(:,:,2),                                                                        & ! (in)
        lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3),                                 & ! (in)
        lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB,                                                                  & ! (in)
        lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D                                     ) ! (in)    
    
      call PROF_rapstart('ATM_PHY_TB_cal_grad', 2)
      call atm_phy_tb_dgm_smg_cal_grad_lc( &
        tb_vars%auxvars(T11_ID)%local(ldomID)%val, tb_vars%auxvars(T12_ID)%local(ldomID)%val, tb_vars%auxvars(T13_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T21_ID)%local(ldomID)%val, tb_vars%auxvars(T22_ID)%local(ldomID)%val, tb_vars%auxvars(T23_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T31_ID)%local(ldomID)%val, tb_vars%auxvars(T32_ID)%local(ldomID)%val, tb_vars%auxvars(T33_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(DF1_ID)%local(ldomID)%val, tb_vars%auxvars(DF2_ID)%local(ldomID)%val, tb_vars%auxvars(DF3_ID)%local(ldomID)%val,  &
        tb_vars%diagvars(TKE_ID)%local(ldomID)%val, tb_vars%diagvars(NU_ID)%local(ldomID)%val, tb_vars%diagvars(KH_ID)%local(ldomID)%val, &
        DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, DRHOT%local(ldomID)%val,         &
        DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val, PRES%local(ldomID)%val, PT%local(ldomID)%val,                             &
        Dx, Dy, Dz, Lift, lcmesh3D, mesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, &
        bnd_info(ldomID)%is_bound )
      call PROF_rapend('ATM_PHY_TB_cal_grad', 2)
    enddo

    !* Exchange halo data
    call PROF_rapstart('ATM_PHY_TB_exchange_prgv', 2)
    call tb_vars%Exchange_auxvars()
    call PROF_rapend('ATM_PHY_TB_exchange_prgv', 2)

    do ldomID=1, mesh3D%LOCAL_MESH_NUM
      lcmesh3D =>  mesh3D%lcmesh_list(ldomID)

      call dyn_bnd%ApplyBC_Grad_TBStress_lc( ldomID, &
        tb_vars%auxvars(T11_ID)%local(ldomID)%val, tb_vars%auxvars(T12_ID)%local(ldomID)%val, tb_vars%auxvars(T13_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T21_ID)%local(ldomID)%val, tb_vars%auxvars(T22_ID)%local(ldomID)%val, tb_vars%auxvars(T23_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T31_ID)%local(ldomID)%val, tb_vars%auxvars(T32_ID)%local(ldomID)%val, tb_vars%auxvars(T33_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(DF1_ID)%local(ldomID)%val, tb_vars%auxvars(DF2_ID)%local(ldomID)%val, tb_vars%auxvars(DF3_ID)%local(ldomID)%val,  &
        lcmesh3D%Gsqrt(:,:), lcmesh3D%GsqrtH(:,:), lcmesh3D%GIJ(:,:,1,1), lcmesh3D%GIJ(:,:,1,2), lcmesh3D%GIJ(:,:,2,2),  & ! (in)
        lcmesh3D%GI3(:,:,1), lcmesh3D%GI3(:,:,2),                                                                        & ! (in)
        lcmesh3D%normal_fn(:,:,1), lcmesh3D%normal_fn(:,:,2), lcmesh3D%normal_fn(:,:,3),                                 & ! (in)
        lcmesh3D%vmapM, lcmesh3D%vmapP, lcmesh3D%vmapB,                                                                  & ! (in)
        lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D                                     ) ! (in)   
           
      call PROF_rapstart('ATM_PHY_TB_cal_tend', 2)
      call atm_phy_tb_common_dgm_cal_tend( &
        tb_vars%tends(TB_MOMX_t_VID)%local(ldomID)%val, tb_vars%tends(TB_MOMY_t_VID)%local(ldomID)%val, tb_vars%tends(TB_MOMZ_t_VID)%local(ldomID)%val,     &
        tb_vars%tends(TB_RHOT_t_VID)%local(ldomID)%val,                                                                                   &
        tb_vars%auxvars(T11_ID)%local(ldomID)%val, tb_vars%auxvars(T12_ID)%local(ldomID)%val, tb_vars%auxvars(T13_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T21_ID)%local(ldomID)%val, tb_vars%auxvars(T22_ID)%local(ldomID)%val, tb_vars%auxvars(T23_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(T31_ID)%local(ldomID)%val, tb_vars%auxvars(T32_ID)%local(ldomID)%val, tb_vars%auxvars(T33_ID)%local(ldomID)%val,  &
        tb_vars%auxvars(DF1_ID)%local(ldomID)%val, tb_vars%auxvars(DF2_ID)%local(ldomID)%val, tb_vars%auxvars(DF3_ID)%local(ldomID)%val,  &
        tb_vars%diagvars(NU_ID)%local(ldomID)%val, tb_vars%diagvars(KH_ID)%local(ldomID)%val,                                             &
        DDENS%local(ldomID)%val, MOMX%local(ldomID)%val, MOMY%local(ldomID)%val, MOMZ%local(ldomID)%val, DRHOT%local(ldomID)%val,         &
        DENS_hyd%local(ldomID)%val, PRES_hyd%local(ldomID)%val, PRES%local(ldomID)%val, PT%local(ldomID)%val,                             &
        Dx, Dy, Dz, Lift, lcmesh3D, mesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D, &
        bnd_info(ldomID)%is_bound )
      call PROF_rapend('ATM_PHY_TB_cal_tend', 2)
    enddo

    cal_grad_flag = .true.
    do iq = 1, QA
      if ( .not. TRACER_ADVC(iq) ) cycle
      
      QTRC => dyn_vars%qtrc_vars(iq)
      tb_RHOQ_t => tb_vars%tends(ATMOS_PHY_TB_TENDS_NUM1 + iq)
      
      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_cal_grad_qtrc', 2)
        call atm_phy_tb_common_dgm_cal_grad_qtrc( &
          tb_vars%auxtrcvars(DFQ1_ID)%local(ldomID)%val, tb_vars%auxtrcvars(DFQ2_ID)%local(ldomID)%val, tb_vars%auxtrcvars(DFQ3_ID)%local(ldomID)%val, & ! (out)
          tb_vars%GRAD_DENS(1)%local(ldomID)%val, tb_vars%GRAD_DENS(2)%local(ldomID)%val, tb_vars%GRAD_DENS(3)%local(ldomID)%val,                      & ! (inout)
          tb_vars%diagvars(KH_ID)%local(ldomID)%val, QTRC%local(ldomID)%val, DDENS%local(ldomID)%val, DENS_hyd%local(ldomID)%val,                      & ! (in) 
          Dx, Dy, Dz,  Lift, lcmesh3D, lcmesh3D%refElem3D,  lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,                                            & ! (in)
          bnd_info(ldomID)%is_bound, cal_grad_flag )                                            ! (in)
        call PROF_rapend('ATM_PHY_TB_cal_grad_qtrc', 2)
      enddo
      cal_grad_flag = .false.
      call tb_vars%Exchange_auxtrcvars()

      do ldomID=1, mesh3D%LOCAL_MESH_NUM
        call PROF_rapstart('ATM_PHY_TB_cal_tend_qtrc', 2)
        call atm_phy_tb_common_dgm_cal_tend_qtrc( tb_RHOQ_t%local(ldomID)%val,  & ! (out)
          tb_vars%auxtrcvars(DFQ1_ID)%local(ldomID)%val, tb_vars%auxtrcvars(DFQ2_ID)%local(ldomID)%val, tb_vars%auxtrcvars(DFQ3_ID)%local(ldomID)%val, & ! (in)
          tb_vars%diagvars(KH_ID)%local(ldomID)%val, DDENS%local(ldomID)%val, DENS_hyd%local(ldomID)%val,                                              & ! (in) 
          Dx, Dy, Dz, Lift, lcmesh3D, lcmesh3D%refElem3D, lcmesh3D%lcmesh2D, lcmesh3D%lcmesh2D%refElem2D,                                              & ! (in)
          bnd_info(ldomID)%is_bound )                                             ! (in)
        call PROF_rapend('ATM_PHY_TB_cal_tend_qtrc', 2)
      enddo        
    enddo

    return
  end subroutine atm_phy_tb_dgm_smg_cal
  
!-------

!> Calculate parameterized stress tensor and eddy heat flux with turbulent model
!!
!OCL SERIAL  
  subroutine atm_phy_tb_dgm_smg_cal_grad_lc( &
    T11, T12, T13, T21, T22, T23, T31, T32, T33,                & ! (out)
    DF1, DF2, DF3,                                              & ! (out)
    TKE, Nu, Kh,                                                & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd,    & ! (in)
    PRES, PT,                                                   & ! (in)
    Dx, Dy, Dz, Lift, lmesh, elem, lmesh2D, elem2D,             & ! (in)
    is_bound                                                    ) ! (in)

    use scale_atmos_phy_tb_common_dgm, only: &
      atm_phy_tb_common_dgm_calc_lambda
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: T11(elem%Np,lmesh%NeA)      !< (1,1) component of stress tensor
    real(RP), intent(out) :: T12(elem%Np,lmesh%NeA)      !< (1,2) component of stress tensor
    real(RP), intent(out) :: T13(elem%Np,lmesh%NeA)      !< (1,3) component of stress tensor
    real(RP), intent(out) :: T21(elem%Np,lmesh%NeA)      !< (2,1) component of stress tensor
    real(RP), intent(out) :: T22(elem%Np,lmesh%NeA)      !< (2,2) component of stress tensor
    real(RP), intent(out) :: T23(elem%Np,lmesh%NeA)      !< (2,3) component of stress tensor
    real(RP), intent(out) :: T31(elem%Np,lmesh%NeA)      !< (3,1) component of stress tensor
    real(RP), intent(out) :: T32(elem%Np,lmesh%NeA)      !< (3,2) component of stress tensor
    real(RP), intent(out) :: T33(elem%Np,lmesh%NeA)      !< (3,3) component of stress tensor
    real(RP), intent(out)  :: DF1(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x1 direction / density
    real(RP), intent(out)  :: DF2(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x2 direction / density
    real(RP), intent(out)  :: DF3(elem%Np,lmesh%NeA)     !< Diffusive heat flux in x3 direction / density
    real(RP), intent(out) :: TKE(elem%Np,lmesh%NeA)      !< Parameterized turbulent kinetic energy
    real(RP), intent(out) :: Nu(elem%Np,lmesh%NeA)       !< Eddy viscosity
    real(RP), intent(out) :: Kh(elem%Np,lmesh%NeA)       !< Eddy diffusivity
    real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)   !< Density perturbation
    real(RP), intent(in)  :: MOMX_ (elem%Np,lmesh%NeA)   !< Momentum in x1 direction
    real(RP), intent(in)  :: MOMY_ (elem%Np,lmesh%NeA)   !< Momentum in x2 direction
    real(RP), intent(in)  :: MOMZ_ (elem%Np,lmesh%NeA)   !< Momentum in x3 direction
    real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)   !< Density x potential temperature perturbation
    real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA) !< Reference pressure in hydrostatic balance
    real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA) !< Reference density in hydrostatic balance
    real(RP), intent(in)  :: PRES(elem%Np,lmesh%NeA)     !< Pressure
    real(RP), intent(in)  :: PT(elem%Np,lmesh%NeA)       !< Potential temperature
    type(SparseMat), intent(in) :: Dx, Dy, Dz             !< Differential matrix managed by sparse matrix type
    type(SparseMat), intent(in) :: Lift                   !< Lifting matrix managed by sparse matrix type
    logical, intent(in) :: is_bound(elem%NfpTot,lmesh%Ne) !< Flag whether nodes are located at domain boundaries

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DENS(elem%Np), RDENS(elem%Np), RHOT(elem%Np), Q(elem%Np)
    real(RP) :: DdensDxi(elem%Np,3)
    real(RP) :: DVelDxi(elem%Np,3,3)
    real(RP) :: del_flux_rho (elem%NfpTot,lmesh%Ne,3)
    real(RP) :: del_flux_mom (elem%NfpTot,lmesh%Ne,3,3)
    real(RP) :: del_flux_rhot(elem%NfpTot,lmesh%Ne,3)

    real(RP) :: S11(elem%Np), S12(elem%Np), S22(elem%Np), S23(elem%Np), S31(elem%Np), S33(elem%Np)
    real(RP) :: SkkOvThree
    real(RP) :: TKEMulTwoOvThree
    real(RP) :: coef

    real(RP) :: Ri ! local gradient Richardson number
    real(RP) :: S2 ! (2SijSij)^1/2
    real(RP) :: fm ! factor in eddy viscosity which represents the stability dependence of 
                   ! the Brown et al (1994)'s subgrid model 
    real(RP) :: Pr ! Parandtl number (=Nu/Kh= fm/fh)

    real(RP) :: lambda  (elem%Np,lmesh%Ne) ! basic mixing length
    real(RP) :: lambda_r(elem%Np)          ! characteristic subgrid length scale 
    real(RP) :: E(elem%Np)  ! subgrid kinetic energy 
    real(RP) :: C1(elem%Np) ! factor in the relation with energy disspation rate, lambda_r, and E

    integer :: ke
    integer :: p
    !--------------------------------------------------------------------

    call cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,        & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT,            & ! (in)
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, is_bound )                                                   ! (in)

    call atm_phy_tb_common_dgm_calc_lambda( lambda, & ! (out)
      Cs, filter_fac, lmesh, elem, lmesh2D, elem2D  ) ! (in)
  
    !$omp parallel do private( &
    !$omp Fx, Fy, Fz, LiftDelFlx,                  &
    !$omp DENS, RHOT, RDENS, Q, DdensDxi, DVelDxi, &
    !$omp S11, S12, S22, S23, S31, S33,            &
    !$omp coef, SkkOvThree, TKEMulTwoOvThree,      &
    !$omp p, Ri, S2, fm, Pr, lambda_r, E, C1       )
    do ke=lmesh%NeS, lmesh%NeE
      !---
      DENS (:) = DENS_hyd(:,ke) + DDENS_(:,ke)
      RDENS(:) = 1.0_RP / DENS(:)
      RHOT(:) = DENS(:) * PT(:,ke)

      ! gradient of density
      call sparsemat_matmul( Dx, DENS, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,1), LiftDelFlx )
      DdensDxi(:,1) = lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:)

      call sparsemat_matmul( Dy, DENS, Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,2), LiftDelFlx )
      DdensDxi(:,2) = lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:)

      call sparsemat_matmul( Dz, DENS, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rho(:,ke,3), LiftDelFlx )
      DdensDxi(:,3) = lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:)
      
      ! gradient of u
      Q(:) = MOMX_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMX_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,1), LiftDelFlx )
      DVelDxi(:,1,1) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMX_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,1), LiftDelFlx )
      DVelDxi(:,2,1) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMX_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,1), LiftDelFlx )
      DVelDxi(:,3,1) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of v
      Q(:) = MOMY_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMY_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,2), LiftDelFlx )
      DVelDxi(:,1,2) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMY_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,2), LiftDelFlx )
      DVelDxi(:,2,2) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMY_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,2), LiftDelFlx )
      DVelDxi(:,3,2) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of w
      Q(:) = MOMZ_(:,ke) * RDENS(:)

      call sparsemat_matmul( Dx, MOMZ_(:,ke), Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,1,3), LiftDelFlx )
      DVelDxi(:,1,3) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, MOMZ_(:,ke), Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,2,3), LiftDelFlx )
      DVelDxi(:,2,3) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, MOMZ_(:,ke), Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_mom(:,ke,3,3), LiftDelFlx )
      DVelDxi(:,3,3) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)

      ! gradient of pt
      Q(:) = RHOT(:) * RDENS(:)

      call sparsemat_matmul( Dx, RHOT, Fx )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,1), LiftDelFlx )
      DF1(:,ke) = ( lmesh%Escale(:,ke,1,1) * Fx(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,1) ) * RDENS(:)

      call sparsemat_matmul( Dy, RHOT, Fy )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,2), LiftDelFlx )
      DF2(:,ke) = ( lmesh%Escale(:,ke,2,2) * Fy(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,2) ) * RDENS(:)

      call sparsemat_matmul( Dz, RHOT, Fz )
      call sparsemat_matmul( Lift, lmesh%Fscale(:,ke) * del_flux_rhot(:,ke,3), LiftDelFlx )
      DF3(:,ke) = ( lmesh%Escale(:,ke,3,3) * Fz(:) + LiftDelFlx(:) - Q(:) * DdensDxi(:,3) ) * RDENS(:)


      ! Calculate the component of strain velocity tensor
      S11(:) = DVelDxi(:,1,1)
      S12(:) = 0.5_RP * ( DVelDxi(:,1,2) + DVelDxi(:,2,1) )
      S22(:) = DVelDxi(:,2,2)
      S23(:) = 0.5_RP * ( DVelDxi(:,2,3) + DVelDxi(:,3,2) )
      S31(:) = 0.5_RP * ( DVelDxi(:,1,3) + DVelDxi(:,3,1) )
      S33(:) = DVelDxi(:,3,3)

      ! Caclulate eddy viscosity & eddy diffusivity
      
      do p=1, elem%Np
        S2 = 2.0_RP * ( S11(p)**2 + S22(p)**2 + S33(p)**2 ) &
           + 4.0_RP * ( S31(p)**2 + S12(p)**2 + S23(p)**2 )
        
        Ri = Grav / PT(p,ke) * DF3(p,ke) / max( S2, EPS )

        ! The Stability functions fm and fh are given by the appendix A of Brown et al. (1994). 
        if (Ri < 0.0_RP ) then ! unstable
          fm = sqrt( 1.0_RP - FmC * Ri )
          Nu(p,ke) = lambda(p,ke)**2 * sqrt( S2 ) * fm
          Pr = fm / sqrt( 1.0_RP - FhB * Ri ) * PrN
        else if ( Ri < RiC ) then ! stable
          fm = ( 1.0_RP - Ri * RRiC )**4
          Nu(p,ke) = lambda(p,ke)**2 * sqrt( S2 ) * fm
          Pr = PrN / ( 1.0_RP - OnemPrNovRiC * Ri )
        else ! strongly stable
          fm = 0.0_RP
          Nu(p,ke) = 0.0_RP
          Kh(p,ke) = 0.0_RP
          Pr = 1.0_RP
        end if

        if ( Ri < RiC ) then
          Kh(p,ke) = max( min( Nu(p,ke) / Pr, NU_MAX ), EPS )
          Nu(p,ke) = max( min( Nu(p,ke), NU_MAX ), EPS )
          Pr = Nu(p,ke) / Kh(p,ke)
          lambda_r(p) = lambda(p,ke) * sqrt( fm / sqrt( 1.0_RP - Ri/Pr ) )
        else
          lambda_r(p) = 0.0_RP
        end if
      enddo
      
!      Nu(:,ke) = 0.0_RP

      ! if ( backscatter ) then
      ! else 
        E (:) = Nu(:,ke)**3 / ( lambda_r(:)**4 + EPS )
        C1(:) = C1o
      ! end if

      ! TKE
      TKE(:,ke) = ( E(:) * lambda_r(:) / C1(:) )**twoOverThree

      !---

      do p=1, elem%Np
        TKEMulTwoOvThree = twoOverThree * TKE(p,ke) * tke_fac
        SkkOvThree = ( S11(p) + S22(p) + S33(p) ) * OneOverThree
        coef = 2.0_RP * Nu(p,ke)

        T11(p,ke) = DENS(p) * ( coef * ( S11(p) - SkkOvThree ) - TKEMulTwoOvThree )
        T12(p,ke) = DENS(p) * coef * S12(p)
        T13(p,ke) = DENS(p) * coef * S31(p)

        T21(p,ke) = DENS(p) * coef * S12(p)
        T22(p,ke) = DENS(p) * ( coef * ( S22(p) - SkkOvThree ) - TKEMulTwoOvThree )
        T23(p,ke) = DENS(p) * coef * S23(p)

        T31(p,ke) = DENS(p) * coef * S31(p)
        T32(p,ke) = DENS(p) * coef * S23(p)
        T33(p,ke) = DENS(p) * ( coef * ( S33(p) - SkkOvThree ) - TKEMulTwoOvThree )
      enddo

      DF1(:,ke) = Kh(:,ke) * DF1(:,ke)
      DF2(:,ke) = Kh(:,ke) * DF2(:,ke)
      DF3(:,ke) = Kh(:,ke) * DF3(:,ke)
    enddo
  
    return
  end subroutine atm_phy_tb_dgm_smg_cal_grad_lc

!-- private --------------------------------------------------------

!OCL SERIAL  
  subroutine cal_del_flux_grad( del_flux_rho, del_flux_mom, del_flux_rhot,  & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DENS_hyd, PRES_hyd, PT_,           & ! (in)
    nx, ny, nz, vmapM, vmapP, lmesh, elem, is_bound                         ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    real(RP), intent(out) ::  del_flux_rho(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(out) ::  del_flux_mom(elem%NfpTot*lmesh%Ne,3,3)
    real(RP), intent(out) ::  del_flux_rhot(elem%NfpTot*lmesh%Ne,3)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)
    logical, intent(in) :: is_bound(elem%NfpTot*lmesh%Ne)
    
    integer :: i, iP, iM
    real(RP) :: densM, densP
    real(RP) :: del
    real(RP) :: facx, facy, facz
    !------------------------------------------------------------------------
    
    !$omp parallel do private ( iM, iP,   &
    !$omp densM, densP,                   &
    !$omp del, facx, facy, facz           )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      densM = DDENS_(iM) + DENS_hyd(iM)
      densP = DDENS_(iP) + DENS_hyd(iP)

      if ( is_bound(i) ) then
        facx = 1.0_RP
        facy = 1.0_RP 
        facz = 1.0_RP
      else
        ! facx = 1.0_RP - sign(1.0_RP,nx(i))
        ! facy = 1.0_RP - sign(1.0_RP,ny(i))
        ! facz = 1.0_RP - sign(1.0_RP,nz(i))
        facx = 1.0_RP
        facy = 1.0_RP
        facz = 1.0_RP
      end if

      del = 0.5_RP * ( densP - densM )
      del_flux_rho(i,1) = facx * del * nx(i)
      del_flux_rho(i,2) = facy * del * ny(i)
      del_flux_rho(i,3) = facz * del * nz(i)

      del = 0.5_RP * ( MOMX_(iP) - MOMX_(iM) )
      del_flux_mom(i,1,1) = facx * del * nx(i)
      del_flux_mom(i,2,1) = facy * del * ny(i)
      del_flux_mom(i,3,1) = facz * del * nz(i)

      del = 0.5_RP * ( MOMY_(iP) - MOMY_(iM) )
      del_flux_mom(i,1,2) = facx * del * nx(i)
      del_flux_mom(i,2,2) = facy * del * ny(i)
      del_flux_mom(i,3,2) = facz * del * nz(i)

      del = 0.5_RP * ( MOMZ_(iP) - MOMZ_(iM) )
      del_flux_mom(i,1,3) = facx * del * nx(i)
      del_flux_mom(i,2,3) = facy * del * ny(i)
      del_flux_mom(i,3,3) = facz * del * nz(i)

      del = 0.5_RP * ( densP * PT_(iP) - densM * PT_(iM) )
      del_flux_rhot(i,1) = facx * del * nx(i)
      del_flux_rhot(i,2) = facy * del * ny(i)
      del_flux_rhot(i,3) = facz * del * nz(i)
    enddo

    return
  end subroutine cal_del_flux_grad

end module scale_atmos_phy_tb_smg_dgm
