!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics RK
!!
!! @par Description
!!          HEVE DGM scheme for Atmospheric dynamical process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tstep_short_dgm_heve
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc
#if defined DEBUG || defined QUICKDEBUG
  use scale_debug, only: &
     CHECK
  use scale_const, only: &
     UNDEF  => CONST_UNDEF, &
     IUNDEF => CONST_UNDEF2
#endif
  use scale_const, only: &
    GRAV => CONST_GRAV,  &
    Rdry => CONST_Rdry,  &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
  
  use scale_sparsemat, only: &
    SparseMat, sparsemat_matmul
  use scale_fem_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_fem_localmesh_2d, only: LocalMesh2D
  use scale_fem_localmesh_3d, only: LocalMesh3D

  use scale_atmos_dyn_dgm_vars, only: &
    DENS_VID => PRGVAR_DDENS_ID, RHOT_VID => PRGVAR_DRHOT_ID, &
    MOMX_VID => PRGVAR_MOMX_ID, MOMY_VID => PRGVAR_MOMY_ID,   &
    MOMZ_VID => PRGVAR_MOMZ_ID,                               &
    PRGVAR_NUM  
  use scale_atmos_dyn_dgm_operator, only: &
    IntrpMat_VPOrdM1
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tstep_short_dgm_heve_regist
  public :: ATMOS_DYN_Tstep_short_dgm_heve_setup
  public :: ATMOS_DYN_Tstep_short_dgm_heve

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, parameter :: VA_DGM_HEVE = 0 ! Dummy for driver module

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tstep_short_dgm_heve_regist( &
       ATMOS_DYN_TYPE, &
       VA_out,         &
       VAR_NAME,       &
       VAR_DESC,       &
       VAR_UNIT        )
    use scale_prc, only: &
       PRC_abort
    implicit none

    character(len=*),       intent(in)  :: ATMOS_DYN_TYPE
    integer,                intent(out) :: VA_out         !< number of prognostic variables
    character(len=H_SHORT), intent(out) :: VAR_NAME(:)    !< name   of the variables
    character(len=H_MID),   intent(out) :: VAR_DESC(:)    !< desc.  of the variables
    character(len=H_SHORT), intent(out) :: VAR_UNIT(:)    !< unit   of the variables
    !---------------------------------------------------------------------------

    if ( ATMOS_DYN_TYPE /= 'DGM-HEVE' ) then
       LOG_ERROR("ATMOS_DYN_Tstep_short_dgm_heve_regist",*) 'ATMOS_DYN_TYPE is not DGM-HEVE. Check!'
       call PRC_abort
    endif

    VA_out      = VA_DGM_HEVE
    VAR_NAME(:) = ""
    VAR_DESC(:) = ""
    VAR_UNIT(:) = ""

    LOG_NEWLINE
    LOG_INFO("ATMOS_DYN_Tstep_short_dgm_heve_regist",*) 'Register additional prognostic variables (HEVE)'
    if ( VA_out < 1 ) then
       LOG_INFO_CONT(*) '=> nothing.'
    endif

    return
  end subroutine ATMOS_DYN_Tstep_short_dgm_heve_regist

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tstep_short_dgm_heve_setup
    return
  end subroutine ATMOS_DYN_Tstep_short_dgm_heve_setup

  !-----------------------------------------------------------------------------
  subroutine ATMOS_DYN_Tstep_short_dgm_heve( &
      DENS_dt, MOMX_dt, MOMY_dt, MOMZ_dt, RHOT_dt,               & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_,               & ! (in) 
      DENS_hyd, PRES_hyd, PRES_hyd_ref, CORIOLIS,                & ! (in)
      Rtot, CVtot, CPtot,                                        & ! (in)
      Dx, Dy, Dz, Lift, lmesh, elem, lmesh2D, elem2D             ) ! (in)

      implicit none

      class(LocalMesh3D), intent(in) :: lmesh
      class(ElementBase3D), intent(in) :: elem
      class(LocalMesh2D), intent(in) :: lmesh2D
      class(ElementBase2D), intent(in) :: elem2D
      type(SparseMat), intent(in) :: Dx, Dy, Dz, Lift
      real(RP), intent(out) :: DENS_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMX_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMY_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: MOMZ_dt(elem%Np,lmesh%NeA)
      real(RP), intent(out) :: RHOT_dt(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DDENS_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMX_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMY_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: MOMZ_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DRHOT_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DPRES_(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: DENS_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: PRES_hyd_ref(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CORIOLIS(elem2D%Np,lmesh2D%NeA)
      real(RP), intent(in)  :: Rtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CVtot(elem%Np,lmesh%NeA)
      real(RP), intent(in)  :: CPtot(elem%Np,lmesh%NeA)
#ifdef HIST_TEND
    use scale_file_history, only: &
       FILE_HISTORY_in
#endif

    real(RP) :: Fx(elem%Np), Fy(elem%Np), Fz(elem%Np), LiftDelFlx(elem%Np)
    real(RP) :: DPRES_hyd(elem%Np), GradPhyd_x(elem%Np), GradPhyd_y(elem%Np)
    real(RP) :: del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP) :: del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP) :: RHOT_(elem%Np)
    real(RP) :: rdens_(elem%Np), u_(elem%Np), v_(elem%Np), w_(elem%Np), wt_(elem%Np)
    real(RP) :: Cori(elem%Np)
    real(RP) :: drho(elem%Np)
    real(RP) :: GsqrtV(elem%Np), RGsqrtV(elem%Np)

    integer :: ke, ke2d

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR
    !-----------------------------------------------------------------------


    call PROF_rapstart('cal_dyn_tend_bndflux', 3)
    call get_ebnd_flux( &
      del_flux, del_flux_hyd,                                                 & ! (out)
      DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd,        & ! (in)
      Rtot, CVtot, CPtot,                                                     & ! (in)
      lmesh%Gsqrt, lmesh%GI3(:,:,1), lmesh%GI3(:,:,2),                        & ! (in)    
      lmesh%normal_fn(:,:,1), lmesh%normal_fn(:,:,2), lmesh%normal_fn(:,:,3), & ! (in)
      lmesh%vmapM, lmesh%vmapP,                                               & ! (in)
      lmesh, elem, lmesh2D, elem2D )                                            ! (in)
    call PROF_rapend('cal_dyn_tend_bndflux', 3)
 
    !-----
    call PROF_rapstart('cal_dyn_tend_interior', 3)
    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( ke, ke2d, Cori,     &
    !$omp RHOT_, rdens_, u_, v_, w_, wt_,          &
    !$omp drho, DPRES_hyd, GradPhyd_x, GradPhyd_y, &
    !$omp GsqrtV, RGsqrtV,                         &
    !$omp Fx, Fy, Fz, LiftDelFlx )
    do ke = lmesh%NeS, lmesh%NeE
      !--
      ke2d = lmesh%EMap3Dto2D(ke)
      Cori(:) = CORIOLIS(elem%IndexH2Dto3D(:),ke2d)

      GsqrtV(:)  = lmesh%Gsqrt(:,ke) / lmesh%GsqrtH(elem%IndexH2Dto3D,ke2d)
      RGsqrtV(:) = 1.0_RP / GsqrtV(:)

      !--
      RHOT_(:) = P0ovR * (PRES_hyd(:,ke) * rP0)**rgamm + DRHOT_(:,ke)
      ! DPRES_(:) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT_(:) )**( CPtot(:,ke) / CVtot(:,ke) ) &
      !           - PRES_hyd(:,ke)

      rdens_(:) = 1.0_RP / (DDENS_(:,ke) + DENS_hyd(:,ke))
      u_ (:) = MOMX_(:,ke) * rdens_(:)
      v_ (:) = MOMY_(:,ke) * rdens_(:)
      w_ (:) = MOMZ_(:,ke) * rdens_(:)
      wt_(:) = w_(:) * RGsqrtV(:) + lmesh%GI3(:,ke,1) * u_(:) + lmesh%GI3(:,ke,2) * v_(:) 

      drho(:) = matmul(IntrpMat_VPOrdM1, DDENS_(:,ke))

      !-- Gradient hydrostatic pressure
      
      DPRES_hyd(:) = PRES_hyd(:,ke) - PRES_hyd_ref(:,ke)

      call sparsemat_matmul(Dx, GsqrtV(:) * DPRES_hyd(:), Fx)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,1) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,1), LiftDelFlx)
      GradPhyd_x(:) = lmesh%Escale(:,ke,1,1) * Fx(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      call sparsemat_matmul(Dy, GsqrtV(:) * DPRES_hyd(:), Fy)
      call sparsemat_matmul(Dz, GsqrtV(:) * lmesh%GI3(:,ke,2) * DPRES_hyd(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux_hyd(:,ke,2), LiftDelFlx)
      GradPhyd_y(:) = lmesh%Escale(:,ke,2,2) * Fy(:) &
                    + lmesh%Escale(:,ke,3,3) * Fz(:) &
                    + LiftDelFlx(:)

      !-- DENS
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * MOMX_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * MOMY_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( DDENS_(:,ke) + DENS_hyd(:,ke) ) * wt_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,DENS_VID), LiftDelFlx)

      DENS_dt(:,ke) = - ( &
            lmesh%Escale(:,ke,1,1) * Fx(:)    &
          + lmesh%Escale(:,ke,2,2) * Fy(:)    &
          + lmesh%Escale(:,ke,3,3) * Fz(:)    &
          + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)
      
      !-- MOMX
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * (  u_(:) * MOMX_(:,ke) + DPRES_(:,ke) ), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *    v_(:) * MOMX_(:,ke)                 , Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMX_(:,ke)                     &
                                                    + lmesh%GI3(:,ke,1) * DPRES_(:,ke)    ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMX_VID), LiftDelFlx)

      MOMX_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_x(:) * RGsqrtV(:)            &
          + Cori(:) * MOMY_(:,ke)

      !-- MOMY
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *    u_(:) * MOMY_(:,ke)                 , Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * (  v_(:) * MOMY_(:,ke) + DPRES_(:,ke) ), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMY_(:,ke)                     &
                                                    + lmesh%GI3(:,ke,2) * DPRES_(:,ke)    ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMY_VID), LiftDelFlx)

      MOMY_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) &
          - GradPhyd_y(:) * RGsqrtV(:)            &
          - Cori(:) * MOMX_(:,ke)

      !-- MOMZ
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) *   u_(:) * MOMZ_(:,ke), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) *   v_(:) * MOMZ_(:,ke), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * ( wt_(:) * MOMZ_(:,ke)           &
                                                    + RGsqrtV(:) * DPRES_(:,ke) ), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,MOMZ_VID), LiftDelFlx)
      
      MOMZ_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)       &
            + lmesh%Escale(:,ke,2,2) * Fy(:)       &
            + lmesh%Escale(:,ke,3,3) * Fz(:)       &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke)  &
          - Grav * drho(:)

      !-- RHOT
      call sparsemat_matmul(Dx, lmesh%Gsqrt(:,ke) * u_ (:) * RHOT_(:), Fx)
      call sparsemat_matmul(Dy, lmesh%Gsqrt(:,ke) * v_ (:) * RHOT_(:), Fy)
      call sparsemat_matmul(Dz, lmesh%Gsqrt(:,ke) * wt_(:) * RHOT_(:), Fz)
      call sparsemat_matmul(Lift, lmesh%Fscale(:,ke) * del_flux(:,ke,RHOT_VID), LiftDelFlx)
      
      RHOT_dt(:,ke) = &
          - ( lmesh%Escale(:,ke,1,1) * Fx(:)      &
            + lmesh%Escale(:,ke,2,2) * Fy(:)      &
            + lmesh%Escale(:,ke,3,3) * Fz(:)      &
            + LiftDelFlx(:) ) / lmesh%Gsqrt(:,ke) 

    end do

    call PROF_rapend('cal_dyn_tend_interior', 3)

    return
  end subroutine ATMOS_DYN_Tstep_short_dgm_heve

!-- private

!OCL SERIAL
  subroutine get_ebnd_flux( &
    del_flux, del_flux_hyd,                                          & ! (out)
    DDENS_, MOMX_, MOMY_, MOMZ_, DRHOT_, DPRES_, DENS_hyd, PRES_hyd, & ! (in)
    Rtot, CVtot, CPtot,                                              & ! (in)
    Gsqrt, G13, G23, nx, ny, nz,                                     & ! (in)
    vmapM, vmapP, lmesh, elem, lmesh2D, elem2D                       ) ! (in)

    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem  
    class(LocalMesh2D), intent(in) :: lmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) ::  del_flux(elem%NfpTot,lmesh%Ne,PRGVAR_NUM)
    real(RP), intent(out) ::  del_flux_hyd(elem%NfpTot,lmesh%Ne,2)
    real(RP), intent(in) ::  DDENS_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  MOMX_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMY_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  MOMZ_(elem%Np*lmesh%NeA)  
    real(RP), intent(in) ::  DRHOT_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DPRES_(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  DENS_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  PRES_hyd(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Rtot (elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CVtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  CPtot(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  Gsqrt(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G13(elem%Np*lmesh%NeA)
    real(RP), intent(in) ::  G23(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot,lmesh%Ne)
    real(RP), intent(in) :: nz(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot,lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot,lmesh%Ne)
    
    integer :: ke, i, iP(elem%NfpTot), iM(elem%NfpTot)
    integer :: ke2D
    real(RP) :: VelP(elem%NfpTot), VelM(elem%NfpTot), alpha(elem%NfpTot)
    real(RP) :: dpresP(elem%NfpTot), dpresM(elem%NfpTot)
    real(RP) :: GsqrtDensM(elem%NfpTot), GsqrtDensP(elem%NfpTot)
    real(RP) :: GsqrtRhotM(elem%NfpTot), GsqrtRhotP(elem%NfpTot)
    real(RP) :: GsqrtDDENS_P(elem%NfpTot), GsqrtDDENS_M(elem%NfpTot)
    real(RP) :: GsqrtMOMX_P(elem%NfpTot), GsqrtMOMX_M(elem%NfpTot)
    real(RP) :: GsqrtMOMY_P(elem%NfpTot), GsqrtMOMY_M(elem%NfpTot)
    real(RP) :: GsqrtMOMZ_P(elem%NfpTot), GsqrtMOMZ_M(elem%NfpTot)
    real(RP) :: GsqrtDRHOT_P(elem%NfpTot), GsqrtDRHOT_M(elem%NfpTot)
    real(RP) :: Phyd_P(elem%NfpTot), Phyd_M(elem%NfpTot)
    real(RP) :: Gsqrt_P(elem%NfpTot), Gsqrt_M(elem%NfpTot)
    real(RP) :: GsqrtV_P(elem%NfpTot), GsqrtV_M(elem%NfpTot)
    real(RP) :: G13_M(elem%NfpTot), G13_P(elem%NfpTot)
    real(RP) :: G23_M(elem%NfpTot), G23_P(elem%NfpTot)
    real(RP) :: Gnn_M(elem%NfpTot), Gnn_P(elem%NfpTot)

    real(RP) :: gamm, rgamm    
    real(RP) :: rP0
    real(RP) :: RovP0, P0ovR     
    !------------------------------------------------------------------------

    gamm  = CPDry / CvDry
    rgamm = CvDry / CpDry
    rP0   = 1.0_RP / PRES00
    RovP0 = Rdry * rP0
    P0ovR = PRES00 / Rdry

    !$omp parallel do private( &
    !$omp ke, iM, iP, ke2D,                                                             &
    !$omp alpha, VelM, VelP,                                                            &
    !$omp dpresM, dpresP, GsqrtDensM, GsqrtDensP, GsqrtRhotM, GsqrtRhotP,               &
    !$omp GsqrtMOMX_M, GsqrtMOMX_P, GsqrtMOMY_M, GsqrtMOMY_P, GsqrtMOMZ_M, GsqrtMOMZ_P, &
    !$omp GsqrtDDENS_M, GsqrtDDENS_P, GsqrtDRHOT_M, GsqrtDRHOT_P,                       &
    !$omp Phyd_M, Phyd_P,                                                               &
    !$omp Gsqrt_P, Gsqrt_M, GsqrtV_P, GsqrtV_M, G13_P, G13_M, G23_P, G23_M,             &
    !$omp Gnn_P, Gnn_M                                                                  )
    do ke=lmesh%NeS, lmesh%NeE
      iM(:) = vmapM(:,ke); iP(:) = vmapP(:,ke)
      ke2D = lmesh%EMap3Dto2D(ke)

      Gsqrt_M(:) = Gsqrt(iM)
      Gsqrt_P(:) = Gsqrt(iP)
      GsqrtV_M(:) = Gsqrt_M(:)
      GsqrtV_P(:) = Gsqrt_P(:)

      G13_M(:) = G13(iM)
      G13_P(:) = G13(iP)
      G23_M(:) = G23(iM)
      G23_P(:) = G23(iP)

      GsqrtDDENS_M(:) = Gsqrt_M(:) * DDENS_(iM)
      GsqrtDDENS_P(:) = Gsqrt_P(:) * DDENS_(iP)
      GsqrtMOMX_M (:) = Gsqrt_M(:) * MOMX_ (iM)
      GsqrtMOMX_P (:) = Gsqrt_P(:) * MOMX_ (iP)
      GsqrtMOMY_M (:) = Gsqrt_M(:) * MOMY_ (iM)
      GsqrtMOMY_P (:) = Gsqrt_P(:) * MOMY_ (iP)
      GsqrtMOMZ_M (:) = Gsqrt_M(:) * MOMZ_ (iM)
      GsqrtMOMZ_P (:) = Gsqrt_P(:) * MOMZ_ (iP)
      GsqrtDRHOT_M(:) = Gsqrt_M(:) * DRHOT_(iM)
      GsqrtDRHOT_P(:) = Gsqrt_P(:) * DRHOT_(iP)
      Phyd_M(:) = PRES_hyd(iM)
      Phyd_P(:) = PRES_hyd(iP)

      Gnn_M(:) = abs( nx(:,ke) ) + abs( ny(:,ke) ) &
               + ( 1.0_RP / GsqrtV_M(:)**2 + G13_M(:)**2 + G23_M(:)**2 ) * abs( nz(:,ke) )
      Gnn_P(:) = abs( nx(:,ke) ) + abs( ny(:,ke) ) &
               + ( 1.0_RP / GsqrtV_P(:)**2 + G13_P(:)**2 + G23_P(:)**2 ) * abs( nz(:,ke) )

      GsqrtDensM(:) = GsqrtDDENS_M(:) + Gsqrt_M(:) * DENS_hyd(iM)
      GsqrtDensP(:) = GsqrtDDENS_P(:) + Gsqrt_P(:) * DENS_hyd(iP)

      GsqrtRhotM(:) = Gsqrt_M(:) * P0ovR * (Phyd_M(:) * rP0)**rgamm + GsqrtDRHOT_M(:)
      GsqrtRhotP(:) = Gsqrt_P(:) * P0ovR * (Phyd_P(:) * rP0)**rgamm + GsqrtDRHOT_P(:)

      VelM(:) = ( GsqrtMOMX_M(:) * nx(:,ke) + GsqrtMOMY_M(:) * ny(:,ke)                    &
                + ( ( GsqrtMOMZ_M(:) / GsqrtV_M(:)                                         &
                    + G13_M(:) * GsqrtMOMX_M(:) + G23_M(:) * GsqrtMOMY_M(:) ) * nz(:,ke) ) &
                ) / GsqrtDensM(:)
      VelP(:) = ( GsqrtMOMX_P(:) * nx(:,ke) + GsqrtMOMY_P(:) * ny(:,ke)                    &
                + ( ( GsqrtMOMZ_P(:) / GsqrtV_P(:)                                         &
                    + G13_P(:) * GsqrtMOMX_P(:) + G23_P(:) * GsqrtMOMY_P(:) ) * nz(:,ke) ) &
                ) / GsqrtDensP(:)
        

      ! dpresM(:) = PRES00 * ( Rtot(iM) * rP0 * GsqrtRhotM(:) / Gsqrt_M(:) )**( CPtot(iM) / CVtot(iM) ) &
      !           - Phyd_M(:)
      ! dpresP(:) = PRES00 * ( Rtot(iP) * rP0 * GsqrtRhotP(:) / Gsqrt_P(:) )**( CPtot(iP) / CVtot(iP) ) &
      !           - Phyd_P(:)
      dpresM(:) = DPRES_(iM)
      dpresP(:) = DPRES_(iP)

      alpha(:) = max( sqrt( Gnn_M(:) * gamm * ( Phyd_M(:) + dpresM(:) ) * Gsqrt_M(:) / GsqrtDensM(:) ) + abs(VelM(:)), &
                      sqrt( Gnn_P(:) * gamm * ( Phyd_P(:) + dpresP(:) ) * Gsqrt_P(:) / GsqrtDensP(:) ) + abs(VelP(:))  )
      
      del_flux(:,ke,DENS_VID) = 0.5_RP * ( &
                    ( GsqrtDensP(:) * VelP(:) - GsqrtDensM(:) * VelM(:) )  &
                    - alpha(:) * ( GsqrtDDENS_P(:) - GsqrtDDENS_M(:) )     )

      del_flux(:,ke,MOMX_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMX_P(:) * VelP(:) - GsqrtMOMX_M(:) * VelM(:) )           &
                    + (  Gsqrt_P(:) * ( nx(:,ke) + G13_P(:) * nz(:,ke)) * dpresP(:)   &
                       - Gsqrt_M(:) * ( nx(:,ke) + G13_M(:) * nz(:,ke)) * dpresM(:) ) &
                    - alpha(:) * ( GsqrtMOMX_P(:) - GsqrtMOMX_M(:) )                  )

      del_flux(:,ke,MOMY_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMY_P(:) * VelP(:) - GsqrtMOMY_M(:) * VelM(:) ) &
                    + (  Gsqrt_P(:) * ( ny(:,ke) + G23_P(:) * nz(:,ke)) * dpresP(:)   &
                       - Gsqrt_M(:) * ( ny(:,ke) + G23_M(:) * nz(:,ke)) * dpresM(:) ) &
                    - alpha(:) * ( GsqrtMOMY_P(:) - GsqrtMOMY_M(:) )        )

      del_flux(:,ke,MOMZ_VID ) = 0.5_RP * ( &
                    ( GsqrtMOMZ_P(:) * VelP(:) - GsqrtMOMZ_M(:) * VelM(:) ) &
                    + (  Gsqrt_P(:) * dpresP(:) / GsqrtV_P(:)               &
                       - Gsqrt_M(:) * dpresM(:) / GsqrtV_M(:) ) * nz(:,ke)  &
                    - alpha(:) * ( GsqrtMOMZ_P(:) - GsqrtMOMZ_M(:) )        )
                    
      del_flux(:,ke,RHOT_VID) = 0.5_RP * ( &
                    ( GsqrtRhotP(:) * VelP(:) - GsqrtRhotM(:) * VelM(:) )   &
                    - alpha(:) * ( GsqrtDRHOT_P(:) - GsqrtDRHOT_M(:) )      )

      del_flux_hyd(:,ke,1) = 0.5_RP * ( &
          GsqrtV_P(:) * ( nx(:,ke) + G13_P(:) * nz(:,ke) ) * Phyd_P(:) &
        - GsqrtV_M(:) * ( nx(:,ke) + G13_M(:) * nz(:,ke) ) * Phyd_M(:) )
      
      del_flux_hyd(:,ke,2) = 0.5_RP * ( &
          GsqrtV_P(:) * ( ny(:,ke) + G23_P(:) * nz(:,ke) ) * Phyd_P(:) &
        - GsqrtV_M(:) * ( ny(:,ke) + G23_M(:) * nz(:,ke) ) * Phyd_M(:) )
    end do

    return
  end subroutine get_ebnd_flux

end module scale_atmos_dyn_tstep_short_dgm_heve
