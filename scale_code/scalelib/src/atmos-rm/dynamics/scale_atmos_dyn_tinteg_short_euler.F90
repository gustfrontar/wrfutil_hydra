!-------------------------------------------------------------------------------
!> module Atmosphere / Dyn Tinteg
!!
!! @par Description
!!          Temporal integration in Dynamical core for Atmospheric process
!!          Euler scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_short_euler
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
#if defined DEBUG || defined QUICKDEBUG
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
  public :: ATMOS_DYN_Tinteg_short_euler_setup
  public :: ATMOS_DYN_Tinteg_short_euler_finalize
  public :: ATMOS_DYN_Tinteg_short_euler
  public :: ATMOS_DYN_Tinteg_short_euler_split_explicit

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

  ! for communication

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_DYN_Tinteg_short_euler_setup( &
       tinteg_type )
    use scale_prc, only: &
       PRC_abort
    use scale_const, only: &
       UNDEF => CONST_UNDEF
    use scale_comm_cartesC, only: &
       COMM_vars8_init
    implicit none

    character(len=*), intent(in) :: tinteg_type
    !---------------------------------------------------------------------------
    return
  end subroutine ATMOS_DYN_Tinteg_short_euler_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_DYN_Tinteg_short_euler_finalize
    implicit none
    return
  end subroutine ATMOS_DYN_Tinteg_short_euler_finalize
  !-----------------------------------------------------------------------------
  !> Euler
  subroutine ATMOS_DYN_Tinteg_short_euler( &
       DENS, MOMZ, MOMX, MOMY, RHOT, PROG,      &
       mflx_hi,  tflx_hi,                       &
       DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,  &
       DPRES0, RT2P, CORIOLI,                   &
       num_diff, wdamp_coef, divdmp_coef, DDIV, &
       FLAG_FCT_MOMENTUM, FLAG_FCT_T,           &
       FLAG_FCT_ALONG_STREAM,                   &
       CDZ, FDZ, FDX, FDY,                      &
       RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,      &
       PHI, GSQRT, J13G, J23G, J33G, MAPF,      &
       REF_pres, REF_dens,                      &
       BND_W, BND_E, BND_S, BND_N, TwoD,        &
       dt                                       )
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait
    use scale_atmos_dyn_tstep_short, only: &
       ATMOS_DYN_tstep => ATMOS_DYN_Tstep_short
    use scale_atmos_dyn_common, only: &
       ATMOS_DYN_Copy_boundary
    implicit none

    real(RP), intent(inout) :: DENS(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)
    real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3)
    real(RP), intent(out  ) :: tflx_hi(KA,IA,JA,3)

    real(RP), intent(in)    :: DENS_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMZ_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMX_t(KA,IA,JA)
    real(RP), intent(in)    :: MOMY_t(KA,IA,JA)
    real(RP), intent(in)    :: RHOT_t(KA,IA,JA)

    real(RP), intent(in)    :: DPRES0(KA,IA,JA)
    real(RP), intent(in)    :: RT2P(KA,IA,JA)
    real(RP), intent(in)    :: CORIOLI(IA,JA)
    real(RP), intent(in)    :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)    :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef
    real(RP), intent(in)    :: DDIV(KA,IA,JA)

    logical,  intent(in)    :: FLAG_FCT_MOMENTUM
    logical,  intent(in)    :: FLAG_FCT_T
    logical,  intent(in)    :: FLAG_FCT_ALONG_STREAM

    real(RP), intent(in)    :: CDZ (KA)
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

    real(RP), intent(in)    :: REF_pres(KA,IA,JA)   !< reference pressure
    real(RP), intent(in)    :: REF_dens(KA,IA,JA)

    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical,  intent(in)    :: TwoD

    real(RP), intent(in)    :: dt

    real(RP) :: DENS0(KA,IA,JA)
    real(RP) :: MOMZ0(KA,IA,JA)
    real(RP) :: MOMX0(KA,IA,JA)
    real(RP) :: MOMY0(KA,IA,JA)
    real(RP) :: RHOT0(KA,IA,JA)
    real(RP) :: PROG0(KA,IA,JA,VA)

    real(RP) :: dtrk

    integer  :: i, j, k, iv, n
    !---------------------------------------------------------------------------

    call PROF_rapstart("DYN_Euler_Prep",3)

#if defined QUICKDEBUG || defined DEBUG
    mflx_hi(   1:KS-1,:,:,:) = UNDEF
    mflx_hi(KE+1:KA  ,:,:,:) = UNDEF
#endif

!OCL XFILL
    DENS0 = DENS
!OCL XFILL
    MOMZ0 = MOMZ
!OCL XFILL
    MOMX0 = MOMX
!OCL XFILL
    MOMY0 = MOMY
!OCL XFILL
    RHOT0 = RHOT
!OCL XFILL
    if ( VA > 0 ) PROG0 = PROG
    call PROF_rapend  ("DYN_Euler_Prep",3)

    !------------------------------------------------------------------------
    ! Start RK
    !------------------------------------------------------------------------

    !##### Euler : PROG0,PROG->PROG_RK1 #####

    call PROF_rapstart("DYN_Euler",3)

    dtrk = dt

    call ATMOS_DYN_tstep( DENS, MOMZ, MOMX, MOMY, RHOT,                     & ! [OUT]
                          PROG,                                             & ! [OUT]
                          mflx_hi, tflx_hi,                                 & ! [INOUT,OUT]
                          DENS0,    MOMZ0,    MOMX0,    MOMY0,    RHOT0,    & ! [IN]
                          DENS,     MOMZ,     MOMX,     MOMY,     RHOT,     & ! [IN]
                          DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! [IN]
                          PROG0, PROG,                                      & ! [IN]
                          DPRES0, RT2P, CORIOLI,                            & ! [IN]
                          num_diff, wdamp_coef, divdmp_coef, DDIV,          & ! [IN]
                          FLAG_FCT_MOMENTUM, FLAG_FCT_T,                    & ! [IN]
                          FLAG_FCT_ALONG_STREAM,                            & ! [IN]
                          CDZ, FDZ, FDX, FDY,                               & ! [IN]
                          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! [IN]
                          PHI, GSQRT, J13G, J23G, J33G, MAPF,               & ! [IN]
                          REF_pres, REF_dens,                               & ! [IN]
                          BND_W, BND_E, BND_S, BND_N, TwoD,                 & ! [IN]
                          dtrk, .true.                                      ) ! [IN]

    call PROF_rapend  ("DYN_Euler",3)

    return
  end subroutine ATMOS_DYN_Tinteg_short_euler

  !-----------------------------------------------------------------------------
  !> Euler (Split-Explicit)
  subroutine ATMOS_DYN_Tinteg_short_euler_split_explicit( &
    DDENS, DMOMZ, DMOMX, DMOMY, DRHOT,                & ! (inout)
    mflx_hi,  tflx_hi,                                & ! (out)
    DENS0, RHOT0,                                     & ! (in)
    DENS_t,   MOMZ_t,   MOMX_t,   MOMY_t,   RHOT_t,   & ! (in)
    DENS_t_se, MOMZ_t_se, MOMX_t_se, MOMY_t_se, RHOT_t_se, &  ! (in)
    RT2P,                                             & ! (in)
    num_diff, wdamp_coef, divdmp_coef, DDIV,          & ! (in)
    CDZ, FDZ, FDX, FDY,                               & ! (in)
    RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,               & ! (in)
    GSQRT, J13G, J23G, J33G, MAPF,                    & ! (in)
    BND_W, BND_E, BND_S, BND_N, TwoD,                 & ! (in)
    dtrk, last                                        ) ! (in)

    use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait
    use scale_atmos_dyn_tstep_short, only: &
      ATMOS_DYN_tstep => ATMOS_DYN_Tstep_short_split_explicit
    use scale_atmos_dyn_common, only: &
      ATMOS_DYN_Copy_boundary
    implicit none

    real(RP), intent(inout) :: DDENS(KA,IA,JA)   ! prognostic variables
    real(RP), intent(inout) :: DMOMZ(KA,IA,JA)   !
    real(RP), intent(inout) :: DMOMX(KA,IA,JA)   !
    real(RP), intent(inout) :: DMOMY(KA,IA,JA)   !
    real(RP), intent(inout) :: DRHOT(KA,IA,JA)   !

    real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3) ! mass flux
    real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3) ! internal energy flux

    real(RP), intent(in) :: DENS0(KA,IA,JA)
    real(RP), intent(in) :: RHOT0(KA,IA,JA)

    real(RP), intent(in)  :: DENS_t(KA,IA,JA)    ! tendency
    real(RP), intent(in)  :: MOMZ_t(KA,IA,JA)    !
    real(RP), intent(in)  :: MOMX_t(KA,IA,JA)    !
    real(RP), intent(in)  :: MOMY_t(KA,IA,JA)    !
    real(RP), intent(in)  :: RHOT_t(KA,IA,JA)    !

    real(RP), intent(in) :: DENS_t_se(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_t_se(KA,IA,JA)
    real(RP), intent(in) :: MOMX_t_se(KA,IA,JA)
    real(RP), intent(in) :: MOMY_t_se(KA,IA,JA)
    real(RP), intent(in) :: RHOT_t_se(KA,IA,JA)       

    real(RP), intent(in)  :: RT2P    (KA,IA,JA)
    real(RP), intent(in) :: num_diff(KA,IA,JA,5,3)
    real(RP), intent(in)  :: wdamp_coef(KA)
    real(RP), intent(in)    :: divdmp_coef
    real(RP), intent(in)    :: DDIV(KA,IA,JA)

    real(RP), intent(in)  :: CDZ (KA)
    real(RP), intent(in)  :: FDZ (KA-1)
    real(RP), intent(in)  :: FDX (IA-1)
    real(RP), intent(in)  :: FDY (JA-1)
    real(RP), intent(in)  :: RCDZ(KA)
    real(RP), intent(in)  :: RCDX(IA)
    real(RP), intent(in)  :: RCDY(JA)
    real(RP), intent(in)  :: RFDZ(KA-1)
    real(RP), intent(in)  :: RFDX(IA-1)
    real(RP), intent(in)  :: RFDY(JA-1)

    real(RP), intent(in)  :: GSQRT   (KA,IA,JA,7) !< vertical metrics {G}^1/2
    real(RP), intent(in)  :: J13G    (KA,IA,JA,7) !< (1,3) element of Jacobian matrix
    real(RP), intent(in)  :: J23G    (KA,IA,JA,7) !< (2,3) element of Jacobian matrix
    real(RP), intent(in)  :: J33G                 !< (3,3) element of Jacobian matrix
    real(RP), intent(in)  :: MAPF    (IA,JA,2,4)  !< map factor

    logical,  intent(in)  :: BND_W
    logical,  intent(in)  :: BND_E
    logical,  intent(in)  :: BND_S
    logical,  intent(in)  :: BND_N
    logical,  intent(in)  :: TwoD

    real(RP), intent(in)  :: dtrk
    logical,  intent(in)  :: last

    !----------------------------------------------------------

#if defined QUICKDEBUG || defined DEBUG
    mflx_hi(   1:KS-1,:,:,:) = UNDEF
    mflx_hi(KE+1:KA  ,:,:,:) = UNDEF
#endif

    call PROF_rapstart("DYN_SE_FBEuler",3)

    call ATMOS_DYN_tstep( &
      DDENS, DMOMZ, DMOMX, DMOMY, DRHOT,                      &
      mflx_hi, tflx_hi,                                       &
      DENS0, RHOT0,                                           & 
      DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,                 &
      DENS_t_se, MOMZ_t_se, MOMX_t_se, MOMY_t_se, RHOT_t_se,  &
      RT2P, num_diff, wdamp_coef, divdmp_coef, DDIV,          &
      CDZ, FDZ, FDX, FDY, RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY, &
      GSQRT, J13G, J23G, J33G, MAPF,                          &
      BND_W, BND_E, BND_S, BND_N, TwoD,                       &
      dtrk, last                                              )

     call PROF_rapend("DYN_SE_FBEuler",3)

     return
   end subroutine ATMOS_DYN_Tinteg_short_euler_split_explicit

end module scale_atmos_dyn_tinteg_short_euler
