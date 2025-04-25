!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics Temporal integration
!!
!! @par Description
!!          Temporal integration scheme selecter for dynamical short time step
!!
!! @author Team SCALE
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_tinteg_short
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
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_Tinteg_short_setup

  !> Runge-Kutta loop
  abstract interface
     subroutine short( &
          DENS, MOMZ, MOMX, MOMY, RHOT, PROG,       & ! (inout)
          mflx_hi, tflx_hi,                         & ! (inout, out)
          DENS_t, MOMZ_t, MOMX_t, MOMY_t, RHOT_t,   & ! (in)
          DPRES0, RT2P, CORIOLI,                    & ! (in)
          num_diff, wdamp_coef, divdmp_coef, DDIV,  & ! (in)
          FLAG_FCT_MOMENTUM, FLAG_FCT_T,            & ! (in)
          FLAG_FCT_ALONG_STREAM,                    & ! (in)
          CDZ, FDZ, FDX, FDY,                       & ! (in)
          RCDZ, RCDX, RCDY, RFDZ, RFDX, RFDY,       & ! (in)
          PHI, GSQRT, J13G, J23G, J33G, MAPF,       & ! (in)
          REF_pres, REF_dens,                       & ! (in)
          BND_W, BND_E, BND_S, BND_N, TwoD,         & ! (in)
          dt                                        ) ! (in)
       use scale_precision
       use scale_atmos_grid_cartesC_index
       use scale_index
       real(RP), intent(inout) :: DENS(KA,IA,JA)
       real(RP), intent(inout) :: MOMZ(KA,IA,JA)
       real(RP), intent(inout) :: MOMX(KA,IA,JA)
       real(RP), intent(inout) :: MOMY(KA,IA,JA)
       real(RP), intent(inout) :: RHOT(KA,IA,JA)
       real(RP), intent(inout) :: PROG(KA,IA,JA,VA)

       real(RP), intent(inout) :: mflx_hi(KA,IA,JA,3)
       real(RP), intent(out)   :: tflx_hi(KA,IA,JA,3)

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
     end subroutine short

     subroutine short_split_explicit( DDENS, DMOMZ, DMOMX, DMOMY, DRHOT, & ! (inout)
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
       use scale_precision
       use scale_atmos_grid_cartesC_index
       use scale_index
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
       real(RP), intent(in)  :: divdmp_coef
       real(RP), intent(in)  :: DDIV(KA,IA,JA)

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
     end subroutine short_split_explicit
          
     subroutine finalize
     end subroutine finalize
  end interface
  procedure(short), pointer :: ATMOS_DYN_Tinteg_short => NULL()
  public :: ATMOS_DYN_Tinteg_short
  procedure(short_split_explicit), pointer :: ATMOS_DYN_Tinteg_short_split_explicit => NULL()
  public :: ATMOS_DYN_Tinteg_short_split_explicit
  procedure(finalize), pointer :: ATMOS_DYN_Tinteg_short_finalize => NULL()
  public :: ATMOS_DYN_Tinteg_short_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: check_availability

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Register
  subroutine ATMOS_DYN_Tinteg_short_setup( &
       ATMOS_DYN_Tinteg_short_TYPE, ATMOS_DYN_Tstep_short_TYPE )

    use scale_precision
    use scale_atmos_grid_cartesC_index
    use scale_index
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_dyn_tinteg_short_euler, only: &
       ATMOS_DYN_Tinteg_short_euler_setup,    &
       ATMOS_DYN_Tinteg_short_euler_finalize, &
       ATMOS_DYN_Tinteg_short_euler,          &
       ATMOS_DYN_Tinteg_short_euler_split_explicit
    use scale_atmos_dyn_tinteg_short_rk3, only: &
       ATMOS_DYN_Tinteg_short_rk3_setup,    &
       ATMOS_DYN_Tinteg_short_rk3_finalize, &
       ATMOS_DYN_Tinteg_short_rk3
    use scale_atmos_dyn_tinteg_short_rk4, only: &
       ATMOS_DYN_Tinteg_short_rk4_setup,    &
       ATMOS_DYN_Tinteg_short_rk4_finalize, &
       ATMOS_DYN_Tinteg_short_rk4
    use scale_atmos_dyn_tinteg_short_rk7s6o, only: &
       ATMOS_DYN_Tinteg_short_rk7s6o_setup, &
       ATMOS_DYN_Tinteg_short_rk7s6o
    use scale_atmos_dyn_tinteg_short_rk11s8o, only: &
       ATMOS_DYN_Tinteg_short_rk11s8o_setup, &
       ATMOS_DYN_Tinteg_short_rk11s8o         
    implicit none

    character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_short_TYPE
    character(len=*), intent(in)  :: ATMOS_DYN_Tstep_short_TYPE

    logical :: is_SplitExplicit
    !---------------------------------------------------------------------------

    call check_availability( ATMOS_DYN_Tinteg_short_TYPE, ATMOS_DYN_Tstep_short_TYPE, &
      is_SplitExplicit  )

    ATMOS_DYN_Tinteg_short_finalize => DYN_Tinteg_short_finalize

    select case( ATMOS_DYN_Tinteg_short_TYPE )
    case( 'EULER' )
      call ATMOS_DYN_Tinteg_short_euler_setup( &
           ATMOS_DYN_Tinteg_short_TYPE )
      ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_euler
      if ( is_SplitExplicit ) ATMOS_DYN_Tinteg_short_split_explicit => ATMOS_DYN_Tinteg_short_euler_split_explicit
      ATMOS_DYN_Tinteg_short_finalize => ATMOS_DYN_Tinteg_short_euler_finalize
    case( 'RK3', 'RK3WS2002' )
       call ATMOS_DYN_Tinteg_short_rk3_setup( &
            ATMOS_DYN_Tinteg_short_TYPE )
       ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_rk3
       ATMOS_DYN_Tinteg_short_finalize => ATMOS_DYN_Tinteg_short_rk3_finalize
    case( 'RK4' )
       call ATMOS_DYN_Tinteg_short_rk4_setup( &
            ATMOS_DYN_Tinteg_short_TYPE )
       ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_rk4
       ATMOS_DYN_Tinteg_short_finalize => ATMOS_DYN_Tinteg_short_rk4_finalize
    case( 'RK7s6o', 'RK7s6oLawson1967', 'RK7s6oButcher1964' )
       call ATMOS_DYN_Tinteg_short_rk7s6o_setup( &
              ATMOS_DYN_Tinteg_short_TYPE )
         ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_rk7s6o
   case( 'RK11s8o', 'RK11s8oCooperVerner1972' )
      call ATMOS_DYN_Tinteg_short_rk11s8o_setup( &
               ATMOS_DYN_Tinteg_short_TYPE )
         ATMOS_DYN_Tinteg_short => ATMOS_DYN_Tinteg_short_rk11s8o 
    case( 'OFF', 'NONE' )
       ! do nothing
    case default
       LOG_ERROR("ATMOS_DYN_Tinteg_short_setup",*) 'ATMOS_DYN_TINTEG_SHORT_TYPE is invalid: ', ATMOS_DYN_Tinteg_short_TYPE
       call PRC_abort
    end select

    return
  end subroutine ATMOS_DYN_Tinteg_short_setup

  subroutine DYN_Tinteg_short_finalize

    return
  end subroutine DYN_Tinteg_short_finalize

!---- private
  subroutine check_availability( ATMOS_DYN_Tinteg_short_TYPE, ATMOS_DYN_Tstep_short_TYPE, &
   is_SplitExplicit )
   implicit none
   character(len=*), intent(in)  :: ATMOS_DYN_Tinteg_short_TYPE
   character(len=*), intent(in)  :: ATMOS_DYN_Tstep_short_TYPE
   logical, intent(out) :: is_SplitExplicit
   !----------------------------------------------

   ! RK7s6o, RK11s8o are available only for FVM-HEVE
   select case( ATMOS_DYN_Tinteg_short_TYPE )
   case( 'RK7s6o', 'RK7s6oLawson1967', 'RK7s6oButcher1964', 'RK11s8o', 'RK11s8oCooperVerner1972' )
      if ( .not. (ATMOS_DYN_Tstep_short_TYPE == 'HEVE' .or. ATMOS_DYN_Tstep_short_TYPE == 'FVM-HEVE') ) then
         LOG_ERROR("ATMOS_DYN_Tinteg_short_setup",*) "The specified ATMOS_DYN_TINTEG_SHORT_TYPE is now supported only for 'HEVE',", ATMOS_DYN_Tinteg_short_TYPE
      end if
   end select

   ! For Split-explicit integration, Euler scheme is only available
   if ( ATMOS_DYN_Tstep_short_TYPE == 'FVM-HEVI-SplitExplicit' .or. ATMOS_DYN_Tstep_short_TYPE == 'HEVI-SplitExplicit' ) then
      select case( ATMOS_DYN_Tinteg_short_TYPE )
      case ('EULER')
      case default
         LOG_ERROR("ATMOS_DYN_Tinteg_short_setup",*) "The specified ATMOS_DYN_TINTEG_SHORT_TYPE is not supported for 'HEVI-SplitExplicit',", ATMOS_DYN_Tinteg_short_TYPE
      end select
      is_SplitExplicit = .true.
   else
      is_SplitExplicit = .false.
   end if

   return
  end subroutine check_availability

end module scale_atmos_dyn_tinteg_short
