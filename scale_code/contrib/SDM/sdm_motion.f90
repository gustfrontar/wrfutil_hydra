!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Module to control the motion of super-droplets (advection, sedmentation, precipitation)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!! @par History
!! @li      2014-07-11 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-11 (S.Shima) [rev] Fixed the bug of interpolation in sdm_getvel
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitoring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2017-02-08 (S.Shima) [mod] sdm_getvz is renamed to sdm_getvz_liq
!! @li      2017-02-08 (S.Shima) [add] sdm_getvz_ice
!! @li      2017-09-12 (S.Shima) [add] New formula by S.Shima to evaluate the area of the particle projected to the flow direction is introduced
!! @li      2017-09-18 (S.Shima) [fix] Apply the new formula of the particle area also to coalescence vz
!! @li      2017-10-02 (S.Shima) [mod] turbulence correction (Mitchell and Heymsfield 2005) is applied
!! @li      2017-10-02 (S.Shima) [mod] var_k formula is modfied to exp(-1.5phi)
!! @li      2017-10-15 (S.Shima) [add] avoid negative Reynolds number of ice. This can happen for unnatural, small and very distorted ice.
!! @li      2017-10-26 (S.Shima) [mod] now var_k_coef is defined in sdm_common.f90
!! @li      2017-11-11 (S.Shima) [mod] terminal velocity formula of Heymsfield and Westbrook (2010) is separated as a subroutine
!! @li      2017-11-11 (S.Shima) [add] terminal velocity formula of Bohm (1992)
!! @li      2018-02-28 (S.Shima) [add] simple-linear interpolation of wind velocity to preserve the divergence (Grabowski et al. 2018)
!! @li      2018-03-31 (S.Shima) [add] set an upper-limit 7mm to the droplet diameter when using Beard formula (1976)
!! @li      2019-01-10 (S.Shima) [mod] no interpolation of scalar variables when evaluating terminal velocity
!! @li      2019-01-12 (S.Shima) [mod] dry air density -> moist air density
!! @li      2020-03-23 (S.Shima) [fix] terminal velocity formula of Bohm (1992)
!! @li      2022-05-21 (S.Shima) [mod] support terrain following and vertically stretched coordinate
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module m_sdm_motion
  use scale_precision
  use scale_io

  implicit none
  private
  public :: sdm_getvz_ice, sdm_getvz_liq, sdm_getvel, sdm_move

  interface sdm_getvel

     module procedure sdm_getvel_simple_linear
!     module procedure sdm_getvel_trilinear

  end interface

contains
  subroutine sdm_move(IA,IS,IE,JA,JS,JE,sdm_dtadv,                      &
       sd_num,sd_u,sd_v,sd_vz,sd_x,sd_y,sd_ri,sd_rj,sd_rk)
    use scale_atmos_grid_cartesC, only: &
         GRID_CX => ATMOS_GRID_CARTESC_CX, &
         GRID_CY => ATMOS_GRID_CARTESC_CY
    use scale_atmos_grid_cartesC_real, only: &
         REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use m_sdm_common, only: &
         VALID2INVALID
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: sdm_dtadv ! time step of {motion of super-droplets} process
    integer, intent(in) :: sd_num ! number of super-droplets
    real(RP), intent(in) :: sd_u(1:sd_num) ! x-direction velocity of super-droplets
    real(RP), intent(in) :: sd_v(1:sd_num) ! y-direction velocity of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num) ! z velocity of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(inout) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(out)   :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
    real(RP), intent(out)   :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    ! Work variables
    integer :: n           ! index

    integer :: ic,icp,icm,jc,jcp,jcm,kc
    real(RP) :: RCDZ           !### Jacobian_33 = 1.0/(pz(i,j,k)/pk) (reciprocal of the center grid cell)
    real(RP) :: jacob31_33     !### Jacobian_31/Jacobian_33 = (pz(i,j,k)/pi)/(dx/di)
    real(RP) :: pz_d2i, dx_d2i
    real(RP) :: jacob32_33     !### Jacobian_32/Jacobian_33 = (pz(i,j,k)/pj)/(dy/dj)
    real(RP) :: pz_d2j, dy_d2j

    !---------------------------------------------------------------------
    ! The advection process of Super-Droplets.
    ! now only support uniform grid

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_move",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(ic,jc,kc,RCDZ) &
    !$omp& private(icp,icm,pz_d2i,dx_d2i,jacob31_33) &
    !$omp& private(jcp,jcm,pz_d2j,dy_d2j,jacob32_33)
    do n=1,sd_num

       !### skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle

       !### get the center grid index of super-droplets ###!
       ic = floor(sd_ri(n))+1
       jc = floor(sd_rj(n))+1
       kc = floor(sd_rk(n))+1

       !### Jacobian_33 = 1.0/(pz(i,j,k)/pk)
       !### reciprocal of the center grid cell
       RCDZ = 1.0_RP / (REAL_FZ(kc,ic,jc) - REAL_FZ(kc-1,ic,jc))

       !### Jacobian_31/Jacobian_33 = (pz(i,j,k)/pi)/(dx/di)
       icp = min(ic+1,IA)
       icm = max(ic-1, 1)
       pz_d2i = REAL_FZ(kc,icp,jc) - REAL_FZ(kc,icm,jc)
       dx_d2i = GRID_CX(   icp   ) - GRID_CX(   icm   )
       jacob31_33 = - pz_d2i/dx_d2i

       !### Jacobian_32/Jacobian_33 = (pz(i,j,k)/pj)/(dy/dj)
       jcp = min(jc+1,JA)
       jcm = max(jc-1, 1)
       pz_d2j = REAL_FZ(kc,ic,jcp) - REAL_FZ(kc,ic,jcm)
       dy_d2j = GRID_CY(      jcp) - GRID_CY(      jcm)
       jacob32_33 = - pz_d2j/dy_d2j

       !### move super-droplets (explicit) ###!
       sd_x(n)  = sd_x(n)  + sd_u(n)  * real(sdm_dtadv,kind=RP)
       sd_y(n)  = sd_y(n)  + sd_v(n)  * real(sdm_dtadv,kind=RP)
       sd_rk(n) = sd_rk(n) + (   jacob31_33*sd_u(n)    &
                               + jacob32_33*sd_v(n)    &
                               + sd_vz(n)           )  &
                             * RCDZ                    &
                             * real(sdm_dtadv,kind=RP)
    enddo
    !$omp end parallel do simd

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_move",1,1)
#endif
    return
  end subroutine sdm_move
  !----------------------------------------------------------------------------
  subroutine sdm_getvel_trilinear(    &
       KA,KS,KE, IA,IS,IE, JA,JS,JE, &
       u_scale,v_scale,w_scale,              &
       sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vzw)
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID

    ! Input variables
    integer, intent(in) :: KA,KS,KE
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: u_scale(KA,IA,JA) ! x-component velocity
    real(RP), intent(in) :: v_scale(KA,IA,JA) ! y-component velocity
    real(RP), intent(in) :: w_scale(KA,IA,JA) ! z-component velocity
    integer, intent(in) :: sd_num  ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! face index[k/real] of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_vzw(1:sd_num)  ! terminal velocity [in] / z velocity [out] of super-droplets
    ! Output variables
    real(RP), intent(out) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
    real(RP), intent(out) :: sd_u(1:sd_num) ! x-direction velocity of super-droplets
    real(RP), intent(out) :: sd_v(1:sd_num) ! y-direction velocity of super-droplets
    ! Work variables
    real(RP) :: sd_vz ! terminal velocity of super-droplets
    real(RP) :: sd_w  ! z velocity of super-droplets
    real(RP) :: ri    ! real index [i] of super-droplets
    real(RP) :: rj    ! real index [j] of super-droplets
    real(RP) :: rk    ! real index [k] of super-droplets
    real(RP) :: sXm   ! variable for inteporation
    real(RP) :: sXp   ! variable for inteporation
    real(RP) :: sYm   ! variable for inteporation
    real(RP) :: sYp   ! variable for inteporation
    real(RP) :: sZm   ! variable for inteporation
    real(RP) :: sZp   ! variable for inteporation
    integer :: iXm    ! index for inteporation
    integer :: iXp    ! index for inteporation
    integer :: iYm    ! index for inteporation
    integer :: iYp    ! index for inteporation
    integer :: iZm    ! index for inteporation
    integer :: iZp    ! index for inteporation
    integer :: n      ! index
    !-------------------------------------------------------------------
    ! The advection process of Super-Droplets.

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_getvel",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(ri,rj,rk) &
    !$omp& private(iXm,iXp,sXm,sXp) &
    !$omp& private(iYm,iYp,sYm,sYp) &
    !$omp& private(iZm,iZp,sZm,sZp) &
    !$omp& private(sd_w,sd_vz)
    do n=1,sd_num

       !### skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle

       !### get the real index of super-droplets ###!

       ri = sd_ri(n)
       rj = sd_rj(n)
       rk = sd_rk(n)

       !### x-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = min(max(floor(ri), 1), IA-1)
       iXp = iXm + 1
       sXm = min(max(ri - real(iXm,kind=RP), 0.0_RP), 1.0_RP)
       sXp = 1.0_RP - sXm

       iYm = min(max(floor(rj+0.5_RP), 1), JA-1)
       iYp = iYm + 1
       sYm = min(max((rj+0.5_RP) - real(iYm,kind=RP), 0.0_RP), 1.0_RP)
       sYp = 1.0_RP - sYm

       iZm = min(max(floor(rk+0.5_RP), KS), KE-1)
       iZp = iZm + 1
       sZm = min(max((rk+0.5_RP) - real(iZm,kind=RP), 0.0_RP), 1.0_RP)
       sZp = 1.0_RP - sZm

       sd_u(n) = u_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )             &
            + u_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )             &
            + u_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )             &
            + u_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )             &
            + u_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )             &
            + u_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )             &
            + u_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )             &
            + u_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       !### y-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = min(max(floor(ri+0.5_RP), 1), IA-1)
       iXp = iXm + 1
       sXm = min(max((ri+0.5_RP) - real(iXm,kind=RP), 0.0_RP), 1.0_RP)
       sXp = 1.0_RP - sXm

       iYm = min(max(floor(rj), 1), JA-1)
       iYp = iYm + 1
       sYm = min(max(rj - real(iYm,kind=RP), 0.0_RP), 1.0_RP)
       sYp = 1.0_RP - sYm

       iZm = min(max(floor(rk+0.5_RP), KS), KE-1)
       iZp = iZm + 1
       sZm = min(max((rk+0.5_RP) - real(iZm,kind=RP), 0.0_RP), 1.0_RP)
       sZp = 1.0_RP - sZm

       sd_v(n) = v_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )             &
            + v_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )             &
            + v_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )             &
            + v_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )             &
            + v_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )             &
            + v_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )             &
            + v_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )             &
            + v_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       !### z velocity ###!
       !### at the position of super-droplets         ###!
       !! interpolation in face grid
       iXm = min(max(floor(ri+0.5_RP), 1), IA-1)
       iXp = iXm + 1
       sXm = min(max((ri+0.5_RP) - real(iXm,kind=RP), 0.0_RP), 1.0_RP)
       sXp = 1.0_RP - sXm

       iYm = min(max(floor(rj+0.5_RP), 1), JA-1)
       iYp = iYm + 1
       sYm = min(max((rj+0.5_RP) - real(iYm,kind=RP), 0.0_RP), 1.0_RP)
       sYp = 1.0_RP - sYm

       iZm = min(max(floor(rk), KS-1), KE-1)
       iZp = iZm + 1
       sZm = min(max(rk - real(iZm,kind=RP), 0.0_RP), 1.0_RP)
       sZp = 1.0_RP - sZm

       sd_w = w_scale(iZm,iXm,iYm) * ( SXp * SYp * SZp )              &
            + w_scale(iZm,iXp,iYm) * ( SXm * SYp * SZp )              &
            + w_scale(iZm,iXm,iYp) * ( SXp * SYm * SZp )              &
            + w_scale(iZm,iXp,iYp) * ( SXm * SYm * SZp )              &
            + w_scale(iZp,iXm,iYm) * ( SXp * SYp * SZm )              &
            + w_scale(iZp,iXp,iYm) * ( SXm * SYp * SZm )              &
            + w_scale(iZp,iXm,iYp) * ( SXp * SYm * SZm )              &
            + w_scale(iZp,iXp,iYp) * ( SXm * SYm * SZm )

       ! the impact of terminal veloicty

       sd_vz = sd_vzw(n)

       sd_vzw(n) = sd_w - sd_vz
    end do
    !$omp end parallel do simd

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvel",1,1)
#endif
    return
  end subroutine sdm_getvel_trilinear
  !----------------------------------------------------------------------------
  subroutine sdm_getvel_simple_linear(    &
       KA,KS,KE, IA,IS,IE, JA,JS,JE, &
       u_scale,v_scale,w_scale,              &
       sd_num,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_u,sd_v,sd_vzw)
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID

    ! Input variables
    integer, intent(in) :: KA,KS,KE
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: u_scale(KA,IA,JA) ! x-component velocity
    real(RP), intent(in) :: v_scale(KA,IA,JA) ! y-component velocity
    real(RP), intent(in) :: w_scale(KA,IA,JA) ! z-component velocity
    integer, intent(in) :: sd_num  ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! face index[k/real] of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_vzw(1:sd_num)  ! terminal velocity [in] / z velocity [out] of super-droplets
    ! Output variables
    real(RP), intent(out) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets
    real(RP), intent(out) :: sd_u(1:sd_num) ! x-direction velocity of super-droplets
    real(RP), intent(out) :: sd_v(1:sd_num) ! y-direction velocity of super-droplets
    ! Work variables
    real(RP) :: sd_vz ! terminal velocity of super-droplets
    real(RP) :: sd_w  ! z velocity of super-droplets
    real(RP) :: ri    ! real index [i] of super-droplets
    real(RP) :: rj    ! real index [j] of super-droplets
    real(RP) :: rk    ! real index [k] of super-droplets
    real(RP) :: sXm   ! variable for inteporation
    real(RP) :: sXp   ! variable for inteporation
    real(RP) :: sYm   ! variable for inteporation
    real(RP) :: sYp   ! variable for inteporation
    real(RP) :: sZm   ! variable for inteporation
    real(RP) :: sZp   ! variable for inteporation
    integer :: iXm    ! index for inteporation
    integer :: iXp    ! index for inteporation
    integer :: iYm    ! index for inteporation
    integer :: iYp    ! index for inteporation
    integer :: iZm    ! index for inteporation
    integer :: iZp    ! index for inteporation
    integer :: n      ! index
    !-------------------------------------------------------------------
    ! The advection process of Super-Droplets.

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_getvel",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(ri,rj,rk) &
    !$omp& private(iXm,iXp,sXm,sXp) &
    !$omp& private(iYm,iYp,sYm,sYp) &
    !$omp& private(iZm,iZp,sZm,sZp) &
    !$omp& private(sd_w,sd_vz)
    do n=1,sd_num

       !### skip invalid super-droplets ###!

       if( sd_rk(n)<VALID2INVALID ) cycle

       !### get the real index of super-droplets ###!

       ri = sd_ri(n)
       rj = sd_rj(n)
       rk = sd_rk(n)

       !### x-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = min(max(floor(ri), 1), IA-1)
       iXp = iXm + 1
       sXm = min(max(ri - real(iXm,kind=RP), 0.0_RP), 1.0_RP)
       sXp = 1.0_RP - sXm

       iYm = min(max(nint(rj+0.5_RP), 1), JA)

       iZm = min(max(nint(rk+0.5_RP), KS), KE)

       sd_u(n) = u_scale(iZm,iXm,iYm) * SXp &
               + u_scale(iZm,iXp,iYm) * SXm

       !### y-velocity at the position of super-droplets ###!
       !! interpolation in face grid
       iXm = min(max(nint(ri+0.5_RP), 1), IA)

       iYm = min(max(floor(rj), 1), JA-1)
       iYp = iYm + 1
       sYm = min(max(rj - real(iYm,kind=RP), 0.0_RP), 1.0_RP)
       sYp = 1.0_RP - sYm

       iZm = min(max(nint(rk+0.5_RP), KS), KE)

       sd_v(n) = v_scale(iZm,iXm,iYm) * SYp &
               + v_scale(iZm,iXm,iYp) * SYm

       !### z velocity ###!
       !### at the position of super-droplets         ###!
       !! interpolation in face grid
       iXm = min(max(nint(ri+0.5_RP), 1), IA)

       iYm = min(max(nint(rj+0.5_RP), 1), JA)

       iZm = min(max(floor(rk), KS-1), KE-1)
       iZp = iZm + 1
       sZm = min(max(rk - real(iZm,kind=RP), 0.0_RP), 1.0_RP)
       sZp = 1.0_RP - sZm

       sd_w = w_scale(iZm,iXm,iYm) * SZp &
            + w_scale(iZp,iXm,iYm) * SZm

       ! the impact of terminal veloicty

       sd_vz = sd_vzw(n)

       sd_vzw(n) = sd_w - sd_vz
    end do
    !$omp end parallel do simd

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvel",1,1)
#endif
    return
  end subroutine sdm_getvel_simple_linear
  !----------------------------------------------------------------------------
  subroutine sdm_getvz_liq(    &
       KA, IA,IS,IE, JA,JS,JE, &
       pres,rhom,temp,            &
       sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sd_r,sd_vz,   &
       ilist_s,ilist_m,ilist_l,ptype        )
    ! evaluate the terminal velocities
    use scale_const, only: &
         rhow_mks => CONST_DWATR, &
         p0 => CONST_PRE00, &
         t0 => CONST_TEM00, &
         cp => CONST_CPdry, &
         rd => CONST_Rdry
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID, O_SIX, vz_b, vz_c, dosedimentation, i2, STAT_LIQ
    use scale_prc, only:  &
         PRC_abort

    ! Input variables
    integer, intent(in) :: KA
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: pres(KA,IA,JA)  ! Pressure
    real(RP), intent(in) :: rhom(KA,IA,JA)  ! moist air density
    real(RP), intent(in) :: temp(KA,IA,JA)  ! temperature
    integer, intent(in) :: sd_num           ! number of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
    real(RP), intent(in) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)  ! equivalent radius of super-droplets
    character(*), intent(in) :: ptype   ! process type : 'motion process' or 'stochastic coalescence process'
    ! Output variables
    real(RP), intent(inout) :: sd_vz(1:sd_num)! terminal velocity of super-droplets in real space
    integer, intent(out) :: ilist_s(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_m(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_l(1:sd_num)  ! buffer for list vectorization

    ! Local parameters
    real(RP), parameter :: grav = 9.80665E+2_RP       ! gravity acceleration [cm/s^2]
    real(RP), parameter :: rhow = 1.0_RP              ! density of water [g/cm^3]
    real(RP), parameter :: lb4l = 6.62E-6_RP          ! base length[cm] for mean free path  of air molecules
    real(RP), parameter :: pb4l = 1013.25_RP          ! base pressure[hPa] for mean free path of air molecules
    real(RP), parameter :: tb4l = 293.15_RP           ! base temperarure[20C] for mean free path of air molecules
    real(RP), parameter :: visb4l = 1.818E-4_RP       ! base dynamic viscosity[g/(cm*s)] for mean free path of air molecules
    real(RP), parameter :: m2cm = 1.E+2_RP            ! converter of unit [m] => [cm]
    real(RP), parameter :: cm2m = 1.E-2_RP            ! converter of unit [cm] => [m]
    real(RP), parameter :: pa2hpa = 1.E-2_RP          ! converter of unit [Pa] => [hPa]
    real(RP), parameter :: kgm2gcm = 1.E-3_RP         ! converter of unit [kg/m^3] => [g/cm^3]

    ! Work variables
    real(RP) :: sd_dia    ! diameter[cm] of super-droplets
    real(RP) :: sd_l      ! mean free path[cm] of air molecules at the location of super-droplets
    real(RP) :: sd_p      ! pressure[hPa] at the location of super-droplets
    real(RP) :: sd_rd     ! air density[g/cm^3] at the location of super-droplets
    real(RP) :: sd_sigma  ! water tension[N/m*10^3=g/s^2] of super-droplets
    real(RP) :: sd_t      ! temperature[K] at the location of super-droplets
    real(RP) :: sd_visco  ! dynamic viscosity[g/(cm*s)] at the location of super-droplets
    real(RP) :: bond      ! modified Bond number
    real(RP) :: csc       ! slip correction factor
    real(RP) :: gxdrow    ! = grav * (rhow-sd_rd)
    real(RP) :: nda       ! Davies number
    real(RP) :: npp       ! physical property number
    real(RP) :: nre       ! Reynolds number
    real(RP) :: p0iv      ! inverse of p0
    real(RP) :: rddvcp    ! rd / cp

    real(RP) :: c1        ! temporary
    real(RP) :: c2        ! temporary
    real(RP) :: c3        ! temporary
    real(RP) :: x1        ! temporary
    real(RP) :: x2        ! temporary
    real(RP) :: x3        ! temporary
    real(RP) :: y         ! temporary
    real(RP) :: tau       ! temporary
    real(RP) :: tdeg      ! temporary

    real(RP) :: ri        ! real index [i] of super-droplets
    real(RP) :: rj        ! real index [j] of super-droplets
    real(RP) :: rk        ! real index [k] of super-droplets

    real(RP) :: sXm       ! variable for inteporation
    real(RP) :: sXp       ! variable for inteporation
    real(RP) :: sYm       ! variable for inteporation
    real(RP) :: sYp       ! variable for inteporation
    real(RP) :: sZm       ! variable for inteporation
    real(RP) :: sZp       ! variable for inteporation

    integer :: iXm        ! index for inteporation
    integer :: iXp        ! index for inteporation
    integer :: iYm        ! index for inteporation
    integer :: iYp        ! index for inteporation
    integer :: iZm        ! index for inteporation
    integer :: iZp        ! index for inteporation

    integer :: tlist_s         ! total list number for small
    integer :: tlist_m         ! total list number for middle
    integer :: tlist_l         ! total list number for large

    integer :: scnt            ! counter for small
    integer :: mcnt            ! counter for middle
    integer :: lcnt            ! counter for large

    integer :: i, j, k, n, m   ! index
    !---------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_getvz_liq",1,1)
#endif
    if( .not. dosedimentation ) then
       sd_vz(:) = 0.0_RP
       return
    endif

    if(( ptype .ne. 'interpolation' ) .and. ( ptype .ne. 'no_interpolation' ))then
       LOG_ERROR("sdm_getvz_liq",*) 'last argument of sdm_getvz_liq has to be "no_interpolation" or "interpolation"!'
       call PRC_abort
    end if

    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize

    tlist_s = 0
    tlist_m = 0
    tlist_l = 0

    p0iv   = 1.0_RP/real(p0*pa2hpa,kind=RP)
    rddvcp = real(rd,kind=RP)/real(cp,kind=RP)

    ! Get index list for compressing buffer.
    scnt = 0
    mcnt = 0
    lcnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n) .ne. STAT_LIQ ) cycle

       ri = sd_ri(n)
       rj = sd_rj(n)

       i = floor(ri)+1
       j = floor(rj)+1

       if ( i < IS .or. i > IE .or. j < JS .or. j > JE ) then ! in the boundary
          sd_vz(n) = 0.0_RP
          cycle
       end if

       sd_dia = 2.0_RP * sd_r(n) * m2cm   !! diameter [m=>cm]

       if( sd_dia<=1.9E-3_RP ) then

          !== small cloud droplets ==!

          scnt = scnt + 1
          ilist_s(scnt) = n

       else if( sd_dia>1.9E-3_RP .and. sd_dia<=1.07E-1_RP ) then

          !== large cloud droplets and small raindrops ==!

          mcnt = mcnt + 1
          ilist_m(mcnt) = n

       else

          !== raindrops ==!

          lcnt = lcnt + 1
          ilist_l(lcnt) = n

       end if

    end do

    tlist_s = scnt
    tlist_m = mcnt
    tlist_l = lcnt

    ! Calculate terminal velocity of super-droplets

    if( ptype .eq. 'interpolation' ) then

       !$omp  parallel 

       !### small cloud droplets ###!

       if( tlist_s>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(iXm,iXp,sXm,sXp) &
          !$omp& private(iYm,iYp,sYm,sYp) &
          !$omp& private(iZm,iZp,sZm,sZp) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,sd_l,c1,csc)
          do m=1,tlist_s
             n = ilist_s(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhom(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhom(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhom(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhom(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhom(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhom(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhom(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhom(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                  + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * sqrt(sd_t/tb4l)

             !== terminal vecocity ==!

             c1  = gxdrow / ( 18.0_RP * sd_visco )
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             sd_vz(n) = c1 * csc * sd_dia * sd_dia

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do nowait

       end if

       !### large cloud droplets and small raindrops ###!

       if( tlist_m>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(iXm,iXp,sXm,sXp) &
          !$omp& private(iYm,iYp,sYm,sYp) &
          !$omp& private(iZm,iZp,sZm,sZp) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,sd_l,c2,nda,csc) &
          !$omp& private(x1,x2,x3,y,nre)
          do m=1,tlist_m
             n = ilist_m(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!
             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhom(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhom(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhom(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhom(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhom(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhom(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhom(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhom(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                           + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * sqrt(sd_t/tb4l)

             !== Davies number ==!

             c2 = sd_rd * (4.0_RP*gxdrow)/(3.0_RP*sd_visco*sd_visco)

             nda = c2 * sd_dia * sd_dia * sd_dia

             !== Reynolds number ==!

             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             x1 = log(nda)
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_b(1)                                       &
                  + vz_b(2) * (x1)                                &
                  + vz_b(3) * (x2)                                &
                  + vz_b(4) * (x3)                                &
                  + vz_b(5) * (x2*x2)                             &
                  + vz_b(6) * (x3*x2)                             &
                  + vz_b(7) * (x3*x3)

             nre = csc * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do nowait

       end if

       !### raindrops ###!

       if( tlist_l>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(iXm,iXp,sXm,sXp) &
          !$omp& private(iYm,iYp,sYm,sYp) &
          !$omp& private(iZm,iZp,sZm,sZp) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,tau,sd_sigma,c3,bond,npp) &
          !$omp& private(x1,x2,x3,y,nre)
          do m=1,tlist_l
             n = ilist_l(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm
             sd_dia = min(sd_dia,7.0E-1_RP)     !! Assume d=7mm if the droplet is too large. Beard formula (1976) is only applicable to d<7mm

             !== density of air ==!

             sd_rd = rhom(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhom(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhom(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhom(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhom(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhom(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhom(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhom(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_rd = sd_rd * kgm2gcm
             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p  = pres(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + pres(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + pres(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + pres(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + pres(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + pres(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + pres(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + pres(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             sd_p  = sd_p * pa2hpa

             !== temperarure ==!

             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                  + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== surface tension [N/m*10^3=g/s^2]   ==!
             !==               (Holten et al. 2005) ==!

             tau = 1.0_RP - sd_t/647.0960_RP
             sd_sigma = 0.2358_RP * exp( 1.2560_RP*log(tau) )      &
                  * ( 1.0_RP - 0.6250_RP*tau )

             if( sd_t<267.5_RP ) then

                tau = tanh( (sd_t-243.9_RP)/35.35_RP )

                sd_sigma = sd_sigma - 2.854E-3_RP * tau + 1.666E-3_RP

             end if

             sd_sigma = sd_sigma * 1.E+3_RP  !! N/m => g/s^2

             !== modified Bond number ==!

             c3 = (4.0_RP*gxdrow)/(3.0_RP*sd_sigma)

             bond = c3 * sd_dia * sd_dia

             !== physical property number ==!

             npp = (sd_rd*sd_rd) * (sd_sigma*sd_sigma*sd_sigma) &
                  / (gxdrow*sd_visco*sd_visco*sd_visco*sd_visco)

             npp = exp( O_SIX * log(npp) )

             !== Reynolds number ==!

             x1 = log( bond * npp )
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_c(1)                                       &
                  + vz_c(2) * (x1)                                &
                  + vz_c(3) * (x2)                                &
                  + vz_c(4) * (x3)                                &
                  + vz_c(5) * (x2*x2)                             &
                  + vz_c(6) * (x3*x2)

             nre = npp * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do

       end if

       !$omp end parallel 

    else if( ptype .eq. 'no_interpolation' ) then

       !$omp parallel 

       !### small cloud droplets ###!

       if( tlist_s>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(i,j,k) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,sd_l,c1,csc)
          do m=1,tlist_s
             n = ilist_s(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhom(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * sqrt(sd_t/tb4l)

             !== terminal vecocity ==!

             c1  = gxdrow / ( 18.0_RP * sd_visco )
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             sd_vz(n) = c1 * csc * sd_dia * sd_dia

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do nowait

       end if

       !### large cloud droplets and small raindrops ###!

       if( tlist_m>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(i,j,k) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,sd_l,c2,nda,csc) &
          !$omp& private(x1,x2,x3,y,nre)
          do m=1,tlist_m
             n = ilist_m(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm

             !== density of air ==!

             sd_rd = rhom(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== mean free path of air molecules ==!

             sd_l = lb4l * (sd_visco/visb4l) * (pb4l/sd_p)      &
                  * sqrt(sd_t/tb4l)

             !== Davies number ==!

             c2 = sd_rd * (4.0_RP*gxdrow)/(3.0_RP*sd_visco*sd_visco)

             nda = c2 * sd_dia * sd_dia * sd_dia

             !== Reynolds number ==!
             csc = 1.0_RP + 2.510_RP * ( sd_l / sd_dia )

             x1 = log(nda)
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_b(1)                                       &
                  + vz_b(2) * (x1)                                &
                  + vz_b(3) * (x2)                                &
                  + vz_b(4) * (x3)                                &
                  + vz_b(5) * (x2*x2)                             &
                  + vz_b(6) * (x3*x2)                             &
                  + vz_b(7) * (x3*x3)

             nre = csc * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do nowait

       end if

       !### raindrops ###!

       if( tlist_l>0 ) then

!OCL NORECURRENCE
          !$omp  do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(i,j,k) &
          !$omp& private(sd_dia,sd_rd,gxdrow,sd_p,sd_t) &
          !$omp& private(tdeg,sd_visco,tau,sd_sigma,c3,bond,npp) &
          !$omp& private(x1,x2,x3,y,nre)
          do m=1,tlist_l
             n = ilist_l(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             !== diameter[cm] of droplets ==!

             sd_dia = 2.0_RP * sd_r(n) * m2cm
             sd_dia = min(sd_dia,7.0E-1_RP)     !! Assume d=7mm if the droplet is too large. Beard formula (1976) is only applicable to d<7mm

             !== density of air ==!

             sd_rd = rhom(k,i,j) * kgm2gcm

             gxdrow = grav * ( rhow - sd_rd )

             !== pressure ==!

             sd_p = pres(k,i,j) * pa2hpa

             !== temperarure ==!

             sd_t = temp(k,i,j)

             !== dynamic viscosity (Pruppacher & Klett,1997) ==!

             tdeg = sd_t - t0     !! [K] => [degC]

             if( tdeg>=0.0_RP ) then

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-4_RP

             else

                sd_visco = ( 1.7180_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-4_RP
             end if

             !== surface tension [N/m*10^3=g/s^2]   ==!
             !==               (Holten et al. 2005) ==!

             tau = 1.0_RP - sd_t/647.0960_RP

             sd_sigma = 0.23580_RP * exp( 1.2560_RP*log(tau) )      &
                  * ( 1.00_RP - 0.6250_RP*tau )

             if( sd_t<267.50_RP ) then

                tau = tanh( (sd_t-243.90_RP)/35.350_RP )

                sd_sigma = sd_sigma - 2.854E-3_RP * tau + 1.666E-3_RP

             end if

             sd_sigma = sd_sigma * 1.E+3_RP  !! N/m => g/s^2

             !== modified Bond number ==!

             c3 = (4.0_RP*gxdrow)/(3.0_RP*sd_sigma)

             bond = c3 * sd_dia * sd_dia

             !== physical property number ==!

             npp = (sd_rd*sd_rd) * (sd_sigma*sd_sigma*sd_sigma) &
                  / (gxdrow*sd_visco*sd_visco*sd_visco*sd_visco)

             npp = exp( O_SIX * log(npp) )

             !== Reynolds number ==!

             x1 = log( bond * npp )
             x2 = x1 * x1
             x3 = x1 * x2

             y  = vz_c(1)                                       &
                  + vz_c(2) * (x1)                                &
                  + vz_c(3) * (x2)                                &
                  + vz_c(4) * (x3)                                &
                  + vz_c(5) * (x2*x2)                             &
                  + vz_c(6) * (x3*x2)

             nre = npp * exp(y)

             !== terminal vecocity ==!

             sd_vz(n) = (sd_visco*nre)/(sd_rd*sd_dia)

             sd_vz(n) = sd_vz(n) * cm2m

          end do
          !$omp end do

       end if

       !$omp end parallel 

    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvz_liq",1,1)
#endif
    return
  end subroutine sdm_getvz_liq
  !----------------------------------------------------------------------------
  subroutine sdm_getvz_ice(    &
       KA, IA,IS,IE, JA,JS,JE, &
       rhom,temp,            &
       sd_num,sd_liqice,sd_x,sd_y,sd_ri,sd_rj,sd_rk,sdi,sd_vz,   &
       ilist,ptype)
    ! evaluate the terminal velocities of ice particles
    ! The formula of Heymsfield and Westbrook (2010) is used
    use scale_const, only: &
         rhoi_mks => CONST_DICE, &   ! density of ice [kg/m^3]
         t0 => CONST_TEM00           ! 0 degC in [K]
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         VALID2INVALID, dosedimentation, i2, STAT_ICE, sdicedef, &
         F_THRD, ONE_PI, var_k_coef
    use scale_prc, only:  &
         PRC_abort

    ! Input variables
    integer, intent(in) :: KA
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: rhom(KA,IA,JA)  ! moist air density [kg/m^3]
    real(RP), intent(in) :: temp(KA,IA,JA)  ! temperature [K]
    integer, intent(in) :: sd_num           ! number of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
    real(RP), intent(in) :: sd_x(1:sd_num)  ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)  ! y-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    type(sdicedef), intent(in) :: sdi       ! ice phase super-droplets
    character(*), intent(in) :: ptype   ! process type : 'motion process' or 'stochastic coalescence process'
    ! Output variables
    real(RP), intent(inout) :: sd_vz(1:sd_num)! terminal velocity of super-droplets in real space (overwrite only for ice perticles)
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization

    ! Local parameters
    real(RP), parameter :: C0 = 0.35_RP      ! Fitting parameter of Re(X*)
    real(RP), parameter :: delta0 = 8.0_RP   ! Fitting parameter of Re(X*)

    ! Work variables
    real(RP) :: sd_maxD   ! maximum dimension of the particle [m]
    real(RP) :: sd_dvisc  ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
    real(RP) :: sd_rd     ! air density [kg/m^3] at the location of super-droplets
    real(RP) :: sd_t      ! temperature [K] at the location of super-droplets
    real(RP) :: sd_mass   ! mass [kg] of the particle
    real(RP) :: sd_area   ! A: area of the particle projected to the flow direction [m^2]
    real(RP) :: sd_area_ratio ! sd_area_ratio Ar: A divided by the area of a circumscribing disk
    real(RP) :: sd_area_ratio_q ! q: A divided by the area of a circumscribed ellipse
    real(RP) :: sd_phi    ! sd_rp/sd_re: aspect ratio
    real(RP) :: var_k     ! k:=exp(-c/a). Needed to evaluate A.

    real(RP) :: tdeg      ! temporary

    real(RP) :: ri        ! real index [i] of super-droplets
    real(RP) :: rj        ! real index [j] of super-droplets
    real(RP) :: rk        ! real index [k] of super-droplets

    integer :: tlist      ! total list number for ice
    integer :: icnt       ! counter for ice

    real(RP) :: sXm       ! variable for inteporation
    real(RP) :: sXp       ! variable for inteporation
    real(RP) :: sYm       ! variable for inteporation
    real(RP) :: sYp       ! variable for inteporation
    real(RP) :: sZm       ! variable for inteporation
    real(RP) :: sZp       ! variable for inteporation

    integer :: iXm        ! index for inteporation
    integer :: iXp        ! index for inteporation
    integer :: iYm        ! index for inteporation
    integer :: iYp        ! index for inteporation
    integer :: iZm        ! index for inteporation
    integer :: iZp        ! index for inteporation

    integer :: i, j, k, n, m   ! index
   !---------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_getvz_ice",1,1)
#endif
    if( .not. dosedimentation ) then
       sd_vz(:) = 0.0_RP
       return
    endif

    if(( ptype .ne. 'interpolation' ) .and. ( ptype .ne. 'no_interpolation' ))then
       LOG_ERROR("sdm_getvz_ice",*) 'last argument of sdm_getvz_ice has to be "no_interpolation" or "interpolation"!'
       call PRC_abort
    end if

    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    tlist = 0

    ! Get index list for compressing buffer.
    icnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n) .ne. STAT_ICE ) cycle

       ri = sd_ri(n)
       rj = sd_rj(n)

       i = floor(ri)+1
       j = floor(rj)+1

       if ( i < IS .or. i > IE .or. j < JS .or. j > JE ) then ! in the boundary
          sd_vz(n) = 0.0_RP
          cycle
       end if

       !== ice particles ==!
       icnt = icnt + 1
       ilist(icnt) = n

    end do

    tlist = icnt

    ! Calculate terminal velocity of ice super-droplets

    if( ptype .eq. 'interpolation' ) then

       if( tlist>0 ) then
!OCL NORECURRENCE
          !$omp  parallel do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(iXm,iXp,sXm,sXp) &
          !$omp& private(iYm,iYp,sYm,sYp) &
          !$omp& private(iZm,iZp,sZm,sZp) &
          !$omp& private(sd_rd,sd_t,sd_maxD) &
          !$omp& private(tdeg,sd_dvisc,sd_mass,sd_phi) &
          !$omp& private(var_k,sd_area,sd_area_ratio,sd_area_ratio_q)
          do m=1,tlist
             n = ilist(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! Interpolation in certer grid
             iXm = floor(ri+0.5_RP)
             iXp = iXm + 1
             sXm = (ri+0.5_RP) - real(iXm,kind=RP)
             sXp = 1.0_RP - sXm

             iYm = floor(rj+0.5_RP)
             iYp = iYm + 1
             sYm = (rj+0.5_RP) - real(iYm,kind=RP)
             sYp = 1.0_RP - sYm

             iZm = floor(rk+0.50_RP)
             iZp = iZm + 1
             sZm = (rk+0.5_RP) - real(iZm,kind=RP)
             sZp = 1.0_RP - sZm

             ! moist air density [kg/m^3] at the location of super-droplets
             sd_rd = rhom(iZm,iXm,iYm) * ( SXp * SYp * SZp )    &
                  + rhom(iZm,iXp,iYm) * ( SXm * SYp * SZp )    &
                  + rhom(iZm,iXm,iYp) * ( SXp * SYm * SZp )    &
                  + rhom(iZm,iXp,iYp) * ( SXm * SYm * SZp )    &
                  + rhom(iZp,iXm,iYm) * ( SXp * SYp * SZm )    &
                  + rhom(iZp,iXp,iYm) * ( SXm * SYp * SZm )    &
                  + rhom(iZp,iXm,iYp) * ( SXp * SYm * SZm )    &
                  + rhom(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             ! temperature[K] at the location of super-droplets
             sd_t =  temp(iZm,iXm,iYm) * ( SXp * SYp * SZp )                        &
                  + temp(iZm,iXp,iYm) * ( SXm * SYp * SZp )                        &
                  + temp(iZm,iXm,iYp) * ( SXp * SYm * SZp )                        &
                  + temp(iZm,iXp,iYp) * ( SXm * SYm * SZp )                        &
                  + temp(iZp,iXm,iYm) * ( SXp * SYp * SZm )                        &
                  + temp(iZp,iXp,iYm) * ( SXm * SYp * SZm )                        &
                  + temp(iZp,iXm,iYp) * ( SXp * SYm * SZm )                        &
                  + temp(iZp,iXp,iYp) * ( SXm * SYm * SZm )

             ! maximum dimension of the particle [m]
             sd_maxD = 2.0_RP * max(sdi%re(n),sdi%rp(n))

             ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
             !== (Pruppacher & Klett,1997) ==!
             tdeg = sd_t - t0     !! [K] => [degC]
             if( tdeg>=0.0_RP ) then
                sd_dvisc = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                sd_dvisc = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             ! mass of the particle [kg]
             sd_mass =  F_THRD * ONE_PI * sdi%re(n)**2 * sdi%rp(n) * sdi%rho(n)

             ! aspect ratio
             sd_phi = sdi%rp(n)/sdi%re(n)

             ! A: area of the particle projected to the flow direction [m^2]
!             var_k = exp( -sdi%rp(n)/sdi%re(n) )
             var_k = exp( -var_k_coef*sd_phi)
             sd_area = ONE_PI * sdi%re(n) * (sd_maxD/2.0_RP) * (sdi%rho(n)/rhoi_mks)**var_k

             ! Ar: A divided by the area of a circumscribing disk
             sd_area_ratio = sd_area / (ONE_PI*sd_maxD**2/4.0_RP)

!!$             ! terminal velocity of ice particles by Heymsfield and Westbrook (2010)
!!$             call terminal_vel_Heymsfield_Westbrook_2010(sd_vz(n),sd_rd,sd_dvisc,sd_mass,sd_area_ratio,sd_maxD)

             ! terminal velocity of ice particles by Bohm (1992)
             !! q: A divided by the area of a circumscribed ellipse
             sd_area_ratio_q = sd_area / (ONE_PI * sdi%re(n) * (sd_maxD/2.0_RP))
             call terminal_vel_Bohm_1992(sd_vz(n),sd_rd,sd_dvisc,sd_mass,sd_area_ratio_q,2.0_RP*sdi%re(n),sd_phi)

          end do
          !$omp end parallel do

       end if

    else if( ptype .eq. 'no_interpolation' ) then

       if( tlist>0 ) then
!OCL NORECURRENCE
          !$omp  parallel do schedule(static,256) &
          !$omp& private(m,n,ri,rj,rk) &
          !$omp& private(i,j,k) &
          !$omp& private(sd_rd,sd_t,sd_maxD) &
          !$omp& private(tdeg,sd_dvisc,sd_mass,sd_phi) &
          !$omp& private(var_k,sd_area,sd_area_ratio,sd_area_ratio_q)
          do m=1,tlist
             n = ilist(m)

             ri = sd_ri(n)
             rj = sd_rj(n)
             rk = sd_rk(n)

             !! coversion to center grid index
             !! No interpolation
             i = floor(ri)+1
             j = floor(rj)+1
             k = floor(rk)+1

             ! air density [kg/m^3] at the location of super-droplets
             sd_rd = rhom(k,i,j)

             ! temperature[K] at the location of super-droplets
             sd_t = temp(k,i,j)

             ! maximum dimension of the particle [m]
             sd_maxD = 2.0_RP * max(sdi%re(n),sdi%rp(n))

             ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
             !== (Pruppacher & Klett,1997) ==!
             tdeg = sd_t - t0     !! [K] => [degC]
             if( tdeg>=0.0_RP ) then
                sd_dvisc = ( 1.7180_RP + 4.9E-3_RP*tdeg ) * 1.E-5_RP
             else
                sd_dvisc = ( 1.718_RP + 4.9E-3_RP*tdeg              &
                     - 1.2E-5_RP*tdeg*tdeg ) * 1.E-5_RP
             end if

             ! mass of the particle [kg]
             sd_mass =  F_THRD * ONE_PI * sdi%re(n)**2 * sdi%rp(n) * sdi%rho(n)

             ! aspect ratio
             sd_phi = sdi%rp(n)/sdi%re(n)

             ! A: area of the particle projected to the flow direction [m^2]
!             var_k = exp( -sdi%rp(n)/sdi%re(n) )
             var_k = exp( -var_k_coef*sd_phi)
             sd_area = ONE_PI * sdi%re(n) * (sd_maxD/2.0_RP) * (sdi%rho(n)/rhoi_mks)**var_k

             ! Ar: A divided by the area of a circumscribing disk
             sd_area_ratio = sd_area / (ONE_PI*sd_maxD**2/4.0_RP)

!!$             ! terminal velocity of ice particles by Heymsfield and Westbrook (2010)
!!$             call terminal_vel_Heymsfield_Westbrook_2010(sd_vz(n),sd_rd,sd_dvisc,sd_mass,sd_area_ratio,sd_maxD)

             ! terminal velocity of ice particles by Bohm (1992)
             !! q: A divided by the area of a circumscribed ellipse
             sd_area_ratio_q = sd_area / (ONE_PI * sdi%re(n) * (sd_maxD/2.0_RP))
             call terminal_vel_Bohm_1992(sd_vz(n),sd_rd,sd_dvisc,sd_mass,sd_area_ratio_q,2.0_RP*sdi%re(n),sd_phi)

          end do
          !$omp end parallel do
       end if
    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_getvz_ice",1,1)
#endif
    return
  end subroutine sdm_getvz_ice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine terminal_vel_Heymsfield_Westbrook_2010(sd_vz,sd_rd,sd_dvisc,sd_mass,sd_area_ratio,sd_maxD)
    ! terminal velocity of ice particles by Heymsfield and Westbrook (2010)
    use scale_const, only: &
         grav_mks => CONST_GRAV  ! gravitational constant [m/s^2]
    use m_sdm_common, only: &
         ONE_PI

    implicit none

    real(RP),intent(out) :: sd_vz ! terminal velocity [m/s]

    real(RP),intent(in) :: sd_rd  ! air density [kg/m^3] at the location of super-droplets
    real(RP),intent(in) :: sd_dvisc  ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
    real(RP),intent(in) :: sd_mass   ! mass [kg] of the particle
    real(RP),intent(in) :: sd_area_ratio ! sd_area_ratio Ar: A divided by the area of a circumscribing disk
    real(RP),intent(in) :: sd_maxD   ! maximum dimension of the particle [m]

    real(RP) :: n_modX    ! modified Best (Davies) number
    real(RP) :: nre       ! Reynolds number

    ! Local parameters
    real(RP), parameter :: C0 = 0.35_RP      ! Fitting parameter of Re(X*)
    real(RP), parameter :: delta0 = 8.0_RP   ! Fitting parameter of Re(X*)

    ! modified Best (Davies) number of the particle
    n_modX = (sd_rd/sd_dvisc**2)*(8.0_RP*sd_mass*grav_mks)/ONE_PI/sqrt(sd_area_ratio)

    ! Reynolds number of the particle
    nre = (delta0**2 / 4.0_RP) * (                            &
         &  sqrt(1.0_RP + 4.0_RP*sqrt(n_modX/C0)/(delta0**2)) &
         &  - 1.0_RP                                          &
         &                        )**2

!!$             ! turbulence correction (Mitchell and Heymsfield 2005)
!!$             if(nre>=1.0d3) then
!!$                n_orgX = n_modX/sqrt(sd_area_ratio)
!!$                nre = nre - 1.7e-3_RP*min(n_orgX,1.0d8)**0.8_RP
!!$             end if
!!$             !! Avoid negative Reynolds number. This can happen for small and very distorted ice.
!!$             !! This is unnatural, and has to be examined.
!!$             nre = max(nre,0.0_RP)

!!$             !! this would be better?
!!$             nre = nre - 1.7e-3_RP*min(n_modX,1.0e8_RP)**0.8_RP

    ! terminal velocity [m/s]
    sd_vz = sd_dvisc * nre / sd_rd / sd_maxD

  end subroutine terminal_vel_Heymsfield_Westbrook_2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine terminal_vel_Bohm_1992(sd_vz,sd_rd,sd_dvisc,sd_mass,sd_area_ratio_q,sd_chard,sd_phi)
    ! terminal velocity of ice particles by Bohm (1989,1992,1999)
    use scale_const, only: &
         grav_mks => CONST_GRAV  ! gravitational constant [m/s^2]
    use m_sdm_common, only: &
         ONE_PI

    implicit none

    real(RP),intent(out) :: sd_vz ! terminal velocity [m/s]

    real(RP),intent(in) :: sd_rd  ! air density [kg/m^3] at the location of super-droplets
    real(RP),intent(in) :: sd_dvisc  ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
    real(RP),intent(in) :: sd_mass   ! mass [kg] of the particle
    real(RP),intent(in) :: sd_area_ratio_q ! q: A divided by the area of a circumscribed ellipse
    real(RP),intent(in) :: sd_chard  ! 2.0*(equatorial radius), NOT the maximum dimension of the particle [m]
    real(RP),intent(in) :: sd_phi    ! sd_rp/sd_re: aspect ratio

    real(RP) :: n_X    ! Best (Davies) number
    real(RP) :: sf_k   ! viscous shape factor
    real(RP) :: drag_gamma ! Gamma for pressure drag coefficient
    real(RP) :: C_DPS,C_DP  ! pressure drag coefficients
    real(RP) :: C_DO        ! Oseen drag coefficients
    real(RP) :: var_beta, var_gamma ! auxiliary variables
    real(RP) :: nre       ! Reynolds number

    real(RP),parameter :: X0 = 2.8d6 ! parameter for turbulence correction

    ! Best (Davies) number of the particle
    n_X = (sd_rd/sd_dvisc**2)*(8.0_RP*sd_mass*grav_mks)/ONE_PI/max(sd_phi,1.0_RP)/max(sd_area_ratio_q**0.25_RP,sd_area_ratio_q)

    ! turbulence correction
    n_X = n_X*(1.0_RP+(n_X/X0)**2)/(1.0_RP+1.6_RP*(n_X/X0)**2)

    ! viscous shape factor
    sf_k = min( &
         & max(0.82_RP+0.18_RP*sd_phi,0.85_RP),      &
         & 0.37_RP+0.63_RP/sqrt(sd_phi),            &
         & 1.33_RP/(max(log(sd_phi),0.0_RP)+1.19_RP) &
         &    )

    ! Gamma for pressure drag coefficient
    drag_gamma = max(1.0_RP, min(1.98_RP, 3.76_RP - 8.41_RP*sd_phi + 9.18_RP*sd_phi**2 - 3.53_RP*sd_phi**3) )

    ! pressure drag coefficients
    C_DPS = max( 0.292_RP*sf_k*drag_gamma, 0.492_RP-0.2_RP/sqrt(sd_phi) )
    C_DP  = max( 1.0_RP, sd_area_ratio_q*(1.46_RP*sd_area_ratio_q-0.46_RP) )*C_DPS
    !!! NOTE
    !!! C_DP  = max( 1.0_RP, (1.46_RP*sd_area_ratio_q-0.46_RP) )*C_DPS, could be the correct equation.

    ! Oseen drag coefficients
    C_DO = 4.5_RP*sf_k**2*max(sd_phi,1.0_RP)

    ! auxiliary variables
    var_beta = sqrt(1.0_RP + (C_DP/(6.0_RP*sf_k))*sqrt(n_X/C_DP)) - 1.0_RP
    var_gamma = (C_DO-C_DP)/(4.0_RP*C_DP)

    ! Reynolds number
    nre = (6.0_RP*sf_k/C_DP)*var_beta**2
    !! matching to low Reynolds number regime
    nre = nre*( 1.0_RP + 2.0_RP*var_beta*exp(-var_beta*var_gamma)/(2.0_RP+var_beta)/(1.0_RP+var_beta))

    ! terminal velocity [m/s]
    sd_vz = sd_dvisc * nre / sd_rd / sd_chard

  end subroutine terminal_vel_Bohm_1992
end module m_sdm_motion
