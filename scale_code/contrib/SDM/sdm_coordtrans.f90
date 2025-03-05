!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coordinate transformation of the SDM variables
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
!! @li      2014-06-26 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-06-27 (S.Shima) [rev] sdm_rk2z, sdm_x2ri, sdm_y2rj added
!! @li      2014-07-04 (S.Shima) [rev] sdm_x2ri, sdm_y2rj modified to reduce dependency.
!! @li      2014-07-05 (S.Shima) [rev] sdm_x2ri, sdm_y2rj optimized for uniform grid temporarily. Bugfix of sdm_getrklu.
!! @li      2018-02-28 (S.Shima) [mod] range check in sdm_x2ri and sdm_y2rj removed for predictor-corrector of SD motion
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module m_sdm_coordtrans

  implicit none
  private
  public :: sdm_getrklu, sdm_x2ri, sdm_y2rj, sdm_z2rk, sdm_rk2z

contains
  subroutine sdm_rk2z(  &
       IS, IE, JS, JE, &
       sd_num,sd_x,sd_y,sd_rk,sd_z,sd_ri,sd_rj)
    use scale_precision
    use scale_atmos_grid_cartesC_real, only: &
         REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
    use m_sdm_common, only: &
         VALID2INVALID, INVALID

    implicit none

    integer, intent(in) :: IS,IE
    integer, intent(in) :: JS,JE
    integer, intent(in) :: sd_num    ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! face index-k(real) of super-droplets
    real(RP), intent(out) :: sd_z(1:sd_num) ! z-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! face index-i(real) of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! face index-j(real) of super-droplets

    ! Work variables
    real(RP) :: ri      ! real face index [i] of super-droplets
    real(RP) :: rj      ! real face index [j] of super-droplets
    real(RP) :: sXm     ! variable for inteporation
    real(RP) :: sXp     ! variable for inteporation
    real(RP) :: sYm     ! variable for inteporation
    real(RP) :: sYp     ! variable for inteporation
    real(RP) :: zph_l   ! "zph_crs" at lower[k]
    real(RP) :: zph_u   ! "zph_crs" at upper[k+1]

    integer :: iXm           ! index for inteporation
    integer :: iXp           ! index for inteporation
    integer :: iYm           ! index for inteporation
    integer :: iYp           ! index for inteporation

    integer :: k, n, i, j    ! index

    !### get horizontal face index(real) of super-droplets ###!
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(k,ri,rj,iXm,iXp,sXm,sXp,iYm,iYp,sYm,sYp,zph_u,zph_l)   
    do n=1,sd_num

       if( sd_rk(n)>VALID2INVALID ) then

          !### get vertical face index of super-droplets ###!
          k  = floor(sd_rk(n))
          
          !### get horizontal face index (real) of super-droplets ###!
!!$          iloop: do i = IS, IE
!!$             if( sd_x(n) < GRID_FX(i) ) then
!!$                ri = real(i-1,kind=RP) + ( sd_x(n)-GRID_FX(i-1) ) / ( GRID_FX(i)-GRID_FX(i-1) )
!!$                exit iloop
!!$             endif
!!$          enddo iloop
!!$          jloop: do j = JS, JE
!!$             if( sd_y(n) < GRID_FY(j) ) then
!!$                rj = real(j-1,kind=RP) + ( sd_y(n)-GRID_FY(j-1) ) / ( GRID_FY(j)-GRID_FY(j-1) ) 
!!$                exit jloop
!!$             endif
!!$          enddo jloop
          ri = sd_ri(n)
          rj = sd_rj(n)

          ! Face to Center conversion
          iXm = floor(ri+0.5d0) 
          iXp = iXm + 1
          sXm = (ri+0.5d0) - real(iXm,kind=RP)
          sXp = 1.d0 - sXm

          ! Face to Center conversion
          iYm = floor(rj+0.5d0)
          iYp = iYm + 1
          sYm = (rj+0.5d0) - real(iYm,kind=RP)
          sYp = 1.d0 - sYm

          zph_u = real(REAL_FZ(k+1,iXm,iYm),kind=RP) * ( SXp * SYp )  &
               + real(REAL_FZ(k+1,iXp,iYm),kind=RP) * ( SXm * SYp )  &
               + real(REAL_FZ(k+1,iXm,iYp),kind=RP) * ( SXp * SYm )  &
               + real(REAL_FZ(k+1,iXp,iYp),kind=RP) * ( SXm * SYm )
          
          zph_l = real(REAL_FZ(k,iXm,iYm),kind=RP) * ( SXp * SYp )    &
               + real(REAL_FZ(k,iXp,iYm),kind=RP) * ( SXm * SYp )    &
               + real(REAL_FZ(k,iXm,iYp),kind=RP) * ( SXp * SYm )    &
               + real(REAL_FZ(k,iXp,iYp),kind=RP) * ( SXm * SYm )
          
          !### get z-coordinate of super-droplets ###!
          sd_z(n) = zph_l+(zph_u-zph_l)*(sd_rk(n) - real(k,kind=RP))

       else

          sd_z(n) = INVALID 
          
       end if

    end do
    !$omp end parallel do simd
  end subroutine sdm_rk2z
  !-----------------------------------------------------------------------------
  subroutine sdm_getrklu(  &
                         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                         sdm_zlower,sdm_zupper,             &
                         sd_rkl,sd_rku)

   use scale_atmos_grid_cartesC_real, only : &
        REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
   use scale_precision
   use m_sdm_common, only : &
        rkumax, rkumin, rklmax, knum_sdm

   integer, intent(in) :: KA, KS, KE
   integer, intent(in) :: IA, IS, IE
   integer, intent(in) :: JA, JS, JE
   real(RP),intent(in) :: sdm_zlower  ! sdm_zlower+surface height is the lower limitaion of SD's position
   real(RP),intent(in) :: sdm_zupper  ! Upper limitaion of SD's position
   real(RP),intent(out) :: sd_rkl(IA,JA)    ! index-k(real) at "sdm_zlower"
   real(RP),intent(out) :: sd_rku(IA,JA)    ! index-k(real) at "sdm_zupper"

   ! Work variables
   real(RP) :: zph_l         ! REAL_FZ(k,:,:)
   real(RP) :: zph_u         ! REAL_FZ(k+1,:,:)
   real(RP) :: zph_s         ! REAL_FZ(KS-1,:,:), i.e., surface height

   integer :: i, j, k        ! index
   !---------------------------------------------------------------------

      rkumax = 0.0_RP
      rkumin = real(KE+1,kind=RP) + 1.0_RP
      rklmax = 0.0_RP

      !$omp parallel
      !$omp do private(j,i) collapse(2)
      do j=1,JA
      do i=1,IA
         sd_rkl(i,j) = 0.0_RP
         sd_rku(i,j) = 0.0_RP
      end do
      end do
      !$omp end do

      ! Get vertical index[k/real] at 'sdm_zlower', 'sdm_zupper'
      !$omp do private(j,i,zph_u,zph_l,zph_s) collapse(2)
      do j=JS,JE
      do i=IS,IE
      do k=KS-1,KE

         zph_u = REAL_FZ(k+1,i,j)
         zph_l = REAL_FZ(k,i,j)
         zph_s = REAL_FZ(KS-1,i,j)

         !### get vertical index of 'sdm_zlower' ###!

         if( sdm_zlower<(zph_u-zph_s) .and.                             &
                       sdm_zlower>=(zph_l-zph_s) ) then

            sd_rkl(i,j) = real(k,kind=RP)                               &
                        + real((sdm_zlower+zph_s-zph_l)/(zph_u-zph_l),kind=RP)
         end if

         !### get vertical index of 'sdm_zupper' ###!

         if( (sdm_zupper>=zph_l) .and. (sdm_zupper<zph_u) ) then

            sd_rku(i,j) = real(k,kind=RP)                               &
                        + real((sdm_zupper-zph_l)/(zph_u-zph_l),kind=RP)

         end if

      end do
      end do
      end do
      !$omp end do 

     ! -----
     ! Copy data

      !$omp do simd
      do i=1,IA

         sd_rkl(i,1:(JS-1))  = sd_rkl(i,JS)
         sd_rkl(i,(JE+1):JA) = sd_rkl(i,JE)

         sd_rku(i,1:(JS-1))  = sd_rku(i,JS)
         sd_rku(i,(JE+1):JA) = sd_rku(i,JE)

      end do
      !$omp end do simd
      !$omp do simd
      do j=1, JA

         sd_rkl(1:(IS-1) ,j) = sd_rkl(IS,j)
         sd_rkl((IE+1):IA,j) = sd_rkl(IE,j)

         sd_rku(1:(IS-1) ,j) = sd_rku(IS,j)
         sd_rku((IE+1):IA,j) = sd_rku(IE,j)

      end do
      !$omp end do simd

     ! Get maximum grid number, maximum index and minimum index in vertical

      !$omp do private(j,i) reduction(max:rkumax) collapse(2)
      do j=JS,JE
      do i=IS,IE

         if( sd_rku(i,j)>rkumax ) then
            rkumax = sd_rku(i,j)
         end if

      end do
      end do
      !$omp end do

      !$omp do private(j,i) reduction(min:rkumin) collapse(2)
      do j=JS,JE
      do i=IS,IE

         if( sd_rku(i,j)<rkumin ) then
            rkumin = sd_rku(i,j)
         end if

      end do
      end do
      !$omp end do

      !$omp do private(j,i) reduction(max:rklmax) collapse(2)
      do j=JS,JE
      do i=IS,IE

         if( sd_rkl(i,j)>rklmax ) then
            rklmax = sd_rkl(i,j)
         end if

      end do
      end do
      !$omp end do

      !$omp end parallel

      knum_sdm = min( KE-KS+1, (floor(rkumax)+1)-KS+1 )

    return
  end subroutine sdm_getrklu
  !-----------------------------------------------------------------------------
  subroutine sdm_z2rk(KS,KE,IS,IE,JS,JE,sd_num,sd_x,sd_y,sd_z,sd_ri,sd_rj,sd_rk)
      use scale_atmos_grid_cartesC_real, only: &
        REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ
      use scale_precision
      use m_sdm_common, only: &
           INVALID,sdm_zlower,sdm_zupper

      integer, intent(in) :: KS,KE
      integer, intent(in) :: IS,IE
      integer, intent(in) :: JS,JE
      integer, intent(in) :: sd_num       ! number of super-droplets
      real(RP), intent(in) :: sd_x(1:sd_num)   ! x-coordinate of super-droplets
      real(RP), intent(in) :: sd_y(1:sd_num)   ! y-coordinate of super-droplets
      real(RP), intent(inout) :: sd_z(1:sd_num)! z-coordinate of super-droplets
      real(RP), intent(out) :: sd_ri(1:sd_num) ! index-i(real) of super-droplets
      real(RP), intent(out) :: sd_rj(1:sd_num) ! index-j(real) of super-droplets
      real(RP), intent(inout) :: sd_rk(1:sd_num) ! index-k(real) of super-droplets

      ! Work variables
      real(RP) :: ri      ! real index [i] of super-droplets
      real(RP) :: rj      ! real index [j] of super-droplets
      real(RP) :: sXm     ! variable for inteporation
      real(RP) :: sXp     ! variable for inteporation
      real(RP) :: sYm     ! variable for inteporation
      real(RP) :: sYp     ! variable for inteporation
      real(RP) :: zph_l   ! "zph_crs" at lower[k]
      real(RP) :: zph_u   ! "zph_crs" at upper[k+1]
      real(RP) :: zph_s   ! "zph_crs" at surface[k=2]

      integer :: iXm           ! index for inteporation
      integer :: iXp           ! index for inteporation
      integer :: iYm           ! index for inteporation
      integer :: iYp           ! index for inteporation

      integer :: k, n, i, j    ! index
   !-------------------------------------------------------------------
      ! -----
      !### get horizontal face index(real) of super-droplets ###!
      call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
      call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

      ! Get vertical index[k/real] of super-droplets
      !$omp  parallel do simd schedule(static,256) &
      !$omp& private(ri,rj,iXm,iXp,sXm,sXp,iYm,iYp,sYm,sYp,zph_s) &   
      !$omp& private(k,zph_u,zph_l)
      do n=1,sd_num

         !### initialize ###!

         sd_rk(n) = INVALID

         !### convert to invalid super-droplets ###!

         if( sd_z(n)>=real(sdm_zupper,kind=RP) ) then
            sd_z(n) = INVALID
            cycle
         end if

!         ri = sd_x(n) * real(dxiv_sdm,kind=RP) + 2.d0
!         rj = sd_y(n) * real(dyiv_sdm,kind=RP) + 2.d0
!!$         iloop: do i = IS, IE
!!$           if( sd_x(n) < GRID_FX(i) ) then
!!$            ri = real(i-1,kind=RP) + ( sd_x(n)-GRID_FX(i-1) ) / ( GRID_FX(i)-GRID_FX(i-1) )
!!$            exit iloop
!!$           endif
!!$         enddo iloop
         ri = sd_ri(n)
!!$         jloop: do j = JS, JE
!!$           if( sd_y(n) < GRID_FY(j) ) then
!!$            rj = real(j-1,kind=RP) + ( sd_y(n)-GRID_FY(j-1) ) / ( GRID_FY(j)-GRID_FY(j-1) ) 
!!$            exit jloop
!!$           endif
!!$         enddo jloop
         rj = sd_rj(n)

         ! Interpolation in center grid
!         iXm = floor(ri-0.5d0)
         iXm = floor(ri+0.5d0) ! Face to Center conversion
         iXp = iXm + 1
!         sXm = (ri-0.5d0) - real(iXm,kind=RP)
         sXm = (ri+0.5d0) - real(iXm,kind=RP)
         sXp = 1.d0 - sXm

!         iYm = floor(rj-0.5d0)
         iYm = floor(rj+0.5d0)
         iYp = iYm + 1
!         sYm = (rj-0.5d0) - real(iYm,kind=RP)
         sYm = (rj+0.5d0) - real(iYm,kind=RP)
         sYp = 1.d0 - sYm

         zph_s = real(REAL_FZ(KS-1,iXm,iYm),kind=RP) * ( SXp * SYp )       &
               + real(REAL_FZ(KS-1,iXp,iYm),kind=RP) * ( SXm * SYp )       &
               + real(REAL_FZ(KS-1,iXm,iYp),kind=RP) * ( SXp * SYm )       &
               + real(REAL_FZ(KS-1,iXp,iYp),kind=RP) * ( SXm * SYm )

         if( sd_z(n)<real(zph_s+sdm_zlower,kind=RP) ) then
            sd_z(n) = INVALID
            cycle
         end if

         !### get vertical index of super-droplets ###!

!         do k=2,nk-1
         do k=KS-1,KE

            zph_u = real(REAL_FZ(k+1,iXm,iYm),kind=RP) * ( SXp * SYp )  &
                  + real(REAL_FZ(k+1,iXp,iYm),kind=RP) * ( SXm * SYp )  &
                  + real(REAL_FZ(k+1,iXm,iYp),kind=RP) * ( SXp * SYm )  &
                  + real(REAL_FZ(k+1,iXp,iYp),kind=RP) * ( SXm * SYm )

            zph_l = real(REAL_FZ(k,iXm,iYm),kind=RP) * ( SXp * SYp )    &
                  + real(REAL_FZ(k,iXp,iYm),kind=RP) * ( SXm * SYp )    &
                  + real(REAL_FZ(k,iXm,iYp),kind=RP) * ( SXp * SYm )    &
                  + real(REAL_FZ(k,iXp,iYp),kind=RP) * ( SXm * SYm )

            if( sd_z(n)>=zph_l .and. sd_z(n)<zph_u ) then

               sd_rk(n) = real(k,kind=RP)                               &
                        + (sd_z(n)-zph_l)/(zph_u-zph_l)

            end if
         end do

      end do
      !$omp end parallel do simd

    return
  end subroutine sdm_z2rk
  !-----------------------------------------------------------------------------
  subroutine sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    use scale_atmos_grid_cartesC, only: &
         GRID_FX => ATMOS_GRID_CARTESC_FX, &
         DX
    use scale_precision
    use scale_io
    use scale_prc, only: &
         PRC_abort
    use m_sdm_common, only: &
         VALID2INVALID

    integer, intent(in) :: IS,IE
    integer, intent(in) :: sd_num       ! number of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num)   ! x-coordinate of super-droplets
    real(RP), intent(out) :: sd_ri(1:sd_num) ! index-i(real) of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! index-k(real) of super-droplets
    
    ! Work variables
    integer :: i,n    ! index
    integer :: error

    ! only for uniform grid. temporarily used
    real(RP) :: dx_inv
    dx_inv=1.0_RP/DX

    ! Get horizontal index[i/real] of super-droplets
    error=0
    !$omp  parallel do simd schedule(static,256) reduction(min:error)
    do n=1,sd_num
       
       if( sd_rk(n)<VALID2INVALID ) cycle
                
!!$       if((sd_x(n)<GRID_FX(IS-1)).or.(sd_x(n)>GRID_FX(IE))) then
!!$          error=n
!!$       end if
       
       !### get horizontal face index(real) of super-droplets ###!
! not used temporarily
!!$       iloop: do i = IS, IE
!!$          if( sd_x(n) < GRID_FX(i) ) then
!!$             sd_ri(n) = real(i-1,kind=RP) + ( sd_x(n)-GRID_FX(i-1) ) / ( GRID_FX(i)-GRID_FX(i-1) )
!!$             exit iloop
!!$          endif
!!$       enddo iloop

       ! only for uniform grid. temporarily used
       sd_ri(n) = real((IS-1),kind=RP) + (sd_x(n)-GRID_FX(IS-1))*dx_inv

    end do
    !$omp end parallel do simd

    if(error>0)then
       LOG_ERROR("sdm_x2ri",*) 'sd_x out of region', n, sd_x(n), GRID_FX(IS-1), GRID_FX(IE)
       call PRC_abort
    end if

    return
  end subroutine sdm_x2ri
  !-----------------------------------------------------------------------------
  subroutine sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)
    use scale_atmos_grid_cartesC, only: &
         GRID_FY => ATMOS_GRID_CARTESC_FY, &
         DY
    use scale_precision
    use scale_io
    use scale_prc, only: &
         PRC_abort
    use m_sdm_common, only: &
         VALID2INVALID

    integer, intent(in) :: JS,JE
    integer, intent(in) :: sd_num       ! number of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)   ! y-coordinate of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num) ! index-j(real) of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num) ! index-k(real) of super-droplets
    
    ! Work variables
    integer :: j,n    ! index
    integer :: error

    ! only for uniform grid. temporarily used
    real(RP) :: dy_inv
    dy_inv=1.0_RP/DY

    ! Get horizontal index[i/real] of super-droplets
    error=0
    !$omp  parallel do simd schedule(static,256) reduction(min:error)
    do n=1,sd_num
       
       if( sd_rk(n)<VALID2INVALID ) cycle

!!$       if((sd_y(n)<GRID_FY(JS-1)).or.(sd_y(n)>GRID_FY(JE))) then
!!$          error=n
!!$       end if
       
       !### get horizontal face index(real) of super-droplets ###!
! not used temporarily
!!$       jloop: do j = JS, JE
!!$          if( sd_y(n) < GRID_FY(j) ) then
!!$             sd_rj(n) = real(j-1,kind=RP) + ( sd_y(n)-GRID_FY(j-1) ) / ( GRID_FY(j)-GRID_FY(j-1) )
!!$             exit jloop
!!$          endif
!!$       enddo jloop
!!$    end do

       ! only for uniform grid. temporarily used
       sd_rj(n) = real((JS-1),kind=RP) + (sd_y(n)-GRID_FY(JS-1))*dy_inv

    end do
    !$omp end parallel do simd

    if(error>0)then
       LOG_ERROR("sdm_y2rj",*) 'sd_y out of region', n, sd_y(n), GRID_FY(JS-1), GRID_FY(JE)
       call PRC_abort
    end if

    return
  end subroutine sdm_y2rj
end module m_sdm_coordtrans
