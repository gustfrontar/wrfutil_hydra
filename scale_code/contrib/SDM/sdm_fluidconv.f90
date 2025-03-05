!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Coversion between atmospheric fluid variables
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
!! @li      2014-07-14 (S.Shima) [rev] sdm_rhot_qtrc2cpexnr added
!! @li      2014-07-24 (Y.Sato)  [mod] Modify bugs accessing upper/lower boundary
!! @li      2014-07-25 (Y.Sato)  [mod] Modify a bug in sdm_rho_mom2uvw
!! @li      2018-06-05 (S.Shima) [fix] Exner function in sdm_rhot_qtrc2cpexnr
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_fluidconv

  implicit none
  private
  public :: sdm_rhot_qtrc2p_t, sdm_rho_rhot2pt, sdm_rho_mom2uvw, sdm_rho_qtrc2rhod,sdm_rhot_qtrc2cpexnr

contains
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_qtrc2rhod(    &
       KA, IA, JA, &
       dens,qtrc,rhod)
    use scale_precision
    use m_sdm_common, only: &
         I_QV_sdm,QA_MP_sdm
    ! Input variables
    integer, intent(in) :: KA, IA, JA
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: qtrc(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    ! Output variables
    real(RP), intent(out) :: rhod(KA,IA,JA)  ! Density of dry air [kg/m3]

    real(RP) :: qdry
    integer :: i, j, k,iq  ! index
    !---------------------------------------------------------------------

    !$omp parallel
    !$omp do private(j,i,k,qdry) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       ! rho = rho_d + rho_v  when the SDM is used
       qdry        = 1.0_RP - qtrc(k,i,j,I_QV_sdm)
       rhod(k,i,j) = dens(k,i,j)*qdry
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_rho_qtrc2rhod
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_mom2uvw(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dens,momx,momy,momz,u,v,w)
    use scale_precision
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in)  :: dens(KA,IA,JA) ! Density [kg/m3]
    real(RP), intent(in)  :: momx(KA,IA,JA) ! Momentum [kg/s/m2]
    real(RP), intent(in)  :: momy(KA,IA,JA) 
    real(RP), intent(in)  :: momz(KA,IA,JA) 
    ! Output variables
    real(RP), intent(out) :: u(KA,IA,JA)  ! wind velocity in face grid [m/s]
    real(RP), intent(out) :: v(KA,IA,JA)  ! 
    real(RP), intent(out) :: w(KA,IA,JA)  ! 

    integer :: i, j, k  ! index
    !---------------------------------------------------------------------

#ifdef QUICKDEBUG
    u(:,:,:) = -999.0_RP
    v(:,:,:) = -999.0_RP
    w(:,:,:) = -999.0_RP
#endif

    !$omp parallel

    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = max(IS-1,1), IE
    do k = KS, KE          
       u(k,i,j) = 2.0_RP * momx(k,i,j) / ( dens(k,i,j)+dens(k,i+1,j) )
    enddo
    enddo
    enddo
    !$omp end do

    !! fill upper and lower HALO
    !$omp do private(j,i)
    do j = JS, JE
    do i = max(IS-1,1), IE
       u(1:KS-1,i,j)  = u(KS,i,j)       
       u(KE+1:KA,i,j) = u(KE,i,j)
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(2)
    do j = max(JS-1,1), JE
    do i = IS, IE
    do k = KS, KE          
       v(k,i,j) = 2.0_RP * momy(k,i,j) / ( dens(k,i,j)+dens(k,i,j+1) )
    enddo
    enddo
    enddo
    !$omp end do

    !! fill upper and lower HALO
    !$omp do private(j,i)
    do j = max(JS-1,1), JE
    do i = IS, IE
       v(1:KS-1,i,j)  = v(KS,i,j)
       v(KE+1:KA,i,j) = v(KE,i,j)       
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS-1, KE
       w(k,i,j) = 2.0_RP * momz(k,i,j) / ( dens(k,i,j)+dens(k+1,i,j) )
    enddo
    enddo
    enddo
    !$omp end do

    !! fill upper and lower HALO with 0
    !! also set w at k=KS-1 and KE as 0
    !$omp do private(j,i)
    do j = JS, JE
    do i = IS, IE
       w(1:KS-1,i,j) = 0.0_RP
       w(KE:KA,i,j)  = 0.0_RP
    enddo
    enddo
    !$omp end do

    !$omp end parallel
    
    return
  end subroutine sdm_rho_mom2uvw
  !-----------------------------------------------------------------------------
  subroutine sdm_rho_rhot2pt(    &
       KA, IA, JA, &
       dens,rhot,pt)
    use scale_precision
    ! Input variables
    integer, intent(in) :: KA, IA, JA
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: rhot(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    ! Output variables
    real(RP), intent(out) :: pt(KA,IA,JA)         !  Potential temperature [K]
    
    integer :: i, j, k  ! index
    !---------------------------------------------------------------------

    !$omp parallel do private(j,i,k) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA          
       pt(k,i,j) = RHOT(k,i,j) / DENS(k,i,j)
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine sdm_rho_rhot2pt
  !-----------------------------------------------------------------------------
  subroutine sdm_rhot_qtrc2p_t(    &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       rhot,qtrc,dens,p,t)
    use scale_precision
    use scale_const, only: &
         P00 => CONST_PRE00, &
         Rdry => CONST_Rdry, &
         Rvap => CONST_Rvap, &
         CPdry => CONST_CPdry, &
         CPw => CONST_CPvap
    use m_sdm_common, only: &
         I_QV_sdm,QA_MP_sdm
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: rhot(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: qtrc(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    ! Output variables
    real(RP), intent(out) :: p(KA,IA,JA)          !  Pressure [Pa]
    real(RP), intent(out) :: t(KA,IA,JA)          !  Temperature [K]

    real(RP) :: qdry
    real(RP) :: rtot
    real(RP) :: cptot
    real(RP) :: cpovcv
    integer :: i, j, k, iq  ! index
    !---------------------------------------------------------------------

    !$omp parallel
    !$omp do private(j,i,k,qdry,rtot,cptot,cpovcv) collapse(3)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
!!$      do iq = QQS, QQE
!!$        qdry(k,i,j) = qdry(k,i,j) - qtrc(k,i,j,iq)
!!$      enddo
      ! rho = rho_d + rho_v  when the SDM is used
      qdry  = 1.0_RP - qtrc(k,i,j,I_QV_sdm)
      rtot  = Rdry * qdry + Rvap * qtrc(k,i,j,I_QV_sdm)
      cptot = CPdry * qdry
!!$      do iq = QQS, QQE
!!$        cptot = cptot + qtrc(k,i,j,iq) * CPw(iq)
!!$      enddo
      cptot  = cptot + qtrc(k,i,j,I_QV_sdm) * CPw
      cpovcv = cptot / ( cptot - rtot )

      p(k,i,j) = P00 * ( rhot(k,i,j) * rtot / P00 )**cpovcv
      t(k,i,j) = (rhot(k,i,j)/dens(k,i,j)) * (p(k,i,j)/P00)**(rtot/cptot)
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_rhot_qtrc2p_t
  !---------------------------------------------------------------------
  subroutine sdm_rhot_qtrc2cpexnr(    &
       KA, IA, JA, &
       rhot,qtrc,dens,cpexnr)
    use scale_precision
    use scale_const, only: &
         P00 => CONST_PRE00, &
         Rdry => CONST_Rdry, &
         Rvap => CONST_Rvap, &
         CPdry => CONST_CPdry, &
         CPw => CONST_CPvap
    use m_sdm_common, only: &
         I_QV_sdm,QA_MP_sdm
    ! Input variables
    integer, intent(in) :: KA, IA, JA
    real(RP), intent(in) :: rhot(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: qtrc(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(in) :: dens(KA,IA,JA)        !! Density [kg/m3]
    ! Output variables
    real(RP), intent(out) :: cpexnr(KA,IA,JA)     !  cp*exner_function

    real(RP) :: qdry
    real(RP) :: rtot
    real(RP) :: cptot
    integer :: i, j, k, iq  ! index
    !---------------------------------------------------------------------

    !$omp parallel
    !$omp do private(j,i,k,qdry,rtot,cptot) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA          
       qdry         = 1.0_RP - qtrc(k,i,j,I_QV_sdm)
       rtot         = Rdry * qdry + Rvap * qtrc(k,i,j,I_QV_sdm)
       cptot = CPdry * qdry
       cptot = cptot + qtrc(k,i,j,I_QV_sdm) * CPw
       cpexnr(k,i,j) = cptot*(rhot(k,i,j)*rtot/P00)**(rtot/(cptot-rtot))
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_rhot_qtrc2cpexnr
end module m_sdm_fluidconv
