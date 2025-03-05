!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Covert super-droplets to fluid variables
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
!! @li      2014-07-14 (S.Shima) [new] Separated from scale_atmos_phy_mp_sdm.F90
!! @li      2014-07-14 (S.Shima) [rev] sdm_sd2prec added
!! @li      2014-07-22 (Y.Sato ) [mod] Modify the definition of kl and ku for calculating drate
!! @li      2014-07-24 (Y.Sato ) [mod] Modify bugs accessing upper/lower boundary
!! @li      2015-06-27 (S.Shima) [add] Add fapp_start/stop calls to monitor the performance
!! @li      2015-06-27 (S.Shima) [mod] Working arrays are introduced to parallelize sd2rhow and sd2rhocr 
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2016-07-16 (S.Shima) [add] sdm_sd2qiqsqg added
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2qcqr, sdm_sd2rhocr, sdm_sd2rhow modified
!! @li      2016-07-16 (S.Shima) [add] sdm_sd2rhosol added
!! @li      2016-07-16 (S.Shima) [mod] sdm_sd2prec modified
!! @li      2016-07-22 (S.Shima) [fix] sdm_sd2rhoisg is fixed not to overcount graupel and snow 
!! @li      2018-07-01 (S.Shima) [mod] categorization of ice, snow, graupel, based on mrime and nmono
!! @li      2019-01-09 (S.Shima) [mod] difinition of ice, snow, graupel
!! @li      2019-07-05 (S.Shima) [fix] sdm_sd2prec to reset rain/snow rate everytime
!! @li      2019-07-05 (S.Shima) [fix] for reverting sdm_sd2prec
!! @li      2019-07-28 (S.Shima) [add] sdm_sd2rhodropmom
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_sd2fluid
  use scale_precision

  implicit none
  private
  public :: sdm_sd2prec, sdm_sd2rhow, sdm_sd2rhocr, sdm_sd2qcqr, sdm_sd2qiqsqg, sdm_sd2rhosol, sdm_sd2rhosd, sdm_sd2rhodropmom

contains
  subroutine sdm_sd2prec(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       dtb_crs,                        &
       prec,sd_num,sd_n,sd_liqice,sd_x,sd_y,     &
       sd_r,sdi,sd_ri,sd_rj,sd_rk,ilist_s,ilist_r,pr_sdm)
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, PREC2INVALID,INVALID, &
         sdm_cold, i2, sdicedef, STAT_LIQ, STAT_ICE
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: dtb_crs
    integer, intent(in) :: sd_num  ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num) ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    ! Input and output variables
    real(RP), intent(out) :: sd_ri(1:sd_num)   ! index-i(real) of super-droplets
    real(RP), intent(out) :: sd_rj(1:sd_num)   ! index-j(real) of super-droplets
    real(RP), intent(inout) :: sd_rk(1:sd_num) ! index[k/real] of super-droplets
    real(RP), intent(inout) :: prec(IA,JA,1:6) ! precipitation rate [m/s] and accumlation [m]
                                               ! 1:rain rate, 2:rain accumulation
                                               ! 3:snow rate, 4:snow accumulation
                                               ! 5:precipitation rate, 6:precipitation accumulation
    ! Output variables
    real(RP), intent(out) :: pr_sdm(1:IA,1:JA) ! temporary buffer of CReSS dimension
    integer, intent(out) :: ilist_r(1:sd_num) ! buffer for list vectorization
    integer, intent(out) :: ilist_s(1:sd_num) ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(IA,JA) ! coef.
    real(RP) :: dtmp,dtmp2   ! temporary variables
    real(RP) :: dtbiv        ! 1.e0 / time step

    real(RP) :: rk_tmp(1:sd_num)

    integer :: tlist_r, tlist_s              ! total list number
    integer :: cnt_r, cnt_s                ! counter

    integer :: i, j, m, n ! index
    !-------------------------------------------------------------------

    ! Initialize
    dtbiv = 1.0_RP / dtb_crs
    !$omp parallel
    !$omp workshare
    dcoef(1:IA,1:JA)=0.0_RP
    !$omp end workshare

    !$omp do private(j,i)
    do j = JS, JE
    do i = IS, IE
       dcoef(i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j)
    enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! Get index list for compressing buffer.
    cnt_r=0
    cnt_s=0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID .and. &
            sd_rk(n)>PREC2INVALID ) then
          if( sd_liqice(n) .eq. STAT_LIQ) then
             cnt_r = cnt_r + 1
             ilist_r(cnt_r) = n
          else if( sd_liqice(n) .eq. STAT_ICE) then
             cnt_s = cnt_s + 1
             ilist_s(cnt_s) = n
          end if
          rk_tmp(n) = 0.0 ! VALID
       else
          rk_tmp(n) = INVALID
       end if
    end do

    !### get horizontal face index(real) of super-droplets ###!
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,rk_tmp)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,rk_tmp)

    tlist_r = cnt_r
    tlist_s = cnt_s

    ! Get rain rate and accumulation
    if( tlist_r>0 ) then

       !$omp parallel
       !$omp do private(j,i) collapse(2)
       do j=1,JA
       do i=1,IA
          pr_sdm(i,j) = 0.0_RP
       end do
       end do
       !$omp end do

       !$omp do private(n,i,j,dtmp2) reduction(+:pr_sdm)
       do m=1,tlist_r
          n=ilist_r(m)
          i=floor(sd_ri(n))+1
          j=floor(sd_rj(n))+1

          dtmp2 = sd_r(n) * sd_r(n) * sd_r(n)           &
                  * real(sd_n(n),kind=RP)
          pr_sdm(i,j) = pr_sdm(i,j) + dtmp2

          sd_rk(n) = INVALID     !! convert to invalid
       end do
       !$omp end do
       
       !### convert super-droplets to precipitation ###!
       !$omp do private(j,i,dtmp) collapse(2)
       do j=1,JA
       do i=1,IA

          dtmp = pr_sdm(i,j) * dcoef(i,j)

          !! rain fall rate
          prec(i,j,1) = dtmp * dtbiv

          !! rain accumulation
          prec(i,j,2) = prec(i,j,2) + dtmp
          
       end do
       end do
       !$omp end do
       !$omp end parallel
    else
       !$omp parallel workshare
       prec(:,:,1) = 0.0_RP
       !$omp end parallel workshare
    end if

    ! Get snow rate and accumulation
    if( sdm_cold .and. (tlist_s>0) ) then

       !$omp parallel
       !$omp do private(j,i) collapse(2)
       do j=1,JA
       do i=1,IA
          pr_sdm(i,j) = 0.0_RP
       end do
       end do
       !$omp end do

       !$omp single
       do m=1,tlist_s
          n=ilist_s(m)
          i=floor(sd_ri(n))+1
          j=floor(sd_rj(n))+1

          dtmp2 = sdi%re(n) * sdi%re(n) * sdi%rp(n) * sdi%rho(n) / rw &
               * real(sd_n(n),kind=RP)

          pr_sdm(i,j) = pr_sdm(i,j) + dtmp2

          sd_rk(n) = INVALID     !! convert to invalid
       end do
       !$omp end single

       !### convert super-droplets to precipitation ###!
       !$omp do private(j,i,dtmp) collapse(2)
       do j=1,JA
       do i=1,IA

          dtmp = real( pr_sdm(i,j) * dcoef(i,j) )

          !! snow fall rate
          prec(i,j,3) = dtmp * dtbiv

          !! snow accumulation
          prec(i,j,4) = prec(i,j,4) + dtmp
          
       end do
       end do
       !$omp end do
       !$omp end parallel
    else
       !$omp parallel workshare
       prec(:,:,3) = 0.0_RP
       !$omp end parallel workshare
    end if

    ! Get precipitation rate and accumulation
    !$omp parallel
    !$omp do private(j,i) collapse(2)
    do j=1,JA
    do i=1,IA

       !! precipitation rate
       prec(i,j,5) = prec(i,j,1) + prec(i,j,3) 

       !! precipitation accumulation
       prec(i,j,6) = prec(i,j,2) + prec(i,j,4) 
          
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_sd2prec
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2qcqr(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,QC,QR,        &
       zph_crs,           &
       sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                           &
       rhoc_sdm,rhor_sdm,              &
       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2    )
    use m_sdm_common, only: &
         sdm_rqc2qr,i2
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(out):: QC(KA,IA,JA)   ! rhoc/rho
    real(RP), intent(out):: QR(KA,IA,JA)   ! rhor/rho
    real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)      ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)  ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)  ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rhoc_sdm(KA,IA,JA)   ! density of cloud water
    real(RP), intent(out) :: rhor_sdm(KA,IA,JA)   ! density of rain water
    real(RP), intent(out) :: crs_dtmp1(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp2(KA,IA,JA)    ! temporary buffer of CReSS dimension
    integer, intent(out) :: sd_itmp1(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp2(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    ! Work variables
    integer :: i, j, k   ! index

    !-------------------------------------------------------------------

    ! Get density of water hydrometeor

    !### the case updating water hydrometeor by SDM ###!
    !! convert super-droplets to density of cloud water
    !! and rain water.
    
    call sdm_sd2rhocr(  &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         sdm_rqc2qr,                          &
         zph_crs,rhoc_sdm,rhor_sdm,           &
         sd_num,sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,    &
         sd_rkl,sd_rku,                       &
         crs_dtmp1,crs_dtmp2,sd_itmp1,sd_itmp2)

    ! Convert water hydrometeor density to mixing ratio
    !$omp parallel do private(j,i,k) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       QC(k,i,j) = rhoc_sdm(k,i,j)/DENS(k,i,j)
       QR(k,i,j) = rhor_sdm(k,i,j)/DENS(k,i,j)
    end do
    end do
    end do
    !$omp end parallel do

    return
  end subroutine sdm_sd2qcqr
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2qiqsqg(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,QI,QS,QG,        &
       zph_crs,           &
       sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                           &
       rhoi_sdm,rhos_sdm,rhog_sdm,              &
       sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3    )
    use m_sdm_common, only: &
         sdm_rqi2qs, i2, sdicedef
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: DENS(KA,IA,JA) ! Density     [kg/m3]
    real(RP), intent(out):: QI(KA,IA,JA)   ! rhoi/rho
    real(RP), intent(out):: QS(KA,IA,JA)   ! rhos/rho
    real(RP), intent(out):: QG(KA,IA,JA)   ! rhog/rho
    real(RP), intent(in) :: zph_crs(KA,IA,JA)  ! z physical coordinate
    integer, intent(in) :: sd_num    ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)      ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)      ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)     ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)     ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)     ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)  ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)  ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(inout) :: rhoi_sdm(KA,IA,JA)   ! density of cloud ice
    real(RP), intent(inout) :: rhos_sdm(KA,IA,JA)   ! density of snow
    real(RP), intent(inout) :: rhog_sdm(KA,IA,JA)   ! density of graupel
    real(RP), intent(out) :: crs_dtmp1(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp2(KA,IA,JA)    ! temporary buffer of CReSS dimension
    real(RP), intent(out) :: crs_dtmp3(KA,IA,JA)    ! temporary buffer of CReSS dimension
    integer, intent(out) :: sd_itmp1(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp2(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    integer, intent(out) :: sd_itmp3(1:sd_num)    ! temporary array of the size of the number of super-droplets.
    ! Work variables
    integer :: i, j, k   ! index

    !-------------------------------------------------------------------

    ! Get density of solid hydrometeor
    call sdm_sd2rhoisg(  &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         sdm_rqi2qs,                        &
         zph_crs,rhoi_sdm,rhos_sdm,rhog_sdm,              &
         sd_num,sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
         sd_rkl,sd_rku,                                   &
         crs_dtmp1,crs_dtmp2,crs_dtmp3,sd_itmp1,sd_itmp2,sd_itmp3)

    ! Convert water hydrometeor density to mixing ratio
    !$omp parallel do private(j,i,k) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       QI(k,i,j) = rhoi_sdm(k,i,j)/DENS(k,i,j)
       QS(k,i,j) = rhos_sdm(k,i,j)/DENS(k,i,j)
       QG(k,i,j) = rhog_sdm(k,i,j)/DENS(k,i,j)
    end do
    end do
    end do
    !$omp end parallel do

    return
  end subroutine sdm_sd2qiqsqg
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhocr(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       sdm_rqc2qr,                        &
       zph_crs,rhoc_sdm,rhor_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       liqc_sdm,liqr_sdm,ilist_c,ilist_r)
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, i2, STAT_LIQ
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP),intent(in) :: sdm_rqc2qr          ! Threshould between qc and qr [m]
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num)     ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhoc_sdm(KA,IA,JA)  ! densitiy of cloud water
    real(RP), intent(out) :: rhor_sdm(KA,IA,JA)  ! densitiy of rain water
    real(RP), intent(out) :: liqc_sdm(KA,IA,JA)  ! liquid cloud water of super-droplets
    real(RP), intent(out) :: liqr_sdm(KA,IA,JA)  ! liquid rain water of super-droplets
    integer, intent(out) :: ilist_c(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_r(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate,dtmp        ! temporary
    integer :: tlist_c            ! total list number for cloud
    integer :: tlist_r            ! total list number for rain
    integer :: ccnt               ! counter
    integer :: rcnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhocr",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)
    
    ! Initialize
    !$omp parallel
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA
    do k=1,KA
       liqc_sdm(k,i,j) = 0.0_RP
       liqr_sdm(k,i,j) = 0.0_RP
       rhoc_sdm(k,i,j) = 0.0_RP
       rhor_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    tlist_c = 0
    tlist_r = 0
    ccnt = 0
    rcnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_LIQ ) cycle

       if( sd_r(n)<real(sdm_rqc2qr,kind=RP) ) then

          !### cloud-water ###!

          ccnt = ccnt + 1
          ilist_c(ccnt) = n

       else
          
          !### rain-water ###!
          rcnt = rcnt + 1
          ilist_r(rcnt) = n

       end if

    end do

    tlist_c = ccnt
    tlist_r = rcnt
    !$omp end single

    ! Get density of cloud-water and rain-water.

    !### cloud-water ###!

    if( tlist_c>0 ) then

       !$omp single
       do m=1,tlist_c
          n = ilist_c(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = rw * sd_r(n) * sd_r(n) * sd_r(n) * real(sd_n(n),kind=RP)

          liqc_sdm(k,i,j) = liqc_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust cloud-water in verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             liqc_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             liqc_sdm(kl,i,j) = liqc_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             liqc_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             liqc_sdm(ku,i,j) = liqc_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do

       !=== convert super-droplets to density of cloud-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhoc_sdm(k,i,j) = liqc_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if

    !### rain-water ###!

    if( tlist_r>0 ) then

       !$omp single
       do m=1,tlist_r
          n = ilist_r(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = rw * sd_r(n) * sd_r(n) * sd_r(n) * real(sd_n(n),kind=RP)

          liqr_sdm(k,i,j) = liqr_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust rain-water in verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)

          if( drate<0.5_RP ) then
             liqr_sdm(kl,i,j) = 0.e0           !! <50% in share
          else
             liqr_sdm(kl,i,j) = liqr_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.5_RP ) then
             liqr_sdm(ku,i,j) = 0.e0           !! <50% in share
          else
             liqr_sdm(ku,i,j) = liqr_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do

       !=== convert super-droplets to density of rain-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhor_sdm(k,i,j) = liqr_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if
    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhocr",1,1)
#endif
    return
  end subroutine sdm_sd2rhocr
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhow(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       zph_crs,rhow_sdm,sd_num,sd_n,sd_liqice,  &
       sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,sd_rkl,sd_rku,      &
       liqw_sdm,ilist)
    use scale_const, only: &
         rw => CONST_DWATR
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, i2, STAT_LIQ
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP),intent(in) :: zph_crs(KA,IA,JA)    ! z physical coordinate
    integer, intent(in) :: sd_num               ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA) ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA) ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rhow_sdm(KA,IA,JA) ! densitiy of liquid water
    real(RP), intent(out) :: liqw_sdm(KA,IA,JA) ! liquid water of super-droplets
    integer, intent(out) :: ilist(1:sd_num)   ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)    ! coef.
    real(RP) :: drate        ! temporary
    real(RP) :: dtmp        ! temporary
    integer :: cnt                ! counter
    integer :: i, j, k, kl, ku, m, n    ! index
    integer :: tlist
    !--------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhow",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    !$omp parallel
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA
    do k=1,KA
       liqw_sdm(k,i,j) = 0.0_RP
       rhow_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    cnt =0
    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_LIQ ) cycle

       cnt = cnt + 1
       ilist(cnt) = n
    end do

    tlist = cnt

    !$omp end single

    ! Get density of liquid-water.
    !### count voulme of super-droplets ###!
    !$omp do private(n,i,j,k,dtmp) reduction(+:liqw_sdm)
    do m=1,tlist
       n = ilist(m)

       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       dtmp = rw * sd_r(n) * sd_r(n) * sd_r(n) * real(sd_n(n),kind=RP)

       liqw_sdm(k,i,j) = liqw_sdm(k,i,j) + dtmp

    end do
!    !$omp end single

    ! Adjust liquid-water in verical boundary.
    !$omp do private(j,i,kl,drate,ku)
    do j=JS,JE
    do i=IS,IE
       !! at lower boundary
       kl    = floor(sd_rkl(i,j))+1
       drate = real(kl,kind=RP) - sd_rkl(i,j)
       if( drate<0.50_RP ) then
          liqw_sdm(kl,i,j) = 0.0_RP           !! <50% in share
       else
          liqw_sdm(kl,i,j) = liqw_sdm(kl,i,j)/drate
       end if

       !! at upper boundary
       ku    = floor(sd_rku(i,j))+1
       drate = sd_rku(i,j) - real(ku-1,kind=RP)
       if( drate<0.50_RP ) then
          liqw_sdm(ku,i,j) = 0.0_RP           !! <50% in share
       else
          liqw_sdm(ku,i,j) = liqw_sdm(ku,i,j)/drate
       end if
    end do
    end do
    !$omp end do

    ! Convert super-droplets to density of liquid-water.
    !$omp do private(j,i,k) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       rhow_sdm(k,i,j) = liqw_sdm(k,i,j) * dcoef(k,i,j)
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhow",1,1)
#endif
    return
  end subroutine sdm_sd2rhow
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhoisg(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       sdm_rqi2qs,                        &
       zph_crs,rhoi_sdm,rhos_sdm,rhog_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       soli_sdm,sols_sdm,solg_sdm,ilist_i,ilist_s,ilist_g)
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, &
         i2, STAT_ICE, sdicedef 
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP),intent(in) :: sdm_rqi2qs          ! Threshould between qi and qs [m]
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhoi_sdm(KA,IA,JA)  ! density of cloud ice
    real(RP), intent(out) :: rhos_sdm(KA,IA,JA)  ! density of snow
    real(RP), intent(out) :: rhog_sdm(KA,IA,JA)  ! density of graupel
    real(RP), intent(out) :: soli_sdm(KA,IA,JA)  ! cloud ice of super-droplets
    real(RP), intent(out) :: sols_sdm(KA,IA,JA)  ! snow of super-droplets
    real(RP), intent(out) :: solg_sdm(KA,IA,JA)  ! graupel of super-droplets
    integer, intent(out) :: ilist_i(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_s(1:sd_num)  ! buffer for list vectorization
    integer, intent(out) :: ilist_g(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate,dtmp        ! temporary
    integer :: tlist_i            ! total list number for cloud
    integer :: tlist_s            ! total list number for rain
    integer :: tlist_g            ! total list number for rain
    integer :: icnt               ! counter
    integer :: scnt               ! counter
    integer :: gcnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    real(RP) :: sd_vol, sd_mass, sd_mrime
    integer  :: sd_nmono
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhoisg",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)
    
    !$omp parallel

    ! Initialize
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA
    do k=1,KA
       soli_sdm(k,i,j) = 0.0_RP
       sols_sdm(k,i,j) = 0.0_RP
       solg_sdm(k,i,j) = 0.0_RP
       rhoi_sdm(k,i,j) = 0.0_RP
       rhos_sdm(k,i,j) = 0.0_RP
       rhog_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    tlist_i = 0
    tlist_s = 0
    tlist_g = 0
    icnt = 0
    scnt = 0
    gcnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_ICE ) cycle

       sd_vol   = F_THRD * ONE_PI * sdi%re(n)**2 * sdi%rp(n)
       sd_mass  = sd_vol * sdi%rho(n) ! mass of the particle [kg]
       sd_mrime = sdi%mrime(n)
       sd_nmono = sdi%nmono(n)

       if( (sd_mass*0.3_RP) < sd_mrime )then
          !### graupel ###!
          gcnt = gcnt + 1
          ilist_g(gcnt) = n
       else if(sd_nmono>10)then
          !### snow flake ###!
          scnt = scnt + 1
          ilist_s(scnt) = n
       else
          !### ice crystal ###!
          icnt = icnt + 1
          ilist_i(icnt) = n
       end if

    end do

    tlist_i = icnt
    tlist_s = scnt
    tlist_g = gcnt
    !$omp end single

    ! Get density of cloud ice, snow, and graupel.

    !### cloud ice ###!

    if( tlist_i>0 ) then

       !$omp single
       do m=1,tlist_i
          n = ilist_i(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = sdi%rho(n) * sdi%re(n) * sdi%re(n) * sdi%rp(n) * real(sd_n(n),kind=RP)

          soli_sdm(k,i,j) = soli_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust at verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             soli_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             soli_sdm(kl,i,j) = soli_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             soli_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             soli_sdm(ku,i,j) = soli_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do

       !=== convert super-droplets to density of cloud-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhoi_sdm(k,i,j) = soli_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if

    !### snow ###!

    if( tlist_s>0 ) then

       !$omp single
       do m=1,tlist_s
          n = ilist_s(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = sdi%rho(n) * sdi%re(n) * sdi%re(n) * sdi%rp(n) * real(sd_n(n),kind=RP)

          sols_sdm(k,i,j) = sols_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust at verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             sols_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sols_sdm(kl,i,j) = sols_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sols_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sols_sdm(ku,i,j) = sols_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do
       
       !=== convert super-droplets to density of cloud-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhos_sdm(k,i,j) = sols_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if

    !### graupel ###!

    if( tlist_g>0 ) then

       !$omp single
       do m=1,tlist_g
          n = ilist_g(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = sdi%rho(n) * sdi%re(n) * sdi%re(n) * sdi%rp(n) * real(sd_n(n),kind=RP)

          solg_sdm(k,i,j) = solg_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust at verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             solg_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             solg_sdm(kl,i,j) = solg_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             solg_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             solg_sdm(ku,i,j) = solg_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do

       !=== convert super-droplets to density of cloud-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhog_sdm(k,i,j) = solg_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if

    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhoisg",1,1)
#endif
    return
  end subroutine sdm_sd2rhoisg
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhosol(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       zph_crs,rhosol_sdm,         &
       sd_num, sd_n,sd_liqice,sd_x,sd_y,sdi,sd_ri,sd_rj,sd_rk, &
       sd_rkl,sd_rku,                     &
       sol_sdm,ilist)
    use m_sdm_common, only: &
         F_THRD, ONE_PI, dxiv_sdm, dyiv_sdm, VALID2INVALID, &
         i2, STAT_ICE, sdicedef 
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP),intent(in) :: zph_crs(KA,IA,JA)   ! z physical coordinate
    integer, intent(in) :: sd_num              ! Number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)  ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num)     ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num)     ! y-coordinate of super-droplets
    type(sdicedef), intent(in) :: sdi   ! ice phase super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num) ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num) ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)    ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA)      ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA)      ! upper boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(out) :: rhosol_sdm(KA,IA,JA)  ! density of cloud ice
    real(RP), intent(out) :: sol_sdm(KA,IA,JA)  ! cloud ice of super-droplets
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)   ! coef.
    real(RP) :: drate,dtmp        ! temporary
    integer :: tlist              ! total list number for cloud
    integer :: cnt               ! counter
    integer :: i                  ! index
    integer :: j                  ! index
    integer :: k                  ! index
    integer :: kl                 ! index
    integer :: ku                 ! index
    integer :: m                  ! index
    integer :: n                  ! index
    !-----------------------------------------------------------------------------
    
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhosol",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)
    
    !$omp parallel

    ! Initialize
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = F_THRD * ONE_PI * dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA
    do k=1,KA
       sol_sdm(k,i,j) = 0.0_RP
       rhosol_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    tlist = 0
    cnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_ICE ) cycle

       cnt = cnt + 1
       ilist(cnt) = n

    end do

    tlist = cnt
    !$omp end single

    ! Get density of solid phase water

    if( tlist>0 ) then

       !$omp single
       do m=1,tlist
          n = ilist(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = sdi%rho(n) * sdi%re(n) * sdi%re(n) * sdi%rp(n) * real(sd_n(n),kind=RP)

          sol_sdm(k,i,j) = sol_sdm(k,i,j) + dtmp
       end do
       !$omp end single

       !=== adjust at verical boundary. ===!
       !$omp do private(j,i,kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary
          
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             sol_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             sol_sdm(kl,i,j) = sol_sdm(kl,i,j)/drate
          end if

          !! at upper boundary

          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sol_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             sol_sdm(ku,i,j) = sol_sdm(ku,i,j)/drate
          end if

       end do
       end do
       !$omp end do

       !=== convert super-droplets to density of cloud-water. ===!
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rhosol_sdm(k,i,j) = sol_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if
    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhosol",1,1)
#endif
    return
  end subroutine sdm_sd2rhosol
  !-----------------------------------------------------------------------------
  subroutine sdm_sd2rhosd(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       zph_crs,rho_sdm,sd_num,sd_n,  &
       sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,sd_rkl,sd_rku,      &
       ilist)
    use m_sdm_common, only: &
         dxiv_sdm, dyiv_sdm, VALID2INVALID
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP),intent(in) :: zph_crs(KA,IA,JA)    ! z physical coordinate
    integer, intent(in) :: sd_num               ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA) ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA) ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rho_sdm(KA,IA,JA) ! densitiy of liquid water
    integer, intent(out) :: ilist(1:sd_num)   ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)    ! coef.
    real(RP) :: drate        ! temporary
    integer :: cnt                ! counter
    integer :: i, j, k, kl, ku, m, n    ! index
    integer :: tlist
    !--------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhosd",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp parallel

    ! Initialize
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA       
    do k=1,KA
       rho_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    cnt =0
    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle

       cnt = cnt + 1
       ilist(cnt) = n
    end do

    tlist = cnt

    ! Get density of liquid-water.
    !### count voulme of super-droplets ###!
   
    do m=1,tlist
       n = ilist(m)

       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       rho_sdm(k,i,j)=rho_sdm(k,i,j) + 1.0_RP
    end do
    !$omp end single

    ! Adjust liquid-water in verical boundary.
    !$omp do private(j,i,kl,drate,ku)
    do j=JS,JE
    do i=IS,IE
       !! at lower boundary
       kl    = floor(sd_rkl(i,j))+1
       drate = real(kl,kind=RP) - sd_rkl(i,j)
       if( drate<0.50_RP ) then
          rho_sdm(kl,i,j) = 0.0_RP           !! <50% in share
       else
          rho_sdm(kl,i,j) = rho_sdm(kl,i,j)/drate
       end if

       !! at upper boundary
       ku    = floor(sd_rku(i,j))+1
       drate = sd_rku(i,j) - real(ku-1,kind=RP)
       if( drate<0.50_RP ) then
          rho_sdm(ku,i,j) = 0.0_RP           !! <50% in share
       else
          rho_sdm(ku,i,j) = rho_sdm(ku,i,j)/drate
       end if
    end do
    end do
    !$omp end do

    ! Convert super-droplets to density of liquid-water.
    !$omp do private(j,i,k) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       rho_sdm(k,i,j) = rho_sdm(k,i,j) * dcoef(k,i,j)
    end do
    end do
    end do
    !$omp end do

    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhosd",1,1)
#endif
    return
  end subroutine sdm_sd2rhosd
  !---------------------------------------------------------------------------------------------------------------------------------
  subroutine sdm_sd2rhodropmom(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       order_n,zph_crs,rho_sdm,sd_num,sd_n,sd_liqice, &
       sd_x,sd_y,sd_r,sd_ri,sd_rj,sd_rk,sd_rkl,sd_rku,      &
       ilist)
    use m_sdm_common, only: &
         dxiv_sdm, dyiv_sdm, VALID2INVALID, STAT_LIQ, i2
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    integer,intent(in)  :: order_n               ! order of the droplet moment
    real(RP),intent(in) :: zph_crs(KA,IA,JA)    ! z physical coordinate
    integer, intent(in) :: sd_num               ! number of super-droplets
    integer(DP), intent(in) :: sd_n(1:sd_num)   ! multiplicity of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_r(1:sd_num) ! equivalent radius of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_rkl(IA,JA) ! lower boundary index[k/real] at scalar point in SDM calculation area
    real(RP), intent(in) :: sd_rku(IA,JA) ! upper boundary index[k/real] at scalar point in SDM calculation area
    ! Output variables
    real(RP), intent(out) :: rho_sdm(KA,IA,JA) ! densitiy of the n-th order droplet moment
    integer, intent(out) :: ilist(1:sd_num)   ! buffer for list vectorization
    ! Work variables
    real(RP) :: dcoef(KA,IA,JA)    ! coef.
    real(RP) :: drate,dtmp        ! temporary
    integer :: cnt                ! counter
    integer :: i, j, k, kl, ku, m, n    ! index
    integer :: tlist
    !--------------------------------------------------------------------

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_sd2rhodropmom",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    !$omp parallel

    ! Initialize
    !$omp do private(j,i,k) collapse(2)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       dcoef(k,i,j) = dxiv_sdm(i) * dyiv_sdm(j) &
            / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do private(j,i,k) collapse(3)
    do j=1,JA
    do i=1,IA       
    do k=1,KA
       rho_sdm(k,i,j) = 0.0_RP
    end do
    end do
    end do
    !$omp end do

    ! Get index list for compressing buffer.
    !$omp single
    cnt =0
    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n)/=STAT_LIQ ) cycle

       cnt = cnt + 1
       ilist(cnt) = n
    end do

    tlist = cnt
    !$omp end single

    ! Get density of the n-th droplet moment
    !### count the n-th droplet moment in each grid box ###!

    if(tlist>0) then

       !$omp single
       do m=1,tlist
          n = ilist(m)

          i = floor(sd_ri(n))+1
          j = floor(sd_rj(n))+1
          k = floor(sd_rk(n))+1

          dtmp = sd_r(n)**order_n * real(sd_n(n),kind=RP)

          rho_sdm(k,i,j)=rho_sdm(k,i,j) + dtmp
       end do
       !$omp end single
    
       ! Adjust the moment at verical boundary.
       !$omp do private(j,i,kl,drate,ku)
       do j=JS,JE
       do i=IS,IE
          !! at lower boundary
          kl    = floor(sd_rkl(i,j))+1
          drate = real(kl,kind=RP) - sd_rkl(i,j)
          if( drate<0.50_RP ) then
             rho_sdm(kl,i,j) = 0.0_RP           !! <50% in share
          else
             rho_sdm(kl,i,j) = rho_sdm(kl,i,j)/drate
          end if

          !! at upper boundary
          ku    = floor(sd_rku(i,j))+1
          drate = sd_rku(i,j) - real(ku-1,kind=RP)
          if( drate<0.50_RP ) then
             rho_sdm(ku,i,j) = 0.0_RP           !! <50% in share
          else
             rho_sdm(ku,i,j) = rho_sdm(ku,i,j)/drate
          end if
       end do
       end do
       !$omp end do

       ! Convert super-droplets to density of the n-th moment
       !$omp do private(j,i,k) collapse(2)
       do j=JS,JE
       do i=IS,IE
       do k=KS,KE
          rho_sdm(k,i,j) = rho_sdm(k,i,j) * dcoef(k,i,j)
       end do
       end do
       end do
       !$omp end do

    end if

    !$omp end parallel

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_sd2rhodropmom",1,1)
#endif
    return
  end subroutine sdm_sd2rhodropmom
  !---------------------------------------------------------------------------------------------------------------------------------
end module m_sdm_sd2fluid
