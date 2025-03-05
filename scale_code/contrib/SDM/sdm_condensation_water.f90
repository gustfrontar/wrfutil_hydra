!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Condensation of water to the particles
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
!! @li      2014-07-18 (Y.Sato)  [rev] Change Teten's formula to SATURATION library of SCALE Library 
!! @li      2014-07-24 (Y.Sato)  [mod] Modify a bug relating to the revision on 2014/07/18
!! @li      2014-12-12 (Y.Sato)  [mod] modify epsva, LatHet, and DNS_RL as those used in SCALE Library
!! @li      2014-12-19 (Y.Sato)  [mod] modify the location for defining LatHet
!! @li      2014-12-24 (Y.Sato)  [mod] Modify the Latent Heat for SCALE library 
!! @li      2015-06-25 (S.Shima) [add] fapp_start/stop added for performance monitoring
!! @li      2015-06-26 (S.Shima) [add] OCL added for Auto parallelization on K/FX10
!! @li      2015-07-30 (Y.Sato)  [add] Add "ifdef" for fapp and fipp module
!! @li      2016-07-18 (S.Shima) [mod] Only for liquid water phase particles 
!! @li      2019-10-05 (S.Shima) [mod] limiter of super saturation for the 1st 1h is introduced for DYCOMSII(RF02) (Ackerman et al., 2009)
!! @li      2019-10-07 (S.Shima) [add] aslset=5 for DYCOMSII(RF02) (Ackerman et al. 2009)
!!
!<
!-------------------------------------------------------------------------------
module m_sdm_condensation_water
  use scale_precision  
  implicit none
  private
  public :: sdm_condevp,sdm_condevp_updatefluid

contains
  subroutine sdm_condevp(    &
       KA,KS,KE, IA,IS,IE, JA,JS,JE, &
       sdm_aslset,                     &
       sdm_aslmw,sdm_aslion,sdm_dtevl, &
       pres_scale,t_scale,qv_scale, &
       sd_num,sd_numasl,sd_liqice,sd_x,sd_y,     &
       sd_r,sd_asl,sd_ri,sd_rj,sd_rk              )
    
    use scale_const, only: &
         cp   => CONST_CPdry, &
         p0   => CONST_PRE00, &
         t0   => CONST_TEM00, &
         es0  => CONST_PSAT0, &
         epsva=> CONST_EPSvap, &
         rd   => CONST_Rdry
    use scale_atmos_saturation, only: &
        ATMOS_SATURATION_pres2qsat_liq, &
        ATMOS_SATURATION_psat_liq
    use scale_time, only:     &
         TIME_NOWSEC
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         sdm_wbc,sdm_ebc,sdm_sbc,sdm_nbc, &
         mass_amsul, ion_amsul, mass_nacl, ion_nacl, &
         VALID2INVALID, &
         ! RLRv_D, LatGas, L_RL_K, CurveF, ASL_FF, STAT_LIQ, i2
         LatGas, CurveF, ASL_FF, STAT_LIQ, i2, LatHet, DNS_RL, GasV_C
    ! Input variables
    integer, intent(in) :: KA,KS,KE
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    integer,  intent(in) :: sdm_aslset
    real(RP), intent(in) :: sdm_aslmw(20)
    real(RP), intent(in) :: sdm_aslion(20)
    real(RP), intent(in) :: sdm_dtevl  ! tims step of {condensation/evaporation} process
    real(RP), intent(in) :: pres_scale(KA,IA,JA) ! Pressure
    real(RP), intent(in) :: t_scale(KA,IA,JA) ! Temperature
    real(RP), intent(in) :: qv_scale(KA,IA,JA)   ! Water vapor mixing ratio
    integer,  intent(in) :: sd_num      ! number of super-droplets
    integer,  intent(in) :: sd_numasl   ! number of kind of chemical material contained as water-soluble aerosol in super droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(in) :: sd_asl(1:sd_num,1:sd_numasl) ! aerosol mass of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
    ! Input and output variables
    real(RP), intent(inout) :: sd_r(1:sd_num) ! equivalent radius of super-droplets

    ! Internal shared variables
    real(RP):: sd_aslmw(1:22) ! Molecular mass of chemical material contained as water-soluble aerosol in super droplets (default+20)
    real(RP):: sd_aslion(1:22) ! Degree of ion dissociation of chemical material contained as water-soluble aerosol in super droplets (default+20)

    ! Work variables
    real(RP) :: dmask(1:22)       ! mask for vactorization
    real(RP) :: ss_sd(KA,IA,JA)   ! Degree of super-saturation
    real(RP) :: Fac_dd            ! Fd in growth EQ.(R.R.Rogers)
    real(RP) :: Fac_kk            ! Fk in growth EQ.(R.R.Rogers)
    real(RP) :: Rc                ! Rc
    real(RP) :: ivt_sd            ! 1.d0 / t_scale
    real(RP) :: dtivfdk(KA,IA,JA) ! dt / ( Fk + Td )
    real(RP) :: crd2              ! sd_r^2 before interation
    real(RP) :: rdi               ! sd_r after interation
    real(RP) :: a                 ! temporary
    real(RP) :: a3                ! temporary
    real(RP) :: b                 ! temporary
    real(RP) :: eq_a              ! temporary for growth EQ.
    real(RP) :: eq_b              ! temporary for growth EQ.
    real(RP) :: eq_c              ! temporary for growth EQ.
    real(RP) :: new_rd            ! temporary for interation
    real(RP) :: rd2               ! temporary for interation
    real(RP) :: rd5               ! temporary for interation
    real(RP) :: rd7               ! temporary for interation
    real(RP) :: dterm1            ! temporary for interation
    real(RP) :: dterm2            ! temporary for interation
    real(RP) :: dterm3            ! temporary for interation
    real(RP) :: dtmp              ! temporary for interation
    real(RP) :: diff              ! temporary for interation
    real(RP) :: rddvcp            ! rd / cp
    integer  :: idx_nasl(1:22)    ! index for vactorization
    real(RP) :: qd_scale(KA,IA,JA)
    real(RP) :: qvs_scale(KA,IA,JA)
    real(RP) :: es_scale(KA,IA,JA)

    ! list
    integer :: list(sd_num)
    integer :: cnt

    integer :: i, j, k, n, m, s, t, it               ! index
    ! Parameters
    integer, parameter :: itr_max = 25   ! iteration number
!!$    real(RP), parameter :: epsva = 0.622_RP   ! Molecular weight ratio of vapor/air

    real(RP) :: L_RL_K_new
    real(RP) :: RLRv_D_new
    real(RP) :: Diff_C_new
    real(RP) :: Heat_C_new

    real(RP) :: tdeg
    real(RP) :: Heat_dry
    real(RP) :: Heat_vapor

    !---------------------------------------------------------------------
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_condevp",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    rddvcp = rd / cp

    !### aerosol type ###!
    
    if( abs(mod(sdm_aslset,10))==1 ) then

       !### numasl=1 @ init+rest : (NH4)2SO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    else if( abs(mod(sdm_aslset,10))==2 ) then

       if( abs(sdm_aslset)==2 ) then
          
          !### numasl=1 @ init : NaCl ###!

          sd_aslmw(1)  = mass_nacl
          sd_aslion(1) = ion_nacl

       else if( abs(sdm_aslset)==12 ) then

          !### numasl=2 @ init : NaCl, rest : (NH4)2SO4 ###!

          sd_aslmw(1) = mass_amsul
          sd_aslmw(2) = mass_nacl
          sd_aslion(1) = ion_amsul
          sd_aslion(2) = ion_nacl

       end if

    else if( abs(mod(sdm_aslset,10))==3 ) then

       !### numasl>=2 @ init+rest : (NH4)2SO4, NaCl, ... ###!

       sd_aslmw(1) = mass_amsul
       sd_aslmw(2) = mass_nacl

       sd_aslion(1) = ion_amsul
       sd_aslion(2) = ion_nacl

       !! Must be a Bug. This cannot be simply commented out
       do n=1,20
       !            call getrname( id_sdm_aslmw  + (n-1), sd_aslmw(n+2)  )
       !            call getrname( id_sdm_aslion + (n-1), sd_aslion(n+2) )
       end do

    else if( abs(mod(sdm_aslset,10))==5 ) then

       !### numasl=1 @ init+rest : (NH4)HSO4 ###!

       sd_aslmw(1)  = mass_amsul
       sd_aslion(1) = ion_amsul

    end if

    do n=1,22

       if( n<=sd_numasl ) then
          idx_nasl(n) = n
          dmask(n) = 1.0_RP
       else
          idx_nasl(n) = sd_numasl
          dmask(n) = 0.0_RP
       end if

    end do


    !$omp workshare
    qd_scale(:,:,:) = 1.0_RP-qv_scale(:,:,:)
    !$omp end workshare

    call ATMOS_SATURATION_pres2qsat_liq( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         t_scale(:,:,:),pres_scale(:,:,:),qd_scale(:,:,:),qvs_scale(:,:,:))
    call ATMOS_SATURATION_psat_liq( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         t_scale(:,:,:), es_scale(:,:,:) )

    ! Condensation of water
!OCL INDEPENDENT
    !$omp  parallel do collapse(2) &
    !$omp& private(i,j,k) &
    !$omp& private(ivt_sd)   &
    !$omp& private(Diff_C_new,tdeg,Heat_dry,Heat_vapor,Heat_C_new,L_RL_K_new,RLRv_D_new) &
    !$omp& private(Fac_dd,Fac_kk)
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE

       !### Calculate degree of super-saturation ###!

       ivt_sd = 1.0_RP / t_scale(k,i,j)
       ss_sd(k,i,j) = qv_scale(k,i,j) / qvs_scale(k,i,j) !! degree of super-saturation


#ifdef DYCOMS2_RF02_SDM
       !! limiter of super saturation for DYCOMSII(RF02)
       if(TIME_NOWSEC<3600.0_RP)then
          ss_sd(k,i,j) = min(ss_sd(k,i,j),1.01_RP)
       end if
#endif

       !### Set parameters for growth EQ. of the radius ###!
       
       !!!! varibales defined in sdm_common.f90
       !  LatGas = LatHet / GasV_C
       !  L_RL_K = LatHet * DNS_RL / Heat_C
       !  RLRv_D = DNS_RL * GasV_C / Diff_C
       !  ASL_RR = ASL_FF * ION_asl / WGTasl

       !!!!!! obsolete formula ignoring the T and P dependencies
       !! Fac_dd  = RLRv_D * t_scale(k,i,j) / es_scale(k,i,j)
       !! Fac_kk  = ( LatGas*ivt_sd - 1.0_RP ) * L_RL_K * ivt_sd

       !!!!!!! Diffusion coefficient [m^2/s] (Beard and Pruppacher, 1971) !!!!!!!
       Diff_C_new = ( 0.211_RP * (( t_scale(k,i,j) / t0 ) ** 1.94_RP )*( 1013.25E2_RP /pres_scale(k,i,j))) * 1.0E-4_RP

       !!!!!! Thermal conductivity [J/(m s K)] (Beard and Pruppacher, 1971) !!!!!!!!!
       tdeg = t_scale(k,i,j) - t0
       Heat_dry   = (( 5.69_RP + 1.68e-2_RP * tdeg ) * 1.0e-5_RP )  ! Thermal conductivity of dry air [cal/(cm*s*C)]
       Heat_vapor = (( 3.73_RP + 2.00e-2_RP * tdeg ) * 1.0e-5_RP )  ! Thermal conductivity of vapor [cal/(cm*s*C)]
       ! Thermal conductivity of moist air [J/(m*s*K)]
       ! It is not clear which definition of cal was used in Beard and Pruppacher (1971)
       ! We adopted the thermochemical calorie: 1 cal_th = 4.184 J
       Heat_C_new = ( Heat_dry * ( 1.0_RP - ( 1.17_RP - 1.02_RP * ( Heat_vapor / Heat_dry )) * qv_scale(k,i,j) )) * 4.184_RP * 1.0e2_RP
       
       L_RL_K_new = LatHet * DNS_RL / Heat_C_new
       RLRv_D_new = DNS_RL * GasV_C / Diff_C_new

       Fac_dd  = RLRv_D_new * t_scale(k,i,j) / es_scale(k,i,j)
       Fac_kk  = ( LatGas*ivt_sd - 1.0_RP ) * L_RL_K_new * ivt_sd

       dtivFdk(k,i,j) = sdm_dtevl / ( Fac_dd + Fac_kk )

    end do
    end do
    end do

    cnt = 0
    do n = 1, sd_num
       ! only liquid droplets
       if( sd_rk(n) > VALID2INVALID .and. sd_liqice(n) == STAT_LIQ ) then
          cnt = cnt + 1
          list(cnt) = n
       end if
    end do

    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(i,j,k,s,t,it,n) &
    !$omp& private(a,b,a3,eq_a,eq_b,eq_c)                                &
    !$omp& private(crd2,new_rd,rdi,rd2,rd5,rd7,dtmp,dterm1,dterm2,dterm3,diff) &
    !$omp& private(Rc)
    do m = 1, cnt

       n = list(m)

       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       eq_a = CurveF / t_scale(k,i,j)

       eq_b = 0.0_RP

!OCL FULLUNROLL_PRE_SIMD
       do t=1,22

          s = idx_nasl(t)

          dtmp = sd_asl(n,s) * ( sd_aslion(s) / sd_aslmw(s) )
          eq_b = eq_b + dmask(t) * dtmp

       end do

       eq_b = eq_b * ASL_FF

       Rc    = sqrt(eq_b/eq_a)
       eq_c  = ss_sd(k,i,j) - 1.0_RP

       a  = eq_a / eq_c
       b  = eq_b / eq_c
       a3 = a * a * a

       !### advance the particle radius ###!

       crd2   = sd_r(n)**2
       new_rd = sd_r(n)

       if( ( ss_sd(k,i,j) > 1.0_RP ) .and. ( a3 < b*(27.0_RP/4.0_RP) ) ) then
          new_rd = 1.0E-3_RP
       end if

       new_rd = max(new_rd, Rc)

       !== iteration ( newton-raphson method ) ==!

       rdi = new_rd
       rd2 = rdi**2
!OCL FULLUNROLL_PRE_SIMD
       do it=1,itr_max

          ! Considering Kohler-curve (R.R.Rogers)

          rd5 = rd2 * rd2 * rdi
          rd7 = rd5 * rd2

          dterm1 = dtivFdk(k,i,j) * ( rd5*eq_c + (eq_b-eq_a*rd2)*rd2 )
          dterm2 = rd5 * crd2
          dterm3 = rd5 - dtivFdk(k,i,j) * ( -3.0_RP*eq_b + eq_a*rd2 )

          diff = ( rd7 - 2.0_RP*dterm1 - dterm2 ) / dterm3
          ! if ( abs(diff) < rd2 * 1e-8_RP ) exit

          rd2 = max(rd2 - diff, rd2 * 1e-4_RP)
          rdi = sqrt(rd2)

       end do

       sd_r(n) = rdi   !! particle radius at future

    end do
    !$omp end parallel do simd

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_condevp",1,1)
#endif
    return
  end subroutine sdm_condevp
  !-----------------------------------------------------------------------------
  subroutine sdm_condevp_updatefluid(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       RHOT,QTRC,DENS,rhowp_sdm,rhowc_sdm)
    use scale_atmos_hydrometeor, only: &
        ATMOS_HYDROMETEOR_LHV
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2cpexnr, &
         sdm_rhot_qtrc2p_t
    use m_sdm_common, only: &
         LatHet,I_QV_sdm,QA_MP_sdm
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: rhowp_sdm(KA,IA,JA) ! density of liquid water at preveous step
    real(RP), intent(in) :: rhowc_sdm(KA,IA,JA) ! density of liquid water at current step
    ! Output variables
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    ! Work variables
    real(RP) :: rhov(KA,IA,JA),cpexnr(KA,IA,JA)
    real(RP) :: delta_rhow, dens_old
    integer  :: i, j, k ! index
    real(RP) :: lhv(KA,IA,JA), pre(KA,IA,JA), temp(KA,IA,JA)
    !-------------------------------------------------------------------7--
     
    ! exchange vapor
    !$omp parallel
    !$omp do private(j,i,k) collapse(3)
    do j = 1, JA
    do i = 1, IA
    do k = 1, KA
       rhov(k,i,j)=DENS(k,i,j)*QTRC(k,i,j,I_QV_sdm)
    end do
    end do
    end do
    !$omp end do
    !$omp do private(j,i,k,delta_rhow,dens_old) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       delta_rhow = rhowc_sdm(k,i,j)-rhowp_sdm(k,i,j)
       dens_old = DENS(k,i,j)
       
       DENS(k,i,j) = DENS(k,i,j) - delta_rhow
       rhov(k,i,j) = rhov(k,i,j) - delta_rhow
       QTRC(k,i,j,I_QV_sdm) = rhov(k,i,j) / DENS(k,i,j)
       RHOT(k,i,j) = RHOT(k,i,j) * (DENS(k,i,j)/dens_old)  
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    call sdm_rhot_qtrc2p_t( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         RHOT,QTRC,DENS,pre,temp)
    call ATMOS_HYDROMETEOR_LHV( &
         KA, KS, KE, &
         IA, IS, IE, &
         JA, JS, JE, &
         temp,lhv)

    ! exchange heat
    call sdm_rhot_qtrc2cpexnr( &
         KA, IA, JA, &
         RHOT,QTRC,DENS,cpexnr)

    !$omp parallel
    !$omp do private(j,i,k,delta_rhow) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       delta_rhow = rhowc_sdm(k,i,j)-rhowp_sdm(k,i,j)
       
!!$       RHOT(k,i,j) = RHOT(k,i,j) + LatHet*delta_rhow / cpexnr(k,i,j)  
       RHOT(k,i,j) = RHOT(k,i,j) + lhv(k,i,j)*delta_rhow / cpexnr(k,i,j)  
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_condevp_updatefluid
end module m_sdm_condensation_water
