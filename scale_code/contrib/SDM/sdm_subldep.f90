!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics / SDM
!!
!! @par Description
!!          Sublimation and deposition water vapor to ice particles
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
!! @li      2016-07-19 (S.Shima) [new] Newly created
!! @li      2016-07-20 (S.Shima) [fix] Correct a wrong equation: ** -> *
!! @li      2017-02-15 (S.Shima) [mod] implement sublimation/deposition growth model of Chen and Lamb (1994)
!! @li      2017-02-15 (S.Shima) [mod] list vector is introduced to improve performance
!! @li      2017-02-15 (S.Shima) [fix] depsition density formula is used only for deposition, not for sublimation
!! @li      2017-05-04 (S.Shima) [add] sdm_tdeg2growthratio
!! @li      2017-09-27 (S.Shima) [fix] serious bug of sdm_tdeg2growthratio
!! @li      2017-10-01 (S.Shima) [mod] calculate Reynolds number from the terminal velocity
!! @li      2017-10-15 (S.Shima) [mod] limiter to avoid unphysical small ice particle
!! @li      2017-11-28 (S.Shima) [mod] When calculating rho_dep, rhov_infty is limited by water saturation, following Miller and Young (1979).
!! @li      2018-03-29 (S.Shima) [mod] deposition density of Jensen and Harrington (2015)
!! @li      2018-03-29 (S.Shima) [mod] set inherent growth ratio = 1 if D < 10um
!! @li      2018-05-25 (S.Shima) [fix] deposition density of Jensen and Harrington (2015)
!! @li      2018-06-30 (S.Shima) [add] calculation of rime mass
!! @li      2019-01-10 (S.Shima) [mod] deposition density of Chen and Lamb (1994) but use rho_bulk for small plate (Jensen and Harrington, 2015)
!! @li      2019-01-12 (S.Shima) [mod] dry air density -> moist air density
!! @li      2019-06-08 (S.Shima) [mod] set inherent growth ratio = 1 for sublimation
!! @li      2020-05-24 (S.Shima) [mod] if already too long, limit Gamma_star<=1 for deposition
!!
!< 
!-------------------------------------------------------------------------------
module m_sdm_subldep
  use scale_precision  
  implicit none
  private
  public :: sdm_subldep,sdm_subldep_updatefluid

contains
  subroutine sdm_subldep(KA,KS,KE, IA,IS,IE, JA,JS,JE, &
       sdm_dtsbl, &
       pres_scale,t_scale,qv_scale,rhom_scale, &
       sd_num,sd_liqice,sd_x,sd_y,     &
       sdi,sd_vz,sd_ri,sd_rj,sd_rk,ilist        )
    
    use scale_const, only: &
         rhoi_mks => CONST_DICE, & ! density of ice [kg/m^3]
         t0       => CONST_TEM00   ! 0 degC in [K]
    use scale_atmos_saturation, only: &
        ATMOS_SATURATION_pres2qsat_liq, &
        ATMOS_SATURATION_pres2qsat_ice, &
        ATMOS_SATURATION_psat_ice
    use m_sdm_coordtrans, only: &
         sdm_x2ri, sdm_y2rj
    use m_sdm_common, only: &
         F_THRD, ONE_PI, &
         VALID2INVALID, STAT_ICE, i2, sdicedef, &
         ! Diff_C, & ! diffusion constant of vapor [m^2/s]  
         GasV_C, & ! Gas Constant of vapor [J/K/kg]
         ! Heat_C, & ! thermal conductivity at 293K, 100kPa [J/m s K]
         LatHet_S  ! latent heat of sublimation at 0C [J/kg]
    ! Input variables
    integer, intent(in) :: KA,KS,KE
    integer, intent(in) :: IA,IS,IE
    integer, intent(in) :: JA,JS,JE
    real(RP), intent(in) :: sdm_dtsbl  ! tims step of {sublimation/deposition} process [s]
    real(RP), intent(in) :: pres_scale(KA,IA,JA) ! Pressure [Pa]
    real(RP), intent(in) :: t_scale(KA,IA,JA) ! Temperature [K]
    real(RP), intent(in) :: qv_scale(KA,IA,JA)   ! Water vapor mixing ratio [kg/kg]
    real(RP), intent(in) :: rhom_scale(KA,IA,JA)  ! moist air density [kg/m^3]
    integer,  intent(in) :: sd_num      ! number of super-droplets
    integer(i2), intent(in) :: sd_liqice(1:sd_num)
                       ! status of super-droplets (liquid/ice)
                       ! 01 = all liquid, 10 = all ice
                       ! 11 = mixture of ice and liquid
    real(RP), intent(in) :: sd_x(1:sd_num) ! x-coordinate of super-droplets
    real(RP), intent(in) :: sd_y(1:sd_num) ! y-coordinate of super-droplets
    real(RP), intent(inout) :: sd_ri(1:sd_num)   ! index[i/real] of super-droplets
    real(RP), intent(inout) :: sd_rj(1:sd_num)   ! index[j/real] of super-droplets
    real(RP), intent(in) :: sd_rk(1:sd_num)   ! index[k/real] of super-droplets
    real(RP), intent(in) :: sd_vz(1:sd_num)! terminal velocity of super-droplets [m/s]
    ! Input and output variables
    type(sdicedef), intent(inout) :: sdi   ! ice phase super-droplets
    ! Output variables
    integer, intent(out) :: ilist(1:sd_num)  ! buffer for list vectorization

    ! Work variables
    real(RP) :: p_sd      ! pressure of the grid contained the SD [Pa]
    real(RP) :: t_sd      ! temperature of the grid contained the SD [K]
    real(RP) :: qd_sd     ! dry air mass ratio at the grid contained the SD [kg/kg]
    real(RP) :: qv_sd     ! water-vapor of the grid contained the SD [kg/kg]
    real(RP) :: qvsi_sd   ! Saturation mixing ratio over ice of the grid contained the SD [kg/kg]
    real(RP) :: esi_sd    ! Saturation vapor pressure over ice of the grid contained the SD [kg/kg]
    real(RP) :: satri_sd  ! Degree of super-saturation over ice of the grid contained the SD 
    real(RP) :: rhom_sd   ! Moist air density at the position of the SD [kg/m3]
    real(RP) :: dvisc_sd  ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
    real(RP) :: barf      ! mean ventilation effect
    real(RP) :: fratio    ! ventilation effect ratio fc/fa
    real(RP) :: ssi       ! super saturation
    real(RP) :: Fac_ice_dd_invrhoi ! Fd/rhoi
    real(RP) :: Fac_ice_kk_invrhoi ! Fk/rhoi
    real(RP) :: dtrhoiivFdk_ice ! dt*rhoi/(Fd+Fk)
    real(RP) :: ivt_sd    ! 1.d0 / t_sd
    real(RP) :: capaci    ! electric capacitance of the SD
    real(RP) :: ecent     ! eccentricity as a spheroid of the SD 
    real(RP) :: qv_sd_max ! for calculating rho_dep (max of vapor is limited by water saturation)
    real(RP) :: satri_sd_max ! for calculating rho_dep (max of satration ratio over ice is limited by water saturation)
    real(RP) :: ssi_max   ! for calculating rho_dep (max of super saturation over ice is limited by water saturation)

    real(RP) :: sd_ra     ! current equatorial radius of the SD [m]
    real(RP) :: sd_rc     ! current polar radius of the SD [m]
    real(RP) :: sd_rho    ! current density of the SD [kg/m^3]
    real(RP) :: sd_tvel   ! current terminal velocity of the SD [m/s^2]
    real(RP) :: sd_mass   ! current mass of the SD [kg]
    real(RP) :: sd_mass_two_thrd   ! (sd_mass)^(2/3)
    real(RP) :: sd_phi    ! sd_rc/sd_ra: aspect ratio
    real(RP) :: sd_vol    ! current volume of the SD [m^3]
    real(RP) :: sd_maxD   ! maximum dimension of the particle [m]
    real(RP) :: sd_mrime  ! current rime mass of the particle [kg]
    real(RP) :: new_sd_vol    ! next volume of the SD
    real(RP) :: new_sd_ra     ! next equatorial radius of the SD
    real(RP) :: new_sd_rc     ! next polar radius of the SD
    real(RP) :: new_sd_rho    ! next density of the SD
    real(RP) :: new_sd_mass   ! next mass of the SD
    real(RP) :: new_sd_mass_two_thrd   ! (new_sd_mass)^(2/3)
    real(RP) :: new_sd_mrime   ! next rime mass of the SD
    real(RP) :: delta_mass   ! new_sd_mass - sd_mass
    real(RP) :: delta_vol   ! new_sd_vol - sd_vol
    real(RP) :: delta_logvol ! log(new_sd_vol) - log(sd_vol)
    real(RP) :: delta_logra  ! log(new_sd_ra) - log(sd_ra)
    real(RP) :: rho_dep   ! deposition density [kg/m^3]
    real(RP) :: growth_ratio ! inherent growth ratio Gamma(T)
    real(RP) :: growth_ratio_star ! ventilation corrected inherent growth ratio Gamma*(T) = Gamma(T)*(fc/fa)
    real(RP) :: tdeg_sd ! temperature in degree C
    real(RP) :: excess_vap_dens ! excess vapor density [kg/m^3]
    real(RP) :: n_ventX ! X = Nsc^(1/3)*Nre^(1/2) for calculating the ventilation effect
    real(RP) :: n_re    ! Reynolds number
    real(RP) :: n_sc    ! Schmidt number
    real(RP) :: qd_scale(KA,IA,JA), qvsi_scale(KA,IA,JA), esi_scale(KA,IA,JA), qvsl_scale(KA,IA,JA)

    integer :: i, j, k, n, m             ! index

    integer :: tlist      ! total list number for ice
    integer :: icnt       ! counter for ice

    real(RP) :: Diff_C_new
    real(RP) :: Heat_C_new

    real(RP) :: Heat_dry
    real(RP) :: Heat_vapor

    !---------------------------------------------------------------------
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_subldep",1,1)
#endif
    call sdm_x2ri(IS,IE,sd_num,sd_x,sd_ri,sd_rk)
    call sdm_y2rj(JS,JE,sd_num,sd_y,sd_rj,sd_rk)

    ! Initialize
    tlist = 0

    ! Get index list for compressing buffer.
    icnt = 0

    do n=1,sd_num
       if( sd_rk(n)<VALID2INVALID ) cycle
       if( sd_liqice(n) .ne. STAT_ICE ) cycle

       !== ice particles ==!
       icnt = icnt + 1
       ilist(icnt) = n
       
    end do

    tlist = icnt

    ! Deposition and sublimation of water
    if( tlist>0 ) then
    
    !$omp parallel workshare
    qd_scale(:,:,:) = 1.0_RP-qv_scale(:,:,:)
    !$omp end parallel workshare

    call ATMOS_SATURATION_pres2qsat_ice( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         t_scale(:,:,:),pres_scale(:,:,:),qd_scale(:,:,:),qvsi_scale(:,:,:))
    call ATMOS_SATURATION_psat_ice( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         t_scale(:,:,:),esi_scale(:,:,:))
    call ATMOS_SATURATION_pres2qsat_liq( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         t_scale(:,:,:),pres_scale(:,:,:),qd_scale(:,:,:),qvsl_scale(:,:,:))

!OCL NORECURRENCE,INDEPENDENT
    !$omp  parallel do simd schedule(static,256) &
    !$omp& private(n,i,j,k) &
    !$omp& private(p_sd,t_sd,qv_sd,rhom_sd,tdeg_sd) &
    !$omp& private(sd_ra,sd_rc,sd_rho,sd_mrime,sd_tvel,sd_phi,sd_vol,sd_mass,sd_maxD) &
    !$omp& private(growth_ratio,ecent,capaci,dvisc_sd,Diff_C_new,Heat_dry,Heat_vapor,Heat_C_new) &
    !$omp& private(n_re,n_sc,n_ventX) &
    !$omp& private(qd_sd,qvsi_sd,satri_sd,ssi) &
    !$omp& private(ivt_sd,Fac_ice_kk_invrhoi) &
    !$omp& private(esi_sd,Fac_ice_dd_invrhoi,barf,dtrhoiivFdk_ice) &
    !$omp& private(sd_mass_two_thrd,new_sd_mass_two_thrd,new_sd_mass,delta_mass) &
    !$omp& private(rho_dep,qv_sd_max,satri_sd_max,ssi_max,excess_vap_dens,delta_vol,new_sd_vol) &
    !$omp& private(new_sd_rho,fratio,growth_ratio_star,delta_logvol,delta_logra,new_sd_ra,new_sd_rc) &
    !$omp& private(new_sd_mrime)
    do m=1,tlist
       n = ilist(m)

       !### Atmospheric condition at the location of the super-droplet ###!
       
       i = floor(sd_ri(n))+1
       j = floor(sd_rj(n))+1
       k = floor(sd_rk(n))+1

       p_sd  = pres_scale(k,i,j)
       t_sd  = t_scale(k,i,j)
       qv_sd = qv_scale(k,i,j)
       rhom_sd = rhom_scale(k,i,j)
       qvsi_sd = qvsi_scale(k,i,j)
       esi_sd = esi_scale(k,i,j)

       tdeg_sd = t_sd - t0     !! [K] => [degC]

       !### Save the current state of the super-droplet ###!
       sd_ra  = sdi%re(n)
       sd_rc  = sdi%rp(n)
       sd_rho = sdi%rho(n)
       sd_mrime = sdi%mrime(n)
       sd_tvel = sd_vz(n)

       sd_phi  = sd_rc/sd_ra
       sd_vol  = F_THRD * ONE_PI * sd_ra**2 * sd_rc
       sd_mass = sd_vol * sd_rho ! mass of the particle [kg]
       sd_maxD = 2.0d0 * max(sd_ra,sd_rc) ! maximum dimension of the particle [m]

       !### Inherent growth ratio (Chen and Lamb 1994) ###!
       call sdm_tdeg2growthratio(tdeg_sd,growth_ratio)
       if(sd_maxD<1.0e-5) growth_ratio=1.0d0 !! set inherent growth ratio = 1 if D < 10um (Shima)

       !###### capacitance ######!
       if( sd_phi >= 1.001d0 ) then ! prolate
          ecent = sqrt(1.0d0 - 1.0d0/(sd_phi**2))
          capaci = ecent*sd_rc / log((1.0d0+ecent)*sd_phi)            
       else if( sd_phi <= 0.999d0 ) then ! oblate         
          ecent = sqrt(1.0d0 - sd_phi**2)
          capaci = ecent*sd_ra / asin(ecent)                  
       else ! approximation around sd_phi=1.0
          capaci = (2.0d0*sd_ra +sd_rc)/3.0d0
       end if

       !### Calculate X = Nsc^(1/3)*Nre^(1/2) for ventilation ###!
       ! dynamic viscosity [kg/(m*s)] at the location of super-droplets
       !== (Pruppacher & Klett,1997) ==!
       if( tdeg_sd>=0.0d0 ) then
          dvisc_sd = ( 1.7180_RP + 4.9E-3_RP*tdeg_sd ) * 1.E-5_RP
       else
          dvisc_sd = ( 1.718_RP + 4.9E-3_RP*tdeg_sd              &
               - 1.2E-5_RP*tdeg_sd*tdeg_sd ) * 1.E-5_RP
       end if

       !!!!!!! Diffusion coefficient [m^2/s] (Beard and Pruppacher, 1971) !!!!!!!
       Diff_C_new = ( 0.211_RP * (( t_sd / t0 ) ** 1.94_RP )*( 1013.25E2_RP /p_sd)) * 1.0E-4_RP

       !!!!!! Thermal conductivity [J/(m s K)] (Beard and Pruppacher, 1971) !!!!!!!!!
       Heat_dry   = (( 5.69d0 + 1.68d-2 * tdeg_sd ) * 1.0d-5 )  ! Thermal conductivity of dry air [cal/(cm*s*C)]
       Heat_vapor = (( 3.73d0 + 2.00d-2 * tdeg_sd ) * 1.0d-5 )  ! Thermal conductivity of vapor [cal/(cm*s*C)]
       ! Thermal conductivity of moist air [J/(m*s*K)]
       ! It is not clear which definition of cal was used in Beard and Pruppacher (1971)
       ! We adopted the thermochemical calorie: 1 cal_th = 4.184 J
       Heat_C_new = ( Heat_dry * ( 1.0d0 - ( 1.17d0 - 1.02d0 * ( Heat_vapor / Heat_dry )) * qv_sd )) * 4.184d0 * 1.0d+2

       ! Reynolds number of the particle
       n_re = rhom_sd * sd_maxD * sd_tvel / dvisc_sd

       ! Schmidt number of the particle
       n_sc = dvisc_sd / rhom_sd / Diff_C_new

       ! Nsc^(1/3)*Nre^(1/2)
       n_ventX = n_sc**(1.0d0/3.0d0) * sqrt(n_re)

       !#################################################!
       !### Mass change ###!
       !###### ice saturation ratio ######!
       satri_sd  = qv_sd / qvsi_sd
       ssi  = satri_sd - 1.0_RP

       !###### kinetic correction ######!
       ivt_sd = 1.0_RP / t_sd
       Fac_ice_kk_invrhoi = ( ivt_sd * LatHet_S / GasV_C - 1.0_RP ) * (LatHet_S/Heat_C_new) * ivt_sd

       !###### psychrometric correction ######!
       Fac_ice_dd_invrhoi = GasV_C / Diff_C_new * t_sd / esi_sd

       !###### ventilation effect ######!
       if(n_ventX <= 1.0d0) then
          barf = 1.0d0 + 0.14d0*n_ventX**2  
       else
          barf = 0.86d0 + 0.28d0*n_ventX  
       end if

       !###### delta mass ######!
       dtrhoiivFdk_ice = real( sdm_dtsbl, kind=RP ) / ( Fac_ice_dd_invrhoi + Fac_ice_kk_invrhoi )
       
       sd_mass_two_thrd = sd_mass**(2.0d0/3.0d0)
       
       ! explicit Euler for m^(2/3)
       new_sd_mass_two_thrd = sd_mass_two_thrd &
            + (2.0d0/3.0d0) * (sd_mass**(-1.0d0/3.0d0)) * 4.0d0*ONE_PI*capaci*barf * dtrhoiivFdk_ice * ssi

       ! avoid negative mass. (18/(6.02*10^23)/1000)^(2/3)=10^17
       new_sd_mass_two_thrd = max(new_sd_mass_two_thrd, 1.E-17_RP)

       new_sd_mass = new_sd_mass_two_thrd**(3.0d0/2.0d0) 

       delta_mass = new_sd_mass - sd_mass

       !#################################################!
       !### Volume change ###!
       if(delta_mass > 0.0d0) then ! deposition

          ! deposition density [kg/m^3] (Chen and Lamb 1994, Fukuta 1969)
          !! Assume true density of ice if the planer ice is small (Jensen and Harrington 2015)
          if((growth_ratio<1.0_RP).and.(sd_ra<1.0e-4_RP))then
             rho_dep = rhoi_mks
          else
             !! max of ssi ( = super saturation over ice at water saturation )
             qv_sd_max = qvsl_scale(k,i,j)
             satri_sd_max  = qv_sd_max / qvsi_sd
             ssi_max  = satri_sd_max - 1.0_RP
             !! excess vapor density [g/m^3]: rhov_infty - rhov_surface, but limited by water saturation
             excess_vap_dens = 1000.0d0 * min(ssi,ssi_max) / ( Fac_ice_dd_invrhoi + Fac_ice_kk_invrhoi ) / Diff_C_new
             rho_dep = 0.91d3 * exp(-3.0d0*max(excess_vap_dens-0.05d0,0.0d0)/growth_ratio)
          end if

!!$          ! deposition density [kg/m^3] (Jensen and Harrington 2015)
!!$          if(growth_ratio<1.0_RP)then
!!$             if(sd_ra<1.0e-4_RP)then
!!$                rho_dep = rhoi_mks
!!$             else
!!$                rho_dep = rhoi_mks * growth_ratio
!!$             end if
!!$          else
!!$             rho_dep = rhoi_mks / growth_ratio
!!$          end if

       else ! sublimation
          rho_dep = sd_rho
          growth_ratio=1.0d0 !! set inherent growth ratio = 1 for sublimation
       end if
       
       delta_vol = delta_mass / rho_dep
       new_sd_vol = sd_vol + delta_vol

       ! avoid negative volume and too dense particle
       new_sd_vol = max(new_sd_vol, new_sd_mass/rhoi_mks)
       delta_vol = new_sd_vol - sd_vol

       !#################################################!
       !### Density change ###!
       new_sd_rho = new_sd_mass / new_sd_vol

       !#################################################!
       !### Radius change ###!
       !###### ventilation effect ######!
       if(n_ventX <= 1.0d0) then
          fratio =  (1.0d0 + 0.14d0*n_ventX**2*sqrt(sd_rc/capaci)) &
               &  / (1.0d0 + 0.14d0*n_ventX**2*sqrt(sd_ra/capaci)) 
       else
          fratio =  (0.86d0 + 0.28d0*n_ventX*sqrt(sd_rc/capaci)) &  
               &  / (0.86d0 + 0.28d0*n_ventX*sqrt(sd_ra/capaci))
       end if

       growth_ratio_star = growth_ratio * fratio

       !! if already too long, limit Gamma_star<=1 for deposition
       if( (delta_mass > 0.0d0) .and. (sd_phi>40.0_RP) )then
          growth_ratio_star = min(growth_ratio_star,1.0_RP)
       end if

       !###### Radius change ######!
       delta_logvol = log(new_sd_vol/sd_vol)

       delta_logra = delta_logvol/(growth_ratio_star + 2.0d0) 
       new_sd_ra = sd_ra*exp(delta_logra)
       new_sd_rc = sd_rc*exp(growth_ratio_star*delta_logra)

       !#################################################!
       !### limiter to avoid unphysical ice particle ###!
       if( min(new_sd_ra,new_sd_rc) <= 1.0d-6 )then ! Regard it is spherical with bulk ice density if smaller than 1um
          new_sd_rho = rhoi_mks
          new_sd_vol = new_sd_mass / new_sd_rho
          new_sd_ra = (new_sd_vol/F_THRD/ONE_PI)**(1.0d0/3.0d0) 
          new_sd_rc = new_sd_ra
       end if

       if( min(new_sd_ra,new_sd_rc) <= 1.0d-9 )then ! Stop sublimation if smaller than 1nm
          new_sd_ra = 1.0d-9
          new_sd_rc = 1.0d-9
          new_sd_rho = rhoi_mks
       end if

       !#################################################!
       !### rime mas change ###!
       if(delta_mass > 0.0d0) then ! deposition
          new_sd_mrime = sd_mrime
       else ! sublimation
          new_sd_mrime = sd_mrime * new_sd_mass/sd_mass
       end if

       !#################################################!
       !### update the super-droplet ###!
       sdi%re(n)  = new_sd_ra
       sdi%rp(n)  = new_sd_rc
       sdi%rho(n) = new_sd_rho
       sdi%mrime(n) = new_sd_mrime

    end do
    !$omp end parallel do simd
    end if

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_subldep",1,1)
#endif
    return
  end subroutine sdm_subldep
  !-----------------------------------------------------------------------------
  subroutine sdm_subldep_updatefluid(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       RHOT,QTRC,DENS,rhosolp_sdm,rhosolc_sdm)
    use scale_atmos_hydrometeor, only: &
        ATMOS_HYDROMETEOR_LHS
    use m_sdm_fluidconv, only: &
         sdm_rhot_qtrc2cpexnr, &
         sdm_rhot_qtrc2p_t
    use m_sdm_common, only: &
         LatHet,I_QV_sdm,QA_MP_sdm
    ! Input variables
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in) :: rhosolp_sdm(KA,IA,JA) ! density of liquid water at preveous step
    real(RP), intent(in) :: rhosolc_sdm(KA,IA,JA) ! density of liquid water at current step
    ! Output variables
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    ! Work variables
    real(RP) :: rhov(KA,IA,JA),cpexnr(KA,IA,JA)
    real(RP) :: delta_rhosol, dens_old
    integer  :: i, j, k ! index
    real(RP) :: lhs(KA,IA,JA), pre(KA,IA,JA), temp(KA,IA,JA)
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

    !$omp do private(j,i,k,delta_rhosol,dens_old) collapse(2)
    do j=JS,JE
    do i=IS,IE
    do k=KS,KE
       delta_rhosol = rhosolc_sdm(k,i,j)-rhosolp_sdm(k,i,j)
       dens_old = DENS(k,i,j)
       
       DENS(k,i,j) = DENS(k,i,j) - delta_rhosol
       rhov(k,i,j) = rhov(k,i,j) - delta_rhosol
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
    call ATMOS_HYDROMETEOR_LHS( &
         KA, KS, KE, &
         IA, IS, IE, &
         JA, JS, JE, &
         temp,lhs)

    ! exchange heat
    call sdm_rhot_qtrc2cpexnr(&
         KA, IA, JA, &
         RHOT,QTRC,DENS,cpexnr)

    !$omp parallel
    !$omp do private(j,i,k,delta_rhosol) collapse(2)
    do k=KS,KE
    do j=JS,JE
    do i=IS,IE
       delta_rhosol = rhosolc_sdm(k,i,j)-rhosolp_sdm(k,i,j)
       
!!$       RHOT(k,i,j) = RHOT(k,i,j) + LatHet*delta_rhosol / cpexnr(k,i,j)  
       RHOT(k,i,j) = RHOT(k,i,j) + lhs(k,i,j)*delta_rhosol / cpexnr(k,i,j)  
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine sdm_subldep_updatefluid
  !-----------------------------------------------------------------------------
  subroutine sdm_tdeg2growthratio(tdeg,growth_ratio)
    use m_sdm_common, only: tb_habit
    ! Input variables
    real(RP), intent(in) :: tdeg ! temperature in degree C
    ! Output variables
    real(RP), intent(out) :: growth_ratio ! inherent growth ratio Gamma(T)

    integer n

    n = min(max(nint(1.375-4.0d0*tdeg),1),121)
    growth_ratio = tb_habit(n)

  end subroutine sdm_tdeg2growthratio
end module m_sdm_subldep
