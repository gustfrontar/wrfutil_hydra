!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Aerosol Microphysics
!!
!! @par Description
!!          Aerosol Microphysics driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_ae_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_driver_tracer_setup
  public :: ATMOS_PHY_AE_driver_setup
  public :: ATMOS_PHY_AE_driver_finalize
  public :: ATMOS_PHY_AE_driver_adjustment
  public :: ATMOS_PHY_AE_driver_calc_tendency

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
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_driver_tracer_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_tracer, only: &
       TRACER_regist
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_tracer_setup, &
       ATMOS_PHY_AE_kajino13_NAME, &
       ATMOS_PHY_AE_kajino13_DESC, &
       ATMOS_PHY_AE_kajino13_UNIT
    use scale_atmos_phy_ae_SPRINTARS, only: &
       ATMOS_PHY_AE_SPRINTARS_tracer_setup, &
       ATMOS_PHY_AE_SPRINTARS_NAME,  &
       ATMOS_PHY_AE_SPRINTARS_DESC,  &
       ATMOS_PHY_AE_SPRINTARS_UNIT,  &
       ATMOS_PHY_AE_SPRINTARS_CV,    &
       ATMOS_PHY_AE_SPRINTARS_CP,    &
       ATMOS_PHY_AE_SPRINTARS_R,     &
       ATMOS_PHY_AE_SPRINTARS_ENGI0, &
       ATMOS_PHY_AE_SPRINTARS_ADVC,  &
       ATMOS_PHY_AE_SPRINTARS_MASS
    use scale_prc, only: &
       PRC_abort
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    implicit none

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_driver_tracer_setup",*) 'Setup'

    if ( ATMOS_sw_phy_ae ) then
       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'OFF', 'NONE' )
          LOG_INFO("ATMOS_PHY_AE_driver_tracer_setup",*) 'this component is never called.'
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_tracer_setup( QA_AE ) ! [OUT]

          call TRACER_regist( QS_AE,                         & ! [OUT]
                              QA_AE,                         & ! [IN]
                              ATMOS_PHY_AE_kajino13_NAME(:), & ! [IN]
                              ATMOS_PHY_AE_kajino13_DESC(:), & ! [IN]
                              ATMOS_PHY_AE_kajino13_UNIT(:)  ) ! [IN]

       case ( 'SPRINTARS' )
          call ATMOS_PHY_AE_SPRINTARS_tracer_setup( QA_AE ) ! [OUT]

          call TRACER_regist( QS_AE,                                   & ! [OUT]
                              QA_AE,                                   & ! [IN]
                              ATMOS_PHY_AE_SPRINTARS_NAME(:),          & ! [IN]
                              ATMOS_PHY_AE_SPRINTARS_DESC(:),          & ! [IN]
                              ATMOS_PHY_AE_SPRINTARS_UNIT(:),          & ! [IN]
                              CV    = ATMOS_PHY_AE_SPRINTARS_CV   (:), & ! [IN]
                              CP    = ATMOS_PHY_AE_SPRINTARS_CP   (:), & ! [IN]
                              R     = ATMOS_PHY_AE_SPRINTARS_R    (:), & ! [IN]
                              ENGI0 = ATMOS_PHY_AE_SPRINTARS_ENGI0(:), & ! [IN]
                              ADVC  = ATMOS_PHY_AE_SPRINTARS_ADVC (:), & ! [IN]
                              MASS  = ATMOS_PHY_AE_SPRINTARS_MASS (:)  ) ! [IN]

       case ( 'OFFLINE' )
          LOG_INFO("ATMOS_PHY_AE_driver_tracer_setup",*) 'offline aerosol module has no tracers'
       case default
          LOG_ERROR("ATMOS_PHY_AE_driver_tracer_setup",*) 'invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_abort
       end select

       QE_AE = QS_AE + QA_AE - 1

    else
       QA_AE = 0
       QS_AE = -1
       QE_AE = -2
    end if

    return
  end subroutine ATMOS_PHY_AE_driver_tracer_setup

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_driver_setup
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
        ATMOS_PHY_AE_kajino13_setup
    use scale_atmos_phy_ae_offline, only: &
        ATMOS_PHY_AE_offline_setup
    use scale_atmos_phy_ae_SPRINTARS, only: &
        ATMOS_PHY_AE_SPRINTARS_setup
    use scale_prc, only: &
       PRC_abort
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_ae ) then

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_setup
       case ( 'OFFLINE' )
          call ATMOS_PHY_AE_offline_setup
       case ( 'SPRINTARS' )
          call ATMOS_PHY_AE_SPRINTARS_setup
       case default
          LOG_ERROR("ATMOS_PHY_AE_driver_setup",*) 'invalid aerosol type(', ATMOS_PHY_AE_TYPE, '). CHECK!'
          call PRC_abort
       end select

    endif

    return
  end subroutine ATMOS_PHY_AE_driver_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_AE_driver_finalize
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_finalize
    use scale_atmos_phy_ae_sprintars, only: &
       ATMOS_PHY_AE_SPRINTARS_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_ae ) then
       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_finalize
       case ( 'SPRINTARS' )
          call ATMOS_PHY_AE_SPRINTARS_finalize
       case ( 'OFFLINE' )
       end select
    endif

    return
  end subroutine ATMOS_PHY_AE_driver_finalize

  !-----------------------------------------------------------------------------
  !> adjustment
  subroutine ATMOS_PHY_AE_driver_adjustment
    use mod_atmos_vars, only: &
       QTRC
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE
    use mod_atmos_admin, only: &
       ATMOS_PHY_AE_TYPE, &
       ATMOS_sw_phy_ae
    use scale_atmos_phy_ae_kajino13, only: &
        ATMOS_PHY_AE_kajino13_negative_fixer
    use scale_atmos_phy_ae_SPRINTARS, only: &
        ATMOS_PHY_AE_SPRINTARS_negative_fixer
    implicit none

    if ( ATMOS_sw_phy_ae ) then

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_negative_fixer( KA, KS, KE, IA, IS, IE, JA, JS, JE, QA_AE, &
                                                     QTRC(:,:,:,QS_AE:QE_AE) ) ! [INOUT]

       case ( 'SPRINTARS' )
          call ATMOS_PHY_AE_SPRINTARS_negative_fixer( QA_AE, &                  ! [IN]
                                                      QTRC(:,:,:,QS_AE:QE_AE) ) ! [INOUT]
       end select

    end if

    return
  end subroutine ATMOS_PHY_AE_driver_adjustment

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_AE_driver_calc_tendency( update_flag )
    use scale_tracer, only: &
       TRACER_NAME
    use scale_prc, only: &
       PRC_abort
    use scale_time, only: &
       dt_AE => TIME_DTSEC_ATMOS_PHY_AE, &
       TIME_NOWDAYSEC
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_atmos_grid_cartesC_real, only: &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_file_history, only: &
       FILE_HISTORY_in
    use mod_atmos_vars, only: &
       DENS => DENS_av, &
       QTRC => QTRC_av, &
       QDRY, &
       PRES, &
       TEMP, &
       W,    &
       QV,   &
       QC,   &
       QI,   &
       RHOQ_t => RHOQ_tp, &
       THETA => POTT
    use mod_atmos_phy_sf_vars, only: &
         SFC_PRES => ATMOS_PHY_SF_SFC_PRES, &
         T2       => ATMOS_PHY_SF_T2,       &
         U10      => ATMOS_PHY_SF_U10,      &
         V10      => ATMOS_PHY_SF_V10,      &
         USTAR    => ATMOS_PHY_SF_Ustar,    &
         WSTAR    => ATMOS_PHY_SF_Wstar,    &
         QSTAR    => ATMOS_PHY_SF_Qstar,    &
         SFC_TEMP => ATMOS_PHY_SF_SFC_TEMP
    use mod_land_vars, only: &
         U10L => LAND_U10, &
         V10L => LAND_V10, &
         LAND_WATER,       &
         LAND_PROPERTY,    &
         I_WaterLimit
    use mod_ocean_vars, only: &
         OCEAN_ICE_FRAC
    use mod_atmos_phy_ae_vars, only: &
       QA_AE, &
       QS_AE, &
       QE_AE, &
       RHOQ_t_AE => ATMOS_PHY_AE_RHOQ_t, &
       CCN       => ATMOS_PHY_AE_CCN,    &
       CCN_t     => ATMOS_PHY_AE_CCN_t,  &
       AE_EMIT   => ATMOS_PHY_AE_EMIT,   &
       ATMOS_PHY_AE_Re, &
       ATMOS_PHY_AE_Qe
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic,                   &
       QA_MP,                                              &
       QS_MP,                                              &
       QE_MP,                                              &
       EVAPORATE           => ATMOS_PHY_MP_EVAPORATE,      &
       flux_qr_mp          => ATMOS_PHY_MP_mflux_qr_total, & ! rain flux of microphysics scheme for SPRINTARS
       flux_qs_mp          => ATMOS_PHY_MP_mflux_qs_total, & ! snow flux of microphysics scheme for SPRINTARS
       VTERM_MP_RAIN       => ATMOS_PHY_MP_vterm_rain,     & ! terminal velocity of rain [m/s]  for SPRINTARS
       VTERM_MP_SNOW       => ATMOS_PHY_MP_vterm_snow,     & ! terminal velocity of snow [m/s]  for SPRINTARS
       RHOQ_t_MP_SPRINTARS => ATMOS_PHY_MP_RHOQ_t_SPRINTARS  ! tendency rho*QTRC         [kg/m3/s]  for SPRINTARS
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HR, &
       I_HS
    use mod_atmos_admin, only: &
       ATMOS_PHY_CP_TYPE, &
       ATMOS_PHY_MP_TYPE, &
       ATMOS_PHY_AE_TYPE
    use scale_atmos_phy_ae_kajino13, only: &
       ATMOS_PHY_AE_kajino13_tendency
    use scale_atmos_phy_ae_offline, only: &
       ATMOS_PHY_AE_offline_tendency
    use scale_atmos_phy_ae_SPRINTARS, only: &
       ATMOS_PHY_AE_SPRINTARS_tendency
    use mod_atmos_phy_cp_vars, only: &
       cldfrac_cp1 => ATMOS_PHY_CP_cldfrac_dp,     & ! cloud fraction (deep convection) (0-1)
       flux_qr_cp1 => ATMOS_PHY_CP_flux_qr_total,  & ! rain flux of cumulus parameterization for SPRINTARS
       flux_qs_cp1 => ATMOS_PHY_CP_flux_qs_total,  & ! snow flux of cumulus parameterization for SPRINTARS
       RHOHYD_t_CP => ATMOS_PHY_CP_RHOHYD_t          ! tendency rho*QHYD [kg/m3/s]
    use mod_atmos_phy_rd_vars, only: &
       dtau_s_sprintars => ATMOS_PHY_RD_dtau_s,                   & ! 0.67 micron cloud optical depth
       SFCFLX_SW_dn_net => ATMOS_PHY_RD_SFLX_SW_dn_net_sprintars, & ! net downward SW flux for SPRINTARS
       RCOSZ            => ATMOS_PHY_RD_cosSZA                      ! cos(solar zenith angle) (0-1)
    use mod_atmos_phy_bl_vars, only: &
       QS_BL => QS

    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: CN    (KA,IA,JA)
    real(RP) :: CN_CC (KA,IA,JA)
    real(RP) :: CCN_CC(KA,IA,JA)
    real(RP) :: NREG  (KA,IA,JA)

    real(RP) :: cldfrac_cp(KA,  IA,JA) ! cloud fraction (deep convection) (0-1)
    real(RP) :: flux_qr_cp(0:KA,IA,JA) ! rain flux of cumulus parameterization for SPRINTARS
    real(RP) :: flux_qs_cp(0:KA,IA,JA) ! snow flux of cumulus parameterization for SPRINTARS
    real(RP) :: cldfrac_mp(KA,IA,JA)   ! cloud fraction of microphysics scheme [0-1]
    real(RP) :: Re(KA,IA,JA,N_HYD)     ! effective radius [cm]
    real(RP) :: RAIN_MP_Re(KA,IA,JA)   ! effective radius [m]
    real(RP) :: SNOW_MP_Re(KA,IA,JA)   ! effective radius [m]
    real(RP) :: RHOQ_t_MP_RAIN(KA,IA,JA) ! tendency rho*QR (MP) [kg/m3/s]
    real(RP) :: RHOQ_t_MP_SNOW(KA,IA,JA) ! tendency rho*QS (MP) [kg/m3/s]
    real(RP) :: RHOQ_t_CP_PREC(KA,IA,JA) ! tendency rho*QS (CP) [kg/m3/s]
    real(RP) :: NUMDEN_MP_CLDW(KA,IA,JA) ! number density (MP;cloud water) [num/m3]
    real(RP) :: NUMDEN_MP_RAIN(KA,IA,JA) ! number density (MP;rain       ) [num/m3]
    real(RP) :: NUMDEN_MP_CLDI(KA,IA,JA) ! number density (MP;cloud ice  ) [num/m3]
    real(RP) :: NUMDEN_MP_SNOW(KA,IA,JA) ! number density (MP;snow       ) [num/m3]

    integer  :: k, i, j, iq
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       !$acc data create(CN, CN_CC, CCN_CC, NREG)

       !$acc kernels
!OCL XFILL
       CCN      (:,:,:)   = 0.0_RP ! reset
!OCL XFILL
       RHOQ_t_AE(:,:,:,:) = 0.0_RP ! reset
       !$acc end kernels

       !$acc kernels
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          NREG(k,i,j) = EVAPORATE(k,i,j) * dt_AE
       enddo
       enddo
       enddo
       !$acc end kernels

       select case ( ATMOS_PHY_AE_TYPE )
       case ( 'KAJINO13' )
          call ATMOS_PHY_AE_kajino13_tendency( KA, KS, KE, IA, IS, IE, JA, JS, JE, QA_AE, &
                                               TEMP     (:,:,:),             & ! [IN]
                                               PRES     (:,:,:),             & ! [IN]
                                               QDRY     (:,:,:),             & ! [IN]
                                               NREG     (:,:,:),             & ! [IN]
                                               DENS     (:,:,:),             & ! [IN]
                                               QV       (:,:,:),             & ! [IN]
                                               QTRC     (:,:,:,QS_AE:QE_AE), & ! [IN]
                                               AE_EMIT  (:,:,:,QS_AE:QE_AE), & ! [IN]
                                               dt_AE,                        & ! [IN]
                                               RHOQ_t_AE(:,:,:,QS_AE:QE_AE), & ! [OUT]
                                               CN       (:,:,:),             & ! [OUT]
                                               CCN      (:,:,:)              ) ! [OUT]
       case ( 'OFFLINE' )
          call ATMOS_PHY_AE_offline_tendency ( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                               TIME_NOWDAYSEC,               & ! [IN]
                                               CCN      (:,:,:)              ) ! [OUT]
          !$acc kernels
          CN(:,:,:) = 0.0_RP ! not supported
          !$acc end kernels


       case ( 'SPRINTARS' )

          ! cloud fraction (CP)
          if ( ATMOS_PHY_CP_TYPE == 'KF' ) then
             !$acc update host(cldfrac_cp1,flux_qr_cp1,flux_qs_cp1,RHOHYD_t_CP)
             cldfrac_cp(1:KA,:,:) = cldfrac_cp1(1:KA,:,:)
             flux_qr_cp(0:KA,:,:) = flux_qr_cp1(0:KA,:,:)
             flux_qs_cp(0:KA,:,:) = flux_qs_cp1(0:KA,:,:)
             do j = JS, JE
             do i = IS, IE
                do k = KS, KE
                   RHOQ_t_CP_PREC(k,i,j) = RHOHYD_t_CP(k,i,j,I_HR) + RHOHYD_t_CP(k,i,j,I_HS)
                enddo
             enddo
             enddo
          else
             cldfrac_cp(1:KA,:,:) = 0.0_RP
             flux_qr_cp(0:KA,:,:) = 0.0_RP
             flux_qs_cp(0:KA,:,:) = 0.0_RP
             do j = JS, JE
             do i = IS, IE
                do k = KS, KE
                   RHOQ_t_CP_PREC(k,i,j) = 0.0_RP
                enddo
             enddo
             enddo
          endif

          ! cloud fraction (MP)
          call ATMOS_PHY_MP_vars_get_diagnostic( DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
                                                 CLDFRAC = cldfrac_mp,                    & ! [OUT]
                                                 Re = Re                                  ) ! [OUT]

          ! ! cloud fraction (Goto et al., 2020, NICAM)
          ! do j = JS, JE
          ! do i = IS, IE
          !    do k = KS, KE
          !       if ( (QC(k,i,j)+QI(k,i,j))*DENS(k,i,j) >= 1.0E-3_RP ) then
          !          cldfrac_mp(k,i,j) = 1.0_RP
          !       else
          !          cldfrac_mp(k,i,j) = 0.0_RP
          !       endif
          !    enddo
          ! enddo
          ! enddo
          ! !effective radius of hedrometors Re [cm]
          ! call ATMOS_PHY_MP_vars_get_diagnostic( DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
          !                                        Re = Re                                  ) ! [OUT]

          !$acc update host(Re,RHOQ_t_MP_SPRINTARS)
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                RAIN_MP_Re    (k,i,j) = Re(k,i,j,I_HR) * 0.01_RP ! [cm] -> [m]
                SNOW_MP_Re    (k,i,j) = Re(k,i,j,I_HS) * 0.01_RP ! [cm] -> [m]
                RHOQ_t_MP_RAIN(k,i,j) = RHOQ_t_MP_SPRINTARS(k,i,j,QS_MP+I_HR)
                RHOQ_t_MP_SNOW(k,i,j) = RHOQ_t_MP_SPRINTARS(k,i,j,QS_MP+I_HS)
             enddo
          enddo
          enddo

          !Hydrometeor Number Density
          if ( ATMOS_PHY_MP_TYPE == 'SN14' ) then
             !$acc update host(QTRC,DENS)
             do j = JS, JE
             do i = IS, IE
                do k = KS, KE
                   NUMDEN_MP_CLDW(k,i,j) = QTRC(k,i,j,QS_MP-1+ 7) * DENS(k,i,j) ![num/kg] -> [num/m3]
                   NUMDEN_MP_RAIN(k,i,j) = QTRC(k,i,j,QS_MP-1+ 8) * DENS(k,i,j) ![num/kg] -> [num/m3]
                   NUMDEN_MP_CLDI(k,i,j) = QTRC(k,i,j,QS_MP-1+ 9) * DENS(k,i,j) ![num/kg] -> [num/m3]
                   NUMDEN_MP_SNOW(k,i,j) = QTRC(k,i,j,QS_MP-1+10) * DENS(k,i,j) ![num/kg] -> [num/m3]
                enddo
             enddo
             enddo
          else
             NUMDEN_MP_CLDW(:,:,:) = 0.0_RP ![num/m3]
             NUMDEN_MP_RAIN(:,:,:) = 0.0_RP ![num/m3]
             NUMDEN_MP_CLDI(:,:,:) = 0.0_RP ![num/m3]
             NUMDEN_MP_SNOW(:,:,:) = 0.0_RP ![num/m3]
          endif

          !$acc update host(PRES,TEMP,DENS,W,QTRC,QDRY,QV,QC,QI,THETA,flux_qr_mp,flux_qs_mp,VTERM_MP_RAIN,VTERM_MP_SNOW,dtau_s_sprintars)
          !$acc update host(SFC_PRES,T2,U10,V10,U10L,V10L,USTAR,WSTAR,QSTAR,SFC_TEMP,LAND_WATER,LAND_PROPERTY,OCEAN_ICE_FRAC)
          !$acc update host(SFCFLX_SW_dn_net,RCOSZ)
          call ATMOS_PHY_AE_SPRINTARS_tendency( QA_AE,                           & ! [IN]
                                                PRES      (:,:,:),               & ! [IN] 3D
                                                TEMP      (:,:,:),               & ! [IN]
                                                DENS      (:,:,:),               & ! [IN]
                                                W         (:,:,:),               & ! [IN]
                                                QTRC      (:,:,:,QS_BL),         & ! [IN]
                                                QDRY      (:,:,:),               & ! [IN]
                                                QV        (:,:,:),               & ! [IN]
                                                QC        (:,:,:),               & ! [IN]
                                                QI        (:,:,:),               & ! [IN]
                                                THETA     (:,:,:),               & ! [IN]
                                                RAIN_MP_Re(:,:,:),               & ! [IN]
                                                SNOW_MP_Re(:,:,:),               & ! [IN]
                                                VTERM_MP_RAIN(:,:,:),            & ! [IN]
                                                VTERM_MP_SNOW(:,:,:),            & ! [IN]
                                                RHOQ_t_MP_RAIN(:,:,:),           & ! [IN]
                                                RHOQ_t_MP_SNOW(:,:,:),           & ! [IN]
                                                RHOQ_t_CP_PREC(:,:,:),           & ! [IN]
                                                flux_qr_cp(0:KA,:,:),            & ! [IN]
                                                flux_qs_cp(0:KA,:,:),            & ! [IN]
                                                cldfrac_cp(:,:,:),               & ! [IN]
                                                flux_qr_mp(0:KA,:,:),            & ! [IN]
                                                flux_qs_mp(0:KA,:,:),            & ! [IN]
                                                cldfrac_mp(:,:,:),               & ! [IN]
                                                ATMOS_PHY_MP_TYPE,               & ! [IN]
                                                NUMDEN_MP_CLDW(:,:,:),           & ! [IN]
                                                NUMDEN_MP_RAIN(:,:,:),           & ! [IN]
                                                NUMDEN_MP_CLDI(:,:,:),           & ! [IN]
                                                NUMDEN_MP_SNOW(:,:,:),           & ! [IN]
                                                dtau_s_sprintars(:,:,:),         & ! [IN]
                                                SFC_PRES(:,:),                   & ! [IN] 2D
                                                T2      (:,:),                   & ! [IN]
                                                U10     (:,:),                   & ! [IN]
                                                V10     (:,:),                   & ! [IN]
                                                U10L    (:,:),                   & ! [IN]
                                                V10L    (:,:),                   & ! [IN]
                                                USTAR   (:,:),                   & ! [IN]
                                                WSTAR   (:,:),                   & ! [IN]
                                                QSTAR   (:,:),                   & ! [IN]
                                                SFC_TEMP(:,:),                   & ! [IN]
                                                LAND_WATER (1,:,:),              & ! [IN]
                                                LAND_PROPERTY(:,:,I_WaterLimit), & ! [IN]
                                                OCEAN_ICE_FRAC(:,:),             & ! [IN]
                                                SFCFLX_SW_dn_net(:,:),           & ! [IN]
                                                RCOSZ           (:,:),           & ! [IN]
                                                dt_AE,                           & ! [IN]
                                                QTRC     (:,:,:,QS_AE:QE_AE),    & ! [IN]
                                                RHOQ_t_AE(:,:,:,QS_AE:QE_AE),    & ! [INOUT]
                                                CCN(:,:,:),                      & ! [OUT]
                                                CN (:,:,:),                      & ! [OUT]
                                                ATMOS_PHY_AE_Re(:,:,:,:),        & ! [OUT]
                                                ATMOS_PHY_AE_Qe(:,:,:,:)         ) ! [OUT]
          !$acc update device(QTRC,RHOQ_t_AE,CCN,CN,ATMOS_PHY_AE_Re,ATMOS_PHY_AE_Qe)

       end select

       !$acc kernels
       CCN_t(:,:,:) = CCN(:,:,:) / dt_AE
       !$acc end kernels

       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             CN_CC (k,i,j) = CN (k,i,j) * 1.0E-6_RP
             CCN_CC(k,i,j) = CCN(k,i,j) * 1.0E-6_RP
          enddo
       enddo
       enddo
       !$acc end kernels
       call FILE_HISTORY_in( CN_CC (:,:,:), 'CN',  'condensation nucrei',       'num/cc', dim_type = 'ZXY' )
!       call FILE_HISTORY_in( cldfrac_mp(:,:,:), 'CN',  'condensation nucrei',       'num/cc', dim_type = 'ZXY' )
       call FILE_HISTORY_in( CCN_CC(:,:,:), 'CCN', 'cloud condensation nucrei', 'num/cc', dim_type = 'ZXY' )

!       call FILE_HISTORY_in( CN (:,:,:)*1.E-6_RP, 'CN',  'condensation nucrei',       'num/cc' )
!       call FILE_HISTORY_in( CCN(:,:,:)*1.E-6_RP, 'CCN', 'cloud condensation nucrei', 'num/cc' )

       !$acc end data

    endif

    !$omp parallel do private(i,j,k,iq) OMP_SCHEDULE_
    !$acc kernels
    do iq = QS_AE, QE_AE
       do j  = JS, JE
       do i  = IS, IE
       do k  = KS, KE
          RHOQ_t(k,i,j,iq) = RHOQ_t(k,i,j,iq) + RHOQ_t_AE(k,i,j,iq)
       enddo
       enddo
       enddo
    enddo
    !$acc end kernels

    if ( STATISTICS_checktotal ) then
       do iq = QS_AE, QE_AE
          call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                                 RHOQ_t_AE(:,:,:,iq), trim(TRACER_NAME(iq))//'_t_AE', &
                                 ATMOS_GRID_CARTESC_REAL_VOL(:,:,:),                  &
                                 ATMOS_GRID_CARTESC_REAL_TOTVOL                       )
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_driver_calc_tendency

end module mod_atmos_phy_ae_driver
