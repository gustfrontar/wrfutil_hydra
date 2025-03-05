!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiation
!!
!! @par Description
!!          Atmospheric radiation transfer process driver
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_atmos_phy_rd_driver
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_atmos_grid_cartesC_index
  use scale_tracer
  use scale_cpl_sfc_index
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_RD_driver_setup
  public :: ATMOS_PHY_RD_driver_finalize
  public :: ATMOS_PHY_RD_driver_calc_tendency

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
  logical, private :: RD_use_PBL_cloud = .false.

  integer, private :: HIST_id_SOLINS
  integer, private :: HIST_id_COSZ
  integer, private :: HIST_id_SFLX_LW_up_c
  integer, private :: HIST_id_SFLX_LW_dn_c
  integer, private :: HIST_id_SFLX_SW_up_c
  integer, private :: HIST_id_SFLX_SW_dn_c
  integer, private :: HIST_id_SFLX_LW_up
  integer, private :: HIST_id_SFLX_LW_dn
  integer, private :: HIST_id_SFLX_SW_up
  integer, private :: HIST_id_SFLX_SW_dn
  integer, private :: HIST_id_SFLX_LW_net
  integer, private :: HIST_id_SFLX_SW_net
  integer, private :: HIST_id_TOAFLX_LW_up_c
  integer, private :: HIST_id_TOAFLX_LW_dn_c
  integer, private :: HIST_id_TOAFLX_SW_up_c
  integer, private :: HIST_id_TOAFLX_SW_dn_c
  integer, private :: HIST_id_TOAFLX_LW_up
  integer, private :: HIST_id_TOAFLX_LW_dn
  integer, private :: HIST_id_TOAFLX_SW_up
  integer, private :: HIST_id_TOAFLX_SW_dn
  integer, private :: HIST_id_TOAFLX_LW_net
  integer, private :: HIST_id_TOAFLX_SW_net
  integer, private :: HIST_id_TOMFLX_LW_up_c
  integer, private :: HIST_id_TOMFLX_LW_dn_c
  integer, private :: HIST_id_TOMFLX_SW_up_c
  integer, private :: HIST_id_TOMFLX_SW_dn_c
  integer, private :: HIST_id_TOMFLX_LW_up
  integer, private :: HIST_id_TOMFLX_LW_dn
  integer, private :: HIST_id_TOMFLX_SW_up
  integer, private :: HIST_id_TOMFLX_SW_dn
  integer, private :: HIST_id_TOMFLX_LW_net
  integer, private :: HIST_id_TOMFLX_SW_net
  integer, private :: HIST_id_SLR
  integer, private :: HIST_id_SSR
  integer, private :: HIST_id_OLR
  integer, private :: HIST_id_OSR
  integer, private :: HIST_id_RADFLUX_LWUP
  integer, private :: HIST_id_RADFLUX_LWDN
  integer, private :: HIST_id_RADFLUX_LW
  integer, private :: HIST_id_RADFLUX_SWUP
  integer, private :: HIST_id_RADFLUX_SWDN
  integer, private :: HIST_id_RADFLUX_SW
  integer, private :: HIST_id_SFLX_IR_dn_dir
  integer, private :: HIST_id_SFLX_IR_dn_dif
  integer, private :: HIST_id_SFLX_NIR_dn_dir
  integer, private :: HIST_id_SFLX_NIR_dn_dif
  integer, private :: HIST_id_SFLX_VIS_dn_dir
  integer, private :: HIST_id_SFLX_VIS_dn_dif
  integer, private :: HIST_id_TEMP_t_rd_LW
  integer, private :: HIST_id_TEMP_t_rd_SW
  integer, private :: HIST_id_TEMP_t_rd
  integer, private :: HIST_id_RHOH_RD
  integer, private :: HIST_id_RFLX_LW_up
  integer, private :: HIST_id_RFLX_LW_dn
  integer, private :: HIST_id_RFLX_SW_up
  integer, private :: HIST_id_RFLX_SW_dn
  integer, private :: HIST_id_dtau_s
  integer, private :: HIST_id_dem_s

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_RD_driver_setup
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_mstrnx_setup
    use scale_atmos_phy_rd_offline, only: &
       ATMOS_PHY_RD_offline_setup
    use scale_atmos_grid_cartesC, only: &
       CZ => ATMOS_GRID_CARTESC_CZ, &
       FZ => ATMOS_GRID_CARTESC_FZ
    use scale_file_history, only: &
       FILE_HISTORY_reg
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_sw_phy_rd
    use mod_atmos_phy_rd_vars, only: &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       TOMFLX_LW_up => ATMOS_PHY_RD_TOMFLX_LW_up, &
       TOMFLX_LW_dn => ATMOS_PHY_RD_TOMFLX_LW_dn, &
       TOMFLX_SW_up => ATMOS_PHY_RD_TOMFLX_SW_up, &
       TOMFLX_SW_dn => ATMOS_PHY_RD_TOMFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    implicit none
    
    namelist / PARAM_ATMOS_PHY_RD / &
         RD_use_PBL_cloud

    integer :: ierr
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'Setup'

    if ( ATMOS_sw_phy_rd ) then

       !--- read namelist
       rewind(IO_FID_CONF)
       read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_RD,iostat=ierr)
       if( ierr < 0 ) then !--- missing
          LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'Not found namelist. Default used.'
       elseif( ierr > 0 ) then !--- fatal error
          LOG_ERROR("ATMOS_PHY_RD_driver_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_RD. Check!'
          call PRC_abort
       endif
       LOG_NML(PARAM_ATMOS_PHY_RD)

       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )
          call ATMOS_PHY_RD_MSTRNX_setup( KA, KS, KE, CZ(:), FZ(:) )
       case ( "OFFLINE" )
          call ATMOS_PHY_RD_offline_setup
       case default
          LOG_ERROR("ATMOS_PHY_RD_driver_setup",*) 'invalid Radiation type(', trim(ATMOS_PHY_RD_TYPE), '). CHECK!'
          call PRC_abort
       end select

    else

       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'this component is never called.'
       LOG_INFO("ATMOS_PHY_RD_driver_setup",*) 'radiation fluxes are set to zero.'
       !$acc kernels
       SFCFLX_LW_up(:,:)     = 0.0_RP
       SFCFLX_LW_dn(:,:)     = 0.0_RP
       SFCFLX_SW_up(:,:)     = 0.0_RP
       SFCFLX_SW_dn(:,:)     = 0.0_RP
       TOMFLX_LW_up(:,:)     = 0.0_RP
       TOMFLX_LW_dn(:,:)     = 0.0_RP
       TOMFLX_SW_up(:,:)     = 0.0_RP
       TOMFLX_SW_dn(:,:)     = 0.0_RP
       SFLX_rad_dn (:,:,:,:) = 0.0_RP
       solins      (:,:)     = 0.0_RP
       cosSZA      (:,:)     = 0.0_RP
       !$acc end kernels

    endif


    ! history
    call FILE_HISTORY_reg( 'SOLINS', 'solar insolation',        'W/m2', HIST_id_SOLINS, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'COSZ',   'cos(solar zenith angle)', '1',    HIST_id_COSZ,   ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'SFLX_LW_up_c',   'SFC upward   longwave  radiation flux (clr)', 'W/m2', HIST_id_SFLX_LW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_LW_dn_c',   'SFC downward longwave  radiation flux (clr)', 'W/m2', HIST_id_SFLX_LW_dn_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_SW_up_c',   'SFC upward   shortwave radiation flux (clr)', 'W/m2', HIST_id_SFLX_SW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_SW_dn_c',   'SFC downward shortwave radiation flux (clr)', 'W/m2', HIST_id_SFLX_SW_dn_c, ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'SFLX_LW_up',     'SFC upward   longwave  radiation flux',       'W/m2', HIST_id_SFLX_LW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_LW_dn',     'SFC downward longwave  radiation flux',       'W/m2', HIST_id_SFLX_LW_dn,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_SW_up',     'SFC upward   shortwave radiation flux',       'W/m2', HIST_id_SFLX_SW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_SW_dn',     'SFC downward shortwave radiation flux',       'W/m2', HIST_id_SFLX_SW_dn,   ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'SFLX_LW_net',    'SFC net      longwave  radiation flux',       'W/m2', HIST_id_SFLX_LW_net,  ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_SW_net',    'SFC net      shortwave radiation flux',       'W/m2', HIST_id_SFLX_SW_net,  ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOAFLX_LW_up_c', 'TOA upward   longwave  radiation flux (clr)', 'W/m2', HIST_id_TOAFLX_LW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_LW_dn_c', 'TOA downward longwave  radiation flux (clr)', 'W/m2', HIST_id_TOAFLX_LW_dn_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_SW_up_c', 'TOA upward   shortwave radiation flux (clr)', 'W/m2', HIST_id_TOAFLX_SW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_SW_dn_c', 'TOA downward shortwave radiation flux (clr)', 'W/m2', HIST_id_TOAFLX_SW_dn_c, ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOAFLX_LW_up',   'TOA upward   longwave  radiation flux',       'W/m2', HIST_id_TOAFLX_LW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_LW_dn',   'TOA downward longwave  radiation flux',       'W/m2', HIST_id_TOAFLX_LW_dn,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_SW_up',   'TOA upward   shortwave radiation flux',       'W/m2', HIST_id_TOAFLX_SW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_SW_dn',   'TOA downward shortwave radiation flux',       'W/m2', HIST_id_TOAFLX_SW_dn,   ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOAFLX_LW_net',  'TOA net      longwave  radiation flux',       'W/m2', HIST_id_TOAFLX_LW_net,  ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOAFLX_SW_net',  'TOA net      shortwave radiation flux',       'W/m2', HIST_id_TOAFLX_SW_net,  ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOMFLX_LW_up_c', 'TOM upward   longwave  radiation flux (clr)', 'W/m2', HIST_id_TOMFLX_LW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_LW_dn_c', 'TOM downward longwave  radiation flux (clr)', 'W/m2', HIST_id_TOMFLX_LW_dn_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_SW_up_c', 'TOM upward   shortwave radiation flux (clr)', 'W/m2', HIST_id_TOMFLX_SW_up_c, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_SW_dn_c', 'TOM downward shortwave radiation flux (clr)', 'W/m2', HIST_id_TOMFLX_SW_dn_c, ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOMFLX_LW_up',   'TOM upward   longwave  radiation flux',       'W/m2', HIST_id_TOMFLX_LW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_LW_dn',   'TOM downward longwave  radiation flux',       'W/m2', HIST_id_TOMFLX_LW_dn,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_SW_up',   'TOM upward   shortwave radiation flux',       'W/m2', HIST_id_TOMFLX_SW_up,   ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_SW_dn',   'TOM downward shortwave radiation flux',       'W/m2', HIST_id_TOMFLX_SW_dn,   ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TOMFLX_LW_net',  'TOM net      longwave  radiation flux',       'W/m2', HIST_id_TOMFLX_LW_net,  ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TOMFLX_SW_net',  'TOM net      shortwave radiation flux',       'W/m2', HIST_id_TOMFLX_SW_net,  ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'SLR',          'SFC net longwave  radiation flux',  'W/m2', HIST_id_SLR, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SSR',          'SFC net shortwave radiation flux',  'W/m2', HIST_id_SSR, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'OLR',          'outgoing longwave  radiation flux', 'W/m2', HIST_id_OLR, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'OSR',          'outgoing shortwave radiation flux', 'W/m2', HIST_id_OSR, ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'RADFLUX_LWUP', 'upward   longwave  radiation flux', 'W/m2', HIST_id_RADFLUX_LWUP, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RADFLUX_LWDN', 'downward longwave  radiation flux', 'W/m2', HIST_id_RADFLUX_LWDN, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RADFLUX_LW',   'net      longwave  radiation flux', 'W/m2', HIST_id_RADFLUX_LW,   fill_halo=.true. )
    call FILE_HISTORY_reg( 'RADFLUX_SWUP', 'upward   shortwave radiation flux', 'W/m2', HIST_id_RADFLUX_SWUP, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RADFLUX_SWDN', 'downward shortwave radiation flux', 'W/m2', HIST_id_RADFLUX_SWDN, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RADFLUX_SW',   'net      shortwave radiation flux', 'W/m2', HIST_id_RADFLUX_SW,   fill_halo=.true. )

    call FILE_HISTORY_reg( 'SFLX_IR_dn_dir',  'SFC downward radiation flux (direct; IR)',   'W/m2', HIST_id_SFLX_IR_dn_dir,  ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_IR_dn_dif',  'SFC downward radiation flux (diffuse; IR)',  'W/m2', HIST_id_SFLX_IR_dn_dif,  ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_NIR_dn_dir', 'SFC downward radiation flux (direct; NIR)',  'W/m2', HIST_id_SFLX_NIR_dn_dir, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_NIR_dn_dif', 'SFC downward radiation flux (diffuse; NIR)', 'W/m2', HIST_id_SFLX_NIR_dn_dif, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_VIS_dn_dir', 'SFC downward radiation flux (direct; VIS)',  'W/m2', HIST_id_SFLX_VIS_dn_dir, ndims=2, fill_halo=.true. )
    call FILE_HISTORY_reg( 'SFLX_VIS_dn_dif', 'SFC downward radiation flux (diffuse; VIS)', 'W/m2', HIST_id_SFLX_VIS_dn_dif, ndims=2, fill_halo=.true. )

    call FILE_HISTORY_reg( 'TEMP_t_rd_LW', 'tendency of temp in rd(LW)',  'K/day',  HIST_id_TEMP_t_rd_LW, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TEMP_t_rd_SW', 'tendency of temp in rd(SW)',  'K/day',  HIST_id_TEMP_t_rd_SW, fill_halo=.true. )
    call FILE_HISTORY_reg( 'TEMP_t_rd',    'tendency of temp in rd',      'K/day',  HIST_id_TEMP_t_rd,    fill_halo=.true. )
    call FILE_HISTORY_reg( 'RHOH_RD',      'diabatic heating rate in rd', 'J/m3/s', HIST_id_RHOH_Rd,      fill_halo=.true. )

    ! output of raw data, for offline output
    call FILE_HISTORY_reg( 'RFLX_LW_up', 'upward   longwave  radiation flux (cell face)', 'W/m2', HIST_id_RFLX_LW_up, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RFLX_LW_dn', 'downward longwave  radiation flux (cell face)', 'W/m2', HIST_id_RFLX_LW_dn, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RFLX_SW_up', 'upward   shortwave radiation flux (cell face)', 'W/m2', HIST_id_RFLX_SW_up, fill_halo=.true. )
    call FILE_HISTORY_reg( 'RFLX_SW_dn', 'downward shortwave radiation flux (cell face)', 'W/m2', HIST_id_RFLX_SW_dn, fill_halo=.true. )

    call FILE_HISTORY_reg( 'dtau_s', '0.67 micron cloud optical depth', '1', HIST_id_dtau_s, fill_halo=.true. )
    call FILE_HISTORY_reg( 'dem_s',  '10.5 micron cloud emissivity',    '1', HIST_id_dem_s,  fill_halo=.true. )

    return
  end subroutine ATMOS_PHY_RD_driver_setup

  !-----------------------------------------------------------------------------
  !> finalize
  subroutine ATMOS_PHY_RD_driver_finalize
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_MSTRNX_finalize
    use scale_atmos_phy_rd_profile, only: &
       ATMOS_PHY_RD_PROFILE_finalize
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE, &
       ATMOS_sw_phy_rd
    implicit none
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_RD_driver_finalize",*) 'Finalize'

    if ( ATMOS_sw_phy_rd ) then
       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )
          call ATMOS_PHY_RD_MSTRNX_finalize
       case ( "OFFLINE" )
       end select
    end if

    call ATMOS_PHY_RD_PROFILE_finalize

    return
  end subroutine ATMOS_PHY_RD_driver_finalize

  !-----------------------------------------------------------------------------
  !> Driver
  subroutine ATMOS_PHY_RD_driver_calc_tendency( update_flag )
    use scale_const, only: &
       EPS => CONST_EPS
    use scale_atmos_grid_cartesC_real, only: &
       REAL_CZ  => ATMOS_GRID_CARTESC_REAL_CZ,  &
       REAL_FZ  => ATMOS_GRID_CARTESC_REAL_FZ,  &
       REAL_LON => ATMOS_GRID_CARTESC_REAL_LON, &
       REAL_LAT => ATMOS_GRID_CARTESC_REAL_LAT, &
       ATMOS_GRID_CARTESC_REAL_VOL, &
       ATMOS_GRID_CARTESC_REAL_TOTVOL
    use scale_landuse, only: &
       fact_ocean => LANDUSE_fact_ocean, &
       fact_land  => LANDUSE_fact_land,  &
       fact_urban => LANDUSE_fact_urban
    use scale_time, only: &
       TIME_NOWDATE,                     &
       TIME_NOWDAYSEC,                   &
       TIME_OFFSET_YEAR
    use scale_statistics, only: &
       STATISTICS_checktotal, &
       STATISTICS_total
    use scale_file_history, only: &
       FILE_HISTORY_put, &
       FILE_HISTORY_query_f
    use scale_atmos_hydrometeor, only: &
       N_HYD, &
       I_HC, &
       I_HI
    use scale_atmos_aerosol, only: &
       N_AE
    use mod_atmos_admin, only: &
       ATMOS_PHY_RD_TYPE
    use scale_atmos_solarins, only: &
       SOLARINS_insolation => ATMOS_SOLARINS_insolation
    use scale_atmos_phy_rd_mstrnx, only: &
       ATMOS_PHY_RD_MSTRNX_flux
    use scale_atmos_phy_rd_offline, only: &
       ATMOS_PHY_RD_OFFLINE_flux
    use scale_atmos_phy_rd_common, only: &
       ATMOS_PHY_RD_calc_heating, &
       I_SW, &
       I_LW, &
       I_dn, &
       I_up, &
       I_Cloud, &
       I_ClearSky
    use mod_atmos_vars, only: &
       TEMP,              &
       PRES,              &
       QV,                &
       CVtot,             &
       DENS   => DENS_av, &
       QTRC   => QTRC_av, &
       RHOH   => RHOH_p
    use mod_atmos_phy_sf_vars, only: &
       SFC_TEMP    => ATMOS_PHY_SF_SFC_TEMP,   &
       SFC_albedo  => ATMOS_PHY_SF_SFC_albedo
    use mod_atmos_phy_rd_vars, only: &
       RHOH_RD      => ATMOS_PHY_RD_RHOH,         &
       SFCFLX_LW_up => ATMOS_PHY_RD_SFLX_LW_up,   &
       SFCFLX_LW_dn => ATMOS_PHY_RD_SFLX_LW_dn,   &
       SFCFLX_SW_up => ATMOS_PHY_RD_SFLX_SW_up,   &
       SFCFLX_SW_dn => ATMOS_PHY_RD_SFLX_SW_dn,   &
       SFCFLX_SW_dn_net_sprintars => ATMOS_PHY_RD_SFLX_SW_dn_net_sprintars, & ! ### for SPRINTARS
       dtau_s_sprintars => ATMOS_PHY_RD_dtau_s,   & ! ### for SPRINTARS
       TOMFLX_LW_up => ATMOS_PHY_RD_TOMFLX_LW_up, &
       TOMFLX_LW_dn => ATMOS_PHY_RD_TOMFLX_LW_dn, &
       TOMFLX_SW_up => ATMOS_PHY_RD_TOMFLX_SW_up, &
       TOMFLX_SW_dn => ATMOS_PHY_RD_TOMFLX_SW_dn, &
       SFLX_rad_dn  => ATMOS_PHY_RD_SFLX_down,    &
       solins       => ATMOS_PHY_RD_solins,       &
       cosSZA       => ATMOS_PHY_RD_cosSZA
    use mod_atmos_phy_bl_vars, only: &
       ATMOS_PHY_BL_Zi,     &
       ATMOS_PHY_BL_QL,     &
       ATMOS_PHY_BL_cldfrac
    use mod_atmos_vars, only: &
       ATMOS_vars_get_diagnostic
    use mod_atmos_phy_mp_vars, only: &
       ATMOS_PHY_MP_vars_get_diagnostic
    use mod_atmos_phy_ae_vars, only: &
       ATMOS_PHY_AE_vars_get_diagnostic
    implicit none

    logical, intent(in) :: update_flag

    real(RP) :: TEMP_t      (KA,IA,JA,3)
    real(RP) :: flux_rad    (KA,IA,JA,2,2,2)
    real(RP) :: flux_rad_top(   IA,JA,2,2,2)

    real(RP) :: dtau_s      (KA,IA,JA) ! 0.67 micron cloud optical depth
    real(RP) :: dem_s       (KA,IA,JA) ! 10.5 micron cloud emissivity

    real(RP) :: flux_up     (KA,IA,JA,2)
    real(RP) :: flux_dn     (KA,IA,JA,2)
    real(RP) :: flux_net    (KA,IA,JA,2)
    real(RP) :: flux_net_sfc(   IA,JA,2)
    real(RP) :: flux_net_toa(   IA,JA,2)
    real(RP) :: flux_net_tom(   IA,JA,2)

    real(RP) :: TOAFLX_LW_up(IA,JA)
    real(RP) :: TOAFLX_LW_dn(IA,JA)
    real(RP) :: TOAFLX_SW_up(IA,JA)
    real(RP) :: TOAFLX_SW_dn(IA,JA)

    real(RP) :: SFCFLX_LW_up_c(IA,JA)
    real(RP) :: SFCFLX_LW_dn_c(IA,JA)
    real(RP) :: SFCFLX_SW_up_c(IA,JA)
    real(RP) :: SFCFLX_SW_dn_c(IA,JA)
    real(RP) :: TOAFLX_LW_up_c(IA,JA)
    real(RP) :: TOAFLX_LW_dn_c(IA,JA)
    real(RP) :: TOAFLX_SW_up_c(IA,JA)
    real(RP) :: TOAFLX_SW_dn_c(IA,JA)
    real(RP) :: TOMFLX_LW_up_c(IA,JA)
    real(RP) :: TOMFLX_LW_dn_c(IA,JA)
    real(RP) :: TOMFLX_SW_up_c(IA,JA)
    real(RP) :: TOMFLX_SW_dn_c(IA,JA)

    real(RP) :: CLDFRAC(KA,IA,JA)
    real(RP) :: MP_Re  (KA,IA,JA,N_HYD)
    real(RP) :: MP_Qe  (KA,IA,JA,N_HYD)
    real(RP) :: AE_Re  (KA,IA,JA,N_AE)
    real(RP) :: AE_Qe  (KA,IA,JA,N_AE)

    real(RP) :: RH(KA,IA,JA)

    logical :: clear_sky

    integer  :: k, i, j
    !---------------------------------------------------------------------------

    if ( update_flag ) then

       !$acc data create(TEMP_t,flux_rad,flux_rad_top,dtau_s,dem_s,flux_up,flux_dn,flux_net,flux_net_sfc,flux_net_toa,flux_net_tom,TOAFLX_LW_up,TOAFLX_LW_dn,TOAFLX_SW_up,TOAFLX_SW_dn,SFCFLX_LW_up_c,SFCFLX_LW_dn_c,SFCFLX_SW_up_c,SFCFLX_SW_dn_c,TOAFLX_LW_up_c,TOAFLX_LW_dn_c,TOAFLX_SW_up_c,TOAFLX_SW_dn_c,TOMFLX_LW_up_c,TOMFLX_LW_dn_c,TOMFLX_SW_up_c,TOMFLX_SW_dn_c,CLDFRAC,MP_Re,MP_Qe,AE_Re,AE_Qe,RH)

       call SOLARINS_insolation( IA, IS, IE, JA, JS, JE, &
                                 REAL_LON(:,:), REAL_LAT(:,:),      & ! [IN]
                                 TIME_NOWDATE(:), TIME_OFFSET_YEAR, & ! [IN]
                                 solins(:,:), cosSZA  (:,:)         ) ! [OUT]

       call ATMOS_PHY_MP_vars_get_diagnostic( &
            DENS(:,:,:), TEMP(:,:,:), QTRC(:,:,:,:), & ! [IN]
            CLDFRAC=CLDFRAC, Re=MP_Re, Qe=MP_Qe      ) ! [OUT]

       if ( RD_use_PBL_cloud ) then
          !$omp parallel do
          !$acc kernels
          do j = JS, JE
          do i = IS, IE
          do k = KS, KE
             if ( REAL_CZ(k,i,j) < ATMOS_PHY_BL_Zi(i,j) + REAL_FZ(KS-1,i,j) ) then
                CLDFRAC(k,i,j) = ATMOS_PHY_BL_cldfrac(k,i,j)
                if ( ATMOS_PHY_BL_cldfrac(k,i,j) > EPS ) then
                   MP_Qe(k,i,j,I_HC) = ATMOS_PHY_BL_QL(k,i,j) / ATMOS_PHY_BL_cldfrac(k,i,j)
                else
                   MP_Qe(k,i,j,I_HC) = 0.0_RP
                end if
                MP_Qe(k,i,j,I_HI) = 0.0_RP
                MP_Re(k,i,j,I_HC) = 8.0E-4
             end if
          end do
          end do
          end do
          !$acc end kernels
          
       end if

       call ATMOS_vars_get_diagnostic( "RH", RH )
       call ATMOS_PHY_AE_vars_get_diagnostic( &
            QTRC(:,:,:,:), RH(:,:,:), & ! [IN]
            Re=AE_Re, Qe=AE_Qe        ) ! [IN]


       select case ( ATMOS_PHY_RD_TYPE )
       case ( "MSTRNX" )

          if (      FILE_HISTORY_query_f(HIST_id_SFLX_LW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_SFLX_LW_dn_c) &
               .or. FILE_HISTORY_query_f(HIST_id_SFLX_SW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_SFLX_SW_dn_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOAFLX_LW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOAFLX_LW_dn_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOAFLX_SW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOAFLX_SW_dn_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOMFLX_LW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOMFLX_LW_dn_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOMFLX_SW_up_c) &
               .or. FILE_HISTORY_query_f(HIST_id_TOMFLX_SW_dn_c) ) then
             clear_sky = .true.
          else
             clear_sky = .false.
          end if

          call ATMOS_PHY_RD_MSTRNX_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               DENS(:,:,:), TEMP(:,:,:), PRES(:,:,:), QV(:,:,:), & ! [IN]
               REAL_CZ(:,:,:), REAL_FZ(:,:,:),                   & ! [IN]
               fact_ocean(:,:), fact_land(:,:), fact_urban(:,:), & ! [IN]
               SFC_TEMP(:,:), SFC_albedo(:,:,:,:),               & ! [IN]
               solins(:,:), cosSZA(:,:),                         & ! [IN]
               CLDFRAC(:,:,:), MP_Re(:,:,:,:), MP_Qe(:,:,:,:),   & ! [IN]
               AE_Re(:,:,:,:), AE_Qe(:,:,:,:),                   & ! [IN]
               flux_rad(:,:,:,:,:,:),                            & ! [OUT]
               flux_rad_top(:,:,:,:,:), SFLX_rad_dn(:,:,:,:),    & ! [OUT]
               clear_sky = clear_sky,                            & ! [IN, optional]
               dtau_s = dtau_s(:,:,:), dem_s = dem_s(:,:,:)      ) ! [OUT, optional]

       case ( "OFFLINE" )

          call ATMOS_PHY_RD_offline_flux( &
               KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               TIME_NOWDAYSEC,        & ! [IN]
               flux_rad(:,:,:,:,:,2), & ! [OUT]
               SFLX_rad_dn(:,:,:,:)   ) ! [OUT]
          !$acc kernels
          clear_sky = .true.
          flux_rad(:,:,:,:,:,I_ClearSky)   = 0.0_RP ! clear sky
          flux_rad_top(:,:,:,:,:) = 0.0_RP
          dtau_s(:,:,:) = 0.0_RP
          dem_s(:,:,:) = 0.0_RP
          !$acc end kernels

       end select

       ! ### for SPRINTARS

       !$omp parallel do collapse(2)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          dtau_s_sprintars(k,i,j) = dtau_s(k,i,j)
       end do
       end do
       end do
       !$acc end kernels

       ! surface
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          SFCFLX_SW_dn_net_sprintars(i,j) = 0.0_RP
       end do
       end do
       !$acc end kernels

!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          if ( clear_sky ) then
             SFCFLX_LW_up_c(i,j)    = flux_rad(KS-1,i,j,I_LW,I_up,I_ClearSky)
             SFCFLX_LW_dn_c(i,j)    = flux_rad(KS-1,i,j,I_LW,I_dn,I_ClearSky)
             SFCFLX_SW_up_c(i,j)    = flux_rad(KS-1,i,j,I_SW,I_up,I_ClearSky)
             SFCFLX_SW_dn_c(i,j)    = flux_rad(KS-1,i,j,I_SW,I_dn,I_ClearSky)
          end if
          ! for all-sky
          SFCFLX_LW_up  (i,j)    = flux_rad(KS-1,i,j,I_LW,I_up,I_Cloud)
          SFCFLX_LW_dn  (i,j)    = flux_rad(KS-1,i,j,I_LW,I_dn,I_Cloud)
          SFCFLX_SW_up  (i,j)    = flux_rad(KS-1,i,j,I_SW,I_up,I_Cloud)
          SFCFLX_SW_dn  (i,j)    = flux_rad(KS-1,i,j,I_SW,I_dn,I_Cloud)

          flux_net_sfc(i,j,I_LW) = SFCFLX_LW_up(i,j) - SFCFLX_LW_dn(i,j)
          flux_net_sfc(i,j,I_SW) = SFCFLX_SW_up(i,j) - SFCFLX_SW_dn(i,j)

          SFCFLX_SW_dn_net_sprintars(i,j) = - flux_net_sfc(i,j,I_SW) ! ### for SPRINTARS

       enddo
       enddo
       !$acc end kernels

       ! top of the atmosphere
!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          if ( clear_sky ) then
             TOAFLX_LW_up_c(i,j)    = flux_rad_top(i,j,I_LW,I_up,I_ClearSky)
             TOAFLX_LW_dn_c(i,j)    = flux_rad_top(i,j,I_LW,I_dn,I_ClearSky)
             TOAFLX_SW_up_c(i,j)    = flux_rad_top(i,j,I_SW,I_up,I_ClearSky)
             TOAFLX_SW_dn_c(i,j)    = flux_rad_top(i,j,I_SW,I_dn,I_ClearSky)
          end if
          ! for all-sky
          TOAFLX_LW_up  (i,j)    = flux_rad_top(i,j,I_LW,I_up,I_Cloud)
          TOAFLX_LW_dn  (i,j)    = flux_rad_top(i,j,I_LW,I_dn,I_Cloud)
          TOAFLX_SW_up  (i,j)    = flux_rad_top(i,j,I_SW,I_up,I_Cloud)
          TOAFLX_SW_dn  (i,j)    = flux_rad_top(i,j,I_SW,I_dn,I_Cloud)

          flux_net_toa(i,j,I_LW) = TOAFLX_LW_up(i,j) - TOAFLX_LW_dn(i,j)
          flux_net_toa(i,j,I_SW) = TOAFLX_SW_up(i,j) - TOAFLX_SW_dn(i,j)
       enddo
       enddo
       !$acc end kernels

       ! top of the model
!OCL XFILL
       !$omp parallel do
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
          ! for clear-sky
          if ( clear_sky ) then
             TOMFLX_LW_up_c(i,j)    = flux_rad(KE,i,j,I_LW,I_up,I_ClearSky)
             TOMFLX_LW_dn_c(i,j)    = flux_rad(KE,i,j,I_LW,I_dn,I_ClearSky)
             TOMFLX_SW_up_c(i,j)    = flux_rad(KE,i,j,I_SW,I_up,I_ClearSky)
             TOMFLX_SW_dn_c(i,j)    = flux_rad(KE,i,j,I_SW,I_dn,I_ClearSky)
          end if
          ! for all-sky
          TOMFLX_LW_up  (i,j)    = flux_rad(KE,i,j,I_LW,I_up,I_Cloud)
          TOMFLX_LW_dn  (i,j)    = flux_rad(KE,i,j,I_LW,I_dn,I_Cloud)
          TOMFLX_SW_up  (i,j)    = flux_rad(KE,i,j,I_SW,I_up,I_Cloud)
          TOMFLX_SW_dn  (i,j)    = flux_rad(KE,i,j,I_SW,I_dn,I_Cloud)

          flux_net_tom(i,j,I_LW) = TOMFLX_LW_up(i,j) - TOMFLX_LW_dn(i,j)
          flux_net_tom(i,j,I_SW) = TOMFLX_SW_up(i,j) - TOMFLX_SW_dn(i,j)
       enddo
       enddo
       !$acc end kernels

!OCL XFILL
       !$omp parallel do collapse(2)
       !$acc kernels
       do j = JS, JE
       do i = IS, IE
       do k = KS, KE
          flux_up (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_up,I_Cloud) + flux_rad(k,i,j,I_LW,I_up,I_Cloud) )
          flux_dn (k,i,j,I_LW) = 0.5_RP * ( flux_rad(k-1,i,j,I_LW,I_dn,I_Cloud) + flux_rad(k,i,j,I_LW,I_dn,I_Cloud) )
          flux_up (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_up,I_Cloud) + flux_rad(k,i,j,I_SW,I_up,I_Cloud) )
          flux_dn (k,i,j,I_SW) = 0.5_RP * ( flux_rad(k-1,i,j,I_SW,I_dn,I_Cloud) + flux_rad(k,i,j,I_SW,I_dn,I_Cloud) )

          flux_net(k,i,j,I_LW) = flux_up(k,i,j,I_LW) - flux_dn(k,i,j,I_LW)
          flux_net(k,i,j,I_SW) = flux_up(k,i,j,I_SW) - flux_dn(k,i,j,I_SW)
       enddo
       enddo
       enddo
       !$acc end kernels

       ! apply radiative flux convergence -> heating rate
       call ATMOS_PHY_RD_calc_heating( &
            KA, KS, KE, IA, IS, IE, JA, JS, JE, &
            flux_rad(:,:,:,:,:,I_Cloud),     & ! [IN]
            DENS(:,:,:), TEMP(:,:,:),  & ! [IN]
            CVtot(:,:,:),              & ! [IN]
            REAL_FZ(:,:,:),            & ! [IN]
            RHOH_RD(:,:,:),            & ! [OUT]
            temp_t = TEMP_t(:,:,:,:)   ) ! [OUT]



       call FILE_HISTORY_put( HIST_id_SOLINS, solins(:,:) )
       call FILE_HISTORY_put( HIST_id_COSZ,   cosSZA(:,:) )

       call FILE_HISTORY_put( HIST_id_SFLX_LW_up_c, SFCFLX_LW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_LW_dn_c, SFCFLX_LW_dn_c(:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_SW_up_c, SFCFLX_SW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_SW_dn_c, SFCFLX_SW_dn_c(:,:) )

       call FILE_HISTORY_put( HIST_id_SFLX_LW_up,   SFCFLX_LW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_LW_dn,   SFCFLX_LW_dn  (:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_SW_up,   SFCFLX_SW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_SFLX_SW_dn,   SFCFLX_SW_dn  (:,:) )

       call FILE_HISTORY_put( HIST_id_SFLX_LW_net,  flux_net_sfc  (:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_SFLX_SW_net,  flux_net_sfc  (:,:,I_SW) )

       call FILE_HISTORY_put( HIST_id_TOAFLX_LW_up_c, TOAFLX_LW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_LW_dn_c, TOAFLX_LW_dn_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_SW_up_c, TOAFLX_SW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_SW_dn_c, TOAFLX_SW_dn_c(:,:) )

       call FILE_HISTORY_put( HIST_id_TOAFLX_LW_up,   TOAFLX_LW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_LW_dn,   TOAFLX_LW_dn  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_SW_up,   TOAFLX_SW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_SW_dn,   TOAFLX_SW_dn  (:,:) )

       call FILE_HISTORY_put( HIST_id_TOAFLX_LW_net,  flux_net_toa  (:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_TOAFLX_SW_net,  flux_net_toa  (:,:,I_SW) )

       call FILE_HISTORY_put( HIST_id_TOMFLX_LW_up_c, TOMFLX_LW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_LW_dn_c, TOMFLX_LW_dn_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_SW_up_c, TOMFLX_SW_up_c(:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_SW_dn_c, TOMFLX_SW_dn_c(:,:) )

       call FILE_HISTORY_put( HIST_id_TOMFLX_LW_up,   TOMFLX_LW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_LW_dn,   TOMFLX_LW_dn  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_SW_up,   TOMFLX_SW_up  (:,:) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_SW_dn,   TOMFLX_SW_dn  (:,:) )

       call FILE_HISTORY_put( HIST_id_TOMFLX_LW_net,  flux_net_tom  (:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_TOMFLX_SW_net,  flux_net_tom  (:,:,I_SW) )

       call FILE_HISTORY_put( HIST_id_SLR, flux_net_sfc(:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_SSR, flux_net_sfc(:,:,I_SW) )
       call FILE_HISTORY_put( HIST_id_OLR, TOAFLX_LW_up(:,:) )
       call FILE_HISTORY_put( HIST_id_OSR, TOAFLX_SW_up(:,:) )

       call FILE_HISTORY_put( HIST_id_RADFLUX_LWUP, flux_up (:,:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_RADFLUX_LWDN, flux_dn (:,:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_RADFLUX_LW,   flux_net(:,:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_RADFLUX_SWUP, flux_up (:,:,:,I_SW) )
       call FILE_HISTORY_put( HIST_id_RADFLUX_SWDN, flux_dn (:,:,:,I_SW) )
       call FILE_HISTORY_put( HIST_id_RADFLUX_SW,   flux_net(:,:,:,I_SW) )

       call FILE_HISTORY_put( HIST_id_SFLX_IR_dn_dir,  SFLX_rad_dn(:,:,I_R_direct ,I_R_IR ) )
       call FILE_HISTORY_put( HIST_id_SFLX_IR_dn_dif,  SFLX_rad_dn(:,:,I_R_diffuse,I_R_IR ) )
       call FILE_HISTORY_put( HIST_id_SFLX_NIR_dn_dir, SFLX_rad_dn(:,:,I_R_direct ,I_R_NIR) )
       call FILE_HISTORY_put( HIST_id_SFLX_NIR_dn_dif, SFLX_rad_dn(:,:,I_R_diffuse,I_R_NIR) )
       call FILE_HISTORY_put( HIST_id_SFLX_VIS_dn_dir, SFLX_rad_dn(:,:,I_R_direct ,I_R_VIS) )
       call FILE_HISTORY_put( HIST_id_SFLX_VIS_dn_dif, SFLX_rad_dn(:,:,I_R_diffuse,I_R_VIS) )

       call FILE_HISTORY_put( HIST_id_TEMP_t_rd_LW, TEMP_t (:,:,:,I_LW) )
       call FILE_HISTORY_put( HIST_id_TEMP_t_rd_SW, TEMP_t (:,:,:,I_SW) )
       call FILE_HISTORY_put( HIST_id_TEMP_t_rd,    TEMP_t (:,:,:,3   ) )
       call FILE_HISTORY_put( HIST_id_RHOH_RD,      RHOH_RD(:,:,:) )

       ! output of raw data, for offline output
       call FILE_HISTORY_put( HIST_id_RFLX_LW_up, flux_rad(:,:,:,I_LW,I_up,I_Cloud) )
       call FILE_HISTORY_put( HIST_id_RFLX_LW_dn, flux_rad(:,:,:,I_LW,I_dn,I_Cloud) )
       call FILE_HISTORY_put( HIST_id_RFLX_SW_up, flux_rad(:,:,:,I_SW,I_up,I_Cloud) )
       call FILE_HISTORY_put( HIST_id_RFLX_SW_dn, flux_rad(:,:,:,I_SW,I_dn,I_Cloud) )

       call FILE_HISTORY_put( HIST_id_dtau_s, dtau_s(:,:,:) )
       call FILE_HISTORY_put( HIST_id_dem_s,  dem_s (:,:,:) )

       !$acc end data

    endif

    !$acc kernels
    do j = JS, JE
    do i = IS, IE
    do k = KS, KE
       RHOH(k,i,j) = RHOH(k,i,j) + RHOH_RD(k,i,j)
    enddo
    enddo
    enddo
    !$acc end kernels

    if ( STATISTICS_checktotal ) then
       call STATISTICS_total( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                              RHOH_RD(:,:,:), 'RHOH_RD',          &
                              ATMOS_GRID_CARTESC_REAL_VOL(:,:,:), &
                              ATMOS_GRID_CARTESC_REAL_TOTVOL      )
    endif

    return
  end subroutine ATMOS_PHY_RD_driver_calc_tendency

end module mod_atmos_phy_rd_driver
