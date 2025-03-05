!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - admin
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_admin
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
       PRC_myrank,   &
       PRC_IsMaster, &
       PRC_abort
  use scale_statistics, only: &
       STATISTICS_detail
  use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE
  use scale_const, only: &
       CONST_UNDEF2,            &
       CONST_UNDEF,             &
       RAIR   => CONST_Rdry,    & !< specific gas constant (dry air)          [J/kg/K]
       CP     => CONST_CPdry,   & !< specific heat (dry air, const. pressure) [J/kg/K]
       EPSTV  => CONST_EPSTvap, & !< 1 / epsilon - 1
       Rvap   => CONST_Rvap,    & !< specific gas constant (water vapor)      [J/kg/K]
       EPSvap => CONST_EPSvap,  & !< Rdry / Rvap
       TMELT  => CONST_TMELT,   & !< melting point of water
       EPS    => CONST_EPS,     & !< small number = 1.E-16_RP
       PI     => CONST_PI,      & !< pi
       GRAV   => CONST_GRAV,    & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       EL     => CONST_LHV0       !< latent heat of vaporizaion at 0C [J/kg] = 2.501E+6_RP
  use scale_landuse, only: &
       FRAC_PFT            => LANDUSE_frac_PFT,   & !< fraction of PFT for each mosaic
       INDEX_PFT           => LANDUSE_index_PFT,  & !< index of PFT for each mosaic
       FACT_LAND           => LANDUSE_fact_land,  & !< land  factor : fact_lake + fact_urban + fact_land + fact_ocean = 1
       FACT_OCEAN          => LANDUSE_fact_ocean, & !< ocean factor : fact_lake + fact_urban + fact_land + fact_ocean = 1
       FACT_URBAN          => LANDUSE_fact_urban, & !< urban factor : fact_lake + fact_urban + fact_land + fact_ocean = 1
       FACT_LAKE           => LANDUSE_fact_lake,  & !< lake  factor : fact_lake + fact_urban + fact_land + fact_ocean = 1
       FLND                => LANDUSE_frac_land,  & !< frac_land = ( fact_lake + fact_urban + fact_land ) / ( fact_lake + fact_urban + fact_land + fact_ocean)
       FLAKE               => LANDUSE_frac_lake,  & !< frac_lake = fact_lake  / ( fact_lake + fact_urban + fact_land )
       FACT_LAKE_SPRINTARS => LANDUSE_fact_lake_SPRINTARS
  use scale_topography, only: & 
       GDZS => TOPOGRAPHY_Zsfc  !< absolute ground height [m]
  use scale_atmos_grid_cartesC_real, only: &
       GDZ      => ATMOS_GRID_CARTESC_REAL_CZ,  & !< GDZ (KA,  IA,JA) : altitude [m]
       GDZM     => ATMOS_GRID_CARTESC_REAL_FZ,  & !< GDZM(0:KA,IA,JA) : altitude (half level) [m]
       LON_REAL => ATMOS_GRID_CARTESC_REAL_LON, & !< longitude [rad,0-2pi]
       LAT_REAL => ATMOS_GRID_CARTESC_REAL_LAT    !< latitude  [rad,-pi,pi]
  use scale_atmos_saturation, only:   &
       ATMOS_SATURATION_psat_all,     &
       ATMOS_SATURATION_dens2qsat_all
  use scale_file_history, only: &
       FILE_HISTORY_in
  use scale_file_cartesC, only: &
       FILE_CARTESC_open, &
       FILE_CARTESC_read, &
       FILE_CARTESC_close
  use scale_time, only: & 
       TIME_NOWDATE             !< current time [YYYY MM DD HH MM SS]: TIME_NOWDATE(6)
  use scale_atmos_aerosol, only: &
       N_AE                     !< number of aerosol types
                                ! 1->soil dust, 2->carbonaceous(BC/OC=0.30), 3->carbonaceous(BC/OC=0.15), 4->carbonaceous(BC/OC=0.00),
                                ! 5->black carbon, 6->sulfate, 7->sea salt
  use scale_atmos_phy_ae_SPRINTARS_aerosulf, only: &
       ATMOS_PHY_AE_SPRINTARS_AEROSU_setup, &
       ATMOS_PHY_AE_SPRINTARS_AEROSU_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AEROSU
  use scale_atmos_phy_ae_SPRINTARS_aerocarb, only: &
       ATMOS_PHY_AE_SPRINTARS_AEROCA_setup, &
       ATMOS_PHY_AE_SPRINTARS_AEROCA_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AEROCA
  use scale_atmos_phy_ae_SPRINTARS_aerodust, only: &
       ATMOS_PHY_AE_SPRINTARS_AERODU_setup, &
       ATMOS_PHY_AE_SPRINTARS_AERODU_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AERODU
  use scale_atmos_phy_ae_SPRINTARS_aerosalt, only: &
       ATMOS_PHY_AE_SPRINTARS_AEROSA_setup, &
       ATMOS_PHY_AE_SPRINTARS_AEROSA_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AEROSA
  use scale_atmos_phy_ae_SPRINTARS_aerotrac, only:   &
       ATMOS_PHY_AE_SPRINTARS_AEROAT_setup,          &
       ATMOS_PHY_AE_SPRINTARS_AEROAT_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AEROAT,                &
       ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY2
  use scale_atmos_phy_ae_SPRINTARS_ccn, only: &
       ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_setup, &
       ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_finalize, &
       ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN
  use scale_atmos_phy_ae_SPRINTARS_common, only: & 
       ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA, &
       ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY, &
       NSU, NCA, NDU, NSA, NAT,                 &
       QA_AE_SU, QS_AE_SU, QE_AE_SU,            &
       QA_AE_CA, QS_AE_CA, QE_AE_CA,            &
       QA_AE_DU, QS_AE_DU, QE_AE_DU,            &
       QA_AE_SA, QS_AE_SA, QE_AE_SA,            &
       QA_AE_AT, QS_AE_AT, QE_AE_AT,            &
       IWA,                                     &
       !VMISS,                                   &
       VTR,                                     &
       RADR,                                    &
       THRESRe,                                 &
       DENSR,                                   & !< rain droplet density     [kg/m3]
       DENSI,                                   & !< density of ice           [kg/m3]
       THRESVTR

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSL_setup
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSL_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSL
  
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
  ! parameter
  integer,  private, parameter   :: AIRLEV = 25                                     !< original vertical level for airfuel
  integer,  private, parameter   :: KCPCL  = 2                                      !< cloud particle
  real(RP), private, parameter   :: WAVEL(IWA) = (/ 0.550_RP, 0.440_RP, 0.870_RP /) !< wavelengths for optical thickness
  real(RP), private, parameter   :: CDVEMIN = 4.0E-5_RP                             !< for CDVE
  real(RP), private, parameter   :: CDVEMAX = 1.0E+3_RP                             !< for CDVE

  ! variables
  real(RP), private, allocatable :: OUTQLTDU(:,:,:,:) !(KA,IA,JA,NDU)   !< dust for radiation 
  real(RP), private, allocatable :: OUTQLTSA(:,:,:,:) !(KA,IA,JA,NSA)   !< sea salt for radiation
  real(RP), private, allocatable :: OUTQLTCA(:,:,:,:) !(KA,IA,JA,NCA)   !< aerosol for radiation
  real(RP), private, allocatable :: NUMCON  (:,:,:,:) !(KA,IA,JA,  3)   !< aerosol number concentration
  real(RP), private, allocatable :: EMITSA    (:,:,:) !   (IA,JA,NSA)   !< sea salt emission flux

  real(RP), private, allocatable, save :: ACTI(:,:,:,:) !(KA,IA,JA,6)   !< fraction of activated aerosols

  real(RP), private, allocatable, save :: FRAIN_MP_RAIN1(:,:,:) !(0:KA,IA,JA) !< FRAIN_MP_RAIN on 1 step before
  real(RP), private, allocatable, save :: FRAIN_MP_SNOW1(:,:,:) !(0:KA,IA,JA) !< FRAIN_MP_SNOW on 1 step before
  real(RP), private, allocatable, save :: DENS1         (:,:,:) !(KA,  IA,JA) !< air density   on 1 step before

  real(RP), private, allocatable, save ::  GRLAI_DATA   (:,:,:)   !   (IA,JA,12)    !< Leaf Area Index (climatology;1981-2015;ORNL)[m2/m2]
  real(RP), private, allocatable, save ::  GRSNW_DATA   (:,:,:,:) !   (IA,JA,12,31) !< snow depth (climatology;2000-2013;JRA55)[m]
  real(RP), private, allocatable, save ::  OHRAD_DATA (:,:,:,:)   !(KA,IA,JA,12)    !< OH concentration [mol/cm3]
  real(RP), private, allocatable, save :: DAYTIM_DATA   (:,:,:)   !   (IA,JA,12)    !< daytime ratio [1]
  real(RP), private, allocatable, save ::   NO3G_DATA (:,:,:,:)   !(KA,IA,JA,12)    !< NO3(gas) mixing ratio: [pptv]
  real(RP), private, allocatable, save ::  O3ORI_DATA (:,:,:,:)   !(KA,IA,JA,12)    !< O3       mixing ratio: [ppmv]
  real(RP), private, allocatable, save :: AIRBCT_DATA (:,:,:,:)   !(KA,IA,JA,day)   !< BC  emission from aircrafts
  real(RP), private, allocatable, save :: AIROCT_DATA (:,:,:,:)   !(KA,IA,JA,day)   !< OC  emission from aircrafts
  real(RP), private, allocatable, save :: AIRSO2T_DATA(:,:,:,:)   !(KA,IA,JA,day)   !< SO2 emission from aircrafts
  
  !namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL /
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_SU  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_CA  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_DU  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_SA  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_NI  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_AT  = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_CCN = .false.              ! default is .false. ( not treated )
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_ARI = .false.              ! default does not interact with radiation scheme
  logical,  private, save :: ATMOS_SW_PHY_AE_SPRINTARS_ACI = .false.              ! default does not interact with microphysics scheme

  !namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA / 
  character(len=H_LONG), private, save :: fname_aeroat = '141115_Suppl_acp_JAEAsourceterm_UTC_hgt.csv' !< file name (sulf)
  character(len=H_LONG), private, save :: fname_emission = 'sprintars_emission_inventory' !< file name (land)
  integer :: fileid

  !namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN /
  real(RP), private, save :: CCMAX  = 0.2_RP                               !
  real(RP), private, save :: RHFCT  = 1.0_RP                               ! 
  logical,  private, save :: ODCA   = .false.                              ! 
  real(RP), private, save :: FLNDMX = 0.99_RP                              !< threshold of land sea mask
  real(RP), private, save :: SBB    = 3000.0_RP                            !< emission height of biomass burning from surface [m] 
  real(RP), private, save :: WDEP_CLD = 1.0E-3_RP                          !< threshold of hydrometeor's mass for triggering wet deposition
  real(RP), private, save :: WDEP_TAU = 0.2_RP                             !< threshold of hydrometeor's optical thickness for triggering wet deposition
  real(RP), private, save :: THRE_HGT_DUST_Tibet = 1800.0_RP               !< threshold of upper limit of dust emission over Tibet area [m]
  
  ! !namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL /
  ! real(RP), private, save :: RCMIN(KCPCL) = (/  1.0E-6_RP,   1.0E-6_RP /)  !< minimum of radius (DUMMY) [m]
  ! real(RP), private, save :: RCMAX(KCPCL) = (/ 30.0E-6_RP, 200.0E-6_RP /)  !< maximum of radius (DUMMY)
  ! real(RP), private, save :: REFCT(KCPCL) = (/  1.1E+0_RP,   1.0E+0_RP /)  !< Re/Rvol           (DUMMY)
  ! real(RP), private, save :: TFLIQ =  13.50_RP                             !< parameter for FLIQ
  ! real(RP), private, save :: SFLIQ = 269.91_RP                             !< parameter for FLIQ
  ! real(RP), private, save :: MFLIQ = 235.15_RP                             !< parameter for FLIQ
  ! real(RP), private, save :: FLIQMIN = 2.0E-2_RP                           !< minimum liquid fraction (for CCN generation)
  ! real(RP), private       :: UAMIN                                         !< minimum aerosol number    (DUMMY)
  ! real(RP), private       :: CCFCT                                         !< CDCN/CCN                  (DUMMY)
  ! real(RP), private, save :: UCMIN(KCPCL) = (/  2.5E+7_RP,   1.0E+1_RP /)  !< minimum cloud number concentration (DUMMY)
  ! real(RP), private       :: UCMAX                                         !< maximum cloud drop number (DUMMY)
  ! logical,  private, save :: OCCN2CLD = .true.                             !<
  ! logical,  private, save :: OCLD2CCN = .true.                             !< use prognoced ice fraction

  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL_setup
    implicit none

    integer :: k, i, j, it, ip
    integer :: ierr
    character(len=H_MID) :: varname
    real(RP) :: work2d(IA,JA), work3d(KA,IA,JA)

    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL / &
         ATMOS_SW_PHY_AE_SPRINTARS_SU,  &
         ATMOS_SW_PHY_AE_SPRINTARS_CA,  &
         ATMOS_SW_PHY_AE_SPRINTARS_DU,  &
         ATMOS_SW_PHY_AE_SPRINTARS_SA,  &
         ATMOS_SW_PHY_AE_SPRINTARS_NI,  &
         ATMOS_SW_PHY_AE_SPRINTARS_AT,  &
         ATMOS_SW_PHY_AE_SPRINTARS_CCN, &
         ATMOS_SW_PHY_AE_SPRINTARS_ARI, &
         ATMOS_SW_PHY_AE_SPRINTARS_ACI
    
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA / &
         fname_aeroat, &
         fname_emission
    
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN / &
         CCMAX,  &
         RHFCT,  &
         ODCA,   &
         FLNDMX, &
         SBB,    &
         WDEP_TAU, &
         WDEP_CLD, &
         THRE_HGT_DUST_Tibet
    
    ! namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL / &
    !      RCMIN,    &
    !      RCMAX,    &
    !      REFCT,    &
    !      TFLIQ,    &
    !      SFLIQ,    &
    !      MFLIQ,    &
    !      FLIQMIN,  &
    !      UAMIN,    &
    !      CCFCT,    &
    !      UCMIN,    &
    !      UCMAX,    &
    !      OCCN2CLD, &
    !      OCLD2CCN
  
    !----------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Setup'

    ! Setup parameter
    !================================================

    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Setup parameter'

    !--- read namelist
    !PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL
    !-----------------
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSL)
    
    !PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA
    !-----------------
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA)

    !PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN
    !-----------------
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_NMAEMN)

    ! !PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL
    ! !-----------------
    ! rewind( IO_FID_CONF )
    ! read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL, iostat=ierr )
    ! if( ierr < 0 ) then !--- missing
    !    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL. Default used.'
    ! elseif( ierr > 0 ) then !--- fatal error
    !    LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL, Check!'
    !    call PRC_abort
    ! endif
    ! LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_NMAECL)
    !-----------------

    call FILE_CARTESC_open( fname_emission, fileid )

    !--- Setup
    !Setup each aerosol type
    !-----------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SU  ) &
         call ATMOS_PHY_AE_SPRINTARS_AEROSU_setup &
!         ( trim(adjustl(fname_chem)), trim(adjustl(fname_sulf)), &
         ( LON_REAL, LAT_REAL,                                   &
           GDZ, GDZM, QA_AE_SU,                                  &
           fileid                                                )
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA  ) &
         call ATMOS_PHY_AE_SPRINTARS_AEROCA_setup &
!         ( trim(adjustl(fname_carb)), LON_REAL, LAT_REAL, fileid )
         ( LON_REAL, LAT_REAL, fileid )
    if ( ATMOS_SW_PHY_AE_SPRINTARS_DU  ) call ATMOS_PHY_AE_SPRINTARS_AERODU_setup( QA_AE_DU )
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SA  ) call ATMOS_PHY_AE_SPRINTARS_AEROSA_setup( QA_AE_SA )
    if ( ATMOS_SW_PHY_AE_SPRINTARS_AT  ) call ATMOS_PHY_AE_SPRINTARS_AEROAT_setup( LON_REAL, LAT_REAL, GDZM, fname_aeroat, QA_AE_AT )
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CCN ) call ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_setup
    !-----------------
    
    !--- read file
    !read parameter (land) file
    !-----------------
    !Setup LAI (NETCDF)
    allocate( GRLAI_DATA(IA,JA,12) )
    GRLAI_DATA(:,:,:) = 0.0_RP
    write(varname,'(A)') 'GRLAI'
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                               work2d(:,:),                 &
                               ip                           )
       GRLAI_DATA(:,:,ip) = work2d(:,:)
    enddo
    !Setup Snow Depth (NETCDF)
    allocate( GRSNW_DATA(IA,JA,12,31) )
    GRSNW_DATA(:,:,:,:) = 0.0_RP
    do it = 1, 31
       write(varname,'(A,I0)') 'GRSNW_',it
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),                 &
                                  ip                           )
          GRSNW_DATA(:,:,ip,it) = work2d(:,:)
       enddo
    enddo
    !-----------------
    
    !read parameter (chemistry) file
    !-----------------
    !Setup chemistry data (OHRAD_DATA, NO3G_DATA, O3ORI_DATA, DAYTIM_DATA) (NETCDF)
    allocate(  OHRAD_DATA(KA,IA,JA,12) )
    allocate(   NO3G_DATA(KA,IA,JA,12) )
    allocate(  O3ORI_DATA(KA,IA,JA,12) )
    allocate( DAYTIM_DATA   (IA,JA,12) )
    OHRAD_DATA (:,:,:,:) = 0.0_RP
    NO3G_DATA  (:,:,:,:) = 0.0_RP
    O3ORI_DATA (:,:,:,:) = 0.0_RP
    DAYTIM_DATA  (:,:,:) = 0.0_RP

    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA .or. ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       write(varname,'(A)') 'OHRAD'
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'ZXY', &
                                  work3d(:,:,:),                &
                                  ip                            )
          OHRAD_DATA(:,:,:,ip) = work3d(:,:,:)
       enddo
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA ) then
       write(varname,'(A)') 'NO3G'
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'ZXY', &
                                  work3d(:,:,:),                &
                                  ip                            )
          NO3G_DATA(:,:,:,ip) = work3d(:,:,:)
       enddo
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA .or. ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       write(varname,'(A)') 'O3ORI'
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'ZXY', &
                                  work3d(:,:,:),                &
                                  ip                            )
          O3ORI_DATA(:,:,:,ip) = work3d(:,:,:)
       enddo
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA .or. ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       write(varname,'(A)') 'DAYTIM'
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY',  &
                                  work2d(:,:),                  &
                                  ip                            )
          DAYTIM_DATA(:,:,ip) = work2d(:,:)
       enddo
    endif
    do it = 1, 12
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             OHRAD_DATA(k,i,j,it) = max( OHRAD_DATA(k,i,j,it), 0.0_RP )
              NO3G_DATA(k,i,j,it) = max(  NO3G_DATA(k,i,j,it), 0.0_RP )
             O3ORI_DATA(k,i,j,it) = max( O3ORI_DATA(k,i,j,it), 0.0_RP )
          enddo
          DAYTIM_DATA(i,j,it) = max( DAYTIM_DATA(i,j,it), 0.0_RP )
       enddo
       enddo
    enddo

    call FILE_CARTESC_close( fileid )
    !-----------------
    
    !read parameter (emission) file
    !-----------------
    !Setup emission data from aircrafts (AIRBCT_DATA, AIROCT_DATA, AIRSO2T_DATA)
    !-----------------   
    !================================================
    !end Setup parameter

    
    ! Allocate arrays
    !================================================
    allocate( FRAIN_MP_RAIN1(0:KA,IA,JA) )
    allocate( FRAIN_MP_SNOW1(0:KA,IA,JA) )
    allocate( DENS1         (KA,  IA,JA) )
    allocate( EMITSA        (IA,JA,NSA)  )
    allocate( NUMCON        (KA,IA,JA,3) )
    FRAIN_MP_RAIN1(0:KA,:,:) = 0.0_RP
    FRAIN_MP_SNOW1(0:KA,:,:) = 0.0_RP
    DENS1         (:,   :,:) = 0.0_RP
    EMITSA        (:,:,:   ) = 0.0_RP
    NUMCON         (:,:,:,:) = 0.0_RP
    
    allocate( ACTI(KA,IA,JA,6) )
    ACTI(:,:,:,:) = 0.0_RP
    
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA ) then
       allocate( OUTQLTCA(KA,IA,JA,NCA) )
       OUTQLTCA(:,:,:,:) = 0.0_RP
    endif
    
    if ( ATMOS_SW_PHY_AE_SPRINTARS_DU ) then
       allocate( OUTQLTDU(KA,IA,JA,NDU) )
       OUTQLTDU(:,:,:,:) = 0.0_RP
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_SA ) then
       allocate( OUTQLTSA(KA,IA,JA,NSA) )
       OUTQLTSA(:,:,:,:) = 0.0_RP
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_NI ) then
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_AT ) then
    endif
    
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CCN ) then
    endif
    !================================================
    !end Allocate arrays

    ! Initialize history value
    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP

    call FILE_HISTORY_in( zero_3d(:,:,:), 'VT_SPRINTARS',           'virtual temperature',               'K',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'VT_HL_SPRINTARS',        'virtual temperature at half level', 'K',       dim_type = 'ZHXY' ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'PRES_HL_SPRINTARS',      'pressure at half level',            'Pa',      dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RH_SPRINTARS',           'ratio of qv to qsat',               'kg/kg',   dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'FLUX_QR_CP_SPRINTARS',   'rain flux (CP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'FLUX_QS_CP_SPRINTARS',   'snow flux (CP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'FLUX_QR_MP_SPRINTARS',   'rain flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'FLUX_QS_MP_SPRINTARS',   'snow flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'PREC_SPRINTARS',         'precipitation flux',                'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'PREC_CP_SPRINTARS',      'precipitation flux (CP)',           'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RAIN_MP_SPRINTARS',      'rain flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'SNOW_MP_SPRINTARS',      'snow flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'VTR_MP_RAIN_SPRINTARS',  'terminal velocity (MP:rain)',       'm/s',     dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'VTR_MP_SNOW_SPRINTARS',  'terminal velocity (MP:snow)',       'm/s',     dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CLDFRAC_CP_SPRINTARS',   'cloud fraction (CP)',               '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CLDFRAC_MP_SPRINTARS',   'cloud fraction (MP)',               '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'TCLDF_SPRINTARS',        'TCLDF_SPRINTARS',                   '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RNWR_SPRINTARS',         'RNWR_SPRINTARS',                    '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RNWR_CP_SPRINTARS',      'RNWR_CP_SPRINTARS',                 '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RNWR_MP_RAIN_SPRINTARS', 'RNWR_MP_RAIN_SPRINTARS',            '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( zero_3d(:,:,:), 'RNWR_MP_SNOW_SPRINTARS', 'RNWR_MP_SNOW_SPRINTARS',            '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( zero_3d(:,:,:), 'AIRFEL_SPRINTARS',       'AIRFEL_SPRINTARS',                  '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'AIRFELOC_SPRINTARS',     'AIRFELOC_SPRINTARS',                '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'AIRFELSO2_SPRINTARS',    'AIRFELSO2_SPRINTARS',               '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'OHRDC_SPRINTARS',        'OH concentration',                  '1/cm3',   dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'OHRAD_SPRINTARS',        'OH concentration',                  '1/cm3',   dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'O3ORI_SPRINTARS',        'O3 mixing ratio',                   'ppmv',    dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NO3G_SPRINTARS',         'NO3(gas) mixing ratio',             'pptv',    dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( zero_2d(:,:  ), 'Uabs10T_SPRINTARS',      '10m absolute wind',                 'm/s',     dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'CDVE_SPRINTARS',         'CDVE_SPRINTARS',                    '1',       dim_type = 'XY'   ) !
    call FILE_HISTORY_in( zero_2d(:,:  ), 'GRWG_SPRINTARS',         'GRWG_SPRINTARS',                    'm3/m3',   dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'CCOVMR_SPRINTARS',       'CCOVMR_SPRINTARS',                  '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'GRLAI_SPRINTARS',        'leaf area index',                   'm2/m2',   dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'GRSNW_SPRINTARS',        'snow depth',                        'm',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QDRY_SPRINTARS',         'QDRY_SPRINTARS',                    'kg/kg',   dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QV_SPRINTARS',           'QV_SPRINTARS',                      'kg/kg',   dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QC_SPRINTARS',           'QC_SPRINTARS',                      'kg/kg',   dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QI_SPRINTARS',           'QI_SPRINTARS',                      'kg/kg',   dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'GDW_SPRINTARS',          'vertical velocity',                 'm/s',     dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_3d(:,:,:), 'TKE_BL_SPRINTARS',       'TKE(BL)',                           'm2/s2',   dim_type = 'ZXY'  ) !
    call FILE_HISTORY_in( zero_2d(:,:),   'LAND_WATER_SPRINTARS',   'LAND_WATER_SPRINTARS',              'm3/m3',   dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:),   'LAND_PROPERTY_SPRINTARS','LAND_PROPERTY_SPRINTARS',           'm3/m3',   dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:),   'GOICE_SPRINTARS',        'GOICE_SPRINTARS',                   '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:),   'SFCFLX_SW_net_SPRINTARS','SFCFLX_SW_net_SPRINTARS',           'W/m2',    dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:),   'RCOSZ_SPRINTARS',        'RCOSZ_SPRINTARS',                   '1',       dim_type = 'XY'   ) !
    call FILE_HISTORY_in( zero_2d(:,:),   'FACT_OCEAN1_SPRINTARS',  'FACT_OCEAN1_SPRINTARS',             '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:),   'FACT_LAKE_SPRINTARS',    'FACT_LAKE_SPRINTARS',               '1',       dim_type = 'XY'   ) !
    call FILE_HISTORY_in( zero_2d(:,:  ), 'INDEX_PFT1_SPRINTARS',   'INDEX_PFT1_SPRINTARS',              '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'INDEX_PFT2_SPRINTARS',   'INDEX_PFT2_SPRINTARS',              '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'FRAC_PFT1_SPRINTARS',    'FRAC_PFT1_SPRINTARS',               '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_2d(:,:  ), 'FRAC_PFT2_SPRINTARS',    'FRAC_PFT2_SPRINTARS',               '1',       dim_type = 'XY'   ) !check OK
    call FILE_HISTORY_in( zero_3d(:,:,:), 'PM25_SPRINTARS',         'mass concentration (PM2.5)',        'kg/m3',   dim_type = 'ZXY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QPM25_SPRINTARS',        'mass mixing ratio (PM2.5)',         'kg/kg',   dim_type = 'ZXY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'WATRAE_SPRINTARS',       'aerosol water',                     'kg/kg',   dim_type = 'ZXY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUT1_SPRINTARS',        'AOT (total;550nm;allsky)',                '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTA1_SPRINTARS',       'AOT (total;550nm;allsky;abs)',            '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTF1_SPRINTARS',       'AOT (total;550nm;allsky;fine)',           '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'ALFA_SPRINTARS',         'angstrom exponent (440-870nm;allsky)',    '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'SSA_SPRINTARS',          'single scatter albedo (550nm;allsky)',    '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUTD_SPRINTARS',        'AOT (total;550nm;allsky;dry)',            '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUT2_SPRINTARS',        'AOT (total;440nm;allsky)',                '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTA2_SPRINTARS',       'AOT (total;440nm;allsky;abs)',            '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUT3_SPRINTARS',        'AOT (total;870nm;allsky)',                '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTA3_SPRINTARS',       'AOT (total;870nm;allsky;abs)',            '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTC1_SPRINTARS',       'AOT (total;550nm;clearsky)',              '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTAC1_SPRINTARS',      'AOT (total;550nm;clearsky;abs)',          '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTFC1_SPRINTARS',      'AOT (total;550nm;clearsky;fine)',         '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'ALFAC_SPRINTARS',        'angstrom exponent (440-870nm;clearsky)',  '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'SSAC_SPRINTARS',         'single scatter albedo (550nm;clearsky)',  '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTC2_SPRINTARS',       'AOT (total;440nm;clearsky)',              '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTAC2_SPRINTARS',      'AOT (total;440nm;clearsky;abs)',          '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTC3_SPRINTARS',       'AOT (total;870nm;clearsky)',              '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUTAC3_SPRINTARS',      'AOT (total;870nm;clearsky;abs)',          '1',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTT1_SPRINTARS',        'extinction  (total;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABST1_SPRINTARS',        'absorption  (total;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKT1_SPRINTARS',        'backscatter (total;532nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTF1_SPRINTARS',       'extinction  (total;532nm;allsky;fine)',   '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTD1_SPRINTARS',       'extinction  (total;532nm;allsky;dry)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTT2_SPRINTARS',        'extinction  (total;1064nm;allsky)',       '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABST2_SPRINTARS',        'absorption  (total;1064nm;allsky)',       '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKT2_SPRINTARS',        'backscatter (total;1064nm;allsky)',       '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTT3_SPRINTARS',        'extinction  (total;355nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABST3_SPRINTARS',        'absorption  (total;355nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKT3_SPRINTARS',        'backscatter (total;355nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTC1_SPRINTARS',       'extinction  (total;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSTC1_SPRINTARS',       'absorption  (total;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKTC1_SPRINTARS',       'backscatter (total;532nm;clearsky)',      '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTFC1_SPRINTARS',      'extinction  (total;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTDC1_SPRINTARS',      'extinction  (total;532nm;clearsky;dry)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTC2_SPRINTARS',       'extinction  (total;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSTC2_SPRINTARS',       'absorption  (total;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKTC2_SPRINTARS',       'backscatter (total;1064nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTTC3_SPRINTARS',       'extinction  (total;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSTC3_SPRINTARS',       'absorption  (total;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKTC3_SPRINTARS',       'backscatter (total;355nm;clearsky)',      '1/m/sr', dim_type = 'ZXY' )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL_setup
  !-----------------------------------------------------------------------------

  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL_finalize
    implicit none

    deallocate( GRLAI_DATA )
    deallocate( GRSNW_DATA )

    deallocate(  OHRAD_DATA )
    deallocate(   NO3G_DATA )
    deallocate(  O3ORI_DATA )
    deallocate( DAYTIM_DATA )

    deallocate( FRAIN_MP_RAIN1 )
    deallocate( FRAIN_MP_SNOW1 )
    deallocate( DENS1          )
    deallocate( ACTI )
    deallocate( EMITSA   )
    deallocate( NUMCON )

    if ( ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROSU_finalize
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROCA_finalize
       deallocate( OUTQLTCA )
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_DU ) then
       call ATMOS_PHY_AE_SPRINTARS_AERODU_finalize
       deallocate( OUTQLTDU )
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_SA ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROSA_finalize
       deallocate( OUTQLTSA )
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_AT ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROAT_finalize
    endif

    if ( ATMOS_SW_PHY_AE_SPRINTARS_CCN ) then
       call ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_finalize
    endif

    deallocate( zero_2d )
    deallocate( zero_3d )

  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL_finalize
    
  !-----------------------------------------------------------------------------
  !> SPRINTARS Admin
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL( &
       QTRC,             & ! [MODIFIED]
       UAPCL,            & ! [OUT]
       NUMAET,           & ! [OUT]
       AE_Re,            & ! [OUT]
       AE_Qe,            & ! [OUT]
       GDP,              & ! [IN] 3D
       GDT,              & ! [IN]
       DENS,             & ! [IN]
       GDW,              & ! [IN]
       TKE_BL,           & ! [IN]
       QDRY,             & ! [IN]
       QV,               & ! [IN]
       QC,               & ! [IN]
       QI,               & ! [IN]
       THETA,            & ! [IN]
       RAIN_MP_Re,       & ! [IN]
       SNOW_MP_Re,       & ! [IN]
       VTERM_MP_RAIN,    & ! [IN]
       VTERM_MP_SNOW,    & ! [IN]
       RHOQ_t_MP_RAIN,   & ! [IN]
       RHOQ_t_MP_SNOW,   & ! [IN]
       RHOQ_t_CP_PREC,   & ! [IN]
       flux_qr_cp,       & ! [IN]
       flux_qs_cp,       & ! [IN]
       cldfrac_cp,       & ! [IN]
       flux_qr_mp,       & ! [IN]
       flux_qs_mp,       & ! [IN]
       cldfrac_mp,       & ! [IN]
       ATMOS_PHY_MP_TYPE, & ! [IN]
       NUMDEN_MP_CLDW,   & ! [IN]
       NUMDEN_MP_RAIN,   & ! [IN]
       NUMDEN_MP_CLDI,   & ! [IN]
       NUMDEN_MP_SNOW,   & ! [IN]
       dtau_s,           & ! [IN]
       GDPS,             & ! [IN] 2D
       T2,               & ! [IN]
       U10,              & ! [IN]
       V10,              & ! [IN]
       U10L,             & ! [IN]
       V10L,             & ! [IN]
       USTAR,            & ! [IN]
       WSTAR,            & ! [IN]
       QSTAR,            & ! [IN]
       SFC_TEMP,         & ! [IN]
       LAND_WATER,       & ! [IN]
       LAND_PROPERTY,    & ! [IN]
       OCEAN_ICE_FRAC,   & ! [IN]
       SFCFLX_SW_dn_net, & ! [IN]
       RCOSZ,            & ! [IN]
       TRCRNAME,         & ! [IN]
       QA_AE,            & ! [IN] parameter
       dt_AE             ) ! [IN]
   
    implicit none

    ! input
    !-------------------------
    !3D
    real(RP), intent(in)  :: GDP           (KA,  IA,JA) !< pressure                                         [Pa]
    real(RP), intent(in)  :: GDT           (KA,  IA,JA) !< temperature                                      [K]
    real(RP), intent(in)  :: DENS          (KA,  IA,JA) !< air density                                      [kg/m3]
    real(RP), intent(in)  :: GDW           (KA,  IA,JA) !< vertical velocity                                [m/s]
    real(RP), intent(in)  :: TKE_BL        (KA,  IA,JA) !< turbulent kinetic energy                         [m2/s2]
    real(RP), intent(in)  :: QDRY          (KA,  IA,JA) !< dry air                                          [kg/kg]
    real(RP), intent(in)  :: QV            (KA,  IA,JA) !< water vapor                                      [kg/kg]
    real(RP), intent(in)  :: QC            (KA,  IA,JA) !< cloud water                                      [kg/kg]
    real(RP), intent(in)  :: QI            (KA,  IA,JA) !< ice water cloud                                  [kg/kg]
    real(RP), intent(in)  :: THETA         (KA,  IA,JA) !< potential temperature                            [K]
    real(RP), intent(in)  :: RAIN_MP_Re    (KA,  IA,JA) !< effective radius of rain (microphysics scheme)   [m]
    real(RP), intent(in)  :: SNOW_MP_Re    (KA,  IA,JA) !< effective radius of snow (microphysics scheme)   [m]
    real(RP), intent(in)  :: VTERM_MP_RAIN (KA,  IA,JA) !< terminal velocity of rain (microphysics scheme)  [m/s]
    real(RP), intent(in)  :: VTERM_MP_SNOW (KA,  IA,JA) !< terminal velocity of snow (microphysics scheme)  [m/s]
    real(RP), intent(in)  :: RHOQ_t_MP_RAIN(KA,  IA,JA) !< tendency rho*QR (MP)                             [kg/m3/s]
    real(RP), intent(in)  :: RHOQ_t_MP_SNOW(KA,  IA,JA) !< tendency rho*QS (MP)                             [kg/m3/s]
    real(RP), intent(in)  :: RHOQ_t_CP_PREC(KA,  IA,JA) !< tendency rho*(QR+QS) (CP)                        [kg/m3/s]
    real(RP), intent(in)  :: flux_qr_cp    (0:KA,IA,JA) !< rain flux calcurated by cumulus parameterization [kg/m2/s]
    real(RP), intent(in)  :: flux_qs_cp    (0:KA,IA,JA) !< snow flux calcurated by cumulus parameterization [kg/m2/s]
    real(RP), intent(in)  :: cldfrac_cp    (KA,  IA,JA) !< cloud fraction (deep convection) (0-1)           [1]
    real(RP), intent(in)  :: flux_qr_mp    (0:KA,IA,JA) !< rain flux calcurated by microphysics scheme      [kg/m2/s]
    real(RP), intent(in)  :: flux_qs_mp    (0:KA,IA,JA) !< snow flux calcurated by microphysics scheme      [kg/m2/s]
    real(RP), intent(in)  :: cldfrac_mp    (KA,  IA,JA) !< cloud fraction from microphysics scheme [0-1]    [1]
    character(len=H_SHORT), intent(in) :: ATMOS_PHY_MP_TYPE !< Cloud Microphysics scheme
    real(RP), intent(in)  :: NUMDEN_MP_CLDW(KA,  IA,JA) !< number density (MP;cloud water) [num/m3]
    real(RP), intent(in)  :: NUMDEN_MP_RAIN(KA,  IA,JA) !< number density (MP;rain       ) [num/m3]
    real(RP), intent(in)  :: NUMDEN_MP_CLDI(KA,  IA,JA) !< number density (MP;cloud ice  ) [num/m3]
    real(RP), intent(in)  :: NUMDEN_MP_SNOW(KA,  IA,JA) !< number density (MP;snow       ) [num/m3]
    real(RP), intent(in)  :: dtau_s        (KA,  IA,JA) !< 0.67 micron cloud optical depth                  [1]
    
    !2D
    real(RP), intent(in)  :: GDPS            (IA,JA) !< surface pressure Ps             [Pa]
    real(RP), intent(in)  :: T2              (IA,JA) !<  2m T                           [K]
    real(RP), intent(in)  :: U10             (IA,JA) !< 10m U                           [m/s]
    real(RP), intent(in)  :: V10             (IA,JA) !< 10m V                           [m/s]
    real(RP), intent(in)  :: U10L            (IA,JA) !< 10m U over land from land model [m/s]
    real(RP), intent(in)  :: V10L            (IA,JA) !< 10m V over land from land model [m/s]
    real(RP), intent(in)  :: USTAR           (IA,JA) !< friction velocity               [m/s]
    real(RP), intent(in)  :: WSTAR           (IA,JA) !< convective velocity scale       [m/s]
    real(RP), intent(in)  :: QSTAR           (IA,JA) !< moisture scale                  [kg/kg]
    real(RP), intent(in)  :: SFC_TEMP        (IA,JA) !< surface skin temperature        [K]
    real(RP), intent(in)  :: LAND_WATER      (IA,JA) !< moisture of top soil layer      [m3/m3]
    real(RP), intent(in)  :: LAND_PROPERTY   (IA,JA) !< maximum soil moisture           [m3/m3]
    real(RP), intent(in)  :: OCEAN_ICE_FRAC  (IA,JA) !< area fraction of sea ice        [1]
    real(RP), intent(in)  :: SFCFLX_SW_dn_net(IA,JA) !< net downward short wave flux    [W/m2]
    real(RP), intent(in)  :: RCOSZ           (IA,JA) !< cos(solar zenith angle) (0-1)   [1]

    integer,  intent(in)  :: QA_AE                   !< number of aerosol tracers       [1]
    real(DP), intent(in)  :: dt_AE                   !< time step (aerosol)             [s]
    
    character(len=H_SHORT), intent(in)  :: TRCRNAME (QA_AE) !< tracer name

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE)   !< Ratio of aerosol mass to total mass [kg/kg]

    ! output
    !-------------------------
    real(RP) :: OUTQLD(KA,IA,JA,QA_AE)             !< aerosol for radiation              []
    real(RP), intent(out) :: UAPCL (KA,IA,JA)      !< cloud condensation nuclei          [1/m3]
    real(RP), intent(out) :: NUMAET(KA,IA,JA)      !< total aerosol number concentration [1/m3]
    real(RP), intent(out) :: AE_Re (KA,IA,JA,N_AE) !< effective radius for radiation     [cm]
    real(RP), intent(out) :: AE_Qe (KA,IA,JA,N_AE) !< mass ratio for radiation           [kg/kg]

    ! internal work
    !-------------------------
    !set parameter
    logical, save :: OFIRST = .true.
    real(RP) :: AKAPPA                    !< Rair/CP = 0.286
    real(RP) :: GDTV         (KA,  IA,JA) !< virtual temperature Tv
    real(RP) :: GDTVM        (0:KA,IA,JA) !< virtual temperature (half level)
    real(RP) :: GDPM         (0:KA,IA,JA) !< pressure (half level)
    real(RP) :: QSAT         (KA)         !< saturated water vapor mixing ratio [kg/kg]
    real(RP) :: RH           (KA,  IA,JA) !< relative humidity (liquid&ice)
    real(RP) :: DELP         (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP) :: DELZ         (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
    real(RP) :: DENSP        (KA,  IA,JA) !< mass of air [kg/m2]
    real(RP) :: GRWG              (IA,JA) !< subsurf. moist. saturat.[0-1]
    real(RP) :: CDVE              (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
    real(RP) :: SFC_PSAT                  !< saturatad water vapor pressure [Pa] for CDVE
    real(RP) :: SFC_QSAT                  !< for CDVE
    real(RP) :: deltaQV                   !< for CDVE
    real(RP) :: FRAIN        (0:KA,IA,JA) !< precipitation flux:cum+LSC
    real(RP) :: FRAIN_CP     (0:KA,IA,JA) !< precipitation flux:cum
    real(RP) :: FRAIN_MP_RAIN(0:KA,IA,JA) !< rain flux:LSC
    real(RP) :: FRAIN_MP_SNOW(0:KA,IA,JA) !< snow flux:LSC
    real(RP) :: VTR_MP_RAIN  (KA,  IA,JA)
    real(RP) :: VTR_MP_SNOW  (KA,  IA,JA)
    real(RP) :: RADR_MP_RAIN (KA,  IA,JA) !< effective radius of rain (microphysics scheme)    [m]
    real(RP) :: RADR_MP_SNOW (KA,  IA,JA) !< effective radius of snow (microphysics scheme)    [m]
    real(RP) :: RAINN_CP     (KA,  IA,JA) !< rain number concentration (CP)                    [1/m3]
    real(RP) :: RAINN_MP_RAIN(KA,  IA,JA) !< rain number concentration (MP;RAIN)               [1/m3]
    real(RP) :: RAINN_MP_SNOW(KA,  IA,JA) !< rain number concentration (MP;SNOW)               [1/m3]
    real(RP) :: TCLDF        (KA,  IA,JA) !< cloud fraction (total)
    real(RP) :: TCLDF2       (KA,  IA,JA) !< cloud fraction (total) for CCOVMR
    real(RP) :: CCOVMR            (IA,JA) !< cloud cover (max-ran)
    real(RP) :: TCLDW        (KA,  IA,JA) !< cloud water & ice (QC+QI)
    real(RP) :: FLIQ         (KA,  IA,JA) !< liquid fraction
    real(RP) :: DENSC             (IA,JA) !< vertically integrated mass of air
    real(RP) :: COLWV             (IA,JA) !< water vapor column loading
    ! real(RP) :: U10O              (IA,JA) !< 10m u over ocean
    ! real(RP) :: V10O              (IA,JA) !< 10m v over ocean
    ! real(RP) :: VABSL             (IA,JA) !< wind at 10m (land      ) [m/s]
    ! real(RP) :: VABSO             (IA,JA) !< wind at 10m (     ocean) [m/s]
    real(RP) :: VABST             (IA,JA) !< wind at 10m (land+ocean) [m/s]
    real(RP) :: GRLAI             (IA,JA) !< LAI [m2/m2]
    real(RP) :: GRSNW             (IA,JA) !< snow depth [m]
    real(RP) :: GOICE             (IA,JA) !< ocean ice fraction [1]
    real(RP) :: FACT_OCEAN1       (IA,JA) !< fact_ocean - fact_lake

    integer  :: KBB               (IA,JA) !< emission level of biomass burning
    real(RP) :: DPKBB             (IA,JA) !< uplift level total delta P for biomass burning [Pa]
    integer  :: KUP               (IA,JA) !< uplift level
    integer  :: KUPMAX                    !< maximum uplift level 
    real(RP) :: DPKUP             (IA,JA) !< uplift level total delta P [Pa]
    real(RP) :: RKUP              (IA,JA) !< uplift level
    
    real(RP) :: WSTE                      !< use for instability
    real(RP) :: WSATSE                    !< use for instability
    real(RP) :: UPSATQ                    !< use for instability
    real(RP) :: SATQ                      !< use for instability
    real(RP) :: UPT                       !< use for instability
    logical  :: ONEWQ        (KA-1,IA,JA) !< use for instability
    
    real(RP) :: NWRAIN                    !< mixing ratio of new rain formation
    real(RP) :: NWRAIN_CP                 !< mixing ratio of new rain formation (CP)
    real(RP) :: NWRAIN_MP_RAIN            !< mixing ratio of new rain formation (MP:rain)
    real(RP) :: NWRAIN_MP_SNOW            !< mixing ratio of new rain formation (MP:snow)
    real(RP) :: RNWR          (KA, IA,JA) !< ratio of rain to cloud
    real(RP) :: RNWR_CP       (KA, IA,JA) !< ratio of rain to cloud (CP)
    real(RP) :: RNWR_MP_RAIN  (KA, IA,JA) !< ratio of rain to cloud (MP:rain)
    real(RP) :: RNWR_MP_SNOW  (KA, IA,JA) !< ratio of rain to cloud (MP:snow)

    real(RP) :: OHRAD        (KA,  IA,JA)  !< OH concentration                            [molecules/cm3]
    real(RP) :: OHRDC        (KA,  IA,JA)  !< OH concentration = OHRAD_DATA / DAYTIM_DATA [molecules/cm3]   
    real(RP) :: NO3G         (KA,  IA,JA)  !< NO3(gas) mixing ratio:                      [pptv]
    real(RP) :: O3ORI        (KA,  IA,JA)  !< O3       mixing ratio:                      [ppmv]
    real(RP) :: DAYTIM            (IA,JA)  !< daytime ratio                               [1]
    real(RP) :: AIRBCT      (IA,JA,AIRLEV) !< BC  emission from aircrafts                 []
    real(RP) :: AIROCT      (IA,JA,AIRLEV) !< OC  emission from aircrafts                 []
    real(RP) :: AIRSO2T     (IA,JA,AIRLEV) !< SO2 emission from aircrafts                 []
    real(RP) :: AIRFEL   (KA,IA,JA)        !< fuel consumption by aircrafts               []
    real(RP) :: AIRFELOC (KA,IA,JA)        !< fuel consumption by aircrafts               []
    real(RP) :: AIRFELSO2(KA,IA,JA)        !< fuel consumption by aircrafts               []

    real(RP) :: NUMAE    (KA,IA,JA,6)      !< aerosol number density (1->dust, 2-4->carbon, 5->sulfate, 6->salt)
    real(RP) :: RADAE    (KA,IA,JA,6)      !< arosol mode radius

    real(RP) :: QLTSA(KA,IA,JA)            !< total sea salt aerosol [kg/kg]
    
    !optical parameter
    real(RP) :: TAUT  (IA,JA,IWA,2)    !< AOT (total;   total     )
    real(RP) :: TAUTA (IA,JA,IWA,2)    !< AOT (total;   absorption)
    real(RP) :: TAUTF (IA,JA,    2)    !< AOT (total;   fine mode )
    real(RP) :: TAUTD (IA,JA)          !< AOT (total;   dry       )
    real(RP) :: TAUSU (IA,JA,IWA,2)    !< AOT (sulfate; total     )
    real(RP) :: TAUSUD(IA,JA)          !< AOT (sulfate; dry       )
    real(RP) :: TAUCA (IA,JA,IWA,2)    !< AOT (carbon;  total     )
    real(RP) :: TAUCAA(IA,JA,IWA,2)    !< AOT (carbon;  absorption)
    real(RP) :: TAUCAD(IA,JA)          !< AOT (carbon;  dry       )
    real(RP) :: TAUDU (IA,JA,IWA,2)    !< AOT (dust;    total     )
    real(RP) :: TAUDUA(IA,JA,IWA,2)    !< AOT (dust;    absorption)
    real(RP) :: TAUDUF(IA,JA,    2)    !< AOT (dust;    fine mode )
    real(RP) :: TAUNN (IA,JA,IWA,2)    !< AOT (salt;    total     )
    real(RP) :: TAUNNF(IA,JA,    2)    !< AOT (salt;    fine mode )
    real(RP) :: TAUNND(IA,JA)          !< AOT (salt;    dry       )

    real(RP) :: CEXTT (KA,IA,JA,IWA,2) !< extinction  coefficient (total;   total    )
    real(RP) :: CABST (KA,IA,JA,IWA,2) !< absorption  coefficient (total;   total    )
    real(RP) :: CBAKT (KA,IA,JA,IWA,2) !< backscatter coefficient (total;   total    )
    real(RP) :: CEXTTF(KA,IA,JA,    2) !< extinction  coefficient (total;   fine mode)
    real(RP) :: CEXTTD(KA,IA,JA,    2) !< extinction  coefficient (total;   dry      )
    real(RP) :: CEXTSU(KA,IA,JA,IWA,2) !< extinction  coefficient (sulfate; total    )
    real(RP) :: CBAKSU(KA,IA,JA,IWA,2) !< backscatter coefficient (sulfate; total    )
    real(RP) :: CEXTSD(KA,IA,JA)       !< extinction  coefficient (sulfate; dry      )
    real(RP) :: CEXTCA(KA,IA,JA,IWA,2) !< extinction  coefficient (carbon;  total    )
    real(RP) :: CABSCA(KA,IA,JA,IWA,2) !< absorption  coefficient (carbon;  total    )
    real(RP) :: CBAKCA(KA,IA,JA,IWA,2) !< backscatter coefficient (carbon;  total    )
    real(RP) :: CEXTCD(KA,IA,JA)       !< extinction  coefficient (carbon;  dry      )
    real(RP) :: CEXTDU(KA,IA,JA,IWA,2) !< extinction  coefficient (dust;    total    )
    real(RP) :: CABSDU(KA,IA,JA,IWA,2) !< absorption  coefficient (dust;    total    )
    real(RP) :: CBAKDU(KA,IA,JA,IWA,2) !< backscatter coefficient (dust;    total    )
    real(RP) :: CEXTDF(KA,IA,JA,    2) !< extinction  coefficient (dust;    fine mode)
    real(RP) :: CEXTNN(KA,IA,JA,IWA,2) !< extinction  coefficient (salt;    total    )
    real(RP) :: CBAKNN(KA,IA,JA,IWA,2) !< backscatter coefficient (salt;    total    )
    real(RP) :: CEXTNF(KA,IA,JA,    2) !< extinction  coefficient (salt;    fine mode)
    real(RP) :: CEXTND(KA,IA,JA)       !< extinction  coefficient (salt;    dry      )

    real(RP) :: ALFA     (IA,JA,2)     !< Angstrom exponent
    real(RP) :: SSA      (IA,JA,2)     !< single scattering albedo

    real(RP) :: PM25  (KA,IA,JA)       !< submicron mass concentration
    real(RP) :: PM25CA(KA,IA,JA)       !< submicron mass concentration (carbon)
    real(RP) :: PM25SU(KA,IA,JA)       !< submicron mass concentration (sulfate)
    real(RP) :: PM25DU(KA,IA,JA)       !< submicron mass concentration (dust)
    real(RP) :: PM25SA(KA,IA,JA)       !< submicron mass concentration (salt)
    real(RP) :: QPM25 (KA,IA,JA)       !< submicron mass mixing ratio

    real(RP) :: WATRAE(KA,IA,JA)       !< aerosol water (total)
    real(RP) :: WATROC(KA,IA,JA)       !< aerosol water (OC)
    real(RP) :: WATRSU(KA,IA,JA)       !< aerosol water (sulfate)
    real(RP) :: WATRSA(KA,IA,JA)       !< aerosol water (salt)

    real(RP) :: QCQR_mass(KA,IA,JA)
    ! others
    !-------------------------
    integer  :: k, i, j, iq, iw
    real(RP) :: index_pft_real(IA,JA,2)
    
    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / aerosol / SPRINTARS - admin'
    LOG_PROGRESS(*) 'atmosphere / physics / aerosol / SPRINTARS '
    
    !SPRINTARS is not used at first step
    !================================================
    if ( OFIRST ) then
       OFIRST = .false.
       UAPCL (:,:,:)   = 0.0_RP
       NUMAET(:,:,:)   = 0.0_RP
       AE_Re (:,:,:,1) = 2.0000E-6_RP * 1.0E+2_RP ! soil dust                [m] -> [cm]
       AE_Re (:,:,:,2) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.30 [m] -> [cm]
       AE_Re (:,:,:,3) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.15 [m] -> [cm]
       AE_Re (:,:,:,4) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.00 [m] -> [cm]
       AE_Re (:,:,:,5) = 0.0118E-6_RP * 1.0E+2_RP ! ext. BC                  [m] -> [cm]
       AE_Re (:,:,:,6) = 0.0695E-6_RP * 1.0E+2_RP ! sulfate                  [m] -> [cm]
       AE_Re (:,:,:,7) = 2.0000E-6_RP * 1.0E+2_RP ! sea salt                 [m] -> [cm]
       AE_Qe (:,:,:,:) = 0.0_RP                   ! [kg/kg]
       FRAIN_MP_RAIN1(0:KA,:,:) = 0.0_RP
       FRAIN_MP_SNOW1(0:KA,:,:) = 0.0_RP
       DENS1         (:,   :,:) = DENS(:,:,:)
       !negative fixer
       do iq = 1, QA_AE, 1
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                if ( QTRC(k,i,j,iq) < 0.0_RP ) then
                   write(*,'("Negative value exists! fixed. subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL(ONCE) ", a10, " [myrank,iq,k,i,j]=[", 5i5, 2x, e12.5, "]")') & 
                        trim(adjustl(TRCRNAME(iq))), PRC_myrank, iq, k, i, j, QTRC(k,i,j,iq)
                   QTRC(k,i,j,iq) = 1.0E-20_RP
                endif
             enddo
          enddo
          enddo
       enddo
       return
    endif
    !================================================
    !end first step

        
    !setup parameter
    !================================================
    !kappa
    AKAPPA = RAIR / CP       ! Rair/CP = 0.286

    do j = JS, JE
    do i = IS, IE

       !GDTV (virtual temperature)
       do k = KS, KE
          GDTV(k,i,j) = GDT(k,i,j) * ( 1.0_RP + EPSTV * QV(k,i,j) )
       enddo

       !GDTVM (virtual temperature at half level)
       do k = KS, KE-1
          GDTVM(k,i,j) = ( GDTV(k+1,i,j)-GDTV(k,i,j) ) / ( GDZ(k+1,i,j)-GDZ(k,i,j) ) &
                       * ( GDZM(k,i,j)-GDZ(k,i,j) ) + GDTV(k,i,j)
       enddo
       GDTVM(KE,  i,j) = GDT(KE,i,j)
       GDTVM(KS-1,i,j) = GDT(KS,i,j)

       !GDPM (pressure at half level)
       GDPM(0:KS-2,i,j) = CONST_UNDEF
       if ( GDPS(i,j) == CONST_UNDEF ) then
          GDPM(KS-1,i,j) = exp( (log(GDP(KS,i,j))-log(GDP(KS+1,i,j))) &
                                  / (GDZ(KS,i,j) -    GDZ(KS+1,i,j) ) * (GDZM(KS-1,i,j)-GDZ(KS+1,i,j)) + log(GDP(KS+1,i,j)) )
       else
          GDPM(KS-1,i,j) = GDPS(i,j)
       endif
       do k = KS, KE-1
          ! GDPM(k,i,j) = exp( (log(GDP(k,i,j))-log(GDP(k+1,i,j))) &
          !                      / (GDZ(k,i,j) -    GDZ(k+1,i,j) ) * (GDZM(k,i,j)-GDZ(k+1,i,j)) + log(GDP(k+1,i,j)) )
          GDPM(k,i,j) = ( (GDP(k,i,j)**(AKAPPA+1.0_RP)-GDP(k+1,i,j)**(AKAPPA+1.0_RP))/(GDP(k,i,j)-GDP(k+1,i,j))/(1.0_RP+AKAPPA) )**(1.0_RP/AKAPPA)
       enddo
       if ( GDP(KE+1,i,j) == CONST_UNDEF ) then
          GDPM(KE,i,j) = exp( (log(GDP(KE-1,i,j))-log(GDP(KE,i,j))) &
                                / (GDZ(KE-1,i,j) -    GDZ(KE,i,j) ) * (GDZM(KE,i,j)-GDZ(KE,i,j)) + log(GDP(KE,i,j)) )
       else
          ! GDPM(KE,i,j) = exp( (log(GDP(KE,i,j))-log(GDP(KE+1,i,j))) &
          !                       / (GDZ(KE,i,j) -    GDZ(KE+1,i,j) ) * (GDZM(KE,i,j)-GDZ(KE+1,i,j)) + log(GDP(KE+1,i,j)) )
          GDPM(KE,i,j) = ( (GDP(KE,i,j)**(AKAPPA+1.0_RP)-GDP(KE+1,i,j)**(AKAPPA+1.0_RP))/(GDP(KE,i,j)-GDP(KE+1,i,j))/(1.0_RP+AKAPPA) )**(1.0_RP/AKAPPA)
       endif
       GDPM(KE+1:KA,i,j) = CONST_UNDEF

       !RH (relative humidity (liquid & ice)) [kg/kg]
       call ATMOS_SATURATION_dens2qsat_all( KA, KS, KE, GDT(:,i,j), DENS(:,i,j), QSAT(:) )
       do k = KS, KE
          RH(k,i,j) = QV(k,i,j) / QSAT(k)
          RH(k,i,j) = max( min( RH(k,i,j), 1.0_RP ), 0.0_RP )
       enddo

       !DELP (delta pressure), DELZ (delta z), DENSP (mass of air)
       do k = KS, KE
          DELP(k,i,j)  = GDPM(k-1,i,j) - GDPM(k,  i,j)
          DELZ(k,i,j)  = GDZM(k,  i,j) - GDZM(k-1,i,j)
          DENSP(k,i,j) = DENS(k,i,j) * DELZ(k,i,j)
       enddo

       !GRWG (subsurf. moist. saturat.[0-1])
       GRWG(i,j) = LAND_WATER(i,j) / LAND_PROPERTY(i,j)
       GRWG(i,j) = max( min( GRWG(i,j), 1.0_RP ), 0.0_RP )
       
       !CDVE (bulk coefficient CeU) (cf. scale_atmos_phy_sf_bulk.F90)
       call ATMOS_SATURATION_psat_all( SFC_TEMP(i,j), SFC_PSAT )
       SFC_QSAT = EPSvap * SFC_PSAT / ( GDPS(i,j) - ( 1.0_RP-EPSvap ) * SFC_PSAT )
       deltaQV = sign( max( abs(QV(KS,i,j)-SFC_QSAT), EPS ), QV(KS,i,j)-SFC_QSAT )
       CDVE(i,j) = max( min( USTAR(i,j) * QSTAR(i,j) / deltaQV, CDVEMAX ), CDVEMIN )

       !FRAIN (precipitation flux : rain (cp) + snow(cp) + rain(mp) + snow(mp))
       do k = KS-1, KE-1
          FRAIN_CP     (k,i,j) = max( flux_qr_cp(k,i,j), 0.0_RP ) + max( flux_qs_cp(k,i,j), 0.0_RP )
          FRAIN_MP_RAIN(k,i,j) = max( flux_qr_mp(k,i,j), 0.0_RP )
          FRAIN_MP_SNOW(k,i,j) = max( flux_qs_mp(k,i,j), 0.0_RP )
          FRAIN        (k,i,j) = FRAIN_CP(k,i,j) + FRAIN_MP_RAIN(k,i,j) + FRAIN_MP_SNOW(k,i,j)
       enddo
       FRAIN        (0 :KS-2,i,j) = 0.0_RP
       FRAIN        (KE:KA,  i,j) = 0.0_RP
       FRAIN_CP     (0 :KS-2,i,j) = 0.0_RP
       FRAIN_CP     (KE:KA,  i,j) = 0.0_RP
       FRAIN_MP_RAIN(0 :KS-2,i,j) = 0.0_RP
       FRAIN_MP_RAIN(KE:KA,  i,j) = 0.0_RP
       FRAIN_MP_SNOW(0 :KS-2,i,j) = 0.0_RP
       FRAIN_MP_SNOW(KE:KA,  i,j) = 0.0_RP
       
       !VTR_MP_RAIN & VTR_MP_SNOW (terminal velocity : RAIN & SNOW (microphysics scheme))
       do k = KS, KE
          VTR_MP_RAIN(k,i,j) = max( -VTERM_MP_RAIN(k,i,j), THRESVTR )
          VTR_MP_SNOW(k,i,j) = max( -VTERM_MP_SNOW(k,i,j), THRESVTR )
       enddo

       !RADR_MP_RAIN, RADR_MP_SNOW
       do k = KS, KE
          RADR_MP_RAIN(k,i,j) = max( RAIN_MP_Re(k,i,j), THRESRe )
          RADR_MP_SNOW(k,i,j) = max( SNOW_MP_Re(k,i,j), THRESRe )
       enddo

       !RAINN_CP, RAINN_MP_RAIN, RAINN_MP_SNOW
       if ( ATMOS_PHY_MP_TYPE == 'SN14' ) then
          do k = KS, KE
             RAINN_CP     (k,i,j) = FRAIN_CP      (k,i,j) / VTR                / ( DENSR * 4.0_RP / 3.0_RP * PI * RADR**3                )
             RAINN_MP_RAIN(k,i,j) = NUMDEN_MP_RAIN(k,i,j)
             RAINN_MP_SNOW(k,i,j) = NUMDEN_MP_SNOW(k,i,j)
          enddo
       else
          do k = KS, KE
             RAINN_CP     (k,i,j) = FRAIN_CP     (k,i,j) / VTR                / ( DENSR * 4.0_RP / 3.0_RP * PI * RADR**3                )
             RAINN_MP_RAIN(k,i,j) = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / ( DENSR * 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
             RAINN_MP_SNOW(k,i,j) = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / ( DENSI * 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
          enddo
       endif
       
       !TCLDF (cloud fraction (total))
       do k = KS, KE
          QCQR_mass(k,i,j) = ( QC(k,i,j)+QI(k,i,j) )*DENS(k,i,j) 
          if( QCQR_mass(k,i,j) > WDEP_CLD .or. dtau_s(k,i,j) > WDEP_TAU ) then
             !TCLDF(k,i,j) = 1.0_RP - ( 1.0_RP - cldfrac_mp(k,i,j) ) * ( 1.0_RP - cldfrac_cp(k,i,j) )
             !TCLDF(k,i,j) = 1.0_RP - ( 1.0_RP - 1.0_RP ) * ( 1.0_RP - cldfrac_cp(k,i,j) )
             TCLDF(k,i,j) = 1.0_RP 
          else
             !TCLDF(k,i,j) = 1.0_RP - ( 1.0_RP - 0.0_RP ) * ( 1.0_RP - cldfrac_cp(k,i,j) )
             TCLDF(k,i,j) = 0.0_RP
          endif
          TCLDF(k,i,j) = min( max( TCLDF(k,i,j), 0.0_RP ), 1.0_RP )
       enddo

       !TCLDF2 (cloud fraction (total) for CCOVMR) (Goto et al., 2020, NICAM)
       do k = KS, KE
          if ( ( QC(k,i,j)+QI(k,i,j) ) * DENS(k,i,j) >= 1.0E-3_RP .or. dtau_s(k,i,j) > 0.2_RP ) then
             TCLDF2(k,i,j) = 1.0_RP - ( 1.0_RP - 1.0_RP ) * ( 1.0_RP - cldfrac_cp(k,i,j) )
          else
             TCLDF2(k,i,j) = 1.0_RP - ( 1.0_RP - 0.0_RP ) * ( 1.0_RP - cldfrac_cp(k,i,j) )
          endif
          TCLDF2(k,i,j) = min( max( TCLDF2(k,i,j), 0.0_RP ), 1.0_RP )
       enddo

       !CCOVMR
       !maximum-random over lapping (Morcrette's formula) : modified MIROC nonstd/pradmX.F (using TCLDF)
       CCOVMR(i,j) = 1.0_RP
       do k = KE, KS, -1
          if ( k == KE ) then
             CCOVMR(i,j) = CCOVMR(i,j) * ( 1.0_RP - TCLDF2(k,i,j) )
          else
             if ( TCLDF2(k+1,i,j) >= 1.0_RP ) then
                CCOVMR(i,j) = CCOVMR(i,j)
             else
                CCOVMR(i,j) = CCOVMR(i,j) * ( 1.0_RP - max(TCLDF2(k,i,j),TCLDF2(k+1,i,j)) ) / ( 1.0_RP - TCLDF2(k+1,i,j) )
             endif
          endif
       enddo
       CCOVMR(i,j) = 1.0_RP - CCOVMR(i,j)

       !TCLDW, FLIQ
       do k = KS, KE
          TCLDW(k,i,j) = QC(k,i,j) + QI(k,i,j)
          FLIQ (k,i,j) = 0.0_RP
          if ( TCLDW(k,i,j) > 0.0_RP ) then
             FLIQ(k,i,j) = QC(k,i,j) / TCLDW(k,i,j)
          endif
          FLIQ (k,i,j) = min( max( FLIQ(k,i,j), 0.0_RP ), 1.0_RP )
       enddo

       !DENSC, COLWV
       DENSC(i,j) = 0.0_RP
       COLWV(i,j) = 0.0_RP
       do k = KS, KE
          DENSC(i,j) = DENSC(i,j) + DENSP(k,i,j)
          COLWV(i,j) = COLWV(i,j) + QV(k,i,j) * DENSP(k,i,j)
       enddo
       
       !U10O, V10O, VABSL, VABSO, VABST (10m wind speed) 
       ! if ( FLND(i,j) < FLNDMX ) then
       !    U10O(i,j) = ( U10(i,j) - U10L(i,j) * FLND(i,j) ) / ( 1.0_RP - FLND(i,j) )
       !    V10O(i,j) = ( V10(i,j) - V10L(i,j) * FLND(i,j) ) / ( 1.0_RP - FLND(i,j) )
       ! else
       !    U10O(i,j) = 0.0_RP
       !    V10O(i,j) = 0.0_RP
       ! endif
       ! VABSL(i,j) = sqrt( U10L(i,j)**2 + V10L(i,j)**2 + WSTAR(i,j)**2 )
       ! VABSO(i,j) = sqrt( U10O(i,j)**2 + V10O(i,j)**2 + WSTAR(i,j)**2 )
       ! VABST(i,j) = FLND(i,j) * VABSL(i,j) + ( 1.0_RP - FLND(i,j) ) * VABSO(i,j)
       ! U10O (i,j) = 0.0_RP
       ! V10O (i,j) = 0.0_RP
       ! VABSL(i,j) = 0.0_RP
       ! VABSO(i,j) = 0.0_RP
       VABST(i,j) = sqrt( U10(i,j)**2 + V10(i,j)**2 + (min(WSTAR(i,j),4.0_RP))**2 )
       
    enddo
    enddo
 
    !GRLAI, GRSNW, GOICE
    !-------------------------------
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( GRLAI(:,:),                                               & ! [OUT]
           GRLAI_DATA(:,:,12), GRLAI_DATA(:,:,:), GRLAI_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)         ) ! [IN]
    
    do j = JS, JE
    do i = IS, IE
       GRLAI(i,j) = max( GRLAI(i,j), 0.0_RP )
       GRSNW(i,j) = max( GRSNW_DATA(i,j,TIME_NOWDATE(2),TIME_NOWDATE(3)), 0.0_RP )
       GOICE(i,j) = max( min( OCEAN_ICE_FRAC(i,j), 1.0_RP ), 0.0_RP )
    enddo
    enddo

    !FACT_OCEAN1
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       if ( FACT_LAKE_SPRINTARS(i,j) > 0.0_RP ) then
          FACT_OCEAN1(i,j) = max( min( FACT_OCEAN(i,j)-FACT_LAKE_SPRINTARS(i,j), 1.0_RP ), 0.0_RP )
       else
          FACT_OCEAN1(i,j) = FACT_OCEAN(i,j)
       endif
    enddo
    enddo
    
    !for emission : KUP, KUPMAX, DPKUP, DPKBB
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       KBB(i,j) = KS
       do k = KS, KE
          if ( GDZ(k,i,j) >= SBB + GDZM(KS-1,i,j) ) then
             KBB(i,j) = k
             exit
          endif
       enddo
    enddo
    enddo
    
    do j = JS, JE
    do i = IS, IE
       KUP(i,j) = KS
       KUPMAX = KS
       KUPMAX = max( KUPMAX, KUP(i,j) )
       DPKUP(i,j) = GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
       DPKBB(i,j) = GDPM(KS-1,i,j) - GDPM(KBB(i,j),i,j)
    enddo
    enddo
       
    ! !for instability : SATQ, WSTE, WSATSE, UPT, UPSATQ, ONEWQ
    ! !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    call ATMOS_SATURATION_dens2qsat_all( KA, KS, KE, GDT(:,i,j), DENS(:,i,j), QSAT(:) )
    !    do k = KS, KE-1
    !       WSTE   = CP * GDT(k  ,i,j) + GRAV * GDZ(k  ,i,j) + EL * QV(k,i,j)
    !       WSATSE = CP * GDT(k+1,i,j) + GRAV * GDZ(k+1,i,j) + EL * QSAT(k+1)
    !       UPT    = ( WSTE - GRAV * GDZ(k+1,i,j) - EL * QV(k,i,j) ) / CP
    !       call ATMOS_SATURATION_dens2qsat_all( UPT, DENS(k+1,i,j), UPSATQ )
    !       if (       WSATSE < WSTE      &
    !            .and. UPSATQ < QV(k,i,j) &
    !            .and. ODCA               ) then
    !          ONEWQ(k,i,j) = .true.
    !       else
    !          ONEWQ(k,i,j) = .false.
    !       endif
    !    enddo
    ! enddo
    ! enddo
    ONEWQ(:,:,:) = .false.

    !for cloud scavenging (rain out) : NWRAIN, RNWR (MPについてはQR_t,QS_tをつかったほうがいい?)
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          ! NWRAIN         = ( FRAIN        (k-1,i,j) - FRAIN        (k,i,j) ) * GRAV / DELP(k,i,j) * real(dt_AE,kind=RP)
          ! NWRAIN_CP      = ( FRAIN_CP     (k-1,i,j) - FRAIN_CP     (k,i,j) ) * GRAV / DELP(k,i,j) * real(dt_AE,kind=RP)
          ! NWRAIN_MP_RAIN = ( FRAIN_MP_RAIN(k-1,i,j) - FRAIN_MP_RAIN(k,i,j) ) * GRAV / DELP(k,i,j) * real(dt_AE,kind=RP)
          ! NWRAIN_MP_SNOW = ( FRAIN_MP_SNOW(k-1,i,j) - FRAIN_MP_SNOW(k,i,j) ) * GRAV / DELP(k,i,j) * real(dt_AE,kind=RP)
          
          NWRAIN_CP = ( FRAIN_CP(k-1,i,j) - FRAIN_CP(k,i,j) ) * GRAV / DELP(k,i,j) * real(dt_AE,kind=RP)
          ! NWRAIN_CP = RHOQ_t_CP_PREC(k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP)
          if ( GDT(k,i,j) > TMELT ) then
             NWRAIN_MP_RAIN = ( RHOQ_t_MP_RAIN(k,i,j) + RHOQ_t_MP_SNOW(k,i,j) ) / DENS(k,i,j) * real(dt_AE,kind=RP)
             NWRAIN_MP_SNOW = 0.0_RP
          else
             NWRAIN_MP_RAIN = RHOQ_t_MP_RAIN(k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP)
             NWRAIN_MP_SNOW = RHOQ_t_MP_SNOW(k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP)
          endif
          
          NWRAIN = 0.0_RP
          if ( NWRAIN_CP      > EPS ) NWRAIN = NWRAIN + NWRAIN_CP
          if ( NWRAIN_MP_RAIN > EPS ) NWRAIN = NWRAIN + NWRAIN_MP_RAIN
          if ( NWRAIN_MP_SNOW > EPS ) NWRAIN = NWRAIN + NWRAIN_MP_SNOW
          
          if ( NWRAIN > EPS ) then
             RNWR(k,i,j) = NWRAIN / ( TCLDW(k,i,j) + NWRAIN )
             if ( NWRAIN_CP      > EPS ) then
                RNWR_CP     (k,i,j) = NWRAIN_CP      / ( TCLDW(k,i,j) + NWRAIN )
             else
                RNWR_CP     (k,i,j) = 0.0_RP
             endif
             if ( NWRAIN_MP_RAIN > EPS ) then
                RNWR_MP_RAIN(k,i,j) = NWRAIN_MP_RAIN / ( TCLDW(k,i,j) + NWRAIN )
             else
                RNWR_MP_RAIN(k,i,j) = 0.0_RP
             endif
             if ( NWRAIN_MP_SNOW > EPS ) then
                RNWR_MP_SNOW(k,i,j) = NWRAIN_MP_SNOW / ( TCLDW(k,i,j) + NWRAIN )
             else
                RNWR_MP_SNOW(k,i,j) = 0.0_RP
             endif
          else
             RNWR        (k,i,j) = 0.0_RP
             RNWR_CP     (k,i,j) = 0.0_RP
             RNWR_MP_RAIN(k,i,j) = 0.0_RP
             RNWR_MP_SNOW(k,i,j) = 0.0_RP
          endif
          RNWR        (k,i,j) = max( min( RNWR        (k,i,j), 1.0_RP ), 0.0_RP )
          RNWR_CP     (k,i,j) = max( min( RNWR_CP     (k,i,j), 1.0_RP ), 0.0_RP )
          RNWR_MP_RAIN(k,i,j) = max( min( RNWR_MP_RAIN(k,i,j), 1.0_RP ), 0.0_RP )
          RNWR_MP_SNOW(k,i,j) = max( min( RNWR_MP_SNOW(k,i,j), 1.0_RP ), 0.0_RP )
          
       enddo
    enddo
    enddo

    !OHRDC, NO3G, O3ORI
    !-------------------------------
    do k = KS, KE
       call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
            ( OHRAD(k,:,:),                                                      & ! [OUT]
              OHRAD_DATA(k,:,:,12), OHRAD_DATA(k,:,:,1:12), OHRAD_DATA(k,:,:,1), & ! [IN]
              TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                  ) ! [IN]

       call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
            ( NO3G (k,:,:),                                                      & ! [OUT]
              NO3G_DATA (k,:,:,12), NO3G_DATA (k,:,:,1:12), NO3G_DATA (k,:,:,1), & ! [IN]
              TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                  ) ! [IN]

       call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
            ( O3ORI(k,:,:),                                                      & ! [OUT]
              O3ORI_DATA(k,:,:,12), O3ORI_DATA(k,:,:,1:12), O3ORI_DATA(k,:,:,1), & ! [IN]
              TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                  ) ! [IN]
    enddo

    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( DAYTIM(:,:),                                                    & ! [OUT]
           DAYTIM_DATA(:,:,12), DAYTIM_DATA(:,:,1:12), DAYTIM_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)               ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       DAYTIM(i,j) = max( DAYTIM(i,j), 0.0_RP )
       do k = KS, KE
          OHRAD(k,i,j) = max( OHRAD(k,i,j), 0.0_RP )
          NO3G (k,i,j) = max( NO3G (k,i,j), 0.0_RP )
          O3ORI(k,i,j) = max( O3ORI(k,i,j), 0.0_RP )
          if ( DAYTIM(i,j) > 0.05_RP ) then
             if ( RCOSZ(i,j) >= 0.0_RP ) then
                OHRDC(k,i,j) = OHRAD(k,i,j) / DAYTIM(i,j)
             else
                OHRDC(k,i,j) = 0.0_RP
             endif
          else
             OHRDC(k,i,j) = OHRAD(k,i,j)
          endif
       enddo
    enddo
    enddo

    !AIRBCT, AIROCT, AIRSO2T, AIRFEL, AIRFELOC, AIRFELSO2
    !-------------------------------
    
    !================================================
    !end setup parameter


    !aerosol transport
    !================================================
    
    !setup
    !-------------------------------
    OUTQLD(:,:,:,:) = 0.0_RP
    
    NUMAE(:,:,:,:) = 0.0_RP
    RADAE(:,:,:,1) = 2.0000E-6_RP ! soil dust
    RADAE(:,:,:,2) = 0.0118E-6_RP ! carbon (ext. BC)
    RADAE(:,:,:,3) = 0.1000E-6_RP ! carbon (int. OC & BC)
    RADAE(:,:,:,4) = 0.1000E-6_RP ! carbon (SOC)
    RADAE(:,:,:,5) = 0.0695E-6_RP ! sulfate
    RADAE(:,:,:,6) = 2.0000E-6_RP ! salt

    AE_Re(:,:,:,1) = 2.0000E-6_RP * 1.0E+2_RP ! soil dust                [m] -> [cm]
    AE_Re(:,:,:,2) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.30 [m] -> [cm]
    AE_Re(:,:,:,3) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.15 [m] -> [cm]
    AE_Re(:,:,:,4) = 0.1000E-6_RP * 1.0E+2_RP ! carbonaceous(BC/OC)=0.00 [m] -> [cm]
    AE_Re(:,:,:,5) = 0.0118E-6_RP * 1.0E+2_RP ! ext. BC                  [m] -> [cm]
    AE_Re(:,:,:,6) = 0.0695E-6_RP * 1.0E+2_RP ! sulfate                  [m] -> [cm]
    AE_Re(:,:,:,7) = 2.0000E-6_RP * 1.0E+2_RP ! sea salt                 [m] -> [cm]
    AE_Qe(:,:,:,:) = 0.0_RP                   ! [kg/kg]


    !soil dust
    !-------------------------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_DU ) then
       call ATMOS_PHY_AE_SPRINTARS_AERODU &
            ( QTRC(:,:,:,QS_AE_DU:QE_AE_DU), & ! [MODIFIED]
              OUTQLTDU(:,:,:,:),             & ! [OUT]
              NUMAE(:,:,:,1),                & ! [OUT]
              RADAE(:,:,:,1),                & ! [OUT]
              PM25DU(:,:,:),                 & ! [OUT]
              TAUDU (:,:,:,:),               & ! [OUT]
              TAUDUA(:,:,:,:),               & ! [OUT]
              TAUDUF(:,:,  :),               & ! [OUT]
              CEXTDU(:,:,:,:,:),             & ! [OUT]
              CABSDU(:,:,:,:,:),             & ! [OUT]
              CBAKDU(:,:,:,:,:),             & ! [OUT]
              CEXTDF(:,:,:,  :),             & ! [OUT]
              CDVE(:,:),                     & ! [IN] 
              CCOVMR(:,:),                   & ! [IN] OK
              CCMAX,                         & ! [IN] OK
              VABST(:,:),                    & ! [IN] OK
              GRWG(:,:),                     & ! [IN] check OK
              GRSNW(:,:),                    & ! [IN] fix
              GRLAI(:,:),                    & ! [IN] fix
              DENS(:,:,:),                   & ! [IN] check OK
              DENS1(:,:,:),                  & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:),      & ! [IN]
              FRAIN_MP_SNOW1(0:KA,:,:),      & ! [IN]
              FRAIN(0:KA,:,:),               & ! [IN] OK
              FRAIN_CP(0:KA,:,:),            & ! [IN]
              FRAIN_MP_RAIN(0:KA,:,:),       & ! [IN]
              FRAIN_MP_SNOW(0:KA,:,:),       & ! [IN]
              VTR_MP_RAIN(:,:,:),            & ! [IN]
              VTR_MP_SNOW(:,:,:),            & ! [IN]
              RADR_MP_RAIN(:,:,:),           & ! [IN]
              RADR_MP_SNOW(:,:,:),           & ! [IN]
              RAINN_CP(:,:,:),               & ! [IN]
              RAINN_MP_RAIN(:,:,:),          & ! [IN]
              RAINN_MP_SNOW(:,:,:),          & ! [IN]
              TCLDF(:,:,:),                  & ! [IN] OK
              TCLDF2(:,:,:),                 & ! [IN]
              KUP(:,:),                      & ! [IN] OK
              KUPMAX,                        & ! [IN] OK
              DPKUP(:,:),                    & ! [IN] OK
              GDPM(0:KA,:,:),                & ! [IN] check OK
              DELP(:,:,:),                   & ! [IN] OK
              GDZM(0:KA,:,:),                & ! [IN] OK
              DELZ(:,:,:),                   & ! [IN] OK
              ONEWQ(1:KA-1,:,:),             & ! [IN] OK
              RNWR(:,:,:),                   & ! [IN] OK
              RNWR_CP(:,:,:),                & ! [IN] OK
              RNWR_MP_RAIN(:,:,:),           & ! [IN] OK
              RNWR_MP_SNOW(:,:,:),           & ! [IN] OK
              LON_REAL(:,:),                 & ! [IN] OK
              LAT_REAL(:,:),                 & ! [IN] OK
              INDEX_PFT(:,:,:),              & ! [IN] OK
              FRAC_PFT(:,:,:),               & ! [IN] OK
              FACT_LAND(:,:),                & ! [IN] OK
              GDZS(:,:),                     & ! [IN] OK
              THRE_HGT_DUST_Tibet,           & ! [IN] 
              ACTI(:,:,:,1),                 & ! [IN]
              ATMOS_SW_PHY_AE_SPRINTARS_CCN, & ! [IN] OK
              dt_AE,                         & ! [IN] OK
              QA_AE_DU                       ) ! [IN] OK
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             AE_Re(k,i,j,1) = RADAE(k,i,j,1) * 1.0E+2_RP           ! soil dust [m] -> [cm]
             AE_Qe(k,i,j,1) = sum( QTRC(k,i,j,QS_AE_DU:QE_AE_DU) ) ! soil dust [kg/kg]
          enddo
       enddo
       enddo
       
    else
       !      k i j w 2
       PM25DU(:,:,:)     = 0.0_RP
       TAUDU   (:,:,:,:) = 0.0_RP
       TAUDUA  (:,:,:,:) = 0.0_RP
       TAUDUF  (:,:,  :) = 0.0_RP
       CEXTDU(:,:,:,:,:) = 0.0_RP
       CABSDU(:,:,:,:,:) = 0.0_RP
       CBAKDU(:,:,:,:,:) = 0.0_RP
       CEXTDF(:,:,:,  :) = 0.0_RP
    endif

    
    !sea salt
    !-------------------------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SA ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROSA &
            ( QTRC(:,:,:,QS_AE_SA:QE_AE_SA), & ! [MODIFIED]
              OUTQLTSA(:,:,:,:),             & ! [OUT]
              NUMAE(:,:,:,6),                & ! [OUT]
              RADAE(:,:,:,6),                & ! [OUT]
              PM25SA(:,:,:),                 & ! [OUT]
              WATRSA(:,:,:),                 & ! [OUT] 
              EMITSA(:,:,1:QA_AE_SA),        & ! [OUT] 
              TAUNN (:,:,:,:),               & ! [OUT]
              TAUNNF(:,:,  :),               & ! [OUT]
              CEXTNN(:,:,:,:,:),             & ! [OUT]
              CBAKNN(:,:,:,:,:),             & ! [OUT]
              CEXTNF(:,:,:,  :),             & ! [OUT]
              TAUNND(:,:),                   & ! [OUT]
              CEXTND(:,:,:),                 & ! [OUT]
              CDVE(:,:),                     & ! [IN]
              CCOVMR(:,:),                   & ! [IN] OK
              CCMAX,                         & ! [IN] OK
              VABST(:,:),                    & ! [IN] OK
              GOICE(:,:),                    & ! [IN] OK
              RH(:,:,:),                     & ! [IN] OK
              DENS(:,:,:),                   & ! [IN] OK
              DENS1(:,:,:),                  & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:),      & ! [IN]
              FRAIN_MP_SNOW1(0:KA,:,:),      & ! [IN]
              FRAIN(0:KA,:,:),               & ! [IN] OK
              FRAIN_CP(0:KA,:,:),            & ! [IN]
              FRAIN_MP_RAIN(0:KA,:,:),       & ! [IN]
              FRAIN_MP_SNOW(0:KA,:,:),       & ! [IN]
              VTR_MP_RAIN(:,:,:),            & ! [IN]
              VTR_MP_SNOW(:,:,:),            & ! [IN]
              RADR_MP_RAIN(:,:,:),           & ! [IN]
              RADR_MP_SNOW(:,:,:),           & ! [IN]
              RAINN_CP(:,:,:),               & ! [IN]
              RAINN_MP_RAIN(:,:,:),          & ! [IN]
              RAINN_MP_SNOW(:,:,:),          & ! [IN]
              TCLDF(:,:,:),                  & ! [IN] OK
              TCLDF2(:,:,:),                 & ! [IN]
              KUP(:,:),                      & ! [IN] OK
              KUPMAX,                        & ! [IN] OK
              DPKUP(:,:),                    & ! [IN] OK
              GDPM(0:KA,:,:),                & ! [IN] OK
              DELP(:,:,:),                   & ! [IN] OK
              GDZM(0:KA,:,:),                & ! [IN] OK
              DELZ(:,:,:),                   & ! [IN] OK
              ONEWQ(1:KA-1,:,:),             & ! [IN] OK
              RNWR(:,:,:),                   & ! [IN] OK
              RNWR_CP(:,:,:),                & ! [IN] OK
              RNWR_MP_RAIN(:,:,:),           & ! [IN] OK
              RNWR_MP_SNOW(:,:,:),           & ! [IN] OK
              FACT_OCEAN1(:,:),              & ! [IN] OK
              ACTI(:,:,:,6),                 & ! [IN]
              ATMOS_SW_PHY_AE_SPRINTARS_CCN, & ! [IN] OK
              dt_AE,                         & ! [IN] OK
              QA_AE_SA                       ) ! [IN] OK
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             AE_Re(k,i,j,7) = RADAE(k,i,j,6) * 1.0E+2_RP           ! sea salt [m] -> [cm]
             AE_Qe(k,i,j,7) = sum( QTRC(k,i,j,QS_AE_SA:QE_AE_SA) ) ! sea salt [kg/kg]
          enddo
       enddo
       enddo

       QLTSA(:,:,:) = 0.0_RP
       do iq = 1, QA_AE_SA
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QLTSA(k,i,j) = QLTSA(k,i,j) + OUTQLD(k,i,j,iq)
             enddo
          enddo
          enddo
       enddo
       
    else
       !      i j q
       EMITSA(:,:,:)     = 0.0_RP
       !      k i j w 2
       PM25SA(:,:,:)     = 0.0_RP
       WATRSA(:,:,:)     = 0.0_RP
       TAUNN   (:,:,:,:) = 0.0_RP
       TAUNNF  (:,:,  :) = 0.0_RP
       CEXTNN(:,:,:,:,:) = 0.0_RP
       CBAKNN(:,:,:,:,:) = 0.0_RP
       CEXTNF(:,:,:,  :) = 0.0_RP
       TAUNND  (:,:)     = 0.0_RP
       CEXTND(:,:,:)     = 0.0_RP
       QLTSA (:,:,:)     = 0.0_RP
    endif

    
    !carbon
    !-------------------------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA ) then
       AIRFEL  (:,:,:) = 0.0_RP    !8888888888888888888888888
       AIRFELOC(:,:,:) = 0.0_RP    !8888888888888888888888888
       call ATMOS_PHY_AE_SPRINTARS_AEROCA &
            ( QTRC(:,:,:,QS_AE_CA:QE_AE_CA), & ! [MODIFIED]
              OUTQLTCA(:,:,:,:),             & ! [OUT]
              NUMCON(:,:,:,:),               & ! [OUT]
              PM25CA(:,:,:),                 & ! [OUT]
              WATROC(:,:,:),                 & ! [OUT]
              TAUCA (:,:,:,:),               & ! [OUT]
              TAUCAA(:,:,:,:),               & ! [OUT]
              CEXTCA(:,:,:,:,:),             & ! [OUT]
              CABSCA(:,:,:,:,:),             & ! [OUT]
              CBAKCA(:,:,:,:,:),             & ! [OUT]
              TAUCAD(:,:),                   & ! [OUT]
              CEXTCD(:,:,:),                 & ! [OUT]
              CDVE(:,:),                     & ! [IN]
              CCOVMR(:,:),                   & ! [IN] OK
              CCMAX,                         & ! [IN] OK
              VABST(:,:),                    & ! [IN] OK
              RH(:,:,:),                     & ! [IN] OK
              DENS(:,:,:),                   & ! [IN] OK
              DENS1(:,:,:),                  & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:),      & ! [IN]
              FRAIN_MP_SNOW1(0:KA,:,:),      & ! [IN]
              FRAIN(0:KA,:,:),               & ! [IN] OK
              FRAIN_CP(0:KA,:,:),            & ! [IN]
              FRAIN_MP_RAIN(0:KA,:,:),       & ! [IN]
              FRAIN_MP_SNOW(0:KA,:,:),       & ! [IN]
              VTR_MP_RAIN(:,:,:),            & ! [IN]
              VTR_MP_SNOW(:,:,:),            & ! [IN]
              RADR_MP_RAIN(:,:,:),           & ! [IN]
              RADR_MP_SNOW(:,:,:),           & ! [IN]
              RAINN_CP(:,:,:),               & ! [IN]
              RAINN_MP_RAIN(:,:,:),          & ! [IN]
              RAINN_MP_SNOW(:,:,:),          & ! [IN]
              TCLDF(:,:,:),                  & ! [IN] OK
              TCLDF2(:,:,:),                 & ! [IN] OK
              KUP(:,:),                      & ! [IN] OK
              KUPMAX,                        & ! [IN] OK
              DPKUP(:,:),                    & ! [IN] OK
              KBB(:,:),                      & ! [IN] OK
              DPKBB(:,:),                    & ! [IN] OK
              GDPM(0:KA,:,:),                & ! [IN] OK
              DELP(:,:,:),                   & ! [IN] OK
              GDZM(0:KA,:,:),                & ! [IN] OK
              DELZ(:,:,:),                   & ! [IN] OK
              ONEWQ(1:KA-1,:,:),             & ! [IN] OK
              RNWR(:,:,:),                   & ! [IN] OK
              RNWR_CP(:,:,:),                & ! [IN] OK
              RNWR_MP_RAIN(:,:,:),           & ! [IN] OK
              RNWR_MP_SNOW(:,:,:),           & ! [IN] OK
              INDEX_PFT(:,:,:),              & ! [IN] OK
              FRAC_PFT(:,:,:),               & ! [IN] OK
              FACT_LAND(:,:),                & ! [IN] OK
              FACT_OCEAN1(:,:),              & ! [IN] OK
              GOICE(:,:),                    & ! [IN] OK
              AIRFEL(:,:,:),                 & ! [IN] fix
              AIRFELOC(:,:,:),               & ! [IN] fix
              EMITSA(:,:,:),                 & ! [IN] OK
              SFCFLX_SW_dn_net(:,:),         & ! [IN] OK
              GDZ(:,:,:),                    & ! [IN] OK
              OHRDC(:,:,:),                  & ! [IN] fix
              O3ORI(:,:,:),                  & ! [IN] fix
              NO3G(:,:,:),                   & ! [IN] fix
              GDT(:,:,:),                    & ! [IN] OK
              ACTI(:,:,:,2),                 & ! [IN]
              ACTI(:,:,:,3),                 & ! [IN]
              ACTI(:,:,:,4),                 & ! [IN]
              ATMOS_SW_PHY_AE_SPRINTARS_CCN, & ! [IN] OK
              TIME_NOWDATE(:),               & ! [IN] OK
              dt_AE,                         & ! [IN] OK
              QA_AE_CA                       ) ! [IN] OK
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             AE_Re(k,i,j,2) = 0.1000E-6_RP * 1.0E+2_RP !carbonaceous(BC/OC)=0.30 [m] -> [cm] 
             AE_Re(k,i,j,3) = 0.1000E-6_RP * 1.0E+2_RP !carbonaceous(BC/OC)=0.15 [m] -> [cm] 
             AE_Re(k,i,j,4) = 0.1000E-6_RP * 1.0E+2_RP !carbonaceous(BC/OC)=0.00 [m] -> [cm] 
             AE_Re(k,i,j,5) = 0.0118E-6_RP * 1.0E+2_RP !ext. BC                  [m] -> [cm] 
             AE_Qe(k,i,j,2) = OUTQLTCA(k,i,j,1)        !carbonaceous(BC/OC)=0.30 [kg/kg]
             AE_Qe(k,i,j,3) = OUTQLTCA(k,i,j,2)        !carbonaceous(BC/OC)=0.15 [kg/kg]
             AE_Qe(k,i,j,4) = OUTQLTCA(k,i,j,3)        !carbonaceous(BC/OC)=0.00 [kg/kg]
             AE_Qe(k,i,j,5) = OUTQLTCA(k,i,j,4)        !ext. BC                  [kg/kg]
          enddo
       enddo
       enddo
       NUMAE(:,:,:,3) = NUMCON(:,:,:,1) !carbon;internal mixed OM(OC) & BC
       NUMAE(:,:,:,4) = NUMCON(:,:,:,2) !carbon;SOC from NVOC
       NUMAE(:,:,:,2) = NUMCON(:,:,:,3) !carbon;external mixed BC
    
    else
       !      k i j w 2 
       PM25CA(:,:,:)     = 0.0_RP
       WATROC(:,:,:)     = 0.0_RP
       TAUCA   (:,:,:,:) = 0.0_RP
       TAUCAA  (:,:,:,:) = 0.0_RP
       CEXTCA(:,:,:,:,:) = 0.0_RP
       CABSCA(:,:,:,:,:) = 0.0_RP
       CBAKCA(:,:,:,:,:) = 0.0_RP
       TAUCAD  (:,:)     = 0.0_RP
       CEXTCD(:,:,:)     = 0.0_RP
       !      k i j 3
       NUMCON(:,:,:,:)   = 0.0_RP
    endif
    

    !sulfate
    !-------------------------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       AIRFELSO2(:,:,:) = 0.0_RP    
       call ATMOS_PHY_AE_SPRINTARS_AEROSU &
            ( QTRC(:,:,:,QS_AE_SU:QE_AE_SU),   & ! [MODIFIED]
              OUTQLD(:,:,:,QS_AE_SU:QE_AE_SU), & ! [OUT]
              NUMAE(:,:,:,5),                  & ! [OUT]
              PM25SU(:,:,:),                   & ! [OUT]
              WATRSU(:,:,:),                   & ! [OUT]
              TAUSU(:,:,:,:),                  & ! [OUT]
              CEXTSU(:,:,:,:,:),               & ! [OUT]
              CBAKSU(:,:,:,:,:),               & ! [OUT]
              TAUSUD(:,:),                     & ! [OUT]
              CEXTSD(:,:,:),                   & ! [OUT]
              CDVE(:,:),                       & ! [IN] 
              CCOVMR(:,:),                     & ! [IN] OK
              CCMAX,                           & ! [IN] OK
              GDT(:,:,:),                      & ! [IN] OK
              GOICE(:,:),                      & ! [IN] OK
              SFCFLX_SW_dn_net(:,:),           & ! [IN] OK
              RCOSZ(:,:),                      & ! [IN] 
              GRLAI(:,:),                      & ! [IN] fix
              KBB(:,:),                        & ! [IN] OK
              DPKBB(:,:),                      & ! [IN] OK
              THETA(:,:,:),                    & ! [IN] OK
              RH(:,:,:),                       & ! [IN] OK
              DENS(:,:,:),                     & ! [IN] OK
              DENS1(:,:,:),                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:),        & ! [IN]
              FRAIN_MP_SNOW1(0:KA,:,:),        & ! [IN]
              FRAIN(0:KA,:,:),                 & ! [IN] OK
              FRAIN_CP(0:KA,:,:),              & ! [IN] 
              FRAIN_MP_RAIN(0:KA,:,:),         & ! [IN] 
              FRAIN_MP_SNOW(0:KA,:,:),         & ! [IN] 
              VTR_MP_RAIN(:,:,:),              & ! [IN] 
              VTR_MP_SNOW(:,:,:),              & ! [IN] 
              RADR_MP_RAIN(:,:,:),             & ! [IN] 
              RADR_MP_SNOW(:,:,:),             & ! [IN]
              RAINN_CP(:,:,:),                 & ! [IN]
              RAINN_MP_RAIN(:,:,:),            & ! [IN]
              RAINN_MP_SNOW(:,:,:),            & ! [IN]
              QC(:,:,:),                       & ! [IN] 
              TCLDF(:,:,:),                    & ! [IN] OK
              TCLDF2(:,:,:),                   & ! [IN] OK
              KUP(:,:),                        & ! [IN] OK
              KUPMAX,                          & ! [IN] OK
              DPKUP(:,:),                      & ! [IN] OK
              FLIQ(:,:,:),                     & ! [IN] 
              OHRDC(:,:,:),                    & ! [IN] fix
              O3ORI(:,:,:),                    & ! [IN] fix
              GDZS(:,:),                       & ! [IN] fix
              GDPM(0:KA,:,:),                  & ! [IN] OK
              DELP(:,:,:),                     & ! [IN] OK
              GDZM(0:KA,:,:),                  & ! [IN] OK
              DELZ(:,:,:),                     & ! [IN] OK
              ONEWQ(1:KA-1,:,:),               & ! [IN] OK
              RNWR(:,:,:),                     & ! [IN] OK
              RNWR_CP(:,:,:),                  & ! [IN] 
              RNWR_MP_RAIN(:,:,:),             & ! [IN] 
              RNWR_MP_SNOW(:,:,:),             & ! [IN] 
              FACT_LAND(:,:),                  & ! [IN] OK
              FACT_OCEAN1(:,:),                & ! [IN] OK
              AIRFELSO2(:,:,:),                & ! [IN] OK
              ACTI(:,:,:,5),                   & ! [IN]
              ATMOS_SW_PHY_AE_SPRINTARS_CCN,   & ! [IN] OK
              TIME_NOWDATE(:),                 & ! [IN] OK
              dt_AE,                           & ! [IN] OK
              QA_AE_SU                         ) ! [IN] OK
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             AE_Re(k,i,j,6) = 0.0695E-6_RP * 1.0E+2_RP ! sulfate [m] -> [cm] 
             AE_Qe(k,i,j,6) = QTRC(k,i,j,QS_AE_SU)     ! sulfate [kg/kg]
          enddo
       enddo
       enddo
       
    else
       !      k i j w 2
       PM25SU(:,:,:)     = 0.0_RP
       WATRSU(:,:,:)     = 0.0_RP
       TAUSU   (:,:,:,:) = 0.0_RP
       CEXTSU(:,:,:,:,:) = 0.0_RP
       CBAKSU(:,:,:,:,:) = 0.0_RP
       TAUSUD  (:,:)     = 0.0_RP
       CEXTSD(:,:,:)     = 0.0_RP
    endif

    !passive tracer
    !-------------------------------
    if ( ATMOS_SW_PHY_AE_SPRINTARS_AT ) then
       call ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY2 &
            ( QTRC(:,:,:,QS_AE_AT:QE_AE_AT),   & ! [MODIFIED]
              CDVE(:,:),                       & ! [IN] 
              RH(:,:,:),                       & ! [IN] OK
              DENS(:,:,:),                     & ! [IN] OK
              DENS1(:,:,:),                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:),        & ! [IN]
              FRAIN_MP_SNOW1(0:KA,:,:),        & ! [IN]
              FRAIN(0:KA,:,:),                 & ! [IN] OK
              FRAIN_CP(0:KA,:,:),              & ! [IN] 
              FRAIN_MP_RAIN(0:KA,:,:),         & ! [IN] 
              FRAIN_MP_SNOW(0:KA,:,:),         & ! [IN] 
              VTR_MP_RAIN(:,:,:),              & ! [IN] 
              VTR_MP_SNOW(:,:,:),              & ! [IN] 
              RADR_MP_RAIN(:,:,:),             & ! [IN] 
              RADR_MP_SNOW(:,:,:),             & ! [IN]
              RAINN_CP(:,:,:),                 & ! [IN]
              RAINN_MP_RAIN(:,:,:),            & ! [IN]
              RAINN_MP_SNOW(:,:,:),            & ! [IN]
              TCLDF(:,:,:),                    & ! [IN] OK
              KUP(:,:),                        & ! [IN] OK
              KUPMAX,                          & ! [IN] OK
              DPKUP(:,:),                      & ! [IN] OK
              GDPM(0:KA,:,:),                  & ! [IN] OK
              DELP(:,:,:),                     & ! [IN] OK
              GDZM(0:KA,:,:),                  & ! [IN] OK
              DELZ(:,:,:),                     & ! [IN] OK
              RNWR(:,:,:),                     & ! [IN] OK
              RNWR_CP(:,:,:),                  & ! [IN] 
              RNWR_MP_RAIN(:,:,:),             & ! [IN] 
              RNWR_MP_SNOW(:,:,:),             & ! [IN] 
              LON_REAL(:,:),                   & ! [IN] OK
              LAT_REAL(:,:),                   & ! [IN] OK
              ACTI(:,:,:,5),                   & ! [IN]
              ATMOS_SW_PHY_AE_SPRINTARS_CCN,   & ! [IN] OK SW_CCN
              TIME_NOWDATE(:),                 & ! [IN] OK
              dt_AE,                           & ! [IN] OK
              QA_AE_AT                         ) ! [IN] OK
    endif

    
    !negative fixer
    !-------------------------------
    do iq = 1, QA_AE, 1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( QTRC(k,i,j,iq) < 0.0E+0_RP .or. QTRC(k,i,j,iq) > 1.0E-3_RP ) then
                write(*,'(a,i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)') 'Time: ', TIME_NOWDATE(1), '/', TIME_NOWDATE(2), &
                                '/', TIME_NOWDATE(3), ' ', TIME_NOWDATE(4), ':', TIME_NOWDATE(5), ':', TIME_NOWDATE(6)
                write(*,'(a,a10,a,i4,a,i4,a,i4,a,i4,a,i4,a,e13.5,a)') &
                     'Invalid value! ATMOS_PHY_AE_SPRINTARS_AEROSOL : ', trim(adjustl(TRCRNAME(iq))), &
                     ' [rank,iq,k,i,j]=[', PRC_myrank, ',', iq, ',', k, ',', i, ',', j, ',', QTRC(k,i,j,iq), ']'
                QTRC(k,i,j,iq) = max( min( QTRC(k,i,j,iq), 1.0E-3_RP ), 0.0E+0_RP )
             endif
          enddo
       enddo
       enddo
    enddo
    !-------------------------------
    !end negative fixer

    
    !check anomalous values
    !-------------------------------
#if defined(DEBUG) || defined(QUICKDEBUG)
    if ( ATMOS_SW_PHY_AE_SPRINTARS_DU ) then
       if (      minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) ) < 0.0E+0_RP &
            .or. maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) ) > 1.0E-3_RP ) then
          write(*,*) 'anomalous values DU', minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) ), maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) ), &
                                            minloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) ), maxloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_DU:QE_AE_DU) )
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'SPRINTARS DU has anomalous values'
          call PRC_abort
       endif
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SA ) then
       if (      minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) ) < 0.0E+0_RP &
            .or. maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) ) > 1.0E-3_RP ) then
          write(*,*) 'anomalous values SA', minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) ), maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) ), &
                                            minloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) ), maxloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SA:QE_AE_SA) )
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'SPRINTARS SA has anomalous values'
          call PRC_abort
       endif
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CA ) then
       if (      minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) ) < 0.0E+0_RP &
            .or. maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) ) > 1.0E-3_RP ) then
          write(*,*) 'anomalous values CA', minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) ), maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) ), &
                                            minloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) ), maxloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_CA:QE_AE_CA) )
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'SPRINTARS CA has anomalous values'
          call PRC_abort
       endif
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_SU ) then
       if (      minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) ) < 0.0E+0_RP &
            .or. maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) ) > 1.0E-3_RP ) then
          write(*,*) 'anomalous values SU', minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) ), maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) ), &
                                            minloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) ), maxloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_SU:QE_AE_SU) )
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'SPRINTARS SU has anomalous values'
          call PRC_abort
       endif
    endif
    if ( ATMOS_SW_PHY_AE_SPRINTARS_AT ) then
       if (      minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) ) < 0.0E+0_RP &
            .or. maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) ) > 1.0E-3_RP ) then
          write(*,*) 'anomalous values AT', minval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) ), maxval( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) ), &
                                            minloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) ), maxloc( QTRC(KS:KE,IS:IE,JS:JE,QS_AE_AT:QE_AE_AT) )
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_admin",*) 'SPRINTARS AT has anomalous values'
          call PRC_abort
       endif
    endif
#endif
    
    !================================================
    !end aerosol transport
    
        
    !optical thickness
    !================================================
    !total
    !-------------------------------
    do iw = 1, IWA
       do j = JS, JE
       do i = IS, IE
          !total            !soil dust         !carbon            !sulfate           !sea salt
          TAUT (i,j,iw,1) = TAUDU (i,j,iw,1) + TAUCA (i,j,iw,1) + TAUSU (i,j,iw,1) + TAUNN (i,j,iw,1)
          TAUTA(i,j,iw,1) = TAUDUA(i,j,iw,1) + TAUCAA(i,j,iw,1)
       enddo
       enddo
    enddo

    !fine mode, dry mode
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       !total         !soil dust        !carbon           !sulfate          !sea salt    
       TAUTF(i,j,1) = TAUDUF(i,j,  1) + TAUCA (i,j,1,1) + TAUSU (i,j,1,1) + TAUNNF(i,j,  1)
       TAUTD(i,j)   = TAUDU (i,j,1,1) + TAUCAD(i,j    ) + TAUSUD(i,j    ) + TAUNND(i,j    )
    enddo
    enddo

    !modification caused by cloud cover
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
          
       if ( CCOVMR(i,j) >= CCMAX ) then
          TAUT (i,j,1:IWA,2) = CONST_UNDEF !VMISS
          TAUTA(i,j,1:IWA,2) = CONST_UNDEF !VMISS
          TAUTF(i,j,      2) = CONST_UNDEF !VMISS
       else
          TAUT (i,j,1:IWA,2) = TAUT (i,j,1:IWA,1)
          TAUTA(i,j,1:IWA,2) = TAUTA(i,j,1:IWA,1)
          TAUTF(i,j,      2) = TAUTF(i,j,      1)
       endif
          
    enddo
    enddo
    !================================================
    !end optical thickness

    
    !Angstrom exponent (ALFA), single scattering albedo (SSA)
    !================================================
    do j = JS, JE
    do i = IS, IE
          
       if (       TAUT(i,j,1,1) > 0.0_RP &
            .and. TAUT(i,j,2,1) > 0.0_RP &
            .and. TAUT(i,j,3,1) > 0.0_RP ) then
          ALFA(i,j,1) = - log( TAUT(i,j,3,1)/TAUT(i,j,2,1) ) / log( WAVEL(3)/WAVEL(2) )
          SSA (i,j,1) = 1.0_RP - TAUTA(i,j,1,1) / TAUT(i,j,1,1)
       else
          ALFA(i,j,1) = CONST_UNDEF !VMISS
          SSA (i,j,1) = CONST_UNDEF !VMISS
       endif
       
       if (       TAUT(i,j,1,2) > 0.0_RP &
            .and. TAUT(i,j,2,2) > 0.0_RP &
            .and. TAUT(i,j,3,2) > 0.0_RP ) then
          ALFA(i,j,2) = - log( TAUT(i,j,3,2)/TAUT(i,j,2,2) ) / log( WAVEL(3)/WAVEL(2) )
          SSA (i,j,2) = 1.0_RP - TAUTA(i,j,1,2) / TAUT(i,j,1,2)
       else
          ALFA(i,j,2) = CONST_UNDEF !VMISS
          SSA (i,j,2) = CONST_UNDEF !VMISS
       endif

    enddo
    enddo
    !================================================
    !end Angstrom exponent (ALFA), single scattering albedo (SSA)

    
    !extinction, absorption & backscatter coefficient
    !================================================
    ! total
    !-------------------------------
    do iw = 1, IWA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
                
             ! extinction  coefficient (total)
             CEXTT(k,i,j,iw,1) = CEXTDU(k,i,j,iw,1) + CEXTCA(k,i,j,iw,1) & ! soil dust & carbon
                               + CEXTSU(k,i,j,iw,1) + CEXTNN(k,i,j,iw,1)   ! sulfate   & sea salt

             ! absorption  coefficient (total)
             CABST(k,i,j,iw,1) = CABSDU(k,i,j,iw,1) + CABSCA(k,i,j,iw,1)   ! soil dust & carbon

             ! backscatter coefficient (total)
             CBAKT(k,i,j,iw,1) = CBAKDU(k,i,j,iw,1) + CBAKCA(k,i,j,iw,1) & ! soil dust & carbon
                               + CBAKSU(k,i,j,iw,1) + CBAKNN(k,i,j,iw,1)   ! sulfate   & sea salt
             
          enddo
       enddo
       enddo
    enddo
    
    ! fine mode, dry mode, water
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          ! extinction  coefficient (total;fine mode)
          CEXTTF(k,i,j,1) = CEXTDF(k,i,j,  1) + CEXTCA(k,i,j,1,1) &     ! soil dust & carbon
                          + CEXTSU(k,i,j,1,1) + CEXTNF(k,i,j,  1)       ! sulfate   & sea salt

          ! extinction  coefficient (total;dry)
          CEXTTD(k,i,j,1) = CEXTDU(k,i,j,1,1) + CEXTCD(k,i,j    ) &     ! soil dust & carbon
                          + CEXTSD(k,i,j    ) + CEXTND(k,i,j    )       ! sulfate   & sea salt
             
          ! aerosol water (total)
          WATRAE(k,i,j) = WATROC(k,i,j) + WATRSU(k,i,j) + WATRSA(k,i,j) ! carbon, sulfate & sea salt

       enddo
    enddo
    enddo
 
    !modification caused by cloud cover
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
             
          if ( TCLDF2(k,i,j) >= CCMAX ) then
             CEXTT (k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CABST (k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CBAKT (k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CEXTTF(k,i,j,      2) = CONST_UNDEF !VMISS
             CEXTTD(k,i,j,      2) = CONST_UNDEF !VMISS
          else
             CEXTT (k,i,j,1:IWA,2) = CEXTT (k,i,j,1:IWA,1)
             CABST (k,i,j,1:IWA,2) = CABST (k,i,j,1:IWA,1)
             CBAKT (k,i,j,1:IWA,2) = CBAKT (k,i,j,1:IWA,1)
             CEXTTF(k,i,j,      2) = CEXTTF(k,i,j,      1)
             CEXTTD(k,i,j,      2) = CEXTTD(k,i,j,      1)
          endif

       enddo
    enddo
    enddo
    !================================================
    !end extinction, absorption & backscatter coefficient

    
    !PM2.5
    !================================================
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          
          ! submicron mass concentration
          PM25(k,i,j) = PM25DU(k,i,j) + PM25CA(k,i,j) & ! soil dust & carbon
                      + PM25SU(k,i,j) + PM25SA(k,i,j)   ! sulfate   & sea salt

          ! submicron mass mixing ratio
          QPM25(k,i,j) = PM25(k,i,j) / DENS (k,i,j)

       enddo
    enddo
    enddo
    !================================================
    !end PM2.5

    
    !CCN
    !================================================
    if ( ATMOS_SW_PHY_AE_SPRINTARS_CCN ) then
       call ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN &
            ( UAPCL(:,:,:),   & ! [OUT]
              NUMAET(:,:,:),  & ! [OUT]
              ACTI(:,:,:,:),  & ! [OUT]
              TCLDF(:,:,:),   & ! [IN] OK
              TCLDF2(:,:,:),  & ! [IN] OK
              QC(:,:,:),      & ! [IN] OK
              NUMAE(:,:,:,:), & ! [INOUT] OK
              RADAE(:,:,:,:), & ! [IN] OK
              DELZ(:,:,:),    & ! [IN] OK
              DENS(:,:,:),    & ! [IN] OK
              GDW(:,:,:),     & ! [IN]
              TKE_BL(:,:,:),  & ! [IN]
              GDT(:,:,:),     & ! [IN] OK
              GDP(:,:,:)      ) ! [IN] OK
    else
       UAPCL(:,:,:)   = 0.0_RP
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             do iq = 1, 6
                NUMAE(k,i,j,iq) = max( NUMAE(k,i,j,iq), 0.0_RP )
             enddo
             NUMAET(k,i,j) = sum( NUMAE(k,i,j,1:6) )
          enddo
       enddo
       enddo
       ACTI (:,:,:,:) = 0.0_RP
    endif
    !================================================
    !end CCN

    
    !Aerosol-Radiation Interaction
    !================================================
    if ( .not. ATMOS_SW_PHY_AE_SPRINTARS_ARI ) then
       AE_Re(:,:,:,:) = 0.0_RP
       AE_Qe(:,:,:,:) = 0.0_RP
    endif
    !================================================
    !end Aerosol-Radiation Interaction

    
    !Aerosol-Cloud Interaction
    !================================================
    if ( .not. ATMOS_SW_PHY_AE_SPRINTARS_ACI ) then
       UAPCL(:,:,:) = 0.0_RP
    endif
    !================================================
    !end Aerosol-Cloud Interaction


    !save
    !================================================
    FRAIN_MP_RAIN1(0:KA,:,:) = FRAIN_MP_RAIN(0:KA,:,:)
    FRAIN_MP_SNOW1(0:KA,:,:) = FRAIN_MP_SNOW(0:KA,:,:)
    DENS1         (:,   :,:) = DENS         (:   ,:,:)
    !================================================
    !end save

    
    !OUTPUT
    !================================================
    !SPRINTARS original
    !-------------------------------
    call FILE_HISTORY_in( GDTV         (1:KA,:,:), 'VT_SPRINTARS',           'virtual temperature',               'K',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( GDTVM        (1:KA,:,:), 'VT_HL_SPRINTARS',        'virtual temperature at half level', 'K',       dim_type = 'ZHXY' ) !
    call FILE_HISTORY_in( GDPM         (1:KA,:,:), 'PRES_HL_SPRINTARS',      'pressure at half level',            'Pa',      dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( RH           (1:KA,:,:), 'RH_SPRINTARS',           'ratio of qv to qsat',               'kg/kg',   dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( flux_qr_cp   (1:KA,:,:), 'FLUX_QR_CP_SPRINTARS',   'rain flux (CP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( flux_qs_cp   (1:KA,:,:), 'FLUX_QS_CP_SPRINTARS',   'snow flux (CP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( flux_qr_mp   (1:KA,:,:), 'FLUX_QR_MP_SPRINTARS',   'rain flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( flux_qs_mp   (1:KA,:,:), 'FLUX_QS_MP_SPRINTARS',   'snow flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( FRAIN        (1:KA,:,:), 'PREC_SPRINTARS',         'precipitation flux',                'kg/m2/s', dim_type = 'ZHXY' ) !check OK
    call FILE_HISTORY_in( FRAIN_CP     (1:KA,:,:), 'PREC_CP_SPRINTARS',      'precipitation flux (CP)',           'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( FRAIN_MP_RAIN(1:KA,:,:), 'RAIN_MP_SPRINTARS',      'rain flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( FRAIN_MP_SNOW(1:KA,:,:), 'SNOW_MP_SPRINTARS',      'snow flux (MP)',                    'kg/m2/s', dim_type = 'ZHXY' ) !check
    call FILE_HISTORY_in( VTR_MP_RAIN  (1:KA,:,:), 'VTR_MP_RAIN_SPRINTARS',  'terminal velocity (MP:rain)',       'm/s',     dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( VTR_MP_SNOW  (1:KA,:,:), 'VTR_MP_SNOW_SPRINTARS',  'terminal velocity (MP:snow)',       'm/s',     dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( cldfrac_cp   (1:KA,:,:), 'CLDFRAC_CP_SPRINTARS',   'cloud fraction (CP)',               '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( cldfrac_mp   (1:KA,:,:), 'CLDFRAC_MP_SPRINTARS',   'cloud fraction (MP)',               '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( TCLDF        (1:KA,:,:), 'TCLDF_SPRINTARS',        'TCLDF_SPRINTARS',                   '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( RNWR         (1:KA,:,:), 'RNWR_SPRINTARS',         'RNWR_SPRINTARS',                    '1',       dim_type = 'ZXY'  ) !check OK
    call FILE_HISTORY_in( RNWR_CP      (1:KA,:,:), 'RNWR_CP_SPRINTARS',      'RNWR_CP_SPRINTARS',                 '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( RNWR_MP_RAIN (1:KA,:,:), 'RNWR_MP_RAIN_SPRINTARS', 'RNWR_MP_RAIN_SPRINTARS',            '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( RNWR_MP_SNOW (1:KA,:,:), 'RNWR_MP_SNOW_SPRINTARS', 'RNWR_MP_SNOW_SPRINTARS',            '1',       dim_type = 'ZXY'  ) !check 
    call FILE_HISTORY_in( AIRFEL       (1:KA,:,:), 'AIRFEL_SPRINTARS',       'AIRFEL_SPRINTARS',                  '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( AIRFELOC     (1:KA,:,:), 'AIRFELOC_SPRINTARS',     'AIRFELOC_SPRINTARS',                '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( AIRFELSO2    (1:KA,:,:), 'AIRFELSO2_SPRINTARS',    'AIRFELSO2_SPRINTARS',               '1',       dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( OHRDC        (1:KA,:,:), 'OHRDC_SPRINTARS',        'OH concentration',                  '1/cm3',   dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( OHRAD        (1:KA,:,:), 'OHRAD_SPRINTARS',        'OH concentration',                  '1/cm3',   dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( O3ORI        (1:KA,:,:), 'O3ORI_SPRINTARS',        'O3 mixing ratio',                   'ppmv',    dim_type = 'ZXY'  ) !check
    call FILE_HISTORY_in( NO3G         (1:KA,:,:), 'NO3G_SPRINTARS',         'NO3(gas) mixing ratio',             'pptv',    dim_type = 'ZXY'  ) !check
    
    !2D
    !                            i j
    ! call FILE_HISTORY_in( U10O  (:,:), 'U10O_SPRINTARS',    'ocean 10m x-wind',           'm/s',   dim_type = 'XY' ) !check OK
    ! call FILE_HISTORY_in( V10O  (:,:), 'V10O_SPRINTARS',    'ocean 10m y-wind',           'm/s',   dim_type = 'XY' ) !check OK
    ! call FILE_HISTORY_in( VABSL (:,:), 'Uabs10L_SPRINTARS', '10m absolute wind on land',  'm/s',   dim_type = 'XY' ) !check OK
    ! call FILE_HISTORY_in( VABSO (:,:), 'Uabs10O_SPRINTARS', '10m absolute wind on ocean', 'm/s',   dim_type = 'XY' ) !check OK
    call FILE_HISTORY_in( VABST (:,:), 'Uabs10T_SPRINTARS', '10m absolute wind',          'm/s',   dim_type = 'XY' ) !check OK
    call FILE_HISTORY_in( CDVE  (:,:), 'CDVE_SPRINTARS',    'CDVE_SPRINTARS',             '1',     dim_type = 'XY' ) !
    call FILE_HISTORY_in( GRWG  (:,:), 'GRWG_SPRINTARS',    'GRWG_SPRINTARS',             'm3/m3', dim_type = 'XY' ) !check OK
    call FILE_HISTORY_in( CCOVMR(:,:), 'CCOVMR_SPRINTARS',  'CCOVMR_SPRINTARS',           '1',     dim_type = 'XY' ) !check OK
    call FILE_HISTORY_in( GRLAI (:,:), 'GRLAI_SPRINTARS',   'leaf area index',            'm2/m2', dim_type = 'XY' ) !check OK
    call FILE_HISTORY_in( GRSNW (:,:), 'GRSNW_SPRINTARS',   'snow depth',                 'm',     dim_type = 'XY' ) !check OK
    
    !reference
    !-------------------------------
    !3D
    !                            k i j 
    call FILE_HISTORY_in( QDRY  (:,:,:), 'QDRY_SPRINTARS',   'QDRY_SPRINTARS',    'kg/kg', dim_type = 'ZXY' ) !
    call FILE_HISTORY_in( QV    (:,:,:), 'QV_SPRINTARS',     'QV_SPRINTARS',      'kg/kg', dim_type = 'ZXY' ) !
    call FILE_HISTORY_in( QC    (:,:,:), 'QC_SPRINTARS',     'QC_SPRINTARS',      'kg/kg', dim_type = 'ZXY' ) !
    call FILE_HISTORY_in( QI    (:,:,:), 'QI_SPRINTARS',     'QI_SPRINTARS',      'kg/kg', dim_type = 'ZXY' ) !
    call FILE_HISTORY_in( GDW   (:,:,:), 'GDW_SPRINTARS',    'vertical velocity', 'm/s',   dim_type = 'ZXY' ) !
    call FILE_HISTORY_in( TKE_BL(:,:,:), 'TKE_BL_SPRINTARS', 'TKE(BL)',           'm2/s2', dim_type = 'ZXY' ) !
    
    !2D
    !                                         i j
    call FILE_HISTORY_in( LAND_WATER         (:,:), 'LAND_WATER_SPRINTARS',     'LAND_WATER_SPRINTARS',     'm3/m3', dim_type='XY' ) !check OK
    call FILE_HISTORY_in( LAND_PROPERTY      (:,:), 'LAND_PROPERTY_SPRINTARS',  'LAND_PROPERTY_SPRINTARS',  'm3/m3', dim_type='XY' ) !check OK
    call FILE_HISTORY_in( GOICE              (:,:), 'GOICE_SPRINTARS',          'GOICE_SPRINTARS',          '1',     dim_type='XY' ) !check OK
    call FILE_HISTORY_in( SFCFLX_SW_dn_net   (:,:), 'SFCFLX_SW_net_SPRINTARS',  'SFCFLX_SW_net_SPRINTARS',  'W/m2',  dim_type='XY' ) !check OK
    call FILE_HISTORY_in( RCOSZ              (:,:), 'RCOSZ_SPRINTARS',          'RCOSZ_SPRINTARS',          '1',     dim_type='XY' ) !
    call FILE_HISTORY_in( FACT_OCEAN1        (:,:), 'FACT_OCEAN1_SPRINTARS',    'FACT_OCEAN1_SPRINTARS',    '1',     dim_type='XY' ) !check OK
    call FILE_HISTORY_in( FACT_LAKE_SPRINTARS(:,:), 'FACT_LAKE_SPRINTARS',      'FACT_LAKE_SPRINTARS',      '1',     dim_type='XY' ) !
    
    index_pft_real(:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       index_pft_real(i,j,1) = real( INDEX_PFT(i,j,1), kind=RP )
       index_pft_real(i,j,2) = real( INDEX_PFT(i,j,2), kind=RP )
    enddo
    enddo
    !                                    i j 2
    call FILE_HISTORY_in( index_pft_real(:,:,1), 'INDEX_PFT1_SPRINTARS', 'INDEX_PFT1_SPRINTARS', '1', dim_type='XY' ) !check OK
    call FILE_HISTORY_in( index_pft_real(:,:,2), 'INDEX_PFT2_SPRINTARS', 'INDEX_PFT2_SPRINTARS', '1', dim_type='XY' ) !check OK
    call FILE_HISTORY_in( FRAC_PFT      (:,:,1), 'FRAC_PFT1_SPRINTARS',  'FRAC_PFT1_SPRINTARS',  '1', dim_type='XY' ) !check OK
    call FILE_HISTORY_in( FRAC_PFT      (:,:,2), 'FRAC_PFT2_SPRINTARS',  'FRAC_PFT2_SPRINTARS',  '1', dim_type='XY' ) !check OK


    !flux
    !-------------------------------
    !aerosol
    !-------------------------------
    !                            k i j
    call FILE_HISTORY_in( PM25  (:,:,:), 'PM25_SPRINTARS',   'mass concentration (PM2.5)', 'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QPM25 (:,:,:), 'QPM25_SPRINTARS',  'mass mixing ratio (PM2.5)',  'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( WATRAE(:,:,:), 'WATRAE_SPRINTARS', 'aerosol water',              'kg/kg', dim_type = 'ZXY' )

    !optical parameter
    !-------------------------------
    !                           i j w 2
    call FILE_HISTORY_in( TAUT (:,:,1,1), 'TAUT1_SPRINTARS',   'AOT (total;550nm;allsky)',               '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,1,1), 'TAUTA1_SPRINTARS',  'AOT (total;550nm;allsky;abs)',           '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTF(:,:,  1), 'TAUTF1_SPRINTARS',  'AOT (total;550nm;allsky;fine)',          '1', dim_type = 'XY' )
    call FILE_HISTORY_in( ALFA (:,:,  1), 'ALFA_SPRINTARS',    'angstrom exponent (440-870nm;allsky)',   '1', dim_type = 'XY' )
    call FILE_HISTORY_in( SSA  (:,:,  1), 'SSA_SPRINTARS',     'single scatter albedo (550nm;allsky)',   '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTD(:,:),     'TAUTD_SPRINTARS',   'AOT (total;550nm;allsky;dry)',           '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUT (:,:,2,1), 'TAUT2_SPRINTARS',   'AOT (total;440nm;allsky)',               '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,2,1), 'TAUTA2_SPRINTARS',  'AOT (total;440nm;allsky;abs)',           '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUT (:,:,3,1), 'TAUT3_SPRINTARS',   'AOT (total;870nm;allsky)',               '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,3,1), 'TAUTA3_SPRINTARS',  'AOT (total;870nm;allsky;abs)',           '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUT (:,:,1,2), 'TAUTC1_SPRINTARS',  'AOT (total;550nm;clearsky)',             '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,1,2), 'TAUTAC1_SPRINTARS', 'AOT (total;550nm;clearsky;abs)',         '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTF(:,:,  2), 'TAUTFC1_SPRINTARS', 'AOT (total;550nm;clearsky;fine)',        '1', dim_type = 'XY' )
    call FILE_HISTORY_in( ALFA (:,:,  2), 'ALFAC_SPRINTARS',   'angstrom exponent (440-870nm;clearsky)', '1', dim_type = 'XY' )
    call FILE_HISTORY_in( SSA  (:,:,  2), 'SSAC_SPRINTARS',    'single scatter albedo (550nm;clearsky)', '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUT (:,:,2,2), 'TAUTC2_SPRINTARS',  'AOT (total;440nm;clearsky)',             '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,2,2), 'TAUTAC2_SPRINTARS', 'AOT (total;440nm;clearsky;abs)',         '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUT (:,:,3,2), 'TAUTC3_SPRINTARS',  'AOT (total;870nm;clearsky)',             '1', dim_type = 'XY' )
    call FILE_HISTORY_in( TAUTA(:,:,3,2), 'TAUTAC3_SPRINTARS', 'AOT (total;870nm;clearsky;abs)',         '1', dim_type = 'XY' )

    !                            k i j w 2
    call FILE_HISTORY_in( CEXTT (:,:,:,1,1), 'EXTT1_SPRINTARS',   'extinction  (total;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,1,1), 'ABST1_SPRINTARS',   'absorption  (total;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,1,1), 'BAKT1_SPRINTARS',   'backscatter (total;532nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTTF(:,:,:,  1), 'EXTTF1_SPRINTARS',  'extinction  (total;532nm;allsky;fine)',   '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTTD(:,:,:,  1), 'EXTTD1_SPRINTARS',  'extinction  (total;532nm;allsky;dry)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTT (:,:,:,2,1), 'EXTT2_SPRINTARS',   'extinction  (total;1064nm;allsky)',       '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,2,1), 'ABST2_SPRINTARS',   'absorption  (total;1064nm;allsky)',       '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,2,1), 'BAKT2_SPRINTARS',   'backscatter (total;1064nm;allsky)',       '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTT (:,:,:,3,1), 'EXTT3_SPRINTARS',   'extinction  (total;355nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,3,1), 'ABST3_SPRINTARS',   'absorption  (total;355nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,3,1), 'BAKT3_SPRINTARS',   'backscatter (total;355nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTT (:,:,:,1,2), 'EXTTC1_SPRINTARS',  'extinction  (total;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,1,2), 'ABSTC1_SPRINTARS',  'absorption  (total;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,1,2), 'BAKTC1_SPRINTARS',  'backscatter (total;532nm;clearsky)',      '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTTF(:,:,:,  2), 'EXTTFC1_SPRINTARS', 'extinction  (total;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTTD(:,:,:,  2), 'EXTTDC1_SPRINTARS', 'extinction  (total;532nm;clearsky;dry)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTT (:,:,:,2,2), 'EXTTC2_SPRINTARS',  'extinction  (total;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,2,2), 'ABSTC2_SPRINTARS',  'absorption  (total;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,2,2), 'BAKTC2_SPRINTARS',  'backscatter (total;1064nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTT (:,:,:,3,2), 'EXTTC3_SPRINTARS',  'extinction  (total;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABST (:,:,:,3,2), 'ABSTC3_SPRINTARS',  'absorption  (total;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKT (:,:,:,3,2), 'BAKTC3_SPRINTARS',  'backscatter (total;355nm;clearsky)',      '1/m/sr', dim_type = 'ZXY' )
    !================================================
    !end OUTPUT

    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSL
  !-----------------------------------------------------------------------------


  
end module scale_atmos_phy_ae_SPRINTARS_admin
