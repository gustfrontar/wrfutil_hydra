!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - dust
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_aerodust
  !-----------------------------------------------------------------------------
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
       PRC_abort
  
  use scale_statistics, only: &
       STATISTICS_detail

  use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE

  use scale_const, only: &
       PI   => CONST_PI,   & !< pi
       GRAV => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       CONST_D2R,          &
       CONST_UNDEF
  
  use scale_file_history, only: &
       FILE_HISTORY_in
  
  use scale_atmos_phy_ae_SPRINTARS_common, only: &
       ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL,        &
       ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z, &
       ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC,          &
       NDU,                                           & !< total number of dust tracers
       IWA,                                           & !< number of wavelengths
       !VMISS,                                         & !< missing value
       NU,                                            & !< air viscosity [kg/m/s] = [Pa*s]
       RADR,                                          & !< rain droplet radius      [m]
       VTR,                                           & !< rain terminal velocity   [m/s]
       DENSR,                                         & !< rain droplet density     [kg/m3]
       THRESP                                           !< precipitation threshold  [kg/m2/s]
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_AERODU_setup
  public :: ATMOS_PHY_AE_SPRINTARS_AERODU_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_AERODU
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! parameter
  ! variables
  real(RP), private, save :: FINCDU = 0.10E+0_RP !< incloud coefficient
  real(RP), private, save :: UTR    = 6.50E+0_RP !< threshold 10m wind speed [m/s]
  real(RP), private, save :: TSNOW  = 0.10E+0_RP !< threshold of snow        [m]
  real(RP), private, save :: &
       WATER(7) &                                !< threshold of soilwater   [m/m]
       = (/ 1.0E-1_RP, 1.0E-1_RP, 1.0E-1_RP, 2.0E-1_RP, 1.0E-1_RP, 1.0E-1_RP, 1.0E-1_RP /)
  real(RP), private, save :: &
       EMFDU(7) &                                !< emission coefficient
       = (/ 1.3E-9_RP, 1.3E-9_RP, 1.3E-9_RP, 1.3E-9_RP, 1.3E-9_RP, 1.3E-9_RP, 1.3E-9_RP /)

  real(RP), private, allocatable, save :: QTRCINP(:,:,:,:) !< mixing ratio in prec(MP)

  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Dust setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AERODU_setup( QA_AE_DU )
    implicit none
    integer, intent(in) :: QA_AE_DU
    integer :: ierr, iq
    character :: cnumber*2
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AERODU_NMAEDU / &
         FINCDU, &
         UTR,    &
         WATER,  &
         TSNOW,  &
         EMFDU
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerodust",*) 'Setup'

    !--- read namelist
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AERODU_NMAEDU, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerodust",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AERODU_NMAEDU. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_aerodust",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AERODU_NMAEDU, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AERODU_NMAEDU)

    allocate( QTRCINP(KA,IA,JA,NDU) )
    QTRCINP(:,:,:,:) = 0.0_RP

    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP
    call FILE_HISTORY_in( zero_2d(:,:), 'EFDU_SPRINTARS',           'emission flux (dust)',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTDU_SPRINTARS',         'wet deposition flux (dust)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTDU_CP_SPRINTARS',      'wet deposition flux (dust;CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTDU_MP_RAIN_SPRINTARS', 'wet deposition flux (dust;MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTDU_MP_SNOW_SPRINTARS', 'wet deposition flux (dust;MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'DRFTDU_SPRINTARS',         'dry deposition flux (dust)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'GRFTDU_SPRINTARS',         'gravitational settling flux (dust)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'TDFTDU_SPRINTARS',         'total deposition flux (dust)',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'EFDU25_SPRINTARS',         'emission flux (dust,PM2.5)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WFDU25_SPRINTARS',         'wet deposition flux (dust,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'DFDU25_SPRINTARS',         'dry deposition flux (dust,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'GFDU25_SPRINTARS',         'gravitational settling flux (dust,PM2.5)', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDDT_SPRINTARS',  'mixing ratio (dust)',               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONDT_SPRINTARS',  'mass concentration (dust)',         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCDT_SPRINTARS',  'number concentration (dust)',       '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:), 'COLDT_SPRINTARS',  'colomn loading (dust)',             'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QPM25D_SPRINTARS', 'mixing ratio (dust,PM2.5)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NPM25D_SPRINTARS', 'number concentration (dust,PM2.5)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:), 'SKUP_SPRINTARS',   'uplift level',                      'level', dim_type = 'XY'  )
    do iq = 1, QA_AE_DU
       write(cnumber,'(i2.2)') iq
       call FILE_HISTORY_in( zero_3d(:,:,:), 'COND'//cnumber//'_SPRINTARS', 'mass concentration (dust'//cnumber//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCD'//cnumber//'_SPRINTARS', 'number concentration (dust'//cnumber//')', '1/m3',  dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_2d  (:,:), 'COLD'//cnumber//'_SPRINTARS', 'column loading (dust'//cnumber//')',       'kg/m2', dim_type = 'XY'  )
    enddo
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDU_SPRINTARS',  'AOT (dust;550nm;allsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDUA_SPRINTARS', 'AOT (dust;550nm;allsky;abs)',         '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDUF_SPRINTARS', 'AOT (dust;550nm;allsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDU1_SPRINTARS', 'extinction (dust;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDU1_SPRINTARS', 'absorption (dust;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDU1_SPRINTARS', 'backscatter (dust;532nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDF1_SPRINTARS', 'extinction (dust;532nm;allsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDU2_SPRINTARS', 'extinction (dust;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDU2_SPRINTARS', 'absorption (dust;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDU2_SPRINTARS', 'backscatter (dust;1064nm;allsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDU3_SPRINTARS', 'extinction (dust;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDU3_SPRINTARS', 'absorption (dust;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDU3_SPRINTARS', 'backscatter (dust;355nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONDTC_SPRINTARS',  'mass concentration (dust;clearsky)',    'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCDTC_SPRINTARS',  'number concentration (dust;clearsky)',  '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLDTC_SPRINTARS',  'column loading (dust;clearsky)',        'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDUC_SPRINTARS',  'AOT (dust;550nm;clearsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDUAC_SPRINTARS', 'AOT (dust;550nm;clearsky;abs)',         '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUDUFC_SPRINTARS', 'AOT (dust;550nm;clearsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDUC1_SPRINTARS', 'extinction (dust;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDUC1_SPRINTARS', 'absorption (dust;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDUC1_SPRINTARS', 'backscatter (dust;532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDFC1_SPRINTARS', 'extinction (dust;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDUC2_SPRINTARS', 'extinction (dust;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDUC2_SPRINTARS', 'absorption (dust;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDUC2_SPRINTARS', 'backscatter (dust;1064nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTDUC3_SPRINTARS', 'extinction (dust;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSDUC3_SPRINTARS', 'absorption (dust;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKDUC3_SPRINTARS', 'backscatter (dust;355nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AERODU_setup
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  !> SPRINTARS Dust finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_AERODU_finalize
    implicit none

    deallocate( QTRCINP )

    deallocate( zero_3d )
    deallocate( zero_2d )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AERODU_finalize
  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Dust
  subroutine ATMOS_PHY_AE_SPRINTARS_AERODU( &
       QTRC,           & ! [MODIFIED]
       QLOAD,          & ! [OUT]   OUTQLTDU
       NUMCNT,         & ! [OUT]   NUMAE
       RADDU,          & ! [OUT]   RADAE
       PM25DU,         & ! [OUT]   PM25DU
       TAUDU,          & ! [OUT]   TAUDU
       TAUDUA,         & ! [OUT]   TAUDUA
       TAUDUF,         & ! [OUT]   TAUDUF
       CEXTDU,         & ! [OUT]   CEXTDU
       CABSDU,         & ! [OUT]   CABSDU
       CBAKDU,         & ! [OUT]   CBAKDU
       CEXTDF,         & ! [OUT]   CEXTDF
       CDVE,           & ! [IN]    CDVE
       CCOVMR,         & ! [IN]    CCOVMR
       CCMAX,          & ! [IN]    CCMAX
       WS10,           & ! [IN]    VABST
       GRWG,           & ! [IN]    GRWG
       GRSNW,          & ! [IN]    GRSNW
       GRLAI,          & ! [IN]    GRLAI
       DENS,           & ! [IN]    DENS
       DENS1,          & ! [IN]
       FRAIN_MP_RAIN1, & ! [IN]
       FRAIN_MP_SNOW1, & ! [IN]
       FRAIN,          & ! [IN]    FRAIN
       FRAIN_CP,       & ! [IN]
       FRAIN_MP_RAIN,  & ! [IN]
       FRAIN_MP_SNOW,  & ! [IN]
       VTR_MP_RAIN,    & ! [IN]
       VTR_MP_SNOW,    & ! [IN]
       RADR_MP_RAIN,   & ! [IN]
       RADR_MP_SNOW,   & ! [IN]
       RAINN_CP,       & ! [IN]
       RAINN_MP_RAIN,  & ! [IN]
       RAINN_MP_SNOW,  & ! [IN]
       TCLDF,          & ! [IN]    TCLDF
       TCLDF2,         & ! [IN]    TCLDF2
       KUP,            & ! [IN]    KUP
       KUPMAX,         & ! [IN]    KUPMAX
       DPKUP,          & ! [IN]    DPKUP
       GDPM,           & ! [IN]    GDPM
       DELP,           & ! [IN]    DELP
       GDZM,           & ! [IN]    GDZM
       DELZ,           & ! [IN]    DELZ
       ONEWQ,          & ! [IN]    ONEWQ
       RNWR,           & ! [IN]    RNWR
       RNWR_CP,        & ! [IN]
       RNWR_MP_RAIN,   & ! [IN]
       RNWR_MP_SNOW,   & ! [IN]
       LON_REAL,       & ! [IN]
       LAT_REAL,       & ! [IN]
       INDEX_PFT,      & ! [IN]
       FRAC_PFT,       & ! [IN]
       FACT_LAND,      & ! [IN]
       Zsfc,           & ! [IN]
       HGT_DUST_Tibet, & ! [IN]
       ACTIDU,         & ! [IN]
       SW_CCN,         & ! [IN]
       dt_AE,          & ! [IN]
       QA_AE_DU        ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE                        !< time step
    integer,  intent(in)    :: QA_AE_DU                     !< total number of dust tracers

    real(RP), intent(in)    :: CDVE               (IA,JA)   !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
    real(RP), intent(in)    :: CCOVMR             (IA,JA)   !< cloud cover (max-ran)
    real(RP), intent(in)    :: CCMAX                        !< maximum cloud fraction for all/clear sky
    real(RP), intent(in)    :: WS10               (IA,JA)   !< wind at 10m (land+ocean) [m/s]
    real(RP), intent(in)    :: GRWG               (IA,JA)   !< subsurf. moist. saturat.[0-1] [m3/m3]
    real(RP), intent(in)    :: GRSNW              (IA,JA)   !< snow depth                      [m]
    real(RP), intent(in)    :: GRLAI              (IA,JA)   !< Leaf Area Index [m2/m2]
    real(RP), intent(in)    :: DENS          (KA,  IA,JA)   !< air density
    real(RP), intent(in)    :: DENS1         (KA,  IA,JA)   !< air density   on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA)   !< precipitation flux (liq    ) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA)   !< precipitation flux (    ice) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN         (0:KA,IA,JA)   !< precipitation flux (liq+ice) : CP + MP
    real(RP), intent(in)    :: FRAIN_CP      (0:KA,IA,JA)   !< precipitation flux (liq+ice) : CP
    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA)   !< precipitation flux (liq    ) : MP
    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA)   !< precipitation flux (    ice) : MP
    real(RP), intent(in)    :: VTR_MP_RAIN   (KA,  IA,JA)   !< terminal velocity of rain (microphysic scheme) [m/s]
    real(RP), intent(in)    :: VTR_MP_SNOW   (KA,  IA,JA)   !< terminal velocity of snow (microphysic scheme) [m/s]
    real(RP), intent(in)    :: RADR_MP_RAIN  (KA,  IA,JA)   !< effective radius of rain (microphysic scheme) [m]
    real(RP), intent(in)    :: RADR_MP_SNOW  (KA,  IA,JA)   !< effective radius of snow (microphysic scheme) [m]
    real(RP), intent(in)    :: RAINN_CP      (KA,  IA,JA)   !< rain number concentration (CP)      [1/m3]
    real(RP), intent(in)    :: RAINN_MP_RAIN (KA,  IA,JA)   !< rain number concentration (MP;RAIN) [1/m3]
    real(RP), intent(in)    :: RAINN_MP_SNOW (KA,  IA,JA)   !< rain number concentration (MP;SNOW) [1/m3]
    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA)   !< cloud fraction (CP+MP)
    real(RP), intent(in)    :: TCLDF2        (KA,  IA,JA)   !< cloud fraction (CP+MP) for CCOVMR
    integer,  intent(in)    :: KUP                (IA,JA)   !< uplift level
    integer,  intent(in)    :: KUPMAX                       !< maximum uplift level
    real(RP), intent(in)    :: DPKUP              (IA,JA)   !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA)   !< pressure at half levels
    real(RP), intent(in)    :: DELP          (KA,  IA,JA)   !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA)   !< height at half levels
    real(RP), intent(in)    :: DELZ          (KA,  IA,JA)   !< GDZM(k,i,j)-GDZM(k-1,i,j)
    logical,  intent(in)    :: ONEWQ         (KA-1,IA,JA)   !< use for instability
    real(RP), intent(in)    :: RNWR          (KA,  IA,JA)   !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA)   !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA)   !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA)   !< ratio of rain to cloud
    real(RP), intent(in)    :: LON_REAL           (IA,JA)   !< longitude
    real(RP), intent(in)    :: LAT_REAL           (IA,JA)   !< latitude
    integer,  intent(in)    :: INDEX_PFT          (IA,JA,2) !< index of PFT for each mosaic
    real(RP), intent(in)    :: FRAC_PFT           (IA,JA,2) !< fraction of PFT for each mosaic
    real(RP), intent(in)    :: FACT_LAND          (IA,JA)   !< land factor
    real(RP), intent(in)    :: Zsfc               (IA,JA)   !< height of surface [m]
    real(RP), intent(in)    :: HGT_DUST_Tibet               !< Upper limit of dust emission eleveation at Tibet area
    real(RP), intent(in)    :: ACTIDU        (KA,  IA,JA)   !< ratio of activated aerosols [0-1]
    logical,  intent(in)    :: SW_CCN                       !< switch for CCN scheme (true->use ACTI)
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_DU)   !< ratio of dust tracer mass to total mass
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: QLOAD (KA,IA,JA,QA_AE_DU) !< QTRC for radiation
    real(RP), intent(out)   :: NUMCNT(KA,IA,JA)          !< number concentration (total)
    real(RP), intent(out)   :: RADDU (KA,IA,JA)          !< mode radius of dust tracers
    real(RP), intent(out)   :: PM25DU(KA,IA,JA)          !< submicron mass concentration
    real(RP), intent(out)   :: TAUDU    (IA,JA,IWA,2)    !< AOT
    real(RP), intent(out)   :: TAUDUA   (IA,JA,IWA,2)    !< AOT (absorption)
    real(RP), intent(out)   :: TAUDUF   (IA,JA,    2)    !< AOT (fine mode)
    real(RP), intent(out)   :: CEXTDU(KA,IA,JA,IWA,2)    !< extinction  coefficient
    real(RP), intent(out)   :: CABSDU(KA,IA,JA,IWA,2)    !< absorption  coefficient
    real(RP), intent(out)   :: CBAKDU(KA,IA,JA,IWA,2)    !< backscatter coefficient
    real(RP), intent(out)   :: CEXTDF(KA,IA,JA,    2)    !< extinction  coefficient (fine mode)
    
    ! internal work
    !-------------------------
    !parameter
    real(RP) :: VTER(QA_AE_DU)             !< terminal velocity                      [m/s]
    real(RP) :: RADN(QA_AE_DU)             !< aerosol radius of number concentration [m]
    real(RP) :: MD  (QA_AE_DU)             !< dust particle mass                     [kg]

    !emission
    real(RP) :: CC    (IA,JA)              !< emission coefficient
    real(RP) :: QEMIT (IA,JA,QA_AE_DU)     !< emission mixing ratio [kg/kg]
    real(RP) :: EMITF (IA,JA)              !< emission flux         [kg/m2/s]
    real(RP) :: SKUP  (IA,JA)              !< uplift level
    real(RP) :: REGION(IA,JA)              !< region
    real(RP) :: FPAR  (IA,JA)              !< absorbed photosyn. rad.
    real(RP) :: AEFF  (IA,JA)              !< effective surf. frac.
    real(RP) :: CEMFDU                     !< use for emission
    real(RP) :: TWATER                     !< use for emission

    !incloud
    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_DU) !< outside cloud water [kg/kg]
    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_DU) !< inside  cloud water [kg/kg]

    !subcloud scavenging (wash out)
    ! real(RP) :: RAINN         (KA,IA,JA)          !< rain number concentration     [1/m3]
    real(RP) :: DUSTN                             !< number density of dust tracer [1/m3]
    real(RP) :: DUSTRN        (KA,IA,JA,QA_AE_DU) !< number density of dust tracer entering into raindrops [1/s/m3] <=> collision frequency between dust and rain
    real(RP) :: DUSTRN_CP     (KA,IA,JA,QA_AE_DU)
    real(RP) :: DUSTRN_MP_RAIN(KA,IA,JA,QA_AE_DU)
    real(RP) :: DUSTRN_MP_SNOW(KA,IA,JA,QA_AE_DU)
    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_DU) !< mixing ratio by wet deposition [kg/kg]
    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_DU)
    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_DU)
    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_DU)
    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_DU) !< dust tracers that do not cllide with rain particles [kg/kg]

    !incloud scavenging (rain out)
    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_DU) !< mixing ratio moving cloud to rain particles [kg/kg]
    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_DU)
    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_DU)
    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_DU)
    
    !re-emission from rain
    real(RP) :: WETF           (IA,JA,QA_AE_DU) !< wet deposition flux
    real(RP) :: WETF_CP        (IA,JA,QA_AE_DU)
    real(RP) :: WETF_MP_RAIN   (IA,JA,QA_AE_DU)
    real(RP) :: WETF_MP_SNOW   (IA,JA,QA_AE_DU)
    real(RP) :: QTRCINR     (KA,IA,JA,QA_AE_DU) !< mixing ratio in rain(MP)
    real(RP) :: QTRCINS     (KA,IA,JA,QA_AE_DU) !< mixing ratio in snow(MP)
    
    !dry deposition
    real(RP) :: QMTX(KA,IA,JA,-1:1)     !< implicit matrix of QTRC
    real(RP) :: DRYF   (IA,JA,QA_AE_DU) !< dry deposition flux

    !gravitational settling
    real(RP) :: VTERG(KA,IA,JA)          !< terminal velocity for gravitational settling
    real(RP) :: GRAVF   (IA,JA,QA_AE_DU) !< gravitational settling flux

    !file output (flux)
    real(RP) :: WETFT        (IA,JA)       !< wet deposition flux (total)
    real(RP) :: WETFT_CP     (IA,JA)       !< wet deposition flux (total;CP)
    real(RP) :: WETFT_MP_RAIN(IA,JA)       !< wet deposition flux (total;MP(RAIN))
    real(RP) :: WETFT_MP_SNOW(IA,JA)       !< wet deposition flux (total;MP(SNOW))
    real(RP) :: DRYFT        (IA,JA)       !< dry deposition flux (total)
    real(RP) :: GRAVFT       (IA,JA)       !< gravitational settling flux (total)
    real(RP) :: DEPFT        (IA,JA)       !< total deposition flux
    real(RP) :: EMIF25       (IA,JA)       !< emission flux (PM2.5)
    real(RP) :: WETF25       (IA,JA)       !< wet deposition flux (PM2.5)
    real(RP) :: DRYF25       (IA,JA)       !< dry deposition flux (PM2.5)
    real(RP) :: GRVF25       (IA,JA)       !< gravitational settling flux (PM2.5)
    
    !file output (aerosol)
    real(RP) :: QLDDT  (KA,IA,JA)          !< mixing ratio (total) [kg/kg]
    real(RP) :: COLMST    (IA,JA)          !< column mass loading (total) [kg/m2]
    real(RP) :: DCON   (KA,IA,JA,QA_AE_DU) !< mass concentration
    real(RP) :: DCONT  (KA,IA,JA)          !< mass concentraion (total)
    real(RP) :: NUMCON (KA,IA,JA,QA_AE_DU) !< number concentration
    real(RP) :: COLMS     (IA,JA,QA_AE_DU) !< column mass loading [kg/m2]
    real(RP) :: QPM25D (KA,IA,JA)          !< mass mixing ratio (PM2.5) [kg/kg]
    real(RP) :: NPM25D (KA,IA,JA)          !< number concentration (PM2.5) [1/m3]
    
    !file output (optical parameter)
    real(RP) :: DCONTC (KA,IA,JA)          !< mass concentraion (total,clear sky) [kg/m3]
    real(RP) :: NUMCTC (KA,IA,JA)          !< number concentration (total,clear sky) [1/m3]
    real(RP) :: COLMTC    (IA,JA)          !< column mass loading (total,clear sky) [kg/m2]
    real(RP) :: TUSE                       !< use for AOT
    
    ! parameter
    !-------------------------
    real(RP), parameter :: &
         CEXTDP(NDU,IWA)   & !<  extinction cross section (550/ 440/870nm)
         !             dust01       dust02       dust03       dust04       dust05       dust06 
         = reshape( (/ 1.917E+3_RP, 3.596E+3_RP, 8.987E+2_RP, 5.204E+2_RP, 1.979E+2_RP, 7.545E+1_RP, & !  550nm
                       3.633E+3_RP, 3.275E+3_RP, 8.608E+2_RP, 5.240E+2_RP, 1.936E+2_RP, 7.502E+1_RP, & !  440nm
                       4.330E+2_RP, 2.208E+3_RP, 9.567E+2_RP, 5.592E+2_RP, 2.026E+2_RP, 7.689E+1_RP /), & !  870nm
                    (/ NDU, IWA /) )
    real(RP), parameter :: &
         CABSDP(NDU,IWA)   & !<  absorption cross section (550/ 440/870nm)
         !             dust01       dust02       dust03       dust04       dust05       dust06 
         = reshape( (/ 2.539E+1_RP, 3.470E+1_RP, 3.632E+1_RP, 3.122E+1_RP, 2.355E+1_RP, 1.730E+1_RP, & !  550nm
                       3.699E+1_RP, 4.454E+1_RP, 4.372E+1_RP, 3.555E+1_RP, 2.730E+1_RP, 1.959E+1_RP, & !  440nm
                       1.172E+1_RP, 1.954E+1_RP, 2.258E+1_RP, 2.271E+1_RP, 1.696E+1_RP, 1.287E+1_RP /), & !  870nm
                    (/ NDU, IWA /) )
    real(RP), parameter :: &
         CEXTDA(NDU,IWA)   & !<  extinction cross section (532/1064/355nm)
         !             dust01       dust02       dust03       dust04       dust05       dust06 
         = reshape( (/ 2.123E+3_RP, 3.611E+3_RP, 9.138E+2_RP, 5.213E+2_RP, 1.972E+2_RP, 7.539E+1_RP, & !  532nm
                       2.074E+2_RP, 1.533E+3_RP, 1.301E+3_RP, 4.739E+2_RP, 2.071E+2_RP, 7.718E+1_RP, & ! 1064nm
                       5.453E+3_RP, 2.348E+3_RP, 8.061E+2_RP, 5.099E+2_RP, 1.927E+2_RP, 7.461E+1_RP /), & !  355nm
                    (/ NDU, IWA /) )
    real(RP), parameter :: &
         CABSDA(NDU,IWA)   & !<  extinction cross section (532/1064/355nm)
         !             dust01       dust02       dust03       dust04       dust05       dust06 
         = reshape( (/ 2.714E+1_RP, 3.621E+1_RP, 3.754E+1_RP, 3.170E+1_RP, 2.409E+1_RP, 1.765E+1_RP, & !  532nm
                       8.766E+0_RP, 1.535E+1_RP, 1.823E+1_RP, 1.868E+1_RP, 1.452E+1_RP, 1.123E+1_RP, & ! 1064nm
                       4.905E+1_RP, 5.752E+1_RP, 4.956E+1_RP, 4.300E+1_RP, 3.250E+1_RP, 2.241E+1_RP /), & !  355nm
                    (/ NDU, IWA /) )
    real(RP), parameter :: &
         CBAKDA(NDU,IWA)   & !< backscatter cross section (532/1064/355nm)
         !             dust01       dust02       dust03       dust04       dust05       dust06 
         = reshape( (/ 2.462E+1_RP, 8.279E+1_RP, 1.432E+2_RP, 2.579E+1_RP, 5.226E+0_RP, 1.523E+0_RP, & !  532nm
                       1.704E+1_RP, 1.788E+1_RP, 5.007E+1_RP, 7.980E+1_RP, 6.879E+0_RP, 3.069E+0_RP, & ! 1064nm
                       7.463E+1_RP, 2.268E+2_RP, 4.619E+1_RP, 2.659E+1_RP, 8.428E+0_RP, 1.025E+0_RP /), & !  355nm
                    (/ NDU, IWA /) )
    real(RP), parameter :: &
         WSRC(NDU)         & !< source strength
         !    dust01        dust02        dust03        dust04        dust05        dust06 
         = (/ 0.0045E+0_RP, 0.0290E+0_RP, 0.1766E+0_RP, 0.2633E+0_RP, 0.2633E+0_RP, 0.2633E+0_RP /)
    real(RP), parameter :: &
         AA  (NDU)       & !< collision efficiency
         !    dust01       dust02       dust03       dust04       dust05       dust06 
         = (/ 1.000E-4_RP, 1.000E-4_RP, 1.000E-3_RP, 1.000E-2_RP, 0.631E+0_RP, 1.000E+0_RP /)
    real(RP), parameter :: &
         RADD(NDU)         & !< dust particle radius [um]
         !    dust01      dust02      dust03      dust04      dust05      dust06 
         = (/ 0.13E-6_RP, 0.33E-6_RP, 0.82E-6_RP, 1.27E-6_RP, 3.20E-6_RP, 8.02E-6_RP /)
    real(RP), parameter :: DENSD  = 2.60E+3_RP !< dust particle density    [kg/m3]
    real(RP), parameter :: SIGMA  = 1.10E+0_RP !< GSD of log-normal
    
    ! others 
    !-------------------------
    integer   :: k, i, j, iq, iw
    character :: cnumber*2
    real(RP) :: check_2d_data   (IA,JA,100)
    real(RP) :: check_3d_data(KA,IA,JA,100)

    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - dust'

    !setup parameter 
    !================================================

    !RADN, MD, VTER
    !-------------------------------
    do iq = 1, QA_AE_DU
       RADN(iq) = RADD(iq) / exp( 3.0_RP * log(SIGMA)**2 )
       MD  (iq) = 4.0_RP / 3.0_RP * PI * RADN(iq)**3 * DENSD * exp( 4.5_RP * log(SIGMA)**2 )
       VTER(iq) = 2.0_RP * DENSD * RADN(iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA)**2 )
    enddo

    !REGION
    !-------------------------------

    !================================================
    !end setup parameter

    
    !aerosol transport 
    !================================================
    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_DU
       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
    enddo

    
    !emission
    !-------------------------------
    CC(:,:) = 0.0_RP
    SKUP(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       CEMFDU = EMFDU(1)        !88888888888888
       TWATER = WATER(1)        !88888888888888

       if (       WS10 (i,j) >= UTR    &
            .and. GRWG (i,j) >  0.0_RP &
            .and. GRWG (i,j) <= TWATER &
            .and. GRSNW(i,j) >= 0.0_RP &
            .and. GRSNW(i,j) <= TSNOW  ) then ! snow_depth is shorter than 0.1m

          if (       INDEX_PFT(i,j,1) >= 1 .and. INDEX_PFT(i,j,1) <= 9 &
                .and. ( .not. ( Zsfc(i,j) >= HGT_DUST_Tibet &  ! exclude Tibetan Plateau
                          .and. LON_REAL(i,j) >= 65.0_RP * CONST_D2R .and. LON_REAL(i,j) <= 110.0_RP * CONST_D2R &
                          .and. LAT_REAL(i,j) >= 24.0_RP * CONST_D2R .and. LAT_REAL(i,j) <=  41.3_RP * CONST_D2R ) ) ) then
             CC(i,j) = CC(i,j) + FACT_LAND(i,j) * FRAC_PFT(i,j,1) * CEMFDU * ( TWATER -  GRWG(i,j) ) / TWATER &
                                                                           * ( TSNOW  - GRSNW(i,j) ) / TSNOW 
          endif
          if (       INDEX_PFT(i,j,2) >= 1 .and. INDEX_PFT(i,j,2) <= 9 &
                .and. ( .not. ( Zsfc(i,j) >= HGT_DUST_Tibet &  ! exclude Tibetan Plateau
                          .and. LON_REAL(i,j) >= 65.0_RP * CONST_D2R .and. LON_REAL(i,j) <= 110.0_RP * CONST_D2R &
                          .and. LAT_REAL(i,j) >= 24.0_RP * CONST_D2R .and. LAT_REAL(i,j) <=  41.3_RP * CONST_D2R ) ) ) then
             CC(i,j) = CC(i,j) + FACT_LAND(i,j) * FRAC_PFT(i,j,2) * CEMFDU * ( TWATER -  GRWG(i,j) ) / TWATER &
                                                                           * ( TSNOW  - GRSNW(i,j) ) / TSNOW 
          endif

          if ( GRLAI(i,j) >= 0.0_RP ) then
             FPAR(i,j) = 1.0_RP - exp( -0.5_RP * GRLAI(i,j) ) 
          else
             FPAR(i,j) = 1.0_RP
          endif
          if ( FPAR(i,j) <= 0.25_RP ) then
             AEFF(i,j) = 1.0_RP - FPAR(i,j)
          else
             AEFF(i,j) = 0.0_RP
          endif
          
          EMITF(i,j) = AEFF(i,j) * CC(i,j) * ( WS10(i,j) - UTR ) * WS10(i,j)**2
          if ( CC(i,j) > 0.0_RP ) SKUP(i,j) = real( KUP(i,j), kind=RP )

       else
          
          EMITF(i,j) = 0.0_RP
          SKUP (i,j) = 0.0_RP

       endif
    enddo
    enddo
    
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          QEMIT(i,j,iq) = WSRC(iq) * EMITF(i,j) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
          do k = KS, KUPMAX
             if ( k <= KUP(i,j) ) then
                QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(i,j,iq)
             endif
          enddo
       enddo
       enddo
    enddo

    
    ! !instability 
    ! !-------------------------------
    ! do iq = 1, QA_AE_DU
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB( QTRC(:,:,:,iq),                                & ! [MODIFIED]
    !                                               ONEWQ(1:KA-1,:,:), DELP(:,:,:), GDPM(0:KA,:,:) ) ! [IN]
    ! enddo

    
    !incloud
    !-------------------------------
    !inside cloud
    if ( SW_CCN ) then
       do iq = 1, QA_AE_DU
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTIDU(k,i,j) * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    else
       do iq = 1, QA_AE_DU
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCDU * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    endif

    !outside cloud
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCTM1(k,i,j,iq) = QTRC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    
    ! !subcloud scavenging (wash out)
    ! !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    do k = KS, KE
    !       RAINN(k,i,j) = FRAIN(k,i,j) / VTR / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3 )
    !    enddo
    ! enddo
    ! enddo

    ! do iq = 1, QA_AE_DU
    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          DUSTN             = QTRCTM1(k,i,j,iq) / MD(iq) * DENS(k,i,j)
    !          DUSTRN (k,i,j,iq) = AA(iq) * PI * ( RADR + RADN(iq) )**2 * ( VTR - VTER(iq) ) * DUSTN * RAINN(k,i,j)
    !          QWTDEP (k,i,j,iq) = DUSTRN(k,i,j,iq) * MD(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
    !          QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP(k,i,j,iq)

    !          if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
    !             QTRCTMP(k,i,j,iq) = 0.0_RP
    !             DUSTRN (k,i,j,iq) = DUSTN / real(dt_AE,kind=RP)
    !             QWTDEP (k,i,j,iq) = QTRCTM1(k,i,j,iq)
    !          endif
    !       enddo
    !    enddo
    !    enddo
    ! enddo

    
    ! !incloud  scavenging (rain out)
    ! !-------------------------------
    ! do iq = 1, QA_AE_DU
    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          CLTORN (k,i,j,iq) = RNWR(k,i,j) * QTRCTM2(k,i,j,iq)
    !          QTRCTM2(k,i,j,iq) = QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq)
    !          QWTDEP (k,i,j,iq) = QWTDEP (k,i,j,iq) + CLTORN(k,i,j,iq)
    !          QTRCTM1(k,i,j,iq) = QTRCTMP(k,i,j,iq)
    !       enddo
    !    enddo
    !    enddo
    ! enddo

    
    ! !re-emission from rain
    ! !-------------------------------
    ! do iq = 1, QA_AE_DU
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER                         &
    !         ( QTRCTM1(:,:,:,iq), QWTDEP(:,:,:,iq),                       & ! [MODIFIED]
    !           WETF(:,:,iq),                                              & ! [OUT]
    !           DUSTRN(:,:,:,iq), MD(iq), CLTORN(:,:,:,iq),                & ! [IN]
    !           DENS(:,:,:), FRAIN(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          QTRC(k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
    !       enddo
    !    enddo
    !    enddo

    ! enddo
             
    
    ! !dry deposition
    ! !-------------------------------
    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !      ( QMTX,                   & ! [OUT]
    !        DENS, DELP, CDVE, dt_AE ) ! [IN]

    ! do iq = 1, QA_AE_DU
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF  (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    ! enddo


    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, QA_AE_DU
    !    VTERG(:,:,:) = VTER(iq)
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !         ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
    !           GDPM, DELP, DENS, VTERG, dt_AE  ) ! [IN]
    ! enddo

    
    !subcloud scavenging (wash out)
    !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    do k = KS, KE
    !       RAINN_CP     (k,i,j) = FRAIN_CP     (k,i,j) / VTR                / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3                )
    !       RAINN_MP_RAIN(k,i,j) = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
    !       RAINN_MP_SNOW(k,i,j) = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
    !    enddo
    ! enddo
    ! enddo

    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             DUSTN                    = QTRCTM1(k,i,j,iq) / MD(iq) * DENS(k,i,j)
             DUSTRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADN(iq) )**2 * abs( VTR                - VTER(iq) ) * DUSTN * RAINN_CP     (k,i,j)
             DUSTRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADN(iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(iq) ) * DUSTN * RAINN_MP_RAIN(k,i,j)
             DUSTRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADN(iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(iq) ) * DUSTN * RAINN_MP_SNOW(k,i,j)
             QWTDEP_CP     (k,i,j,iq) = DUSTRN_CP     (k,i,j,iq) * MD(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_RAIN(k,i,j,iq) = DUSTRN_MP_RAIN(k,i,j,iq) * MD(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_SNOW(k,i,j,iq) = DUSTRN_MP_SNOW(k,i,j,iq) * MD(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
                QTRCTMP(k,i,j,iq) = 0.0_RP
                DUSTRN (k,i,j,iq) = DUSTRN_CP(k,i,j,iq) + DUSTRN_MP_RAIN(k,i,j,iq) + DUSTRN_MP_SNOW(k,i,j,iq)
                DUSTRN_CP     (k,i,j,iq) = DUSTN / real(dt_AE,kind=RP) * DUSTRN_CP     (k,i,j,iq) / DUSTRN(k,i,j,iq)
                DUSTRN_MP_RAIN(k,i,j,iq) = DUSTN / real(dt_AE,kind=RP) * DUSTRN_MP_RAIN(k,i,j,iq) / DUSTRN(k,i,j,iq)
                DUSTRN_MP_SNOW(k,i,j,iq) = DUSTN / real(dt_AE,kind=RP) * DUSTRN_MP_SNOW(k,i,j,iq) / DUSTRN(k,i,j,iq)
                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * DUSTRN_CP     (k,i,j,iq) / DUSTRN(k,i,j,iq)
                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * DUSTRN_MP_RAIN(k,i,j,iq) / DUSTRN(k,i,j,iq)
                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * DUSTRN_MP_SNOW(k,i,j,iq) / DUSTRN(k,i,j,iq)
             endif
          enddo
       enddo
       enddo
    enddo

    
    !incloud  scavenging (rain out)
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN        (k,i,j,iq) = CLTORN_CP(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
             if ( QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq) < 0.0_RP ) then
                CLTORN_CP     (k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_CP     (k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
                CLTORN_MP_RAIN(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_RAIN(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
                CLTORN_MP_SNOW(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_SNOW(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
                QTRCTM2(k,i,j,iq) = 0.0_RP
             else
                QTRCTM2(k,i,j,iq) = QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq)
             endif
             QWTDEP_CP     (k,i,j,iq) = QWTDEP_CP     (k,i,j,iq) + CLTORN_CP     (k,i,j,iq)
             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq)
             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
             QTRCTM1       (k,i,j,iq) = QTRCTMP       (k,i,j,iq)
          enddo
       enddo
       enddo
    enddo


    !dry deposition
    !-------------------------------
    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
         ( QMTX,                   & ! [OUT]
           DENS, DELP, CDVE, dt_AE ) ! [IN]

    do iq = 1, QA_AE_DU
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
            ( QTRCTM1(:,:,:,iq),            & ! [MODIFIED]
              DRYF(:,:,iq),                 & ! [OUT]
              QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    enddo


    !gravitational settling
    !-------------------------------
    do iq = 1, QA_AE_DU
       VTERG(:,:,:) = VTER(iq)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
            ( QTRCTM1(:,:,:,iq),              & ! [MODIFIED]
              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
              GDPM, DELP, DENS, VTERG, dt_AE  ) ! [IN]
    enddo

    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + QTRCINR(k,i,j,iq)
             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + QTRCINS(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo
    

    !re-emission from rain : CP
    !-------------------------------
    do iq = 1, QA_AE_DU
       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
              WETF_CP(:,:,iq),                                              & ! [OUT]
              DUSTRN_CP(:,:,:,iq), MD(iq), CLTORN_CP(:,:,:,iq),             & ! [IN]
              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
    enddo

    
    !re-emission from rain & snow : MP
    !-------------------------------
    do iq = 1, QA_AE_DU
       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
       !        WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
       !        QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
       !        GDPM, DELP, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
       !        WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
       !        QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
       !        GDPM, DELP, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
       
       QTRCINR(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-rain(MP)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
              WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
              QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
              GDZM, DELZ, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCINR(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINR(k,i,j,iq), 0.0_RP )
          enddo
       enddo
       enddo
          
       QTRCINS(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-snow(MP)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
              WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
              QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
              GDZM, DELZ, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCINS(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINS(k,i,j,iq), 0.0_RP )
          enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCINP(k,i,j,iq) = QTRCINR(k,i,j,iq) + QTRCINS(k,i,j,iq)
          enddo
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          WETF(i,j,iq) = WETF_CP(i,j,iq) + WETF_MP_RAIN(i,j,iq) + WETF_MP_SNOW(i,j,iq)
       enddo
       enddo
    enddo

    
    !inside cloud + outside cloud 
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRC (k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
             QLOAD(k,i,j,iq) = QTRC   (k,i,j,iq)
          enddo
       enddo
       enddo
    enddo
             
    
    ! !dry deposition
    ! !-------------------------------
    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !      ( QMTX,                   & ! [OUT]
    !        DENS, DELP, CDVE, dt_AE ) ! [IN]

    ! do iq = 1, QA_AE_DU
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF  (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    ! enddo


    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, QA_AE_DU
    !    VTERG(:,:,:) = VTER(iq)
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !         ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
    !           GDPM, DELP, DENS, VTERG, dt_AE  ) ! [IN]
    ! enddo

    !================================================
    !end aerosol transport

    
    !flux 
    !================================================
    !Reset 
    !-------------------------------
    !             i j
    WETFT        (:,:) = 0.0_RP
    WETFT_CP     (:,:) = 0.0_RP
    WETFT_MP_RAIN(:,:) = 0.0_RP
    WETFT_MP_SNOW(:,:) = 0.0_RP
    DRYFT        (:,:) = 0.0_RP
    GRAVFT       (:,:) = 0.0_RP
    DEPFT        (:,:) = 0.0_RP
    EMIF25       (:,:) = 0.0_RP
    WETF25       (:,:) = 0.0_RP
    DRYF25       (:,:) = 0.0_RP
    GRVF25       (:,:) = 0.0_RP
    
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          WETFT        (i,j) = WETFT        (i,j) + WETF        (i,j,iq)
          WETFT_CP     (i,j) = WETFT_CP     (i,j) + WETF_CP     (i,j,iq)
          WETFT_MP_RAIN(i,j) = WETFT_MP_RAIN(i,j) + WETF_MP_RAIN(i,j,iq)
          WETFT_MP_SNOW(i,j) = WETFT_MP_SNOW(i,j) + WETF_MP_SNOW(i,j,iq)
          DRYFT        (i,j) = DRYFT        (i,j) + DRYF        (i,j,iq)
          GRAVFT       (i,j) = GRAVFT       (i,j) + GRAVF       (i,j,iq)
          DEPFT        (i,j) = DEPFT(i,j) + WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
          if ( iq <= 3 ) then
             EMIF25(i,j) = EMIF25(i,j) + EMITF(i,j) * WSRC(iq)
             WETF25(i,j) = WETF25(i,j) + WETF (i,j,iq)
             DRYF25(i,j) = DRYF25(i,j) + DRYF (i,j,iq)
             GRVF25(i,j) = GRVF25(i,j) + GRAVF(i,j,iq)
          endif
       enddo
       enddo
    enddo
    !================================================
    !end flux 


    !aerosol 
    !================================================
    !Reset 
    !-------------------------------
    !      k i j q
    QLDDT (:,:,:)   = 0.0_RP
    DCONT (:,:,:)   = 0.0_RP
    NUMCNT(:,:,:)   = 0.0_RP
    COLMS   (:,:,:) = 0.0_RP
    COLMST  (:,:)   = 0.0_RP
    RADDU (:,:,:)   = RADD(1)
    
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QLDDT (k,i,j)    = QLDDT(k,i,j) + QLOAD(k,i,j,iq)
             DCON  (k,i,j,iq) = QLOAD(k,i,j,iq) * DENS(k,i,j)
             NUMCON(k,i,j,iq) = DCON(k,i,j,iq) / MD(iq)
             DCONT (k,i,j)    = DCONT(k,i,j) + DCON(k,i,j,iq)
             NUMCNT(k,i,j)    = NUMCNT(k,i,j) + NUMCON(k,i,j,iq)
             COLMS   (i,j,iq) = COLMS(i,j,iq) + QLOAD(k,i,j,iq) * DELP(k,i,j) / GRAV
             if ( iq >= 2 ) then
                if ( NUMCON(k,i,j,iq) > NUMCON(k,i,j,iq-1) ) RADDU(k,i,j) = RADD(iq)
             endif
          enddo
          COLMST(i,j) = COLMST(i,j) + COLMS(i,j,iq)
       enddo
       enddo
       
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          PM25DU(k,i,j) = sum( DCON  (k,i,j,1:3) )
          NPM25D(k,i,j) = sum( NUMCON(k,i,j,1:3) )
          QPM25D(k,i,j) = PM25DU(k,i,j) / DENS(k,i,j)
       enddo
    enddo
    enddo
    !================================================
    !end aerosol     

    
    !optical parameters
    !================================================
    !      k i j w 2
    TAUDU   (:,:,:,:) = 0.0_RP
    TAUDUA  (:,:,:,:) = 0.0_RP
    TAUDUF  (:,:,  :) = 0.0_RP
    CEXTDU(:,:,:,:,:) = 0.0_RP
    CABSDU(:,:,:,:,:) = 0.0_RP
    CBAKDU(:,:,:,:,:) = 0.0_RP
    CEXTDF(:,:,:,:)   = 0.0_RP
    
    ! all sky
    !-------------------------------
    do iq = 1, QA_AE_DU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             TUSE = QLOAD(k,i,j,iq) * DELP(k,i,j) / GRAV
             do iw = 1, IWA
                TAUDU   (i,j,iw,1) = TAUDU   (i,j,iw,1) + TUSE * CEXTDP(iq,iw)
                TAUDUA  (i,j,iw,1) = TAUDUA  (i,j,iw,1) + TUSE * CABSDP(iq,iw)
                CEXTDU(k,i,j,iw,1) = CEXTDU(k,i,j,iw,1) + TUSE * CEXTDA(iq,iw) / DELZ(k,i,j)
                CABSDU(k,i,j,iw,1) = CABSDU(k,i,j,iw,1) + TUSE * CABSDA(iq,iw) / DELZ(k,i,j)
                CBAKDU(k,i,j,iw,1) = CBAKDU(k,i,j,iw,1) + TUSE * CBAKDA(iq,iw) / DELZ(k,i,j)
             enddo
             if ( iq <= 3 ) then
                TAUDUF  (i,j,1) = TAUDUF  (i,j,1) + TUSE * CEXTDP(iq,1)
                CEXTDF(k,i,j,1) = CEXTDF(k,i,j,1) + TUSE * CEXTDA(iq,1) / DELZ(k,i,j)
             endif
          enddo
       enddo
       enddo
    enddo

   
    !clear sky diagnoses
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       
       do k = KS, KE
          if ( TCLDF2(k,i,j) >= CCMAX ) then
             DCONTC(k,i,j)         = CONST_UNDEF !VMISS
             NUMCTC(k,i,j)         = CONST_UNDEF !VMISS
             CEXTDF(k,i,j,      2) = CONST_UNDEF !VMISS
             CEXTDU(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CABSDU(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CBAKDU(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
          else
             DCONTC(k,i,j)         = DCONT (k,i,j)
             NUMCTC(k,i,j)         = NUMCNT(k,i,j)
             CEXTDF(k,i,j,      2) = CEXTDF(k,i,j,      1)
             CEXTDU(k,i,j,1:IWA,2) = CEXTDU(k,i,j,1:IWA,1)
             CABSDU(k,i,j,1:IWA,2) = CABSDU(k,i,j,1:IWA,1)
             CBAKDU(k,i,j,1:IWA,2) = CBAKDU(k,i,j,1:IWA,1)
          endif
       enddo

       if ( CCOVMR(i,j) >= CCMAX ) then
          COLMTC(i,j)         = CONST_UNDEF !VMISS
          TAUDUF(i,j,      2) = CONST_UNDEF !VMISS
          TAUDU (i,j,1:IWA,2) = CONST_UNDEF !VMISS
          TAUDUA(i,j,1:IWA,2) = CONST_UNDEF !VMISS
       else
          COLMTC(i,j)         = COLMST(i,j)
          TAUDUF(i,j,      2) = TAUDUF(i,j,      1)
          TAUDU (i,j,1:IWA,2) = TAUDU (i,j,1:IWA,1)
          TAUDUA(i,j,1:IWA,2) = TAUDUA(i,j,1:IWA,1)
       endif
       
    enddo
    enddo
    !================================================
    !end optical parameters


#if defined(DEBUG) || defined(QUICKDEBUG)
    !checl OUTPUT
    !================================================
    check_2d_data(:,:, 1) = EMITF        (:,:)
    check_2d_data(:,:, 2) = WETFT        (:,:)
    check_2d_data(:,:, 3) = WETFT_CP     (:,:)
    check_2d_data(:,:, 4) = WETFT_MP_RAIN(:,:)
    check_2d_data(:,:, 5) = WETFT_MP_SNOW(:,:)
    check_2d_data(:,:, 6) = DRYFT        (:,:)
    check_2d_data(:,:, 7) = GRAVFT       (:,:)
    check_2d_data(:,:, 8) = DEPFT        (:,:)
    check_2d_data(:,:, 9) = EMIF25       (:,:)
    check_2d_data(:,:,10) = WETF25       (:,:)
    check_2d_data(:,:,11) = DRYF25       (:,:)
    check_2d_data(:,:,12) = GRVF25       (:,:)
    check_2d_data(:,:,13) = COLMST  (:,:)
    check_2d_data(:,:,14) = SKUP    (:,:)
    check_2d_data(:,:,15) = COLMS  (:,:,1)
    check_2d_data(:,:,16) = COLMS  (:,:,2)
    check_2d_data(:,:,17) = COLMS  (:,:,3)
    check_2d_data(:,:,18) = COLMS  (:,:,4)
    check_2d_data(:,:,19) = COLMS  (:,:,5)
    check_2d_data(:,:,20) = COLMS  (:,:,6)
    check_2d_data(:,:,21) = TAUDU   (:,:,1,1)
    check_2d_data(:,:,22) = TAUDUA  (:,:,1,1)
    check_2d_data(:,:,23) = TAUDUF  (:,:,  1)
    check_2d_data(:,:,24) = COLMTC  (:,:)
    check_2d_data(:,:,25) = TAUDU   (:,:,1,2)
    check_2d_data(:,:,26) = TAUDUA  (:,:,1,2)
    check_2d_data(:,:,27) = TAUDUF  (:,:,  2)

    call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 27, &
                            (/ 'EMITF(:,:)        ', 'WETFT(:,:)        ', 'WETFT_CP(:,:)     ', 'WETFT_MP_RAIN(:,:)', 'WETFT_MP_SNOW(:,:)', &
                               'DRYFT(:,:)        ', 'GRAVFT(:,:)       ', 'DEPFT(:,:)        ', 'EMIF25(:,:)       ', 'WETF25(:,:)       ', &
                               'DRYF25(:,:)       ', 'GRVF25(:,:)       ', 'COLMST(:,:)       ', 'SKUP(:,:)         ', 'COLMS(:,:,1)      ', &
                               'COLMS(:,:,2)      ', 'COLMS(:,:,3)      ', 'COLMS(:,:,4)      ', 'COLMS(:,:,5)      ', 'COLMS(:,:,6)      ', &
                               'TAUDU(:,:,1,1)    ', 'TAUDUA(:,:,1,1)   ', 'TAUDUF(:,:,1)     ', 'COLMTC(:,:)       ', 'TAUDU(:,:,1,2)    ', &
                               'TAUDUA(:,:,1,2)   ', 'TAUDUF(:,:,  2)   '/), &
                               check_2d_data(:,:,1:27) )

    check_3d_data(:,:,:, 1) = QLDDT (:,:,:)
    check_3d_data(:,:,:, 2) = DCONT (:,:,:)
    check_3d_data(:,:,:, 3) = NUMCNT(:,:,:)
    check_3d_data(:,:,:, 4) = QPM25D(:,:,:)
    check_3d_data(:,:,:, 5) = NPM25D(:,:,:)
    check_3d_data(:,:,:, 6) = DCON  (:,:,:,1)
    check_3d_data(:,:,:, 7) = DCON  (:,:,:,2)
    check_3d_data(:,:,:, 8) = DCON  (:,:,:,3)
    check_3d_data(:,:,:, 9) = DCON  (:,:,:,4)
    check_3d_data(:,:,:,10) = DCON  (:,:,:,5)
    check_3d_data(:,:,:,11) = DCON  (:,:,:,6)
    check_3d_data(:,:,:,12) = NUMCON(:,:,:,1)
    check_3d_data(:,:,:,13) = NUMCON(:,:,:,2)
    check_3d_data(:,:,:,14) = NUMCON(:,:,:,3)
    check_3d_data(:,:,:,15) = NUMCON(:,:,:,4)
    check_3d_data(:,:,:,16) = NUMCON(:,:,:,5)
    check_3d_data(:,:,:,17) = NUMCON(:,:,:,6)
    check_3d_data(:,:,:,18) = CEXTDU(:,:,:,1,1)
    check_3d_data(:,:,:,19) = CABSDU(:,:,:,1,1)
    check_3d_data(:,:,:,20) = CBAKDU(:,:,:,1,1)
    check_3d_data(:,:,:,21) = CEXTDF(:,:,:,  1)
    check_3d_data(:,:,:,22) = CEXTDU(:,:,:,2,1)
    check_3d_data(:,:,:,23) = CABSDU(:,:,:,2,1)
    check_3d_data(:,:,:,24) = CBAKDU(:,:,:,2,1)
    check_3d_data(:,:,:,25) = CEXTDU(:,:,:,3,1)
    check_3d_data(:,:,:,26) = CABSDU(:,:,:,3,1)
    check_3d_data(:,:,:,27) = CBAKDU(:,:,:,3,1)
    check_3d_data(:,:,:,28) = DCONTC(:,:,:)
    check_3d_data(:,:,:,29) = NUMCTC(:,:,:)
    check_3d_data(:,:,:,30) = CEXTDU(:,:,:,1,2)
    check_3d_data(:,:,:,31) = CABSDU(:,:,:,1,2)
    check_3d_data(:,:,:,32) = CBAKDU(:,:,:,1,2)
    check_3d_data(:,:,:,33) = CEXTDF(:,:,:,  2)
    check_3d_data(:,:,:,34) = CEXTDU(:,:,:,2,2)
    check_3d_data(:,:,:,35) = CABSDU(:,:,:,2,2)
    check_3d_data(:,:,:,36) = CBAKDU(:,:,:,2,2)
    check_3d_data(:,:,:,37) = CEXTDU(:,:,:,3,2)
    check_3d_data(:,:,:,38) = CABSDU(:,:,:,3,2)
    check_3d_data(:,:,:,39) = CBAKDU(:,:,:,3,2)

    call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 39, &
                            (/ 'QLDDT(:,:,:)     ', 'DCONT(:,:,:)     ', 'NUMCNT(:,:,:)    ', 'QPM25D(:,:,:)    ', 'NPM25D(:,:,:)    ', &
                               'DCON1            ', 'DCON2            ', 'DCON3            ', 'DCON4            ', 'DCON5            ', &
                               'DCON6            ', 'NUMCON1          ', 'NUMCON2          ', 'NUMCON3          ', 'NUMCON4          ', &
                               'NUMCON5          ', 'NUMCON6          ', 'CEXTDU(:,:,:,1,1)', 'CABSDU(:,:,:,1,1)', 'CBAKDU(:,:,:,1,1)', &
                               'CEXTDF(:,:,:,  1)', 'CEXTDU(:,:,:,2,1)', 'CABSDU(:,:,:,2,1)', 'CBAKDU(:,:,:,2,1)', 'CEXTDU(:,:,:,3,1)', &
                               'CABSDU(:,:,:,3,1)', 'CBAKDU(:,:,:,3,1)', 'DCONTC(:,:,:)    ', 'NUMCTC(:,:,:)    ', 'CEXTDU(:,:,:,1,2)', &
                               'CABSDU(:,:,:,1,2)', 'CBAKDU(:,:,:,1,2)', 'CEXTDF(:,:,:,  2)', 'CEXTDU(:,:,:,2,2)', 'CABSDU(:,:,:,2,2)', &
                               'CBAKDU(:,:,:,2,2)', 'CEXTDU(:,:,:,3,2)', 'CABSDU(:,:,:,3,2)', 'CBAKDU(:,:,:,3,2)' /), &
                            check_3d_data(:,:,:,1:39) )
    !================================================
    !end check OUTPUT
#endif

    
    !OUTPUT
    !================================================
    !flux
    !-------------------------------
    !                                   i j
    call FILE_HISTORY_in( EMITF        (:,:), 'EFDU_SPRINTARS',           'emission flux (dust)',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT        (:,:), 'WEFTDU_SPRINTARS',         'wet deposition flux (dust)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_CP     (:,:), 'WEFTDU_CP_SPRINTARS',      'wet deposition flux (dust;CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_RAIN(:,:), 'WEFTDU_MP_RAIN_SPRINTARS', 'wet deposition flux (dust;MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_SNOW(:,:), 'WEFTDU_MP_SNOW_SPRINTARS', 'wet deposition flux (dust;MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYFT        (:,:), 'DRFTDU_SPRINTARS',         'dry deposition flux (dust)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRAVFT       (:,:), 'GRFTDU_SPRINTARS',         'gravitational settling flux (dust)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DEPFT        (:,:), 'TDFTDU_SPRINTARS',         'total deposition flux (dust)',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMIF25       (:,:), 'EFDU25_SPRINTARS',         'emission flux (dust,PM2.5)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF25       (:,:), 'WFDU25_SPRINTARS',         'wet deposition flux (dust,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYF25       (:,:), 'DFDU25_SPRINTARS',         'dry deposition flux (dust,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRVF25       (:,:), 'GFDU25_SPRINTARS',         'gravitational settling flux (dust,PM2.5)', 'kg/m2/s', dim_type = 'XY'  )
   
    !aerosol
    !-------------------------------
    !                            k i j
    call FILE_HISTORY_in( QLDDT (:,:,:), 'QLDDT_SPRINTARS',  'mixing ratio (dust)',               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( DCONT (:,:,:), 'CONDT_SPRINTARS',  'mass concentration (dust)',         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCNT(:,:,:), 'NMCDT_SPRINTARS',  'number concentration (dust)',       '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMST  (:,:), 'COLDT_SPRINTARS',  'colomn loading (dust)',             'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( QPM25D(:,:,:), 'QPM25D_SPRINTARS', 'mixing ratio (dust,PM2.5)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NPM25D(:,:,:), 'NPM25D_SPRINTARS', 'number concentration (dust,PM2.5)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( SKUP    (:,:), 'SKUP_SPRINTARS',   'uplift level',                      'level', dim_type = 'XY'  )
    do iq = 1, QA_AE_DU
       write(cnumber,'(i2.2)') iq
       !                            k i j q
       call FILE_HISTORY_in( DCON  (:,:,:,iq), 'COND'//cnumber//'_SPRINTARS', 'mass concentration (dust'//cnumber//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( NUMCON(:,:,:,iq), 'NMCD'//cnumber//'_SPRINTARS', 'number concentration (dust'//cnumber//')', '1/m3',  dim_type = 'ZXY' )
       call FILE_HISTORY_in(  COLMS  (:,:,iq), 'COLD'//cnumber//'_SPRINTARS', 'column loading (dust'//cnumber//')',       'kg/m2', dim_type = 'XY'  )
    enddo
    
    !optical properties
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( TAUDU   (:,:,1,1), 'TAUDU_SPRINTARS',  'AOT (dust;550nm;allsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUDUA  (:,:,1,1), 'TAUDUA_SPRINTARS', 'AOT (dust;550nm;allsky;abs)',         '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUDUF  (:,:,  1), 'TAUDUF_SPRINTARS', 'AOT (dust;550nm;allsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,1,1), 'EXTDU1_SPRINTARS', 'extinction (dust;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CABSDU(:,:,:,1,1), 'ABSDU1_SPRINTARS', 'absorption (dust;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,1,1), 'BAKDU1_SPRINTARS', 'backscatter (dust;532nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTDF(:,:,:,  1), 'EXTDF1_SPRINTARS', 'extinction (dust;532nm;allsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,2,1), 'EXTDU2_SPRINTARS', 'extinction (dust;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CABSDU(:,:,:,2,1), 'ABSDU2_SPRINTARS', 'absorption (dust;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,2,1), 'BAKDU2_SPRINTARS', 'backscatter (dust;1064nm;allsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,3,1), 'EXTDU3_SPRINTARS', 'extinction (dust;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CABSDU(:,:,:,3,1), 'ABSDU3_SPRINTARS', 'absorption (dust;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,3,1), 'BAKDU3_SPRINTARS', 'backscatter (dust;355nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm

    !clear sky diagnoses
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( DCONTC(:,:,:),     'CONDTC_SPRINTARS',  'mass concentration (dust;clearsky)',    'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCTC(:,:,:),     'NMCDTC_SPRINTARS',  'number concentration (dust;clearsky)',  '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMTC  (:,:),     'COLDTC_SPRINTARS',  'column loading (dust;clearsky)',        'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUDU   (:,:,1,2), 'TAUDUC_SPRINTARS',  'AOT (dust;550nm;clearsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUDUA  (:,:,1,2), 'TAUDUAC_SPRINTARS', 'AOT (dust;550nm;clearsky;abs)',         '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUDUF  (:,:,  2), 'TAUDUFC_SPRINTARS', 'AOT (dust;550nm;clearsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,1,2), 'EXTDUC1_SPRINTARS', 'extinction (dust;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CABSDU(:,:,:,1,2), 'ABSDUC1_SPRINTARS', 'absorption (dust;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,1,2), 'BAKDUC1_SPRINTARS', 'backscatter (dust;532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTDF(:,:,:,  2), 'EXTDFC1_SPRINTARS', 'extinction (dust;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,2,2), 'EXTDUC2_SPRINTARS', 'extinction (dust;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CABSDU(:,:,:,2,2), 'ABSDUC2_SPRINTARS', 'absorption (dust;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,2,2), 'BAKDUC2_SPRINTARS', 'backscatter (dust;1064nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CEXTDU(:,:,:,3,2), 'EXTDUC3_SPRINTARS', 'extinction (dust;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CABSDU(:,:,:,3,2), 'ABSDUC3_SPRINTARS', 'absorption (dust;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CBAKDU(:,:,:,3,2), 'BAKDUC3_SPRINTARS', 'backscatter (dust;355nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm
    
    !================================================
    !end OUTPUT
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AERODU
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_aerodust
