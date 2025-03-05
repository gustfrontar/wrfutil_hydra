!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - sea salt
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_aerosalt
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
       PI   => CONST_PI,   & !< 
       GRAV => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
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
       NSA,                                           & !< total number of sea salt tracers
       IWA,                                           & !< number of wavelengths
       NRH,                                           & !< number of RH bin
       BRH,                                           & !< RH bin
       !VMISS,                                         & !< missing value
       NU,                                            & !< air viscosity [kg/m/s] = [Pa*s]
       RADR,                                          & !< rain droplet radius      [m]
       VTR,                                           & !< rain terminal velocity   [m/s]
       DENSR,                                         & !< rain droplet density     [kg/m3]
       DENSW,                                         & !< density of water         [kg/m3]
       THRESP                                           !< precipitation threshold  [kg/m2/s]
       
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public  :: ATMOS_PHY_AE_SPRINTARS_AEROSA_setup
  public  :: ATMOS_PHY_AE_SPRINTARS_AEROSA_finalize
  public  :: ATMOS_PHY_AE_SPRINTARS_AEROSA
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86 !< Sea salt aerosol source function : Monahan et al. (1986)
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! parameter
  ! variables
  real(RP), private, save :: FINCSA = 0.600_RP !< incloud coefficient
  real(RP), private, save :: DENSSF = 1.373_RP !< factor for sea salt density
  
  real(RP), private, allocatable, save :: QTRCINP(:,:,:,:) !< mixing ratio in prec(MP)

  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Sea salt setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA_setup( QA_AE_SA )
    implicit none
    integer, intent(in) :: QA_AE_SA
    integer :: ierr, iq
    character :: cnumber*2
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSA_NMAESA / &
         FINCSA, &
         DENSSF
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosalt",*) 'Setup'

    !--- read namelist
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSA_NMAESA, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosalt",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSA_NMAESA. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosalt",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSA_NMAESA, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSA_NMAESA)
    
    allocate( QTRCINP(KA,IA,JA,NSA) )
    QTRCINP(:,:,:,:) = 0.0_RP

    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP
    call FILE_HISTORY_in( zero_2d(:,:), 'EFSA_SPRINTARS',           'emission flux (salt)',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTSA_SPRINTARS',         'wet deposition flux (salt)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTSA_CP_SPRINTARS',      'wet deposition flux (salt;CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTSA_MP_RAIN_SPRINTARS', 'wet deposition flux (salt;MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WEFTSA_MP_SNOW_SPRINTARS', 'wet deposition flux (salt;MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'DRFTSA_SPRINTARS',         'dry deposition flux (salt)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'GRFTSA_SPRINTARS',         'gravitational settling flux (salt)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'TDFTSA_SPRINTARS',         'total deposition flux (salt)',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'EFSA25_SPRINTARS',         'emission flux (salt,PM2.5)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'WFSA25_SPRINTARS',         'wet deposition flux (salt,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'DFSA25_SPRINTARS',         'dry deposition flux (salt,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:), 'GFSA25_SPRINTARS',         'gravitational settling flux (salt,PM2.5)', 'kg/m2/s', dim_type = 'XY'  )
    
    !aerosol
    !-------------------------------
    !                            k i j
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDNT_SPRINTARS',  'mixing ratio (salt)',               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONNT_SPRINTARS',  'mass concentration (salt)',         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCNT_SPRINTARS',  'number concentration (salt)',       '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLNT_SPRINTARS',  'column loading (salt)',             'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QPM25N_SPRINTARS', 'mixing ratio (salt,PM2.5)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NPM25N_SPRINTARS', 'number concentration (salt,PM2.5)', '1/m3',  dim_type = 'ZXY' )
    do iq = 1, QA_AE_SA
       write(cnumber,'(i2.2)') iq
       !                            k i j q
       call FILE_HISTORY_in( zero_3d(:,:,:), 'CONN'//cnumber//'_SPRINTARS', 'mass concentration (salt;salt'//cnumber//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCN'//cnumber//'_SPRINTARS', 'number concentration (salt;salt'//cnumber//')', '1/m3',  dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_2d(:,:),   'COLN'//cnumber//'_SPRINTARS', 'column loading (salt;salt'//cnumber//')',       'kg/m2', dim_type = 'XY'  )
    enddo

    !optical properties
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUNN_SPRINTARS',  'AOT (salt;550nm;allsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUNNF_SPRINTARS', 'AOT (salt;550nm;allsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNN1_SPRINTARS', 'extinction (salt;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNN1_SPRINTARS', 'backscatter (salt;532nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNF1_SPRINTARS', 'extinction (salt;532nm;allsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNN2_SPRINTARS', 'extinction (salt;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNN2_SPRINTARS', 'backscatter (salt;1064nm;allsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNN3_SPRINTARS', 'extinction (salt;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNN3_SPRINTARS', 'backscatter (salt;355nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm

    !clear sky diagnoses
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONNTC_SPRINTARS',  'mass concentration (salt;clearsky)',    'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCNTC_SPRINTARS',  'number concentration (salt;clearsky)',  '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'COLNTC_SPRINTARS',  'column loading (salt;clearsky)',        'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUNNC_SPRINTARS',  'AOT (salt;550nm;clearsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_2d(:,:  ), 'TAUNNFC_SPRINTARS', 'AOT (salt;550nm;clearsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNNC1_SPRINTARS', 'extinction (salt;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNNC1_SPRINTARS', 'backscatter (salt;532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNFC1_SPRINTARS', 'extinction (salt;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNNC2_SPRINTARS', 'extinction (salt;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNNC2_SPRINTARS', 'backscatter (salt;1064nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTNNC3_SPRINTARS', 'extinction (salt;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKNNC3_SPRINTARS', 'backscatter (salt;355nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA_setup
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> SPRINTARS Sea salt finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA_finalize
    implicit none
  
    deallocate( QTRCINP )
    deallocate( zero_3d )
    deallocate( zero_2d )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA_finalize

  !-----------------------------------------------------------------------------
  !> SPRINTARS Sea salt aerosol source function : Monahan et al. (1986) only bubble (not spume) 
  function ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86( RW, FRH ) result( SF )
    real(RP), intent(in) :: RW  !< sea salt particle radius [um]
    real(RP), intent(in) :: FRH !< growth factor
    real(RP)             :: SF  !< source function
    real(RP)             :: B   !< working area
    !----------------------------------------------------
    ! B = ( 0.38_RP - log10(RW) ) / 0.65_RP 
    ! SF = ( 1.0_RP + 0.057_RP * RW**1.05_RP ) * 10.0_RP**( 1.19_RP * exp(-B**2) ) / FRH**3
    B = ( 0.38_RP - log10(RW*FRH) ) / 0.65_RP 
    SF = ( 1.0_RP + 0.057_RP * (RW*FRH)**1.05_RP ) * 10.0_RP**( 1.19_RP * exp(-B**2) ) / FRH**3
  end function ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Sea salt
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA &
       ( QTRC,           & ! [MODIFIED]
         OUTQLD,         & ! [OUT]
         NUMCNT,         & ! [OUT]
         RADSA,          & ! [OUT]
         PM25SA,         & ! [OUT]
         WATRSA,         & ! [OUT]
         EMITF,          & ! [OUT]
         TAUNN,          & ! [OUT]
         TAUNNF,         & ! [OUT]
         CEXTNN,         & ! [OUT]
         CBAKNN,         & ! [OUT]
         CEXTNF,         & ! [OUT]
         TAUNND,         & ! [OUT]
         CEXTND,         & ! [OUT]
         CDVE,           & ! [IN]
         CCOVMR,         & ! [IN]
         CCMAX,          & ! [IN]
         WS10,           & ! [IN]
         GOICE,          & ! [IN]
         RH,             & ! [IN]
         DENS,           & ! [IN]
         DENS1,          & ! [IN]
         FRAIN_MP_RAIN1, & ! [IN]
         FRAIN_MP_SNOW1, & ! [IN]
         FRAIN,          & ! [IN]
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
         TCLDF,          & ! [IN]
         TCLDF2,         & ! [IN]
         KUP,            & ! [IN]
         KUPMAX,         & ! [IN]
         DPKUP,          & ! [IN]
         GDPM,           & ! [IN]
         DELP,           & ! [IN]
         GDZM,           & ! [IN]
         DELZ,           & ! [IN]
         ONEWQ,          & ! [IN]
         RNWR,           & ! [IN]
         RNWR_CP,        & ! [IN]
         RNWR_MP_RAIN,   & ! [IN]
         RNWR_MP_SNOW,   & ! [IN]
         FACT_OCEAN,     & ! [IN]
         ACTISA,         & ! [IN]
         SW_CCN,         & ! [IN]
         dt_AE,          & ! [IN]
         QA_AE_SA        ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE                      !< time step
    integer,  intent(in)    :: QA_AE_SA                   !< total number of sea salt tracers
    
    real(RP), intent(in)    :: CDVE               (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
    real(RP), intent(in)    :: CCOVMR             (IA,JA) !< cloud cover (max-ran)
    real(RP), intent(in)    :: CCMAX                      !< maximum cloud fraction for all/clear sky
    real(RP), intent(in)    :: WS10               (IA,JA) !< wind at 10m (land+ocean) [m/s]
    real(RP), intent(in)    :: GOICE              (IA,JA) !< sea ice
    real(RP), intent(in)    :: RH            (KA,  IA,JA) !< relative humidity
    real(RP), intent(in)    :: DENS          (KA,  IA,JA) !< air density
    real(RP), intent(in)    :: DENS1         (KA,  IA,JA) !< air density   on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA) !< precipitation flux (liq    ) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA) !< precipitation flux (    ice) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN         (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP + MP
    real(RP), intent(in)    :: FRAIN_CP      (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP
    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA) !< precipitation flux (liq    ) : MP
    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA) !< precipitation flux (    ice) : MP
    real(RP), intent(in)    :: VTR_MP_RAIN   (KA,  IA,JA) !< terminal velocity of rain (microphysic scheme) [m/s]
    real(RP), intent(in)    :: VTR_MP_SNOW   (KA,  IA,JA) !< terminal velocity of snow (microphysic scheme) [m/s]
    real(RP), intent(in)    :: RADR_MP_RAIN  (KA,  IA,JA) !< effective radius of rain (microphysic scheme) [m]
    real(RP), intent(in)    :: RADR_MP_SNOW  (KA,  IA,JA) !< effective radius of snow (microphysic scheme) [m]
    real(RP), intent(in)    :: RAINN_CP      (KA,  IA,JA) !< rain number concentration (CP)      [1/m3]
    real(RP), intent(in)    :: RAINN_MP_RAIN (KA,  IA,JA) !< rain number concentration (MP;RAIN) [1/m3]
    real(RP), intent(in)    :: RAINN_MP_SNOW (KA,  IA,JA) !< rain number concentration (MP;SNOW) [1/m3]
    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA) !< cloud fraction (CP+MP)
    real(RP), intent(in)    :: TCLDF2        (KA,  IA,JA) !< cloud fraction (CP+MP) for CCOVMR
    integer,  intent(in)    :: KUP                (IA,JA) !< uplift level
    integer,  intent(in)    :: KUPMAX                     !< maximum uplift level
    real(RP), intent(in)    :: DPKUP              (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA) !< pressure at half levels
    real(RP), intent(in)    :: DELP          (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA) !< height at half levels
    real(RP), intent(in)    :: DELZ          (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
    logical,  intent(in)    :: ONEWQ         (KA-1,IA,JA) !< use for instability
    real(RP), intent(in)    :: RNWR          (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: FACT_OCEAN         (IA,JA) !< ocean factor
    real(RP), intent(in)    :: ACTISA        (KA,  IA,JA) !< ratio of activated aerosols [0-1]
    logical,  intent(in)    :: SW_CCN                     !< switch for CCN scheme (true->use ACTI)
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC  (KA,IA,JA,QA_AE_SA) !< ratio of sea salt tracer mass to total mass
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: OUTQLD(KA,IA,JA,QA_AE_SA) !< QTRC for radiation
    real(RP), intent(out)   :: NUMCNT(KA,IA,JA)          !< number concentration (total)
    real(RP), intent(out)   :: RADSA (KA,IA,JA)          !< mode radius
    real(RP), intent(out)   :: PM25SA(KA,IA,JA)          !< submicron mass concentration
    real(RP), intent(out)   :: WATRSA(KA,IA,JA)          !< aerosol water
    real(RP), intent(out)   :: EMITF    (IA,JA,QA_AE_SA) !< emission flux [kg/m2/s]
    real(RP), intent(out)   :: TAUNN    (IA,JA,IWA,2)    !< AOT
    real(RP), intent(out)   :: TAUNNF   (IA,JA,    2)    !< AOT (fine mode)
    real(RP), intent(out)   :: CEXTNN(KA,IA,JA,IWA,2)    !< extinction  coefficient
    real(RP), intent(out)   :: CBAKNN(KA,IA,JA,IWA,2)    !< backscatter coefficient
    real(RP), intent(out)   :: CEXTNF(KA,IA,JA,    2)    !< extinction  coefficient (fine mode)
    real(RP), intent(out)   :: TAUNND   (IA,JA)          !< AOT (dry)
    real(RP), intent(out)   :: CEXTND(KA,IA,JA)          !< extinction coefficient (salt;dry)
    
    ! internal work
    !-------------------------
    !parameter 
    real(RP) :: PI2                            !< 4.0_RP / 3.0_RP * PI * 1.0E-18_RP
    real(RP) :: DENSSAX                        !< sea salt density modified by DENSSF
    real(RP) :: RADN            (QA_AE_SA)     !< aerosol radius of number concentration
    real(RP) :: MSA             (QA_AE_SA)     !< sea salt particle mass
    real(RP) :: VTER            (QA_AE_SA)     !< terminal velocity
    real(RP) :: FRH    (KA,IA,JA)              !< growth factor at each grid
    real(RP) :: CEXTP  (KA,IA,JA,QA_AE_SA,IWA) !< extinction  cross section (550/ 440/870nm)
    real(RP) :: CEXTA  (KA,IA,JA,QA_AE_SA,IWA) !< extinction  cross section (532/1064/355nm)
    real(RP) :: CBAKA  (KA,IA,JA,QA_AE_SA,IWA) !< backscatter cross section (532/1064/355nm)
    real(RP) :: BRH1, BRH2, BRH3               !< use for interpolation
    
    !emission
    real(RP) :: SUMO      (IA,JA)              !< use for emission
    real(RP) :: SUME      (IA,JA)              !< use for emission
    real(RP) :: SSE       (IA,JA)              !< use for emission
    real(RP) :: DIVH                           !< use for emission
    real(RP) :: RW                             !< use for emission
    real(RP) :: RW2                            !< use for emission
    real(RP) :: SUM1                           !< use for emission
    real(RP) :: SUM2                           !< use for emission

    !instability
    real(RP) :: NEWQTRC                        !< use for instability
    
    !incloud
    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_SA)     !< outside cloud water [kg/kg]
    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_SA)     !< inside  cloud water [kg/kg]
    
    !subcloud scavenging (wash out)
    ! real(RP) :: RAINN         (KA,IA,JA)          !< rain number concentration         [1/m3]
    ! real(RP) :: RAINN_CP      (KA,IA,JA)
    ! real(RP) :: RAINN_MP_RAIN (KA,IA,JA)
    ! real(RP) :: RAINN_MP_SNOW (KA,IA,JA)
    real(RP) :: SALTN                             !< number density of sea salt tracer [1/m3]
    real(RP) :: SALTRN        (KA,IA,JA,QA_AE_SA) !< number density of sea salt tracer entering into raindrops [1/s/m3] <=> collision frequency between sea salt and rain
    real(RP) :: SALTRN_CP     (KA,IA,JA,QA_AE_SA)
    real(RP) :: SALTRN_MP_RAIN(KA,IA,JA,QA_AE_SA)
    real(RP) :: SALTRN_MP_SNOW(KA,IA,JA,QA_AE_SA)
    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_SA) !< mixing ratio by wet deposition [kg/kg]
    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_SA)
    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_SA)
    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_SA)
    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_SA) !< sea salt tracers that do not cllide with rain particles [kg/kg]
    
    !incloud scavenging (rain out)
    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_SA) !< mixing ratio moving cloud to rain particles [kg/kg]
    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_SA)
    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_SA)
    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_SA)
    
    !re-emission from rain
    real(RP) :: WETF           (IA,JA,QA_AE_SA) !< wet deposition flux
    real(RP) :: WETF_CP        (IA,JA,QA_AE_SA)
    real(RP) :: WETF_MP_RAIN   (IA,JA,QA_AE_SA)
    real(RP) :: WETF_MP_SNOW   (IA,JA,QA_AE_SA)
    real(RP) :: QTRCINR     (KA,IA,JA,QA_AE_SA) !< mixing ratio in rain(MP)
    real(RP) :: QTRCINS     (KA,IA,JA,QA_AE_SA) !< mixing ratio in snow(MP)

    !dry deposition
    real(RP) :: QMTX(KA,IA,JA,-1:1)         !< implicit matrix of QTRC
    real(RP) :: DRYF   (IA,JA,QA_AE_SA)     !< dry deposition flux
    
    !gravitational settling
    real(RP) :: VTERG(KA,IA,JA)             !< terminal velocity for gravitational settling
    real(RP) :: QLOAD(KA,IA,JA,QA_AE_SA)    !< mixing ratio for radiation
    real(RP) :: GRAVF   (IA,JA,QA_AE_SA)    !< gravitational settling flux

    !file output (flux)
    real(RP) :: EMITFT       (IA,JA)           !< emission flux         [kg/m2/s]
    real(RP) :: WETFT        (IA,JA)           !< wet deposition flux (total)
    real(RP) :: WETFT_CP     (IA,JA)           !< wet deposition flux (total;CP)
    real(RP) :: WETFT_MP_RAIN(IA,JA)           !< wet deposition flux (total;MP(RAIN))
    real(RP) :: WETFT_MP_SNOW(IA,JA)           !< wet deposition flux (total;MP(SNOW))
    real(RP) :: DRYFT        (IA,JA)           !< dry deposition flux (total)
    real(RP) :: GRAVFT       (IA,JA)           !< gravitational settling flux (total)
    real(RP) :: DEPFT        (IA,JA)           !< total deposition flux
    real(RP) :: EMIF25       (IA,JA)           !< emission flux (PM2.5) [kg/m2/s]
    real(RP) :: WETF25       (IA,JA)           !< wet deposition flux (PM2.5)
    real(RP) :: DRYF25       (IA,JA)           !< dry deposition flux (PM2.5)
    real(RP) :: GRVF25       (IA,JA)           !< gravitational settling flux (PM2.5)
    
    !file output (aerosol)
    real(RP) :: QLDNT  (KA,IA,JA)              !< mixing ratio (total) [kg/kg]
    real(RP) :: DCON   (KA,IA,JA,QA_AE_SA)     !< mass concentration
    real(RP) :: DCONT  (KA,IA,JA)              !< mass concentration (total)
    real(RP) :: NUMCON (KA,IA,JA,QA_AE_SA)     !< number concentration
    real(RP) :: COLMS     (IA,JA,QA_AE_SA)     !< column mass loading
    real(RP) :: COLMST    (IA,JA)              !< column mass loading (total)
    real(RP) :: QPM25N (KA,IA,JA)              !< mass mixing ratio (PM2.5)
    real(RP) :: NPM25N (KA,IA,JA)              !< number concentration (PM2.5)
    
    !file output (optical parameter)
    real(RP) :: DCONTC (KA,IA,JA)              !< mass concentraion (total,clear sky)
    real(RP) :: NUMCTC (KA,IA,JA)              !< number concentration (total,clear sky)
    real(RP) :: COLMTC    (IA,JA)              !< column mass loading (total,clear sky)
    real(RP) :: TUSE                           !< use for AOT

    ! parameter
    !-------------------------
    real(RP), parameter :: &
         RADD(NSA)         &    !< sea salt particle radius [m]
         !    salt01,     salt02,     salt03,     salt04
         = (/ 1.78E-7_RP, 5.62E-7_RP, 1.78E-6_RP, 5.62E-6_RP /)
    real(RP), parameter :: &
         RADDM(NSA+1)      &    !< sea salt particle radius (median) [m]
         = (/ 1.00E-7_RP, 3.16E-7_RP, 1.00E-6_RP, 3.16E-6_RP, 1.00E-5_RP /)
    real(RP), parameter :: &
         FRHI(NRH)         &    !< growth factor
         = (/ 1.0000E+0_RP, 1.0694E+0_RP, 1.2756E+0_RP, 1.9875E+0_RP, & ! RH = (0.00,0.50,0.70,0.80) 
              2.3769E+0_RP, 2.8788E+0_RP, 3.7650E+0_RP, 4.6906E+0_RP /) ! RH = (0.90,0.95,0.98,0.99)
    real(RP), parameter :: &
         CEXTNP(NSA,NRH,IWA) &  !<  extinction cross section (550/ 440/870nm)
         !             salt01,      salt02,      salt03,      salt04
         = reshape( (/ 3.311E+3_RP, 1.693E+3_RP, 4.402E+2_RP, 1.307E+2_RP, & ! RH = 0.00, IWA =  550nm
                       3.876E+3_RP, 1.945E+3_RP, 4.993E+2_RP, 1.485E+2_RP, & ! RH = 0.50, IWA =  550nm
                       6.116E+3_RP, 2.529E+3_RP, 7.074E+2_RP, 2.097E+2_RP, & ! RH = 0.70, IWA =  550nm
                       2.249E+4_RP, 5.564E+3_RP, 1.668E+3_RP, 5.044E+2_RP, & ! RH = 0.80, IWA =  550nm
                       3.699E+4_RP, 8.061E+3_RP, 2.366E+3_RP, 7.185E+2_RP, & ! RH = 0.90, IWA =  550nm
                       5.818E+4_RP, 1.180E+4_RP, 3.433E+3_RP, 1.050E+3_RP, & ! RH = 0.95, IWA =  550nm
                       8.920E+4_RP, 1.960E+4_RP, 5.797E+3_RP, 1.789E+3_RP, & ! RH = 0.98, IWA =  550nm
                       1.105E+5_RP, 2.992E+4_RP, 8.954E+3_RP, 2.765E+3_RP, & ! RH = 0.99, IWA =  550nm
                       4.946E+3_RP, 1.401E+3_RP, 4.339E+2_RP, 1.294E+2_RP, & ! RH = 0.00, IWA =  440nm
                       5.730E+3_RP, 1.589E+3_RP, 4.948E+2_RP, 1.474E+2_RP, & ! RH = 0.50, IWA =  440nm
                       8.739E+3_RP, 2.207E+3_RP, 6.964E+2_RP, 2.089E+2_RP, & ! RH = 0.70, IWA =  440nm
                       2.728E+4_RP, 5.691E+3_RP, 1.650E+3_RP, 5.017E+2_RP, & ! RH = 0.80, IWA =  440nm
                       3.989E+4_RP, 8.002E+3_RP, 2.333E+3_RP, 7.155E+2_RP, & ! RH = 0.90, IWA =  440nm
                       5.344E+4_RP, 1.151E+4_RP, 3.389E+3_RP, 1.046E+3_RP, & ! RH = 0.95, IWA =  440nm
                       6.982E+4_RP, 1.927E+4_RP, 5.765E+3_RP, 1.782E+3_RP, & ! RH = 0.98, IWA =  440nm
                       9.640E+4_RP, 2.951E+4_RP, 8.889E+3_RP, 2.757E+3_RP, & ! RH = 0.99, IWA =  440nm
                       8.996E+2_RP, 2.358E+3_RP, 4.696E+2_RP, 1.336E+2_RP, & ! RH = 0.00, IWA =  870nm
                       1.112E+3_RP, 2.681E+3_RP, 5.402E+2_RP, 1.523E+2_RP, & ! RH = 0.50, IWA =  870nm
                       2.019E+3_RP, 3.767E+3_RP, 7.486E+2_RP, 2.148E+2_RP, & ! RH = 0.70, IWA =  870nm
                       1.045E+4_RP, 7.265E+3_RP, 1.715E+3_RP, 5.102E+2_RP, & ! RH = 0.80, IWA =  870nm
                       2.006E+4_RP, 8.756E+3_RP, 2.424E+3_RP, 7.277E+2_RP, & ! RH = 0.90, IWA =  870nm
                       3.899E+4_RP, 1.154E+4_RP, 3.514E+3_RP, 1.061E+3_RP, & ! RH = 0.95, IWA =  870nm
                       8.948E+4_RP, 1.991E+4_RP, 5.933E+3_RP, 1.802E+3_RP, & ! RH = 0.98, IWA =  870nm
                       1.536E+5_RP, 3.137E+4_RP, 9.101E+3_RP, 2.787E+3_RP /), & ! RH = 0.99, IWA =  870nm
                    (/ NSA, NRH, IWA /) )
    real(RP), parameter :: &
         CEXTNA(NSA,NRH,IWA) &  !<  extinction cross section (532/1064/355nm)
          !            salt01,      salt02,      salt03,      salt04
         = reshape( (/ 3.547E+3_RP, 1.623E+3_RP, 4.392E+2_RP, 1.305E+2_RP, & ! RH = 0.00, IWA =  532nm
                       4.143E+3_RP, 1.866E+3_RP, 4.989E+2_RP, 1.483E+2_RP, & ! RH = 0.50, IWA =  532nm
                       6.496E+3_RP, 2.434E+3_RP, 7.057E+2_RP, 2.096E+2_RP, & ! RH = 0.70, IWA =  532nm
                       2.336E+4_RP, 5.602E+3_RP, 1.666E+3_RP, 5.040E+2_RP, & ! RH = 0.80, IWA =  532nm
                       3.786E+4_RP, 8.076E+3_RP, 2.362E+3_RP, 7.180E+2_RP, & ! RH = 0.90, IWA =  532nm
                       5.825E+4_RP, 1.177E+4_RP, 3.425E+3_RP, 1.050E+3_RP, & ! RH = 0.95, IWA =  532nm
                       8.635E+4_RP, 1.952E+4_RP, 5.791E+3_RP, 1.788E+3_RP, & ! RH = 0.98, IWA =  532nm
                       1.068E+5_RP, 2.987E+4_RP, 8.946E+3_RP, 2.764E+3_RP, & ! RH = 0.99, IWA =  532nm
                       4.610E+2_RP, 2.031E+3_RP, 4.652E+2_RP, 1.354E+2_RP, & ! RH = 0.00, IWA = 1064nm
                       5.819E+2_RP, 2.329E+3_RP, 5.260E+2_RP, 1.541E+2_RP, & ! RH = 0.50, IWA = 1064nm
                       1.113E+3_RP, 3.454E+3_RP, 7.611E+2_RP, 2.170E+2_RP, & ! RH = 0.70, IWA = 1064nm
                       6.620E+3_RP, 8.614E+3_RP, 1.752E+3_RP, 5.139E+2_RP, & ! RH = 0.80, IWA = 1064nm
                       1.349E+4_RP, 1.092E+4_RP, 2.456E+3_RP, 7.306E+2_RP, & ! RH = 0.90, IWA = 1064nm
                       2.790E+4_RP, 1.325E+4_RP, 3.556E+3_RP, 1.068E+3_RP, & ! RH = 0.95, IWA = 1064nm
                       7.150E+4_RP, 1.950E+4_RP, 5.986E+3_RP, 1.811E+3_RP, & ! RH = 0.98, IWA = 1064nm
                       1.397E+5_RP, 3.091E+4_RP, 9.197E+3_RP, 2.795E+3_RP, & ! RH = 0.99, IWA = 1064nm
                       6.593E+3_RP, 1.489E+3_RP, 4.288E+2_RP, 1.285E+2_RP, & ! RH = 0.00, IWA =  355nm
                       7.604E+3_RP, 1.649E+3_RP, 4.884E+2_RP, 1.467E+2_RP, & ! RH = 0.50, IWA =  355nm
                       1.112E+4_RP, 2.403E+3_RP, 6.873E+2_RP, 2.078E+2_RP, & ! RH = 0.70, IWA =  355nm
                       2.762E+4_RP, 5.552E+3_RP, 1.625E+3_RP, 4.998E+2_RP, & ! RH = 0.80, IWA =  355nm
                       3.497E+4_RP, 7.822E+3_RP, 2.306E+3_RP, 7.130E+2_RP, & ! RH = 0.90, IWA =  355nm
                       4.196E+4_RP, 1.130E+4_RP, 3.371E+3_RP, 1.043E+3_RP, & ! RH = 0.95, IWA =  355nm
                       6.220E+4_RP, 1.902E+4_RP, 5.726E+3_RP, 1.777E+3_RP, & ! RH = 0.98, IWA =  355nm
                       9.868E+4_RP, 2.920E+4_RP, 8.834E+3_RP, 2.751E+3_RP /), & ! RH = 0.99, IWA =  355nm
                    (/ NSA, NRH, IWA /) )
    real(RP), parameter :: &
         CBAKNA(NSA,NRH,IWA) &  !< backscatter cross section (532/1064/355nm)
          !            salt01,      salt02,      salt03,      salt04
         = reshape( (/ 4.495E+1_RP, 1.481E+2_RP, 2.388E+1_RP, 7.126E+0_RP, & ! RH = 0.00, IWA =  532nm
                       5.306E+1_RP, 1.142E+2_RP, 2.646E+1_RP, 1.011E+1_RP, & ! RH = 0.50, IWA =  532nm
                       7.805E+1_RP, 7.670E+1_RP, 4.243E+1_RP, 2.131E+1_RP, & ! RH = 0.70, IWA =  532nm
                       2.306E+2_RP, 3.398E+2_RP, 1.205E+2_RP, 6.109E+1_RP, & ! RH = 0.80, IWA =  532nm
                       3.364E+2_RP, 6.173E+2_RP, 1.765E+2_RP, 8.173E+1_RP, & ! RH = 0.90, IWA =  532nm
                       4.390E+2_RP, 8.679E+2_RP, 2.601E+2_RP, 1.113E+2_RP, & ! RH = 0.95, IWA =  532nm
                       6.202E+2_RP, 1.041E+3_RP, 5.097E+2_RP, 2.074E+2_RP, & ! RH = 0.98, IWA =  532nm
                       1.837E+3_RP, 1.409E+3_RP, 8.587E+2_RP, 3.340E+2_RP, & ! RH = 0.99, IWA =  532nm
                       2.419E+1_RP, 2.860E+1_RP, 4.056E+1_RP, 9.290E+0_RP, & ! RH = 0.00, IWA = 1064nm
                       2.672E+1_RP, 3.132E+1_RP, 3.347E+1_RP, 9.858E+0_RP, & ! RH = 0.50, IWA = 1064nm
                       3.163E+1_RP, 4.212E+1_RP, 5.062E+1_RP, 1.748E+1_RP, & ! RH = 0.70, IWA = 1064nm
                       7.122E+1_RP, 6.471E+1_RP, 1.176E+2_RP, 3.934E+1_RP, & ! RH = 0.80, IWA = 1064nm
                       1.458E+2_RP, 7.867E+1_RP, 1.294E+2_RP, 6.234E+1_RP, & ! RH = 0.90, IWA = 1064nm
                       2.712E+2_RP, 1.865E+2_RP, 1.632E+2_RP, 9.711E+1_RP, & ! RH = 0.95, IWA = 1064nm
                       6.373E+2_RP, 8.955E+2_RP, 3.316E+2_RP, 1.750E+2_RP, & ! RH = 0.98, IWA = 1064nm
                       1.108E+3_RP, 2.107E+3_RP, 6.180E+2_RP, 2.715E+2_RP, & ! RH = 0.99, IWA = 1064nm
                       1.002E+2_RP, 2.046E+2_RP, 3.632E+1_RP, 7.082E+0_RP, & ! RH = 0.00, IWA =  355nm
                       1.145E+2_RP, 1.563E+2_RP, 3.976E+1_RP, 1.033E+1_RP, & ! RH = 0.50, IWA =  355nm
                       1.529E+2_RP, 1.549E+2_RP, 6.098E+1_RP, 2.517E+1_RP, & ! RH = 0.70, IWA =  355nm
                       2.657E+2_RP, 4.143E+2_RP, 1.306E+2_RP, 5.909E+1_RP, & ! RH = 0.80, IWA =  355nm
                       3.115E+2_RP, 4.407E+2_RP, 2.057E+2_RP, 8.673E+1_RP, & ! RH = 0.90, IWA =  355nm
                       5.797E+2_RP, 5.579E+2_RP, 3.365E+2_RP, 1.357E+2_RP, & ! RH = 0.95, IWA =  355nm
                       2.683E+3_RP, 1.160E+3_RP, 6.290E+2_RP, 2.276E+2_RP, & ! RH = 0.98, IWA =  355nm
                       6.999E+3_RP, 2.105E+3_RP, 1.006E+3_RP, 2.915E+2_RP /), & ! RH = 0.99, IWA =  355nm
                    (/ NSA, NRH, IWA /) )
    real(RP), parameter :: &
         AA(NSA)           &    !< collision efficiency
         !    salt01,      salt02,      salt03,      salt04
         = (/ 1.000E-4_RP, 2.050E-4_RP, 0.105E+0_RP, 0.829E+0_RP /)
    real(RP), parameter :: DENSSA =  2.20E+3_RP !< sea salt density
    real(RP), parameter :: SIGMA  =  1.20E+0_RP !< GSD of log-normal distribution
    integer,  parameter :: NDIV   = 20          !< numerical integral method : composite Simpson's rule N = NDIV
    
    ! others
    !-------------------------
    integer   :: k, i, j, iq, iw
    integer   :: IRH
    integer   :: l
    character :: cnumber*2
    real(RP) :: check_2d_data   (IA,JA,100)
    real(RP) :: check_3d_data(KA,IA,JA,100)

    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - sea salt'
        
    !setup parameter 
    !================================================
    !PI2, DENSSAX, RADN, MD, VTER
    !-------------------------------
    PI2 = 4.0_RP / 3.0_RP * PI * 1.0E-18_RP
    DENSSAX = DENSSA * DENSSF
    do iq = 1, QA_AE_SA
       RADN(iq) = RADD(iq) / exp( 3.0_RP * log(SIGMA)**2 )
       MSA (iq) = 4.0_RP / 3.0_RP * PI * RADN(iq)**3 * DENSSA * exp( 4.5_RP * log(SIGMA)**2 )
       VTER(iq) = 2.0_RP * DENSSA * RADN(iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA)**2 )
    enddo

    !CEXTP, CEXTA, CBAKA
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          
          if (                               RH(k,i,j) < BRH(2) ) then ! RH =  0 - 50 %
             IRH = 1
          elseif ( RH(k,i,j) >= BRH(2) .and. RH(k,i,j) < BRH(3) ) then ! RH = 50 - 70 %
             IRH = 2
          elseif ( RH(k,i,j) >= BRH(3) .and. RH(k,i,j) < BRH(4) ) then ! RH = 70 - 80 %
             IRH = 3
          elseif ( RH(k,i,j) >= BRH(4) .and. RH(k,i,j) < BRH(5) ) then ! RH = 80 - 90 %
             IRH = 4
          elseif ( RH(k,i,j) >= BRH(5) .and. RH(k,i,j) < BRH(6) ) then ! RH = 90 - 95 %
             IRH = 5
          elseif ( RH(k,i,j) >= BRH(6) .and. RH(k,i,j) < BRH(7) ) then ! RH = 95 - 98 %
             IRH = 6
          elseif ( RH(k,i,j) >= BRH(7) .and. RH(k,i,j) < BRH(8) ) then ! RH = 98 - 99 %
             IRH = 7
          else                                                         ! RH = 99 -    %
             IRH = 8           
          endif
          
          if ( IRH == 8 ) then
             FRH(k,i,j) = FRHI(IRH)
             CEXTP(k,i,j,1:QA_AE_SA,1:IWA) = CEXTNP(1:QA_AE_SA,IRH,1:IWA)
             CEXTA(k,i,j,1:QA_AE_SA,1:IWA) = CEXTNA(1:QA_AE_SA,IRH,1:IWA)
             CBAKA(k,i,j,1:QA_AE_SA,1:IWA) = CBAKNA(1:QA_AE_SA,IRH,1:IWA)
          else
             BRH1 =  RH(k,i,j) - BRH(IRH)
             BRH2 = BRH(IRH+1) -  RH(k,i,j)
             BRH3 = BRH(IRH+1) - BRH(IRH)
             FRH(k,i,j) = ( BRH1 * FRHI(IRH+1) + BRH2 * FRHI(IRH) ) / BRH3
             do iw = 1, IWA
                do iq = 1, QA_AE_SA
                   CEXTP(k,i,j,iq,iw) = ( BRH1 * CEXTNP(iq,IRH+1,iw) + BRH2 * CEXTNP(iq,IRH,iw) ) / BRH3
                   CEXTA(k,i,j,iq,iw) = ( BRH1 * CEXTNA(iq,IRH+1,iw) + BRH2 * CEXTNA(iq,IRH,iw) ) / BRH3
                   CBAKA(k,i,j,iq,iw) = ( BRH1 * CBAKNA(iq,IRH+1,iw) + BRH2 * CBAKNA(iq,IRH,iw) ) / BRH3
                enddo
             enddo
          endif

       enddo
    enddo
    enddo
    !================================================
    !end setup parameter


    !aerosol transport 
    !================================================
    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_SA
       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
    enddo

    
    !emission (Monahan et al.,1986); numerical integral method : composite Simpson's rule N = NDIV( = 20)
    !-------------------------------
    EMITF(:,:,:) = 0.0_RP

    do iq = 1, QA_AE_SA
       
       SUMO(:,:) = 0.0_RP
       SUME(:,:) = 0.0_RP
       SSE (:,:) = 0.0_RP
       DIVH = ( RADDM(iq+1) - RADDM(iq) ) * 1.0E+6_RP / real( NDIV, kind = RP )
       do l = 1, NDIV-1, 2
          RW  = RADDM(iq) * 1.0E+6_RP + real( l, kind = RP ) * DIVH
          RW2 = ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86( RW, FRHI(4) ) ! RH = 0.80
          do j = JS, JE
          do i = IS, IE
             SUMO(i,j) = SUMO(i,j) + 2.0_RP * PI2 * DENSSAX * WS10(i,j)**3.41_RP * RW2
          enddo
          enddo
       enddo
       
       do l = 2, NDIV-2, 2
          RW  = RADDM(iq) * 1.0E+6_RP + real( l, kind = RP ) * DIVH
          RW2 = ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86( RW, FRHI(4) ) ! RH = 0.80
          do j = JS, JE
          do i = IS, IE
             SUME(i,j) = SUME(i,j) + 2.0_RP * PI2 * DENSSAX * WS10(i,j)**3.41_RP * RW2
          enddo
          enddo
       enddo
       
       RW  = ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86( RADDM(iq  )*1.0E+6_RP, FRHI(4) ) ! RH = 0.80
       RW2 = ATMOS_PHY_AE_SPRINTARS_AEROSA_SSASF_M86( RADDM(iq+1)*1.0E+6_RP, FRHI(4) ) ! RH = 0.80
       do j = JS, JE
       do i = IS, IE
          SUM1 = 2.0_RP * PI2 * DENSSAX * WS10(i,j)**3.41_RP * RW
          SUM2 = 2.0_RP * PI2 * DENSSAX * WS10(i,j)**3.41_RP * RW2
          if ( FACT_OCEAN(i,j) > 0.10_RP ) then
             SSE(i,j) = ( SUM1 + SUM2 + 4.0_RP * SUMO(i,j) + 2.0_RP * SUME(i,j) ) * DIVH / 3.0_RP * ( 1.0_RP - GOICE(i,j) ) * FACT_OCEAN(i,j)
          endif
          EMITF(i,j,iq) = EMITF(i,j,iq) + SSE(i,j)
       enddo
       enddo

       do j = JS, JE
       do i = IS, IE
          do k = KS, KUPMAX
             if ( k <= KUP(i,j) ) then
                QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + SSE(i,j) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
             endif
          enddo
       enddo
       enddo
          
    enddo

    
    ! !instability
    ! !-------------------------------
    ! do iq = 1, QA_AE_SA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB( QTRC(:,:,:,iq),                                & ! [MODIFIED]
    !                                               ONEWQ(1:KA-1,:,:), DELP(:,:,:), GDPM(0:KA,:,:) ) ! [IN]
    ! enddo

    
    !incloud
    !-------------------------------
    !inside cloud
    if ( SW_CCN ) then
       do iq = 1, QA_AE_SA
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTISA(k,i,j) * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    else
       do iq = 1, QA_AE_SA
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCSA * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    endif

    !outside cloud
    do iq = 1, QA_AE_SA
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

    ! do iq = 1, QA_AE_SA
    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          SALTN             = QTRCTM1(k,i,j,iq) / MSA(iq) * DENS(k,i,j)
    !          SALTRN (k,i,j,iq) = AA(iq) * PI * ( RADR + RADN(iq) )**2 * ( VTR - VTER(iq) ) * SALTN * RAINN(k,i,j)
    !          QWTDEP (k,i,j,iq) = SALTRN(k,i,j,iq) * MSA(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
    !          QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP(k,i,j,iq)
             
    !          if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
    !             QTRCTMP(k,i,j,iq) = 0.0_RP
    !             SALTRN (k,i,j,iq) = SALTN / real(dt_AE,kind=RP)
    !             QWTDEP (k,i,j,iq) = QTRCTM1(k,i,j,iq)
    !          endif
    !       enddo
    !    enddo
    !    enddo
    ! enddo
    

    ! !incloud  scavenging (rain out)
    ! !-------------------------------
    ! do iq = 1, QA_AE_SA
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
    ! do iq = 1, QA_AE_SA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER                         &
    !         ( QTRCTM1(:,:,:,iq), QWTDEP(:,:,:,iq),                       & ! [MODIFIED]
    !           WETF(:,:,iq),                                              & ! [OUT]
    !           SALTRN(:,:,:,iq), MSA(iq), CLTORN(:,:,:,iq),               & ! [IN]
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

    ! do iq = 1, QA_AE_SA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF  (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    ! enddo

    
    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, QA_AE_SA
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

    do iq = 1, QA_AE_SA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             SALTN = QTRCTM1(k,i,j,iq) / MSA(iq) * DENS(k,i,j)
             SALTRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADN(iq) )**2 * abs( VTR                - VTER(iq) ) * SALTN * RAINN_CP     (k,i,j)
             SALTRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADN(iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(iq) ) * SALTN * RAINN_MP_RAIN(k,i,j)
             SALTRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADN(iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(iq) ) * SALTN * RAINN_MP_SNOW(k,i,j)
             QWTDEP_CP     (k,i,j,iq) = SALTRN_CP     (k,i,j,iq) * MSA(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_RAIN(k,i,j,iq) = SALTRN_MP_RAIN(k,i,j,iq) * MSA(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_SNOW(k,i,j,iq) = SALTRN_MP_SNOW(k,i,j,iq) * MSA(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
                QTRCTMP(k,i,j,iq) = 0.0_RP
                SALTRN (k,i,j,iq) = SALTRN_CP(k,i,j,iq) + SALTRN_MP_RAIN(k,i,j,iq) + SALTRN_MP_SNOW(k,i,j,iq)
                SALTRN_CP     (k,i,j,iq) = SALTN / real(dt_AE,kind=RP) * SALTRN_CP     (k,i,j,iq) / SALTRN(k,i,j,iq)
                SALTRN_MP_RAIN(k,i,j,iq) = SALTN / real(dt_AE,kind=RP) * SALTRN_MP_RAIN(k,i,j,iq) / SALTRN(k,i,j,iq)
                SALTRN_MP_SNOW(k,i,j,iq) = SALTN / real(dt_AE,kind=RP) * SALTRN_MP_SNOW(k,i,j,iq) / SALTRN(k,i,j,iq)
                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * SALTRN_CP     (k,i,j,iq) / SALTRN(k,i,j,iq)
                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * SALTRN_MP_RAIN(k,i,j,iq) / SALTRN(k,i,j,iq)
                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * SALTRN_MP_SNOW(k,i,j,iq) / SALTRN(k,i,j,iq)
             endif
          enddo
       enddo
       enddo
    enddo

    
    !incloud  scavenging (rain out)
    !-------------------------------
    do iq = 1, QA_AE_SA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
             CLTORN        (k,i,j,iq) = CLTORN_CP(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
             if ( QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq) < 0.0_RP ) then
                CLTORN_CP     (k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_CP     (k,i,j) / (RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
                CLTORN_MP_RAIN(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_RAIN(k,i,j) / (RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
                CLTORN_MP_SNOW(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_SNOW(k,i,j) / (RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
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

    do iq = 1, QA_AE_SA
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
            ( QTRCTM1(:,:,:,iq),            & ! [MODIFIED]
              DRYF(:,:,iq),                 & ! [OUT]
              QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    enddo

    
    !gravitational settling
    !-------------------------------
    do iq = 1, QA_AE_SA
       VTERG(:,:,:) = VTER(iq)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
            ( QTRCTM1(:,:,:,iq),              & ! [MODIFIED]
              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
              GDPM, DELP, DENS, VTERG, dt_AE  ) ! [IN]
    enddo

    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_SA
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
    do iq = 1, QA_AE_SA
       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
              WETF_CP(:,:,iq),                                              & ! [OUT]
              SALTRN_CP(:,:,:,iq), MSA(iq), CLTORN_CP(:,:,:,iq),            & ! [IN]
              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
    enddo

    
    !re-emission from rain & snow : MP
    !-------------------------------
    do iq = 1, QA_AE_SA
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
    do iq = 1, QA_AE_SA
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

    ! do iq = 1, QA_AE_SA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF  (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
    ! enddo

    
    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, QA_AE_SA
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
    EMITFT       (:,:) = 0.0_RP
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
    do iq = 1, QA_AE_SA
       do j = JS, JE
       do i = IS, IE
          EMITFT       (i,j) = EMITFT       (i,j) + EMITF       (i,j,iq)
          WETFT        (i,j) = WETFT        (i,j) + WETF        (i,j,iq)
          WETFT_CP     (i,j) = WETFT_CP     (i,j) + WETF_CP     (i,j,iq)
          WETFT_MP_RAIN(i,j) = WETFT_MP_RAIN(i,j) + WETF_MP_RAIN(i,j,iq)
          WETFT_MP_SNOW(i,j) = WETFT_MP_SNOW(i,j) + WETF_MP_SNOW(i,j,iq)
          DRYFT        (i,j) = DRYFT        (i,j) + DRYF        (i,j,iq)
          GRAVFT       (i,j) = GRAVFT       (i,j) + GRAVF       (i,j,iq)
          DEPFT        (i,j) = DEPFT(i,j) + WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
          if ( iq <= 2 ) then
             EMIF25(i,j) = EMIF25(i,j) + EMITF(i,j,iq)
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
    OUTQLD(:,:,:,:) = 0.0_RP
    QLDNT (:,:,:)   = 0.0_RP
    DCONT (:,:,:)   = 0.0_RP
    NUMCNT(:,:,:)   = 0.0_RP
    COLMS   (:,:,:) = 0.0_RP
    COLMST  (:,:)   = 0.0_RP
    WATRSA(:,:,:)   = 0.0_RP
    RADSA (:,:,:)   = RADD(1)
    
    !-------------------------------
    do iq = 1, QA_AE_SA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             OUTQLD(k,i,j,iq) = OUTQLD(k,i,j,iq) + QLOAD(k,i,j,iq)
             QLDNT (k,i,j)    = QLDNT (k,i,j)    + QLOAD(k,i,j,iq)
             DCON  (k,i,j,iq) = QLOAD (k,i,j,iq) * DENS (k,i,j)
             DCONT (k,i,j)    = DCONT (k,i,j)    + DCON(k,i,j,iq)
             NUMCON(k,i,j,iq) = DCON  (k,i,j,iq) / MSA(iq)
             NUMCNT(k,i,j)    = NUMCNT(k,i,j)    + NUMCON(k,i,j,iq)
             COLMS   (i,j,iq) = COLMS   (i,j,iq) + QLOAD (k,i,j,iq) * DELP(k,i,j) / GRAV
             WATRSA(k,i,j)    = WATRSA(k,i,j)    + ( FRH(k,i,j)**3 - 1.0_RP ) * OUTQLD(k,i,j,iq) * DENSW / DENSSA
             if ( iq >= 2 ) then
                if ( NUMCON(k,i,j,iq) > NUMCON(k,i,j,iq-1) ) RADSA(k,i,j) = RADD(iq)
             endif
          enddo
          COLMST(i,j) = COLMST(i,j) + COLMS(i,j,iq)
       enddo
       enddo
    enddo
    
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          PM25SA(k,i,j) = sum( DCON  (k,i,j,1:2) )
          NPM25N(k,i,j) = sum( NUMCON(k,i,j,1:2) )
          QPM25N(k,i,j) = PM25SA(k,i,j) / DENS(k,i,j)
       enddo
    enddo
    enddo
    !================================================
    !end aerosol     

    
    !optical parameters
    !================================================
    !Reset 
    !-------------------------------
    !      k i j w 2
    TAUNN   (:,:,:,:) = 0.0_RP
    TAUNNF  (:,:,  :) = 0.0_RP
    TAUNND  (:,:)     = 0.0_RP
    CEXTNN(:,:,:,:,:) = 0.0_RP
    CBAKNN(:,:,:,:,:) = 0.0_RP
    CEXTNF(:,:,:,  :) = 0.0_RP
    CEXTND(:,:,:)     = 0.0_RP

    ! all sky
    !-------------------------------
    do iq = 1, QA_AE_SA
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             TUSE = QLOAD(k,i,j,iq) * DELP(k,i,j) / GRAV
             CEXTND(k,i,j) = CEXTND(k,i,j) + TUSE * CEXTNA(iq,1,1) / DELZ(k,i,j)
             do iw = 1, IWA
                TAUNN   (i,j,iw,1) = TAUNN   (i,j,iw,1) + TUSE * CEXTP(k,i,j,iq,iw)
                CEXTNN(k,i,j,iw,1) = CEXTNN(k,i,j,iw,1) + TUSE * CEXTA(k,i,j,iq,iw) / DELZ(k,i,j)
                CBAKNN(k,i,j,iw,1) = CBAKNN(k,i,j,iw,1) + TUSE * CBAKA(k,i,j,iq,iw) / DELZ(k,i,j)
             enddo
             if ( iq <= 2 ) then
                TAUNNF  (i,j,1) = TAUNNF  (i,j,1) + TUSE * CEXTP(k,i,j,iq,1)
                CEXTNF(k,i,j,1) = CEXTNF(k,i,j,1) + TUSE * CEXTA(k,i,j,iq,1) / DELZ(k,i,j)
             endif
             if ( iq <= 4 ) then
                TAUNND(i,j) = TAUNND(i,j) + TUSE * CEXTNP(iq,1,1)
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
             CEXTNF(k,i,j,      2) = CONST_UNDEF !VMISS
             CEXTNN(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CBAKNN(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
          else
             DCONTC(k,i,j)         = DCONT (k,i,j)
             NUMCTC(k,i,j)         = NUMCNT(k,i,j)
             CEXTNF(k,i,j,      2) = CEXTNF(k,i,j,      1)
             CEXTNN(k,i,j,1:IWA,2) = CEXTNN(k,i,j,1:IWA,1)
             CBAKNN(k,i,j,1:IWA,2) = CBAKNN(k,i,j,1:IWA,1)
          endif
       enddo

       if ( CCOVMR(i,j) >= CCMAX ) then
          COLMTC(i,j)         = CONST_UNDEF !VMISS
          TAUNNF(i,j,      2) = CONST_UNDEF !VMISS
          TAUNN (i,j,1:IWA,2) = CONST_UNDEF !VMISS
       else
          COLMTC(i,j)         = COLMST(i,j)
          TAUNNF(i,j,      2) = TAUNNF(i,j,      1)
          TAUNN (i,j,1:IWA,2) = TAUNN (i,j,1:IWA,1)
       endif
       
    enddo
    enddo
    !================================================
    !end optical parameters

    
#if defined(DEBUG) || defined(QUICKDEBUG)
    !checl OUTPUT
    !================================================
    check_2d_data(:,:, 1) = EMITFT       (:,:)
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
    check_2d_data(:,:,14) = COLMS   (:,:,1)
    check_2d_data(:,:,15) = COLMS   (:,:,2)
    check_2d_data(:,:,16) = COLMS   (:,:,3)
    check_2d_data(:,:,17) = COLMS   (:,:,4)
    check_2d_data(:,:,18) = TAUNN   (:,:,1,1)
    check_2d_data(:,:,19) = TAUNNF  (:,:,  1)
    check_2d_data(:,:,20) = COLMTC  (:,:)
    check_2d_data(:,:,21) = TAUNN   (:,:,1,2)
    check_2d_data(:,:,22) = TAUNNF  (:,:,  2)
    
    call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 22, &
                            (/ 'EMITFT       ', 'WETFT        ', 'WETFT_CP     ', 'WETFT_MP_RAIN', 'WETFT_MP_SNOW', &
                               'DRYFT        ', 'GRAVFT       ', 'DEPFT        ', 'EMIF25       ', 'WETF25       ', &
                               'DRYF25       ', 'GRVF25       ', 'COLMST       ', 'COLMS1       ', 'COLMS2       ', &
                               'COLMS3       ', 'COLMS4       ', 'TAUNN11      ', 'TAUNNF1      ', 'COLMTC       ', &
                               'TAUNN12      ', 'TAUNNF2      ' /), &
                               check_2d_data(:,:,1:22) )
    
    check_3d_data(:,:,:, 1) = QLDNT (:,:,:)
    check_3d_data(:,:,:, 2) = DCONT (:,:,:)
    check_3d_data(:,:,:, 3) = NUMCNT(:,:,:)
    check_3d_data(:,:,:, 4) = QPM25N(:,:,:)
    check_3d_data(:,:,:, 5) = NPM25N(:,:,:)
    check_3d_data(:,:,:, 6) = DCON  (:,:,:,1)
    check_3d_data(:,:,:, 7) = DCON  (:,:,:,2)
    check_3d_data(:,:,:, 8) = DCON  (:,:,:,3)
    check_3d_data(:,:,:, 9) = DCON  (:,:,:,4)
    check_3d_data(:,:,:,10) = NUMCON(:,:,:,1)
    check_3d_data(:,:,:,11) = NUMCON(:,:,:,2)
    check_3d_data(:,:,:,12) = NUMCON(:,:,:,3)
    check_3d_data(:,:,:,13) = NUMCON(:,:,:,4)
    check_3d_data(:,:,:,14) = CEXTNN(:,:,:,1,1)
    check_3d_data(:,:,:,15) = CBAKNN(:,:,:,1,1)
    check_3d_data(:,:,:,16) = CEXTNF(:,:,:,  1)
    check_3d_data(:,:,:,17) = CEXTNN(:,:,:,2,1)
    check_3d_data(:,:,:,18) = CBAKNN(:,:,:,2,1)
    check_3d_data(:,:,:,19) = CEXTNN(:,:,:,3,1)
    check_3d_data(:,:,:,20) = CBAKNN(:,:,:,3,1)
    check_3d_data(:,:,:,21) = DCONTC(:,:,:)
    check_3d_data(:,:,:,22) = NUMCTC(:,:,:)
    check_3d_data(:,:,:,23) = CEXTNN(:,:,:,1,2)
    check_3d_data(:,:,:,24) = CBAKNN(:,:,:,1,2)
    check_3d_data(:,:,:,25) = CEXTNF(:,:,:,  2)
    check_3d_data(:,:,:,26) = CEXTNN(:,:,:,2,2)
    check_3d_data(:,:,:,27) = CBAKNN(:,:,:,2,2)
    check_3d_data(:,:,:,28) = CEXTNN(:,:,:,3,2)
    check_3d_data(:,:,:,29) = CBAKNN(:,:,:,3,2)
    
    call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 29, &
                            (/ 'QLDNT            ', 'DCONT            ', 'NUMCNT           ', 'QPM25N           ', 'NPM25N           ', &
                               'DCON1            ', 'DCON2            ', 'DCON3            ', 'DCON4            ', 'NUMCON1          ', &
                               'NUMCON2          ', 'NUMCON3          ', 'NUMCON4          ', 'CEXTNN11         ', 'CBAKNN11         ', &
                               'CEXTNF(:,:,:,1)  ', 'CEXTNN(:,:,:,2,1)', 'CBAKNN(:,:,:,2,1)', 'CEXTNN(:,:,:,3,1)', 'CBAKNN(:,:,:,3,1)', &
                               'DCONTC(:,:,:)    ', 'NUMCTC(:,:,:)    ', 'CEXTNN(:,:,:,1,2)', 'CBAKNN(:,:,:,1,2)', 'CEXTNF(:,:,:,2)  ', &
                               'CEXTNN(:,:,:,2,2)', 'CBAKNN(:,:,:,2,2)', 'CEXTNN(:,:,:,3,2)', 'CBAKNN(:,:,:,3,2)'  /), &
                            check_3d_data(:,:,:,1:29) )
    !================================================
    !end check OUTPUT
#endif

    
    !OUTPUT
    !================================================
    !flux
    !-------------------------------
    !                                   i j
    call FILE_HISTORY_in( EMITFT       (:,:), 'EFSA_SPRINTARS',           'emission flux (salt)',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT        (:,:), 'WEFTSA_SPRINTARS',         'wet deposition flux (salt)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_CP     (:,:), 'WEFTSA_CP_SPRINTARS',      'wet deposition flux (salt;CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_RAIN(:,:), 'WEFTSA_MP_RAIN_SPRINTARS', 'wet deposition flux (salt;MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_SNOW(:,:), 'WEFTSA_MP_SNOW_SPRINTARS', 'wet deposition flux (salt;MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYFT        (:,:), 'DRFTSA_SPRINTARS',         'dry deposition flux (salt)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRAVFT       (:,:), 'GRFTSA_SPRINTARS',         'gravitational settling flux (salt)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DEPFT        (:,:), 'TDFTSA_SPRINTARS',         'total deposition flux (salt)',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMIF25       (:,:), 'EFSA25_SPRINTARS',         'emission flux (salt,PM2.5)',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF25       (:,:), 'WFSA25_SPRINTARS',         'wet deposition flux (salt,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYF25       (:,:), 'DFSA25_SPRINTARS',         'dry deposition flux (salt,PM2.5)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRVF25       (:,:), 'GFSA25_SPRINTARS',         'gravitational settling flux (salt,PM2.5)', 'kg/m2/s', dim_type = 'XY'  )
    
    !aerosol
    !-------------------------------
    !                            k i j
    call FILE_HISTORY_in( QLDNT (:,:,:), 'QLDNT_SPRINTARS',  'mixing ratio (salt)',               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( DCONT (:,:,:), 'CONNT_SPRINTARS',  'mass concentration (salt)',         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCNT(:,:,:), 'NMCNT_SPRINTARS',  'number concentration (salt)',       '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMST  (:,:), 'COLNT_SPRINTARS',  'column loading (salt)',             'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( QPM25N(:,:,:), 'QPM25N_SPRINTARS', 'mixing ratio (salt,PM2.5)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NPM25N(:,:,:), 'NPM25N_SPRINTARS', 'number concentration (salt,PM2.5)', '1/m3',  dim_type = 'ZXY' )
    do iq = 1, QA_AE_SA
       write(cnumber,'(i2.2)') iq
       !                            k i j q
       call FILE_HISTORY_in( DCON  (:,:,:,iq), 'CONN'//cnumber//'_SPRINTARS', 'mass concentration (salt;salt'//cnumber//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( NUMCON(:,:,:,iq), 'NMCN'//cnumber//'_SPRINTARS', 'number concentration (salt;salt'//cnumber//')', '1/m3',  dim_type = 'ZXY' )
       call FILE_HISTORY_in( COLMS   (:,:,iq), 'COLN'//cnumber//'_SPRINTARS', 'column loading (salt;salt'//cnumber//')',       'kg/m2', dim_type = 'XY'  )
    enddo

    !optical properties
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( TAUNN   (:,:,1,1), 'TAUNN_SPRINTARS',  'AOT (salt;550nm;allsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUNNF  (:,:,  1), 'TAUNNF_SPRINTARS', 'AOT (salt;550nm;allsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,1,1), 'EXTNN1_SPRINTARS', 'extinction (salt;532nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,1,1), 'BAKNN1_SPRINTARS', 'backscatter (salt;532nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTNF(:,:,:,  1), 'EXTNF1_SPRINTARS', 'extinction (salt;532nm;allsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,2,1), 'EXTNN2_SPRINTARS', 'extinction (salt;1064nm;allsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,2,1), 'BAKNN2_SPRINTARS', 'backscatter (salt;1064nm;allsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,3,1), 'EXTNN3_SPRINTARS', 'extinction (salt;355nm;allsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,3,1), 'BAKNN3_SPRINTARS', 'backscatter (salt;355nm;allsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm

    !clear sky diagnoses
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( DCONTC(:,:,:),     'CONNTC_SPRINTARS',  'mass concentration (salt;clearsky)',    'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCTC(:,:,:),     'NMCNTC_SPRINTARS',  'number concentration (salt;clearsky)',  '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMTC  (:,:),     'COLNTC_SPRINTARS',  'column loading (salt;clearsky)',        'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUNN   (:,:,1,2), 'TAUNNC_SPRINTARS',  'AOT (salt;550nm;clearsky)',             '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( TAUNNF  (:,:,  2), 'TAUNNFC_SPRINTARS', 'AOT (salt;550nm;clearsky;fine)',        '1',      dim_type = 'XY'  ) ! 550nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,1,2), 'EXTNNC1_SPRINTARS', 'extinction (salt;532nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,1,2), 'BAKNNC1_SPRINTARS', 'backscatter (salt;532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTNF(:,:,:,  2), 'EXTNFC1_SPRINTARS', 'extinction (salt;532nm;clearsky;fine)', '1/m',    dim_type = 'ZXY' ) ! 532nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,2,2), 'EXTNNC2_SPRINTARS', 'extinction (salt;1064nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,2,2), 'BAKNNC2_SPRINTARS', 'backscatter (salt;1064nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !1064nm
    call FILE_HISTORY_in( CEXTNN(:,:,:,3,2), 'EXTNNC3_SPRINTARS', 'extinction (salt;355nm;clearsky)',      '1/m',    dim_type = 'ZXY' ) ! 355nm
    call FILE_HISTORY_in( CBAKNN(:,:,:,3,2), 'BAKNNC3_SPRINTARS', 'backscatter (salt;355nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' ) ! 355nm
    !================================================
    !end OUTPUT

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSA
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_aerosalt
