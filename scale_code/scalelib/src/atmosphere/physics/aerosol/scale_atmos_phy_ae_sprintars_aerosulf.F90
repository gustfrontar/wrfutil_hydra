!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - sulfate
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_aerosulf
  !-----------------------------------------------------------------------------
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
       PRC_abort
  
  use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE

  use scale_const, only: &
       PI     => CONST_PI,   & !< pi
       GRAV   => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       MOLAIR => CONST_Mdry, & !< mass weight (dry air) = 28.966 [g/mol]
       CONST_UNDEF
       
  use scale_file_history, only: &
       FILE_HISTORY_in

  use scale_file_cartesC, only: &
       FILE_CARTESC_read
  
  use scale_atmos_phy_ae_SPRINTARS_common, only: &
       ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA,       &
       ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY,       &
       ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST,          &
       ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL,        &
       ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z, &
       ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC,          &
       NSU,                                           & !< total number of sulfate tracers
       IWA,                                           & !< number of wavelengths
       NRH,                                           & !< number of RH bin
       BRH,                                           & !< RH bin
       !VMISS,                                         & !< missing value
       NU,                                            & !< air viscosity [kg/m/s] = [Pa*s]
       RADR,                                          & !< rain droplet radius      [m]
       VTR,                                           & !< rain terminal velocity   [m/s]
       DENSR,                                         & !< rain droplet density     [kg/m3]
       DENSW,                                         & !< density of water         [kg/m3]
       AVOG,                                          & !< Avogadro's number
       THRESP                                           !< precipitation threshold  [kg/m2/s]
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSU_setup
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSU_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_AEROSU
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_PHY_AE_SPRINTARS_AEROSU_RSO2AQ
  private :: ATMOS_PHY_AE_SPRINTARS_AEROSU_QSSA
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! parameter
  ! variables
  real(RP), private, save :: FINCSU =    0.5E+0_RP !< incloud coefficient
  real(RP), private, save :: DELTCQ =  120.0E+0_RP !< time step for liquid-phase reaction
  real(RP), private, save :: DELTCG = 1200.0E+0_RP !< time step for gas-phase reaction
  real(RP), private, save :: FSO4   =    3.0E-2_RP !<

  real(RP), private, allocatable, save :: SULDIO_DATA  (:,:,:,:) !<    (IA,JA,14,12)
  real(RP), private, allocatable, save ::  BBSO2_DATA  (:,:,:,:) !<    (IA,JA,14,12)
  real(RP), private, allocatable, save :: VOLSO2_DATA  (:,:,  :) !<    (IA,JA,    1)
  real(RP), private, allocatable, save :: VOLALT_DATA  (:,:,  :) !<    (IA,JA,    1) 
  real(RP), private, allocatable, save :: ANSO4T_DATA  (:,:,  :) !<    (IA,JA,    1)
  real(RP), private, allocatable, save ::  HOORI_DATA(:,:,:,  :) !< (KA,IA,JA,   12) 
  ! real(RP), private, allocatable, save :: SVLSO2_DATA() !< 
  ! real(RP), private, allocatable, save :: SVLALT_DATA() !< 
  ! real(RP), private, allocatable, save :: SVLCLD_DATA() !<
  
  real(RP), private, allocatable, save :: QTRCINP(:,:,:,:) !< mixing ratio in prec(MP)

  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  character(len=7), private, parameter :: ctype(3) = (/ 'sulfate', 'SO2    ', 'DMS    ' /)
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Sulfate setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_setup &
       ( LON_REAL, LAT_REAL,        & ! [IN]
         GDZ, GDZM,                 & ! [IN]
         QA_AE_SU,                  & ! [IN]
         fileid                     ) ! [IN]
    implicit none
    real(RP),     intent(in) :: LON_REAL     (IA,JA) !< longitude [rad]
    real(RP),     intent(in) :: LAT_REAL     (IA,JA) !< latitude  [rad]
    real(RP),     intent(in) :: GDZ     (KA,  IA,JA) !< CZ
    real(RP),     intent(in) :: GDZM    (0:KA,IA,JA) !< FZ
    integer,      intent(in) :: QA_AE_SU
    integer,      intent(in) :: fileid
    character :: cnumber*2
    character(len=H_MID) :: varname
    real(RP) :: work2d(IA,JA), work3d(KA,IA,JA)
    integer  :: ierr, iq, ip
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSU_NMAESU / &
         FINCSU, &
         DELTCQ, &
         DELTCG, &
         FSO4
    !----------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosulf",*) 'Setup'

    !--- read namelist
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSU_NMAESU, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosulf",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSU_NMAESU. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_aerosulf",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSU_NMAESU, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROSU_NMAESU)
    
    !--- read emission data (netcdf file)
    allocate( SULDIO_DATA   (IA,JA,22,12) ) ! 14 -> 2000-2013
    allocate(  BBSO2_DATA   (IA,JA,22,12) ) ! 14 -> 2000-2013
    allocate( VOLSO2_DATA   (IA,JA,    1) ) !
    allocate( VOLALT_DATA   (IA,JA,    1) ) !
    allocate( ANSO4T_DATA   (IA,JA,    1) ) !
    allocate(  HOORI_DATA(KA,IA,JA,   12) ) !
    SULDIO_DATA  (:,:,:,:) = 0.0_RP
     BBSO2_DATA  (:,:,:,:) = 0.0_RP
    VOLSO2_DATA  (:,:,  :) = 0.0_RP
    VOLALT_DATA  (:,:,  :) = 0.0_RP
    ANSO4T_DATA  (:,:,  :) = 0.0_RP
     HOORI_DATA(:,:,:,  :) = 0.0_RP

    do iq = 1, 22
       write(varname,'(A,I0)') 'ANTSO2_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          SULDIO_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    do iq = 1, 22
       write(varname,'(A,I0)') 'BBSO2_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          BBSO2_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    call FILE_CARTESC_read( fileid, 'VOLSO2', 'XY', work2d(:,:) )
    VOLSO2_DATA(:,:,1) = work2d(:,:)
    call FILE_CARTESC_read( fileid, 'VOLALT', 'XY', work2d(:,:) )
    VOLALT_DATA(:,:,1) = work2d(:,:)
    call FILE_CARTESC_read( fileid, 'ANTSO4', 'XY', work2d(:,:) )
    ANSO4T_DATA(:,:,1) = work2d(:,:)
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, 'H2O2', 'ZXY',       &
                               work3d(:,:,:),               &
                               ip                           )
       HOORI_DATA(:,:,:,ip) = work3d(:,:,:)
    enddo

    allocate( QTRCINP(KA,IA,JA,2) )
    QTRCINP(:,:,:,:) = 0.0_RP

    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP
    call FILE_HISTORY_in( zero_2d(:,:),     'DMSFLX_SPRINTARS',         'DMS emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'SO2FLX_SPRINTARS',         'SO2 emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'SO4FLX_SPRINTARS',         'SO4 emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'TDFTSU_SPRINTARS',         'total deposition flux (sulfate)',        'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFTSU_SPRINTARS',         'wet deposition flux (sulfate)',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFTSU_CP_SPRINTARS',      'wet deposition flux (sulfate;CP)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFTSU_MP_RAIN_SPRINTARS', 'wet deposition flux (sulfate;MP(RAIN))', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFTSU_MP_SNOW_SPRINTARS', 'wet deposition flux (sulfate;MP(SNOW))', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'DRFTSU_SPRINTARS',         'dry deposition flux (sulfate)',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'GRFTSU_SPRINTARS',         'gravitational settling flux (sulfate)',  'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'TDFSO2_SPRINTARS',         'total deposition flux (SO2)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'CFSO2_SPRINTARS',          'chemical reaction flux (SO2)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'CFSO2W_SPRINTARS',         'chemical reaction flux (SO2;wet)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFSO2_SPRINTARS',         'wet deposition flux (SO2)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFSO2_CP_SPRINTARS',      'wet deposition flux (SO2;CP)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFSO2_MP_RAIN_SPRINTARS', 'wet deposition flux (SO2;MP(RAIN))',     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'WEFSO2_MP_SNOW_SPRINTARS', 'wet deposition flux (SO2;MP(SNOW))',     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'DRFSO2_SPRINTARS',         'dry deposotion flux (SO2)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'TDFDMS_SPRINTARS',         'total deposition flux (DMS)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'CFDMS_SPRINTARS',          'chemical reaction flux (DMS)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),     'DRFDMS_SPRINTARS',         'dry deposition flux (DMS)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:),   'CSO2K_SPRINTARS',          'chemical reaction flux (SO2)',           'kg/m2/s', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:),   'CSO2WK_SPRINTARS',         'chemical reaction flux (SO2;wet)',       'kg/m2/s', dim_type = 'ZXY' )
    do iq = 1, QA_AE_SU
       write(cnumber,'(i2.2)') iq
       call FILE_HISTORY_in( zero_3d(:,:,:), 'CONS'//cnumber//'_SPRINTARS', &
                         'mass concentration ('//trim(adjustl(ctype(iq)))//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_3d(:,:,:), 'SPPT'//cnumber//'_SPRINTARS', &
                         'volume concentration ('//trim(adjustl(ctype(iq)))//')', 'pptv',  dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_2d(:,:),   'COLS'//cnumber//'_SPRINTARS', &
                         'column loading ('//trim(adjustl(ctype(iq)))//')',       'kg/m2', dim_type = 'XY'  )
    enddo
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDST_SPRINTARS',  'mixing ratio (sulfate)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCST_SPRINTARS',  'number concentration (sulfate)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCSO2_SPRINTARS', 'number concentration (SO2)',     '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCDMS_SPRINTARS', 'number concentration (DMS)',     '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'PH_SPRINTARS',     'pH',                             '1',     dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'HOORI_SPRINTARS',  'H2O2',                           'ppm',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUSU_SPRINTARS',  'AOT (sulfate;550nm;allsky)',          '1',      dim_type = 'XY'  ) !  550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSU1_SPRINTARS', 'extinction (sulfate;532nm;allsky)',   '1/m',    dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSU1_SPRINTARS', 'backscatter (sulfate;532nm;allsky)',  '1/m/sr', dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSU2_SPRINTARS', 'extinction (sulfate;1064nm;allsky)',  '1/m',    dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSU2_SPRINTARS', 'backscatter (sulfate;1064nm;allsky)', '1/m/sr', dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSU3_SPRINTARS', 'extinction (sulfate;355nm;allsky)',   '1/m',    dim_type = 'ZXY' ) !  355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSU3_SPRINTARS', 'backscatter (sulfate;355nm;allsky)',  '1/m/sr', dim_type = 'ZXY' ) !  355nm
    do iq = 1, QA_AE_SU
       write(cnumber,'(i2.2)') iq
       call FILE_HISTORY_in( zero_3d(:,:,:), 'CNSC'//cnumber//'_SPRINTARS', &
                        'mass concentration ('//trim(adjustl(ctype(iq)))//';clearsky)', 'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_2d(:,:)  , 'CLSC'//cnumber//'_SPRINTARS', &
                        'column loading ('//trim(adjustl(ctype(iq)))//';clearsky)',     'kg/m2', dim_type = 'XY'  )
    enddo
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCSTC_SPRINTARS',  'number concentration (sulfate;clearsky)', '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUSUC_SPRINTARS',  'AOT (sulfate;550nm;clearsky)',            '1',      dim_type = 'XY'  ) !  550nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSUC1_SPRINTARS', 'extinction (sulfate;532nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSUC1_SPRINTARS', 'backscatter (sulfate;532nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSUC2_SPRINTARS', 'extinction (sulfate;1064nm;clearsky)',    '1/m',    dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSUC2_SPRINTARS', 'backscatter (sulfate;1064nm;clearsky)',   '1/m/sr', dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTSUC3_SPRINTARS', 'extinction (sulfate;355nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !  355nm
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKSUC3_SPRINTARS', 'backscatter (sulfate;355nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !  355nm

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_setup
  !-----------------------------------------------------------------------------
  !> SPRINTARS Sulfate finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_finalize
    implicit none

    deallocate( SULDIO_DATA )
    deallocate(  BBSO2_DATA ) 
    deallocate( VOLSO2_DATA )
    deallocate( VOLALT_DATA )
    deallocate( ANSO4T_DATA )
    deallocate(  HOORI_DATA )
    deallocate( QTRCINP )
    deallocate( zero_3d )
    deallocate( zero_2d )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_finalize
  
  !-----------------------------------------------------------------------------
  !> SPRINTARS RSO2AQ (liquid phase reaction)
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_RSO2AQ &
       ( SO2G,                                & ! [MODIFIED]
         SULA, SO2A, PH,                      & ! [OUT]
         WATER, OXIGX, SOLOXX, OXIGY, SOLOXY, & ! [IN]
         EQU1, EQU2, SOL1, GDT, dt_AE, DELTCQ ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: SO2G(KA,IA,JA) !< existed SO2 [molecules/cm3]
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: SULA(KA,IA,JA) !< formed sulfate [molecules/cm3]
    real(RP), intent(out)   :: SO2A(KA,IA,JA) !< existed SO2 [molecules/cm3]
    real(RP), intent(out)   :: PH  (KA,IA,JA) !< pH
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: WATER (KA,IA,JA) !< 
    real(RP), intent(in)    :: OXIGX (KA,IA,JA) !< existed oxidant [molecules/cm3]
    real(RP), intent(in)    :: SOLOXX(KA,IA,JA) !< Henry constant for SO2
    real(RP), intent(in)    :: OXIGY (KA,IA,JA) !< existed oxidant [molecules/cm3]
    real(RP), intent(in)    :: SOLOXY(KA,IA,JA) !< Henry constant for SO2
    real(RP), intent(in)    :: EQU1  (KA,IA,JA) !< reaction constant
    real(RP), intent(in)    :: EQU2  (KA,IA,JA) !< reaction constant
    real(RP), intent(in)    :: SOL1  (KA,IA,JA) !< reaction constant for SO2
    real(RP), intent(in)    :: GDT   (KA,IA,JA) !< temperature
    real(DP), intent(in)    :: dt_AE            !< time step
    real(RP), intent(in)    :: DELTCQ           !< time step for liquid-phase reaction

    ! internal work
    !-------------------------
    real(RP) :: SO4XT(KA,IA,JA) !< total sulfate concentration formed by SO4X [molecule/cm3]
    real(RP) :: SO4YT(KA,IA,JA) !< total sulfate concentration formed by SO4Y [molecule/cm3]
    real(RP) :: OXIAX(KA,IA,JA) !< existed oxidant [molecule/cm3]
    real(RP) :: OXIAY(KA,IA,JA) !< existed oxidant [molecule/cm3]
    real(RP) :: SO4X            !< formed sulfate  [molecule/cm3]
    real(RP) :: SO4Y            !< formed sulfate  [molecule/cm3]
    real(RP) :: SOLC            !< Henry constant for SO2
    real(RP) :: EXX             !< 
    real(RP) :: EXY             !< 
    real(RP) :: KFORMX          !< reaction rate, which depends on pH
    real(RP) :: KFORMY          !< reaction rate, which depends on pH
    real(RP) :: SOLADD          !< 
    real(RP) :: HION            !< H+ ion concentration [mol/L]
    real(RP) :: DELTC           !< time step for liquid-phase reaction
    
    ! parameter
    !-------------------------
    real(RP), parameter :: F1  = 0.1E+0_RP            !< 
    real(RP), parameter :: EPS = 1.0E-6_RP            !< 
    real(RP), parameter :: HI0 = 10.0_RP**( -5.6_RP ) !< H+ ion concentration [mol/L] (background)

    ! others
    !-------------------------
    integer :: k, i, j
    integer :: NS, NSMAX

    !----------------------------------------------------
    SULA (:,:,:) = 0.0_RP
    SO2A (:,:,:) = 0.0_RP
    SO4XT(:,:,:) = 0.0_RP
    SO4YT(:,:,:) = 0.0_RP
    OXIAX(:,:,:) = 0.0_RP
    OXIAY(:,:,:) = 0.0_RP
    
    NSMAX = int( real(dt_AE,kind=RP) / DELTCQ )
    NSMAX = max( NSMAX, 1 )
    DELTC = real(dt_AE,kind=RP) / real(NSMAX,kind=RP)

    do NS = 1, NSMAX
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( WATER(k,i,j) > EPS ) then

                !pH calculation
                HION = HI0 + F1 * ( 2.0_RP * ( SO4XT(k,i,j) + SO4YT(k,i,j) ) + SO2A(k,i,j) ) / AVOG / WATER(k,i,j) * 1.0E+6_RP
                HION = min( max( HION, 10.0_RP**(-5.6_RP) ), 10.0_RP**(-4.0_RP) )

                !reaction rate, which depends on pH.
                KFORMX = 5.6E+6_RP * HION / ( 1.0E-1_RP + HION )
                KFORMY = ( 2.0E+4_RP + 3.2E+5_RP * EQU1(k,i,j) / HION + 1.0E+9_RP * EQU1(k,i,j) * EQU2(k,i,j) / HION**2 ) &
                       / ( 1.0_RP + ( EQU1(k,i,j)/HION ) * ( 1.0_RP + EQU2(k,i,j)/HION ) )

                !partitioning coefficient
                SOLADD = SOL1(k,i,j) + SOL1(k,i,j) * EQU1(k,i,j) / HION + SOL1(k,i,j) * EQU1(k,i,j) * EQU2(k,i,j) / HION**2
                SOLC   = 1.0_RP + SOLADD * 1.0E+3_RP * 0.0821_RP * GDT(k,i,j) * WATER(k,i,j) * 1.0E-6_RP

                !re-partitioning gas-aqueous phase after aqueous reaction
                SO2G (k,i,j) = ( SO2A(k,i,j) + SO2G(k,i,j) ) / SOLC
                SO2A (k,i,j) = SO2G(k,i,j) * ( SOLC - 1.0_RP )
                OXIAX(k,i,j) = ( OXIAX(k,i,j) + OXIGX(k,i,j) ) / SOLOXX(k,i,j) * ( SOLOXX(k,i,j) - 1.0_RP )

                if ( SO2A(k,i,j) - 100.0_RP * OXIAX(k,i,j) > 0.0_RP ) then
                   SO4X = OXIAX(k,i,j)
                   
                elseif (       SO2A(k,i,j) - 100.0_RP * OXIAX(k,i,j) <= 0.0_RP &
                         .and. SO2A(k,i,j) -            OXIAX(k,i,j) >  0.0_RP ) then
                   EXX = ( SO2A(k,i,j) - OXIAX(k,i,j) ) * KFORMX / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC
                   if ( EXX >= 100.0_RP ) then
                      SO4X = OXIAX(k,i,j)
                   else
                      SO4X = SO2A(k,i,j) * OXIAX(k,i,j) * ( 1.0_RP - exp(EXX) ) / ( OXIAX(k,i,j) - SO2A(k,i,j) * exp(EXX) )
                   endif
                   
                elseif ( SO2A(k,i,j) - OXIAX(k,i,j) == 0.0_RP ) then
                   SO4X = SO2A(k,i,j)**2 * KFORMX / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC &
                         / ( SO2A(k,i,j) * KFORMX / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC + 1.0_RP )

                elseif (       100.0_RP * SO2A(k,i,j) - OXIAX(k,i,j) >= 0.0_RP &
                         .and. SO2A(k,i,j) - OXIAX(k,i,j) < 0.0_RP             ) then
                   EXX = ( SO2A(k,i,j) - OXIAX(k,i,j) ) * KFORMX / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC
                   if ( EXX <= -100.0_RP ) then
                      SO4X = SO2A(k,i,j)
                   else
                      SO4X = SO2A(k,i,j) * OXIAX(k,i,j) * ( 1.0_RP - exp(EXX) ) / ( OXIAX(k,i,j) - SO2A(k,i,j) * exp(EXX) )
                   endif

                else
                   SO4X = SO2A(k,i,j)
                endif
                
                SO4XT(k,i,j) = SO4XT(k,i,j) + SO4X
                SO2A (k,i,j) = SO2A (k,i,j) - SO4X
                OXIAX(k,i,j) = OXIAX(k,i,j) - SO4X
                if ( SO2A (k,i,j) < 0.0_RP ) SO2A (k,i,j) = 0.0_RP
                if ( OXIAX(k,i,j) < 0.0_RP ) OXIAX(k,i,j) = 0.0_RP 

                !re-partitioning gas-aqueous phase after aqueous reaction
                SO2G (k,i,j) = ( SO2A(k,i,j) + SO2G(k,i,j) ) / SOLC
                SO2A (k,i,j) = SO2G(k,i,j) * ( SOLC - 1.0_RP )
                OXIAY(k,i,j) = ( OXIAY(k,i,j) + OXIGY(k,i,j) ) / SOLOXY(k,i,j) * ( SOLOXY(k,i,j) - 1.0_RP )
                
                if ( SO2A(k,i,j) - 100.0_RP * OXIAY(k,i,j) > 0.0_RP ) then
                   SO4Y = OXIAY(k,i,j)
                   
                elseif (       SO2A(k,i,j) - 100.0_RP * OXIAY(k,i,j) <= 0.0_RP &
                         .and. SO2A(k,i,j) -            OXIAY(k,i,j) >  0.0_RP ) then
                   EXY = ( SO2A(k,i,j) - OXIAY(k,i,j) ) * KFORMY / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC
                   if ( EXY >= 100.0_RP ) then
                      SO4Y = OXIAY(k,i,j)
                   else
                      SO4Y = SO2A(k,i,j) * OXIAY(k,i,j) * ( 1.0_RP - exp(EXY) ) / ( OXIAY(k,i,j) - SO2A(k,i,j) * exp(EXY) )
                   endif
                   
                elseif ( SO2A(k,i,j) - OXIAY(k,i,j) == 0.0_RP ) then
                   SO4Y = SO2A(k,i,j)**2 * KFORMY / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC &
                         / ( SO2A(k,i,j) * KFORMY / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC + 1.0_RP )
               
                elseif (       100.0_RP * SO2A(k,i,j) - OXIAY(k,i,j) >= 0.0_RP &
                         .and. SO2A(k,i,j) - OXIAY(k,i,j) < 0.0_RP ) then
                   EXY = ( SO2A(k,i,j) - OXIAY(k,i,j) ) * KFORMY / AVOG / WATER(k,i,j) * 1.0E+6_RP * DELTC
                   if ( EXY <= -100.0_RP ) then
                      SO4Y = SO2A(k,i,j)
                   else
                      SO4Y = SO2A(k,i,j) * OXIAY(k,i,j) * ( 1.0_RP - exp(EXY) ) / ( OXIAY(k,i,j) - SO2A(k,i,j) * exp(EXY) )
                   endif

                else
                   SO4Y = SO2A(k,i,j)
                endif
                
                SO4YT(k,i,j) = SO4YT(k,i,j) + SO4Y
                SO2A (k,i,j) = SO2A (k,i,j) - SO4Y
                OXIAY(k,i,j) = OXIAY(k,i,j) - SO4Y
                if ( SO2A (k,i,j) < 0.0_RP ) SO2A (k,i,j) = 0.0_RP
                if ( OXIAY(k,i,j) < 0.0_RP ) OXIAY(k,i,j) = 0.0_RP

             endif
                
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          if ( WATER(k,i,j) > EPS ) then
             HION = HI0 + F1 * ( 2.0_RP * ( SO4XT(k,i,j) + SO4YT(k,i,j) ) + SO2A(k,i,j) ) / AVOG / WATER(k,i,j) * 1.0E+6_RP
             HION = min( max( HION, 10.0_RP**(-5.6_RP) ), 10.0_RP**(-4.0_RP) )
             PH  (k,i,j) = - log10( HION )
             SULA(k,i,j) = SO4XT(k,i,j) + SO4YT(k,i,j)
          else
             PH  (k,i,j) = CONST_UNDEF !VMISS
          endif
       enddo
    enddo
    enddo
             
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_RSO2AQ
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS QSSA (gas phase reaction)
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_QSSA &
       ( SO2, DMS,                  & ! [MODIFIED]
         C3,                        & ! [OUT]
         REACRA, OHM, dt_AE, DELTCG ) ! [IN]
    implicit none
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: SO2   (KA,IA,JA)   !< SO2 concentration [molecules/cm3]
    real(RP), intent(inout) :: DMS   (KA,IA,JA)   !< DMS concentration [molecules/cm3]
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: C3    (KA,IA,JA)   !< SO4 concentration [molecules/cm3]
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: REACRA(KA,IA,JA,3) !< reaction rate    [cm3/molecule/s]
    real(RP), intent(in)    :: OHM   (KA,IA,JA)   !< OH concentration [molecules/cm3]
    real(DP), intent(in)    :: dt_AE              !< time step
    real(RP), intent(in)    :: DELTCG             !< time step for gas-phase reaction
    
    ! internal work
    !-------------------------
    real(RP) :: C1(KA,IA,JA)                      !< DMS concentration [molecules/cm3]
    real(RP) :: C2(KA,IA,JA)                      !< SO2 concentration [molecules/cm3]
    real(RP) :: B1(KA,IA,JA)                      !< for DMS reaction : ( k1 + k2 ) * [OH] * dt
    real(RP) :: B2(KA,IA,JA)                      !< for SO2 reaction : ( k1 + 0.6k2 + 0.24k2 ) * [OH] * dt
    real(RP) :: B3(KA,IA,JA)                      !< for SO4 reaction : 0.16 * k2 * [OH] * dt
    real(RP) :: B4(KA,IA,JA)                      !< for SO4 reaction : k3 * [OH] * dt
    real(RP) :: C1B                               !< DMS concentration
    real(RP) :: C2B                               !< SO2 concentration
    real(RP) :: DELTC                             !< time step for gas-phase reaction
    real(RP) :: WRK                               !< [OH] * dt
    
    ! others
    !-------------------------
    integer :: k, i, j
    integer :: it
    integer :: TMAX
    
    !----------------------------------------------------

    TMAX  = int( real(dt_AE,kind=RP) / DELTCG )
    TMAX  = max( TMAX, 1 )
    DELTC = real(dt_AE,kind=RP) / real(TMAX,kind=RP)
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          WRK = OHM(k,i,j) * DELTC
          B1(k,i,j) = ( REACRA(k,i,j,1) + REACRA(k,i,j,2) ) * WRK
          B2(k,i,j) = ( REACRA(k,i,j,1) + ( 0.60_RP + 0.24_RP ) * REACRA(k,i,j,2) ) * WRK
          B3(k,i,j) = 0.16_RP * REACRA(k,i,j,2) * WRK
          B4(k,i,j) = REACRA(k,i,j,3) * WRK
          C1(k,i,j) = max( DMS(k,i,j), 0.0_RP )
          C2(k,i,j) = max( SO2(k,i,j), 0.0_RP )
          C3(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo

    do it = 1, TMAX
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             C1B = C1(k,i,j)
             C2B = C2(k,i,j)
             C1(k,i,j) = C1B * ( 1.0_RP - B1(k,i,j) )
             if ( C1(k,i,j) < 0.0_RP ) then
                C1(k,i,j) = 0.0_RP
                C2B       = C2B       + B2(k,i,j) / B1(k,i,j) * C1B
                C3(k,i,j) = C3(k,i,j) + B3(k,i,j) / B1(k,i,j) * C1B
             else
                C2B       = C2B       + B2(k,i,j) * C1B
                C3(k,i,j) = C3(k,i,j) + B3(k,i,j) * C1B
             endif
             C2(k,i,j) = C2B * ( 1.0_RP - B4(k,i,j) )
             if ( C2(k,i,j) < 0.0_RP ) then
                C2(k,i,j) = 0.0_RP
                C3(k,i,j) = C3(k,i,j) + C2B
             else
                C3(k,i,j) = C3(k,i,j) + B4(k,i,j) * C2B
             endif
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          DMS(k,i,j) = C1(k,i,j)
          SO2(k,i,j) = C2(k,i,j)
       enddo
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU_QSSA
  !-----------------------------------------------------------------------------

  
  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Sulfate
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU &
       ( QTRC,           & ! [MODIFIED]
         OUTQLD,         & ! [OUT]
         NUMSUL,         & ! [OUT]
         PM25SU,         & ! [OUT]
         WATRSU,         & ! [OUT]
         TAUSU,          & ! [OUT]
         CEXTSU,         & ! [OUT]
         CBAKSU,         & ! [OUT]
         TAUSUD,         & ! [OUT]
         CEXTSD,         & ! [OUT]
         CDVE,           & ! [IN]
         CCOVMR,         & ! [IN]
         CCMAX,          & ! [IN]
         GDT,            & ! [IN]
         GOICE,          & ! [IN]
         RFLXSD,         & ! [IN]
         RCOSZ,          & ! [IN]
         GRLAI,          & ! [IN]
         KBB,            & ! [IN]
         DPKBB,          & ! [IN]
         THETA,          & ! [IN]
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
         QC,             & ! [IN]
         TCLDF,          & ! [IN]
         TCLDF2,         & ! [IN]
         KUP,            & ! [IN]
         KUPMAX,         & ! [IN]
         DPKUP,          & ! [IN]
         FLIQ,           & ! [IN]
         OHRAD,          & ! [IN]
         O3ORI,          & ! [IN]
         GDZS,           & ! [IN]
         GDPM,           & ! [IN]
         DELP,           & ! [IN]
         GDZM,           & ! [IN]
         DELZ,           & ! [IN]
         ONEWQ,          & ! [IN]
         RNWR,           & ! [IN]
         RNWR_CP,        & ! [IN]
         RNWR_MP_RAIN,   & ! [IN]
         RNWR_MP_SNOW,   & ! [IN]
         FACT_LAND,      & ! [IN]
         FACT_OCEAN,     & ! [IN]
         AIRFEL,         & ! [IN]
         ACTISU,         & ! [IN]
         SW_CCN,         & ! [IN]
         TIME_NOWDATE,   & ! [IN]
         dt_AE,          & ! [IN]
         QA_AE_SU        ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE                      !< time step
    integer,  intent(in)    :: QA_AE_SU                   !< total number of sulfate tracers

    real(RP), intent(in)    :: CDVE               (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
    real(RP), intent(in)    :: CCOVMR             (IA,JA) !< cloud cover (max-ran)
    real(RP), intent(in)    :: CCMAX                      !< maximum cloud fraction for all/clear sky
    real(RP), intent(in)    :: GDT           (KA,  IA,JA) !< temperature
    real(RP), intent(in)    :: GOICE              (IA,JA) !< sea ice
    real(RP), intent(in)    :: RFLXSD             (IA,JA) !< downward short wave radiation
    real(RP), intent(in)    :: RCOSZ              (IA,JA) !< cos (solar angle)
    real(RP), intent(in)    :: GRLAI              (IA,JA) !< leaf area index [m2/m2]
    integer,  intent(in)    :: KBB                (IA,JA) !<
    real(RP), intent(in)    :: DPKBB              (IA,JA) !<
    real(RP), intent(in)    :: THETA         (KA,  IA,JA) !< potential temperature
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
    real(RP), intent(in)    :: QC            (KA,  IA,JA) !< cloud water
    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA) !< cloud fraction (CP+MP)
    real(RP), intent(in)    :: TCLDF2        (KA,  IA,JA) !< cloud fraction (CP+MP) for CCOVMR
    integer,  intent(in)    :: KUP                (IA,JA) !< uplift level
    integer,  intent(in)    :: KUPMAX                     !< maximum uplift level
    real(RP), intent(in)    :: DPKUP              (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
    real(RP), intent(in)    :: FLIQ          (KA,  IA,JA) !< water/ice partition
    real(RP), intent(in)    :: OHRAD         (KA,  IA,JA) !< OH concentration
    real(RP), intent(in)    :: O3ORI         (KA,  IA,JA) !< O3 concentration
    real(RP), intent(in)    :: GDZS               (IA,JA) !< absolute ground height [m]
    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA) !< pressure at half levels
    real(RP), intent(in)    :: DELP          (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA) !< altitude at half levels
    real(RP), intent(in)    :: DELZ          (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
    logical,  intent(in)    :: ONEWQ         (KA-1,IA,JA) !< use for instability
    real(RP), intent(in)    :: RNWR          (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: FACT_LAND          (IA,JA) !< land  factor
    real(RP), intent(in)    :: FACT_OCEAN         (IA,JA) !< ocean factor
    real(RP), intent(in)    :: AIRFEL        (KA,  IA,JA) !< SO2 emission from aircrafts
    real(RP), intent(in)    :: ACTISU        (KA,  IA,JA) !< ratio of activated aerosols [0-1]
    logical,  intent(in)    :: SW_CCN                     !< switch for CCN scheme (true->use ACTI)
    integer,  intent(in)    :: TIME_NOWDATE(6)            !< current time
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_SU) !< ratio of sulfate tracer mass to total mass

    ! output
    !-------------------------
    real(RP), intent(out)   :: OUTQLD(KA,IA,JA)        !< QTRC for radiation
    real(RP), intent(out)   :: NUMSUL(KA,IA,JA)        !< number concentration
    real(RP), intent(out)   :: PM25SU(KA,IA,JA)        !< submicron mass concentration
    real(RP), intent(out)   :: WATRSU(KA,IA,JA)        !< aerosol water
    real(RP), intent(out)   :: TAUSU    (IA,JA,IWA,2)  !< AOT
    real(RP), intent(out)   :: CEXTSU(KA,IA,JA,IWA,2)  !< extinction  coefficient
    real(RP), intent(out)   :: CBAKSU(KA,IA,JA,IWA,2)  !< backscatter coefficient
    real(RP), intent(out)   :: TAUSUD   (IA,JA)        !< AOT (dry)
    real(RP), intent(out)   :: CEXTSD(KA,IA,JA)        !< extinction  coefficient (sulfate;dry)
    
    ! internal work
    !-------------------------
    !set parameter
    real(RP) :: DRYMSU                      !< dry particle mass                 [kg]
    real(RP) :: SULDIO   (IA,JA)            !< emission flux (fossil fuel)       [kg/m2/s]
    real(RP) :: BBSO2    (IA,JA)            !< emission flux (biomass burning)   [kg/m2/s]
    real(RP) :: VOLSO2   (IA,JA)            !< emission flux (volcanic)          [kg/m2/s]
    real(RP) :: VOLALT   (IA,JA)            !< volcanic altitude                 []
    real(RP) :: ANSO4T   (IA,JA)            !< emission flux (SO4)               [kg/m2/s]
    real(RP) :: SVLSO2   (IA,JA)            !< volcanic emission flux (sporadic) [kg/m2/s]
    real(RP) :: SVLALT   (IA,JA)            !< volcanic altitude (sporadic)      []
    real(RP) :: SVLCLD   (IA,JA)            !< volcanic cloud height (sporadic)  [kg/m2/s]
    real(RP) :: HOORI (KA,IA,JA)            !< H2O2 concentration                []
    real(RP) :: RADS  (KA,IA,JA)            !< sulfate particle radius           [m]
    real(RP) :: CEXTP (KA,IA,JA,IWA)        !< extinction  cross section (550/ 440/870nm)
    real(RP) :: CEXTA (KA,IA,JA,IWA)        !< extinction  cross section (532/1064/355nm)
    real(RP) :: CBAKA (KA,IA,JA,IWA)        !< backscatter cross section (532/1064/355nm)
    real(RP) :: BRH1, BRH2, BRH3            !<
    real(RP) :: VTER  (KA,IA,JA)            !< terminal velocity of sulfate aerosol [m/s]
    
    !emission
    integer  :: VOLK  (IA,JA)               !< volcano layer
    integer  :: VOLKS (IA,JA)               !< volcano layer
    integer  :: KUPS  (IA,JA)               !< uplift level
    integer  :: KUPMXS                      !< KUPS maximum
    real(RP) :: DMSFLX(IA,JA)               !< DMS emission flux [kg/m2/s]
    real(RP) :: ANTSO4(IA,JA)               !< emission flux (SO4) [kg/m2/s]
    real(RP) :: DPVOL                       !< GDPM(i,j,VOLK (i,j)-1) - GDPM(i,j,VOLK(i,j))
    real(RP) :: DPSVL                       !< GDPM(i,j,VOLKS(i,j)-1) - GDPM(i,j,KUPS(i,j))
    real(RP) :: QEMIT (IA,JA,6)             !< emission mixing ratio [kg/kg]
    real(RP) :: SO2FLX(IA,JA)               !< SO2 emission flux[kg/m2/s]
    
    !chemical reaction
    real(RP) :: AIRNUM(KA,IA,JA)            !< air number density   [molecules/cm3]
    real(RP) :: WATER (KA,IA,JA)            !< cloud water          [kg/m3]
    real(RP) :: NUMO3 (KA,IA,JA)            !< O3 number density    [molecules/cm3]
    real(RP) :: NUMHO (KA,IA,JA)            !< H2O2 number density  [molecules/cm3]
    real(RP) :: INREAC(4)                   !< use for reaction
    real(RP) :: ALPHA1                      !< use for reaction
    real(RP) :: ALPHA2                      !< use for reaction
    real(RP) :: ALPHA3                      !< reaction rate        [cm3/molecules/s]
    real(RP) :: REACRA(KA,IA,JA,3)          !< reaction rate        [cm3/molecules/s]
    real(RP) :: INEXP1                      !< use for reaction
    real(RP) :: INEXP2                      !< use for reaction
    real(RP) :: INEPO3                      !< use for reaction
    real(RP) :: INEPHO                      !< use for reaction
    real(RP) :: SOL1  (KA,IA,JA)            !< solubility constant  [mol/L/atm]
    real(RP) :: SOLO3 (KA,IA,JA)            !< solubility constant  [mol/L/atm]
    real(RP) :: SOLHO (KA,IA,JA)            !< solubility constant  [mol/L/atm]
    real(RP) :: SOLOX1(KA,IA,JA)            !< use for reaction
    real(RP) :: SOLOX2(KA,IA,JA)            !< use for reaction
    real(RP) :: INEPQ1                      !< use for reaction
    real(RP) :: INEPQ2                      !< use for reaction
    real(RP) :: EQU1  (KA,IA,JA)            !< equilibrium constant [mol/L]
    real(RP) :: EQU2  (KA,IA,JA)            !< equilibrium constant [mol/L]
    real(RP) :: NUMORI(KA,IA,JA,QA_AE_SU)   !< S number density     [molecules/cm3]
    real(RP) :: SO2G1 (KA,IA,JA)            !< use for reaction     [molecules/cm3]
    real(RP) :: SULA  (KA,IA,JA)            !< use for reaction
    real(RP) :: SO2A  (KA,IA,JA)            !< use for reaction
    real(RP) :: PH    (KA,IA,JA)            !< pH

    real(RP) :: CHEDPT(KA,IA,JA,3)          !< 
    real(RP) :: NUMNW2(KA,IA,JA,QA_AE_SU,3) !< use for reaction
    
    real(RP) :: DMS1  (KA,IA,JA)            !< use for reaction [molecules/cm3]
    real(RP) :: DMS2  (KA,IA,JA)            !< use for reaction [molecules/cm3]
    real(RP) :: SO2G2 (KA,IA,JA)            !< use for reaction [molecules/cm3]
    real(RP) :: SULG1 (KA,IA,JA)            !< use for reaction
    real(RP) :: SULG2 (KA,IA,JA)            !< use for reaction
    real(RP) :: NUMNW (KA,IA,JA,QA_AE_SU)   !< 
    real(RP) :: CHEDEP(KA,IA,JA,QA_AE_SU)   !<
    
    real(RP) :: ADDNUM (KA,IA,JA)           !< 
    real(RP) :: CCNNUM (KA,IA,JA)           !< 
    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_SU)  !< 
    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_SU)  !< 
    real(RP) :: SO2WTM    (IA,JA)           !< 
    real(RP) :: SO2C
    real(RP) :: SO2CW
    real(RP) :: SO2CK  (KA,IA,JA)           !<
    real(RP) :: SO2CWK (KA,IA,JA)           !<
    
    real(RP) :: CHEDPM    (IA,JA,QA_AE_SU)  !< 
    real(RP) :: SO2CHE    (IA,JA)           !< SO2 chemical deposition          [kg/m2/s]
    real(RP) :: DMSCHE    (IA,JA)           !< DMS chemical deposition          [kg/m2/s]
    real(RP) :: SO2CHW    (IA,JA)           !< SO2 chemical deposition (aq)     [kg/m2/s]
    
    !subcloud scavenging (wash out)
    ! real(RP) :: RAINN                       !< rain number concentrarion        [1/m3]
    ! real(RP) :: RAINN_CP                    
    ! real(RP) :: RAINN_MP_RAIN               
    ! real(RP) :: RAINN_MP_SNOW           
    real(RP) :: SUN                         !< number density of sulfate tracer [1/m3]
    real(RP) :: SURN          (KA,IA,JA,2)  !< number density of sulfate tracer entering into raindrops [1/s/m3] <=> collision frequency between sulfate and rain
    real(RP) :: SURN_CP       (KA,IA,JA,2)
    real(RP) :: SURN_MP_RAIN  (KA,IA,JA,2)
    real(RP) :: SURN_MP_SNOW  (KA,IA,JA,2)
    real(RP) :: QWTDEP        (KA,IA,JA,2)  !< mixing ratio by wet deposition   [kg/kg] 
    real(RP) :: QWTDEP_CP     (KA,IA,JA,2)
    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,2)
    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,2)
    real(RP) :: QTRCTMP       (KA,IA,JA)    !< sulfate tracers that do not cllide with rain particles [kg/kg]
    
    !incloud scavenging (rain out)
    real(RP) :: CLTORN        (KA,IA,JA,2)  !< mixing ratio moving cloud to rain particles [kg/kg]
    real(RP) :: CLTORN_CP     (KA,IA,JA,2)
    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,2)
    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,2)
    
    !re-emission from rain
    real(RP) :: MSO2                     !< 
    real(RP) :: WETF           (IA,JA,2) !< wet deposition flux
    real(RP) :: WETF_CP        (IA,JA,2) !< wet deposition flux (CP)
    real(RP) :: WETF_MP_RAIN   (IA,JA,2) !< wet deposition flux (MP:rain)
    real(RP) :: WETF_MP_SNOW   (IA,JA,2) !< wet deposition flux (MP:snow)
    real(RP) :: QTRCINR     (KA,IA,JA,2) !< mixing ratio in rain(MP)
    real(RP) :: QTRCINS     (KA,IA,JA,2) !< mixing ratio in snow(MP)

    !dry deposition
    real(RP) :: CDVES   (IA,JA)           !< 
    real(RP) :: QMTX (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
    real(RP) :: DRYF    (IA,JA,QA_AE_SU)  !< dry deposition flux
    
    !gravitational settling
    real(RP) :: QLOAD (KA,IA,JA,QA_AE_SU)  !< mixing ratio for radiation [kg/kg]
    real(RP) :: SULGRV   (IA,JA)           !< gravitational settling flux (sulfate)
    
    !aerosol
    real(RP) :: SCON   (KA,IA,JA,QA_AE_SU)  !< mass concentration  [kg/m3]
    real(RP) :: COLMS     (IA,JA,QA_AE_SU)  !< column mass loading [kg/m2]
    real(RP) :: NUMSO2 (KA,IA,JA)           !< SO2 number concentration
    real(RP) :: NUMDMS (KA,IA,JA)           !< DMS number concentration
    real(RP) :: SPPT   (KA,IA,JA,QA_AE_SU)  !< volume concentration
    
    !flux
    real(RP) :: SULTDP    (IA,JA)           !< sulfate total deposition
    real(RP) :: SO2TDP    (IA,JA)           !< SO2 total deposition
    real(RP) :: DMSTDP    (IA,JA)           !< DMS total deposition
    
    !optical parameters
    real(RP) :: TUSE                        !<
    real(RP) :: NUMSUC (KA,IA,JA)           !< number concentration (clearsky)
    real(RP) :: SCONC  (KA,IA,JA,QA_AE_SU)  !< mass concentration (clearsky)
    real(RP) :: COLMSC    (IA,JA,QA_AE_SU)  !< column mass loading (clearsky)
    
    ! parameter
    !-------------------------
    real(RP),  parameter :: &
         MOLM(NSU) &            !< molecule mass
         = (/ 132.0E+0_RP, 64.0E+0_RP, 62.0E+0_RP /)
    real(RP),  parameter :: &
         DDV(4) &               !< deposition velocity
         = (/ 2.0E-3_RP, 8.0E-3_RP, 6.0E-3_RP, 1.0E-3_RP /)
    real(RP),  parameter :: &
         CEXTSP(NRH,IWA)   &    !< extinction  cross section (550/ 440/870nm)
         !             RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
         = reshape( (/ 1.801E+3_RP, 3.435E+3_RP, 5.092E+3_RP, 6.802E+3_RP, 1.255E+4_RP, 3.043E+4_RP, 6.185E+4_RP, 1.026E+5_RP, & !  550nm
                       3.109E+3_RP, 5.657E+3_RP, 8.145E+3_RP, 1.064E+4_RP, 1.866E+4_RP, 4.155E+4_RP, 7.764E+4_RP, 1.201E+5_RP, & !  440nm
                       4.843E+2_RP, 1.018E+3_RP, 1.598E+3_RP, 2.225E+3_RP, 4.500E+3_RP, 1.268E+4_RP, 2.978E+4_RP, 5.598E+4_RP /), & !  870nm
                    (/ NRH, IWA /) )
    real(RP),  parameter :: &
         CEXTSA(NRH,IWA)   &    !< extinction  cross section (532/1064/355nm)
         !             RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
         = reshape( (/ 1.963E+3_RP, 3.717E+3_RP, 5.486E+3_RP, 7.305E+3_RP, 1.338E+4_RP, 3.205E+4_RP, 6.437E+4_RP, 1.057E+5_RP, & !  532nm
                       2.480E+2_RP, 5.476E+2_RP, 8.827E+2_RP, 1.253E+3_RP, 2.636E+3_RP, 7.882E+3_RP, 1.957E+4_RP, 3.855E+4_RP, & ! 1064nm
                       5.159E+3_RP, 8.794E+3_RP, 1.218E+4_RP, 1.547E+4_RP, 2.553E+4_RP, 5.155E+4_RP, 8.809E+4_RP, 1.271E+5_RP /), & !  355nm
                    (/ NRH, IWA /) )
    real(RP),  parameter :: &
         CBAKSA(NRH,IWA)   &    !< backscatter cross section (532/1064/355nm)
         !             RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
         = reshape( (/ 4.204E+1_RP, 6.108E+1_RP, 7.876E+1_RP, 9.631E+1_RP, 1.530E+2_RP, 3.191E+2_RP, 5.928E+2_RP, 9.387E+2_RP, & !  532nm
                       1.428E+1_RP, 2.524E+1_RP, 3.493E+1_RP, 4.398E+1_RP, 7.046E+1_RP, 1.398E+2_RP, 2.594E+2_RP, 4.343E+2_RP, & ! 1064nm
                       7.736E+1_RP, 1.132E+2_RP, 1.459E+2_RP, 1.770E+2_RP, 2.701E+2_RP, 5.085E+2_RP, 8.886E+2_RP, 1.444E+3_RP /), & !  355nm
                    (/ NRH, IWA /) )
    real(RP),  parameter :: &
         RADSU(NRH)        &    !< sulfate particle radius
         !    RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
         = (/ 0.695E-7_RP, 0.850E-7_RP, 0.950E-7_RP, 0.103E-6_RP, 0.122E-6_RP, 0.157E-6_RP, 0.195E-6_RP, 0.231E-6_RP /)
    real(RP),  parameter :: DENSSU =    1.769E+3_RP    !< sulfate particle density [kg/m3]
    real(RP),  parameter :: AA     =    1.000E-4_RP    !< collision efficiency
    real(RP),  parameter :: MTONSU =    4.210E+16_RP   !< mcon -> ncon
    real(RP),  parameter :: DMSTMX =  308.15E+0_RP     !< max. tmp. of DMS flux
    real(RP),  parameter :: DELTCQ =  120.0E+0_RP      !< time step for liquid-phase reaction [s]
    real(RP),  parameter :: DELTCG = 1200.0E+0_RP      !< time step for gas-phase reaction [s]
    real(RP),  parameter :: EPS    =    1.0E-20_RP     !< 
    real(RP),  parameter :: AFFSO2 =    1.0E+0_RP      !< SO2 (input AIRSO2 in CMIP6)
    real(RP),  parameter :: SIGMA  =    1.526E+0_RP    !< GSD of log-normal distribution
!    character, parameter :: ctype(3)*7 = (/ 'sulfate', 'SO2', 'DMS' /)
!    character(len=7), parameter :: ctype(3) = (/ 'sulfate', 'SO2    ', 'DMS    ' /)
    
    ! others 
    !-------------------------
    integer   :: k, i, j, iq, iw
    integer   :: IRH
    character :: cnumber*2
    integer   :: nyear

    integer   :: nyear_m, nyear_p
    
    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - sulfate'

    !setup parameter 
    !================================================
    !DRYMSU
    !-------------------------------
    DRYMSU = 4.0_RP / 3.0_RP * PI * RADSU(1)**3 * DENSSU * exp( 4.5_RP * log(SIGMA)**2 )
          
    !set emission data
    !-------------------------------
    nyear = TIME_NOWDATE(1)-1998

    !if( nyear == 14 ) then
    !  nyear_p = nyear
    if( nyear >= 22 ) then
      nyear_p = 22 !tentative
      nyear = 22 
    else
      nyear_p = nyear + 1
    endif
    !if( nyear == 1 ) then
    !  nyear_m = nyear
    if( nyear <= 1 ) then
      nyear_m = 1 ! tentative
      nyear = 1 ! tentative
    else
      nyear_m = nyear - 1
    endif

    
    !ANTSO2
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( SULDIO(:,:),                                                                       & ! [OUT]
!           SULDIO_DATA(:,:,nyear-1,12), SULDIO_DATA(:,:,nyear,:), SULDIO_DATA(:,:,nyear+1,1), & ! [IN]
           SULDIO_DATA(:,:,nyear_m,12), SULDIO_DATA(:,:,nyear,:), SULDIO_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                                  ) ! [IN]

    !BBSO2
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( BBSO2(:,:),                                                                     & ! [OUT]
!           BBSO2_DATA(:,:,nyear-1,12), BBSO2_DATA(:,:,nyear,:), BBSO2_DATA(:,:,nyear+1,1), & ! [IN]
           BBSO2_DATA(:,:,nyear_m,12), BBSO2_DATA(:,:,nyear,:), BBSO2_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                               ) ! [IN]

    !VOLSO2, VOLALT
    VOLSO2(:,:) = VOLSO2_DATA(:,:,1)
    VOLALT(:,:) = VOLALT_DATA(:,:,1)

    !ANTSO4
    ANSO4T(:,:) = ANSO4T_DATA(:,:,1)
    
    !SVLSO2, SVLALT, SVLCLD 
    SVLSO2(:,:) = 0.0_RP
    SVLALT(:,:) = 0.0_RP
    SVLCLD(:,:) = 0.0_RP
    
    !HOORI
    do k = KS, KE
       call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
            ( HOORI(k,:,:),                                                      & ! [OUT]
              HOORI_DATA(k,:,:,12), HOORI_DATA(k,:,:,1:12), HOORI_DATA(k,:,:,1), & ! [IN]
              TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                  ) ! [IN]
    enddo
    
    do j = JS, JE
    do i = IS, IE
       SULDIO(i,j) = max( SULDIO(i,j), 0.0_RP )
       BBSO2 (i,j) = max( BBSO2 (i,j), 0.0_RP )
       VOLSO2(i,j) = max( VOLSO2(i,j), 0.0_RP )
       VOLALT(i,j) = max( VOLALT(i,j), 0.0_RP )
       ANSO4T(i,j) = max( ANSO4T(i,j), 0.0_RP )
       SVLSO2(i,j) = max( SVLSO2(i,j), 0.0_RP )
       SVLALT(i,j) = max( SVLALT(i,j), 0.0_RP )
       SVLCLD(i,j) = max( SVLCLD(i,j), 0.0_RP )
       do k = KS, KE
          HOORI(k,i,j) = max( HOORI(k,i,j), 0.0_RP )
       enddo
    enddo
    enddo
 
    !RADS, CEXTP, CEXTA, CBAKA
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
             RADS (k,i,j)       = RADSU (IRH)
             CEXTP(k,i,j,1:IWA) = CEXTSP(IRH,1:IWA)
             CEXTA(k,i,j,1:IWA) = CEXTSA(IRH,1:IWA)
             CBAKA(k,i,j,1:IWA) = CBAKSA(IRH,1:IWA)
          else
             BRH1 =  RH(k,i,j) - BRH(IRH)
             BRH2 = BRH(IRH+1) -  RH(k,i,j)
             BRH3 = BRH(IRH+1) - BRH(IRH)
             RADS(k,i,j) = ( BRH1 * RADSU(IRH+1) + BRH2 * RADSU(IRH) ) / BRH3
             do iw = 1, IWA
                CEXTP(k,i,j,iw) = ( BRH1 * CEXTSP(IRH+1,iw) + BRH2 * CEXTSP(IRH,iw) ) / BRH3
                CEXTA(k,i,j,iw) = ( BRH1 * CEXTSA(IRH+1,iw) + BRH2 * CEXTSA(IRH,iw) ) / BRH3
                CBAKA(k,i,j,iw) = ( BRH1 * CBAKSA(IRH+1,iw) + BRH2 * CBAKSA(IRH,iw) ) / BRH3
             enddo
          endif

       enddo
    enddo
    enddo

    !VTER
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          VTER(k,i,j) = 2.0_RP * (              ( RADSU(1)/RADS(k,i,j) )**3   * DENSSU &
                                   + ( 1.0_RP - ( RADSU(1)/RADS(k,i,j) )**3 ) * DENSW    ) &
                               * RADS(k,i,j)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA)**2 )
       enddo
    enddo
    enddo
  
    !================================================
    !end setup parameter

    
    !aerosol transport 
    !================================================
    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, 2
       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
    enddo

    
    !emission
    !-------------------------------
    VOLK (:,:) = KS
    VOLKS(:,:) = KS
    do j = JS, JE
    do i = IS, IE
       do k = KS+1, KE
          if ( VOLSO2(i,j) > 0.0_RP .and. VOLALT(i,j) > GDZM(k-1,i,j) + GDZS(i,j) .and. VOLK (i,j) == k-1 ) VOLK (i,j) = k
          if ( SVLSO2(i,j) > 0.0_RP .and. SVLALT(i,j) > GDZM(k-1,i,j) + GDZS(i,j) .and. VOLKS(i,j) == k-1 ) VOLKS(i,j) = k
       enddo
    enddo
    enddo

    KUPS(:,:) = VOLKS(:,:)
    do j = JS, JE
    do i = IS, IE
       do k = KS+1, KE
          if (       SVLSO2(i,j) > 0.0_RP                    .and. k         >= VOLKS(i,j)+1 &
               .and. SVLCLD(i,j) > GDZM(k-1,i,j) + GDZS(i,j) .and. KUPS(i,j) == k-1          ) KUPS(i,j) = k
       enddo
    enddo
    enddo
    KUPMXS = maxval( KUPS(IS:IE,JS:JE) )

    !DMS flux (Bates et al., Nature, 1987; Spiro et al., JGR, 1992)
    do j = JS, JE
    do i = IS, IE
       DMSFLX(i,j) = FACT_OCEAN(i,j) * ( 1.0_RP - GOICE(i,j) ) * (  1.08E-14_RP * RFLXSD(i,j) + 3.56E-13_RP ) &
                   + FACT_LAND (i,j) * 1.0E-12_RP * 62.0_RP / 32.0_RP / 60.0_RP &
                                     * (   2.20_RP * GRLAI(i,j) * max( RCOSZ(i,j), 0.0_RP ) * exp( 0.10_RP * (min(GDT(KS,i,j),DMSTMX)-298.0_RP) ) &
                                         + 0.58_RP                                          * exp( 0.14_RP * (min(GDT(KS,i,j),DMSTMX)-298.0_RP) ) )
       if ( DMSFLX(i,j) < EPS ) DMSFLX(i,j) = 0.0_RP
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       DPVOL = GDPM(VOLK (i,j)-1,i,j) - GDPM(VOLK(i,j),i,j)
       DPSVL = GDPM(VOLKS(i,j)-1,i,j) - GDPM(KUPS(i,j),i,j)
       ANTSO4(i,j) = ANSO4T(i,j) + SULDIO(i,j) * FSO4
       
       QEMIT(i,j,1) = DMSFLX(i,j) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP) ! DMS (ocean etc)
       QEMIT(i,j,2) = SULDIO(i,j) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP) ! SO2 (FF)
       QEMIT(i,j,3) = VOLSO2(i,j) / DPVOL      * GRAV * real(dt_AE,kind=RP) ! SO2 (volcanic)
       QEMIT(i,j,4) = ANTSO4(i,j) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP) ! SO4 (FF)
       QEMIT(i,j,5) = SVLSO2(i,j) / DPSVL      * GRAV * real(dt_AE,kind=RP) ! SO2 (volcanic;sporadic)
       QEMIT(i,j,6) =  BBSO2(i,j) / DPKBB(i,j) * GRAV * real(dt_AE,kind=RP) ! SO2 (BB)
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          if ( k <= KUP(i,j) ) then
             QTRC(k,i,j,1) = QTRC(k,i,j,1) + QEMIT(i,j,4) ! SO4 (FF)
             QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j,2) ! SO2 (FF)
             QTRC(k,i,j,3) = QTRC(k,i,j,3) + QEMIT(i,j,1) ! DMS (ocean, canopy and land)
          endif
          if ( k <= KBB  (i,j)                      ) QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j,6) ! SO2 (BB)
          if ( k == VOLK (i,j)                      ) QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j,3) ! SO2 (volcanic)
          if ( k >= VOLKS(i,j) .and. k <= KUPS(i,j) ) QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j,5) ! SO2 (volcanic;sporadic)
          QTRC(k,i,j,2) = QTRC(k,i,j,2) + AIRFEL(k,i,j) * AFFSO2 / DENS(k,i,j) * real(dt_AE,kind=RP)
       enddo
    enddo
    enddo
 
    do j = JS, JE
    do i = IS, IE
       SO2FLX(i,j) = SULDIO(i,j) + BBSO2(i,j) + VOLSO2(i,j) + SVLSO2(i,j)
       do k = KS, KE
          SO2FLX(i,j) = SO2FLX(i,j) + AIRFEL(k,i,j) * AFFSO2 * DELZ(k,i,j)
       enddo
    enddo
    enddo

    
    ! !instability
    ! !-------------------------------
    ! do iq = 1, QA_AE_SU
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB( QTRC(:,:,:,iq),                                & ! [MODIFIED]
    !                                               ONEWQ(1:KA-1,:,:), DELP(:,:,:), GDPM(0:KA,:,:) ) ! [IN]
    ! enddo

    
    !chemical reaction
    !-------------------------------
    !      k i j q 3
    ADDNUM(:,:,:)     = 0.0_RP
    CCNNUM(:,:,:)     = 0.0_RP
    NUMNW2(:,:,:,:,:) = 0.0_RP
    CHEDEP(:,:,:,:)   = 0.0_RP
    CHEDPT(:,:,:,:)   = 0.0_RP
    CHEDPM  (:,:,:)   = 0.0_RP
    SO2WTM  (:,:)     = 0.0_RP
    CLTORN        (:,:,:,:) = 0.0_RP
    CLTORN_CP     (:,:,:,:) = 0.0_RP
    CLTORN_MP_RAIN(:,:,:,:) = 0.0_RP
    CLTORN_MP_SNOW(:,:,:,:) = 0.0_RP
    QWTDEP        (:,:,:,:) = 0.0_RP
    QWTDEP_CP     (:,:,:,:) = 0.0_RP
    QWTDEP_MP_RAIN(:,:,:,:) = 0.0_RP
    QWTDEP_MP_SNOW(:,:,:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          !cloud water, O3, H2O2
          AIRNUM(k,i,j) = DENS(k,i,j) * 1.0E-3_RP * AVOG / MOLAIR ! [molecules/cm3]
          if ( TCLDF(k,i,j) > EPS ) then
             WATER(k,i,j) = QC(k,i,j) / TCLDF(k,i,j) * DENS(k,i,j) ! [kg/m3]
          else
             WATER(k,i,j) = 0.0_RP
          endif
          WATER(k,i,j) = max( WATER(k,i,j), 0.0_RP )
          NUMO3(k,i,j) = O3ORI(k,i,j) * AIRNUM(k,i,j) * 1.0E-6_RP ! [molecules/cm3]
          NUMHO(k,i,j) = HOORI(k,i,j) * AIRNUM(k,i,j) * 1.0E-6_RP ! [molecules/cm3]

          !reaction rate [cm3/molecules/sec]
          INREAC(1) = -234.0_RP / GDT(k,i,j)
          INREAC(2) = 7460.0_RP / GDT(k,i,j)
          INREAC(3) =  350.0_RP / GDT(k,i,j)
          INREAC(4) =  300.0_RP / GDT(k,i,j)
          ALPHA1 = 1.106E-31_RP * exp( INREAC(2) )  * AIRNUM(k,i,j)
          ALPHA2 = 3.000E-31_RP * INREAC(4)**3.3_RP * AIRNUM(k,i,j)
          ALPHA3 = 1.500E-12_RP
          REACRA(k,i,j,1) = 9.60E-12_RP * exp( INREAC(1) )
          REACRA(k,i,j,2) = 3.04E-12_RP * exp( INREAC(3) ) * ALPHA1 / ( 1.0_RP + ALPHA1 )
          REACRA(k,i,j,3) = ALPHA2 / ( 1.0_RP + ALPHA2/ALPHA3 ) & ! SO2-SO4(2-) (OUTSIDE CLOUD)
                          * ( 0.6_RP**( 1.0_RP / (1.0_RP+log10(ALPHA2/ALPHA3)**2) ) )
          
          !solubility constant [mol/L/atm]
          INEXP1 = 1.0_RP / GDT(k,i,j) - 1.0_RP / 298.0_RP
          INEXP2 = 3120.0_RP * INEXP1
          INEPO3 = 2560.0_RP * INEXP1
          INEPHO = 6600.0_RP * INEXP1
          SOL1 (k,i,j) = 1.23E+0_RP * exp( INEXP2 )
          SOLO3(k,i,j) = 1.15E-2_RP * exp( INEPO3 )
          SOLHO(k,i,j) = 9.70E+4_RP * exp( INEPHO )
          SOLOX1(k,i,j) = 1.0_RP + SOLHO(k,i,j) * 0.0821_RP * GDT(k,i,j) * WATER(k,i,j) * 1.0E-3_RP
          SOLOX2(k,i,j) = 1.0_RP + SOLO3(k,i,j) * 0.0821_RP * GDT(k,i,j) * WATER(k,i,j) * 1.0E-3_RP

          !equilibrium constant [mol/L]
          INEPQ1 = 2090.0_RP * INEXP1
          INEPQ2 = 1120.0_RP * INEXP1
          EQU1(k,i,j) = 1.7E-2_RP * exp( INEPQ1 )
          EQU2(k,i,j) = 6.0E-8_RP * exp( INEPQ2 )

          !number density of sulfate tracer
          do iq = 1, QA_AE_SU
             NUMORI(k,i,j,iq) = QTRC(k,i,j,iq) * DENS(k,i,j) * 1.0E-3_RP * AVOG / MOLM(iq) ! [molecules/cm3]
          enddo

       enddo
    enddo
    enddo

    !liquid phase reaction
    SO2G1(:,:,:) = NUMORI(:,:,:,2)
    call ATMOS_PHY_AE_SPRINTARS_AEROSU_RSO2AQ &
         ( SO2G1,                               & ! [MODIFIED]
           SULA, SO2A, PH,                      & ! [OUT]
           WATER, NUMHO, SOLOX1, NUMO3, SOLOX2, & ! [IN]
           EQU1, EQU2, SOL1, GDT, dt_AE, DELTCQ ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! CLTORN(k,i,j,2)   = RNWR(k,i,j) * SO2A(k,i,j) / AVOG * MOLM(2) * 1.0E+3_RP / DENS(k,i,j)
          ! QWTDEP(k,i,j,2)   = CLTORN(k,i,j,2)
          ! CHEDPT(k,i,j,2)   = SULA  (k,i,j)
          ! NUMNW2(k,i,j,1,2) = SULA  (k,i,j)
          ! NUMNW2(k,i,j,2,2) = SO2A  (k,i,j) * ( 1.0_RP - RNWR(k,i,j) )
          CLTORN_CP     (k,i,j,2) = RNWR_CP     (k,i,j) * SO2A(k,i,j) / AVOG * MOLM(2) * 1.0E+3_RP / DENS(k,i,j)
          CLTORN_MP_RAIN(k,i,j,2) = RNWR_MP_RAIN(k,i,j) * SO2A(k,i,j) / AVOG * MOLM(2) * 1.0E+3_RP / DENS(k,i,j)
          CLTORN_MP_SNOW(k,i,j,2) = RNWR_MP_SNOW(k,i,j) * SO2A(k,i,j) / AVOG * MOLM(2) * 1.0E+3_RP / DENS(k,i,j)
          QWTDEP_CP     (k,i,j,2) = CLTORN_CP     (k,i,j,2)
          QWTDEP_MP_RAIN(k,i,j,2) = CLTORN_MP_RAIN(k,i,j,2)
          QWTDEP_MP_SNOW(k,i,j,2) = CLTORN_MP_SNOW(k,i,j,2)
          CHEDPT(k,i,j,2)   = SULA(k,i,j)
          NUMNW2(k,i,j,1,2) = SULA(k,i,j)
          NUMNW2(k,i,j,2,2) = SO2A(k,i,j) * ( 1.0_RP - ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) ) )
       enddo
    enddo
    enddo

    !gas-phase reaction
    DMS1 (:,:,:) = NUMORI(:,:,:,3)
    DMS2 (:,:,:) = NUMORI(:,:,:,3)
    SO2G2(:,:,:) = NUMORI(:,:,:,2)
    call ATMOS_PHY_AE_SPRINTARS_AEROSU_QSSA &
         ( SO2G1, DMS1,                 & ! [MODIFIED]
           SULG1,                       & ! [OUT]
           REACRA, OHRAD, dt_AE, DELTCG ) ! [IN]
    
    call ATMOS_PHY_AE_SPRINTARS_AEROSU_QSSA &
         ( SO2G2, DMS2,                 & ! [MODIFIED]
           SULG2,                       & ! [OUT]
           REACRA, OHRAD, dt_AE, DELTCG ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          NUMNW (k,i,j,3)   = DMS1(k,i,j) * TCLDF(k,i,j) + DMS2(k,i,j) * ( 1.0_RP - TCLDF(k,i,j) )
          CHEDEP(k,i,j,3)   = NUMORI(k,i,j,3) - NUMNW(k,i,j,3)
          CHEDPT(k,i,j,3)   = SULG1(k,i,j)
          CHEDPT(k,i,j,1)   = SULG2(k,i,j)
          NUMNW2(k,i,j,1,3) = SULG1(k,i,j) + NUMORI(k,i,j,1)
          NUMNW2(k,i,j,1,1) = SULG2(k,i,j) + NUMORI(k,i,j,1)
          NUMNW2(k,i,j,2,3) = SO2G1(k,i,j)
          NUMNW2(k,i,j,2,1) = SO2G2(k,i,j)
       enddo
    enddo
    enddo

    !incloud
    !-------------------------------
    if ( SW_CCN ) then
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             ADDNUM(k,i,j) = NUMNW2(k,i,j,1,2) + NUMNW2(k,i,j,1,3)
             CCNNUM(k,i,j) = ADDNUM(k,i,j) * ACTISU(k,i,j)
          enddo
       enddo
       enddo
    else
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             ADDNUM(k,i,j) = NUMNW2(k,i,j,1,2) + NUMNW2(k,i,j,1,3)
             CCNNUM(k,i,j) = ADDNUM(k,i,j) * FINCSU
          enddo
       enddo
       enddo
    endif
    
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          !number density to mixing ratio
          QTRCTM2(k,i,j,1) = CCNNUM(k,i,j)     * TCLDF(k,i,j) / AVOG * MOLM( 1) * 1.0E+3_RP / DENS(k,i,j)
          QTRCTM2(k,i,j,2) = NUMNW2(k,i,j,2,2) * TCLDF(k,i,j) / AVOG * MOLM( 2) * 1.0E+3_RP / DENS(k,i,j)
          QTRCTM2(k,i,j,3) = 0.0_RP
          QTRCTM1(k,i,j,1) = (   NUMNW2(k,i,j,1,1)                 * ( 1.0_RP - TCLDF(k,i,j) ) &
                               + ( ADDNUM(k,i,j) - CCNNUM(k,i,j) ) * TCLDF(k,i,j)              ) &
                                                              / AVOG * MOLM( 1) * 1.0E+3_RP / DENS(k,i,j)
          QTRCTM1(k,i,j,2) = (   NUMNW2(k,i,j,2,1) * ( 1.0_RP - TCLDF(k,i,j) ) &
                               + NUMNW2(k,i,j,2,3) * TCLDF(k,i,j)              ) &
                                                              / AVOG * MOLM( 2) * 1.0E+3_RP / DENS(k,i,j)
          QTRCTM1(k,i,j,3) = NUMNW(k,i,j,3)                   / AVOG * MOLM( 3) * 1.0E+3_RP / DENS(k,i,j)

          !chemical reaction flux
          CHEDEP (k,i,j,2) = CHEDPT(k,i,j,1) * ( 1.0_RP - TCLDF(k,i,j) ) + ( CHEDPT(k,i,j,2) + CHEDPT(k,i,j,3) ) * TCLDF(k,i,j)
          SO2WTM   (i,j)   = SO2WTM(i,j) + CHEDPT(k,i,j,2) * TCLDF(k,i,j) * DELP(k,i,j) &
                                                              / AVOG * MOLM( 2) * 1.0E+3_RP / DENS(k,i,j)
          SO2C  = CHEDEP(k,i,j,2) * DELP(k,i,j) / GRAV / real(dt_AE,kind=RP) &
                                                              / AVOG * MOLM( 2) * 1.0E+3_RP / DENS(k,i,j) 
          SO2CW = CHEDPT(k,i,j,2) * DELP(k,i,j) / GRAV / real(dt_AE,kind=RP) * TCLDF(k,i,j) &
                                                              / AVOG * MOLM( 2) * 1.0E+3_RP / DENS(k,i,j)
          SO2CK (k,i,j) = max( SO2C,  0.0_RP )
          SO2CWK(k,i,j) = max( SO2CW, 0.0_RP )

       enddo
    enddo
    enddo

    do iq = 2, 3
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             CHEDPM(i,j,iq) = CHEDPM(i,j,iq) + CHEDEP(k,i,j,iq) * DELP(k,i,j) &
                                                              / AVOG * MOLM(iq) * 1.0E+3_RP / DENS(k,i,j)
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       SO2CHE(i,j) = max( CHEDPM(i,j,2) / GRAV / real(dt_AE,kind=RP), 0.0_RP )
       DMSCHE(i,j) = max( CHEDPM(i,j,3) / GRAV / real(dt_AE,kind=RP), 0.0_RP )
       SO2CHW(i,j) = max( SO2WTM(i,j)   / GRAV / real(dt_AE,kind=RP), 0.0_RP )
    enddo
    enddo

    
    ! !subcloud scavenging (wash out)
    ! !-------------------------------
    ! SURN(:,:,:,:) = 0.0_RP
    ! do j = JS, JE
    ! do i = IS, IE
    !    do k = KS, KE
    !       RAINN            = FRAIN(k,i,j) / VTR / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3 )
    !       SUN              = QTRCTM1(k,i,j,1) * DENS(k,i,j) / DRYMSU
    !       SURN   (k,i,j,1) = AA * PI * ( RADR + RADS(k,i,j) )**2 * ( VTR - VTER(k,i,j) ) * SUN * RAINN
    !       QWTDEP (k,i,j,1) = SURN(k,i,j,1) * DRYMSU / DENS(k,i,j) * real(dt_AE,kind=RP)
    !       QTRCTMP(k,i,j)   = QTRCTM1(k,i,j,1) - QWTDEP(k,i,j,1)
    !       if ( QTRCTMP(k,i,j) < 0.0_RP ) then
    !          QTRCTMP(k,i,j)   = 0.0_RP
    !          SURN   (k,i,j,1) = SUN / real(dt_AE,kind=RP)
    !          QWTDEP (k,i,j,1) = QTRCTM1(k,i,j,1)
    !       endif
    !    enddo
    ! enddo
    ! enddo

    
    ! !incloud scavenging (rain out)
    ! !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    do k = KS, KE
    !       CLTORN (k,i,j,1) = RNWR(k,i,j) * QTRCTM2(k,i,j,1)
    !       QTRCTM2(k,i,j,1) = QTRCTM2(k,i,j,1) - CLTORN(k,i,j,1)
    !       QWTDEP (k,i,j,1) = QWTDEP (k,i,j,1) + CLTORN(k,i,j,1)
    !       QTRCTM1(k,i,j,1) = QTRCTMP(k,i,j)
    !    enddo
    ! enddo
    ! enddo

    ! !re-emission from rain
    ! !-------------------------------
    ! MSO2 = MOLM(2) / AVOG
    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER                         &
    !      ( QTRCTM1(:,:,:,1), QWTDEP(:,:,:,1),                         & ! [MODIFIED]
    !        WETF(:,:,1),                                               & ! [OUT]
    !        SURN(:,:,:,1), DRYMSU, CLTORN(:,:,:,1),                    & ! [IN]
    !        DENS(:,:,:), FRAIN(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER                         &
    !      ( QTRCTM1(:,:,:,2), QWTDEP(:,:,:,2),                         & ! [MODIFIED]
    !        WETF(:,:,2),                                               & ! [OUT]
    !        SURN(:,:,:,2), MSO2, CLTORN(:,:,:,2),                      & ! [IN]
    !        DENS(:,:,:), FRAIN(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    ! do iq = 1, QA_AE_SU
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
    ! do iq = 1, QA_AE_SU
       
    !    if ( iq == 1 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(1) * CDVE(i,j) / ( DDV(1) + CDVE(i,j) )
    !       enddo
    !       enddo
    !    elseif ( iq == 2 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(2) * CDVE(i,j) / ( DDV(2) + CDVE(i,j) ) * FACT_OCEAN(i,j)              &
    !                     + DDV(3) * CDVE(i,j) / ( DDV(3) + CDVE(i,j) ) * ( 1.0_RP - FACT_OCEAN(i,j) )
    !       enddo
    !       enddo
    !    elseif ( iq == 3 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(4) * CDVE(i,j) / ( DDV(4) + CDVE(i,j) )
    !       enddo
    !       enddo
    !    endif
       
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !         ( QMTX,                    & ! [OUT]
    !           DENS, DELP, CDVES, dt_AE ) ! [IN]
       
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),                & ! [MODIFIED]
    !           DRYF  (:,:,iq),                & ! [OUT]
    !           QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
       
    ! enddo
    

    ! !gravitational settling
    ! !-------------------------------
    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !      ( QTRC(:,:,:,1),                 & ! [MODIFIED]
    !        QLOAD(:,:,:,1), SULGRV(:,:),   & ! [OUT]
    !        GDPM, DELP, DENS, VTER, dt_AE  ) ! [IN]
    ! QLOAD(:,:,:,2) = QTRC(:,:,:,2)

    
    !subcloud scavenging (wash out)
    !-------------------------------
    SURN_CP     (:,:,:,2) = 0.0_RP
    SURN_MP_RAIN(:,:,:,2) = 0.0_RP
    SURN_MP_SNOW(:,:,:,2) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! RAINN_CP      = FRAIN_CP     (k,i,j) / VTR                / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3                )
          ! RAINN_MP_RAIN = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
          ! RAINN_MP_SNOW = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
          SUN           = QTRCTM1(k,i,j,1) * DENS(k,i,j) / DRYMSU
          SURN_CP     (k,i,j,1) = AA * PI * ( RADR                + RADS(k,i,j) )**2 * abs( VTR                - VTER(k,i,j) ) * SUN * RAINN_CP     (k,i,j)
          SURN_MP_RAIN(k,i,j,1) = AA * PI * ( RADR_MP_RAIN(k,i,j) + RADS(k,i,j) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j) ) * SUN * RAINN_MP_RAIN(k,i,j)
          SURN_MP_SNOW(k,i,j,1) = AA * PI * ( RADR_MP_SNOW(k,i,j) + RADS(k,i,j) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j) ) * SUN * RAINN_MP_SNOW(k,i,j)
          QWTDEP_CP     (k,i,j,1) = SURN_CP     (k,i,j,1) * DRYMSU / DENS(k,i,j) * real(dt_AE,kind=RP)
          QWTDEP_MP_RAIN(k,i,j,1) = SURN_MP_RAIN(k,i,j,1) * DRYMSU / DENS(k,i,j) * real(dt_AE,kind=RP)
          QWTDEP_MP_SNOW(k,i,j,1) = SURN_MP_SNOW(k,i,j,1) * DRYMSU / DENS(k,i,j) * real(dt_AE,kind=RP)
          QTRCTMP(k,i,j)   = QTRCTM1(k,i,j,1) - QWTDEP_CP(k,i,j,1) - QWTDEP_MP_RAIN(k,i,j,1) - QWTDEP_MP_SNOW(k,i,j,1)
          if ( QTRCTMP(k,i,j) < 0.0_RP ) then
             QTRCTMP(k,i,j)   = 0.0_RP
             SURN   (k,i,j,1) = SURN_CP(k,i,j,1) + SURN_MP_RAIN(k,i,j,1) + SURN_MP_SNOW(k,i,j,1)
             SURN_CP     (k,i,j,1) = SUN / real(dt_AE,kind=RP) * SURN_CP     (k,i,j,1) / SURN(k,i,j,1)
             SURN_MP_RAIN(k,i,j,1) = SUN / real(dt_AE,kind=RP) * SURN_MP_RAIN(k,i,j,1) / SURN(k,i,j,1)
             SURN_MP_SNOW(k,i,j,1) = SUN / real(dt_AE,kind=RP) * SURN_MP_SNOW(k,i,j,1) / SURN(k,i,j,1)
             QWTDEP_CP     (k,i,j,1) = QTRCTM1(k,i,j,1) * SURN_CP     (k,i,j,1) / SURN(k,i,j,1)
             QWTDEP_MP_RAIN(k,i,j,1) = QTRCTM1(k,i,j,1) * SURN_MP_RAIN(k,i,j,1) / SURN(k,i,j,1)
             QWTDEP_MP_SNOW(k,i,j,1) = QTRCTM1(k,i,j,1) * SURN_MP_SNOW(k,i,j,1) / SURN(k,i,j,1)
          endif
       enddo
    enddo
    enddo

    
    !incloud scavenging (rain out)
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          CLTORN_CP     (k,i,j,1) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,1)
          CLTORN_MP_RAIN(k,i,j,1) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,1)
          CLTORN_MP_SNOW(k,i,j,1) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,1)
          CLTORN        (k,i,j,1) = CLTORN_CP(k,i,j,1) + CLTORN_MP_RAIN(k,i,j,1) + CLTORN_MP_SNOW(k,i,j,1)
          if ( QTRCTM2(k,i,j,1) - CLTORN(k,i,j,1) < 0.0_RP ) then
             CLTORN_CP     (k,i,j,1) = QTRCTM2(k,i,j,1) * RNWR_CP     (k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
             CLTORN_MP_RAIN(k,i,j,1) = QTRCTM2(k,i,j,1) * RNWR_MP_RAIN(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
             CLTORN_MP_SNOW(k,i,j,1) = QTRCTM2(k,i,j,1) * RNWR_MP_SNOW(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
             QTRCTM2(k,i,j,1) = 0.0_RP
          else
             QTRCTM2(k,i,j,1) = QTRCTM2(k,i,j,1) - CLTORN(k,i,j,1)
          endif
          QWTDEP_CP     (k,i,j,1) = QWTDEP_CP     (k,i,j,1) + CLTORN_CP     (k,i,j,1)
          QWTDEP_MP_RAIN(k,i,j,1) = QWTDEP_MP_RAIN(k,i,j,1) + CLTORN_MP_RAIN(k,i,j,1)
          QWTDEP_MP_SNOW(k,i,j,1) = QWTDEP_MP_SNOW(k,i,j,1) + CLTORN_MP_SNOW(k,i,j,1)
          QTRCTM1       (k,i,j,1) = QTRCTMP(k,i,j)
       enddo
    enddo
    enddo

    
    !dry deposition
    !-------------------------------
    do iq = 1, QA_AE_SU
       
       if ( iq == 1 ) then
          do j = JS, JE
          do i = IS, IE
             CDVES(i,j) = DDV(1) * CDVE(i,j) / ( DDV(1) + CDVE(i,j) )
          enddo
          enddo
       elseif ( iq == 2 ) then
          do j = JS, JE
          do i = IS, IE
             CDVES(i,j) = DDV(2) * CDVE(i,j) / ( DDV(2) + CDVE(i,j) ) * FACT_OCEAN(i,j)              &
                        + DDV(3) * CDVE(i,j) / ( DDV(3) + CDVE(i,j) ) * ( 1.0_RP - FACT_OCEAN(i,j) )
          enddo
          enddo
       elseif ( iq == 3 ) then
          do j = JS, JE
          do i = IS, IE
             CDVES(i,j) = DDV(4) * CDVE(i,j) / ( DDV(4) + CDVE(i,j) )
          enddo
          enddo
       endif
       
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
            ( QMTX,                    & ! [OUT]
              DENS, DELP, CDVES, dt_AE ) ! [IN]
       
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
            ( QTRCTM1(:,:,:,iq),             & ! [MODIFIED]
              DRYF(:,:,iq),                  & ! [OUT]
              QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
       
    enddo
    

    !gravitational settling
    !-------------------------------
    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
         ( QTRCTM1(:,:,:,1),              & ! [MODIFIED]
           QLOAD(:,:,:,1), SULGRV(:,:),   & ! [OUT]
           GDPM, DELP, DENS, VTER, dt_AE  ) ! [IN]
    QLOAD(:,:,:,2) = QTRCTM1(:,:,:,2)
    QLOAD(:,:,:,3) = QTRCTM1(:,:,:,3)


    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, 2
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
    MSO2 = MOLM(2) / AVOG
    call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
         ( QTRCTM1(:,:,:,1), QWTDEP_CP(:,:,:,1),                         & ! [MODIFIED]
           WETF_CP(:,:,1),                                               & ! [OUT]
           SURN_CP(:,:,:,1), DRYMSU, CLTORN_CP(:,:,:,1),                 & ! [IN]
           DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
         ( QTRCTM1(:,:,:,2), QWTDEP_CP(:,:,:,2),                         & ! [MODIFIED]
           WETF_CP(:,:,2),                                               & ! [OUT]
           SURN_CP(:,:,:,2), MSO2, CLTORN_CP(:,:,:,2),                   & ! [IN]
           DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    
    !re-emission from rain & snow : MP
    !-------------------------------
    do iq = 1, 2
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
    do iq = 1, QA_AE_SU
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
    ! do iq = 1, QA_AE_SU
       
    !    if ( iq == 1 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(1) * CDVE(i,j) / ( DDV(1) + CDVE(i,j) )
    !       enddo
    !       enddo
    !    elseif ( iq == 2 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(2) * CDVE(i,j) / ( DDV(2) + CDVE(i,j) ) * FACT_OCEAN(i,j)              &
    !                     + DDV(3) * CDVE(i,j) / ( DDV(3) + CDVE(i,j) ) * ( 1.0_RP - FACT_OCEAN(i,j) )
    !       enddo
    !       enddo
    !    elseif ( iq == 3 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(4) * CDVE(i,j) / ( DDV(4) + CDVE(i,j) )
    !       enddo
    !       enddo
    !    endif
       
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !         ( QMTX,                    & ! [OUT]
    !           DENS, DELP, CDVES, dt_AE ) ! [IN]
       
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRC(:,:,:,iq),                & ! [MODIFIED]
    !           DRYF  (:,:,iq),                & ! [OUT]
    !           QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
       
    ! enddo
    

    ! !gravitational settling
    ! !-------------------------------
    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !      ( QTRC(:,:,:,1),                 & ! [MODIFIED]
    !        QLOAD(:,:,:,1), SULGRV(:,:),   & ! [OUT]
    !        GDPM, DELP, DENS, VTER, dt_AE  ) ! [IN]
    ! QLOAD(:,:,:,2) = QTRC(:,:,:,2)
    ! QLOAD(:,:,:,3) = QTRC(:,:,:,3)

    !================================================
    !end aerosol transport

    
    !aerosol
    !================================================
    !Reset
    !-------------------------------
    COLMS(:,:,:) = 0.0_RP

    !-------------------------------
    do iq = 1, QA_AE_SU
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             SCON (k,i,j,iq) = QLOAD(k,i,j,iq) * DENS(k,i,j)
             COLMS  (i,j,iq) = COLMS  (i,j,iq) + QLOAD(k,i,j,iq) * DELP(k,i,j) / GRAV
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          OUTQLD(k,i,j)   = QLOAD (k,i,j,1)
          NUMSUL(k,i,j)   = SCON  (k,i,j,1) * MTONSU
          NUMSO2(k,i,j)   = SCON  (k,i,j,2) * 1.0E+3_RP * AVOG / MOLM(2)
          NUMDMS(k,i,j)   = SCON  (k,i,j,3) * 1.0E+3_RP * AVOG / MOLM(3)
          SPPT  (k,i,j,1) = SCON  (k,i,j,1) * 1.0E+3_RP * AVOG / MOLM(1) / AIRNUM(k,i,j) * 1.0E+6_RP
          SPPT  (k,i,j,2) = NUMSO2(k,i,j) / AIRNUM(k,i,j) * 1.0E+6_RP
          SPPT  (k,i,j,3) = NUMDMS(k,i,j) / AIRNUM(k,i,j) * 1.0E+6_RP
          PM25SU(k,i,j)   = SCON  (k,i,j,1)
          WATRSU(k,i,j)   = ( RADS(k,i,j)**3 / RADSU(1)**3 - 1 ) * QLOAD(k,i,j,1) * DENSW / DENSSU
       enddo
    enddo
    enddo
    !================================================
    !end aerosol

    
    !flux 
    !================================================
    do j = JS, JE
    do i = IS, IE
       SULTDP(i,j) = WETF  (i,j,1) + DRYF(i,j,1) + SULGRV(i,j)
       SO2TDP(i,j) = SO2CHE(i,j)   + WETF(i,j,2) + DRYF  (i,j,2)
       DMSTDP(i,j) = DMSCHE(i,j)   + DRYF(i,j,3)
    enddo
    enddo
    !================================================
    !end flux

    
    !optical parameters
    !================================================
    !Reset
    !-------------------------------
    !      i j w 2
    TAUSU (:,:,:,:) = 0.0_RP
    TAUSUD(:,:)     = 0.0_RP
    
    !all sky
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          TUSE = QLOAD(k,i,j,1) * DELP(k,i,j) / GRAV
          do iw = 1, IWA
             TAUSU   (i,j,iw,1) = TAUSU(i,j,iw,1) + TUSE * CEXTP(k,i,j,iw)
             CEXTSU(k,i,j,iw,1) = TUSE * CEXTA(k,i,j,iw) / DELZ(k,i,j)
             CBAKSU(k,i,j,iw,1) = TUSE * CBAKA(k,i,j,iw) / DELZ(k,i,j)
          enddo
          TAUSUD  (i,j) = TAUSUD(i,j) + TUSE * CEXTSP(1,1)
          CEXTSD(k,i,j) = TUSE * CEXTSA(1,1) / DELZ(k,i,j)
       enddo
    enddo
    enddo

    !clear sky diagnoses
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
          
       do k = KS, KE
          if ( TCLDF2(k,i,j) >= CCMAX ) then
             NUMSUC(k,i,j)            = CONST_UNDEF !VMISS
             SCONC (k,i,j,1:QA_AE_SU) = CONST_UNDEF !VMISS
             CEXTSU(k,i,j,1:IWA,2)    = CONST_UNDEF !VMISS
             CBAKSU(k,i,j,1:IWA,2)    = CONST_UNDEF !VMISS
          else
             NUMSUC(k,i,j)            = NUMSUL(k,i,j)
             SCONC (k,i,j,1:QA_AE_SU) = SCON  (k,i,j,1:QA_AE_SU)
             CEXTSU(k,i,j,1:IWA,2)    = CEXTSU(k,i,j,1:IWA,1)
             CBAKSU(k,i,j,1:IWA,2)    = CBAKSU(k,i,j,1:IWA,1)
          endif
       enddo

       if ( CCOVMR(i,j) >= CCMAX ) then
          COLMSC(i,j,1:QA_AE_SU) = CONST_UNDEF !VMISS
          TAUSU (i,j,1:IWA,2)    = CONST_UNDEF !VMISS
       else
          COLMSC(i,j,1:QA_AE_SU) = COLMS(i,j,1:QA_AE_SU)
          TAUSU (i,j,1:IWA,2)    = TAUSU(i,j,1:IWA,1)
       endif
       
    enddo
    enddo
    !================================================
    !end optical parameters

    
    !OUTPUT
    !================================================
    !flux
    !-------------------------------
    !                                  k i j q
    call FILE_HISTORY_in( DMSFLX        (:,:),   'DMSFLX_SPRINTARS',         'DMS emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SO2FLX        (:,:),   'SO2FLX_SPRINTARS',         'SO2 emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( ANTSO4        (:,:),   'SO4FLX_SPRINTARS',         'SO4 emission flux',                      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SULTDP        (:,:),   'TDFTSU_SPRINTARS',         'total deposition flux (sulfate)',        'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF          (:,:,1), 'WEFTSU_SPRINTARS',         'wet deposition flux (sulfate)',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_CP       (:,:,1), 'WEFTSU_CP_SPRINTARS',      'wet deposition flux (sulfate;CP)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_MP_RAIN  (:,:,1), 'WEFTSU_MP_RAIN_SPRINTARS', 'wet deposition flux (sulfate;MP(RAIN))', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_MP_SNOW  (:,:,1), 'WEFTSU_MP_SNOW_SPRINTARS', 'wet deposition flux (sulfate;MP(SNOW))', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYF          (:,:,1), 'DRFTSU_SPRINTARS',         'dry deposition flux (sulfate)',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SULGRV        (:,:),   'GRFTSU_SPRINTARS',         'gravitational settling flux (sulfate)',  'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SO2TDP        (:,:),   'TDFSO2_SPRINTARS',         'total deposition flux (SO2)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SO2CHE        (:,:),   'CFSO2_SPRINTARS',          'chemical reaction flux (SO2)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SO2CHW        (:,:),   'CFSO2W_SPRINTARS',         'chemical reaction flux (SO2;wet)',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF          (:,:,2), 'WEFSO2_SPRINTARS',         'wet deposition flux (SO2)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_CP       (:,:,2), 'WEFSO2_CP_SPRINTARS',      'wet deposition flux (SO2;CP)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_MP_RAIN  (:,:,2), 'WEFSO2_MP_RAIN_SPRINTARS', 'wet deposition flux (SO2;MP(RAIN))',     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETF_MP_SNOW  (:,:,2), 'WEFSO2_MP_SNOW_SPRINTARS', 'wet deposition flux (SO2;MP(SNOW))',     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYF          (:,:,2), 'DRFSO2_SPRINTARS',         'dry deposotion flux (SO2)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DMSTDP        (:,:),   'TDFDMS_SPRINTARS',         'total deposition flux (DMS)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DMSCHE        (:,:),   'CFDMS_SPRINTARS',          'chemical reaction flux (DMS)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYF          (:,:,3), 'DRFDMS_SPRINTARS',         'dry deposition flux (DMS)',              'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( SO2CK       (:,:,:),   'CSO2K_SPRINTARS',          'chemical reaction flux (SO2)',           'kg/m2/s', dim_type = 'ZXY' )
    call FILE_HISTORY_in( SO2CWK      (:,:,:),   'CSO2WK_SPRINTARS',         'chemical reaction flux (SO2;wet)',       'kg/m2/s', dim_type = 'ZXY' )

    !aerosol
    !-------------------------------
    do iq = 1, QA_AE_SU
       write(cnumber,'(i2.2)') iq
       !                           k i j q
       call FILE_HISTORY_in( SCON (:,:,:,iq), 'CONS'//cnumber//'_SPRINTARS', &
                         'mass concentration ('//trim(adjustl(ctype(iq)))//')',   'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( SPPT (:,:,:,iq), 'SPPT'//cnumber//'_SPRINTARS', &
                         'volume concentration ('//trim(adjustl(ctype(iq)))//')', 'pptv',  dim_type = 'ZXY' )
       call FILE_HISTORY_in( COLMS  (:,:,iq), 'COLS'//cnumber//'_SPRINTARS', &
                         'column loading ('//trim(adjustl(ctype(iq)))//')',       'kg/m2', dim_type = 'XY'  )
    enddo
    !                            k i j q
    call FILE_HISTORY_in( QLOAD (:,:,:,1), 'QLDST_SPRINTARS',  'mixing ratio (sulfate)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMSUL(:,:,:),   'NMCST_SPRINTARS',  'number concentration (sulfate)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMSO2(:,:,:),   'NMCSO2_SPRINTARS', 'number concentration (SO2)',     '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMDMS(:,:,:),   'NMCDMS_SPRINTARS', 'number concentration (DMS)',     '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( PH    (:,:,:),   'PH_SPRINTARS',     'pH',                             '1',     dim_type = 'ZXY' )
    call FILE_HISTORY_in( HOORI (:,:,:),   'HOORI_SPRINTARS',  'H2O2',                           'ppm',   dim_type = 'ZXY' )
    
    !optical properties
    !-------------------------------
    !                            k i j w 2
    call FILE_HISTORY_in( TAUSU   (:,:,1,1), 'TAUSU_SPRINTARS',  'AOT (sulfate;550nm;allsky)',          '1',      dim_type = 'XY'  ) !  550nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,1,1), 'EXTSU1_SPRINTARS', 'extinction (sulfate;532nm;allsky)',   '1/m',    dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,1,1), 'BAKSU1_SPRINTARS', 'backscatter (sulfate;532nm;allsky)',  '1/m/sr', dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,2,1), 'EXTSU2_SPRINTARS', 'extinction (sulfate;1064nm;allsky)',  '1/m',    dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,2,1), 'BAKSU2_SPRINTARS', 'backscatter (sulfate;1064nm;allsky)', '1/m/sr', dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,3,1), 'EXTSU3_SPRINTARS', 'extinction (sulfate;355nm;allsky)',   '1/m',    dim_type = 'ZXY' ) !  355nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,3,1), 'BAKSU3_SPRINTARS', 'backscatter (sulfate;355nm;allsky)',  '1/m/sr', dim_type = 'ZXY' ) !  355nm
    
    !clear sky diagnoses
    !-------------------------------
    do iq = 1, QA_AE_SU
       write(cnumber,'(i2.2)') iq
       !                            k i j q
       call FILE_HISTORY_in( SCONC (:,:,:,iq), 'CNSC'//cnumber//'_SPRINTARS', &
                        'mass concentration ('//trim(adjustl(ctype(iq)))//';clearsky)', 'kg/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( COLMSC  (:,:,iq), 'CLSC'//cnumber//'_SPRINTARS', &
                        'column loading ('//trim(adjustl(ctype(iq)))//';clearsky)',     'kg/m2', dim_type = 'XY'  )
    enddo
    !                            k i j w 2
    call FILE_HISTORY_in( NUMSUC(:,:,:),     'NMCSTC_SPRINTARS',  'number concentration (sulfate;clearsky)', '1/m3',   dim_type = 'ZXY' )
    call FILE_HISTORY_in( TAUSU   (:,:,1,2), 'TAUSUC_SPRINTARS',  'AOT (sulfate;550nm;clearsky)',            '1',      dim_type = 'XY'  ) !  550nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,1,2), 'EXTSUC1_SPRINTARS', 'extinction (sulfate;532nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,1,2), 'BAKSUC1_SPRINTARS', 'backscatter (sulfate;532nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !  532nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,2,2), 'EXTSUC2_SPRINTARS', 'extinction (sulfate;1064nm;clearsky)',    '1/m',    dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,2,2), 'BAKSUC2_SPRINTARS', 'backscatter (sulfate;1064nm;clearsky)',   '1/m/sr', dim_type = 'ZXY' ) ! 1064nm
    call FILE_HISTORY_in( CEXTSU(:,:,:,3,2), 'EXTSUC3_SPRINTARS', 'extinction (sulfate;355nm;clearsky)',     '1/m',    dim_type = 'ZXY' ) !  355nm
    call FILE_HISTORY_in( CBAKSU(:,:,:,3,2), 'BAKSUC3_SPRINTARS', 'backscatter (sulfate;355nm;clearsky)',    '1/m/sr', dim_type = 'ZXY' ) !  355nm
    !================================================
    !end OUTPUT

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROSU
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_aerosulf
