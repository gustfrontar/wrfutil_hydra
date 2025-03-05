!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - carbon
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_aerocarb
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
       PI     => CONST_PI,   & !< pi
       GRAV   => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       RCON   => CONST_R,    & !< universal gas constant = 8.31436 [J/mol/K]
       MOLAIR => CONST_Mdry, & !< mass weight (dry air) = 28.966   [g/mol]
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
       NCA,                                           & !< total number of carbon tracers
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
       AVOG,                                          & !< Avogadro's number
       THRESP                                           !< precipitation threshold  [kg/m2/s]
       
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_AEROCA_setup
  public :: ATMOS_PHY_AE_SPRINTARS_AEROCA_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_AEROCA
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: ATMOS_PHY_AE_SPRINTARS_AEROCA_FSOA !< offline SOA

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  ! parameter
  integer,  private, parameter :: I_intOC = 1 !< internal mixed OC
  integer,  private, parameter :: I_intBC = 2 !< internal mixed BC
  integer,  private, parameter :: I_SOC   = 3 !< secondary organic carbon (SOC) from NVOC
  integer,  private, parameter :: I_extBC = 4 !< external mixed BC
  integer,  private, parameter :: I_terp  = 5 !< terpene
  integer,  private, parameter :: I_isop  = 6 !< isoprene
  
  integer,  private, parameter :: NCA1           = 3 !< 1->internal mixed OC & BC, 2->SOC from NVOC, 3->external mixed BC
  integer,  private, parameter :: NSOA           = 2 !< 1->terpene, 2->isoprene
  integer,  private, parameter :: IQTRCC_intBCOC = 1 !< internal mixed OC & BC
  integer,  private, parameter :: IQTRCC_SOC     = 2 !< SOC from NVOC
  integer,  private, parameter :: IQTRCC_extBC   = 3 !< external mixed BC
  integer,  private, parameter :: IQTRCC_terp    = 4 !< terpene
  integer,  private, parameter :: IQTRCC_isop    = 5 !< isoprene
  integer,  private, parameter :: ISOA_terp      = 1 !< terpene
  integer,  private, parameter :: ISOA_isop      = 2 !< isoprene
    
  integer,  private, parameter :: NRBCOM   = 4 !< BC/OM = 0.30, 0.15, 0.00, & BC
  integer,  private, parameter :: I_BCOM30 = 1 !< BC/OM = 0.30
  integer,  private, parameter :: I_BCOM15 = 2 !< BC/OM = 0.15
  integer,  private, parameter :: I_BCOM00 = 3 !< BC/OM = 0.00
  integer,  private, parameter :: I_BCOMBC = 4 !< BC

  !SOA VOCs
  integer,  parameter :: NVOC = 10 !< No. of SOA VOCs

  ! variables
  real(RP), private, allocatable, save :: BBBC_DATA  (:,:,:,:) !< BC emission (biomass burning)       [kg(BC)/m2/s] (IA,JA,year,month)
  real(RP), private, allocatable, save :: BBOC_DATA  (:,:,:,:) !< OC emission (biomass burning)       [kg(OC)/m2/s] (IA,JA,year,month)
  real(RP), private, allocatable, save :: ANTBC_DATA (:,:,:,:) !< BC emission (fossil fuel)           [kg(BC)/m2/s] (IA,JA,year,month)
  real(RP), private, allocatable, save :: ANTOC_DATA (:,:,:,:) !< OC emission (fossil fuel)           [kg(OC)/m2/s] (IA,JA,year,month)
  real(RP), private, allocatable, save :: ISOP_DATA  (:,:,  :) !< isoprene emission                   [kg(C)/m2/s]  (IA,JA,     month) for OPT_OFFLINE_SOA
  real(RP), private, allocatable, save :: TERP_DATA  (:,:,  :) !< terpene  emission                   [kg(C)/m2/s]  (IA,JA,     month) for OPT_OFFLINE_SOA
  real(RP), private, allocatable, save :: CHLOA_DATA (:,:,  :) !< Chlorophyll-A (CHL-A) concentration [mg/m3]       (IA,JA,     month) for OPT_OFFLINE_SOA
  real(RP), private, allocatable, save :: DIFATT_DATA(:,:,  :) !< diffuse attenuation coef. at 490nm  [1/m]         (IA,JA,     month) for OPT_OFFLINE_SOA

  real(RP), private, allocatable, save :: GX(:,:,:,:,:)        !< for SOA 
  real(RP), private, allocatable, save :: QTRCINP(:,:,:,:)     !< mixing ratio in prec(MP)
  
  ! namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA /
  real(RP), private, save :: FINCCA(NCA1) = (/ 0.3E+0_RP,  0.3E+0_RP,  0.1E+0_RP /) !< incloud coefficient
  
  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Carbon setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_setup &
       ( LON_REAL, LAT_REAL, fileid ) ! [IN]
    implicit none
    real(RP),     intent(in) :: LON_REAL(IA,JA) !< longitude [rad]
    real(RP),     intent(in) :: LAT_REAL(IA,JA) !< latitude  [rad]
    integer,      intent(in) :: fileid
    character(len=H_MID) :: varname
    real(RP) :: work2d(IA,JA), work3d(KA,IA,JA)
    integer :: iq, ip 
    integer :: ierr
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA / &
         FINCCA
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerocarb",*) 'Setup'
    
    !--- read namelist
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerocarb",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_aerocarb",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCA_NMAECA)


    !Setup BBBC, BBOC, ANTBC, ANTOC, ISOP, TERP, CHLOA, DIFATT
    allocate(   BBBC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(   BBOC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(  ANTBC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(  ANTOC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(   ISOP_DATA(IA,JA,   12) )
    allocate(   TERP_DATA(IA,JA,   12) )
    allocate(  CHLOA_DATA(IA,JA,   12) )
    allocate( DIFATT_DATA(IA,JA,   12) )

    do iq = 1, 22
       write(varname,'(A,I0)') 'BBBC_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          BBBC_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    do iq = 1, 22
       write(varname,'(A,I0)') 'BBOC_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          BBOC_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    do iq = 1, 22
       write(varname,'(A,I0)') 'ANTBC_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          ANTBC_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    do iq = 1, 22
       write(varname,'(A,I0)') 'ANTOC_', iq
       do ip = 1, 12
          call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                                  work2d(:,:),       &
                                  ip                           )
          ANTOC_DATA(:,:,iq,ip) = work2d(:,:)
       enddo
    enddo
    write(varname,'(A)') 'ISOP'
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                               work2d(:,:),       &
                               ip                           )
       ISOP_DATA(:,:,ip) = work2d(:,:)
    enddo
    write(varname,'(A)') 'TERP'
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                               work2d(:,:),       &
                               ip                           )
       TERP_DATA(:,:,ip) = work2d(:,:)
    enddo
    write(varname,'(A)') 'CHLOA'
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                               work2d(:,:),       &
                               ip                           )
       CHLOA_DATA(:,:,ip) = work2d(:,:)
    enddo
    write(varname,'(A)') 'DIFATT'
    do ip = 1, 12
       call FILE_CARTESC_read( fileid, trim(varname), 'XY', &
                               work2d(:,:),       &
                               ip                           )
       DIFATT_DATA(:,:,ip) = work2d(:,:)
    enddo

    allocate( GX(KA,IA,JA,NVOC,2) )
    GX(:,:,:,:,:) = 0.0_RP
    
    allocate( QTRCINP(KA,IA,JA,4) )
    QTRCINP(:,:,:,:) = 0.0_RP
 
    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP
    call FILE_HISTORY_in( zero_2d(:,:),   'EFBC_SPRINTARS',           'emission flux (BC)',                         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'EFOC_SPRINTARS',           'emission flux (OM(OC))',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'EFSSA_SPRINTARS',          'emission flux (oceanic OM)',                 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'CHESOA_SPRINTARS',         'secondary organic aerosol (SOA) production', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TERP_SPRINTARS',           'emission flux (terpene) ',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'ISOP_SPRINTARS',           'emission flux (isoprene)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'EMIISO1_SPRINTARS',        'emission flux (terpene;oceanic precursor)',  'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'EMIISO2_SPRINTARS',        'emission flux (isoprene;oceanic precursor)', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTCA_SPRINTARS',         'wet deposition flux (BC+OM(OC))',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTCA_CP_SPRINTARS',      'wet deposition flux (BC+OM(OC);CP)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTCA_MP_RAIN_SPRINTARS', 'wet deposition flux (BC+OM(OC);MP(RAIN))',   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTCA_MP_SNOW_SPRINTARS', 'wet deposition flux (BC+OM(OC);MP(SNOW))',   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'DRFTCA_SPRINTARS',         'dry deposition flux (BC+OM(OC))',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'GRFTCA_SPRINTARS',         'gravitational settling flux (BC+OM(OC))',    'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TDFTCA_SPRINTARS',         'total deposition flux (BC+OM(OC))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTOC_SPRINTARS',         'wet deposition flux (OM(OC))',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTOC_CP_SPRINTARS',      'wet deposition flux (OM(OC);CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTOC_MP_RAIN_SPRINTARS', 'wet deposition flux (OM(OC);MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTOC_MP_SNOW_SPRINTARS', 'wet deposition flux (OM(OC);MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'DRFTOC_SPRINTARS',         'dry deposition flux (OM(OC))',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'GRFTOC_SPRINTARS',         'gravitational settling flux (OM(OC))',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TDFTOC_SPRINTARS',         'total deposition flux (OM(OC))',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTBC_SPRINTARS',         'wet deposition flux (BC)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTBC_CP_SPRINTARS',      'wet deposition flux (BC;CP)',                'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTBC_MP_RAIN_SPRINTARS', 'wet deposition flux (BC;MP(RAIN))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'WEFTBC_MP_SNOW_SPRINTARS', 'wet deposition flux (BC;MP(SNOW))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'DRFTBC_SPRINTARS',         'dry deposition flux (BC)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'GRFTBC_SPRINTARS',         'gravitational settling flux (BC)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TDFTBC_SPRINTARS',         'total deposition flux (BC)',                 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDC01_SPRINTARS', 'mixing ratio (carbon;internal mixed OM(OC) & BC)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDC02_SPRINTARS', 'mixing ratio (carbon;SOC from NVOC)',                      'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDC03_SPRINTARS', 'mixing ratio (carbon;external mixed BC)',                  'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONC01_SPRINTARS', 'mass concentration (carbon;internal mixed OM(OC) & BC)',   'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONC02_SPRINTARS', 'mass concentration (carbon;SOC from NVOC)',                'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONC03_SPRINTARS', 'mass concentration (carbon;external mixed BC)',            'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCC01_SPRINTARS', 'number concentration (carbon;internal mixed OM(OC) & BC)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCC02_SPRINTARS', 'number concentration (carbon;SOC from NVOC)',              '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMCC03_SPRINTARS', 'number concentration (carbon;external mixed BC)',          '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLC01_SPRINTARS', 'column loading (carbon;internal mixed OM(OC) & BC)',       'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLC02_SPRINTARS', 'column loading (carbon;SOC from NVOC)',                    'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLC03_SPRINTARS', 'column loading (carbon;external mixed BC)',                'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'HCPPB1_SPRINTARS', 'terpene',                                                  'ppbv',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'HCPPB2_SPRINTARS', 'isoprene',                                                 'ppbv',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QLDCT_SPRINTARS',  'mass mixing ratio (BC+OM(OC))',                            'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QTRCOC_SPRINTARS', 'mass mixing ratio (OM(OC))',                               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'QTRCBC_SPRINTARS', 'mass mixing ratio (BC)',                                   'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONCT_SPRINTARS',  'mass concentration (BC & OM(OC))',                         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONBT_SPRINTARS',  'mass concentration (BC)',                                  'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CONOT_SPRINTARS',  'mass concentration (OM(OC))',                              'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLCT_SPRINTARS',  'column loading (BC+OM(OC))',                               'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLBT_SPRINTARS',  'column loading (BC)',                                      'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLOT_SPRINTARS',  'column loading (OM(OC))',                                  'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUCA_SPRINTARS',   'AOT (BC+OM(OC);550nm;allsky)',         '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUCAA_SPRINTARS',  'AOT (BC+OM(OC);550nm;allsky;abs)',     '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUOC_SPRINTARS',   'AOT (OM(OC);550nm;allsky)',            '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUOCA_SPRINTARS',  'AOT (OM(OC);550nm;allsky;abs)',        '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUBC_SPRINTARS',   'AOT (BC;550nm;allsky)',                '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUBCA_SPRINTARS',  'AOT (BC;550nm;allsky;abs)',            '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTCA1_SPRINTARS',  'extinction  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSCA1_SPRINTARS',  'absorption  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKCA1_SPRINTARS',  'backscatter (BC+OM(OC);532nm;allsky)', '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTOC1_SPRINTARS',  'extinction  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSOC1_SPRINTARS',  'absorption  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKOC1_SPRINTARS',  'backscatter (OM(OC);532nm;allsky)',    '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTBC1_SPRINTARS',  'extinction  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSBC1_SPRINTARS',  'absorption  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKBC1_SPRINTARS',  'backscatter (BC;532nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'CONCTC_SPRINTARS',  'mass concentration (BC+OM(OC);clearsky)', 'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'CONBTC_SPRINTARS',  'mass concentration (BC);clearsky)',       'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLCTC_SPRINTARS',  'column loading (BC+OM(OC);clearsky)',     'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'COLBTC_SPRINTARS',  'column loading (BC;clearsky)',            'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUCAC_SPRINTARS',  'AOT (BC+OM(OC);550nm;clearsky)',          '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUCACA_SPRINTARS', 'AOT (BC+OM(OC);550nm;clearsky;abs)',      '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUOCC_SPRINTARS',  'AOT (OM(OC);550nm;clearsky)',             '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUOCCA_SPRINTARS', 'AOT (OM(OC);550nm;clearsky;abs)',         '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUBCC_SPRINTARS',  'AOT (BC;550nm;clearsky)',                 '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'TAUBCCA_SPRINTARS', 'AOT (BC;550nm;clearsky;abs)',             '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTCAC1_SPRINTARS', 'extinction  (BC+OM(OC);532nm;clearsky)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSCAC1_SPRINTARS', 'absorption  (BC+OM(OC);532nm;clearsky)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKCAC1_SPRINTARS', 'backscatter (BC+OM(OC);532nm;clearsky)',  '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTOCC1_SPRINTARS', 'extinction  (OM(OC);532nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSOCC1_SPRINTARS', 'absorption  (OM(OC);532nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKOCC1_SPRINTARS', 'backscatter (OM(OC);532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'EXTBCC1_SPRINTARS', 'extinction  (BC;532nm;clearsky)',         '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'ABSBCC1_SPRINTARS', 'absorption  (BC;532nm;clearsky)',         '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'BAKBCC1_SPRINTARS', 'backscatter (BC;532nm;clearsky)',         '1/m/sr', dim_type = 'ZXY' )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_setup
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> SPRINTARS Carbon finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_finalize
    implicit none

    deallocate(   BBBC_DATA )
    deallocate(   BBOC_DATA )
    deallocate(  ANTBC_DATA )
    deallocate(  ANTOC_DATA )
    deallocate(   ISOP_DATA )
    deallocate(   TERP_DATA )
    deallocate(  CHLOA_DATA )
    deallocate( DIFATT_DATA )
    deallocate( GX )
    deallocate( QTRCINP )
    deallocate( zero_3d )
    deallocate( zero_2d )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_finalize

  !-----------------------------------------------------------------------------
  !> SPRINTARS chemical reaction for SOA in aerocarb
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_FSOA &
       ( NUMHC,                       & ! [MODIFIED]
         NEWSOA,                      & ! [OUT]
         NUMOH, NUMO3, NUMNO3, CONOA, & ! [IN]
         GDT, dt_AE                   ) ! [IN]
    
    implicit none
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: NUMHC (KA,IA,JA,2) !< terpene/isoprene [molec/cm3]

    ! output
    !-------------------------
    real(RP), intent(out)   :: NEWSOA(KA,IA,JA)   !< produced SOA [ug/m3]
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: NUMOH (KA,IA,JA)   !< OH  [molec/cm3]
    real(RP), intent(in)    :: NUMO3 (KA,IA,JA)   !< O3  [molec/cm3]
    real(RP), intent(in)    :: NUMNO3(KA,IA,JA)   !< NO3 [molec/cm3]
    real(RP), intent(in)    :: CONOA (KA,IA,JA)   !< OC  [ug/m3]
    real(RP), intent(in)    :: GDT   (KA,IA,JA)   !< temperature T
    real(DP), intent(in)    :: dt_AE              !< time
    
    ! parameter
    !-------------------------
    real(RP), parameter :: KC1(NVOC/2) &       !< reaction rate
                           = (/  444.E+0_RP, & !  alpha-pinene + OH
                                -730.E+0_RP, & !  alpha-pinene + O3
                                -650.E+0_RP, & !  alpha-pinene + NO3
                                 410.E+0_RP, & !  isoprene + OH
                                -446.E+0_RP /) !  isoprene + NO3
    real(RP), parameter :: KC2(NVOC/2) &       !< reaction rate
                           = (/ 1.20E-11_RP, & !  alpha-pinene + OH
                                9.90E-16_RP, & !  alpha-pinene + O3
                                5.60E-11_RP, & !  alpha-pinene + NO3
                                2.45E-11_RP, & !  isoprene + OH 2.54E-11_RP
                                3.03E-12_RP /) !  isoprene + NO3
    real(RP), parameter :: COEF(NVOC)  &                  !< stoichiometric parameters
                           = (/ 2.58E-2_RP, 2.22E-1_RP, & !  alpha-pinene + OH
                                8.50E-2_RP, 6.94E-2_RP, & !  alpha-pinene + O3
                                0.68E+0_RP, 0.00E+0_RP, & !  alpha-pinene + NO3
                                1.21E-1_RP, 1.51E-2_RP, & !  isoprene + OH
                                4.66E-2_RP, 1.06E-1_RP /) !  isoprene + NO3
    real(RP), parameter :: KSOA0(NVOC) &                    !< KSOA @298K
                           = (/ 1.710E-1_RP, 0.040E-1_RP, & !  alpha-pienen + OH
                                0.880E-1_RP, 0.778E-1_RP, & !  alpha-pienen + O3
                                0.163E-1_RP, 0.000E+0_RP, & !  alpha-pienen + NO3
                                0.862E-2_RP, 1.650E+0_RP, & !  isoprene + OH
                                1.820E-1_RP, 0.046E+0_RP /) !  isoprene + NO3
    real(RP), parameter :: DSOAH(NVOC) &                  !< evaporation enthapy from Tsi [J/mol]
                           = (/ 79.0E+3_RP, 79.0E+3_RP, &
                                79.0E+3_RP, 79.0E+3_RP, &
                                79.0E+3_RP, 79.0E+3_RP, &
                                79.0E+3_RP, 79.0E+3_RP, &
                                79.0E+3_RP, 79.0E+3_RP /)
    real(RP), parameter :: MWSOA(NVOC) &                    !< molecular weight [g/mol]
                           = (/ 200.0E+0_RP, 200.0E+0_RP, &
                                200.0E+0_RP, 200.0E+0_RP, &
                                200.0E+0_RP, 200.0E+0_RP, &
                                130.0E+0_RP, 130.0E+0_RP, &
                                130.0E+0_RP, 130.0E+0_RP /)
    real(RP), parameter :: TS    =   6.0E+2_RP      !< time step
    real(RP), parameter :: GDT0  = 298.15E+0_RP     !< standard temperature
    real(RP), parameter :: FGAS  =   0.50E+0_RP     !< gas remaining ratio (50% of the previous step)
    real(RP), parameter :: FAUTO =   0.003E+0_RP    !< auto converstion ratio
    integer,  parameter :: NTERP = 6
    integer,  parameter :: NPMAX = 5
    
    ! internal work
    !-------------------------
    real(RP) :: CT   (KA,IA,JA,NVOC  )                   !<
    real(RP) :: CA   (KA,IA,JA)                          !< 
    ! real(RP) :: INREAC        (NVOC/2)                   !<
    real(RP) :: BETAR(KA,IA,JA,NVOC/2)                   !<
    real(RP) :: BETAX(KA,IA,JA,NVOC  )                   !<
    real(RP) :: KSOA (KA,IA,JA,NVOC  )                   !< 
    real(RP) :: M0   (KA,IA,JA)                          !<
    real(RP) :: NEWSOC                                   !<
    real(RP) :: XX                                       !< 
    real(RP) :: TSC                                      !< time step
    integer  :: NX                                       !< 
    integer  :: NSTEP                                    !< 
    real(RP) :: DSEL                                     !< 
    ! real(RP), pointer                   :: G0(:,:,:,:)   !< 
    ! real(RP), target, allocatable, save :: GX(:,:,:,:,:) !< 
    real(RP), save :: RST                                !< 
    real(RP), save :: KVAA                               !< sum of KSOA0
    integer,  save :: JCNT                               !< 
    logical,  save :: OFIRST = .true.                    !< 
    integer        :: ICNT
    ! others
    !-------------------------
    integer :: k, i, j, iq
    integer :: n
    
    !----------------------------------------------------

    ! set parameter
    !===========================
    if ( OFIRST ) then
       OFIRST = .false.
       ! allocate( GX(KA,IA,JA,NVOC,2) )
       ! GX(:,:,:,:,:) = 0.0_RP
       JCNT = 0
       KVAA = 0.0_RP
       do iq = 1, NVOC
          KVAA = KVAA + KSOA0(iq)
       enddo
       KVAA = max( KVAA, 1.0E-8_RP )
    endif
    
    JCNT = JCNT + 1
    DSEL = real( mod(JCNT,2), kind = RP )
    !===========================
    ! end set parameter

    
    ! offline SOA
    !===========================
    if ( real(dt_AE,kind=RP) < TS ) then
       TSC = real(dt_AE,kind=RP)
       NSTEP = 1
    else
       RST = real(dt_AE,kind=RP) / TS
       NSTEP = int(RST)
       if ( RST - real(NSTEP,kind=RP) == 0.0_RP ) then
          TSC = TS
       else
          NSTEP = NSTEP + 1
          TSC = real(dt_AE,kind=RP) / real(NSTEP,kind=RP)
       endif
    endif

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! INREAC(1) = KC1(1) / GDT(k,i,j)
          ! INREAC(2) = KC1(2) / GDT(k,i,j)
          ! INREAC(3) = KC1(3) / GDT(k,i,j)
          ! INREAC(4) = KC1(4) / GDT(k,i,j)
          ! INREAC(5) = KC1(5) / GDT(k,i,j)
          ! BETAR(k,i,j,1) = KC2(1) * exp(INREAC(1)) * NUMOH (k,i,j)
          ! BETAR(k,i,j,2) = KC2(2) * exp(INREAC(2)) * NUMO3 (k,i,j)
          ! BETAR(k,i,j,3) = KC2(3) * exp(INREAC(3)) * NUMNO3(k,i,j)
          ! BETAR(k,i,j,4) = KC2(4) * exp(INREAC(4)) * NUMOH (k,i,j)
          ! BETAR(k,i,j,5) = KC2(5) * exp(INREAC(5)) * NUMNO3(k,i,j)
          BETAR(k,i,j,1) = KC2(1) * exp( KC1(1) / GDT(k,i,j) ) * NUMOH (k,i,j)
          BETAR(k,i,j,2) = KC2(2) * exp( KC1(2) / GDT(k,i,j) ) * NUMO3 (k,i,j)
          BETAR(k,i,j,3) = KC2(3) * exp( KC1(3) / GDT(k,i,j) ) * NUMNO3(k,i,j)
          BETAR(k,i,j,4) = KC2(4) * exp( KC1(4) / GDT(k,i,j) ) * NUMOH (k,i,j)
          BETAR(k,i,j,5) = KC2(5) * exp( KC1(5) / GDT(k,i,j) ) * NUMNO3(k,i,j)
       enddo
    enddo
    enddo

    do iq = 1, NVOC
       NX = max( int(iq/2), int((iq+1)/2) )
       BETAX(:,:,:,iq) = BETAR(:,:,:,NX)
    enddo

    ! Backward Euler
    CT(:,:,:,:) = 0.0_RP
    do n = 1, NSTEP             
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             
             NUMHC(k,i,j,1) = NUMHC(k,i,j,1) / ( 1.0_RP + sum(BETAR(k,i,j,1:3))*TSC )
             NUMHC(k,i,j,2) = NUMHC(k,i,j,2) / ( 1.0_RP + sum(BETAR(k,i,j,4:5))*TSC )
             do iq = 1, NTERP
                CT(k,i,j,iq) = CT(k,i,j,iq) + BETAX(k,i,j,iq) * COEF(iq) * NUMHC(k,i,j,1) * TSC
             enddo
             do iq = NTERP+1,  NVOC
                CT(k,i,j,iq) = CT(k,i,j,iq) + BETAX(k,i,j,iq) * COEF(iq) * NUMHC(k,i,j,2) * TSC
             enddo
             
          enddo
       enddo
       enddo
    enddo
    !===========================
    ! end offline SOA
    
    CA(:,:,:) = 0.0_RP
    do iq = 1, NVOC
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             KSOA(k,i,j,iq) = KSOA0(iq) * GDT(k,i,j) / GDT0 * exp( DSOAH(iq)/RCON * ( 1.0_RP/GDT(k,i,j) - 1.0_RP/GDT0 ) )
             CT  (k,i,j,iq) = CT(k,i,j,iq) * MWSOA(iq) / AVOG * 1.0E+12_RP & ! [ug/m3]
                            + FGAS * (           DSEL  * GX(k,i,j,iq,1) &
                                       + (1.0_RP-DSEL) * GX(k,i,j,iq,2) )
             CA  (k,i,j)    = CA(k,i,j) + KSOA0(iq) * CT(k,i,j,iq) / KVAA
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          M0(k,i,j) = CONOA(k,i,j) ! initial OA
          M0(k,i,j) = max( M0(k,i,j), FAUTO*CA(k,i,j) )
       enddo
    enddo
    enddo

    ! G0 => GX(:,:,:,:,1+iand(JCNT-1,1))
    ICNT = 1 + iand(JCNT-1,1)
    
    do n = 1, NPMAX

       do j = JS, JE
       do i = IS, IE
          do k = KS, KE

             NEWSOC = 0.0_RP
             do iq = 1, NVOC
                XX = CT(k,i,j,iq) *            KSOA(k,i,j,iq) * M0(k,i,j) &
                                  / ( 1.0_RP + KSOA(k,i,j,iq) * M0(k,i,j) )
                ! NEWSOC = NEWSOC + XX
                NEWSOC = NEWSOC + max( XX, 0.0_RP )
                ! G0(k,i,j,iq) = CT(k,i,j,iq) - XX
                ! GX(k,i,j,iq,ICNT) = CT(k,i,j,iq) - XX
                GX(k,i,j,iq,ICNT) = max( CT(k,i,j,iq) - XX, 0.0_RP )
             enddo
             NEWSOA(k,i,j) = NEWSOC
             M0(k,i,j) = CONOA(k,i,j) + NEWSOA(k,i,j)

          enddo
       enddo
       enddo
       
    enddo
 
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA_FSOA
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Carbon
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA &
       ( QTRC,           & ! [MODIFIED]
         OUTQLD,         & ! [OUT]
         NUMCON,         & ! [OUT]
         CCONT,          & ! [OUT]
         WATROC,         & ! [OUT]
         TAUCA,          & ! [OUT]
         TAUCAA,         & ! [OUT]
         CEXTCA,         & ! [OUT]
         CABSCA,         & ! [OUT]
         CBAKCA,         & ! [OUT]
         TAUCAD,         & ! [OUT]
         CEXTCD,         & ! [OUT]
         CDVE,           & ! [IN]
         CCOVMR,         & ! [IN]
         CCMAX,          & ! [IN]
         WS10,           & ! [IN]
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
         KBB,            & ! [IN]
         DPKBB,          & ! [IN]
         GDPM,           & ! [IN]
         DELP,           & ! [IN]
         GDZM,           & ! [IN]
         DELZ,           & ! [IN]
         ONEWQ,          & ! [IN]
         RNWR,           & ! [IN]
         RNWR_CP,        & ! [IN] OK
         RNWR_MP_RAIN,   & ! [IN] OK
         RNWR_MP_SNOW,   & ! [IN] OK
         INDEX_PFT,      & ! [IN]
         FRAC_PFT,       & ! [IN]
         FACT_LAND,      & ! [IN]
         FACT_OCEAN,     & ! [IN]
         GOICE,          & ! [IN]
         AIRFEL,         & ! [IN]
         AIRFELOC,       & ! [IN]
         EMITSA,         & ! [IN]
         RFLXSD,         & ! [IN]
         GDZ,            & ! [IN]
         NUMOH,          & ! [IN]
         PPMO3,          & ! [IN]
         NO3G,           & ! [IN]
         GDT,            & ! [IN]
         ACTIBC,         & ! [IN]
         ACTICA,         & ! [IN]
         ACTIOC,         & ! [IN]
         SW_CCN,         & ! [IN]
         TIME_NOWDATE,   & ! [IN]
         dt_AE,          & ! [IN]
         QA_AE_CA        ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE                          !< time step
    integer,  intent(in)    :: QA_AE_CA                       !< total number of carbon tracers

    real(RP), intent(in)    :: CDVE               (IA,JA)     !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
    real(RP), intent(in)    :: CCOVMR             (IA,JA)     !< cloud cover (max-ran)
    real(RP), intent(in)    :: CCMAX                          !< maximum cloud fraction for all/clear sky
    real(RP), intent(in)    :: WS10               (IA,JA)     !< wind speed at 10m
    real(RP), intent(in)    :: RH            (KA,  IA,JA)     !< relative humidity
    real(RP), intent(in)    :: DENS          (KA,  IA,JA)     !< air density
    real(RP), intent(in)    :: DENS1         (KA,  IA,JA)     !< air density   on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA)     !< precipitation flux (liq    ) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA)     !< precipitation flux (    ice) : MP on 1 step before
    real(RP), intent(in)    :: FRAIN         (0:KA,IA,JA)     !< precipitation flux
    real(RP), intent(in)    :: FRAIN_CP      (0:KA,IA,JA)     !< precipitation flux (liq+ice) : CP
    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA)     !< precipitation flux (liq    ) : MP
    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA)     !< precipitation flux (    ice) : MP
    real(RP), intent(in)    :: VTR_MP_RAIN   (KA,  IA,JA)     !< terminal velocity of rain (microphysic scheme) [m/s]
    real(RP), intent(in)    :: VTR_MP_SNOW   (KA,  IA,JA)     !< terminal velocity of snow (microphysic scheme) [m/s]
    real(RP), intent(in)    :: RADR_MP_RAIN  (KA,  IA,JA)     !< effective radius of rain (microphysic scheme) [m]
    real(RP), intent(in)    :: RADR_MP_SNOW  (KA,  IA,JA)     !< effective radius of snow (microphysic scheme) [m]
    real(RP), intent(in)    :: RAINN_CP      (KA,  IA,JA)     !< rain number concentration (CP)      [1/m3]
    real(RP), intent(in)    :: RAINN_MP_RAIN (KA,  IA,JA)     !< rain number concentration (MP;RAIN) [1/m3]
    real(RP), intent(in)    :: RAINN_MP_SNOW (KA,  IA,JA)     !< rain number concentration (MP;SNOW) [1/m3]
    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA)     !< cloud fraction (total)
    real(RP), intent(in)    :: TCLDF2        (KA,  IA,JA)     !< cloud fraction (total) for CCOVMR
    integer,  intent(in)    :: KUP                (IA,JA)     !< uplift level
    integer,  intent(in)    :: KUPMAX                         !< maximum uplift level
    real(RP), intent(in)    :: DPKUP              (IA,JA)     !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
    integer,  intent(in)    :: KBB                (IA,JA)     !< uplift level for biomass burning 
    real(RP), intent(in)    :: DPKBB              (IA,JA)     !< GDPM(KS-1,i,j) - GDPM(KBB(i,j),i,j)
    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA)     !< pressure at half levels
    real(RP), intent(in)    :: DELP          (KA,  IA,JA)     !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA)     !< height at half levels
    real(RP), intent(in)    :: DELZ          (KA,  IA,JA)     !< GDZM(k,i,j)-GDZM(k-1,i,j)
    logical,  intent(in)    :: ONEWQ         (KA-1,IA,JA)     !< use for instability
    real(RP), intent(in)    :: RNWR          (KA,  IA,JA)     !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA)     !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA)     !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA)     !< ratio of rain to cloud
    integer,  intent(in)    :: INDEX_PFT          (IA,JA,2)   !< index of PFT for each mosaic
    real(RP), intent(in)    :: FRAC_PFT           (IA,JA,2)   !< fraction of PFT for each mosaic
    real(RP), intent(in)    :: FACT_LAND          (IA,JA)     !< land factor
    real(RP), intent(in)    :: FACT_OCEAN         (IA,JA)     !< ocean factor
    real(RP), intent(in)    :: GOICE              (IA,JA)     !< ocean ice fraction
    real(RP), intent(in)    :: AIRFEL        (KA,  IA,JA)     !< fuel consumption      by aircrafts
    real(RP), intent(in)    :: AIRFELOC      (KA,  IA,JA)     !< fuel consumption (OC) by aircrafts
    real(RP), intent(in)    :: EMITSA             (IA,JA,NSA) !< Sea salt emission flux
    real(RP), intent(in)    :: RFLXSD             (IA,JA)     !< downward short wave radiation [W/m2]
    real(RP), intent(in)    :: GDZ           (KA,  IA,JA)     !< altitude
    real(RP), intent(in)    :: NUMOH         (KA,  IA,JA)     !< OH concentration [molecule/cm3]
    real(RP), intent(in)    :: PPMO3         (KA,  IA,JA)     !< O3       mixing ratio: ppmv
    real(RP), intent(in)    :: NO3G          (KA,  IA,JA)     !< NO3(gas) mixing ratio: pptv
    real(RP), intent(in)    :: GDT           (KA,  IA,JA)     !< temperature T
    real(RP), intent(in)    :: ACTIBC        (KA,  IA,JA)     !< ratio of activated aerosols (BC)         [0-1]
    real(RP), intent(in)    :: ACTICA        (KA,  IA,JA)     !< ratio of activated aerosols (carbon)     [0-1]
    real(RP), intent(in)    :: ACTIOC        (KA,  IA,JA)     !< ratio of activated aerosols (natural OC) [0-1]
    logical,  intent(in)    :: SW_CCN                         !< switch for CCN scheme (true->use ACTI)
    integer,  intent(in)    :: TIME_NOWDATE(6)                !< time
 
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_CA) !< ratio of carbon tracer mass to total mass
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: OUTQLD(KA,IA,JA,NRBCOM) !< QTRC for radiation
    real(RP), intent(out)   :: NUMCON(KA,IA,JA,NCA1)   !< number concentration
    real(RP), intent(out)   :: CCONT (KA,IA,JA)        !< mass concentraion (total)
    real(RP), intent(out)   :: WATROC(KA,IA,JA)        !< aerosol water
    real(RP), intent(out)   :: TAUCA    (IA,JA,IWA,2)  !< AOT
    real(RP), intent(out)   :: TAUCAA   (IA,JA,IWA,2)  !< AOT (absorption)
    real(RP), intent(out)   :: CEXTCA(KA,IA,JA,IWA,2)  !< extinction  coefficient
    real(RP), intent(out)   :: CABSCA(KA,IA,JA,IWA,2)  !< absorption  coefficient
    real(RP), intent(out)   :: CBAKCA(KA,IA,JA,IWA,2)  !< backscatter coefficient
    real(RP), intent(out)   :: TAUCAD   (IA,JA)        !< AOT (dry)
    real(RP), intent(out)   :: CEXTCD(KA,IA,JA)        !< extinction coefficient (carbon;dry)
    
    ! internal work
    !-------------------------
    !set parameter
    real(RP) :: DRYMCB        (NCA1)  !< dry particle mass (BC) [kg]
    real(RP) :: VTERBC                !< BC terminal velocity   [m/s]
    real(RP) :: VTER (KA,IA,JA,NCA1)  !< terminal velocity      [m/s]
    real(RP) :: RADC (KA,IA,JA)       !< carbon particle radius [m]
    real(RP) :: BRH1, BRH2, BRH3      !< 
    real(RP) :: CEXTP(KA,IA,JA,4,IWA) !< extinction  cross section
    real(RP) :: CABSP(KA,IA,JA,4,IWA) !< absorption  cross section
    real(RP) :: CEXTA(KA,IA,JA,4,IWA) !< extinction  cross section
    real(RP) :: CABSA(KA,IA,JA,4,IWA) !< absorption  cross section
    real(RP) :: CBAKA(KA,IA,JA,4,IWA) !< backscatter cross section
 
    !emission
    real(RP) :: BBBC     (IA,JA)     !< BC emission (biomass burning)
    real(RP) :: BBOC     (IA,JA)     !< OC emission (biomass burning)
    real(RP) :: ANTBC    (IA,JA)     !< BC emission (fossil fuel)
    real(RP) :: ANTOC    (IA,JA)     !< OC emission (fossil fuel)
    real(RP) :: TERP     (IA,JA)     !< emission (terpene)
    real(RP) :: ISOP     (IA,JA)     !< emission (isoprene)
    real(RP) :: CHLOA    (IA,JA)     !< CHL-A concentration
    real(RP) :: DIFATT   (IA,JA)     !< diffuse attenuation coef. at 490nm
    real(RP) :: DIFAT    (IA,JA)     !< diffuse attenuation coef. at 490nm
    real(RP) :: EMITBC   (IA,JA)     !< BC emission flux [kg/m2/s]
    real(RP) :: EMITOC   (IA,JA)     !< OM(OC) emission flux [kg/m2/s]
    real(RP) :: QEMIT    (IA,JA,10)  !< emission mixing ratio [kg/kg]
    real(RP) :: OMSSA    (IA,JA,NSA) !< organic mass fraction of sea spray
    real(RP) :: EMITSS   (IA,JA,NSA) !< OM(SSA) emission flux
    real(RP) :: EMISSA   (IA,JA)     !< OM(SSA) emission flux
    real(RP) :: RMAG     (IA,JA)     !< log( downward SW flux )
    real(RP) :: HMAX     (IA,JA)     !< water depth where downward SW flux attenuates until 2.5W/m2
    real(RP) :: DELTSP               !< delta of integration
    real(RP) :: PISO                 !< isoprene production rate
    real(RP) :: EMIISO   (IA,JA,2)   !< emission of isopren                        [kg/m2/s]
    real(RP) :: EMIOBS   (IA,JA,2)   !< emission of isopren (compered observation) [molecules/cm2/s]
    real(RP) :: EMOMPR   (IA,JA,2)   !< emission of oceanic OM precursor
    real(RP) :: AFFBC (KA,IA,JA)     !< BC emission factor for aircrafts

    !chemical reaction for SOA
    real(RP) :: NUMHC (KA,IA,JA,NSOA)      !< 1 -> terpene, 2 -> isoprene [molecules/cm3]
    real(RP) :: NUMO3 (KA,IA,JA)           !< O3                          [molecules/cm3]
    real(RP) :: NUMNO3(KA,IA,JA)           !< NO3                         [molecules/cm3]
    real(RP) :: CONOA (KA,IA,JA)           !< OC [ug/m3]
    real(RP) :: NEWSOA(KA,IA,JA)           !< produced SOA [ug/m3]
    real(RP) :: CHESOA   (IA,JA)           !< produced SOA [kg/m2/s]
    real(RP) :: QTRCC (KA,IA,JA,NCA1+NSOA) !< carbon mixing ratio (1->int.OC&BC,2->SOC,3->ext.BC,4->terpene,5->isoprene) [kg/kg]
    real(RP) :: RCABC (KA,IA,JA)           !< BC/(OM(OC)+BC) (internal mixing)

    !incloud
    real(RP) :: QTRCTM1(KA,IA,JA,NCA1) !< outside cloud water [kg/kg]
    real(RP) :: QTRCTM2(KA,IA,JA,NCA1) !< inside  cloud water [kg/kg]

    !subcloud scavenging (wash out)
    ! real(RP) :: RAINN         (KA,IA,JA)      !< rain number concentration [1/m3]
    ! real(RP) :: RAINN_CP      (KA,IA,JA)
    ! real(RP) :: RAINN_MP_RAIN (KA,IA,JA)
    ! real(RP) :: RAINN_MP_SNOW (KA,IA,JA)
    real(RP) :: CARBN                         !< number density of carbon tracer [1/m3]
    real(RP) :: CBRN          (KA,IA,JA,NCA1) !< number density of carbon tracer entering into raindrops [1/s/m3] <=> collision frequency between carbon tracer and rain
    real(RP) :: CBRN_CP       (KA,IA,JA,NCA1)
    real(RP) :: CBRN_MP_RAIN  (KA,IA,JA,NCA1)
    real(RP) :: CBRN_MP_SNOW  (KA,IA,JA,NCA1)
    real(RP) :: QWTDEP        (KA,IA,JA,NCA1) !< mixing ratio by wet deposition [kg/kg]
    real(RP) :: QWTDEP_CP     (KA,IA,JA,NCA1)
    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,NCA1)
    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,NCA1)
    real(RP) :: QTRCTMP       (KA,IA,JA,NCA1) !< carbon tracers that do not cllide with rain particles [kg/kg]
    
    !incloud scavenging (rain out)
    real(RP) :: CLTORN        (KA,IA,JA,NCA1) !< mixing ratio moving cloud to rain particles [kg/kg]
    real(RP) :: CLTORN_CP     (KA,IA,JA,NCA1)
    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,NCA1)
    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,NCA1)
    
    !re-emission from rain
    real(RP) :: WETF           (IA,JA,NCA1) !< wet deposition flux
    real(RP) :: WETF_CP        (IA,JA,NCA1)
    real(RP) :: WETF_MP_RAIN   (IA,JA,NCA1)
    real(RP) :: WETF_MP_SNOW   (IA,JA,NCA1)
    real(RP) :: QTRCINR     (KA,IA,JA,4)    !< mixing ratio in rain(MP)
    real(RP) :: QTRCINS     (KA,IA,JA,4)    !< mixing ratio in snow(MP)
    real(RP) :: QTRCINRC    (KA,IA,JA,NCA1) !< mixing ratio in rain(MP)
    real(RP) :: QTRCINSC    (KA,IA,JA,NCA1) !< mixing ratio in snow(MP)

    !dry deposition
    real(RP) :: CDVEC   (IA,JA)           !< bulk coefficient
    real(RP) :: QMTX (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
    real(RP) :: DRYF    (IA,JA,NCA1+NSOA) !< dry deposition flux

    !gravitational settling
    real(RP) :: GRAVF(IA,JA,NCA1)      !< gravitational settling flux

    !distribution to BC/OM(OC) = 0.30, 0.15, 0.00 for radiation
    real(RP) :: QLOAD (KA,IA,JA,NCA1) !< QTRCC for radiation
    real(RP) :: QBC   (KA,IA,JA)      !<
    real(RP) :: QOC   (KA,IA,JA)      !<
    real(RP) :: QC                    !< 
    real(RP) :: F                     !< 
    real(RP) :: FF                    !< 
    real(RP) :: QTRCBC(KA,IA,JA)      !< mass mixing ratio (BC)     [kg/kg]
    real(RP) :: QTRCOC(KA,IA,JA)      !< mass mixing ratio (OM(OC)) [kg/kg]

    !aerosol
    real(RP) :: QLOADT(KA,IA,JA)      !< mixing ratio (total)                       [kg/kg]
    real(RP) :: CCON  (KA,IA,JA,NCA1) !< mass concentration                         [kg/m3]
    real(RP) :: COLMS    (IA,JA,NCA1) !< column mass loading                        [kg/m2]
    real(RP) :: COLMST   (IA,JA)      !< column mass loading (total)                [kg/m2]
    real(RP) :: BCCONT(KA,IA,JA)      !< BC mass concentration (total)              [kg/m3]
    real(RP) :: OCCONT(KA,IA,JA)      !< OM(OC) mass concentration (total)          [kg/m3]
    real(RP) :: COLBCT   (IA,JA)      !< BC column mass load (total)                [kg/m2]
    real(RP) :: COLOCT   (IA,JA)      !< OM(OC) column mass load (total)            [kg/m2]
    real(RP) :: COL1BC   (IA,JA)      !< BC column mass load (internal mixing)      [kg/m2]
    real(RP) :: HCPPB (KA,IA,JA,NSOA) !< volume concentration of terpene & isoprene [ppbv]

    !flux
    real(RP) :: WETFT         (IA,JA)   !< wet deposition flux (total)           [kg/m2/s]
    real(RP) :: WETFT_CP      (IA,JA)   !< wet deposition flux (total;CP)        [kg/m2/s]
    real(RP) :: WETFT_MP_RAIN (IA,JA)   !< wet deposition flux (total;MP(RAIN))  [kg/m2/s]
    real(RP) :: WETFT_MP_SNOW (IA,JA)   !< wet deposition flux (total;MP(SNOW))  [kg/m2/s]
    real(RP) :: WETFBC        (IA,JA)   !< wet deposition flux (BC)              [kg/m2/s]
    real(RP) :: WETFBC_CP     (IA,JA)   !< wet deposition flux (BC;CP)           [kg/m2/s]
    real(RP) :: WETFBC_MP_RAIN(IA,JA)   !< wet deposition flux (BC;MP(RAIN))     [kg/m2/s]
    real(RP) :: WETFBC_MP_SNOW(IA,JA)   !< wet deposition flux (BC;MP(SNOW))     [kg/m2/s]
    real(RP) :: WETFOC        (IA,JA)   !< wet deposition flux (OM(OC))          [kg/m2/s]
    real(RP) :: WETFOC_CP     (IA,JA)   !< wet deposition flux (OM(OC);CP)       [kg/m2/s]
    real(RP) :: WETFOC_MP_RAIN(IA,JA)   !< wet deposition flux (OM(OC);MP(RAIN)) [kg/m2/s]
    real(RP) :: WETFOC_MP_SNOW(IA,JA)   !< wet deposition flux (OM(OC);MP(SNOW)) [kg/m2/s]
    real(RP) :: DRYFT         (IA,JA)   !< dry deposition flux (total)           [kg/m2/s]
    real(RP) :: DRYFBC        (IA,JA)   !< dry deposition flux (BC)              [kg/m2/s]
    real(RP) :: DRYFOC        (IA,JA)   !< dry deposition flux (OM(OC))          [kg/m2/s]
    real(RP) :: GRAVFT        (IA,JA)   !< gravtational settling flux (total)    [kg/m2/s]
    real(RP) :: GRVFBC        (IA,JA)   !< gravtational settling flux (BC)       [kg/m2/s]
    real(RP) :: GRVFOC        (IA,JA)   !< gravtational settling flux (OM(OC))   [kg/m2/s]
    real(RP) :: DEPFT         (IA,JA)   !< total deposition flux                 [kg/m2/s]
    real(RP) :: DEPFBC        (IA,JA)   !< total deposition flux (BC)            [kg/m2/s]
    real(RP) :: DEPFOC        (IA,JA)   !< total deposition flux (OM(OC))        [kg/m2/s]
    
    !optical parameter
    real(RP) :: TUSE               !< use for AOT
    real(RP) :: TAUOC    (IA,JA,2) !< optical thickness (OM(OC))
    real(RP) :: TAUOCA   (IA,JA,2) !< optical thickness (OM(OC);abs)
    real(RP) :: TAUBC    (IA,JA,2) !< optical thickness (BC)
    real(RP) :: TAUBCA   (IA,JA,2) !< optical thickness (BC;abs)
    real(RP) :: CEXTBC(KA,IA,JA,2) !< extinction  cross section (BC)
    real(RP) :: CABSBC(KA,IA,JA,2) !< absorption  cross section (BC)
    real(RP) :: CBAKBC(KA,IA,JA,2) !< backscatter cross section (BC)
    real(RP) :: CEXTOC(KA,IA,JA,2) !< extinction  cross section (OM(OC))
    real(RP) :: CABSOC(KA,IA,JA,2) !< absorption  cross section (OM(OC))
    real(RP) :: CBAKOC(KA,IA,JA,2) !< backscatter cross section (OM(OC))
   
    !clear sky diagnoses
    real(RP) :: BCCNTC(KA,IA,JA)      !< mass concentraion (total,clearsky)
    real(RP) :: CCONTC(KA,IA,JA)      !< BC mass concentraion (total,clearsky)
    real(RP) :: NUMCNC(KA,IA,JA,NCA1) !< number concentration (clearsky)
    real(RP) :: COLMTC   (IA,JA)      !< column mass loading (total,clearsky)
    real(RP) :: COLBTC   (IA,JA)      !< BC column mass load (total,clearsky)
    
    ! parameter
    !-------------------------
    real(RP), parameter :: &
         MTONCA(NCA1)      &    !< mcon -> ncon
         !    OC & int.BC   SOC           ext. BC  
         = (/ 2.802E+16_RP, 2.802E+16_RP, 7.271E+18_RP /)
    real(RP), parameter :: &
         SIGMA (NCA1)      &    !< GSD of log-normal distribution
         !    OC & int.BC  SOC          ext. BC  
         = (/ 1.562E+0_RP, 1.562E+0_RP, 2.000E+0_RP /)
    real(RP), parameter :: &
         CEXTBP(IWA)       &    !< BC extinction  cross section (550/ 440/870nm)
         = (/ 3.715E+3_RP, 5.118E+3_RP, 1.976E+3_RP /)
    real(RP), parameter :: &
         CABSBP(IWA)       &    !< BC absorption  cross section (550/ 440/870nm)
         = (/ 3.097E+3_RP, 4.053E+3_RP, 1.807E+3_RP /)
    real(RP), parameter :: &
         CEXTBA(IWA)       &    !< BC extinction  cross section (532/1064/355nm)
         = (/ 3.923E+3_RP, 1.561E+3_RP, 6.914E+3_RP /)
    real(RP), parameter :: &
         CABSBA(IWA)       &    !< BC absorption  cross section (532/1064/355nm)
         = (/ 3.250E+3_RP, 1.471E+3_RP, 5.211E+3_RP /)
    real(RP), parameter :: &
         CBAKBA(IWA)       &    !< BC backscatter cross section (532/1064/355nm)
         = (/ 3.650E+1_RP, 6.923E+0_RP, 7.184E+1_RP /)
    real(RP), parameter :: &
         CEXTCP(3,NRH,IWA) &    !< extinction  cross section (550/ 440/870nm)
         !             BC/OC=0.30   BC/OC=0.15   BC/OC=0.00
         = reshape( (/ 4.385E+3_RP, 4.249E+3_RP, 4.035E+3_RP, & ! RH = 0.00, IWA =  550nm
                       5.054E+3_RP, 4.932E+3_RP, 4.742E+3_RP, & ! RH = 0.50, IWA =  550nm
                       5.242E+3_RP, 5.127E+3_RP, 4.940E+3_RP, & ! RH = 0.70, IWA =  550nm
                       9.583E+3_RP, 9.625E+3_RP, 9.673E+3_RP, & ! RH = 0.80, IWA =  550nm
                       1.449E+4_RP, 1.469E+4_RP, 1.496E+4_RP, & ! RH = 0.90, IWA =  550nm
                       2.147E+4_RP, 2.190E+4_RP, 2.248E+4_RP, & ! RH = 0.95, IWA =  550nm
                       5.056E+4_RP, 5.178E+4_RP, 5.342E+4_RP, & ! RH = 0.98, IWA =  550nm
                       6.819E+4_RP, 6.985E+4_RP, 7.205E+4_RP, & ! RH = 0.99, IWA =  550nm
                       5.530E+3_RP, 5.528E+3_RP, 5.511E+3_RP, & ! RH = 0.00, IWA =  440nm
                       6.394E+3_RP, 6.425E+3_RP, 6.444E+3_RP, & ! RH = 0.50, IWA =  440nm
                       6.636E+3_RP, 6.677E+3_RP, 6.700E+3_RP, & ! RH = 0.70, IWA =  440nm
                       1.199E+4_RP, 1.224E+4_RP, 1.257E+4_RP, & ! RH = 0.80, IWA =  440nm
                       1.769E+4_RP, 1.811E+4_RP, 1.867E+4_RP, & ! RH = 0.90, IWA =  440nm
                       2.528E+4_RP, 2.593E+4_RP, 2.681E+4_RP, & ! RH = 0.95, IWA =  440nm
                       5.350E+4_RP, 5.485E+4_RP, 5.666E+4_RP, & ! RH = 0.98, IWA =  440nm
                       6.914E+4_RP, 7.086E+4_RP, 7.313E+4_RP, & ! RH = 0.99, IWA =  440nm
                       2.142E+3_RP, 1.919E+3_RP, 1.595E+3_RP, & ! RH = 0.00, IWA =  870nm
                       2.466E+3_RP, 2.227E+3_RP, 1.911E+3_RP, & ! RH = 0.50, IWA =  870nm
                       2.556E+3_RP, 2.320E+3_RP, 2.003E+3_RP, & ! RH = 0.70, IWA =  870nm
                       4.758E+3_RP, 4.577E+3_RP, 4.337E+3_RP, & ! RH = 0.80, IWA =  870nm
                       7.488E+3_RP, 7.397E+3_RP, 7.264E+3_RP, & ! RH = 0.90, IWA =  870nm
                       1.183E+4_RP, 1.186E+4_RP, 1.190E+4_RP, & ! RH = 0.95, IWA =  870nm
                       3.416E+4_RP, 3.481E+4_RP, 3.569E+4_RP, & ! RH = 0.98, IWA =  870nm
                       5.063E+4_RP, 5.172E+4_RP, 5.316E+4_RP /), & ! RH = 0.99, IWA =  870nm
                    (/ 3, NRH, IWA /) )
    real(RP), parameter :: &
         CABSCP(3,NRH,IWA) &    !< absorption  cross section (550/ 440/870nm)
         !             BC/OC=0.30   BC/OC=0.15   BC/OC=0.00
         = reshape( (/ 1.313E+3_RP, 8.771E+2_RP, 1.238E+2_RP, & ! RH = 0.00, IWA =  550nm
                       1.337E+3_RP, 8.879E+2_RP, 1.224E+2_RP, & ! RH = 0.50, IWA =  550nm
                       1.344E+3_RP, 8.888E+2_RP, 1.221E+2_RP, & ! RH = 0.70, IWA =  550nm
                       1.455E+3_RP, 9.264E+2_RP, 1.205E+2_RP, & ! RH = 0.80, IWA =  550nm
                       1.521E+3_RP, 9.518E+2_RP, 1.209E+2_RP, & ! RH = 0.90, IWA =  550nm
                       1.569E+3_RP, 9.703E+2_RP, 1.218E+2_RP, & ! RH = 0.95, IWA =  550nm
                       1.668E+3_RP, 1.014E+3_RP, 1.244E+2_RP, & ! RH = 0.98, IWA =  550nm
                       1.693E+3_RP, 1.027E+3_RP, 1.255E+2_RP, & ! RH = 0.99, IWA =  550nm
                       1.631E+3_RP, 1.119E+3_RP, 1.435E+2_RP, & ! RH = 0.00, IWA =  440nm
                       1.666E+3_RP, 1.135E+3_RP, 1.406E+2_RP, & ! RH = 0.50, IWA =  440nm
                       1.677E+3_RP, 1.136E+3_RP, 1.400E+2_RP, & ! RH = 0.70, IWA =  440nm
                       1.840E+3_RP, 1.187E+3_RP, 1.358E+2_RP, & ! RH = 0.80, IWA =  440nm
                       1.936E+3_RP, 1.220E+3_RP, 1.354E+2_RP, & ! RH = 0.90, IWA =  440nm
                       2.002E+3_RP, 1.243E+3_RP, 1.357E+2_RP, & ! RH = 0.95, IWA =  440nm
                       2.138E+3_RP, 1.297E+3_RP, 1.375E+2_RP, & ! RH = 0.98, IWA =  440nm
                       2.169E+3_RP, 1.314E+3_RP, 1.383E+2_RP, & ! RH = 0.99, IWA =  440nm
                       7.988E+2_RP, 5.277E+2_RP, 1.287E+2_RP, & ! RH = 0.00, IWA =  870nm
                       8.177E+2_RP, 5.421E+2_RP, 1.305E+2_RP, & ! RH = 0.50, IWA =  870nm
                       8.228E+2_RP, 5.442E+2_RP, 1.309E+2_RP, & ! RH = 0.70, IWA =  870nm
                       8.968E+2_RP, 5.822E+2_RP, 1.368E+2_RP, & ! RH = 0.80, IWA =  870nm
                       9.404E+2_RP, 6.036E+2_RP, 1.406E+2_RP, & ! RH = 0.90, IWA =  870nm
                       9.743E+2_RP, 6.229E+2_RP, 1.441E+2_RP, & ! RH = 0.95, IWA =  870nm
                       1.046E+3_RP, 6.621E+2_RP, 1.515E+2_RP, & ! RH = 0.98, IWA =  870nm
                       1.067E+3_RP, 6.754E+2_RP, 1.540E+2_RP /), & ! RH = 0.99, IWA =  870nm
                    (/ 3, NRH, IWA /) )
    real(RP), parameter :: &
         CEXTCC(3,NRH,IWA) &    !< extinction  cross section (532/1064/355nm)
         !             BC/OC=0.30   BC/OC=0.15   BC/OC=0.00
         = reshape( (/ 4.563E+3_RP, 4.443E+3_RP, 4.250E+3_RP, & ! RH = 0.00, IWA =  532nm
                       5.261E+3_RP, 5.157E+3_RP, 4.992E+3_RP, & ! RH = 0.50, IWA =  532nm
                       5.457E+3_RP, 5.361E+3_RP, 5.199E+3_RP, & ! RH = 0.70, IWA =  532nm
                       9.961E+3_RP, 1.003E+4_RP, 1.012E+4_RP, & ! RH = 0.80, IWA =  532nm
                       1.501E+4_RP, 1.524E+4_RP, 1.555E+4_RP, & ! RH = 0.90, IWA =  532nm
                       2.212E+4_RP, 2.259E+4_RP, 2.321E+4_RP, & ! RH = 0.95, IWA =  532nm
                       5.128E+4_RP, 5.252E+4_RP, 5.420E+4_RP, & ! RH = 0.98, IWA =  532nm
                       6.867E+4_RP, 7.036E+4_RP, 7.259E+4_RP, & ! RH = 0.99, IWA =  532nm
                       1.481E+3_RP, 1.286E+3_RP, 9.986E+2_RP, & ! RH = 0.00, IWA = 1064nm
                       1.703E+3_RP, 1.486E+3_RP, 1.206E+3_RP, & ! RH = 0.50, IWA = 1064nm
                       1.764E+3_RP, 1.542E+3_RP, 1.265E+3_RP, & ! RH = 0.70, IWA = 1064nm
                       3.257E+3_RP, 3.057E+3_RP, 2.813E+3_RP, & ! RH = 0.80, IWA = 1064nm
                       5.152E+3_RP, 5.026E+3_RP, 4.835E+3_RP, & ! RH = 0.90, IWA = 1064nm
                       8.296E+3_RP, 8.235E+3_RP, 8.159E+3_RP, & ! RH = 0.95, IWA = 1064nm
                       2.571E+4_RP, 2.614E+4_RP, 2.669E+4_RP, & ! RH = 0.98, IWA = 1064nm
                       3.956E+4_RP, 4.033E+4_RP, 4.136E+4_RP, & ! RH = 0.99, IWA = 1064nm
                       6.407E+3_RP, 6.564E+3_RP, 6.805E+3_RP, & ! RH = 0.00, IWA =  355nm
                       7.432E+3_RP, 7.646E+3_RP, 7.929E+3_RP, & ! RH = 0.50, IWA =  355nm
                       7.718E+3_RP, 7.945E+3_RP, 8.232E+3_RP, & ! RH = 0.70, IWA =  355nm
                       1.375E+4_RP, 1.420E+4_RP, 1.482E+4_RP, & ! RH = 0.80, IWA =  355nm
                       1.976E+4_RP, 2.036E+4_RP, 2.117E+4_RP, & ! RH = 0.90, IWA =  355nm
                       2.721E+4_RP, 2.801E+4_RP, 2.909E+4_RP, & ! RH = 0.95, IWA =  355nm
                       5.250E+4_RP, 5.383E+4_RP, 5.561E+4_RP, & ! RH = 0.98, IWA =  355nm
                       6.594E+4_RP, 6.755E+4_RP, 6.968E+4_RP /), & ! RH = 0.99, IWA =  355nm
                    (/ 3, NRH, IWA /) )
    real(RP), parameter :: &
         CABSCC(3,NRH,IWA) &    !< absorption  cross section (532/1064/355nm)
         !             BC/OC=0.30   BC/OC=0.15   BC/OC=0.00
         = reshape( (/ 1.363E+3_RP, 9.118E+2_RP, 1.184E+2_RP, & ! RH = 0.00, IWA =  532nm
                       1.389E+3_RP, 9.232E+2_RP, 1.168E+2_RP, & ! RH = 0.50, IWA =  532nm
                       1.397E+3_RP, 9.240E+2_RP, 1.165E+2_RP, & ! RH = 0.70, IWA =  532nm
                       1.515E+3_RP, 9.629E+2_RP, 1.146E+2_RP, & ! RH = 0.80, IWA =  532nm
                       1.584E+3_RP, 9.891E+2_RP, 1.148E+2_RP, & ! RH = 0.90, IWA =  532nm
                       1.635E+3_RP, 1.008E+3_RP, 1.155E+2_RP, & ! RH = 0.95, IWA =  532nm
                       1.739E+3_RP, 1.053E+3_RP, 1.178E+2_RP, & ! RH = 0.98, IWA =  532nm
                       1.765E+3_RP, 1.067E+3_RP, 1.188E+2_RP, & ! RH = 0.99, IWA =  532nm
                       6.494E+2_RP, 4.340E+2_RP, 1.311E+2_RP, & ! RH = 0.00, IWA = 1064nm
                       6.676E+2_RP, 4.482E+2_RP, 1.340E+2_RP, & ! RH = 0.50, IWA = 1064nm
                       6.722E+2_RP, 4.517E+2_RP, 1.347E+2_RP, & ! RH = 0.70, IWA = 1064nm
                       7.389E+2_RP, 4.900E+2_RP, 1.440E+2_RP, & ! RH = 0.80, IWA = 1064nm
                       7.784E+2_RP, 5.107E+2_RP, 1.495E+2_RP, & ! RH = 0.90, IWA = 1064nm
                       8.094E+2_RP, 5.308E+2_RP, 1.546E+2_RP, & ! RH = 0.95, IWA = 1064nm
                       8.787E+2_RP, 5.713E+2_RP, 1.653E+2_RP, & ! RH = 0.98, IWA = 1064nm
                       9.000E+2_RP, 5.856E+2_RP, 1.691E+2_RP, & ! RH = 0.99, IWA = 1064nm
                       1.955E+3_RP, 1.385E+3_RP, 1.713E+2_RP, & ! RH = 0.00, IWA =  355nm
                       2.009E+3_RP, 1.416E+3_RP, 1.668E+2_RP, & ! RH = 0.50, IWA =  355nm
                       2.025E+3_RP, 1.418E+3_RP, 1.659E+2_RP, & ! RH = 0.70, IWA =  355nm
                       2.277E+3_RP, 1.500E+3_RP, 1.594E+2_RP, & ! RH = 0.80, IWA =  355nm
                       2.425E+3_RP, 1.548E+3_RP, 1.583E+2_RP, & ! RH = 0.90, IWA =  355nm
                       2.519E+3_RP, 1.582E+3_RP, 1.580E+2_RP, & ! RH = 0.95, IWA =  355nm
                       2.717E+3_RP, 1.655E+3_RP, 1.594E+2_RP, & ! RH = 0.98, IWA =  355nm
                       2.756E+3_RP, 1.676E+3_RP, 1.601E+2_RP /), & ! RH = 0.99, IWA =  355nm
                    (/ 3, NRH, IWA /) )
    real(RP), parameter :: &
         CBAKCC(3,NRH,IWA) &    !< backscatter cross section (532/1064/355nm)
         !             BC/OC=0.30   BC/OC=0.15   BC/OC=0.00
         = reshape( (/ 2.818E+1_RP, 3.779E+1_RP, 6.832E+1_RP, & ! RH = 0.00, IWA =  532nm
                       3.250E+1_RP, 4.434E+1_RP, 7.496E+1_RP, & ! RH = 0.50, IWA =  532nm
                       3.397E+1_RP, 4.604E+1_RP, 7.665E+1_RP, & ! RH = 0.70, IWA =  532nm
                       6.908E+1_RP, 8.797E+1_RP, 1.221E+2_RP, & ! RH = 0.80, IWA =  532nm
                       1.121E+2_RP, 1.343E+2_RP, 1.698E+2_RP, & ! RH = 0.90, IWA =  532nm
                       1.696E+2_RP, 1.949E+2_RP, 2.361E+2_RP, & ! RH = 0.95, IWA =  532nm
                       4.334E+2_RP, 4.725E+2_RP, 5.335E+2_RP, & ! RH = 0.98, IWA =  532nm
                       6.325E+2_RP, 6.882E+2_RP, 7.668E+2_RP, & ! RH = 0.99, IWA =  532nm
                       2.372E+1_RP, 2.412E+1_RP, 2.566E+1_RP, & ! RH = 0.00, IWA = 1064nm
                       2.619E+1_RP, 2.709E+1_RP, 2.887E+1_RP, & ! RH = 0.50, IWA = 1064nm
                       2.689E+1_RP, 2.789E+1_RP, 2.974E+1_RP, & ! RH = 0.70, IWA = 1064nm
                       4.236E+1_RP, 4.521E+1_RP, 4.893E+1_RP, & ! RH = 0.80, IWA = 1064nm
                       6.062E+1_RP, 6.407E+1_RP, 6.984E+1_RP, & ! RH = 0.90, IWA = 1064nm
                       8.787E+1_RP, 9.369E+1_RP, 1.017E+2_RP, & ! RH = 0.95, IWA = 1064nm
                       2.366E+2_RP, 2.475E+2_RP, 2.642E+2_RP, & ! RH = 0.98, IWA = 1064nm
                       3.485E+2_RP, 3.656E+2_RP, 3.864E+2_RP, & ! RH = 0.99, IWA = 1064nm
                       2.863E+1_RP, 5.075E+1_RP, 1.508E+2_RP, & ! RH = 0.00, IWA =  355nm
                       3.509E+1_RP, 6.192E+1_RP, 1.527E+2_RP, & ! RH = 0.50, IWA =  355nm
                       3.720E+1_RP, 6.417E+1_RP, 1.540E+2_RP, & ! RH = 0.70, IWA =  355nm
                       8.521E+1_RP, 1.231E+2_RP, 2.018E+2_RP, & ! RH = 0.80, IWA =  355nm
                       1.413E+2_RP, 1.826E+2_RP, 2.639E+2_RP, & ! RH = 0.90, IWA =  355nm
                       2.134E+2_RP, 2.637E+2_RP, 3.587E+2_RP, & ! RH = 0.95, IWA =  355nm
                       6.616E+2_RP, 7.594E+2_RP, 9.395E+2_RP, & ! RH = 0.98, IWA =  355nm
                       1.072E+3_RP, 1.227E+3_RP, 1.467E+3_RP /), & ! RH = 0.99, IWA =  355nm
                    (/ 3, NRH, IWA /) )
    real(RP), parameter :: &
         RADCA(NRH)        &    !< carbon particle radius [m]
         = (/ 0.100E-6_RP, 0.108E-6_RP, 0.110E-6_RP, 0.144E-6_RP, & ! RH = (0.00,0.50,0.70,0.80) 
              0.169E-6_RP, 0.196E-6_RP, 0.274E-6_RP, 0.312E-6_RP /) ! RH = (0.90,0.95,0.98,0.99)
    real(RP), parameter :: &
         RADSA(4)          &    !< sea salt particle radius [um]
         = (/ 3.560E-1_RP,  1.124E+0_RP,  3.560E+0_RP,  11.24E+0_RP /)
    real(RP), parameter :: BCEXT        = 0.500E+0_RP                   !< BC external factor
    real(RP), parameter :: AFOCBC       = 3.000E+0_RP                   !< aircraft emission factor BC/OC
    real(RP), parameter :: FOMFF        = 1.600E+0_RP                   !< conversion rate from OC to OM for fossil fuel
    real(RP), parameter :: FOMBB        = 2.600E+0_RP                   !< conversion rate from OC to OM for biomass burning
    real(RP), parameter :: DENSBC       = 2.300E+3_RP                   !< BC density                                [kg/m3]
    real(RP), parameter :: DENSCA       = 1.800E+3_RP                   !< density of carbonaceous aerosols          [kg/m3]
    real(RP), parameter :: DRYRBC       = 0.118E-7_RP                   !< BC (external mixture) dry particle radius [m]
    real(RP), parameter :: DRYROC       = 0.100E-6_RP                   !< particle radius of natural OC             [m]
    real(RP), parameter :: AA           = 1.000E-4_RP                   !< collision efficiency
    real(RP), parameter :: MOLMC (NSOA) = (/ 120.0E+0_RP, 60.0E+0_RP /) !< C mass of terpene/isoprene                [g/mol]
    real(RP), parameter :: MASSHC(NSOA) = (/ 136.0E+0_RP, 68.1E+0_RP /) !< mass of BVOCs                             [g/mol]
    real(RP), parameter :: EFAISO(NSOA) = (/ 0.30E-3_RP, 0.30E-1_RP /)  !< emission factor of isopren
    real(RP), parameter :: EFRISO       = 0.70E+0_RP                    !< emission fraction of isopren
    
    ! others 
    !-------------------------
    integer :: k, i, j, iq, iw
    integer :: IRH
    integer :: nyear
    integer :: nint
    real(RP) :: check_2d_data   (IA,JA,100)
    real(RP) :: check_3d_data(KA,IA,JA,100)

    integer :: nyear_p, nyear_m
    
    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - carbon'

    !setup parameter 
    !================================================

    !set daily data
    !-------------------------------
    nyear = TIME_NOWDATE(1)-1998

!    if( nyear == 14 ) then
!      nyear_p = nyear 
    if( nyear >= 22 ) then
      nyear_p = 22  ! tentative 
      nyear = 22   ! tentative 
    else
      nyear_p = nyear + 1
    endif
!    if( nyear == 1 ) then
!      nyear_m = nyear
    if( nyear <= 1 ) then
      nyear_m = 1   ! tentative 
      nyear = 1   ! tentative 
    else
      nyear_m = nyear - 1
    endif

    !BBBC
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( BBBC(:,:),                                                                   & ! [OUT]
!           BBBC_DATA(:,:,nyear-1,12), BBBC_DATA(:,:,nyear,:), BBBC_DATA(:,:,nyear+1,1), & ! [IN]
           BBBC_DATA(:,:,nyear_m,12), BBBC_DATA(:,:,nyear,:), BBBC_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                            ) ! [IN]

    !BBOC
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( BBOC(:,:),                                                                   & ! [OUT]
!           BBOC_DATA(:,:,nyear-1,12), BBOC_DATA(:,:,nyear,:), BBOC_DATA(:,:,nyear+1,1), & ! [IN]
           BBOC_DATA(:,:,nyear_m,12), BBOC_DATA(:,:,nyear,:), BBOC_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                            ) ! [IN]

    !ANTBC
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( ANTBC(:,:),                                                                     & ! [OUT]
!           ANTBC_DATA(:,:,nyear-1,12), ANTBC_DATA(:,:,nyear,:), ANTBC_DATA(:,:,nyear+1,1), & ! [IN]
           ANTBC_DATA(:,:,nyear_m,12), ANTBC_DATA(:,:,nyear,:), ANTBC_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                               ) ! [IN]

    !ANTOC
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( ANTOC(:,:),                                                                     & ! [OUT]
!           ANTOC_DATA(:,:,nyear-1,12), ANTOC_DATA(:,:,nyear,:), ANTOC_DATA(:,:,nyear+1,1), & ! [IN]
           ANTOC_DATA(:,:,nyear_m,12), ANTOC_DATA(:,:,nyear,:), ANTOC_DATA(:,:,nyear_p,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)                               ) ! [IN]

    !ISOP
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( ISOP(:,:),                                             & ! [OUT]
           ISOP_DATA(:,:,12), ISOP_DATA(:,:,:), ISOP_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)      ) ! [IN]

    !TERP
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( TERP(:,:),                                             & ! [OUT]
           TERP_DATA(:,:,12), TERP_DATA(:,:,:), TERP_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)      ) ! [IN]

    !CHLOA
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( CHLOA(:,:),                                               & ! [OUT]
           CHLOA_DATA(:,:,12), CHLOA_DATA(:,:,:), CHLOA_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)         ) ! [IN]

    !DIFATT
    call ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
         ( DIFATT(:,:),                                                 & ! [OUT]
           DIFATT_DATA(:,:,12), DIFATT_DATA(:,:,:), DIFATT_DATA(:,:,1), & ! [IN]
           TIME_NOWDATE(1), TIME_NOWDATE(2), TIME_NOWDATE(3)            ) ! [IN]
    
    do j = JS, JE
    do i = IS, IE
         BBBC(i,j) = max(   BBBC(i,j), 0.0_RP )
         BBOC(i,j) = max(   BBOC(i,j), 0.0_RP )
        ANTBC(i,j) = max(  ANTBC(i,j), 0.0_RP )
        ANTOC(i,j) = max(  ANTOC(i,j), 0.0_RP )
         ISOP(i,j) = max(   ISOP(i,j) * MASSHC(2) / MOLMC(2), 0.0_RP ) ! [kgC/m2/s] -> [kg/m2/s]
         TERP(i,j) = max(   TERP(i,j) * MASSHC(1) / MOLMC(1), 0.0_RP ) ! [kgC/m2/s] -> [kg/m2/s]
        CHLOA(i,j) = max(  CHLOA(i,j), 0.0_RP )
       DIFATT(i,j) = max( DIFATT(i,j), 0.0_RP )
    enddo
    enddo

    
    !CEXTP, CABSP, CEXTA, CABSA, CBAKA, RADC
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
             RADC (k,i,j)           = RADCA     (IRH)
             CEXTP(k,i,j,1:3,1:IWA) = CEXTCP(1:3,IRH,1:IWA)
             CABSP(k,i,j,1:3,1:IWA) = CABSCP(1:3,IRH,1:IWA)
             CEXTA(k,i,j,1:3,1:IWA) = CEXTCC(1:3,IRH,1:IWA)
             CABSA(k,i,j,1:3,1:IWA) = CABSCC(1:3,IRH,1:IWA)
             CBAKA(k,i,j,1:3,1:IWA) = CBAKCC(1:3,IRH,1:IWA)
          else
             BRH1 =  RH(k,i,j) - BRH(IRH)
             BRH2 = BRH(IRH+1) -  RH(k,i,j)
             BRH3 = BRH(IRH+1) - BRH(IRH)
             RADC(k,i,j) = ( BRH1 * RADCA(IRH+1) + BRH2 * RADCA(IRH) ) / BRH3
             do iw = 1, IWA
                do iq = 1, 3
                   CEXTP(k,i,j,iq,iw) = ( BRH1 * CEXTCP(iq,IRH+1,iw) + BRH2 * CEXTCP(iq,IRH,iw) ) / BRH3
                   CABSP(k,i,j,iq,iw) = ( BRH1 * CABSCP(iq,IRH+1,iw) + BRH2 * CABSCP(iq,IRH,iw) ) / BRH3
                   CEXTA(k,i,j,iq,iw) = ( BRH1 * CEXTCC(iq,IRH+1,iw) + BRH2 * CEXTCC(iq,IRH,iw) ) / BRH3
                   CABSA(k,i,j,iq,iw) = ( BRH1 * CABSCC(iq,IRH+1,iw) + BRH2 * CABSCC(iq,IRH,iw) ) / BRH3
                   CBAKA(k,i,j,iq,iw) = ( BRH1 * CBAKCC(iq,IRH+1,iw) + BRH2 * CBAKCC(iq,IRH,iw) ) / BRH3
                enddo
             enddo
          endif

          CEXTP(k,i,j,4,1:IWA) = CEXTBP(1:IWA)
          CABSP(k,i,j,4,1:IWA) = CABSBP(1:IWA)
          CEXTA(k,i,j,4,1:IWA) = CEXTBA(1:IWA)
          CABSA(k,i,j,4,1:IWA) = CABSBA(1:IWA)
          CBAKA(k,i,j,4,1:IWA) = CBAKBA(1:IWA)
    
       enddo
    enddo
    enddo

    !DRYMCB, VTERBC, VTER
    !-------------------------------
    DRYMCB( 1 ) = 4.0_RP / 3.0_RP * PI * RADCA(1)**3 * DENSCA * exp( 4.5_RP * log(SIGMA(1))**2 )
    DRYMCB( 2 ) = 4.0_RP / 3.0_RP * PI * RADCA(1)**3 * DENSCA * exp( 4.5_RP * log(SIGMA(2))**2 )
    DRYMCB( 3 ) = 4.0_RP / 3.0_RP * PI *   DRYRBC**3 * DENSBC * exp( 4.5_RP * log(SIGMA(3))**2 )
    VTERBC      = 2.0_RP * DENSBC * DRYRBC**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(3))**2 )

    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( iq /= NCA1 ) then
                VTER(k,i,j,iq) = 2.0_RP * (               ( RADCA(1)/RADC(k,i,j) )**3   * DENSCA   &
                                            +  ( 1.0_RP - ( RADCA(1)/RADC(k,i,j) )**3 ) * DENSW  ) &
                               * RADC(k,i,j)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(1))**2 )
             else
                VTER(k,i,j,iq) = VTERBC
             endif
          enddo
       enddo
       enddo
    enddo

    !================================================
    !end setup parameter


             
    !aerosol transport 
    !================================================
    
    !in rain(MP), in snow(MP)    
    !-------------------------------
    do iq = 1, 4
       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
    enddo

    
    !emission
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       EMITBC(i,j)   =   BBBC(i,j)         + ANTBC(i,j)
       EMITOC(i,j)   =   BBOC(i,j) * FOMBB + ANTOC(i,j) * FOMFF
       QEMIT (i,j,1) =   BBOC(i,j) * FOMBB * GRAV * real(dt_AE,kind=RP) / DPKBB(i,j)                      ! BC(BB)
       QEMIT (i,j,2) =   BBBC(i,j)         * GRAV * real(dt_AE,kind=RP) / DPKBB(i,j)                      ! OC(BB)
       QEMIT (i,j,3) =  ANTOC(i,j) * FOMFF * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j)                      ! OC(FF)
       QEMIT (i,j,4) =  ANTBC(i,j)         * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j) * ( 1.0_RP - BCEXT ) ! BC(FF;int.)
       QEMIT (i,j,5) =  ANTBC(i,j)         * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j) * BCEXT              ! BC(FF;ext.)
       QEMIT (i,j,6) =  TERP (i,j)         * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j)                      ! terpene
       QEMIT (i,j,7) =  ISOP (i,j)         * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j)                      ! isoprene
    enddo
    enddo

    !oceanic primary OM (Gantt et al., ACP, 2011)
    EMISSA(:,:) = 0.0_RP
    do iq = 1, NSA
       do j = JS, JE
       do i = IS, IE
          if ( FACT_OCEAN(i,j) > 0.0_RP ) then
             OMSSA(i,j,iq) = ( 1.0_RP  / ( 1.0_RP + exp( 0.18_RP * WS10(i,j) - 2.63_RP * CHLOA(i,j) ) ) ) &
                           / ( 1.0_RP + 0.03_RP * exp( 6.81_RP * RADSA(iq) ) ) &
                           +   0.03_RP / ( 1.0_RP + exp( 0.18_RP * WS10(i,j) - 2.63_RP * CHLOA(i,j) ) )
          else
             OMSSA(i,j,iq) = 0.0_RP
          endif
          OMSSA (i,j,iq) = min( max( OMSSA(i,j,iq), 0.0_RP ), 1.0_RP )
          EMITSS(i,j,iq) = OMSSA(i,j,iq) * EMITSA(i,j,iq)
          EMISSA(i,j)    = EMISSA(i,j) + EMITSS(i,j,iq)
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       QEMIT (i,j,8) = EMISSA(i,j)         * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j)                      ! oceanic OM
    enddo
    enddo
          
    !oceanic OM precursors (Gantt et al., ACP, 2009)
    do j = JS, JE
    do i = IS, IE
       if (           DIFATT(i,j) > 0.0_RP &
            .and. FACT_OCEAN(i,j) > 0.0_RP & 
            .and.     RFLXSD(i,j) > 0.0_RP ) then
          DIFAT(i,j) = max( DIFATT(i,j), 1.9E-2_RP )
          RMAG (i,j) = log( RFLXSD(i,j) ) 
          HMAX (i,j) = max( (RMAG(i,j)-log(2.5_RP))/DIFAT(i,j), 0.0_RP )
       else
          HMAX(i,j) = 0.0_RP
          RMAG(i,j) = 0.0_RP
          DIFAT(i,j) = 1.9E-2_RP
       endif
    enddo
    enddo

    do iq = 1, 2
       do j = JS, JE
       do i = IS, IE
          ! PISO = EFAISO(iq) * ( 2.0_RP*HMAX(i,j)*(log(2.0_RP)+RMAG(i,j)) - DIFAT(i,j)*HMAX(i,j)**2 )
          PISO = 0.0_RP
          if ( HMAX(i,j) > 0.0_RP ) then
             DELTSP = 1.0_RP / 20.0_RP * HMAX(i,j)
             do nint = 1, 20
                PISO = PISO + DELTSP * EFAISO(iq) * ( max( log(2.0_RP) + RMAG(i,j) - DIFAT(i,j) * nint * DELTSP, 0.0_RP ) )**2
             enddo
          endif

          PISO = max( PISO, 0.0_RP )
          EMIISO(i,j,  iq) = HMAX(i,j) * CHLOA(i,j) * 1.0E-3_RP * EFRISO * PISO * 1.0E-6_RP * MASSHC(iq) * 1.0E-3_RP / 3.6E+3_RP &
                           * FACT_OCEAN(i,j) * ( 1.0_RP - GOICE(i,j) )
          EMIOBS(i,j,  iq) = EMIISO(i,j,iq) * 1.0E+3_RP / MASSHC(iq) / 1.0E+4_RP * AVOG
          QEMIT (i,j,8+iq) = EMIISO(i,j,iq) * GRAV * real(dt_AE,kind=RP) / DPKUP(i,j)                     ! oceanic OM precursors 
       enddo
       enddo
    enddo
    
    do j = JS, JE
    do i = IS, IE
       do k = KS, KBB(i,j)
          QTRC(k,i,j,1) = QTRC(k,i,j,1) + QEMIT(i,j,1)     ! OC : biomass burning (internal mixture)
          QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j,2)     ! BC : biomass burning (internal mixture)
       enddo
       do k = KS, KUPMAX
          if ( k <= KUP(i,j) ) then
             QTRC(k,i,j,1) = QTRC(k,i,j,1) + QEMIT(i,j, 3) ! OC : anthropogenic source (internal mixture)
             QTRC(k,i,j,2) = QTRC(k,i,j,2) + QEMIT(i,j, 4) ! BC : anthropogenic source (internal mixture)
             QTRC(k,i,j,3) = QTRC(k,i,j,3) + QEMIT(i,j, 8) ! oceanic primary OM : secondary organic carbon (SOC) from natural volatile organic compound (NVOC)
             QTRC(k,i,j,4) = QTRC(k,i,j,4) + QEMIT(i,j, 5) ! BC : anthropogenic source (external mixture)
             QTRC(k,i,j,5) = QTRC(k,i,j,5) + QEMIT(i,j, 6) ! terpene
             QTRC(k,i,j,5) = QTRC(k,i,j,5) + QEMIT(i,j, 9) ! oceanic terpene
             QTRC(k,i,j,6) = QTRC(k,i,j,6) + QEMIT(i,j, 7) ! isoprene
             QTRC(k,i,j,6) = QTRC(k,i,j,6) + QEMIT(i,j,10) ! oceanic isoprene
          endif
       enddo
    enddo
    enddo

    AFFBC(:,:,:) = 1.0_RP       ! BC emission factor for aircrafts
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QTRC(k,i,j,1) + AIRFELOC(k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP) * FOMFF
          QTRC(k,i,j,2) = QTRC(k,i,j,2) + AIRFEL  (k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP) * ( 1.0_RP - BCEXT ) * AFFBC(k,i,j)
          QTRC(k,i,j,4) = QTRC(k,i,j,4) + AIRFEL  (k,i,j) / DENS(k,i,j) * real(dt_AE,kind=RP) * BCEXT              * AFFBC(k,i,j)
          EMITBC(i,j)   = EMITBC(i,j)   + AIRFEL  (k,i,j) * DELZ(k,i,j) * AFFBC(k,i,j)
          EMITOC(i,j)   = EMITOC(i,j)   + AIRFELOC(k,i,j) * DELZ(k,i,j) * FOMFF
       enddo
    enddo
    enddo

    
    !chemical reaction for SOA
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          ! NUMHC (k,i,j,1) = QTRC (k,i,j,5) * DENS(k,i,j) * 1.0E-03_RP * AVOG / MASSHC(1) ! terpene  [molecule/cm3]
          ! NUMHC (k,i,j,2) = QTRC (k,i,j,6) * DENS(k,i,j) * 1.0E-03_RP * AVOG / MASSHC(2) ! isoprene [molecule/cm3]
          ! NUMO3 (k,i,j)   = PPMO3(k,i,j)   * DENS(k,i,j) * 1.0E-09_RP * AVOG / MOLAIR    ! ozone    [molecule/cm3]
          ! NUMNO3(k,i,j)   = NO3G (k,i,j)   * DENS(k,i,j) * 1.0E-15_RP * AVOG / MOLAIR    ! NO3      [molecule/cm3]
          ! CONOA (k,i,j)   = ( QTRC(k,i,j,1) + QTRC(k,i,j,3) ) * DENS(k,i,j) * 1.0E+09_RP ! OC       [ug/m3]
          NUMHC (k,i,j,1) = max( QTRC (k,i,j,5) * DENS(k,i,j) * 1.0E-03_RP * AVOG / MASSHC(1), 0.0_RP ) ! terpene  [molecule/cm3]
          NUMHC (k,i,j,2) = max( QTRC (k,i,j,6) * DENS(k,i,j) * 1.0E-03_RP * AVOG / MASSHC(2), 0.0_RP ) ! isoprene [molecule/cm3]
          NUMO3 (k,i,j)   = max( PPMO3(k,i,j)   * DENS(k,i,j) * 1.0E-09_RP * AVOG / MOLAIR,    0.0_RP ) ! ozone    [molecule/cm3]
          NUMNO3(k,i,j)   = max( NO3G (k,i,j)   * DENS(k,i,j) * 1.0E-15_RP * AVOG / MOLAIR,    0.0_RP ) ! NO3      [molecule/cm3]
          CONOA (k,i,j)   = max( ( QTRC(k,i,j,1) + QTRC(k,i,j,3) ) * DENS(k,i,j) * 1.0E+09_RP, 0.0_RP ) ! OC       [ug/m3]
       enddo
    enddo
    enddo

    call ATMOS_PHY_AE_SPRINTARS_AEROCA_FSOA &
         ( NUMHC,                       & ! [MODIFIED]
           NEWSOA,                      & ! [OUT]
           NUMOH, NUMO3, NUMNO3, CONOA, & ! [IN]
           GDT, dt_AE                   ) ! [IN]

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,3) = QTRC (k,i,j,3) + NEWSOA(k,i,j) * 1.0E-9_RP / DENS(k,i,j)
          QTRC(k,i,j,5) = NUMHC(k,i,j,1) / DENS(k,i,j) * 1.0E+3_RP / AVOG * MASSHC(1)
          QTRC(k,i,j,6) = NUMHC(k,i,j,2) / DENS(k,i,j) * 1.0E+3_RP / AVOG * MASSHC(2)
       enddo
    enddo
    enddo

    CHESOA(:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          CHESOA(i,j) = CHESOA(i,j) + NEWSOA(k,i,j) * 1.0E-9_RP * DELZ(k,i,j) / real(dt_AE,kind=RP)
       enddo
    enddo
    enddo

    
    !internal mixed OC & BC
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QTRCC(k,i,j,1) = QTRC(k,i,j,1) + QTRC(k,i,j,2) ! internal mixed OC & BC
          QTRCC(k,i,j,2) = QTRC(k,i,j,3)                 ! secondary organic carbon (SOC) from NVOC
          QTRCC(k,i,j,3) = QTRC(k,i,j,4)                 ! external mixed BC
          QTRCC(k,i,j,4) = QTRC(k,i,j,5)                 ! terpene
          QTRCC(k,i,j,5) = QTRC(k,i,j,6)                 ! isoprene
          if ( QTRCC(k,i,j,1) > 0.0_RP ) then
             RCABC(k,i,j) = QTRC(k,i,j,2) / QTRCC(k,i,j,1)
          else
             RCABC(k,i,j) = 0.0_RP
          endif
       enddo
    enddo
    enddo


    ! !instability
    ! !-------------------------------
    ! do iq = 1, NCA1+NSOA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB( QTRCC(:,:,:,iq),                               & ! [MODIFIED]
    !                                               ONEWQ(1:KA-1,:,:), DELP(:,:,:), GDPM(0:KA,:,:) ) ! [IN]
    ! enddo

    
    !incloud
    !-------------------------------
    !inside cloud
    if ( SW_CCN ) then
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCTM2(k,i,j,1) = QTRCC(k,i,j,1) * ACTICA(k,i,j) * TCLDF(k,i,j)
             QTRCTM2(k,i,j,2) = QTRCC(k,i,j,2) * ACTIOC(k,i,j) * TCLDF(k,i,j)
             QTRCTM2(k,i,j,3) = QTRCC(k,i,j,3) * FINCCA(3)     * TCLDF(k,i,j)
          enddo
       enddo
       enddo
    else
       do iq = 1, NCA1
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRCC(k,i,j,iq) * FINCCA(iq) * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    endif

    !outside cloud
    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCTM1(k,i,j,iq) = QTRCC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
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

    ! do iq = 1, NCA1
    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          CARBN = QTRCTM1(k,i,j,iq) / DRYMCB(iq) * DENS(k,i,j)
    !          if ( iq <= NCA1-1 ) then
    !             CBRN(k,i,j,iq) = AA * PI * ( RADR + RADC(k,i,j) )**2 * ( VTR - VTER(k,i,j,iq) ) * CARBN * RAINN(k,i,j)
    !          elseif ( iq == NCA1 ) then
    !             CBRN(k,i,j,iq) = AA * PI * ( RADR + DRYRBC      )**2 * ( VTR - VTER(k,i,j,iq) ) * CARBN * RAINN(k,i,j)
    !          endif
    !          QWTDEP (k,i,j,iq) = CBRN(k,i,j,iq) * DRYMCB(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
    !          QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP(k,i,j,iq)
    !          if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
    !             QTRCTMP(k,i,j,iq) = 0.0_RP
    !             CBRN   (k,i,j,iq) = CARBN / real(dt_AE,kind=RP)
    !             QWTDEP (k,i,j,iq) = QTRCTM1(k,i,j,iq)
    !          endif
    !       enddo
    !    enddo
    !    enddo
    ! enddo

    
    ! !incloud scavenging (rain out)
    ! !-------------------------------
    ! do iq = 1, NCA1
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
    ! do iq = 1, NCA1
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER                         &
    !         ( QTRCTM1(:,:,:,iq), QWTDEP(:,:,:,iq),                       & ! [MODIFIED]
    !           WETF(:,:,iq),                                              & ! [OUT]
    !           CBRN(:,:,:,iq), DRYMCB(iq), CLTORN(:,:,:,iq),              & ! [IN]
    !           DENS(:,:,:), FRAIN(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]

    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          QTRCC(k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
    !       enddo
    !    enddo
    !    enddo
    
    ! enddo

    
    ! !dry deposition
    ! !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    CDVEC(i,j) = 1.0E-3_RP * CDVE(i,j) / ( 1.0E-3_RP + CDVE(i,j) )
    ! enddo
    ! enddo

    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !      ( QMTX,                    & ! [OUT]
    !        DENS, DELP, CDVEC, dt_AE ) ! [IN]

    ! do iq = 1, NCA1+2
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRCC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF   (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVEC, dt_AE ) ! [IN]
    ! enddo

    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, NCA1
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !         ( QTRCC(:,:,:,iq),                         & ! [MODIFIED]
    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq),          & ! [OUT]
    !           GDPM, DELP, DENS, VTER(:,:,:,iq), dt_AE  ) ! [IN]
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

    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             CARBN = QTRCTM1(k,i,j,iq) / DRYMCB(iq) * DENS(k,i,j)
             if ( iq <= NCA1-1 ) then
                CBRN_CP     (k,i,j,iq) = AA * PI * ( RADR                + RADC(k,i,j) )**2 * abs( VTR                - VTER(k,i,j,iq) ) * CARBN * RAINN_CP     (k,i,j)
                CBRN_MP_RAIN(k,i,j,iq) = AA * PI * ( RADR_MP_RAIN(k,i,j) + RADC(k,i,j) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * CARBN * RAINN_MP_RAIN(k,i,j)
                CBRN_MP_SNOW(k,i,j,iq) = AA * PI * ( RADR_MP_SNOW(k,i,j) + RADC(k,i,j) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * CARBN * RAINN_MP_SNOW(k,i,j)
             elseif ( iq == NCA1 ) then
                CBRN_CP     (k,i,j,iq) = AA * PI * ( RADR                + DRYRBC      )**2 * abs( VTR                - VTER(k,i,j,iq) ) * CARBN * RAINN_CP     (k,i,j)
                CBRN_MP_RAIN(k,i,j,iq) = AA * PI * ( RADR_MP_RAIN(k,i,j) + DRYRBC      )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * CARBN * RAINN_MP_RAIN(k,i,j)
                CBRN_MP_SNOW(k,i,j,iq) = AA * PI * ( RADR_MP_SNOW(k,i,j) + DRYRBC      )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * CARBN * RAINN_MP_SNOW(k,i,j)
             endif
             QWTDEP_CP     (k,i,j,iq) = CBRN_CP     (k,i,j,iq) * DRYMCB(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_RAIN(k,i,j,iq) = CBRN_MP_RAIN(k,i,j,iq) * DRYMCB(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_SNOW(k,i,j,iq) = CBRN_MP_SNOW(k,i,j,iq) * DRYMCB(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
                QTRCTMP(k,i,j,iq) = 0.0_RP
                CBRN        (k,i,j,iq) = CBRN_CP(k,i,j,iq) + CBRN_MP_RAIN(k,i,j,iq) + CBRN_MP_SNOW(k,i,j,iq)
                CBRN_CP     (k,i,j,iq) = CARBN / real(dt_AE,kind=RP) * CBRN_CP     (k,i,j,iq) / CBRN(k,i,j,iq)
                CBRN_MP_RAIN(k,i,j,iq) = CARBN / real(dt_AE,kind=RP) * CBRN_MP_RAIN(k,i,j,iq) / CBRN(k,i,j,iq)
                CBRN_MP_SNOW(k,i,j,iq) = CARBN / real(dt_AE,kind=RP) * CBRN_MP_SNOW(k,i,j,iq) / CBRN(k,i,j,iq)
                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * CBRN_CP     (k,i,j,iq) / CBRN(k,i,j,iq)
                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * CBRN_MP_RAIN(k,i,j,iq) / CBRN(k,i,j,iq)
                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * CBRN_MP_SNOW(k,i,j,iq) / CBRN(k,i,j,iq)
             endif
          enddo
       enddo
       enddo
    enddo

    
    !incloud scavenging (rain out)
    !-------------------------------
    do iq = 1, NCA1
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
    do j = JS, JE
    do i = IS, IE
       CDVEC(i,j) = 1.0E-3_RP * CDVE(i,j) / ( 1.0E-3_RP + CDVE(i,j) )
    enddo
    enddo

    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
         ( QMTX,                    & ! [OUT]
           DENS, DELP, CDVEC, dt_AE ) ! [IN]

    do iq = 1, NCA1
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
            ( QTRCTM1(:,:,:,iq),             & ! [MODIFIED]
              DRYF(:,:,iq),                  & ! [OUT]
              QMTX, DENS, DELP, CDVEC, dt_AE ) ! [IN]
    enddo
    do iq = NCA1+1, NCA1+NSOA
       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
            ( QTRCC(:,:,:,iq),               & ! [MODIFIED]
              DRYF(:,:,iq),                  & ! [OUT]
              QMTX, DENS, DELP, CDVEC, dt_AE ) ! [IN]
    enddo

    
    !gravitational settling
    !-------------------------------
    do iq = 1, NCA1
       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
            ( QTRCTM1(:,:,:,iq),                      & ! [MODIFIED]
              QLOAD(:,:,:,iq), GRAVF(:,:,iq),         & ! [OUT]
              GDPM, DELP, DENS, VTER(:,:,:,iq), dt_AE ) ! [IN]
    enddo

    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QWTDEP_MP_RAIN(k,i,j,1) = QWTDEP_MP_RAIN(k,i,j,1) + QTRCINR(k,i,j,1) + QTRCINR(k,i,j,2)
          QWTDEP_MP_SNOW(k,i,j,1) = QWTDEP_MP_SNOW(k,i,j,1) + QTRCINS(k,i,j,1) + QTRCINS(k,i,j,2)
          QWTDEP_MP_RAIN(k,i,j,2) = QWTDEP_MP_RAIN(k,i,j,2) + QTRCINR(k,i,j,3)
          QWTDEP_MP_SNOW(k,i,j,2) = QWTDEP_MP_SNOW(k,i,j,2) + QTRCINS(k,i,j,3)
          QWTDEP_MP_RAIN(k,i,j,3) = QWTDEP_MP_RAIN(k,i,j,3) + QTRCINR(k,i,j,4)
          QWTDEP_MP_SNOW(k,i,j,3) = QWTDEP_MP_SNOW(k,i,j,3) + QTRCINS(k,i,j,4)
       enddo
    enddo
    enddo
 
    
    !re-emission from rain : CP
    !-------------------------------
    do iq = 1, NCA1
       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
              WETF_CP(:,:,iq),                                              & ! [OUT]
              CBRN_CP(:,:,:,iq), DRYMCB(iq), CLTORN_CP(:,:,:,iq),           & ! [IN]
              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
    enddo

    !re-emission from rain & snow : MP
    !-------------------------------
    do iq = 1, NCA1
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
       
       QTRCINRC(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-rain(MP)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
              WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
              QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
              GDZM, DELZ, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!        WETF_MP_RAIN(:,:,iq) = 0.0_RP
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCINRC(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINRC(k,i,j,iq), 0.0_RP )
          enddo
       enddo
       enddo
          
       QTRCINSC(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-snow(MP)
       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
              WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
              QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
              GDZM, DELZ, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCINSC(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINSC(k,i,j,iq), 0.0_RP )
          enddo
       enddo
       enddo
       
       do j = JS, JE
       do i = IS, IE
          WETF(i,j,iq) = WETF_CP(i,j,iq) + WETF_MP_RAIN(i,j,iq) + WETF_MP_SNOW(i,j,iq)
       enddo
       enddo
    enddo

    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QTRCINP(k,i,j,1) = ( QTRCINRC(k,i,j,1) + QTRCINSC(k,i,j,1) ) * ( 1.0_RP - RCABC(k,i,j) )
          QTRCINP(k,i,j,2) = ( QTRCINRC(k,i,j,1) + QTRCINSC(k,i,j,1) ) *            RCABC(k,i,j)
          QTRCINP(k,i,j,3) = ( QTRCINRC(k,i,j,2) + QTRCINSC(k,i,j,2) )
          QTRCINP(k,i,j,4) = ( QTRCINRC(k,i,j,3) + QTRCINSC(k,i,j,3) )
       enddo
    enddo
    enddo

    
    !inside cloud + outside cloud 
    !-------------------------------
    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCC(k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
             QLOAD(k,i,j,iq) = QTRCC  (k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    
    ! !dry deposition
    ! !-------------------------------
    ! do j = JS, JE
    ! do i = IS, IE
    !    CDVEC(i,j) = 1.0E-3_RP * CDVE(i,j) / ( 1.0E-3_RP + CDVE(i,j) )
    ! enddo
    ! enddo

    ! call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
    !      ( QMTX,                    & ! [OUT]
    !        DENS, DELP, CDVEC, dt_AE ) ! [IN]

    ! do iq = 1, NCA1+NSOA
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
    !         ( QTRCC(:,:,:,iq),               & ! [MODIFIED]
    !           DRYF   (:,:,iq),               & ! [OUT]
    !           QMTX, DENS, DELP, CDVEC, dt_AE ) ! [IN]
    ! enddo

    
    ! !gravitational settling
    ! !-------------------------------
    ! do iq = 1, NCA1
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !         ( QTRCC(:,:,:,iq),                        & ! [MODIFIED]
    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq),         & ! [OUT]
    !           GDPM, DELP, DENS, VTER(:,:,:,iq), dt_AE ) ! [IN]
    ! enddo

    
    !distribution to BC/OM(OC) = 0.30, 0.15, 0.00 for radiation
    !-------------------------------
    OUTQLD(:,:,:,:) = 0.0_RP
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QBC(k,i,j) = QLOAD(k,i,j,1) * RCABC(k,i,j)
          QOC(k,i,j) = QLOAD(k,i,j,1) * ( 1.0_RP - RCABC(k,i,j) ) + QLOAD(k,i,j,2)
          OUTQLD(k,i,j,I_BCOMBC) = QLOAD(k,i,j,NCA1)
          QC = QBC(k,i,j) + QOC(k,i,j)

          if ( QOC(k,i,j) > 0.0_RP ) then
             F = QBC(k,i,j) / QOC(k,i,j)
             if ( F > 0.15_RP ) then
                FF = min( (F-0.15_RP)/(0.30_RP-0.15_RP), 1.0_RP )
                OUTQLD(k,i,j,I_BCOM30) =            FF   * QC
                OUTQLD(k,i,j,I_BCOM15) = ( 1.0_RP - FF ) * QC
             else
                FF = min( (F-0.00_RP)/(0.15_RP-0.00_RP), 1.0_RP )
                OUTQLD(k,i,j,I_BCOM15) =            FF   * QC
                OUTQLD(k,i,j,I_BCOM00) = ( 1.0_RP - FF ) * QC
             endif
          endif

       enddo
    enddo
    enddo


    !internal mixed OC & BC
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          QTRC(k,i,j,1) = QTRCC(k,i,j,1) * ( 1.0_RP - RCABC(k,i,j) ) ! internal mixed OC
          QTRC(k,i,j,2) = QTRCC(k,i,j,1) *            RCABC(k,i,j)   ! internal mixed BC
          QTRC(k,i,j,3) = QTRCC(k,i,j,2)                             ! secondary organic carbon (SOC) from NVOC
          QTRC(k,i,j,4) = QTRCC(k,i,j,3)                             ! external mixed BC
          QTRC(k,i,j,5) = QTRCC(k,i,j,4)                             ! terpene
          QTRC(k,i,j,6) = QTRCC(k,i,j,5)                             ! isoprene
          QTRCOC(k,i,j) = QTRC(k,i,j,1) + QTRC(k,i,j,3)              ! total OC
          QTRCBC(k,i,j) = QTRC(k,i,j,2) + QTRC(k,i,j,4)              ! total BC
          WATROC(k,i,j) = ( RADC(k,i,j)**3 / RADCA(1)**3 - 1 ) * QTRCOC(k,i,j) * DENSW / DENSCA
       enddo
    enddo
    enddo

    !================================================
    !end aerosol transport

    
    !aerosol 
    !================================================
    !Reset 
    !-------------------------------
    !      k i j q
    QLOADT(:,:,:)   = 0.0_RP
    CCONT (:,:,:)   = 0.0_RP
    COLMS   (:,:,:) = 0.0_RP
    COLMST  (:,:)   = 0.0_RP
    COLBCT  (:,:)   = 0.0_RP
    COL1BC  (:,:)   = 0.0_RP

    !-------------------------------
    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QLOADT(k,i,j)    = QLOADT(k,i,j)    + QLOAD(k,i,j,iq)
             CCON  (k,i,j,iq) = QLOAD (k,i,j,iq) * DENS(k,i,j)
             NUMCON(k,i,j,iq) = CCON  (k,i,j,iq) * MTONCA(iq)
             CCONT (k,i,j)    = CCONT (k,i,j)    + CCON(k,i,j,iq)
             COLMS   (i,j,iq) = COLMS   (i,j,iq) + QLOAD(k,i,j,iq) * DELP(k,i,j) / GRAV
          enddo
          COLMST(i,j) = COLMST(i,j) + COLMS(i,j,iq)
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          BCCONT(k,i,j)   = ( QBC(k,i,j) + QLOAD(k,i,j,3) ) * DENS(k,i,j)
          OCCONT(k,i,j)   = max( CCONT(k,i,j)-BCCONT(k,i,j), 0.0_RP )
          NUMCON(k,i,j,2) = NUMCON(k,i,j,2) * ( RADCA(1) / DRYROC )**3
          COLBCT  (i,j)   = COLBCT  (i,j) + ( QBC(k,i,j) + QLOAD(k,i,j,3) ) * DELP(k,i,j) / GRAV
          COL1BC  (i,j)   = COL1BC  (i,j) +   QLOAD(k,i,j,1) * RCABC(k,i,j) * DELP(k,i,j) / GRAV
          HCPPB (k,i,j,1) = QTRCC (k,i,j,4) * MOLAIR / MASSHC(1) * 1.0E+9_RP
          HCPPB (k,i,j,2) = QTRCC (k,i,j,5) * MOLAIR / MASSHC(2) * 1.0E+9_RP
       enddo
       COLOCT(i,j) = max( COLMST(i,j)-COLBCT(i,j), 0.0_RP )
    enddo
    enddo
    !================================================
    !end aerosol     


    !flux 
    !================================================
    !Reset 
    !-------------------------------
    !              i j
    WETFT         (:,:) = 0.0_RP
    WETFT_CP      (:,:) = 0.0_RP
    WETFT_MP_RAIN (:,:) = 0.0_RP
    WETFT_MP_SNOW (:,:) = 0.0_RP
    DRYFT         (:,:) = 0.0_RP
    GRAVFT        (:,:) = 0.0_RP
    WETFOC        (:,:) = 0.0_RP
    WETFOC_CP     (:,:) = 0.0_RP
    WETFOC_MP_RAIN(:,:) = 0.0_RP
    WETFOC_MP_SNOW(:,:) = 0.0_RP
    DRYFOC        (:,:) = 0.0_RP
    GRVFOC        (:,:) = 0.0_RP
    WETFBC        (:,:) = 0.0_RP
    WETFBC_CP     (:,:) = 0.0_RP
    WETFBC_MP_RAIN(:,:) = 0.0_RP
    WETFBC_MP_SNOW(:,:) = 0.0_RP
    DRYFBC        (:,:) = 0.0_RP
    GRVFBC        (:,:) = 0.0_RP

    !-------------------------------
    do iq = 1, NCA1
       do j = JS, JE
       do i = IS, IE
          WETFT        (i,j) = WETFT         (i,j) + WETF        (i,j,iq)
          WETFT_CP     (i,j) = WETFT_CP      (i,j) + WETF_CP     (i,j,iq)
          WETFT_MP_RAIN(i,j) = WETFT_MP_RAIN (i,j) + WETF_MP_RAIN(i,j,iq)
          WETFT_MP_SNOW(i,j) = WETFT_MP_SNOW (i,j) + WETF_MP_SNOW(i,j,iq)
          DRYFT        (i,j) = DRYFT         (i,j) + DRYF        (i,j,iq)
          GRAVFT       (i,j) = GRAVFT        (i,j) + GRAVF       (i,j,iq)
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       if ( COLMS(i,j,1) > 0.0_RP ) then
          WETFOC        (i,j) = WETFOC        (i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * WETF        (i,j,1) + WETF        (i,j,2)
          WETFOC_CP     (i,j) = WETFOC_CP     (i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * WETF_CP     (i,j,1) + WETF_CP     (i,j,2)
          WETFOC_MP_RAIN(i,j) = WETFOC_MP_RAIN(i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * WETF_MP_RAIN(i,j,1) + WETF_MP_RAIN(i,j,2)
          WETFOC_MP_SNOW(i,j) = WETFOC_MP_SNOW(i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * WETF_MP_SNOW(i,j,1) + WETF_MP_SNOW(i,j,2)
          DRYFOC        (i,j) = DRYFOC        (i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * DRYF        (i,j,1) + DRYF        (i,j,2)
          GRVFOC        (i,j) = GRVFOC        (i,j) + ( COLMS(i,j,1) - COL1BC(i,j) ) / COLMS(i,j,1) * GRAVF       (i,j,1) + GRAVF       (i,j,2)
          WETFBC        (i,j) = WETFBC        (i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * WETF        (i,j,1) + WETF        (i,j,3)
          WETFBC_CP     (i,j) = WETFBC_CP     (i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * WETF_CP     (i,j,1) + WETF_CP     (i,j,3)
          WETFBC_MP_RAIN(i,j) = WETFBC_MP_RAIN(i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * WETF_MP_RAIN(i,j,1) + WETF_MP_RAIN(i,j,3)
          WETFBC_MP_SNOW(i,j) = WETFBC_MP_SNOW(i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * WETF_MP_SNOW(i,j,1) + WETF_MP_SNOW(i,j,3)
          DRYFBC        (i,j) = DRYFBC        (i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * DRYF        (i,j,1) + DRYF        (i,j,3)
          GRVFBC        (i,j) = GRVFBC        (i,j) +                  COL1BC(i,j)   / COLMS(i,j,1) * GRAVF       (i,j,1) + GRAVF       (i,j,3)
       endif
       DEPFT (i,j) = WETFT (i,j) + DRYFT (i,j) + GRAVFT(i,j)
       DEPFOC(i,j) = WETFOC(i,j) + DRYFOC(i,j) + GRVFOC(i,j)
       DEPFBC(i,j) = WETFBC(i,j) + DRYFBC(i,j) + GRVFBC(i,j)
    enddo
    enddo
    !================================================
    !end flux 


    !optical parameters
    !================================================
    !      k i j w 2
    TAUCA   (:,:,:,:) = 0.0_RP
    TAUCAA  (:,:,:,:) = 0.0_RP
    TAUOC   (:,:,  :) = 0.0_RP
    TAUOCA  (:,:,  :) = 0.0_RP
    TAUBC   (:,:,  :) = 0.0_RP
    TAUBCA  (:,:,  :) = 0.0_RP
    TAUCAD  (:,:)     = 0.0_RP
    CEXTCA(:,:,:,:,:) = 0.0_RP
    CABSCA(:,:,:,:,:) = 0.0_RP
    CBAKCA(:,:,:,:,:) = 0.0_RP
    CEXTCD(:,:,:)     = 0.0_RP
    
    !all sky
    !-------------------------------
    do iq = 1, 4
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             TUSE = OUTQLD(k,i,j,iq) * DELP(k,i,j) / GRAV
             do iw = 1, IWA
                TAUCA   (i,j,iw,1) = TAUCA   (i,j,iw,1) + TUSE * CEXTP(k,i,j,iq,iw)
                TAUCAA  (i,j,iw,1) = TAUCAA  (i,j,iw,1) + TUSE * CABSP(k,i,j,iq,iw)
                CEXTCA(k,i,j,iw,1) = CEXTCA(k,i,j,iw,1) + TUSE * CEXTA(k,i,j,iq,iw) / DELZ(k,i,j)
                CABSCA(k,i,j,iw,1) = CABSCA(k,i,j,iw,1) + TUSE * CABSA(k,i,j,iq,iw) / DELZ(k,i,j)
                CBAKCA(k,i,j,iw,1) = CBAKCA(k,i,j,iw,1) + TUSE * CBAKA(k,i,j,iq,iw) / DELZ(k,i,j)
             enddo
             if ( iq <= 3 ) then
                TAUCAD  (i,j) = TAUCAD  (i,j) + TUSE * CEXTCP(iq,1,1)
                CEXTCD(k,i,j) = CEXTCD(k,i,j) + TUSE * CEXTCC(iq,1,1) / DELZ(k,i,j)
             else
                TAUCAD  (i,j) = TAUCAD  (i,j) + TUSE * CEXTBP(1)
                CEXTCD(k,i,j) = CEXTCD(k,i,j) + TUSE * CEXTBA(1)      / DELZ(k,i,j)
             endif
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          TAUOC   (i,j,1) = TAUOC (i,j,1) + QTRCOC(k,i,j) * DELP(k,i,j) / GRAV * CEXTP(k,i,j,3,1)
          TAUOCA  (i,j,1) = TAUOCA(i,j,1) + QTRCOC(k,i,j) * DELP(k,i,j) / GRAV * CABSP(k,i,j,3,1)
          TAUBC   (i,j,1) = TAUBC (i,j,1) + QTRCBC(k,i,j) * DELP(k,i,j) / GRAV * CEXTBP(1)
          TAUBCA  (i,j,1) = TAUBCA(i,j,1) + QTRCBC(k,i,j) * DELP(k,i,j) / GRAV * CABSBP(1)
          CEXTOC(k,i,j,1) = QTRCOC(k,i,j) * DELP(k,i,j) / GRAV * CEXTA(k,i,j,3,1) / DELZ(k,i,j)
          CABSOC(k,i,j,1) = QTRCOC(k,i,j) * DELP(k,i,j) / GRAV * CABSA(k,i,j,3,1) / DELZ(k,i,j)
          CBAKOC(k,i,j,1) = QTRCOC(k,i,j) * DELP(k,i,j) / GRAV * CBAKA(k,i,j,3,1) / DELZ(k,i,j)
          CEXTBC(k,i,j,1) = QTRCBC(k,i,j) * DELP(k,i,j) / GRAV * CEXTBA(1)        / DELZ(k,i,j)
          CABSBC(k,i,j,1) = QTRCBC(k,i,j) * DELP(k,i,j) / GRAV * CABSBA(1)        / DELZ(k,i,j)
          CBAKBC(k,i,j,1) = QTRCBC(k,i,j) * DELP(k,i,j) / GRAV * CBAKBA(1)        / DELZ(k,i,j)
       enddo
    enddo
    enddo

    !clear sky diagnoses
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
          
       do k = KS, KE
          if ( TCLDF2(k,i,j) >= CCMAX ) then
             CCONTC(k,i,j)         = CONST_UNDEF !VMISS
             BCCNTC(k,i,j)         = CONST_UNDEF !VMISS
             CEXTOC(k,i,j,      2) = CONST_UNDEF !VMISS
             CABSOC(k,i,j,      2) = CONST_UNDEF !VMISS
             CBAKOC(k,i,j,      2) = CONST_UNDEF !VMISS
             CEXTBC(k,i,j,      2) = CONST_UNDEF !VMISS
             CABSBC(k,i,j,      2) = CONST_UNDEF !VMISS
             CBAKBC(k,i,j,      2) = CONST_UNDEF !VMISS
             CEXTCA(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CABSCA(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             CBAKCA(k,i,j,1:IWA,2) = CONST_UNDEF !VMISS
             NUMCNC(k,i,j,1:NCA1)  = CONST_UNDEF !VMISS
          else
             CCONTC(k,i,j)         = CCONT (k,i,j)
             BCCNTC(k,i,j)         = BCCONT(k,i,j)
             CEXTOC(k,i,j,      2) = CEXTOC(k,i,j,      1)
             CABSOC(k,i,j,      2) = CABSOC(k,i,j,      1)
             CBAKOC(k,i,j,      2) = CBAKOC(k,i,j,      1)
             CEXTBC(k,i,j,      2) = CEXTBC(k,i,j,      1)
             CABSBC(k,i,j,      2) = CABSBC(k,i,j,      1)
             CBAKBC(k,i,j,      2) = CBAKBC(k,i,j,      1)
             CEXTCA(k,i,j,1:IWA,2) = CEXTCA(k,i,j,1:IWA,1)
             CABSCA(k,i,j,1:IWA,2) = CABSCA(k,i,j,1:IWA,1)
             CBAKCA(k,i,j,1:IWA,2) = CBAKCA(k,i,j,1:IWA,1)
             NUMCNC(k,i,j,1:NCA1)  = NUMCON(k,i,j,1:NCA1)
          endif
       enddo

       if ( CCOVMR(i,j) >= CCMAX ) then
          COLMTC(i,j)         = CONST_UNDEF !VMISS
          COLBTC(i,j)         = CONST_UNDEF !VMISS
          TAUOC (i,j,      2) = CONST_UNDEF !VMISS
          TAUOCA(i,j,      2) = CONST_UNDEF !VMISS
          TAUBC (i,j,      2) = CONST_UNDEF !VMISS
          TAUBCA(i,j,      2) = CONST_UNDEF !VMISS
          TAUCA (i,j,1:IWA,2) = CONST_UNDEF !VMISS
          TAUCAA(i,j,1:IWA,2) = CONST_UNDEF !VMISS
       else
          COLMTC(i,j)         = COLMST(i,j)
          COLBTC(i,j)         = COLBCT(i,j)
          TAUOC (i,j,      2) = TAUOC (i,j,      1)
          TAUOCA(i,j,      2) = TAUOCA(i,j,      1)
          TAUBC (i,j,      2) = TAUBC (i,j,      1)
          TAUBCA(i,j,      2) = TAUBCA(i,j,      1) 
          TAUCA (i,j,1:IWA,2) = TAUCA (i,j,1:IWA,1)
          TAUCAA(i,j,1:IWA,2) = TAUCAA(i,j,1:IWA,1)
       endif
       
    enddo
    enddo
    !================================================
    !end optical parameters


#if defined(DEBUG) || defined(QUICKDEBUG)
    ! !checl OUTPUT
    ! !================================================
    ! check_2d_data(:,:, 1) = EMITBC(:,:)
    ! check_2d_data(:,:, 2) = EMITOC(:,:)
    ! check_2d_data(:,:, 3) = EMISSA(:,:)
    ! check_2d_data(:,:, 4) = CHESOA(:,:)
    ! check_2d_data(:,:, 5) = TERP(:,:)
    ! check_2d_data(:,:, 6) = ISOP(:,:)
    ! check_2d_data(:,:, 7) = EMIISO(:,:,1)
    ! check_2d_data(:,:, 8) = EMIISO(:,:,2)
    ! check_2d_data(:,:, 9) = WETFT(:,:)
    ! check_2d_data(:,:,10) = WETFT_CP(:,:)
    ! check_2d_data(:,:,11) = WETFT_MP_RAIN(:,:)
    ! check_2d_data(:,:,12) = WETFT_MP_SNOW(:,:)
    ! check_2d_data(:,:,13) = DRYFT(:,:)
    ! check_2d_data(:,:,14) = GRAVFT(:,:)
    ! check_2d_data(:,:,15) = DEPFT(:,:)
    ! check_2d_data(:,:,16) = WETFOC(:,:)
    ! check_2d_data(:,:,17) = WETFOC_CP(:,:)
    ! check_2d_data(:,:,18) = WETFOC_MP_RAIN(:,:)
    ! check_2d_data(:,:,19) = WETFOC_MP_SNOW(:,:)
    ! check_2d_data(:,:,20) = DRYFOC(:,:)
    ! check_2d_data(:,:,21) = GRVFOC(:,:)
    ! check_2d_data(:,:,22) = DEPFOC(:,:)
    ! check_2d_data(:,:,23) = WETFBC(:,:)
    ! check_2d_data(:,:,24) = WETFBC_CP(:,:)
    ! check_2d_data(:,:,25) = WETFBC_MP_RAIN(:,:)
    ! check_2d_data(:,:,26) = DRYFBC(:,:)
    ! check_2d_data(:,:,27) = GRVFBC(:,:)
    ! check_2d_data(:,:,28) = DEPFBC(:,:)
    ! check_2d_data(:,:,29) = COLMS (:,:,1)
    ! check_2d_data(:,:,30) = COLMS (:,:,2)
    ! check_2d_data(:,:,31) = COLMS (:,:,3)
    ! check_2d_data(:,:,32) = COLMST(:,:)
    ! check_2d_data(:,:,33) = COLBCT(:,:)
    ! check_2d_data(:,:,34) = COLOCT(:,:)
    ! check_2d_data(:,:,35) = TAUCA(:,:,1,1)
    ! check_2d_data(:,:,36) = TAUCAA(:,:,1,1)
    ! check_2d_data(:,:,37) = TAUOC(:,:,  1)
    ! check_2d_data(:,:,38) = TAUOCA(:,:,  1)
    ! check_2d_data(:,:,39) = TAUBC(:,:,  1)
    ! check_2d_data(:,:,40) = TAUBCA(:,:,  1)
    ! check_2d_data(:,:,41) = COLMTC(:,:)
    ! check_2d_data(:,:,42) = COLBTC(:,:)
    ! check_2d_data(:,:,43) = TAUCA  (:,:,1,2)
    ! check_2d_data(:,:,44) = TAUCAA (:,:,1,2)
    ! check_2d_data(:,:,45) = TAUOC  (:,:,  2)
    ! check_2d_data(:,:,46) = TAUOCA (:,:,  2)
    ! check_2d_data(:,:,47) = TAUBC  (:,:,  2)
    ! check_2d_data(:,:,48) = TAUBCA (:,:,  2)
    
    ! call STATISTICS_detail( IA, IS, IE, JA, JS, JE, 48, &
    !                         (/ 'EMITBC', 'EMITOC', 'EMISSA', 'CHESOA', 'TERP', &
    !                            'ISOP', 'EMIISO1', 'EMIISO2', 'WETFT', 'WETFT_CP', &
    !                            'WETFT_MP_RAIN', 'WETFT_MP_SNOW', 'DRYFT', 'GRAVFT', 'DEPFT', &
    !                            'WETFOC', 'WETFOC_CP', 'WETFOC_MP_RAIN', 'WETFOC_MP_SNOW', ' DRYFOC', &
    !                            'GRVFOC', 'DEPFOC', 'WETFBC', 'WETFBC_CP', 'WETFBC_MP_RAIN', &
    !                            'DRYFBC', 'GRVFBC', 'DEPFBC', 'COLMS1', 'COLMS2', &
    !                            'COLMS3', 'COLMST', 'COLBCT', 'COLOCT', 'TAUCA', &
    !                            'TAUCAA', 'TAUOC', 'TAUOCA', 'TAUBC', 'TAUBCA', &
    !                            'COLMTC', 'COLBTC', 'TAUCA12', 'TAUCAA12', 'TAUOC2', &
    !                            'TAUOCA2', 'TAUBC2', 'TAUBCA2' /), &
    !                            check_2d_data(:,:,1:48) )

    ! check_3d_data(:,:,:, 1) = QLOAD (:,:,:,1)
    ! check_3d_data(:,:,:, 2) = QLOAD (:,:,:,2)
    ! check_3d_data(:,:,:, 3) = QLOAD (:,:,:,3)
    ! check_3d_data(:,:,:, 4) = CCON  (:,:,:,1)
    ! check_3d_data(:,:,:, 5) = CCON  (:,:,:,2)
    ! check_3d_data(:,:,:, 6) = CCON  (:,:,:,3)
    ! check_3d_data(:,:,:, 7) = NUMCON(:,:,:,1)
    ! check_3d_data(:,:,:, 8) = NUMCON(:,:,:,2)
    ! check_3d_data(:,:,:, 9) = NUMCON(:,:,:,3)
    ! check_3d_data(:,:,:,10) = HCPPB (:,:,:,1)
    ! check_3d_data(:,:,:,11) = HCPPB (:,:,:,2)
    ! check_3d_data(:,:,:,12) = QLOADT(:,:,:)
    ! check_3d_data(:,:,:,13) = QTRCOC(:,:,:)
    ! check_3d_data(:,:,:,14) = QTRCBC(:,:,:)
    ! check_3d_data(:,:,:,15) = CCONT (:,:,:)
    ! check_3d_data(:,:,:,16) = BCCONT(:,:,:) 
    ! check_3d_data(:,:,:,17) = OCCONT(:,:,:)
    ! check_3d_data(:,:,:,18) = CEXTCA(:,:,:,1,1)
    ! check_3d_data(:,:,:,19) = CABSCA(:,:,:,1,1)
    ! check_3d_data(:,:,:,20) = CBAKCA(:,:,:,1,1)
    ! check_3d_data(:,:,:,21) = CEXTOC(:,:,:,  1)
    ! check_3d_data(:,:,:,22) = CABSOC(:,:,:,  1)
    ! check_3d_data(:,:,:,23) = CBAKOC(:,:,:,  1)
    ! check_3d_data(:,:,:,24) = CEXTBC(:,:,:,  1)
    ! check_3d_data(:,:,:,25) = CABSBC(:,:,:,  1)
    ! check_3d_data(:,:,:,26) = CBAKBC(:,:,:,  1)
    ! check_3d_data(:,:,:,27) = NUMCNC(:,:,:,1)
    ! check_3d_data(:,:,:,28) = NUMCNC(:,:,:,2)
    ! check_3d_data(:,:,:,29) = NUMCNC(:,:,:,3)
    ! check_3d_data(:,:,:,30) = CCONTC(:,:,:)
    ! check_3d_data(:,:,:,31) = BCCNTC(:,:,:)
    ! check_3d_data(:,:,:,32) = CEXTCA(:,:,:,1,2)
    ! check_3d_data(:,:,:,33) = CABSCA(:,:,:,1,2)
    ! check_3d_data(:,:,:,34) = CBAKCA(:,:,:,1,2)
    ! check_3d_data(:,:,:,35) = CEXTOC(:,:,:,  2)
    ! check_3d_data(:,:,:,36) = CABSOC(:,:,:,  2)
    ! check_3d_data(:,:,:,37) = CBAKOC(:,:,:,  2)
    ! check_3d_data(:,:,:,38) = CEXTBC(:,:,:,  2)
    ! check_3d_data(:,:,:,39) = CABSBC(:,:,:,  2)
    ! check_3d_data(:,:,:,40) = CBAKBC(:,:,:,  2)
    ! call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, 40, &
    !                         (/ 'QLOAD1', 'QLOAD2', 'QLOAD3', 'CCON1', 'CCON2', &
    !                            'CCON3', 'NUMCON1', 'NUMCON2', 'NUMCON3', 'HCPPB1', &
    !                            'HCPPB2', 'QLOADT', 'QTRCOC', 'QTRCBC', 'CCONT', &
    !                            'BCCONT', 'OCCONT', 'CEXTCA11', 'CABSCA11', 'CBAKCA11', &
    !                            'CEXTOC1', 'CABSOC1', 'CBAKOC1', 'CEXTBC1', 'CABSBC1', &
    !                            'CBAKBC1', 'NUMCNC1', 'NUMCNC2', 'NUMCNC3', 'CCONTC', &
    !                            'BCCNTC', 'CEXTCA12', 'CABSCA12', 'CBAKCA12', 'CEXTOC2', &
    !                            'CABSOC2', 'CBAKOC2', 'CEXTBC2', 'CABSBC2', 'CBAKBC2' /), &
    !                         check_3d_data(:,:,:,1:40) )
    ! !================================================
    ! !end check OUTPUT
#endif

    
    !OUTPUT
    !================================================
    !flux
    !-------------------------------
    !                                    i j
    call FILE_HISTORY_in( EMITBC        (:,:),   'EFBC_SPRINTARS',           'emission flux (BC)',                         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMITOC        (:,:),   'EFOC_SPRINTARS',           'emission flux (OM(OC))',                     'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMISSA        (:,:),   'EFSSA_SPRINTARS',          'emission flux (oceanic OM)',                 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( CHESOA        (:,:),   'CHESOA_SPRINTARS',         'secondary organic aerosol (SOA) production', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( TERP          (:,:),   'TERP_SPRINTARS',           'emission flux (terpene) ',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( ISOP          (:,:),   'ISOP_SPRINTARS',           'emission flux (isoprene)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMIISO        (:,:,1), 'EMIISO1_SPRINTARS',        'emission flux (terpene;oceanic precursor)',  'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( EMIISO        (:,:,2), 'EMIISO2_SPRINTARS',        'emission flux (isoprene;oceanic precursor)', 'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT         (:,:),   'WEFTCA_SPRINTARS',         'wet deposition flux (BC+OM(OC))',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_CP      (:,:),   'WEFTCA_CP_SPRINTARS',      'wet deposition flux (BC+OM(OC);CP)',         'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_RAIN (:,:),   'WEFTCA_MP_RAIN_SPRINTARS', 'wet deposition flux (BC+OM(OC);MP(RAIN))',   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFT_MP_SNOW (:,:),   'WEFTCA_MP_SNOW_SPRINTARS', 'wet deposition flux (BC+OM(OC);MP(SNOW))',   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYFT         (:,:),   'DRFTCA_SPRINTARS',         'dry deposition flux (BC+OM(OC))',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRAVFT        (:,:),   'GRFTCA_SPRINTARS',         'gravitational settling flux (BC+OM(OC))',    'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DEPFT         (:,:),   'TDFTCA_SPRINTARS',         'total deposition flux (BC+OM(OC))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFOC        (:,:),   'WEFTOC_SPRINTARS',         'wet deposition flux (OM(OC))',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFOC_CP     (:,:),   'WEFTOC_CP_SPRINTARS',      'wet deposition flux (OM(OC);CP)',            'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFOC_MP_RAIN(:,:),   'WEFTOC_MP_RAIN_SPRINTARS', 'wet deposition flux (OM(OC);MP(RAIN))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFOC_MP_SNOW(:,:),   'WEFTOC_MP_SNOW_SPRINTARS', 'wet deposition flux (OM(OC);MP(SNOW))',      'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYFOC        (:,:),   'DRFTOC_SPRINTARS',         'dry deposition flux (OM(OC))',               'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRVFOC        (:,:),   'GRFTOC_SPRINTARS',         'gravitational settling flux (OM(OC))',       'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DEPFOC        (:,:),   'TDFTOC_SPRINTARS',         'total deposition flux (OM(OC))',             'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFBC        (:,:),   'WEFTBC_SPRINTARS',         'wet deposition flux (BC)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFBC_CP     (:,:),   'WEFTBC_CP_SPRINTARS',      'wet deposition flux (BC;CP)',                'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFBC_MP_RAIN(:,:),   'WEFTBC_MP_RAIN_SPRINTARS', 'wet deposition flux (BC;MP(RAIN))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( WETFBC_MP_SNOW(:,:),   'WEFTBC_MP_SNOW_SPRINTARS', 'wet deposition flux (BC;MP(SNOW))',          'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DRYFBC        (:,:),   'DRFTBC_SPRINTARS',         'dry deposition flux (BC)',                   'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( GRVFBC        (:,:),   'GRFTBC_SPRINTARS',         'gravitational settling flux (BC)',           'kg/m2/s', dim_type = 'XY'  )
    call FILE_HISTORY_in( DEPFBC        (:,:),   'TDFTBC_SPRINTARS',         'total deposition flux (BC)',                 'kg/m2/s', dim_type = 'XY'  )
    
    !aerosol
    !-------------------------------
    !                            k i j q
    call FILE_HISTORY_in( QLOAD (:,:,:,1), 'QLDC01_SPRINTARS', 'mixing ratio (carbon;internal mixed OM(OC) & BC)',         'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QLOAD (:,:,:,2), 'QLDC02_SPRINTARS', 'mixing ratio (carbon;SOC from NVOC)',                      'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QLOAD (:,:,:,3), 'QLDC03_SPRINTARS', 'mixing ratio (carbon;external mixed BC)',                  'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CCON  (:,:,:,1), 'CONC01_SPRINTARS', 'mass concentration (carbon;internal mixed OM(OC) & BC)',   'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CCON  (:,:,:,2), 'CONC02_SPRINTARS', 'mass concentration (carbon;SOC from NVOC)',                'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CCON  (:,:,:,3), 'CONC03_SPRINTARS', 'mass concentration (carbon;external mixed BC)',            'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCON(:,:,:,1), 'NMCC01_SPRINTARS', 'number concentration (carbon;internal mixed OM(OC) & BC)', '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCON(:,:,:,2), 'NMCC02_SPRINTARS', 'number concentration (carbon;SOC from NVOC)',              '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( NUMCON(:,:,:,3), 'NMCC03_SPRINTARS', 'number concentration (carbon;external mixed BC)',          '1/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMS   (:,:,1), 'COLC01_SPRINTARS', 'column loading (carbon;internal mixed OM(OC) & BC)',       'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( COLMS   (:,:,2), 'COLC02_SPRINTARS', 'column loading (carbon;SOC from NVOC)',                    'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( COLMS   (:,:,3), 'COLC03_SPRINTARS', 'column loading (carbon;external mixed BC)',                'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( HCPPB (:,:,:,1), 'HCPPB1_SPRINTARS', 'terpene',                                                  'ppbv',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( HCPPB (:,:,:,2), 'HCPPB2_SPRINTARS', 'isoprene',                                                 'ppbv',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( QLOADT(:,:,:),   'QLDCT_SPRINTARS',  'mass mixing ratio (BC+OM(OC))',                            'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QTRCOC(:,:,:),   'QTRCOC_SPRINTARS', 'mass mixing ratio (OM(OC))',                               'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QTRCBC(:,:,:),   'QTRCBC_SPRINTARS', 'mass mixing ratio (BC)',                                   'kg/kg', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CCONT (:,:,:),   'CONCT_SPRINTARS',  'mass concentration (BC & OM(OC))',                         'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( BCCONT(:,:,:),   'CONBT_SPRINTARS',  'mass concentration (BC)',                                  'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( OCCONT(:,:,:),   'CONOT_SPRINTARS',  'mass concentration (OM(OC))',                              'kg/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMST  (:,:),   'COLCT_SPRINTARS',  'column loading (BC+OM(OC))',                               'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( COLBCT  (:,:),   'COLBT_SPRINTARS',  'column loading (BC)',                                      'kg/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( COLOCT  (:,:),   'COLOT_SPRINTARS',  'column loading (OM(OC))',                                  'kg/m2', dim_type = 'XY'  )
    
    !optical properties
    !-------------------------------
    !                              i j w 2
    call FILE_HISTORY_in( TAUCA   (:,:,1,1), 'TAUCA_SPRINTARS',   'AOT (BC+OM(OC);550nm;allsky)',         '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUCAA  (:,:,1,1), 'TAUCAA_SPRINTARS',  'AOT (BC+OM(OC);550nm;allsky;abs)',     '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUOC   (:,:,  1), 'TAUOC_SPRINTARS',   'AOT (OM(OC);550nm;allsky)',            '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUOCA  (:,:,  1), 'TAUOCA_SPRINTARS',  'AOT (OM(OC);550nm;allsky;abs)',        '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUBC   (:,:,  1), 'TAUBC_SPRINTARS',   'AOT (BC;550nm;allsky)',                '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUBCA  (:,:,  1), 'TAUBCA_SPRINTARS',  'AOT (BC;550nm;allsky;abs)',            '1',      dim_type = 'XY'  )
!""    call FILE_HISTORY_in( CEXTCA(:,:,:,1,1), 'EXTCA1_SPRINTARS',  'extinction  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CABSCA(:,:,:,1,1), 'ABSCA1_SPRINTARS',  'absorption  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CBAKCA(:,:,:,1,1), 'BAKCA1_SPRINTARS',  'backscatter (BC+OM(OC);532nm;allsky)', '1/m/sr', dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CEXTOC(:,:,:,  1), 'EXTOC1_SPRINTARS',  'extinction  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CABSOC(:,:,:,  1), 'ABSOC1_SPRINTARS',  'absorption  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CBAKOC(:,:,:,  1), 'BAKOC1_SPRINTARS',  'backscatter (OM(OC);532nm;allsky)',    '1/m/sr', dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CEXTBC(:,:,:,  1), 'EXTBC1_SPRINTARS',  'extinction  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CABSBC(:,:,:,  1), 'ABSBC1_SPRINTARS',  'absorption  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( CBAKBC(:,:,:,  1), 'BAKBC1_SPRINTARS',  'backscatter (BC;532nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )

    call FILE_HISTORY_in( QWTDEP_CP(:,:,:,1), 'EXTCA1_SPRINTARS',  'extinction  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_CP(:,:,:,2), 'ABSCA1_SPRINTARS',  'absorption  (BC+OM(OC);532nm;allsky)', '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_CP(:,:,:,3), 'BAKCA1_SPRINTARS',  'backscatter (BC+OM(OC);532nm;allsky)', '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_RAIN(:,:,:,1), 'EXTOC1_SPRINTARS',  'extinction  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_RAIN(:,:,:,2), 'ABSOC1_SPRINTARS',  'absorption  (OM(OC);532nm;allsky)',    '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_RAIN(:,:,:,3), 'BAKOC1_SPRINTARS',  'backscatter (OM(OC);532nm;allsky)',    '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_SNOW(:,:,:,1), 'EXTBC1_SPRINTARS',  'extinction  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_SNOW(:,:,:,2), 'ABSBC1_SPRINTARS',  'absorption  (BC;532nm;allsky)',        '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( QWTDEP_MP_SNOW(:,:,:,3), 'BAKBC1_SPRINTARS',  'backscatter (BC;532nm;allsky)',        '1/m/sr', dim_type = 'ZXY' )

    !clear sky diagnoses
    !-------------------------------
    !                            k i j q
!""    call FILE_HISTORY_in( NUMCNC(:,:,:,1), 'NMCCC01_SPRINTARS', 'number concentration (carbon;internal mixed OM(OC) & BC;clearsky)', '1/m3', dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( NUMCNC(:,:,:,2), 'NMCCC02_SPRINTARS', 'number concentration (carbon;SOC from NVOC;clearsky)',              '1/m3', dim_type = 'ZXY' )
!""    call FILE_HISTORY_in( NUMCNC(:,:,:,3), 'NMCCC03_SPRINTARS', 'number concentration (carbon;external mixed BC;clearsky)',          '1/m3', dim_type = 'ZXY' )
    !                            k i j w 2
    call FILE_HISTORY_in( CCONTC(:,:,:),     'CONCTC_SPRINTARS',  'mass concentration (BC+OM(OC);clearsky)', 'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( BCCNTC(:,:,:),     'CONBTC_SPRINTARS',  'mass concentration (BC);clearsky)',       'kg/m3',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( COLMTC  (:,:),     'COLCTC_SPRINTARS',  'column loading (BC+OM(OC);clearsky)',     'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in( COLBTC  (:,:),     'COLBTC_SPRINTARS',  'column loading (BC;clearsky)',            'kg/m2',  dim_type = 'XY'  )
    call FILE_HISTORY_in(  TAUCA  (:,:,1,2), 'TAUCAC_SPRINTARS',  'AOT (BC+OM(OC);550nm;clearsky)',          '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUCAA  (:,:,1,2), 'TAUCACA_SPRINTARS', 'AOT (BC+OM(OC);550nm;clearsky;abs)',      '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in(  TAUOC  (:,:,  2), 'TAUOCC_SPRINTARS',  'AOT (OM(OC);550nm;clearsky)',             '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUOCA  (:,:,  2), 'TAUOCCA_SPRINTARS', 'AOT (OM(OC);550nm;clearsky;abs)',         '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in(  TAUBC  (:,:,  2), 'TAUBCC_SPRINTARS',  'AOT (BC;550nm;clearsky)',                 '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( TAUBCA  (:,:,  2), 'TAUBCCA_SPRINTARS', 'AOT (BC;550nm;clearsky;abs)',             '1',      dim_type = 'XY'  )
    call FILE_HISTORY_in( CEXTCA(:,:,:,1,2), 'EXTCAC1_SPRINTARS', 'extinction  (BC+OM(OC);532nm;clearsky)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABSCA(:,:,:,1,2), 'ABSCAC1_SPRINTARS', 'absorption  (BC+OM(OC);532nm;clearsky)',  '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKCA(:,:,:,1,2), 'BAKCAC1_SPRINTARS', 'backscatter (BC+OM(OC);532nm;clearsky)',  '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTOC(:,:,:,  2), 'EXTOCC1_SPRINTARS', 'extinction  (OM(OC);532nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABSOC(:,:,:,  2), 'ABSOCC1_SPRINTARS', 'absorption  (OM(OC);532nm;clearsky)',     '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKOC(:,:,:,  2), 'BAKOCC1_SPRINTARS', 'backscatter (OM(OC);532nm;clearsky)',     '1/m/sr', dim_type = 'ZXY' )
    call FILE_HISTORY_in( CEXTBC(:,:,:,  2), 'EXTBCC1_SPRINTARS', 'extinction  (BC;532nm;clearsky)',         '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CABSBC(:,:,:,  2), 'ABSBCC1_SPRINTARS', 'absorption  (BC;532nm;clearsky)',         '1/m',    dim_type = 'ZXY' )
    call FILE_HISTORY_in( CBAKBC(:,:,:,  2), 'BAKBCC1_SPRINTARS', 'backscatter (BC;532nm;clearsky)',         '1/m/sr', dim_type = 'ZXY' )
    !================================================
    !end OUTPUT

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROCA
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_aerocarb
