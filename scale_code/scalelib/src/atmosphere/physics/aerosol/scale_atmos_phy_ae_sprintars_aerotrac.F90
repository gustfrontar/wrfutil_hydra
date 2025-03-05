!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - tracer
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_aerotrac
  !-----------------------------------------------------------------------------
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
  use scale_atmos_grid_cartesC_real, only: &
       HAREA => ATMOS_GRID_CARTESC_REAL_AREA !< horizontal area ( xy, normal z) [m2]

  use scale_const, only: &
       PI   => CONST_PI,   & !< pi
       GRAV => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       CONST_D2R
  
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
       NAT,                                           & !< total number of air tracers
       IWA,                                           & !< number of wavelengths
       NRH,                                           & !< number of RH bin
       BRH,                                           & !< RH bin
       NU,                                            & !< air viscosity [kg/m/s] = [Pa*s]
       RADR,                                          & !< rain droplet radius      [m]
       VTR,                                           & !< rain terminal velocity   [m/s]
       DENSR,                                         & !< rain droplet density     [kg/m3]
       DENSW,                                         & !< density of water         [kg/m3]
       THRESP                                           !< precipitation threshold 
       
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT_setup
  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT
!  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY
  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY2
!  public :: ATMOS_PHY_AE_SPRINTARS_AEROAT_SILICATE_PROXY
  
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
  real(RP), private, save :: FINCAT = 0.05E+0_RP !< incloud coefficient
  
  real(RP), private, allocatable, save :: QTRCINP(:,:,:,:) !< mixing ratio in prec(MP)

  !for ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY
  integer,  private, save :: I_FDNPP = -1
  integer,  private, save :: J_FDNPP = -1
  integer,  private, save :: myrank_FDNPP = -1
  integer,  private, save :: emit_time (66,2)
  real(RP), private, save :: emit_data (66) !< Katata et al. (2015) [Bq/h]
  integer,  private, save :: emit_kbidx(66) = -1
  integer,  private, save :: emit_ktidx(66) = -1
  real(RP), private, allocatable, save :: emit_coeff(:,:) !(66,KA)
  real(RP), private, allocatable, save :: ACCUM_TOT_DEP        (:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_WET_DEP        (:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_WET_DEP_CP     (:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_WET_DEP_MP_RAIN(:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_WET_DEP_MP_SNOW(:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_DRY_DEP        (:,:,:) ![kg/m2]
  real(RP), private, allocatable, save :: ACCUM_GRV_DEP        (:,:,:) ![kg/m2]

  real(RP), private, allocatable :: zero_2d(:,:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Tracer setup
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_setup( LON_REAL, LAT_REAL, GDZM, fname_aeroat, QA_AE_AT )
    implicit none
    real(RP), intent(in) :: LON_REAL     (IA,JA) !< longitude
    real(RP), intent(in) :: LAT_REAL     (IA,JA) !< latitude
    real(RP), intent(in) :: GDZM    (0:KA,IA,JA) !< GDZM(0:KA,IA,JA) : altitude (half level) [m]
    character(len=H_LONG), intent(in) :: fname_aeroat !< file for emission for aero trc
    integer,  intent(in) :: QA_AE_AT
    integer  :: fi, ierr
    integer  :: n
    integer  :: sdate(66,6), edate(66,6)
    real(RP) :: emit_hgt(66,2)
    real(RP) :: mean_hgt
    integer  :: k, i, j, iq
    character :: cnumber*2
    !----------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_aerotrac",*) 'Setup'
    
    allocate( QTRCINP(KA,IA,JA,NAT) )
    QTRCINP(:,:,:,:) = 0.0_RP

    !for ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY
    !-------------------------
    do j = JS, JE-1
    do i = IS, IE-1
       if (        37.42_RP * CONST_D2R >= (LAT_REAL(i,j)+LAT_REAL(i,  j-1))/2.0_RP .and.  37.42_RP * CONST_D2R < (LAT_REAL(i,j)+LAT_REAL(i,  j+1))/2.0_RP &
            .and. 141.03_RP * CONST_D2R >= (LON_REAL(i,j)+LON_REAL(i-1,j  ))/2.0_RP .and. 141.03_RP * CONST_D2R < (LON_REAL(i,j)+LON_REAL(i+1,j  ))/2.0_RP ) then
          I_FDNPP = i
          J_FDNPP = j
          myrank_FDNPP = PRC_myrank
       endif
    enddo
    enddo


    allocate( emit_coeff(66,KA) )
    emit_coeff(:,:) = 0.0_RP
    
    fi = io_get_available_fid()
    n = 0
    open( fi, file = trim(fname_aeroat), status = 'old', iostat = ierr )
    if ( ierr /= 0 ) then
       write(*,*) 'Cannot open 137Cs emission profile. check!'
       call PRC_abort
    endif
    do
       n = n + 1
       read( fi, *, iostat = ierr ) sdate(n,1:6), edate(n,1:6), emit_data(n), emit_hgt(n,1:2)
       if ( ierr /= 0 .and. ierr /= -1 ) then
          write(*,*) 'Cannot read 137Cs emission profile. check!'
          call PRC_abort
       endif
       
       emit_time(n,1) = sdate(n,3) * 1000000 + sdate(n,4) * 10000 + sdate(n,5) * 100 + sdate(n,6)
       emit_time(n,2) = edate(n,3) * 1000000 + edate(n,4) * 10000 + edate(n,5) * 100 + edate(n,6)


       if ( myrank_FDNPP == -1 ) exit
       do k = KS-1, KE-1
          !emit_kbidx
          if ( emit_hgt(n,1) >= GDZM(k,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP) .and. emit_hgt(n,1) < GDZM(k+1,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP) ) then
             emit_kbidx(n) = k+1
          endif
          !emit_ktidx
          if ( emit_hgt(n,2) >= GDZM(k,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP) .and. emit_hgt(n,2) < GDZM(k+1,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP) ) then
             emit_ktidx(n) = k+1
          endif
       enddo

       do k = KS, KE
          if     ( k == emit_kbidx(n) .and. k == emit_ktidx(n) ) then
             emit_coeff(n,k) = 1.0_RP
             exit
          elseif ( k == emit_kbidx(n) ) then
             emit_coeff(n,k) = max( ( (GDZM(k,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP)) - emit_hgt(n,1)                                          ) / ( emit_hgt(n,2)-emit_hgt(n,1) ), 0.0_RP )
          elseif ( k >  emit_kbidx(n) .and. k <  emit_ktidx(n) ) then
             emit_coeff(n,k) = max( ( (GDZM(k,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP)) - (GDZM(k-1,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP)) ) / ( emit_hgt(n,2)-emit_hgt(n,1) ), 0.0_RP )
          elseif ( k == emit_ktidx(n) ) then
             emit_coeff(n,k) = max( ( emit_hgt(n,2)                                        - (GDZM(k-1,I_FDNPP,J_FDNPP)-GDZM(KS-1,I_FDNPP,J_FDNPP)) ) / ( emit_hgt(n,2)-emit_hgt(n,1) ), 0.0_RP )
          endif
       enddo
             
       if ( emit_kbidx(n) == -1 .or. emit_ktidx(n) == -1 ) then
          ! write(*,*) mean_hgt
          do k = KS-1, KE
             write(*,*) k, GDZM(k,IS,JS)
          enddo
          LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_AEROAT_setup",*) 'emit_kidx is -1. check!'
          call PRC_abort
       endif
       
       write(*,*) 'time:', emit_time(n,1), emit_time(n,2)
       write(*,'(e12.5)') emit_data(n) * 1.0E-14_RP
       write(*,*) emit_kbidx(n), emit_ktidx(n)
       do k = emit_kbidx(n), emit_ktidx(n)
          write(*,'(i2,2x,f5.3,2x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5)') k, emit_coeff(n,k), emit_hgt(n,1), emit_hgt(n,2), GDZM(k-1,I_FDNPP,J_FDNPP), GDZM(k,I_FDNPP,J_FDNPP)
       enddo
       write(*,*) ''
       if ( n == 66 ) then
          exit
       endif
       
    enddo
    close( fi )

    allocate( ACCUM_TOT_DEP        (IA,JA,NAT) )
    allocate( ACCUM_WET_DEP        (IA,JA,NAT) )
    allocate( ACCUM_WET_DEP_CP     (IA,JA,NAT) )
    allocate( ACCUM_WET_DEP_MP_RAIN(IA,JA,NAT) )
    allocate( ACCUM_WET_DEP_MP_SNOW(IA,JA,NAT) )
    allocate( ACCUM_DRY_DEP        (IA,JA,NAT) )
    allocate( ACCUM_GRV_DEP        (IA,JA,NAT) )
    ACCUM_TOT_DEP        (:,:,:) = 0.0_RP
    ACCUM_WET_DEP        (:,:,:) = 0.0_RP 
    ACCUM_WET_DEP_CP     (:,:,:) = 0.0_RP
    ACCUM_WET_DEP_MP_RAIN(:,:,:) = 0.0_RP
    ACCUM_WET_DEP_MP_SNOW(:,:,:) = 0.0_RP 
    ACCUM_DRY_DEP        (:,:,:) = 0.0_RP
    ACCUM_GRV_DEP        (:,:,:) = 0.0_RP
    !-------------------------

    allocate( zero_2d(IA,JA) )
    zero_2d(:,:) = 0.0_RP
    do iq = 1, QA_AE_AT
       write(cnumber,'(i2.2)') iq
       call FILE_HISTORY_in( zero_2d(:,:), 'EFTAT'//cnumber//'_SPRINTARS',          'emission flux (tracer'//cnumber//')',                'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'TDFTAT'//cnumber//'_SPRINTARS',         'total deposition flux (tracer'//cnumber//')',        'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'WEFTAT'//cnumber//'_SPRINTARS',         'wet deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'WEFTAT'//cnumber//'_CP_SPRINTARS',      'wet deposition flux (tracer'//cnumber//';CP)',       'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'WEFTAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(RAIN))', 'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'WEFTAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(SNOW))', 'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'DRFTAT'//cnumber//'_SPRINTARS',         'dry deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'GRFTAT'//cnumber//'_SPRINTARS',         'gravitational settling flux (tracer'//cnumber//')',  'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACTDAT'//cnumber//'_SPRINTARS',         'total deposition (tracer'//cnumber//')',                  'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACWDAT'//cnumber//'_SPRINTARS',         'wet deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACWDAT'//cnumber//'_CP_SPRINTARS',      'wet deposition (tracer'//cnumber//';CP)',                 'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACWDAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(RAIN))',           'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACWDAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(SNOW))',           'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACDDAT'//cnumber//'_SPRINTARS',         'dry deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( zero_2d(:,:), 'ACGDAT'//cnumber//'_SPRINTARS',         'gravitational settling deposition (tracer'//cnumber//')', 'kg/m2', dim_type = 'XY' )
    enddo

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_setup
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> SPRINTARS Tracer finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_finalize
    implicit none

    deallocate( QTRCINP )
    deallocate( emit_coeff )
    deallocate( ACCUM_TOT_DEP         )
    deallocate( ACCUM_WET_DEP         )
    deallocate( ACCUM_WET_DEP_CP      )
    deallocate( ACCUM_WET_DEP_MP_RAIN )
    deallocate( ACCUM_WET_DEP_MP_SNOW )
    deallocate( ACCUM_DRY_DEP         )
    deallocate( ACCUM_GRV_DEP         )

    deallocate( zero_2d )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_finalize
  
  !-----------------------------------------------------------------------------
  !> SPRINTARS Tracer
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT &
       ( QTRC,    & ! [MODIFIED]
         dt_AE,   & ! [IN]
         QA_AE_AT ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE
    integer,  intent(in)    :: QA_AE_AT

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_AT)
    
    ! output
    !-------------------------

    ! internal work
    !-------------------------
    real(RP) :: decay_ratio_tracer01

    ! parameter
    !-------------------------
    real(RP), parameter :: half_life_tracer01 = 240.0_RP * 3600.0_RP ! half life time [s]

    ! others 
    !-------------------------
    integer :: k, i, j, iq
    
    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - tracer'

    !setup parameter 
    !================================================
    !Decay based on half life
    !-------------------------------
    decay_ratio_tracer01 = log(2.0_RP) / half_life_tracer01
    !================================================
    !end setup parameter

    
    !aerosol transport 
    !================================================
    !Decay
    !-------------------------------
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) - QTRC(k,i,j,iq) * decay_ratio_tracer01 * real(dt_AE,kind=RP)
          enddo
       enddo
       enddo
    enddo
    !================================================
    !end aerosol transport

    
    !aerosol
    !================================================
    !================================================
    !end aerosol

    
    !flux 
    !================================================
    !================================================
    !end flux

    
    !optical parameters
    !================================================
    !================================================
    !end optical parameters

    
    !OUTPUT
    !================================================
    !================================================
    !end OUTPUT

    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  !> SPRINTARS Tracer (sulfate proxy)
!  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY &
!       ( QTRC,           & ! [MODIFIED]
!         CDVE,           & ! [IN]
!         RH,             & ! [IN]
!         DENS,           & ! [IN]
!         DENS1,          & ! [IN]
!         FRAIN_MP_RAIN1, & ! [IN]
!         FRAIN_MP_SNOW1, & ! [IN]
!         FRAIN,          & ! [IN]
!         FRAIN_CP,       & ! [IN]
!         FRAIN_MP_RAIN,  & ! [IN]
!         FRAIN_MP_SNOW,  & ! [IN]
!         VTR_MP_RAIN,    & ! [IN]
!         VTR_MP_SNOW,    & ! [IN]
!         RADR_MP_RAIN,   & ! [IN]
!         RADR_MP_SNOW,   & ! [IN]
!         RAINN_CP,       & ! [IN]
!         RAINN_MP_RAIN,  & ! [IN]
!         RAINN_MP_SNOW,  & ! [IN]
!         TCLDF,          & ! [IN]
!         KUP,            & ! [IN]
!         KUPMAX,         & ! [IN]
!         DPKUP,          & ! [IN]
!         GDPM,           & ! [IN]
!         DELP,           & ! [IN]
!         GDZM,           & ! [IN]
!         DELZ,           & ! [IN]
!         RNWR,           & ! [IN]
!         RNWR_CP,        & ! [IN]
!         RNWR_MP_RAIN,   & ! [IN]
!         RNWR_MP_SNOW,   & ! [IN]
!         LON_REAL,       & ! [IN]
!         LAT_REAL,       & ! [IN]
!         ACTIAT,         & ! [IN]
!         SW_CCN,         & ! [IN]
!         TIME_NOWDATE,   & ! [IN]
!         dt_AE,          & ! [IN]
!         QA_AE_AT        ) ! [IN]
!    implicit none
!
!    ! input
!    !-------------------------
!    real(DP), intent(in)    :: dt_AE
!    integer,  intent(in)    :: QA_AE_AT
!    
!    real(RP), intent(in)    :: CDVE               (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
!    real(RP), intent(in)    :: RH            (KA,  IA,JA) !< relative humidity
!    real(RP), intent(in)    :: DENS          (KA,  IA,JA) !< air density
!    real(RP), intent(in)    :: DENS1         (KA,  IA,JA) !< air density   on 1 step before
!    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA) !< precipitation flux (liq    ) : MP on 1 step before
!    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA) !< precipitation flux (    ice) : MP on 1 step before
!    real(RP), intent(in)    :: FRAIN         (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP + MP
!    real(RP), intent(in)    :: FRAIN_CP      (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP
!    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA) !< precipitation flux (liq    ) : MP
!    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA) !< precipitation flux (    ice) : MP
!    real(RP), intent(in)    :: VTR_MP_RAIN   (KA,  IA,JA) !< terminal velocity of rain (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: VTR_MP_SNOW   (KA,  IA,JA) !< terminal velocity of snow (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: RADR_MP_RAIN  (KA,  IA,JA) !< effective radius of rain (microphysic scheme) [m]
!    real(RP), intent(in)    :: RADR_MP_SNOW  (KA,  IA,JA) !< effective radius of snow (microphysic scheme) [m]
!    real(RP), intent(in)    :: RAINN_CP      (KA,  IA,JA) !< rain number concentration (CP)      [1/m3]
!    real(RP), intent(in)    :: RAINN_MP_RAIN (KA,  IA,JA) !< rain number concentration (MP;RAIN) [1/m3]
!    real(RP), intent(in)    :: RAINN_MP_SNOW (KA,  IA,JA) !< rain number concentration (MP;SNOW) [1/m3]
!    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA) !< cloud fraction (CP+MP)
!    integer,  intent(in)    :: KUP                (IA,JA) !< uplift level
!    integer,  intent(in)    :: KUPMAX                     !< maximum uplift level
!    real(RP), intent(in)    :: DPKUP              (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
!    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA) !< pressure at half levels
!    real(RP), intent(in)    :: DELP          (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
!    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA) !< altitude at half levels
!    real(RP), intent(in)    :: DELZ          (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
!    real(RP), intent(in)    :: RNWR          (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: LON_REAL           (IA,JA) !< longitude
!    real(RP), intent(in)    :: LAT_REAL           (IA,JA) !< latitude
!    real(RP), intent(in)    :: ACTIAT        (KA,  IA,JA) !< ratio of activated aerosols [0-1]
!    logical,  intent(in)    :: SW_CCN                     !< switch for CCN scheme (true->use ACTI)
!    integer,  intent(in)    :: TIME_NOWDATE(6)            !< current time
!    
!    ! modified
!    !-------------------------
!    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_AT)
!    
!    ! output
!    !-------------------------
!
!    ! internal work
!    !-------------------------
!    !set parameter
!    real(RP) :: decay_ratio(QA_AE_AT)
!
!    real(RP) :: MAT           (QA_AE_AT)
!    real(RP) :: RADAT(KA,IA,JA,QA_AE_AT) !< tracer particle radius      [m]
!    real(RP) :: BRH1, BRH2, BRH3         !<
!    real(RP) :: VTER (KA,IA,JA,QA_AE_AT) !< terminal velocity of tracer [m/s]
!
!    !emission
!    real(RP) :: QEMIT(KA,IA,JA,QA_AE_AT) !emission mixing ratio [kg/kg]
!    real(RP) :: EMITF   (IA,JA,QA_AE_AT) !tracer emission flux  [kg/m2/s]
!    
!    !incloud
!    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_AT)  !< outside cloud
!    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_AT)  !< inside  cloud
!    
!    !subcloud scavenging (wash out)
!    ! real(RP) :: RAINN        (KA,IA,JA)                     !< rain number concentrarion        [1/m3]
!    ! real(RP) :: RAINN_CP     (KA,IA,JA)               
!    ! real(RP) :: RAINN_MP_RAIN(KA,IA,JA)              
!    ! real(RP) :: RAINN_MP_SNOW(KA,IA,JA)               
!    real(RP) :: TRACN                             !< number density of tracer [1/m3]
!    real(RP) :: TRACRN        (KA,IA,JA,QA_AE_AT) !< number density of tracer entering into raindrops [1/s/m3] <=> collision frequency between tracer and rain
!    real(RP) :: TRACRN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_AT) !< mixing ratio by wet deposition   [kg/kg] 
!    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_AT) !< tracers that do not cllide with rain particles [kg/kg]
!    
!    !incloud scavenging (rain out)
!    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_AT) !< mixing ratio moving cloud to rain particles [kg/kg]
!    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    
!    !re-emission from rain
!    real(RP) :: WETF           (IA,JA,QA_AE_AT) !< wet deposition flux
!    real(RP) :: WETF_CP        (IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_RAIN   (IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_SNOW   (IA,JA,QA_AE_AT)
!    real(RP) :: QTRCINR     (KA,IA,JA,QA_AE_AT) !< mixing ratio in rain(MP)
!    real(RP) :: QTRCINS     (KA,IA,JA,QA_AE_AT) !< mixing ratio in snow(MP)
!
!    !dry deposition
!    real(RP) :: CDVES   (IA,JA)           !< 
!    real(RP) :: QMTX (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
!    real(RP) :: DRYF    (IA,JA,QA_AE_AT)  !< dry deposition flux
!    
!    !gravitational settling
!    real(RP) :: QLOAD(KA,IA,JA,QA_AE_AT)  !< mixing ratio for radiation [kg/kg]
!    real(RP) :: GRAVF   (IA,JA,QA_AE_AT)  !< gravitational settling flux 
!    
!    !flux
!    real(RP) :: TRACDP  (IA,JA,QA_AE_AT)  !< sulfate total deposition
!    
!    ! parameter
!    !-------------------------
!    real(RP), parameter :: hltime(NAT) = (/ 30.2_RP * 365.0_RP * 86400.0_RP /) ! half life time [s] : 137Cs (30.1yrs)
!    real(RP), parameter :: DDV(NAT) = (/ 2.0E-3_RP /)                !< deposition velocity
!    real(RP), parameter :: &
!         RADTRAC(NRH)      &    !< sulfate particle radius
!         !    RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
!         = (/ 0.695E-7_RP, 0.850E-7_RP, 0.950E-7_RP, 0.103E-6_RP, 0.122E-6_RP, 0.157E-6_RP, 0.195E-6_RP, 0.231E-6_RP /)
!    real(RP), parameter :: DENSAT(NAT) = (/ 1.769E+3_RP /) !< tracer particle density [kg/m3]
!    real(RP), parameter :: AA    (NAT) = (/ 1.000E-4_RP /) !< collision efficiency
!    real(RP), parameter :: SIGMA (NAT) = (/ 1.526E+0_RP /) !< GSD of log-normal distribution
!    
!    ! others 
!    !-------------------------
!    integer   :: k, i, j, iq
!    integer   :: IRH
!    character :: cnumber*2
!    integer   :: time           ! ddhhmmss
!    integer   :: ievent, mevent
!    logical, save :: OFIRST = .true. 
!
!    !----------------------------------------------------
!    
!    LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - tracer'
!
!    !once 
!    !================================================
!    if ( OFIRST ) then
!       OFIRST = .false.
!       if ( PRC_myrank == myrank_FDNPP ) then
!          write(*,*) 'FDNPP PRC_myrank =', PRC_myrank
!          write(*,*) 'FDNPP =', I_FDNPP, J_FDNPP, LON_REAL(I_FDNPP,J_FDNPP)/CONST_D2R, LAT_REAL(I_FDNPP,J_FDNPP)/CONST_D2R
!       endif
!    endif
!
!    !================================================
!    !end once 
!    
!    !setup parameter 
!    !================================================
!    !MAT
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       MAT(iq) = 4.0_RP / 3.0_RP * PI * RADTRAC(1)**3 * DENSAT(iq) * exp( 4.5_RP * log(SIGMA(iq))**2 )
!    enddo
!    
!    !RADAT
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!          
!             if (                               RH(k,i,j) < BRH(2) ) then ! RH =  0 - 50 %
!                IRH = 1
!             elseif ( RH(k,i,j) >= BRH(2) .and. RH(k,i,j) < BRH(3) ) then ! RH = 50 - 70 %
!                IRH = 2
!             elseif ( RH(k,i,j) >= BRH(3) .and. RH(k,i,j) < BRH(4) ) then ! RH = 70 - 80 %
!                IRH = 3
!             elseif ( RH(k,i,j) >= BRH(4) .and. RH(k,i,j) < BRH(5) ) then ! RH = 80 - 90 %
!                IRH = 4
!             elseif ( RH(k,i,j) >= BRH(5) .and. RH(k,i,j) < BRH(6) ) then ! RH = 90 - 95 %
!                IRH = 5
!             elseif ( RH(k,i,j) >= BRH(6) .and. RH(k,i,j) < BRH(7) ) then ! RH = 95 - 98 %
!                IRH = 6
!             elseif ( RH(k,i,j) >= BRH(7) .and. RH(k,i,j) < BRH(8) ) then ! RH = 98 - 99 %
!                IRH = 7
!             else                                                         ! RH = 99 -    %
!                IRH = 8           
!             endif
!
!             if ( IRH == 8 ) then
!                RADAT(k,i,j,iq) = RADTRAC(IRH)
!             else
!                BRH1 =  RH(k,i,j) - BRH(IRH)
!                BRH2 = BRH(IRH+1) -  RH(k,i,j)
!                BRH3 = BRH(IRH+1) - BRH(IRH)
!                RADAT(k,i,j,iq) = ( BRH1 * RADTRAC(IRH+1) + BRH2 * RADTRAC(IRH) ) / BRH3
!             endif
!
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !VTER
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             VTER(k,i,j,iq) = 2.0_RP * (              ( RADTRAC(1)/RADAT(k,i,j,iq) )**3   * DENSAT(iq) &
!                                         + ( 1.0_RP - ( RADTRAC(1)/RADAT(k,i,j,iq) )**3 ) * DENSW    ) &
!                                     * RADAT(k,i,j,iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(iq))**2 )
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !Decay based on half life
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       decay_ratio(iq) = log(2.0_RP) / hltime(iq)
!    enddo
!    
!    !================================================
!    !end setup parameter
!
!
!    !aerosol transport 
!    !================================================
!      
!    !in rain(MP), in snow(MP)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
!            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
!              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
!              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
!              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
!              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
!    enddo
!
!    
!    !emission
!    !-------------------------------
!    time = TIME_NOWDATE(3) * 1000000 + TIME_NOWDATE(4) * 10000 + TIME_NOWDATE(5) * 100 + TIME_NOWDATE(6) !ddhhmmss
!    
!    EMITF(:,:,:) = 0.0_RP       ![kg/m2/s]
!    mevent = 1
!    if ( PRC_myrank == myrank_FDNPP ) then
!       do iq = 1, QA_AE_AT
!          ! do j = JS, JE-1
!          ! do i = IS, IE-1
!             
!          ! if (        37.42_RP * CONST_D2R >= (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP .and.  37.42_RP * CONST_D2R < (LAT_REAL(i,j)+LAT_REAL(i,j+1))/2.0_RP &
!          !      .and. 141.03_RP * CONST_D2R >= (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP .and. 141.03_RP * CONST_D2R < (LON_REAL(i,j)+LON_REAL(i+1,j))/2.0_RP ) then
!
!          ! write(*,*) iq, i, j, (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP/CONST_D2R, (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP/CONST_D2R
!          ! write(*,*) PRC_myrank, myrank_FDNPP
!          ! stop
!
!          do ievent = 1, 66
!             if ( time >= emit_time(ievent,1) .and. time < emit_time(ievent,2) ) then
!                EMITF(I_FDNPP,J_FDNPP,iq) = emit_data(ievent) * 1.0E-14_RP
!                mevent = ievent
!             endif
!          enddo
!          
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 3.0864E-11_RP ! [*1e+14*Bq/m2/s] (3000m*3000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.736E-11_RP ! [*1e+14*Bq/m2/s] (4000m*4000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.111E-09_RP ! [*1e+14*Bq/m2/s] (500m*500m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 2.778E-10_RP ! [*1e+14*Bq/m2/s] (1000m*1000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 6.944E-11_RP ! [*1e+14*Bq/m2/s] (2000m*2000m)
!          EMITF(I_FDNPP,J_FDNPP,iq) = EMITF(I_FDNPP,J_FDNPP,iq) / ( HAREA(I_FDNPP,J_FDNPP) * 3600.0_RP )  ! [*1e+14*Bq/m2/s]
!          
!          !   endif
!          ! enddo
!          ! enddo
!       enddo
!    endif
!    
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
!          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DELP(emit_kidx(mevent),i,j) * GRAV * real(dt_AE,kind=RP)
!          do k = KS, KE
!             QEMIT(k,i,j,iq) = emit_coeff(mevent,k) * EMITF(i,j,iq) / DELP(k,i,j) * GRAV * real(dt_AE,kind=RP)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             ! if ( k <= KUP(i,j) ) then
!             ! if ( k == emit_kidx(mevent) ) then
!             !    QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(i,j,iq)
!             ! endif
!             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud
!    !-------------------------------
!    !inside cloud
!    if ( SW_CCN ) then
!       do iq = 1, QA_AE_AT
!          do j = JS, JE
!          do i = IS, IE
!             do k = KS, KE
!                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTIAT(k,i,j) * TCLDF(k,i,j)
!             enddo
!          enddo
!          enddo
!       enddo
!    else
!       do iq = 1, QA_AE_AT
!          do j = JS, JE
!          do i = IS, IE
!             do k = KS, KE
!                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCAT * TCLDF(k,i,j)
!             enddo
!          enddo
!          enddo
!       enddo
!    endif
!
!    !outside cloud
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCTM1(k,i,j,iq) = QTRC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    
!    !subcloud scavenging (wash out)
!    !-------------------------------
!    ! do j = JS, JE
!    ! do i = IS, IE
!    !    do k = KS, KE
!    !       RAINN_CP     (k,i,j) = FRAIN_CP     (k,i,j) / VTR                / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3                )
!    !       RAINN_MP_RAIN(k,i,j) = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
!    !       RAINN_MP_SNOW(k,i,j) = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
!    !    enddo
!    ! enddo
!    ! enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             TRACN                    = QTRCTM1(k,i,j,iq) / MAT(iq) * DENS(k,i,j)
!             TRACRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADAT(k,i,j,iq) )**2 * abs( VTR                - VTER(k,i,j,iq) ) * TRACN * RAINN_CP     (k,i,j)
!             TRACRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_RAIN(k,i,j)
!             TRACRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_SNOW(k,i,j)
!             QWTDEP_CP     (k,i,j,iq) = TRACRN_CP     (k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_RAIN(k,i,j,iq) = TRACRN_MP_RAIN(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_SNOW(k,i,j,iq) = TRACRN_MP_SNOW(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
!             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
!                QTRCTMP(k,i,j,iq) = 0.0_RP
!                TRACRN (k,i,j,iq) = TRACRN_CP(k,i,j,iq) + TRACRN_MP_RAIN(k,i,j,iq) + TRACRN_MP_SNOW(k,i,j,iq)
!                TRACRN_CP     (k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_RAIN(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_SNOW(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!             endif
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud  scavenging (rain out)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN        (k,i,j,iq) = CLTORN_CP(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             if ( QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq) < 0.0_RP ) then
!                CLTORN_CP     (k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_CP     (k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_RAIN(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_RAIN(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_SNOW(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_SNOW(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                QTRCTM2(k,i,j,iq) = 0.0_RP
!             else
!                QTRCTM2(k,i,j,iq) = QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq)
!             endif
!             QWTDEP_CP     (k,i,j,iq) = QWTDEP_CP     (k,i,j,iq) + CLTORN_CP     (k,i,j,iq)
!             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq)
!             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             QTRCTM1       (k,i,j,iq) = QTRCTMP       (k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !dry deposition
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       if ( iq == 1 ) then
!          do j = JS, JE
!          do i = IS, IE
!             CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
!          enddo
!          enddo
!       endif
!       
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
!            ( QMTX,                    & ! [OUT]
!              DENS, DELP, CDVES, dt_AE ) ! [IN]
!
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
!            ( QTRCTM1(:,:,:,iq),             & ! [MODIFIED]
!              DRYF(:,:,iq),                  & ! [OUT]
!              QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
!    enddo
!
!
!    !gravitational settling
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
!            ( QTRCTM1(:,:,:,iq),              & ! [MODIFIED]
!              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
!              GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
!    enddo
!
!    
!    !in rain(MP), in snow(MP)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + QTRCINR(k,i,j,iq)
!             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + QTRCINS(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !re-emission from rain : CP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
!            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
!              WETF_CP(:,:,iq),                                              & ! [OUT]
!              TRACRN_CP(:,:,:,iq), MAT(iq), CLTORN_CP(:,:,:,iq),            & ! [IN]
!              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
!    enddo
!
!    
!    !re-emission from rain & snow : MP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!       !        GDPM, DELP, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!       !        GDPM, DELP, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!       
!       QTRCINR(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-rain(MP)
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!              GDZM, DELZ, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINR(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINR(k,i,j,iq), 0.0_RP )
!          enddo
!       enddo
!       enddo
!          
!       QTRCINS(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-snow(MP)
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!              GDZM, DELZ, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINS(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINS(k,i,j,iq), 0.0_RP )
!          enddo
!       enddo
!       enddo
!
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINP(k,i,j,iq) = QTRCINR(k,i,j,iq) + QTRCINS(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!
!       do j = JS, JE
!       do i = IS, IE
!          WETF(i,j,iq) = WETF_CP(i,j,iq) + WETF_MP_RAIN(i,j,iq) + WETF_MP_SNOW(i,j,iq)
!       enddo
!       enddo
!    enddo
!
!
!    !inside cloud + outside cloud 
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRC (k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
!             QLOAD(k,i,j,iq) = QTRC   (k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!             
!    
!    ! !dry deposition
!    ! !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    if ( iq == 1 ) then
!    !       do j = JS, JE
!    !       do i = IS, IE
!    !          CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
!    !       enddo
!    !       enddo
!    !    endif
!       
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
!    !         ( QMTX,                    & ! [OUT]
!    !           DENS, DELP, CDVES, dt_AE ) ! [IN]
!
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
!    !         ( QTRC(:,:,:,iq),                & ! [MODIFIED]
!    !           DRYF  (:,:,iq),                & ! [OUT]
!    !           QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
!    ! enddo
!
!
!    ! !gravitational settling
!    ! !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
!    !         ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
!    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
!    !           GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
!    ! enddo
!
!    
!    !Decay
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) - QTRC(k,i,j,iq) * decay_ratio(iq) * real(dt_AE,kind=RP)
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !================================================
!    !end aerosol transport
!
!    
!    !aerosol
!    !================================================
!    !================================================
!    !end aerosol
!
!    
!    !flux 
!    !================================================
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          TRACDP(i,j,iq) = WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          ACCUM_TOT_DEP        (i,j,iq) = ACCUM_TOT_DEP        (i,j,iq) + TRACDP      (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP        (i,j,iq) = ACCUM_WET_DEP        (i,j,iq) + WETF        (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_CP     (i,j,iq) = ACCUM_WET_DEP_CP     (i,j,iq) + WETF_CP     (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_MP_RAIN(i,j,iq) = ACCUM_WET_DEP_MP_RAIN(i,j,iq) + WETF_MP_RAIN(i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_MP_SNOW(i,j,iq) = ACCUM_WET_DEP_MP_SNOW(i,j,iq) + WETF_MP_SNOW(i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_DRY_DEP        (i,j,iq) = ACCUM_DRY_DEP        (i,j,iq) + DRYF        (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_GRV_DEP        (i,j,iq) = ACCUM_GRV_DEP        (i,j,iq) + GRAVF       (i,j,iq) * real(dt_AE,kind=RP)
!       enddo
!       enddo
!    enddo
!    !================================================
!    !end flux
!
!    
!    !optical parameters
!    !================================================
!    !================================================
!    !end optical parameters
!
!    
!    !OUTPUT
!    !================================================
!    !flux
!    !-------------------------------    
!    do iq = 1, QA_AE_AT
!       write(cnumber,'(i2.2)') iq
!       !                                  i j q
!       call FILE_HISTORY_in( EMITF       (:,:,iq), 'EFTAT'//cnumber//'_SPRINTARS',          'emission flux (tracer'//cnumber//')',                'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( TRACDP      (:,:,iq), 'TDFTAT'//cnumber//'_SPRINTARS',         'total deposition flux (tracer'//cnumber//')',        'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF        (:,:,iq), 'WEFTAT'//cnumber//'_SPRINTARS',         'wet deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_CP     (:,:,iq), 'WEFTAT'//cnumber//'_CP_SPRINTARS',      'wet deposition flux (tracer'//cnumber//';CP)',       'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_MP_RAIN(:,:,iq), 'WEFTAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(RAIN))', 'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_MP_SNOW(:,:,iq), 'WEFTAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(SNOW))', 'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( DRYF        (:,:,iq), 'DRFTAT'//cnumber//'_SPRINTARS',         'dry deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( GRAVF       (:,:,iq), 'GRFTAT'//cnumber//'_SPRINTARS',         'gravitational settling flux (tracer'//cnumber//')',  'kg/m2/s', dim_type = 'XY' )
!
!       !                                           i j q
!       call FILE_HISTORY_in( ACCUM_TOT_DEP        (:,:,iq), 'ACTDAT'//cnumber//'_SPRINTARS',         'total deposition (tracer'//cnumber//')',                  'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP        (:,:,iq), 'ACWDAT'//cnumber//'_SPRINTARS',         'wet deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_CP     (:,:,iq), 'ACWDAT'//cnumber//'_CP_SPRINTARS',      'wet deposition (tracer'//cnumber//';CP)',                 'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_RAIN(:,:,iq), 'ACWDAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(RAIN))',           'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_SNOW(:,:,iq), 'ACWDAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(SNOW))',           'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_DRY_DEP        (:,:,iq), 'ACDDAT'//cnumber//'_SPRINTARS',         'dry deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_GRV_DEP        (:,:,iq), 'ACGDAT'//cnumber//'_SPRINTARS',         'gravitational settling deposition (tracer'//cnumber//')', 'kg/m2', dim_type = 'XY' )
!    enddo
!       
!    !aerosol
!    !-------------------------------
!
!    !optical properties
!    !-------------------------------
!    
!    !clear sky diagnoses
!    !-------------------------------
!    
!    !in-precipitation
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       write(cnumber,'(i2.2)') iq
!       !                             k i j q
!       call FILE_HISTORY_in( QTRCINP(:,:,:,iq), 'tracer'//cnumber//'_INPREC_SPRINTARS', 'Ratio of tracer'//cnumber//' in PREC', 'kg/kg', dim_type = 'ZXY' )
!       call FILE_HISTORY_in( QTRCINR(:,:,:,iq), 'tracer'//cnumber//'_INRAIN_SPRINTARS', 'Ratio of tracer'//cnumber//' in RAIN', 'kg/kg', dim_type = 'ZXY' )
!       call FILE_HISTORY_in( QTRCINS(:,:,:,iq), 'tracer'//cnumber//'_INSNOW_SPRINTARS', 'Ratio of tracer'//cnumber//' in SNOW', 'kg/kg', dim_type = 'ZXY' )
!    enddo
!    !================================================
!    !end OUTPUT
!
!    return
!  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY
  !-----------------------------------------------------------------------------
  !> SPRINTARS Tracer (sulfate proxy)
  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY2 &
       ( QTRC,           & ! [MODIFIED]
         CDVE,           & ! [IN]
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
         KUP,            & ! [IN]
         KUPMAX,         & ! [IN]
         DPKUP,          & ! [IN]
         GDPM,           & ! [IN]
         DELP,           & ! [IN]
         GDZM,           & ! [IN]
         DELZ,           & ! [IN]
         RNWR,           & ! [IN]
         RNWR_CP,        & ! [IN]
         RNWR_MP_RAIN,   & ! [IN]
         RNWR_MP_SNOW,   & ! [IN]
         LON_REAL,       & ! [IN]
         LAT_REAL,       & ! [IN]
         ACTIAT,         & ! [IN]
         SW_CCN,         & ! [IN]
         TIME_NOWDATE,   & ! [IN]
         dt_AE,          & ! [IN]
         QA_AE_AT        ) ! [IN]
    implicit none

    ! input
    !-------------------------
    real(DP), intent(in)    :: dt_AE
    integer,  intent(in)    :: QA_AE_AT
    
    real(RP), intent(in)    :: CDVE               (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
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
    integer,  intent(in)    :: KUP                (IA,JA) !< uplift level
    integer,  intent(in)    :: KUPMAX                     !< maximum uplift level
    real(RP), intent(in)    :: DPKUP              (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA) !< pressure at half levels
    real(RP), intent(in)    :: DELP          (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA) !< altitude at half levels
    real(RP), intent(in)    :: DELZ          (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
    real(RP), intent(in)    :: RNWR          (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA) !< ratio of rain to cloud
    real(RP), intent(in)    :: LON_REAL           (IA,JA) !< longitude
    real(RP), intent(in)    :: LAT_REAL           (IA,JA) !< latitude
    real(RP), intent(in)    :: ACTIAT        (KA,  IA,JA) !< ratio of activated aerosols [0-1]
    logical,  intent(in)    :: SW_CCN                     !< switch for CCN scheme (true->use ACTI)
    integer,  intent(in)    :: TIME_NOWDATE(6)            !< current time
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_AT)
    
    ! output
    !-------------------------

    ! internal work
    !-------------------------
    !set parameter
    real(RP) :: decay_ratio(QA_AE_AT)

    real(RP) :: MAT           (QA_AE_AT)
    real(RP) :: RADAT(KA,IA,JA,QA_AE_AT) !< tracer particle radius      [m]
    real(RP) :: BRH1, BRH2, BRH3         !<
    real(RP) :: VTER (KA,IA,JA,QA_AE_AT) !< terminal velocity of tracer [m/s]

    !emission
    real(RP) :: QEMIT(KA,IA,JA,QA_AE_AT) !emission mixing ratio [kg/kg]
    real(RP) :: EMITF   (IA,JA,QA_AE_AT) !tracer emission flux  [kg/m2/s]
    
    !incloud
    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_AT)  !< outside cloud
    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_AT)  !< inside  cloud
    
    !subcloud scavenging (wash out)
    ! real(RP) :: RAINN        (KA,IA,JA)                     !< rain number concentrarion        [1/m3]
    ! real(RP) :: RAINN_CP     (KA,IA,JA)               
    ! real(RP) :: RAINN_MP_RAIN(KA,IA,JA)              
    ! real(RP) :: RAINN_MP_SNOW(KA,IA,JA)               
    real(RP) :: TRACN                             !< number density of tracer [1/m3]
    real(RP) :: TRACRN        (KA,IA,JA,QA_AE_AT) !< number density of tracer entering into raindrops [1/s/m3] <=> collision frequency between tracer and rain
    real(RP) :: TRACRN_CP     (KA,IA,JA,QA_AE_AT)
    real(RP) :: TRACRN_MP_RAIN(KA,IA,JA,QA_AE_AT)
    real(RP) :: TRACRN_MP_SNOW(KA,IA,JA,QA_AE_AT)
    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_AT) !< mixing ratio by wet deposition   [kg/kg] 
    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_AT)
    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_AT)
    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_AT)
    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_AT) !< tracers that do not cllide with rain particles [kg/kg]
    
    !incloud scavenging (rain out)
    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_AT) !< mixing ratio moving cloud to rain particles [kg/kg]
    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_AT)
    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_AT)
    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_AT)
    
    !re-emission from rain
    real(RP) :: WETF           (IA,JA,QA_AE_AT) !< wet deposition flux
    real(RP) :: WETF_CP        (IA,JA,QA_AE_AT)
    real(RP) :: WETF_MP_RAIN   (IA,JA,QA_AE_AT)
    real(RP) :: WETF_MP_SNOW   (IA,JA,QA_AE_AT)
    real(RP) :: QTRCINR     (KA,IA,JA,QA_AE_AT) !< mixing ratio in rain(MP)
    real(RP) :: QTRCINS     (KA,IA,JA,QA_AE_AT) !< mixing ratio in snow(MP)

    !dry deposition
    real(RP) :: CDVES   (IA,JA)           !< 
    real(RP) :: QMTX (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
    real(RP) :: DRYF    (IA,JA,QA_AE_AT)  !< dry deposition flux
    
    !gravitational settling
    real(RP) :: QLOAD(KA,IA,JA,QA_AE_AT)  !< mixing ratio for radiation [kg/kg]
    real(RP) :: GRAVF   (IA,JA,QA_AE_AT)  !< gravitational settling flux 
    
    !flux
    real(RP) :: TRACDP  (IA,JA,QA_AE_AT)  !< sulfate total deposition
    
    ! parameter
    !-------------------------
    real(RP), parameter :: hltime(NAT) = (/ 30.2_RP * 365.0_RP * 86400.0_RP /) ! half life time [s] : 137Cs (30.1yrs)
    real(RP), parameter :: DDV(NAT) = (/ 2.0E-3_RP /)                !< deposition velocity
    real(RP), parameter :: &
         RADTRAC(NRH)      &    !< sulfate particle radius
         !    RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
         = (/ 0.695E-7_RP, 0.850E-7_RP, 0.950E-7_RP, 0.103E-6_RP, 0.122E-6_RP, 0.157E-6_RP, 0.195E-6_RP, 0.231E-6_RP /)
    real(RP), parameter :: DENSAT(NAT) = (/ 1.769E+3_RP /) !< tracer particle density [kg/m3]
    real(RP), parameter :: AA    (NAT) = (/ 1.000E-4_RP /) !< collision efficiency
    real(RP), parameter :: SIGMA (NAT) = (/ 1.526E+0_RP /) !< GSD of log-normal distribution

    ! real(RP), parameter :: THRESPAT = 1.39E-4_RP ![kg/m2/s] = 0.50 [mm/h]
    ! real(RP), parameter :: THRESPAT = 5.56E-5_RP ![kg/m2/s] = 0.20 [mm/h]
    real(RP), parameter :: THRESPAT = 2.78E-5_RP ![kg/m2/s] = 0.10 [mm/h]
    ! real(RP), parameter :: THRESPAT = 1.39E-5_RP ![kg/m2/s] = 0.05 [mm/h]
    ! real(RP), parameter :: THRESPAT = 5.56E-6_RP ![kg/m2/s] = 0.02 [mm/h]
    ! real(RP), parameter :: THRESPAT = 2.78E-6_RP ![kg/m2/s] = 0.01 [mm/h]
    
    ! others 
    !-------------------------
    integer   :: k, i, j, iq
    integer   :: IRH
    character :: cnumber*2
    integer   :: time           ! ddhhmmss
    integer   :: ievent, mevent
    logical, save :: OFIRST = .true. 

    !----------------------------------------------------
    
    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - tracer SULFATE_PROXY2'

    !once 
    !================================================
    if ( OFIRST ) then
       OFIRST = .false.
       if ( PRC_myrank == myrank_FDNPP ) then
          write(*,*) 'FDNPP PRC_myrank =', PRC_myrank
          write(*,*) 'FDNPP =', I_FDNPP, J_FDNPP, LON_REAL(I_FDNPP,J_FDNPP)/CONST_D2R, LAT_REAL(I_FDNPP,J_FDNPP)/CONST_D2R
       endif
    endif

    !================================================
    !end once 
    
    !setup parameter 
    !================================================
    !MAT
    !-------------------------------
    do iq = 1, QA_AE_AT
       MAT(iq) = 4.0_RP / 3.0_RP * PI * RADTRAC(1)**3 * DENSAT(iq) * exp( 4.5_RP * log(SIGMA(iq))**2 )
    enddo
    
    !RADAT
    !-------------------------------
    do iq = 1, QA_AE_AT
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
                RADAT(k,i,j,iq) = RADTRAC(IRH)
             else
                BRH1 =  RH(k,i,j) - BRH(IRH)
                BRH2 = BRH(IRH+1) -  RH(k,i,j)
                BRH3 = BRH(IRH+1) - BRH(IRH)
                RADAT(k,i,j,iq) = ( BRH1 * RADTRAC(IRH+1) + BRH2 * RADTRAC(IRH) ) / BRH3
             endif

          enddo
       enddo
       enddo
    enddo
    
    !VTER
    !-------------------------------
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             VTER(k,i,j,iq) = 2.0_RP * (              ( RADTRAC(1)/RADAT(k,i,j,iq) )**3   * DENSAT(iq) &
                                         + ( 1.0_RP - ( RADTRAC(1)/RADAT(k,i,j,iq) )**3 ) * DENSW    ) &
                                     * RADAT(k,i,j,iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(iq))**2 )
          enddo
       enddo
       enddo
    enddo
    
    !Decay based on half life
    !-------------------------------
    do iq = 1, QA_AE_AT
       decay_ratio(iq) = log(2.0_RP) / hltime(iq)
    enddo
    
    !================================================
    !end setup parameter


    !aerosol transport 
    !================================================
      
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_AT
       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
              DENS1(:,:,:), DENS(:,:,:), THRESPAT                                                                   ) ! [IN]
    enddo

    
    !emission
    !-------------------------------
    time = TIME_NOWDATE(3) * 1000000 + TIME_NOWDATE(4) * 10000 + TIME_NOWDATE(5) * 100 + TIME_NOWDATE(6) !ddhhmmss
    
    EMITF(:,:,:) = 0.0_RP       ![kg/m2/s]
    mevent = 1
    if ( PRC_myrank == myrank_FDNPP ) then
       do iq = 1, QA_AE_AT
          ! do j = JS, JE-1
          ! do i = IS, IE-1
             
          ! if (        37.42_RP * CONST_D2R >= (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP .and.  37.42_RP * CONST_D2R < (LAT_REAL(i,j)+LAT_REAL(i,j+1))/2.0_RP &
          !      .and. 141.03_RP * CONST_D2R >= (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP .and. 141.03_RP * CONST_D2R < (LON_REAL(i,j)+LON_REAL(i+1,j))/2.0_RP ) then

          ! write(*,*) iq, i, j, (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP/CONST_D2R, (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP/CONST_D2R
          ! write(*,*) PRC_myrank, myrank_FDNPP
          ! stop

          do ievent = 1, 66
             if ( time >= emit_time(ievent,1) .and. time < emit_time(ievent,2) ) then
                EMITF(I_FDNPP,J_FDNPP,iq) = emit_data(ievent) * 1.0E-14_RP
                mevent = ievent
             endif
          enddo
          
          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 3.0864E-11_RP ! [*1e+14*Bq/m2/s] (3000m*3000m)
          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.736E-11_RP ! [*1e+14*Bq/m2/s] (4000m*4000m)
          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.111E-09_RP ! [*1e+14*Bq/m2/s] (500m*500m)
          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 2.778E-10_RP ! [*1e+14*Bq/m2/s] (1000m*1000m)
          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 6.944E-11_RP ! [*1e+14*Bq/m2/s] (2000m*2000m)
          EMITF(I_FDNPP,J_FDNPP,iq) = EMITF(I_FDNPP,J_FDNPP,iq) / ( HAREA(I_FDNPP,J_FDNPP) * 3600.0_RP )  ! [*1e+14*Bq/m2/s]
          
          !   endif
          ! enddo
          ! enddo
       enddo
    endif
    
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DELP(emit_kidx(mevent),i,j) * GRAV * real(dt_AE,kind=RP)
          do k = KS, KE
             QEMIT(k,i,j,iq) = emit_coeff(mevent,k) * EMITF(i,j,iq) / DELP(k,i,j) * GRAV * real(dt_AE,kind=RP)
          enddo
       enddo
       enddo
    enddo

    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             ! if ( k <= KUP(i,j) ) then
             ! if ( k == emit_kidx(mevent) ) then
             !    QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(i,j,iq)
             ! endif
             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    
    !incloud
    !-------------------------------
    !inside cloud
    if ( SW_CCN ) then
       do iq = 1, QA_AE_AT
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTIAT(k,i,j) * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    else
       do iq = 1, QA_AE_AT
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCAT * TCLDF(k,i,j)
             enddo
          enddo
          enddo
       enddo
    endif

    !outside cloud
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRCTM1(k,i,j,iq) = QTRC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo
    
    
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

    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             TRACN                    = QTRCTM1(k,i,j,iq) / MAT(iq) * DENS(k,i,j)
             if ( FRAIN_CP     (k-1,i,j) >= THRESPAT ) then
                TRACRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADAT(k,i,j,iq) )**2 * abs( VTR                - VTER(k,i,j,iq) ) * TRACN * RAINN_CP     (k,i,j)
             else
                TRACRN_CP     (k,i,j,iq) = 0.0_RP
             endif
             if ( FRAIN_MP_RAIN(k-1,i,j) >= THRESPAT ) then
                TRACRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_RAIN(k,i,j)
             else
                TRACRN_MP_RAIN(k,i,j,iq) = 0.0_RP
             endif
             if ( FRAIN_MP_SNOW(k-1,i,j) >= THRESPAT ) then
                TRACRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_SNOW(k,i,j)
             else
                TRACRN_MP_SNOW(k,i,j,iq) = 0.0_RP
             endif
             QWTDEP_CP     (k,i,j,iq) = TRACRN_CP     (k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_RAIN(k,i,j,iq) = TRACRN_MP_RAIN(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QWTDEP_MP_SNOW(k,i,j,iq) = TRACRN_MP_SNOW(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
                QTRCTMP(k,i,j,iq) = 0.0_RP
                TRACRN (k,i,j,iq) = TRACRN_CP(k,i,j,iq) + TRACRN_MP_RAIN(k,i,j,iq) + TRACRN_MP_SNOW(k,i,j,iq)
                TRACRN_CP     (k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
                TRACRN_MP_RAIN(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
                TRACRN_MP_SNOW(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
             endif
          enddo
       enddo
       enddo
    enddo

    
    !incloud  scavenging (rain out)
    !-------------------------------
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( FRAIN_CP     (k-1,i,j) >= THRESPAT ) then
                CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
             else
                CLTORN_CP     (k,i,j,iq) = 0.0_RP
             endif
             if ( FRAIN_MP_RAIN(k-1,i,j) >= THRESPAT ) then
                CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
             else
                CLTORN_MP_RAIN(k,i,j,iq) = 0.0_RP
             endif
             if ( FRAIN_MP_SNOW(k-1,i,j) >= THRESPAT ) then
                CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
             else
                CLTORN_MP_SNOW(k,i,j,iq) = 0.0_RP
             endif
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
    do iq = 1, QA_AE_AT
       if ( iq == 1 ) then
          do j = JS, JE
          do i = IS, IE
             CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
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
    do iq = 1, QA_AE_AT
       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
            ( QTRCTM1(:,:,:,iq),              & ! [MODIFIED]
              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
              GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
    enddo

    
    !in rain(MP), in snow(MP)
    !-------------------------------
    do iq = 1, QA_AE_AT
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
    do iq = 1, QA_AE_AT
       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
              WETF_CP(:,:,iq),                                              & ! [OUT]
              TRACRN_CP(:,:,:,iq), MAT(iq), CLTORN_CP(:,:,:,iq),            & ! [IN]
              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
    enddo

    
    !re-emission from rain & snow : MP
    !-------------------------------
    do iq = 1, QA_AE_AT
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
    do iq = 1, QA_AE_AT
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
    ! do iq = 1, QA_AE_AT
    !    if ( iq == 1 ) then
    !       do j = JS, JE
    !       do i = IS, IE
    !          CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
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
    ! do iq = 1, QA_AE_AT
    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
    !         ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
    !           GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
    ! enddo

    
    !Decay
    !-------------------------------
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) - QTRC(k,i,j,iq) * decay_ratio(iq) * real(dt_AE,kind=RP)
          enddo
       enddo
       enddo
    enddo
    
    !================================================
    !end aerosol transport

    
    !aerosol
    !================================================
    !================================================
    !end aerosol

    
    !flux 
    !================================================
    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          TRACDP(i,j,iq) = WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
       enddo
       enddo
    enddo

    do iq = 1, QA_AE_AT
       do j = JS, JE
       do i = IS, IE
          ACCUM_TOT_DEP        (i,j,iq) = ACCUM_TOT_DEP        (i,j,iq) + TRACDP      (i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_WET_DEP        (i,j,iq) = ACCUM_WET_DEP        (i,j,iq) + WETF        (i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_WET_DEP_CP     (i,j,iq) = ACCUM_WET_DEP_CP     (i,j,iq) + WETF_CP     (i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_WET_DEP_MP_RAIN(i,j,iq) = ACCUM_WET_DEP_MP_RAIN(i,j,iq) + WETF_MP_RAIN(i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_WET_DEP_MP_SNOW(i,j,iq) = ACCUM_WET_DEP_MP_SNOW(i,j,iq) + WETF_MP_SNOW(i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_DRY_DEP        (i,j,iq) = ACCUM_DRY_DEP        (i,j,iq) + DRYF        (i,j,iq) * real(dt_AE,kind=RP)
          ACCUM_GRV_DEP        (i,j,iq) = ACCUM_GRV_DEP        (i,j,iq) + GRAVF       (i,j,iq) * real(dt_AE,kind=RP)
       enddo
       enddo
    enddo
    !================================================
    !end flux

    
    !optical parameters
    !================================================
    !================================================
    !end optical parameters

    
    !OUTPUT
    !================================================
    !flux
    !-------------------------------    
    do iq = 1, QA_AE_AT
       write(cnumber,'(i2.2)') iq
       !                                  i j q
       call FILE_HISTORY_in( EMITF       (:,:,iq), 'EFTAT'//cnumber//'_SPRINTARS',          'emission flux (tracer'//cnumber//')',                'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( TRACDP      (:,:,iq), 'TDFTAT'//cnumber//'_SPRINTARS',         'total deposition flux (tracer'//cnumber//')',        'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( WETF        (:,:,iq), 'WEFTAT'//cnumber//'_SPRINTARS',         'wet deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( WETF_CP     (:,:,iq), 'WEFTAT'//cnumber//'_CP_SPRINTARS',      'wet deposition flux (tracer'//cnumber//';CP)',       'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( WETF_MP_RAIN(:,:,iq), 'WEFTAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(RAIN))', 'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( WETF_MP_SNOW(:,:,iq), 'WEFTAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(SNOW))', 'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( DRYF        (:,:,iq), 'DRFTAT'//cnumber//'_SPRINTARS',         'dry deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
       call FILE_HISTORY_in( GRAVF       (:,:,iq), 'GRFTAT'//cnumber//'_SPRINTARS',         'gravitational settling flux (tracer'//cnumber//')',  'kg/m2/s', dim_type = 'XY' )

       !                                           i j q
       call FILE_HISTORY_in( ACCUM_TOT_DEP        (:,:,iq), 'ACTDAT'//cnumber//'_SPRINTARS',         'total deposition (tracer'//cnumber//')',                  'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_WET_DEP        (:,:,iq), 'ACWDAT'//cnumber//'_SPRINTARS',         'wet deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_WET_DEP_CP     (:,:,iq), 'ACWDAT'//cnumber//'_CP_SPRINTARS',      'wet deposition (tracer'//cnumber//';CP)',                 'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_RAIN(:,:,iq), 'ACWDAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(RAIN))',           'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_SNOW(:,:,iq), 'ACWDAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(SNOW))',           'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_DRY_DEP        (:,:,iq), 'ACDDAT'//cnumber//'_SPRINTARS',         'dry deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
       call FILE_HISTORY_in( ACCUM_GRV_DEP        (:,:,iq), 'ACGDAT'//cnumber//'_SPRINTARS',         'gravitational settling deposition (tracer'//cnumber//')', 'kg/m2', dim_type = 'XY' )
    enddo
       
    !aerosol
    !-------------------------------

    !optical properties
    !-------------------------------
    
    !clear sky diagnoses
    !-------------------------------
    
    !in-precipitation
    !-------------------------------
    do iq = 1, QA_AE_AT
       write(cnumber,'(i2.2)') iq
       !                             k i j q
       call FILE_HISTORY_in( QTRCINP(:,:,:,iq), 'tracer'//cnumber//'_INPREC_SPRINTARS', 'Ratio of tracer'//cnumber//' in PREC', 'kg/kg', dim_type = 'ZXY' )
       call FILE_HISTORY_in( QTRCINR(:,:,:,iq), 'tracer'//cnumber//'_INRAIN_SPRINTARS', 'Ratio of tracer'//cnumber//' in RAIN', 'kg/kg', dim_type = 'ZXY' )
       call FILE_HISTORY_in( QTRCINS(:,:,:,iq), 'tracer'//cnumber//'_INSNOW_SPRINTARS', 'Ratio of tracer'//cnumber//' in SNOW', 'kg/kg', dim_type = 'ZXY' )
    enddo
    !================================================
    !end OUTPUT

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SULFATE_PROXY2
  !-----------------------------------------------------------------------------
!  !-----------------------------------------------------------------------------
!  !> SPRINTARS Tracer (silicate glass proxy)
!  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SILICATE_PROXY &
!       ( QTRC,           & ! [MODIFIED]
!         CDVE,           & ! [IN]
!         RH,             & ! [IN]
!         DENS,           & ! [IN]
!         DENS1,          & ! [IN]
!         FRAIN_MP_RAIN1, & ! [IN]
!         FRAIN_MP_SNOW1, & ! [IN]
!         FRAIN,          & ! [IN]
!         FRAIN_CP,       & ! [IN]
!         FRAIN_MP_RAIN,  & ! [IN]
!         FRAIN_MP_SNOW,  & ! [IN]
!         VTR_MP_RAIN,    & ! [IN]
!         VTR_MP_SNOW,    & ! [IN]
!         RADR_MP_RAIN,   & ! [IN]
!         RADR_MP_SNOW,   & ! [IN]
!         RAINN_CP,       & ! [IN]
!         RAINN_MP_RAIN,  & ! [IN]
!         RAINN_MP_SNOW,  & ! [IN]
!         TCLDF,          & ! [IN]
!         KUP,            & ! [IN]
!         KUPMAX,         & ! [IN]
!         DPKUP,          & ! [IN]
!         GDPM,           & ! [IN]
!         DELP,           & ! [IN]
!         GDZM,           & ! [IN]
!         DELZ,           & ! [IN]
!         RNWR,           & ! [IN]
!         RNWR_CP,        & ! [IN]
!         RNWR_MP_RAIN,   & ! [IN]
!         RNWR_MP_SNOW,   & ! [IN]
!         LON_REAL,       & ! [IN]
!         LAT_REAL,       & ! [IN]
!         ACTIAT,         & ! [IN]
!         SW_CCN,         & ! [IN]
!         TIME_NOWDATE,   & ! [IN]
!         dt_AE,          & ! [IN]
!         QA_AE_AT        ) ! [IN]
!    implicit none
!
!    ! input
!    !-------------------------
!    real(DP), intent(in)    :: dt_AE
!    integer,  intent(in)    :: QA_AE_AT
!    
!    real(RP), intent(in)    :: CDVE               (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
!    real(RP), intent(in)    :: RH            (KA,  IA,JA) !< relative humidity
!    real(RP), intent(in)    :: DENS          (KA,  IA,JA) !< air density
!    real(RP), intent(in)    :: DENS1         (KA,  IA,JA) !< air density   on 1 step before
!    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA) !< precipitation flux (liq    ) : MP on 1 step before
!    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA) !< precipitation flux (    ice) : MP on 1 step before
!    real(RP), intent(in)    :: FRAIN         (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP + MP
!    real(RP), intent(in)    :: FRAIN_CP      (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP
!    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA) !< precipitation flux (liq    ) : MP
!    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA) !< precipitation flux (    ice) : MP
!    real(RP), intent(in)    :: VTR_MP_RAIN   (KA,  IA,JA) !< terminal velocity of rain (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: VTR_MP_SNOW   (KA,  IA,JA) !< terminal velocity of snow (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: RADR_MP_RAIN  (KA,  IA,JA) !< effective radius of rain (microphysic scheme) [m]
!    real(RP), intent(in)    :: RADR_MP_SNOW  (KA,  IA,JA) !< effective radius of snow (microphysic scheme) [m]
!    real(RP), intent(in)    :: RAINN_CP      (KA,  IA,JA) !< rain number concentration (CP)      [1/m3]
!    real(RP), intent(in)    :: RAINN_MP_RAIN (KA,  IA,JA) !< rain number concentration (MP;RAIN) [1/m3]
!    real(RP), intent(in)    :: RAINN_MP_SNOW (KA,  IA,JA) !< rain number concentration (MP;SNOW) [1/m3]
!    real(RP), intent(in)    :: TCLDF         (KA,  IA,JA) !< cloud fraction (CP+MP)
!    integer,  intent(in)    :: KUP                (IA,JA) !< uplift level
!    integer,  intent(in)    :: KUPMAX                     !< maximum uplift level
!    real(RP), intent(in)    :: DPKUP              (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
!    real(RP), intent(in)    :: GDPM          (0:KA,IA,JA) !< pressure at half levels
!    real(RP), intent(in)    :: DELP          (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
!    real(RP), intent(in)    :: GDZM          (0:KA,IA,JA) !< altitude at half levels
!    real(RP), intent(in)    :: DELZ          (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
!    real(RP), intent(in)    :: RNWR          (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_CP       (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_RAIN  (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_SNOW  (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: LON_REAL           (IA,JA) !< longitude
!    real(RP), intent(in)    :: LAT_REAL           (IA,JA) !< latitude
!    real(RP), intent(in)    :: ACTIAT        (KA,  IA,JA) !< ratio of activated aerosols [0-1]
!    logical,  intent(in)    :: SW_CCN                     !< switch for CCN scheme (true->use ACTI)
!    integer,  intent(in)    :: TIME_NOWDATE(6)            !< current time
!    
!    ! modified
!    !-------------------------
!    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_AT)
!    
!    ! output
!    !-------------------------
!
!    ! internal work
!    !-------------------------
!    !set parameter
!    real(RP) :: decay_ratio(QA_AE_AT)
!
!    real(RP) :: MAT           (QA_AE_AT)
!    real(RP) :: RADAT(KA,IA,JA,QA_AE_AT) !< tracer particle radius      [m]
!    real(RP) :: BRH1, BRH2, BRH3         !<
!    real(RP) :: VTER (KA,IA,JA,QA_AE_AT) !< terminal velocity of tracer [m/s]
!
!    !emission
!    real(RP) :: QEMIT(KA,IA,JA,QA_AE_AT) !emission mixing ratio [kg/kg]
!    real(RP) :: EMITF   (IA,JA,QA_AE_AT) !tracer emission flux  [kg/m2/s]
!    
!    !incloud
!    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_AT)  !< outside cloud
!    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_AT)  !< inside  cloud
!    
!    !subcloud scavenging (wash out)
!    ! real(RP) :: RAINN        (KA,IA,JA)                     !< rain number concentrarion        [1/m3]
!    ! real(RP) :: RAINN_CP     (KA,IA,JA)               
!    ! real(RP) :: RAINN_MP_RAIN(KA,IA,JA)              
!    ! real(RP) :: RAINN_MP_SNOW(KA,IA,JA)               
!    real(RP) :: TRACN                             !< number density of tracer [1/m3]
!    real(RP) :: TRACRN        (KA,IA,JA,QA_AE_AT) !< number density of tracer entering into raindrops [1/s/m3] <=> collision frequency between tracer and rain
!    real(RP) :: TRACRN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_AT) !< mixing ratio by wet deposition   [kg/kg] 
!    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_AT) !< tracers that do not cllide with rain particles [kg/kg]
!    
!    !incloud scavenging (rain out)
!    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_AT) !< mixing ratio moving cloud to rain particles [kg/kg]
!    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    
!    !re-emission from rain
!    real(RP) :: WETF           (IA,JA,QA_AE_AT) !< wet deposition flux
!    real(RP) :: WETF_CP        (IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_RAIN   (IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_SNOW   (IA,JA,QA_AE_AT)
!    real(RP) :: QTRCINR     (KA,IA,JA,QA_AE_AT) !< mixing ratio in rain(MP)
!    real(RP) :: QTRCINS     (KA,IA,JA,QA_AE_AT) !< mixing ratio in snow(MP)
!
!    !dry deposition
!    real(RP) :: CDVES   (IA,JA)           !< 
!    real(RP) :: QMTX (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
!    real(RP) :: DRYF    (IA,JA,QA_AE_AT)  !< dry deposition flux
!    
!    !gravitational settling
!    real(RP) :: QLOAD(KA,IA,JA,QA_AE_AT)  !< mixing ratio for radiation [kg/kg]
!    real(RP) :: GRAVF   (IA,JA,QA_AE_AT)  !< gravitational settling flux 
!    
!    !flux
!    real(RP) :: TRACDP  (IA,JA,QA_AE_AT)  !< sulfate total deposition
!    
!    ! parameter
!    !-------------------------
!    real(RP), parameter :: hltime(NAT) = (/ 30.2_RP * 365.0_RP * 86400.0_RP /) ! half life time [s] : 137Cs (30.1yrs)
!    ! real(RP), parameter :: DDV(NAT) = (/ 2.0E-3_RP /)                !< deposition velocity
!    ! real(RP), parameter :: &
!    !      RADTRAC(NRH)      &    !< sulfate particle radius
!    !      !    RH = 0.00    RH = 0.50    RH = 0.70    RH = 0.80    RH = 0.90    RH = 0.95    RH = 0.98    RH = 0.99
!    !      = (/ 0.695E-7_RP, 0.850E-7_RP, 0.950E-7_RP, 0.103E-6_RP, 0.122E-6_RP, 0.157E-6_RP, 0.195E-6_RP, 0.231E-6_RP /)
!    ! real(RP), parameter :: DENSAT(NAT) = (/ 1.769E+3_RP /) !< tracer particle density [kg/m3]
!    ! real(RP), parameter :: AA    (NAT) = (/ 1.000E-4_RP /) !< collision efficiency
!    ! real(RP), parameter :: SIGMA (NAT) = (/ 1.526E+0_RP /) !< GSD of log-normal distribution
!    real(RP), parameter :: RADTRAC(NAT) = (/ 3.200E-6_RP /) !< silicate glass particle radius [m]
!    real(RP), parameter :: DENSAT (NAT) = (/ 2.200E+3_RP /) !< tracer particle density [kg/m3]
!    real(RP), parameter :: AA     (NAT) = (/ 0.631E+0_RP /) !< collision efficiency
!    real(RP), parameter :: SIGMA  (NAT) = (/ 1.100E+0_RP /) !< GSD of log-normal distribution
!    
!    ! others 
!    !-------------------------
!    integer   :: k, i, j, iq
!    integer   :: IRH
!    character :: cnumber*2
!    integer   :: time           ! ddhhmmss
!    integer   :: ievent, mevent
!    logical, save :: OFIRST = .true. 
!
!    !----------------------------------------------------
!    
!    LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - tracer'
!
!    !once 
!    !================================================
!    if ( OFIRST ) then
!       OFIRST = .false.
!       if ( PRC_myrank == myrank_FDNPP ) then
!          write(*,*) 'FDNPP PRC_myrank =', PRC_myrank
!          write(*,*) 'FDNPP =', I_FDNPP, J_FDNPP, LON_REAL(I_FDNPP,J_FDNPP)/CONST_D2R, LAT_REAL(I_FDNPP,J_FDNPP)/CONST_D2R
!       endif
!    endif
!
!    !================================================
!    !end once 
!    
!    !setup parameter 
!    !================================================
!    !MAT
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       MAT(iq) = 4.0_RP / 3.0_RP * PI * RADTRAC(1)**3 * DENSAT(iq) * exp( 4.5_RP * log(SIGMA(iq))**2 )
!    enddo
!    
!    !RADAT
!    !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    do j = JS, JE
!    !    do i = IS, IE
!    !       do k = KS, KE
!          
!    !          if (                               RH(k,i,j) < BRH(2) ) then ! RH =  0 - 50 %
!    !             IRH = 1
!    !          elseif ( RH(k,i,j) >= BRH(2) .and. RH(k,i,j) < BRH(3) ) then ! RH = 50 - 70 %
!    !             IRH = 2
!    !          elseif ( RH(k,i,j) >= BRH(3) .and. RH(k,i,j) < BRH(4) ) then ! RH = 70 - 80 %
!    !             IRH = 3
!    !          elseif ( RH(k,i,j) >= BRH(4) .and. RH(k,i,j) < BRH(5) ) then ! RH = 80 - 90 %
!    !             IRH = 4
!    !          elseif ( RH(k,i,j) >= BRH(5) .and. RH(k,i,j) < BRH(6) ) then ! RH = 90 - 95 %
!    !             IRH = 5
!    !          elseif ( RH(k,i,j) >= BRH(6) .and. RH(k,i,j) < BRH(7) ) then ! RH = 95 - 98 %
!    !             IRH = 6
!    !          elseif ( RH(k,i,j) >= BRH(7) .and. RH(k,i,j) < BRH(8) ) then ! RH = 98 - 99 %
!    !             IRH = 7
!    !          else                                                         ! RH = 99 -    %
!    !             IRH = 8           
!    !          endif
!
!    !          if ( IRH == 8 ) then
!    !             RADAT(k,i,j,iq) = RADTRAC(IRH)
!    !          else
!    !             BRH1 =  RH(k,i,j) - BRH(IRH)
!    !             BRH2 = BRH(IRH+1) -  RH(k,i,j)
!    !             BRH3 = BRH(IRH+1) - BRH(IRH)
!    !             RADAT(k,i,j,iq) = ( BRH1 * RADTRAC(IRH+1) + BRH2 * RADTRAC(IRH) ) / BRH3
!    !          endif
!
!    !       enddo
!    !    enddo
!    !    enddo
!    ! enddo
!    RADAT(:,:,:,:) = RADTRAC(1) 
!    
!    !VTER
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             ! VTER(k,i,j,iq) = 2.0_RP * (              ( RADTRAC(1)/RADAT(k,i,j,iq) )**3   * DENSAT(iq) &
!             !                             + ( 1.0_RP - ( RADTRAC(1)/RADAT(k,i,j,iq) )**3 ) * DENSW    ) &
!             !                         * RADAT(k,i,j,iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(iq))**2 )
!             VTER(k,i,j,iq) = 2.0_RP * DENSAT(iq) * RADAT(k,i,j,iq)**2 * GRAV / 9.0_RP / NU * exp( 6.0_RP * log(SIGMA(iq))**2 )
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !Decay based on half life
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       decay_ratio(iq) = log(2.0_RP) / hltime(iq)
!    enddo
!    
!    !================================================
!    !end setup parameter
!
!
!    !aerosol transport 
!    !================================================
!      
!    !in rain(MP), in snow(MP)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
!            ( QTRC(:,:,:,iq),                                                                                       & ! [INOUT]
!              QTRCINR(:,:,:,iq), QTRCINS(:,:,:,iq),                                                                 & ! [OUT]
!              QTRCINP(:,:,:,iq),                                                                                    & ! [IN]
!              FRAIN_MP_RAIN1(0:KA,:,:), FRAIN_MP_SNOW1(0:KA,:,:), FRAIN_MP_RAIN(0:KA,:,:), FRAIN_MP_SNOW(0:KA,:,:), & ! [IN]
!              DENS1(:,:,:), DENS(:,:,:), THRESP                                                                     ) ! [IN]
!    enddo
!
!    
!    !emission
!    !-------------------------------
!    time = TIME_NOWDATE(3) * 1000000 + TIME_NOWDATE(4) * 10000 + TIME_NOWDATE(5) * 100 + TIME_NOWDATE(6) !ddhhmmss
!    
!    EMITF(:,:,:) = 0.0_RP       ![kg/m2/s]
!    mevent = 1
!    if ( PRC_myrank == myrank_FDNPP ) then
!       do iq = 1, QA_AE_AT
!          ! do j = JS, JE-1
!          ! do i = IS, IE-1
!             
!          ! if (        37.42_RP * CONST_D2R >= (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP .and.  37.42_RP * CONST_D2R < (LAT_REAL(i,j)+LAT_REAL(i,j+1))/2.0_RP &
!          !      .and. 141.03_RP * CONST_D2R >= (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP .and. 141.03_RP * CONST_D2R < (LON_REAL(i,j)+LON_REAL(i+1,j))/2.0_RP ) then
!
!          ! write(*,*) iq, i, j, (LAT_REAL(i,j)+LAT_REAL(i,j-1))/2.0_RP/CONST_D2R, (LON_REAL(i,j)+LON_REAL(i-1,j))/2.0_RP/CONST_D2R
!          ! write(*,*) PRC_myrank, myrank_FDNPP
!          ! stop
!
!          do ievent = 1, 66
!             if ( time >= emit_time(ievent,1) .and. time < emit_time(ievent,2) ) then
!                EMITF(I_FDNPP,J_FDNPP,iq) = emit_data(ievent) * 1.0E-14_RP
!                mevent = ievent
!             endif
!          enddo
!          
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 3.0864E-11_RP ! [*1e+14*Bq/m2/s] (3000m*3000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.736E-11_RP ! [*1e+14*Bq/m2/s] (4000m*4000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 1.111E-09_RP ! [*1e+14*Bq/m2/s] (500m*500m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 2.778E-10_RP ! [*1e+14*Bq/m2/s] (1000m*1000m)
!          ! EMITF(i,j,iq) = EMITF(i,j,iq) * 6.944E-11_RP ! [*1e+14*Bq/m2/s] (2000m*2000m)
!          EMITF(I_FDNPP,J_FDNPP,iq) = EMITF(I_FDNPP,J_FDNPP,iq) / ( HAREA(I_FDNPP,J_FDNPP) * 3600.0_RP )  ! [*1e+14*Bq/m2/s]
!          
!          !   endif
!          ! enddo
!          ! enddo
!       enddo
!    endif
!    
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
!          ! QEMIT(i,j,iq) = EMITF(i,j,iq) / DELP(emit_kidx(mevent),i,j) * GRAV * real(dt_AE,kind=RP)
!          do k = KS, KE
!             QEMIT(k,i,j,iq) = emit_coeff(mevent,k) * EMITF(i,j,iq) / DELP(k,i,j) * GRAV * real(dt_AE,kind=RP)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             ! if ( k <= KUP(i,j) ) then
!             ! if ( k == emit_kidx(mevent) ) then
!             !    QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(i,j,iq)
!             ! endif
!             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud
!    !-------------------------------
!    !inside cloud
!    if ( SW_CCN ) then
!       do iq = 1, QA_AE_AT
!          do j = JS, JE
!          do i = IS, IE
!             do k = KS, KE
!                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTIAT(k,i,j) * TCLDF(k,i,j)
!             enddo
!          enddo
!          enddo
!       enddo
!    else
!       do iq = 1, QA_AE_AT
!          do j = JS, JE
!          do i = IS, IE
!             do k = KS, KE
!                QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCAT * TCLDF(k,i,j)
!             enddo
!          enddo
!          enddo
!       enddo
!    endif
!
!    !outside cloud
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCTM1(k,i,j,iq) = QTRC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    
!    !subcloud scavenging (wash out)
!    !-------------------------------
!    ! do j = JS, JE
!    ! do i = IS, IE
!    !    do k = KS, KE
!    !       RAINN_CP     (k,i,j) = FRAIN_CP     (k,i,j) / VTR                / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3                )
!    !       RAINN_MP_RAIN(k,i,j) = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
!    !       RAINN_MP_SNOW(k,i,j) = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
!    !    enddo
!    ! enddo
!    ! enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             TRACN                    = QTRCTM1(k,i,j,iq) / MAT(iq) * DENS(k,i,j)
!             TRACRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADAT(k,i,j,iq) )**2 * abs( VTR                - VTER(k,i,j,iq) ) * TRACN * RAINN_CP     (k,i,j)
!             TRACRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_RAIN(k,i,j)
!             TRACRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_SNOW(k,i,j)
!             QWTDEP_CP     (k,i,j,iq) = TRACRN_CP     (k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_RAIN(k,i,j,iq) = TRACRN_MP_RAIN(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_SNOW(k,i,j,iq) = TRACRN_MP_SNOW(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
!             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
!                QTRCTMP(k,i,j,iq) = 0.0_RP
!                TRACRN (k,i,j,iq) = TRACRN_CP(k,i,j,iq) + TRACRN_MP_RAIN(k,i,j,iq) + TRACRN_MP_SNOW(k,i,j,iq)
!                TRACRN_CP     (k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_RAIN(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_SNOW(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!             endif
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud  scavenging (rain out)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN        (k,i,j,iq) = CLTORN_CP(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             if ( QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq) < 0.0_RP ) then
!                CLTORN_CP     (k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_CP     (k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_RAIN(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_RAIN(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_SNOW(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_SNOW(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                QTRCTM2(k,i,j,iq) = 0.0_RP
!             else
!                QTRCTM2(k,i,j,iq) = QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq)
!             endif
!             QWTDEP_CP     (k,i,j,iq) = QWTDEP_CP     (k,i,j,iq) + CLTORN_CP     (k,i,j,iq)
!             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq)
!             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             QTRCTM1       (k,i,j,iq) = QTRCTMP       (k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !dry deposition
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       if ( iq == 1 ) then
!          do j = JS, JE
!          do i = IS, IE
!             ! CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
!             CDVES(i,j) = CDVE(i,j)
!          enddo
!          enddo
!       endif
!       
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
!            ( QMTX,                    & ! [OUT]
!              DENS, DELP, CDVES, dt_AE ) ! [IN]
!
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
!            ( QTRCTM1(:,:,:,iq),             & ! [MODIFIED]
!              DRYF(:,:,iq),                  & ! [OUT]
!              QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
!    enddo
!
!
!    !gravitational settling
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
!            ( QTRCTM1(:,:,:,iq),              & ! [MODIFIED]
!              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
!              GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
!    enddo
!
!    
!    !in rain(MP), in snow(MP)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + QTRCINR(k,i,j,iq)
!             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + QTRCINS(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !re-emission from rain : CP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
!            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
!              WETF_CP(:,:,iq),                                              & ! [OUT]
!              TRACRN_CP(:,:,:,iq), MAT(iq), CLTORN_CP(:,:,:,iq),            & ! [IN]
!              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
!    enddo
!
!    
!    !re-emission from rain & snow : MP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!       !        GDPM, DELP, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!       !        GDPM, DELP, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!       
!       QTRCINR(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-rain(MP)
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!              GDZM, DELZ, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINR(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINR(k,i,j,iq), 0.0_RP )
!          enddo
!       enddo
!       enddo
!          
!       QTRCINS(:,:,:,iq) = QTRCTM1(:,:,:,iq) !for in-snow(MP)
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!              GDZM, DELZ, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINS(k,i,j,iq) = max( QTRCTM1(k,i,j,iq) - QTRCINS(k,i,j,iq), 0.0_RP )
!          enddo
!       enddo
!       enddo
!
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRCINP(k,i,j,iq) = QTRCINR(k,i,j,iq) + QTRCINS(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!
!       do j = JS, JE
!       do i = IS, IE
!          WETF(i,j,iq) = WETF_CP(i,j,iq) + WETF_MP_RAIN(i,j,iq) + WETF_MP_SNOW(i,j,iq)
!       enddo
!       enddo
!    enddo
!
!
!    !inside cloud + outside cloud 
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRC (k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
!             QLOAD(k,i,j,iq) = QTRC   (k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!             
!    
!    ! !dry deposition
!    ! !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    if ( iq == 1 ) then
!    !       do j = JS, JE
!    !       do i = IS, IE
!    !          CDVES(i,j) = DDV(iq) * CDVE(i,j) / ( DDV(iq) + CDVE(i,j) )
!    !       enddo
!    !       enddo
!    !    endif
!       
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
!    !         ( QMTX,                    & ! [OUT]
!    !           DENS, DELP, CDVES, dt_AE ) ! [IN]
!
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
!    !         ( QTRC(:,:,:,iq),                & ! [MODIFIED]
!    !           DRYF  (:,:,iq),                & ! [OUT]
!    !           QMTX, DENS, DELP, CDVES, dt_AE ) ! [IN]
!    ! enddo
!
!
!    ! !gravitational settling
!    ! !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
!    !         ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
!    !           QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
!    !           GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
!    ! enddo
!
!    
!    !Decay
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) - QTRC(k,i,j,iq) * decay_ratio(iq) * real(dt_AE,kind=RP)
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !================================================
!    !end aerosol transport
!
!    
!    !aerosol
!    !================================================
!    !================================================
!    !end aerosol
!
!    
!    !flux 
!    !================================================
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          TRACDP(i,j,iq) = WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          ACCUM_TOT_DEP        (i,j,iq) = ACCUM_TOT_DEP        (i,j,iq) + TRACDP      (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP        (i,j,iq) = ACCUM_WET_DEP        (i,j,iq) + WETF        (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_CP     (i,j,iq) = ACCUM_WET_DEP_CP     (i,j,iq) + WETF_CP     (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_MP_RAIN(i,j,iq) = ACCUM_WET_DEP_MP_RAIN(i,j,iq) + WETF_MP_RAIN(i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_WET_DEP_MP_SNOW(i,j,iq) = ACCUM_WET_DEP_MP_SNOW(i,j,iq) + WETF_MP_SNOW(i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_DRY_DEP        (i,j,iq) = ACCUM_DRY_DEP        (i,j,iq) + DRYF        (i,j,iq) * real(dt_AE,kind=RP)
!          ACCUM_GRV_DEP        (i,j,iq) = ACCUM_GRV_DEP        (i,j,iq) + GRAVF       (i,j,iq) * real(dt_AE,kind=RP)
!       enddo
!       enddo
!    enddo
!    !================================================
!    !end flux
!
!    
!    !optical parameters
!    !================================================
!    !================================================
!    !end optical parameters
!
!    
!    !OUTPUT
!    !================================================
!    !flux
!    !-------------------------------    
!    do iq = 1, QA_AE_AT
!       write(cnumber,'(i2.2)') iq
!       !                                  i j q
!       call FILE_HISTORY_in( EMITF       (:,:,iq), 'EFTAT'//cnumber//'_SPRINTARS',          'emission flux (tracer'//cnumber//')',                'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( TRACDP      (:,:,iq), 'TDFTAT'//cnumber//'_SPRINTARS',         'total deposition flux (tracer'//cnumber//')',        'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF        (:,:,iq), 'WEFTAT'//cnumber//'_SPRINTARS',         'wet deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_CP     (:,:,iq), 'WEFTAT'//cnumber//'_CP_SPRINTARS',      'wet deposition flux (tracer'//cnumber//';CP)',       'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_MP_RAIN(:,:,iq), 'WEFTAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(RAIN))', 'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( WETF_MP_SNOW(:,:,iq), 'WEFTAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition flux (tracer'//cnumber//';MP(SNOW))', 'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( DRYF        (:,:,iq), 'DRFTAT'//cnumber//'_SPRINTARS',         'dry deposition flux (tracer'//cnumber//')',          'kg/m2/s', dim_type = 'XY' )
!       call FILE_HISTORY_in( GRAVF       (:,:,iq), 'GRFTAT'//cnumber//'_SPRINTARS',         'gravitational settling flux (tracer'//cnumber//')',  'kg/m2/s', dim_type = 'XY' )
!
!       !                                           i j q
!       call FILE_HISTORY_in( ACCUM_TOT_DEP        (:,:,iq), 'ACTDAT'//cnumber//'_SPRINTARS',         'total deposition (tracer'//cnumber//')',                  'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP        (:,:,iq), 'ACWDAT'//cnumber//'_SPRINTARS',         'wet deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_CP     (:,:,iq), 'ACWDAT'//cnumber//'_CP_SPRINTARS',      'wet deposition (tracer'//cnumber//';CP)',                 'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_RAIN(:,:,iq), 'ACWDAT'//cnumber//'_MP_RAIN_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(RAIN))',           'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_WET_DEP_MP_SNOW(:,:,iq), 'ACWDAT'//cnumber//'_MP_SNOW_SPRINTARS', 'wet deposition (tracer'//cnumber//';MP(SNOW))',           'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_DRY_DEP        (:,:,iq), 'ACDDAT'//cnumber//'_SPRINTARS',         'dry deposition (tracer'//cnumber//')',                    'kg/m2', dim_type = 'XY' )
!       call FILE_HISTORY_in( ACCUM_GRV_DEP        (:,:,iq), 'ACGDAT'//cnumber//'_SPRINTARS',         'gravitational settling deposition (tracer'//cnumber//')', 'kg/m2', dim_type = 'XY' )
!    enddo
!       
!    !aerosol
!    !-------------------------------
!
!    !optical properties
!    !-------------------------------
!    
!    !clear sky diagnoses
!    !-------------------------------
!
!    !================================================
!    !end OUTPUT
!
!    return
!  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_SILICATE_PROXY
!  !-----------------------------------------------------------------------------
!  !-----------------------------------------------------------------------------
!  !> SPRINTARS Tracer (cedar pollen)
!  subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_CEDAR_POLLEN &
!       ( QTRC,          & ! [MODIFIED]
!         CDVE,          & ! [IN]
!         RH,            & ! [IN]
!         DENS,          & ! [IN]
!         FRAIN,         & ! [IN]
!         FRAIN_CP,      & ! [IN]
!         FRAIN_MP_RAIN, & ! [IN]
!         FRAIN_MP_SNOW, & ! [IN]
!         VTR_MP_RAIN,   & ! [IN]
!         VTR_MP_SNOW,   & ! [IN]
!         RADR_MP_RAIN,  & ! [IN]
!         RADR_MP_SNOW,  & ! [IN]
!         TCLDF,         & ! [IN]
!         KUP,           & ! [IN]
!         KUPMAX,        & ! [IN]
!         DPKUP,         & ! [IN]
!         GDPM,          & ! [IN]
!         DELP,          & ! [IN]
!         GDZM,          & ! [IN]
!         DELZ,          & ! [IN]
!         RNWR,          & ! [IN]
!         RNWR_CP,       & ! [IN]
!         RNWR_MP_RAIN,  & ! [IN]
!         RNWR_MP_SNOW,  & ! [IN]
!         TIME_NOWDATE,  & ! [IN]
!         dt_AE,         & ! [IN]
!         QA_AE_AT       ) ! [IN]
!    implicit none
!
!    ! input
!    !-------------------------
!    real(DP), intent(in)    :: dt_AE
!    integer,  intent(in)    :: QA_AE_AT
!    
!    real(RP), intent(in)    :: CDVE              (IA,JA) !< bulk coef. * Uabs = Ustar * Qstar / (QV(KS)-SFC_QVsat)
!    real(RP), intent(in)    :: RH           (KA,  IA,JA) !< relative humidity
!    real(RP), intent(in)    :: DENS         (KA,  IA,JA) !< air density
!    real(RP), intent(in)    :: FRAIN        (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP + MP
!    real(RP), intent(in)    :: FRAIN_CP     (0:KA,IA,JA) !< precipitation flux (liq+ice) : CP
!    real(RP), intent(in)    :: FRAIN_MP_RAIN(0:KA,IA,JA) !< precipitation flux (liq    ) : MP
!    real(RP), intent(in)    :: FRAIN_MP_SNOW(0:KA,IA,JA) !< precipitation flux (    ice) : MP
!    real(RP), intent(in)    :: VTR_MP_RAIN  (KA,  IA,JA) !< terminal velocity of rain (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: VTR_MP_SNOW  (KA,  IA,JA) !< terminal velocity of snow (microphysic scheme) [m/s]
!    real(RP), intent(in)    :: RADR_MP_RAIN (KA,  IA,JA) !< effective radius of rain (microphysic scheme) [m]
!    real(RP), intent(in)    :: RADR_MP_SNOW (KA,  IA,JA) !< effective radius of snow (microphysic scheme) [m]
!    real(RP), intent(in)    :: TCLDF        (KA,  IA,JA) !< cloud fraction (CP+MP)
!    integer,  intent(in)    :: KUP               (IA,JA) !< uplift level
!    integer,  intent(in)    :: KUPMAX                    !< maximum uplift level
!    real(RP), intent(in)    :: DPKUP             (IA,JA) !< GDPM(KS-1,i,j) - GDPM(KUP(i,j),i,j)
!    real(RP), intent(in)    :: GDPM         (0:KA,IA,JA) !< pressure at half levels
!    real(RP), intent(in)    :: DELP         (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
!    real(RP), intent(in)    :: GDZM         (0:KA,IA,JA) !< altitude at half levels
!    real(RP), intent(in)    :: DELZ         (KA,  IA,JA) !< GDZM(k,i,j)-GDZM(k-1,i,j)
!    real(RP), intent(in)    :: RNWR         (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_CP      (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_RAIN (KA,  IA,JA) !< ratio of rain to cloud
!    real(RP), intent(in)    :: RNWR_MP_SNOW (KA,  IA,JA) !< ratio of rain to cloud
!    integer,  intent(in)    :: TIME_NOWDATE(6)           !< current time
!    
!    ! modified
!    !-------------------------
!    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE_AT)
!    
!    ! output
!    !-------------------------
!
!    ! internal work
!    !-------------------------
!    !set parameter
!    real(RP) :: MAT           (QA_AE_AT)
!    real(RP) :: RADAT(KA,IA,JA,QA_AE_AT) !< tracer particle radius      [m]
!    real(RP) :: BRH1, BRH2, BRH3         !<
!    real(RP) :: VTER (KA,IA,JA,QA_AE_AT) !< terminal velocity of tracer [m/s]
!
!    !emission
!    real(RP) :: QEMIT(IA,JA,QA_AE_AT) !emission mixing ratio [kg/kg]
!    real(RP) :: EMITF(IA,JA,QA_AE_AT) !tracer emission flux  [kg/m2/s]
!    
!    !incloud
!    real(RP) :: QTRCTM1(KA,IA,JA,QA_AE_AT)  !< outside cloud
!    real(RP) :: QTRCTM2(KA,IA,JA,QA_AE_AT)  !< inside  cloud
!
!    !subcloud scavenging (wash out)
!    real(RP) :: RAINN        (KA,IA,JA)     !< rain number concentrarion        [1/m3]
!    real(RP) :: RAINN_CP     (KA,IA,JA)               
!    real(RP) :: RAINN_MP_RAIN(KA,IA,JA)              
!    real(RP) :: RAINN_MP_SNOW(KA,IA,JA)               
!    real(RP) :: TRACN                             !< number density of tracer [1/m3]
!    real(RP) :: TRACRN        (KA,IA,JA,QA_AE_AT) !< number density of tracer entering into raindrops [1/s/m3] <=> collision frequency between tracer and rain
!    real(RP) :: TRACRN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: TRACRN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP        (KA,IA,JA,QA_AE_AT) !< mixing ratio by wet deposition   [kg/kg] 
!    real(RP) :: QWTDEP_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QWTDEP_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    real(RP) :: QTRCTMP       (KA,IA,JA,QA_AE_AT) !< tracers that do not cllide with rain particles [kg/kg]
!    
!    !incloud scavenging (rain out)
!    real(RP) :: CLTORN        (KA,IA,JA,QA_AE_AT) !< mixing ratio moving cloud to rain particles [kg/kg]
!    real(RP) :: CLTORN_CP     (KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_RAIN(KA,IA,JA,QA_AE_AT)
!    real(RP) :: CLTORN_MP_SNOW(KA,IA,JA,QA_AE_AT)
!    
!    !re-emission from rain
!    real(RP) :: WETF        (IA,JA,QA_AE_AT)      !< wet deposition flux
!    real(RP) :: WETF_CP     (IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_RAIN(IA,JA,QA_AE_AT)
!    real(RP) :: WETF_MP_SNOW(IA,JA,QA_AE_AT)
!    
!    !dry deposition
!    real(RP) :: QMTX   (KA,IA,JA,-1:1)      !< implicit matrix of QTRC
!    real(RP) :: DRYF      (IA,JA,QA_AE_AT)  !< dry deposition flux
!    
!    !gravitational settling
!    real(RP) :: QLOAD  (KA,IA,JA,QA_AE_AT)  !< mixing ratio for radiation [kg/kg]
!    real(RP) :: GRAVF     (IA,JA,QA_AE_AT)  !< gravitational settling flux 
!    
!    !flux
!    real(RP) :: TRACDP    (IA,JA,QA_AE_AT)  !< total deposition
!    
!    ! parameter
!    !-------------------------
!    real(RP), parameter :: RADTRAC(NAT) = (/ 15.000E-6_RP /) !< particle radius         [m]
!    real(RP), parameter :: DENSAT (NAT) = (/  1.000E+3_RP /) !< tracer particle density [kg/m3]
!    real(RP), parameter :: AA     (NAT) = (/  1.000E-2_RP /) !< collision efficiency
!    ! real(RP), parameter :: SIGMA  (NAT) = (/  1.526E+0_RP /) !< GSD of log-normal distribution
!    
!    ! others 
!    !-------------------------
!    integer   :: k, i, j, iq
!    integer   :: IRH
!    character :: cnumber*2
!    
!    !----------------------------------------------------
!    
!    LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - tracer'
!
!    !setup parameter 
!    !================================================
!    !RADAT
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       RADAT(:,:,:,iq) = RADTRAC(iq)
!    enddo
!    
!    !MAT
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       MAT(iq) = 4.0_RP / 3.0_RP * PI * RADTRAC(1)**3 * DENSAT(iq)
!    enddo
!    
!    !VTER
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             VTER(k,i,j,iq) = DENSAT(iq) * 2.0_RP * RADAT(k,i,j,iq)**2 * GRAV / 9.0_RP / NU
!          enddo
!       enddo
!       enddo
!    enddo
!    
!    !================================================
!    !end setup parameter
!
!    
!    !aerosol transport 
!    !================================================
!    !emission
!    !-------------------------------
!    EMITF(:,:,:) = 0.0_RP
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          QEMIT(i,j,iq) = EMITF(i,j,iq) / DPKUP(i,j) * GRAV * real(dt_AE,kind=RP)
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          if ( k <= KUP(i,j) ) then
!             QTRC(k,i,j,iq) = QTRC(k,i,j,iq) + QEMIT(i,j,iq)
!          endif
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud
!    !-------------------------------
!    ! do iq = 1, QA_AE_AT
!    !    do j = JS, JE
!    !    do i = IS, IE
!    !       do k = KS, KE
!    !          if ( SW_CCN ) then
!    !             QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * ACTIAT(k,i,j) * TCLDF(k,i,j)
!    !          else
!    !             QTRCTM2(k,i,j,iq) = QTRC(k,i,j,iq) * FINCAT * TCLDF(k,i,j)
!    !          endif
!    !          QTRCTM1(k,i,j,iq) = QTRC(k,i,j,iq) - QTRCTM2(k,i,j,iq)
!    !       enddo
!    !    enddo
!    !    enddo
!    ! enddo
!    QTRCTM2(:,:,:,:) = 0.0_RP
!    QTRCTM1(:,:,:,:) = QTRC(:,:,:,:)
!
!    
!    !subcloud scavenging (wash out)
!    !-------------------------------
!    do j = JS, JE
!    do i = IS, IE
!       do k = KS, KE
!          RAINN_CP     (k,i,j) = FRAIN_CP     (k,i,j) / VTR                / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR**3                )
!          RAINN_MP_RAIN(k,i,j) = FRAIN_MP_RAIN(k,i,j) / VTR_MP_RAIN(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_RAIN(k,i,j)**3 )
!          RAINN_MP_SNOW(k,i,j) = FRAIN_MP_SNOW(k,i,j) / VTR_MP_SNOW(k,i,j) / DENSR / ( 4.0_RP / 3.0_RP * PI * RADR_MP_SNOW(k,i,j)**3 )
!       enddo
!    enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             TRACN                    = QTRCTM1(k,i,j,iq) / MAT(iq) * DENS(k,i,j)
!             TRACRN_CP     (k,i,j,iq) = AA(iq) * PI * ( RADR                + RADAT(k,i,j,iq) )**2 * abs( VTR                - VTER(k,i,j,iq) ) * TRACN * RAINN_CP     (k,i,j)
!             TRACRN_MP_RAIN(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_RAIN(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_RAIN(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_RAIN(k,i,j)
!             TRACRN_MP_SNOW(k,i,j,iq) = AA(iq) * PI * ( RADR_MP_SNOW(k,i,j) + RADAT(k,i,j,iq) )**2 * abs( VTR_MP_SNOW(k,i,j) - VTER(k,i,j,iq) ) * TRACN * RAINN_MP_SNOW(k,i,j)
!             QWTDEP_CP     (k,i,j,iq) = TRACRN_CP     (k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_RAIN(k,i,j,iq) = TRACRN_MP_RAIN(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QWTDEP_MP_SNOW(k,i,j,iq) = TRACRN_MP_SNOW(k,i,j,iq) * MAT(iq) / DENS(k,i,j) * real(dt_AE,kind=RP)
!             QTRCTMP(k,i,j,iq) = QTRCTM1(k,i,j,iq) - QWTDEP_CP(k,i,j,iq) - QWTDEP_MP_RAIN(k,i,j,iq) - QWTDEP_MP_SNOW(k,i,j,iq)
!             if ( QTRCTMP(k,i,j,iq) < 0.0_RP ) then
!                QTRCTMP(k,i,j,iq) = 0.0_RP
!                TRACRN (k,i,j,iq) = TRACRN_CP(k,i,j,iq) + TRACRN_MP_RAIN(k,i,j,iq) + TRACRN_MP_SNOW(k,i,j,iq)
!                TRACRN_CP     (k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_RAIN(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                TRACRN_MP_SNOW(k,i,j,iq) = TRACN / real(dt_AE,kind=RP) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_CP     (k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_CP     (k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_RAIN(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_RAIN(k,i,j,iq) / TRACRN(k,i,j,iq)
!                QWTDEP_MP_SNOW(k,i,j,iq) = QTRCTM1(k,i,j,iq) * TRACRN_MP_SNOW(k,i,j,iq) / TRACRN(k,i,j,iq)
!             endif
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !incloud  scavenging (rain out)
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             CLTORN_CP     (k,i,j,iq) = RNWR_CP     (k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_RAIN(k,i,j,iq) = RNWR_MP_RAIN(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN_MP_SNOW(k,i,j,iq) = RNWR_MP_SNOW(k,i,j) * QTRCTM2(k,i,j,iq)
!             CLTORN        (k,i,j,iq) = CLTORN_CP(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             if ( QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq) < 0.0_RP ) then
!                CLTORN_CP     (k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_CP     (k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_RAIN(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_RAIN(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                CLTORN_MP_SNOW(k,i,j,iq) = QTRCTM2(k,i,j,iq) * RNWR_MP_SNOW(k,i,j) / ( RNWR_CP(k,i,j) + RNWR_MP_RAIN(k,i,j) + RNWR_MP_SNOW(k,i,j) )
!                QTRCTM2(k,i,j,iq) = 0.0_RP
!             else
!                QTRCTM2(k,i,j,iq) = QTRCTM2(k,i,j,iq) - CLTORN(k,i,j,iq)
!             endif
!             QWTDEP_CP     (k,i,j,iq) = QWTDEP_CP     (k,i,j,iq) + CLTORN_CP     (k,i,j,iq)
!             QWTDEP_MP_RAIN(k,i,j,iq) = QWTDEP_MP_RAIN(k,i,j,iq) + CLTORN_MP_RAIN(k,i,j,iq)
!             QWTDEP_MP_SNOW(k,i,j,iq) = QWTDEP_MP_SNOW(k,i,j,iq) + CLTORN_MP_SNOW(k,i,j,iq)
!             QTRCTM1       (k,i,j,iq) = QTRCTMP       (k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!
!    
!    !re-emission from rain : CP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
!            ( QTRCTM1(:,:,:,iq), QWTDEP_CP(:,:,:,iq),                       & ! [MODIFIED]
!              WETF_CP(:,:,iq),                                              & ! [OUT]
!              TRACRN_CP(:,:,:,iq), MAT(iq), CLTORN_CP(:,:,:,iq),            & ! [IN]
!              DENS(:,:,:), FRAIN_CP(:,:,:), DELP(:,:,:), DELZ(:,:,:), dt_AE ) ! [IN]
!    enddo
!
!    !re-emission from rain & snow : MP
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!              GDPM, DELP, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
!            ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!              WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!              QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!              GDPM, DELP, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_RAIN(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_RAIN(:,:,:,iq),                   & ! [IN]
!       !        GDZM, DELZ, DENS, VTR_MP_RAIN(:,:,:), dt_AE ) ! [IN]
!       ! call ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
!       !      ( QTRCTM1(:,:,:,iq),                          & ! [MODIFIED]
!       !        WETF_MP_SNOW(:,:,iq),                       & ! [OUT]
!       !        QWTDEP_MP_SNOW(:,:,:,iq),                   & ! [IN]
!       !        GDZM, DELZ, DENS, VTR_MP_SNOW(:,:,:), dt_AE ) ! [IN]
!       
!       do j = JS, JE
!       do i = IS, IE
!          WETF(i,j,iq) = WETF_CP(i,j,iq) + WETF_MP_RAIN(i,j,iq) + WETF_MP_SNOW(i,j,iq)
!       enddo
!       enddo
!    enddo
!
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          do k = KS, KE
!             QTRC(k,i,j,iq) = QTRCTM1(k,i,j,iq) + QTRCTM2(k,i,j,iq)
!          enddo
!       enddo
!       enddo
!    enddo
!             
!    
!    !dry deposition
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
!            ( QMTX,                   & ! [OUT]
!              DENS, DELP, CDVE, dt_AE ) ! [IN]
!
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
!            ( QTRC(:,:,:,iq),               & ! [MODIFIED]
!              DRYF  (:,:,iq),               & ! [OUT]
!              QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]
!    enddo
!
!
!    !gravitational settling
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
!            ( QTRC(:,:,:,iq),                 & ! [MODIFIED]
!              QLOAD(:,:,:,iq), GRAVF(:,:,iq), & ! [OUT]
!              GDPM, DELP, DENS, VTER, dt_AE   ) ! [IN]
!    enddo
!
!    !================================================
!    !end aerosol transport
!
!    
!    !aerosol
!    !================================================
!    !================================================
!    !end aerosol
!
!    
!    !flux 
!    !================================================
!    do iq = 1, QA_AE_AT
!       do j = JS, JE
!       do i = IS, IE
!          TRACDP(i,j,iq) = WETF(i,j,iq) + DRYF(i,j,iq) + GRAVF(i,j,iq)
!       enddo
!       enddo
!    enddo
!    !================================================
!    !end flux
!
!    
!    !optical parameters
!    !================================================
!    !================================================
!    !end optical parameters
!
!    
!    !OUTPUT
!    !================================================
!    !flux
!    !-------------------------------
!    do iq = 1, QA_AE_AT
!       write(cnumber,'(i2.2)') iq
!       !                            i j q
!       call FILE_HISTORY_in( EMITF (:,:,iq), 'EFTAT'//cnumber//'_SPRINTARS',  'emission flux (tracer'//cnumber//')',               'kg/m2/s', dim_type = 'XY'  )
!       call FILE_HISTORY_in( TRACDP(:,:,iq), 'TDFTAT'//cnumber//'_SPRINTARS', 'total deposition flux (tracer'//cnumber//')',       'kg/m2/s', dim_type = 'XY'  )
!       call FILE_HISTORY_in( WETF  (:,:,iq), 'WEFTAT'//cnumber//'_SPRINTARS', 'wet deposition flux (tracer'//cnumber//')',         'kg/m2/s', dim_type = 'XY'  )
!       call FILE_HISTORY_in( DRYF  (:,:,iq), 'DRFTAT'//cnumber//'_SPRINTARS', 'dry deposition flux (tracer'//cnumber//')',         'kg/m2/s', dim_type = 'XY'  )
!       call FILE_HISTORY_in( GRAVF (:,:,iq), 'GRFTAT'//cnumber//'_SPRINTARS', 'gravitational settling flux (tracer'//cnumber//')', 'kg/m2/s', dim_type = 'XY'  )
!    enddo
!       
!    !aerosol
!    !-------------------------------
!
!    !optical properties
!    !-------------------------------
!    
!    !clear sky diagnoses
!    !-------------------------------
!
!    !================================================
!    !end OUTPUT
!
!    
!    return
!  end subroutine ATMOS_PHY_AE_SPRINTARS_AEROAT_CEDAR_POLLEN
  !-----------------------------------------------------------------------------
end module scale_atmos_phy_ae_SPRINTARS_aerotrac
