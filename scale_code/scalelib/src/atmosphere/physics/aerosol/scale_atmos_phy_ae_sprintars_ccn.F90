!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - CCN
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_ccn
  !-----------------------------------------------------------------------------
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
       PRC_IsMaster, &
       PRC_abort
  
  use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE

  use scale_const, only: &
       PI    => CONST_PI,     & !< pi
       GRAV  => CONST_GRAV,   & !< standard acceleration of gravity          [m/s2] = 9.80665_RP
       P00   => CONST_Pstd,   & !< standard pressure                         [Pa]
       RVAP  => CONST_Rvap,   & !< specific gas constant (water vapor)       [J/kg/K]
       CP    => CONST_CPdry,  & !< specific heat (dry air,constant pressure) [J/kg/K]
       EL    => CONST_LHV0,   & !< latent heat of vaporizaion at 0C          [J/kg] = 2.501E+6_RP
       EPSV  => CONST_EPSvap, & !< Rdry / Rvap
       EMELT => CONST_EMELT,  & !< latent heat of melting
       TMELT => CONST_TMELT,  & !< melting point of water
       DWATR => CONST_DWATR,  & !< density of water [kg/m3]
       PSAT0 => CONST_PSAT0,  & !< saturate pressure of water vapor at 0C [Pa] = 610.78_RP  
       CONST_UNDEF
  
  use scale_file_history, only: &
       FILE_HISTORY_in

  use scale_atmos_saturation, only: &
       ATMOS_SATURATION_dens2qsat_all, &
       ATMOS_SATURATION_psat_all, &
       ATMOS_SATURATION_dqs_dtem_dens_liq, &
       ATMOS_SATURATION_dqs_dtem_dens_ice, &
       ATMOS_SATURATION_alpha, &
       ATMOS_SATURATION_dalphadT

  ! use scale_atmos_phy_ae_SPRINTARS_common, only: &
  !      VMISS                  !< missing value
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_setup
  public :: ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN
  
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
  ! variables
  ! namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN
  real(RP), private, save :: &
       SIGAE(6) &                                    !< standard deviation of aerosol size distribution
       !    soil dust   ext.BC      int.OC&BC   SOC         sulfate     sea salt
       = (/ 2.20E+0_RP, 2.00E+0_RP, 2.24E+0_RP, 2.24E+0_RP, 2.03E+0_RP, 2.03E+0_RP /)
  real(RP), private, save :: &
       HYGRAE(6) &                                   !< aerosol hygroscopic coefficient
       !    soil dust   ext.BC      int.OC&BC   SOC         sulfate     sea salt
       = (/ 0.14E+0_RP, 5.00E-7_RP, 0.14E+0_RP, 0.14E+0_RP, 0.51E+0_RP, 1.16E+0_RP /)
  real(RP), private, save :: SC      =   0.70E+00_RP !< coefficient of turbulent updraft velocity
  real(RP), private, save :: OMGMIN  =   1.00E-01_RP !< minimum updraft velocity [m/s]
  real(RP), private, save :: TKEMIN  =   0.00E+00_RP !< minimum TKE
  real(RP), private, save :: TKEMAX  = 100.00E+00_RP !< maximum TKE
  real(RP), private, save :: NAEMIN  =   1.00E+06_RP !< minimum aerosol number concentration [1/m3]
  real(RP), private, save :: CLWMIN  =   1.00E-06_RP !< minimum cloud water for cloud base diagnosis
  real(RP), private, save :: DFEMIN  =   1.00E-05_RP
  real(RP), private, save :: EPS     =   1.00E-20_RP
  real(RP), private, save :: EPSC    =   1.00E-06_RP
  ! real(RP), private, save :: EPSQ    =   1.00E-08_RP !< minimum for qc/C
  ! real(RP), private, save :: THMNC   = 238.15E+00_RP !< temperature of homogeneous freezing
  ! real(RP), private, save :: TCTFR   = 270.15E+00_RP !< temperature of contact freezing
  ! real(RP), private, save :: &
  !      PRCCOL(3) &                                   !< accretion coefficient
  !      = (/ 1.0E+0_RP, 5.0E-1_RP, 5.0E-2_RP /)
  ! real(RP), private, save :: TNUW    =   1.44E+04_RP !< time scale of cloud droplet nucleated
  ! real(RP), private, save :: IMNIC   =   1.00E-12_RP !< initial mass of a nucleated ice crystal [kg]
  ! real(RP), private, save :: ADIF    =   1.40E-08_RP !< aerosol diffusivity [1/m2/s]
  ! real(RP), private, save :: NAO     =   2.00E+05_RP !< concentration of active contact nuclei [1/m3]
  ! real(RP), private, save :: SFE     =   1.00E+00_RP !< snow formation efficiency
  ! real(RP), private, save :: ICCE    =   0.10E+00_RP !< collision efficiency between ice crystals
  ! real(RP), private, save :: FVSICD  =   0.25E+00_RP !< dispersion of ice crystal fall velocity spectrum
  ! real(RP), private, save :: DICE    =   5.00E+02_RP !< density of ice crystal
  ! real(RP), private, save :: EMPC    =   7.00E+02_RP !< empilical constant for aggregation
  ! real(RP), private, save :: MI0     =   1.00E-12_RP !< mass of newly formed ice crystal
  ! real(RP), private, save :: FBCHY   =   1.00E-01_RP !< ratio of hydrophilic BC
  ! real(RP), private, save :: TWSNOWA = 275.15E+00_RP !< temperature of rain/snow partition
  ! real(RP), private, save :: FHMP    =   3.00E-03_RP

  ! ! namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY
  ! real(RP), private, save :: B1 = 0.01E+00_RP        !< Berry coefficient
  ! real(RP), private, save :: B2 = 0.12E+00_RP        !< Berry coefficient
  ! real(RP), private, save :: B3 = 1.00E-12_RP        !< Berry coefficient
  
  ! ! namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL
  ! real(RP), private, save :: &
  !      RCMIN(KCPCL) &                               !< minimum of radius (DUMMY)
  !      = (/  1.0E-6_RP,   1.0E-6_RP /)
  ! real(RP), private, save :: &
  !      RCMAX(KCPCL) &                               !< maximum of radius (DUMMY)
  !      = (/ 30.0E-6_RP, 200.0E-6_RP /)
  ! real(RP), private, save :: &
  !      REFCT(KCPCL) &                               !< Re/Rvol   (DUMMY)
  !      = (/  1.1E+0_RP,   1.0E+0_RP /)
  ! real(RP), private, save :: TFLIQ   =  13.50E+0_RP !< parameter for FLIQ        (DUMMY)
  ! real(RP), private, save :: SFLIQ   = 269.91E+0_RP !< parameter for FLIQ        (DUMMY)
  ! real(RP), private, save :: MFLIQ   = 235.15E+0_RP !< parameter for FLIQ        (DUMMY)
  ! real(RP), private, save :: FLIQMIN =   2.00E-2_RP !< minimum liquid fraction (for CCN generation)
  ! real(RP), private       :: UAMIN                  !< minimum aerosol number    (DUMMY)
  ! real(RP), private       :: CCFCT                  !< CDCN/CCN                  (DUMMY)
  ! real(RP), private, save :: &
  !      UCMIN(KCPCL) &                               !< minimum cloud number concentration
  !      = (/ 1.0E+07_RP, 1.0E+01_RP /)
  ! real(RP), private, save :: &
  !      UCMAX(KCPCL) &                               !< maximum cloud drop number (DUMMY)
  !      = (/ 1.0E+14_RP, 1.0E+14_RP /)
  ! logical,  private, save :: OCCN2CLD = .true.      !< use prognoced ice fraction
  ! logical,  private, save :: OCLD2CCN = .true.      !< use prognoced ice fraction
  real(RP), private, allocatable :: zero_3d(:,:,:)
  real(RP), private, allocatable :: zero_2d(:,:)
  character(len=7), private :: cae(6)
  data cae / 'DUST', 'EXTBC', 'INTOCBC', 'SOC', 'SULFATE', 'SALT' /
  
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS AEROCCN setup
  subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_setup
    implicit none
    integer :: ierr, iq
    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN / &
         SIGAE,  HYGRAE,  SC,     OMGMIN, TKEMIN, &
         TKEMAX, NAEMIN,  CLWMIN, DFEMIN, EPS,    &
         EPSC
         !EPSQ,    THMNC,  TCTFR,  PRCCOL, &
         !TNUW,   IMNIC,   ADIF,   NAO,    SFE,    &
         !ICCE,   FVSICD,  DICE,   EMPC,   MI0,    &
         !FBCHY,  TWSNOWA, FHMP
    ! namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY / &
    !      B1, B2, B3
    ! namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL / &
    !      RCMIN, RCMAX,    REFCT,   TFLIQ, SFLIQ, &
    !      MFLIQ, FLIQMIN,  UAMIN ,  CCFCT, UCMIN, &
    !      UCMAX, OCCN2CLD, OCLD2CCN
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*) 'Setup'

    !--- read namelist
    !PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMACCN)

    ! !PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY
    ! rewind( IO_FID_CONF )
    ! read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY, iostat=ierr )
    ! if( ierr < 0 ) then !--- missing
    !    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY. Default used.'
    ! elseif( ierr > 0 ) then !--- fatal error
    !    LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY, Check!'
    !    call PRC_abort
    ! endif
    ! LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMBERY)

    ! !PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL
    ! rewind( IO_FID_CONF )
    ! read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL, iostat=ierr )
    ! if( ierr < 0 ) then !--- missing
    !    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL. Default used.'
    ! elseif( ierr > 0 ) then !--- fatal error
    !    LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_ccn",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL, Check!'
    !    call PRC_abort
    ! endif
    ! LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_AEROCCN_NMAECL)

    allocate( zero_3d(KA,IA,JA) )
    allocate( zero_2d(IA,JA) )
    zero_3d(:,:,:) = 0.0_RP
    zero_2d(:,:) = 0.0_RP
    do iq = 1, 6
       !                            k i j
       call FILE_HISTORY_in( zero_3d(:,:,:), 'NMAE_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'aerosol number concentration ('//trim(adjustl(cae(iq)))//')',                 '1/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_3d(:,:,:), 'NCCN_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'activated aerosol (CCN) number concentration ('//trim(adjustl(cae(iq)))//')', '1/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( zero_3d(:,:,:), 'ACTI_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'ratio of activated aerosol (CCN) ('//trim(adjustl(cae(iq)))//')',             '1',    dim_type = 'ZXY' )
    enddo
    !                            k i j
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NMAET_SPRINTARS',  'total aerosol number concentration',           '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'NCCN_SPRINTARS',   'activated aerosol (CCN) number concentration', '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CCN01_SPRINTARS',  'CCN(0.1%)',                                    '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'CCN03_SPRINTARS',  'CCN(0.3%)',                                    '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_2d(:,:),   'CCN01C_SPRINTARS', 'CCN(0.1%) (in column)',                        '1/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'CCN03C_SPRINTARS', 'CCN(0.3%) (in column)',                        '1/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_2d(:,:),   'UAPCB_SPRINTARS',  'CCN number at cloud base',                     '1/m3', dim_type = 'XY'  )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'OMGCCN_SPRINTARS', 'updraft velocity',                             'm/s',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'AVRW_SPRINTARS',   'updraft velocity (average)',                   'm/s',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( zero_3d(:,:,:), 'TURW_SPRINTARS',   'updraft velocity (turbulence)',                'm/s',  dim_type = 'ZXY' )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_setup
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> SPRINTARS AEROCCN finalize
  subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_finalize
    implicit none

    deallocate( zero_3d )
    deallocate( zero_2d )
  
  end subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN_finalize

  !-----------------------------------------------------------------------------
  !> SPRINTARS AEROCCN (prognosis/diagnosis of cloud water/ice number concentration)
  subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN &
       ( UAPCL,  & ! [OUT]
         NUMAET, & ! [OUT]
         ACTI,   & ! [OUT]
         TCLDF,  & ! [IN]
         TCLDF2, & ! [IN]
         QC,     & ! [IN]
         NUMAE,  & ! [INOUT]
         RADAE,  & ! [IN]
         DELZ,   & ! [IN]
         DENS,   & ! [IN]
         GDW,    & ! [IN]
         TKE,    & ! [IN]
         GDT,    & ! [IN]
         GDP     ) ! [IN]
    implicit none

    ! output
    !-------------------------
    real(RP), intent(out)   :: UAPCL (KA,IA,JA)   !< CCN number concentration    [1/m3]
    real(RP), intent(out)   :: NUMAET(KA,IA,JA)   !< total aerosol number concentration [1/m3]
    real(RP), intent(out)   :: ACTI  (KA,IA,JA,6) !< ratio of activated aerosols [1]
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: TCLDF (KA,IA,JA)   !< cloud fraction (total)
    real(RP), intent(in)    :: TCLDF2(KA,IA,JA)   !< cloud fraction (total)
    real(RP), intent(in)    :: QC    (KA,IA,JA)   !< cloud water                  [kg/kg]
    real(RP), intent(inout) :: NUMAE (KA,IA,JA,6) !< aerosol number concentration [1/m3]
    real(RP), intent(in)    :: RADAE (KA,IA,JA,6) !< aerosol radius               [m]
    real(RP), intent(in)    :: DELZ  (KA,IA,JA)   !< GDZM(k,i,j)-GDZM(k-1,i,j)
    real(RP), intent(in)    :: DENS  (KA,IA,JA)   !< air density
    real(RP), intent(in)    :: GDW   (KA,IA,JA)   !< vertical velocity            [m/s]
    real(RP), intent(in)    :: TKE   (KA,IA,JA)   !< turbulent kinetic energy     [m2/s2]
    real(RP), intent(in)    :: GDT   (KA,IA,JA)   !< temperature
    real(RP), intent(in)    :: GDP   (KA,IA,JA)   !< pressure
    
    ! internal work
    !-------------------------
    !setup parameter 
    real(RP) :: F1(6)               !< 
    real(RP) :: F2(6)               !< 
    real(RP) :: SB(6)               !< 
    real(RP) :: CURV   (KA,IA,JA)   !< curvature effect
    real(RP) :: SATQ   (KA,IA,JA)   !< saturation water vapour mixing ratio
    real(RP) :: SATQT               !< d(QSAT)/d(T)
    real(RP) :: ALFA   (KA,IA,JA)   !< coefficient
    real(RP) :: BETA   (KA,IA,JA)   !< coefficient
    real(RP) :: BG     (KA,IA,JA)   !< coefficient
    real(RP) :: DFF                 !<
    real(RP) :: AVRW   (KA,IA,JA)   !< updraft velocity
    real(RP) :: TURW   (KA,IA,JA)   !< updraft velocity (turbulent)
    real(RP) :: OMGCCN (KA,IA,JA)   !< updraft velocity
    real(RP) :: FJ_term1
    real(RP) :: FJ_term2
    real(RP) :: FJ     (KA,IA,JA)   !< 
    real(RP) :: FJ01   (KA,IA,JA)   !< 
    real(RP) :: FJ03   (KA,IA,JA)   !<

    real(RP) :: qsat_liq_, qsat_ice_ 
    real(RP) :: alpha_, dalpha_dT
    real(RP) :: dqsat_dt_liq, dqsat_dT_ice 

    !CCN
    real(RP) :: UAPCL1 (KA,IA,JA)   !< cloud droplet number concentration
    real(RP) :: UAPCL3 (KA,IA,JA)   !< cloud droplet number concentration
    real(RP) :: UAPCL1C   (IA,JA)   !< cloud droplet number concentration
    real(RP) :: UAPCL3C   (IA,JA)   !< cloud droplet number concentration
    real(RP) :: UAPCB     (IA,JA)   !< cloud base droplet number concentration
    real(RP) :: SM     (KA,IA,JA,6) !< 
    real(RP) :: ACTI01 (KA,IA,JA,6) !< ratio of activated aerosols
    real(RP) :: ACTI03 (KA,IA,JA,6) !< ratio of activated aerosols
    real(RP) :: NUMCLD (KA,IA,JA,6) !< 
    real(RP) :: NCLD01 (KA,IA,JA,6) !< 
    real(RP) :: NCLD03 (KA,IA,JA,6) !<
    real(RP) :: TCLDWC1
      
    ! parameter
    !-------------------------
    real(RP),  parameter :: STENSW = 7.562E-2_RP !< surface tension of water against air [kg/s2]
    real(RP),  parameter :: CONDE  = 2.41E-2_RP  !< thermal conductivity [W/m/K]
!    character, parameter :: cae(6)*7 = (/ 'DUST', 'EXTBC', 'INTOCBC', 'SOC', 'SULFATE', 'SALT' /)
!    character(len=7) :: cae(6)
!    data cae / 'DUST', 'EXTBC', 'INTOCBC', 'SOC', 'SULFATE', 'SALT' /

    ! others
    !-------------------------
    integer :: k, i, j, iq
    
    !----------------------------------------------------

    !LOG_PROGRESS(*) 'atmosphere / physics / SPRINTARS - ccn'
    
    !setup parameter 
    !================================================
    !F1, F2, SB
    !-------------------------------
    do iq = 1, 6
       F1(iq) = exp( 3.16_RP * (log(SIGAE(iq)))**1.97_RP )
       F2(iq) = 0.382_RP * log(SIGAE(iq))
       SB(iq) = 4.0_RP / ( 3.0_RP * sqrt(2.0_RP*PI) * log(SIGAE(iq)) )
    enddo

    !NUMAE, NUMAET
    !-------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          do iq = 1, 6
             NUMAE(k,i,j,iq) = max( NUMAE(k,i,j,iq), 0.0_RP )
          enddo
          NUMAET(k,i,j) = sum( NUMAE(k,i,j,1:6) )
          if ( NUMAET(k,i,j) < NAEMIN .and. NUMAET(k,i,j) > EPS ) then
             do iq = 1, 6
                NUMAE(k,i,j,iq) = NAEMIN * NUMAE(k,i,j,iq) / NUMAET(k,i,j)
             enddo
          endif
          if ( NUMAET(k,i,j) < NAEMIN .and. NUMAET(k,i,j) > 0.0_RP ) then
             NUMAET(k,i,j) = NAEMIN
          endif
          
       enddo
    enddo
    enddo
    !================================================
    !end setup parameter 


    !CCN (Ghan et al., JGR, 1997)
    !================================================
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          CURV(k,i,j) = 2.0_RP * STENSW / DWATR / RVAP / GDT(k,i,j)
          call ATMOS_SATURATION_dens2qsat_all( GDT(k,i,j), DENS(k,i,j), SATQ(k,i,j) )
!!          SATQ(k,i,j) = EPSV * PSAT0 / GDP(k,i,j) * exp( ( EL+EMELT/2.0_RP*(1.0_RP-sign(1.0_RP,GDT(k,i,j)-TMELT)) ) &
!!                                                         / RVAP * ( 1.0_RP/TMELT - 1.0_RP/GDT(k,i,j) )              )
!!          SATQT       = ( EL + EMELT / 2.0_RP * ( 1.0_RP - sign(1.0_RP,GDT(k,i,j)-TMELT) ) ) * SATQ(k,i,j) / ( RVAP * (GDT(k,i,j)**2) )
          call ATMOS_SATURATION_dqs_dtem_dens_liq( GDT(k,i,j), DENS(k,i,j), & ! [IN]
                                                   dqsat_dT_liq, qsat_liq_  ) ! [OUT]
          call ATMOS_SATURATION_dqs_dtem_dens_ice( GDT(k,i,j), DENS(k,i,j), & ! [IN]
                                                   dqsat_dT_ice, qsat_ice_  ) ! [OUT]
          call ATMOS_SATURATION_alpha   ( GDT(k,i,j), alpha_    ) ! [IN], [OUT]
          call ATMOS_SATURATION_dalphadT( GDT(k,i,j), dalpha_dT ) ! [IN], [OUT]
      
          SATQT  = qsat_liq_ * dalpha_dT + dqsat_dT_liq * (        alpha_ ) &
                 - qsat_ice_ * dalpha_dT + dqsat_dT_ice * ( 1.0_RP-alpha_ )

          ALFA(k,i,j) = GRAV * ( SATQT / SATQ(k,i,j) / CP - DENS(k,i,j) / GDP(k,i,j) )
          if ( k > KS ) then
             DFF = 2.11E-5_RP * ( GDT(k,i,j)/TMELT )** 1.94_RP * ( P00/GDP(k,i,j) )
             BG(k,i,j) = 1.0_RP / ( DWATR * RVAP * GDT(k,i,j) / max( DFF, DFEMIN ) / ( SATQ(k,i,j) * GDP(k,i,j) / EPSV ) &
                                    + EL * DWATR / CONDE / GDT(k,i,j) * ( EL / RVAP / GDT(k,i,j) - 1.0_RP ) )
             AVRW(k,i,j) = GDW(k,i,j)
             !TURW(k,i,j) = SC * sqrt( min( max( TKE(k,i,j)*0.5_RP, TKEMIN ), TKEMAX ) )
             TURW(k,i,j) = SC * sqrt( min( max( TKE(k,i,j), TKEMIN ), TKEMAX ) )
          else
             BG(k,i,j) = 1.0_RP / ( EL * DWATR / CONDE / GDT(k,i,j) * ( EL / RVAP / GDT(k,i,j) - 1.0_RP ) )
             AVRW(k,i,j) = 0.0_RP
             TURW(k,i,j) = 0.0_RP
          endif
          BETA(k,i,j) = 4.0_RP * PI * DWATR * BG(k,i,j) / DENS(k,i,j) / SATQ(k,i,j) * ( 1.0_RP + EL / CP * SATQT )
          OMGCCN(k,i,j) = AVRW(k,i,j) + TURW(k,i,j)
          OMGCCN(k,i,j) = max( OMGCCN(k,i,j), OMGMIN )
       enddo
    enddo
    enddo

    FJ    (:,:,:) = 0.0_RP
    FJ01  (:,:,:) = 0.0_RP
    FJ03  (:,:,:) = 0.0_RP
    UAPCL (:,:,:) = 0.0_RP
    UAPCL1(:,:,:) = 0.0_RP
    UAPCL3(:,:,:) = 0.0_RP
    UAPCL1C (:,:) = 0.0_RP
    UAPCL3C (:,:) = 0.0_RP
    UAPCB   (:,:) = CONST_UNDEF !VMISS
    
    do iq = 1, 6
       if ( iq /= 2 ) then
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                SM(k,i,j,iq) = 2.0_RP / sqrt( HYGRAE(iq) ) * ( CURV(k,i,j) / 3.0_RP / RADAE(k,i,j,iq) )**1.5_RP
                FJ_term1 = F1(iq) * ( (CURV(k,i,j)*NUMAE(k,i,j,iq)*BETA(k,i,j))/(3.0_RP*ALFA(k,i,j)*OMGCCN(k,i,j)) )**2
                FJ_term2 = F2(iq) * ( 2.0_RP * CURV(k,i,j)**3 * NUMAE(k,i,j,iq) * BETA(k,i,j) * sqrt( BG(k,i,j) ) ) &
                                  / ( 27.0_RP * HYGRAE(iq) * RADAE(k,i,j,iq)**3 * ( ALFA(k,i,j) * OMGCCN(k,i,j) )**1.5_RP )
                FJ  (k,i,j) = FJ  (k,i,j) + ( FJ_term1 + FJ_term2 )**3 / SM(k,i,j,iq)**2
                FJ01(k,i,j) = FJ01(k,i,j) + ( FJ_term1 + FJ_term2 )**3 / 1.0E-3_RP**2
                FJ03(k,i,j) = FJ03(k,i,j) + ( FJ_term1 + FJ_term2 )**3 / 3.0E-3_RP**2
             enddo
          enddo
          enddo
       else
          SM(KS:KE,IS:IE,JS:JE,iq) = 0.0_RP
       endif
    enddo
    
    do iq = 1, 6
       if ( iq /= 2 ) then
          do j = JS, JE
          do i = IS, IE
             do k = KS, KE
                ACTI  (k,i,j,iq) = 1.0_RP / ( 1.0_RP + ( SM(k,i,j,iq)**2 * FJ  (k,i,j) )**(SB(iq)/3.0_RP) )!1/3???88888888888
                ACTI01(k,i,j,iq) = 1.0_RP / ( 1.0_RP + ( SM(k,i,j,iq)**2 * FJ01(k,i,j) )**(SB(iq)/3.0_RP) )!1/3???88888888888
                ACTI03(k,i,j,iq) = 1.0_RP / ( 1.0_RP + ( SM(k,i,j,iq)**2 * FJ03(k,i,j) )**(SB(iq)/3.0_RP) )!1/3???88888888888
             enddo
          enddo
          enddo
       else
          ACTI  (KS:KE,IS:IE,JS:JE,iq) = 0.0_RP
          ACTI01(KS:KE,IS:IE,JS:JE,iq) = 0.0_RP
          ACTI03(KS:KE,IS:IE,JS:JE,iq) = 0.0_RP
       endif
    enddo

    do iq = 1, 6
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             NUMCLD(k,i,j,iq) = NUMAE(k,i,j,iq) * ACTI  (k,i,j,iq)
             NCLD01(k,i,j,iq) = NUMAE(k,i,j,iq) * ACTI01(k,i,j,iq)
             NCLD03(k,i,j,iq) = NUMAE(k,i,j,iq) * ACTI03(k,i,j,iq)
             UAPCL (k,i,j)    = UAPCL (k,i,j) + NUMCLD(k,i,j,iq)
             UAPCL1(k,i,j)    = UAPCL1(k,i,j) + NCLD01(k,i,j,iq)
             UAPCL3(k,i,j)    = UAPCL3(k,i,j) + NCLD03(k,i,j,iq)
          enddo
       enddo
       enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          UAPCL1C(i,j) = UAPCL1C(i,j) + UAPCL1(k,i,j) * DELZ(k,i,j)
          UAPCL3C(i,j) = UAPCL3C(i,j) + UAPCL3(k,i,j) * DELZ(k,i,j)
       enddo
    enddo
    enddo
    
    do j = JS, JE
    do i = IS, IE
       do k = KE, KS+1, -1
          if ( TCLDF2(k,i,j) > EPSC ) then
             TCLDWC1 = QC(k,i,j)
          else
             TCLDWC1 = 0.0_RP
          endif
          if ( TCLDWC1 > CLWMIN ) then
             UAPCB(i,j) = UAPCL(k,i,j)
          endif
       enddo
    enddo
    enddo
    
    !================================================
    !end CCN

    
    !OUTPUT
    !================================================
    do iq = 1, 6
       !                            k i j
       call FILE_HISTORY_in( NUMAE (:,:,:,iq), 'NMAE_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'aerosol number concentration ('//trim(adjustl(cae(iq)))//')',                 '1/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( NUMCLD(:,:,:,iq), 'NCCN_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'activated aerosol (CCN) number concentration ('//trim(adjustl(cae(iq)))//')', '1/m3', dim_type = 'ZXY' )
       call FILE_HISTORY_in( ACTI  (:,:,:,iq), 'ACTI_'//trim(adjustl(cae(iq)))//'_SPRINTARS', 'ratio of activated aerosol (CCN) ('//trim(adjustl(cae(iq)))//')',             '1',    dim_type = 'ZXY' )
    enddo
    !                            k i j
    call FILE_HISTORY_in( NUMAET(:,:,:), 'NMAET_SPRINTARS',  'total aerosol number concentration',           '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( UAPCL (:,:,:), 'NCCN_SPRINTARS',   'activated aerosol (CCN) number concentration', '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( UAPCL1(:,:,:), 'CCN01_SPRINTARS',  'CCN(0.1%)',                                    '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( UAPCL3(:,:,:), 'CCN03_SPRINTARS',  'CCN(0.3%)',                                    '1/m3', dim_type = 'ZXY' )
    call FILE_HISTORY_in( UAPCL1C (:,:), 'CCN01C_SPRINTARS', 'CCN(0.1%) (in column)',                        '1/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( UAPCL3C (:,:), 'CCN03C_SPRINTARS', 'CCN(0.3%) (in column)',                        '1/m2', dim_type = 'XY'  )
    call FILE_HISTORY_in( UAPCB   (:,:), 'UAPCB_SPRINTARS',  'CCN number at cloud base',                     '1/m3', dim_type = 'XY'  )
    call FILE_HISTORY_in( OMGCCN(:,:,:), 'OMGCCN_SPRINTARS', 'updraft velocity',                             'm/s',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( AVRW  (:,:,:), 'AVRW_SPRINTARS',   'updraft velocity (average)',                   'm/s',  dim_type = 'ZXY' )
    call FILE_HISTORY_in( TURW  (:,:,:), 'TURW_SPRINTARS',   'updraft velocity (turbulence)',                'm/s',  dim_type = 'ZXY' )
    !================================================
    !end OUTPUT
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_CCN_AEROCCN
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_ccn

