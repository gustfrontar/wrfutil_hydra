!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          Cloud Microphysics by Super Droplet Method (SDM)
!!
!! - Reference
!!  - Shima et al., 2009:
!!    The super-droplet method for the numerical simulation of clouds and precipitation:
!!    A particle-based and probabilistic microphysics model coupled with a non-hydrostatic model.
!!    Quart. J. Roy. Meteorol. Soc., 135: 1307-1320
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_mp_sdm
  !-----------------------------------------------------------------------------
  !
  !++ used modules ! For encapsulation, reduce the use of modules here as far as possible.
  !   Modules should be called inside the subroutine here.
  !
  use mpi
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only:  &
     mype => PRC_myrank,    &
     PRC_nprocs,            &
     PRC_IsMaster,          &
     PRC_abort
  use scale_time, only:     &
     TIME_NOWSEC

  use scale_atmos_hydrometeor, only: &
    N_HYD, &
    I_HC,  &
    I_HR,  &
    I_HI,  &
    I_HS,  &
    I_HG,  &
    I_HH

#ifndef DISABLE_SDM
  use m_sdm_common, only: &
    QA_MP_sdm
#endif

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  public :: ATMOS_PHY_MP_sdm_tracer_setup
  public :: ATMOS_PHY_MP_sdm_setup
  public :: ATMOS_PHY_MP_sdm_adjustment
  public :: ATMOS_PHY_MP_sdm_Cloud_Fraction
  public :: ATMOS_PHY_MP_sdm_Effective_Radius
  public :: ATMOS_PHY_MP_sdm_qtrc2qhyd
  public :: ATMOS_PHY_MP_sdm_qtrc2nhyd
  public :: ATMOS_PHY_MP_sdm_restart_read
  public :: ATMOS_PHY_MP_sdm_restart_write
  public :: ATMOS_PHY_MP_sdm_random_setup
  public :: ATMOS_PHY_MP_sdm_additional_output

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------

  integer, public :: ATMOS_PHY_MP_sdm_ntracers
  integer, public :: ATMOS_PHY_MP_sdm_nwaters
  integer, public :: ATMOS_PHY_MP_sdm_nices
  integer, public :: ATMOS_PHY_MP_sdm_nccn
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_MP_sdm_tracer_names(:)
  character(len=H_MID)  , public, allocatable :: ATMOS_PHY_MP_sdm_tracer_descriptions(:)
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_MP_sdm_tracer_units(:)
  logical, public, allocatable :: ATMOS_PHY_MP_sdm_tracer_advc(:)
  logical,  public, save   :: sd_rest_flg_out = .false. ! restart flg of Super D

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer :: histitemid(5)
  logical :: initialized

#ifdef DISABLE_SDM
  integer, parameter, private :: QA_MP_sdm = 1
#endif
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Configure
  subroutine ATMOS_PHY_MP_sdm_tracer_setup
#ifdef DISABLE_SDM
    use scale_prc, only: &
         PRC_abort
#else
    use m_sdm_common, only: &
         I_QV_sdm
#endif
    implicit none

    !---------------------------------------------------------------------------

#ifdef DISABLE_SDM
    LOG_ERROR("ATMOS_PHY_MP_sdm_tracer_setup",*) "Scalelib was compiled without SDM"
    LOG_ERROR_CONT(*) "Recompile scalelib with SCALE_DISABLE_SDM=F"
    call PRC_abort
#else

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_sdm_tracer_setup",*) 'Tracer setup for Super Droplet Method'

    ! regist only water vapor
    QA_MP_sdm = 1
    ATMOS_PHY_MP_sdm_ntracers = QA_MP_sdm
    ATMOS_PHY_MP_sdm_nwaters = 0
    ATMOS_PHY_MP_sdm_nices = 0

    allocate( ATMOS_PHY_MP_sdm_tracer_names       (QA_MP_sdm) )
    allocate( ATMOS_PHY_MP_sdm_tracer_descriptions(QA_MP_sdm) )
    allocate( ATMOS_PHY_MP_sdm_tracer_units       (QA_MP_sdm) )
    allocate( ATMOS_PHY_MP_sdm_tracer_advc        (QA_MP_sdm) )

    ATMOS_PHY_MP_sdm_tracer_names(1:1) = (/ &
         'QV' /)

    ATMOS_PHY_MP_sdm_tracer_descriptions(1:1) = (/ &
         'Ratio of Water Vapor mass to total mass (Specific humidity)' /)

    ATMOS_PHY_MP_sdm_tracer_units(1:1) = (/ &
         'kg/kg' /)

    ATMOS_PHY_MP_sdm_tracer_advc(1:1) = (/ &
         .true. /)

    I_QV_sdm = 1

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_tracer_setup
  !-----------------------------------------------------------------------------
  !> Setup Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_setup( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       REAL_FZ, FX, FY, &
       CDX, CDY, RCDX, RCDY, DX, DY, &
       dt, &
       do_precipitation )
#ifndef DISABLE_SDM
    use scale_prc_cartesC, only: &
       PRC_next,    &
       PRC_NUM_X,   &
       PRC_NUM_Y,   &
       PRC_W,       &
       PRC_E,       &
       PRC_S,       &
       PRC_N
    use scale_const, only: &
       dens_w => CONST_DWATR, &
       dens_i => CONST_DICE
!    use scale_specfunc, only: &
!       SF_gamma
!    use scale_comm_cartesC, only: &
!       COMM_horizontal_mean
    use scale_prc, only: &
         PRC_abort
    use scale_file_history, only: &
         FILE_HISTORY_reg
    use m_sdm_coordtrans, only: &
       sdm_getrklu
    use m_sdm_ptlformation, only: &
         sdm_numset
    use m_sdm_memmgr, only: &
         sdm_allocinit
    use m_sdm_common
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: REAL_FZ(0:KA,IA,JA)
    real(RP), intent(in) :: FX(0:IA)
    real(RP), intent(in) :: FY(0:JA)
    real(RP), intent(in) :: CDX(IA)
    real(RP), intent(in) :: CDY(JA)
    real(RP), intent(in) :: RCDX(IA)
    real(RP), intent(in) :: RCDY(JA)
    real(RP), intent(in) :: DX
    real(RP), intent(in) :: DY
    real(DP), intent(in) :: dt
    logical,  intent(in) :: do_precipitation

#ifndef DISABLE_SDM
!    real(RP) :: dtevl
!    real(RP) :: n0, dry_r
!    real(RP) :: delta1, delta2 !, sdn_tmp
    real(RP) :: buffact
    integer :: ierr
    integer :: i, j
    integer :: ilcm, igcd
    character(len=17) :: fmt1="(A, '.', A, I*.*)"
    character(:), allocatable :: tmp_value
    integer len, status

! namelist variables
  NAMELIST / PARAM_ATMOS_PHY_MP_SDM / &
       RANDOM_IN_BASENAME,  &
       RANDOM_OUT_BASENAME, &
       SD_IN_BASENAME,      &
       SD_OUT_BASENAME,     &
       docondensation,      &
       docollision,    &
       domovement,          &
       sdm_cold,            &
       domeltfreeze,   &
       dosublimation,   &
       donegative_fixer,    &
       sdm_dtcmph,          &
       sdm_rdnc,            &
       sdm_sdnmlvol,        &
       sdm_inisdnc,         &
       sdm_aslset,          &
       sdm_fctr2multi,      &
       sdm_aslmw,           &
       sdm_zlower,          &
       sdm_zupper,          &
       sdm_extbuf,          &
       sdm_rqc2qr,          &
       sdm_rqi2qs,          &
       sdm_aslfmrate,       &
       sdm_aslfmdt,         &
       sdm_aslfmsdnc,       &
       sdm_colbrwn,         &
       sdm_colkrnl,         &
       sdm_mvexchg,         &
       sdm_nadjdt,          &
       sdm_nadjvar,         &
       sdm_dmpvar,          &
       sdm_dmpitva,         &
       sdm_dmpnskip,        &
       sdm_dmpitvb,         &
       sdm_dmpitvl,         &
       sdm_dmpsdsiz,        &
       sdm_wbc,sdm_ebc,sdm_sbc,sdm_nbc

    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_MP_sdm_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_ATMOS_PHY_MP_SDM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sdm_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sdm_setup",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_MP_SDM. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_MP_SDM)

    dosedimentation = do_precipitation

    ! [Bug] Horizontally stretched coordinate is not supported yet, but the check below is not working any more
    ! [Bug] Must be fixed in the future when introducing generalized coordinate
    buffact = 0.0_RP
    do j = JS, JE
      buffact = max( buffact,CDY(j)/DY )
    enddo
    do i = IS, IE
      buffact = max( buffact,CDX(i)/DX )
    enddo

    if( buffact < 1.0_RP .or. buffact > 1.0_RP ) then
       sthopt = 1
    else
       sthopt = 0
    endif

    if(sthopt==1) then
       write(*,*) 'ERROR: horizontally stretched coordinate is not yet supported!', buffact
       call PRC_abort
    end if
    ! [Bug] Horizontally stretched coordinate is not supported yet, but the check above is not working any more

    nsub  = max( PRC_nprocs,1 )

    !--- set process next to the process
    dstw_sub = PRC_next(PRC_W)
    dste_sub = PRC_next(PRC_E)
    srcw_sub = PRC_next(PRC_W)
    srce_sub = PRC_next(PRC_E)

    dsts_sub = PRC_next(PRC_S)
    dstn_sub = PRC_next(PRC_N)
    srcs_sub = PRC_next(PRC_S)
    srcn_sub = PRC_next(PRC_N)

    tag  = 0

    if( .not. sdm_cold ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_setup",*) 'Warm SDM is used'
    else if ( sdm_cold ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_setup",*) 'Cold SDM is used'
    endif

    if( (.not. sdm_cold ) .and. domeltfreeze ) then
       LOG_WARN("ATMOS_PHY_MP_sdm_setup",*) 'domeltfreeze option does not work for warm SDM'
    endif

    if( RANDOM_IN_BASENAME == '' ) then
       LOG_ERROR("ATMOS_PHY_MP_sdm_setup",*) 'Random number generator file not specified! stop'
       LOG_ERROR_CONT(*) 'To create random number generator file, run scale_init first'
       call PRC_abort
    else
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(RANDOM_IN_BASENAME,fmt1) trim(RANDOM_IN_BASENAME),'pe',mype
    endif
    fid_random_i = IO_get_available_fid()

    if( SD_IN_BASENAME == '' ) then
       sd_rest_flg_in = .false.
       LOG_INFO("ATMOS_PHY_MP_sdm_setup",*) 'Not found S.D. file, generate from random number.'
    else
       sd_rest_flg_in = .true.
       fid_sd_i = IO_get_available_fid()
       write(fmt1(14:14),'(I1)') 6
       write(fmt1(16:16),'(I1)') 6
       write(SD_IN_BASENAME,fmt1) trim(SD_IN_BASENAME),'pe',mype
    endif

    if( SD_OUT_BASENAME == '' ) then
       sd_rest_flg_out = .false.
    else
       sd_rest_flg_out = .true.
    endif

    if( docondensation )   sdm_calvar(1) = .true.
    if( docollision ) sdm_calvar(2) = .true.
    if( domovement )       sdm_calvar(3) = .true.
    if( domeltfreeze )     sdm_calvar(4) = .true.
    if( dosublimation )    sdm_calvar(5) = .true.

! zlower+surface height is the lower boundary of SDs.
!!$     if( sdm_zlower < CZ(KS) ) then
!!$      if( mype == PRC_master )  write(*,*) "sdm_zlower was set to CZ(KS) because zlower < CZ(KS)"
!!$      sdm_zlower = CZ(KS)
!!$     endif
!

!     if( sdm_zupper > CZ(KE) ) then
!      if( mype == PRC_master )  write(*,*) "sdm_zupper was set to CZ(KE) because zupper > CZ(KE)"
!      sdm_zupper = CZ(KE)
!     endif

    if( sdm_zupper > minval(REAL_FZ(KE,IS:IE,JS:JE)) ) then
!       if( mype == PRC_master )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE)) because zupper > minval(REAL_FZ(KE))"
       if( PRC_IsMaster )  write(*,*) "sdm_zupper was set to minval(REAL_FZ(KE)) because zupper > minval(REAL_FZ(KE))"
       sdm_zupper = minval(REAL_FZ(KE,IS:IE,JS:JE))
    endif

    ! check whether sdm_dtcmph(1:3) > 0
     if(  ( (sdm_dtcmph(1) <= 0.0_RP) .and. docondensation   )   .or. &
         ( (sdm_dtcmph(2) <= 0.0_RP) .and. docollision )   .or. &
         ( (sdm_dtcmph(3) <= 0.0_RP) .and. domovement       )   .or. &
         ( (sdm_dtcmph(4) <= 0.0_RP) .and. domeltfreeze     )   .or. &
         ( (sdm_dtcmph(5) <= 0.0_RP) .and. dosublimation    )        ) then
       write(*,*) 'ERROR: sdm_dtcmph(1:5) have to be positive'
       call PRC_abort
     end if

    ! check whether dt (dt of mp) is larger than sdm_dtcmph(i).
    if( (dt < minval(sdm_dtcmph(:))) ) then
       write(*,*) 'ERROR: For now, sdm_dtcmph should be smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_abort
    end if

    ! aerosol nucleation and sd number adjustment functions are not supported yet
    if ( (abs(sdm_aslset) >= 10) .or. (sdm_nadjvar /= 0)) then
       write(*,*) 'ERROR: aerosol nucleation and sd number adjustment functions are not supported yet'
       write(*,*) 'ERROR: set sdm_aslset < 10 and sdm_nadjvar =0'
       call PRC_abort
    end if

    ! rigorous momentum exchange function is not supported yet
    if ( sdm_mvexchg >= 2) then
       write(*,*) 'ERROR: Rigorous momentum exchange not yet supported. set sdm_mvexchg = 0 (none) or 1 (rhow only)'
       call PRC_abort
    end if

    if( docondensation ) then
       nclstp(1)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(1)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(1)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_abort
       end if
    else
       nclstp(1) = 1
    end if

    if( docollision ) then
       nclstp(2)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(2)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(2)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_abort
       end if
    else
       nclstp(2) = 1
    end if

    if( domovement ) then
       nclstp(3)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(3)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(3)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_abort
       end if
    else
       nclstp(3) = 1
    end if

    if( domeltfreeze ) then
       nclstp(4)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(4)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(4)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_abort
       end if
    else
       nclstp(4) = 1
    end if

    if( dosublimation ) then
       nclstp(5)=10*int(1.E+5_RP*(dt+1.E-6_RP))            &
            /int(1.E+6_RP*(sdm_dtcmph(5)+1.E-7_RP))
       if(mod(10*int(1.E+5_RP*(dt+1.E-6_RP)),int(1.E+6_RP*(sdm_dtcmph(5)+1.E-7_RP))) /= 0) then
          write(*,*) 'ERROR: sdm_dtcmph should be submultiple of TIME_DTSEC_ATMOS_PHY_MP'
          call PRC_abort
       end if
    else
       nclstp(5) = 1
    end if

    nclstp(0)=min(nclstp(1),nclstp(2),nclstp(3),nclstp(4),nclstp(5))

    ilcm=0
    iterate_2: do
       ilcm=ilcm+1
       nclstp(0)=nclstp(0)*ilcm
       if( mod(nclstp(0),nclstp(1)) == 0 .and.           &
           mod(nclstp(0),nclstp(2)) == 0 .and.           &
           mod(nclstp(0),nclstp(3)) == 0 .and.           &
           mod(nclstp(0),nclstp(4)) == 0 .and.           &
           mod(nclstp(0),nclstp(5)) == 0   ) then
           exit iterate_2
       else
          nclstp(0)=nclstp(0)/ilcm
       end if
    end do iterate_2

    ! check whether dt is the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt
    !! find the smallest nclstp(1:3) that satisfies sdm_dtcmph(i)<= dt
    igcd=maxval(nclstp(1:5))
    if((sdm_dtcmph(1)<= dt).and.docondensation)   igcd=min(nclstp(1),igcd)
    if((sdm_dtcmph(2)<= dt).and.docollision) igcd=min(nclstp(2),igcd)
    if((sdm_dtcmph(3)<= dt).and.domovement)       igcd=min(nclstp(3),igcd)
    if((sdm_dtcmph(4)<= dt).and.domeltfreeze)     igcd=min(nclstp(4),igcd)
    if((sdm_dtcmph(5)<= dt).and.dosublimation)    igcd=min(nclstp(5),igcd)
    !! find the greatest common divisor of nclstp(1:5) that satisfies sdm_dtcmph(i)<= dt
    do
       if( igcd == 1) exit
       if( mod(nclstp(1),igcd) == 0 .and.           &
           mod(nclstp(2),igcd) == 0 .and.           &
           mod(nclstp(3),igcd) == 0 .and.           &
           mod(nclstp(4),igcd) == 0 .and.           &
           mod(nclstp(5),igcd) == 0   ) then
           exit
       end if
       igcd = igcd - 1
    end do
    !! if igcd>1, it meanst dt is not the least common multiple of sdm_dtcmph(i) that satisfies sdm_dtcmph(i)<= dt
    if( igcd > 1)then
       write(*,*) 'ERROR: TIME_DTSEC_ATMOS_PHY_MP should be the least comon multiple of sdm_dtcmph(1:5)', &
                  ' that are smaller than TIME_DTSEC_ATMOS_PHY_MP'
       call PRC_abort
    end if

    xmax_sdm = FX(IE)-FX(IS-1)
    ymax_sdm = FY(JE)-FY(JS-1)

    allocate( zph_crs(KA,IA,JA) )

    !--- set number of super droplet etc...
    call sdm_numset(              &
      KA, KS, KE, IA, IS, IE, JA, JS, JE, &
      sdm_extbuf,                 &
      sdm_aslset, sdm_sdnmlvol,   &
      sdm_inisdnc,                &
      sdm_aslfmsdnc,              &
      sdm_zlower, sdm_zupper,     &
      minzph, sdininum_s2c,       &
      sdfmnum_s2c, sdnum_s2c,     &
      ni_s2c, nj_s2c, nk_s2c, zph_crs )

    call sdm_allocinit( KA, IA, JA )

    !### Get index[k/real] at "sdm_zlower" and "sdm_zupper"  ###!
    call sdm_getrklu( &
         KA, KS, KE, IA, IS, IE, JA, JS, JE, &
         sdm_zlower,sdm_zupper,      &
         sdrkl_s2c,sdrku_s2c)

    dx_sdm(1:IA) = CDX(1:IA)
    dy_sdm(1:JA) = CDY(1:JA)
    dxiv_sdm(1:IA) = RCDX(1:IA)
    dyiv_sdm(1:JA) = RCDY(1:JA)


    ! history
    call FILE_HISTORY_reg( 'SNC', 'SD number density', 'num/m3', histitemid(1), ndims=3)
    call FILE_HISTORY_reg( 'NC', 'droplet number density', 'num/m3', histitemid(2), ndims=3)
    call FILE_HISTORY_reg( 'DMOM1', 'density of the 1st droplet moment', 'm/m3', histitemid(3), ndims=3)
    call FILE_HISTORY_reg( 'DMOM2', 'density of the 2nd droplet moment', 'm2/m3', histitemid(4), ndims=3)
    call FILE_HISTORY_reg( 'DMOM3', 'density of the 3rd droplet moment', 'm3/m3', histitemid(5), ndims=3)

    ! reset the array to record surface precipitation
    prr_crs(1:IA,1:JA,1:6)=0.0_RP

    initialized = .false.

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_setup
  !-----------------------------------------------------------------------------
  !> Cloud Microphysics
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_adjustment(  &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS,      &
       MOMZ,      &
       MOMX,      &
       MOMY,      &
       RHOT,      &
       QTRC,      &
       CCN,       &
       EVAPORATE, &
       SFLX_rain, &
       SFLX_snow, &
       dt         )
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_in,    &
       FILE_HISTORY_put
    use scale_const, only: &
       rw => CONST_DWATR, &
       Rvap  => CONST_Rvap
#ifndef DISABLE_SDM
    use m_sdm_calc, only: &
       sdm_calc
    use m_sdm_coordtrans, only: &
       sdm_rk2z
    use m_sdm_motion, only: &
       sdm_getvz_liq, sdm_getvz_ice
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr,sdm_sd2qiqsqg,sdm_sd2rhosd,sdm_sd2rhodropmom
    use m_sdm_idutil, only: &
         sdm_copy_selected_sd
    use m_sdm_ptlformation, only: &
         sdm_iniset
    use m_sdm_common
#endif
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(inout) :: MOMZ(KA,IA,JA)        !! Momentum [kg/s/m2]
    real(RP), intent(inout) :: MOMX(KA,IA,JA)
    real(RP), intent(inout) :: MOMY(KA,IA,JA)
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]
    real(RP), intent(out)   :: EVAPORATE(KA,IA,JA)   !---- evaporated cloud number concentration [/m3]
    real(RP), intent(in)    :: CCN(KA,IA,JA)         !! cloud condensation nuclei calculated by aerosol module
    real(RP), intent(out)   :: SFLX_rain(IA,JA)      !! rain flux at surface (used in land and urban model)
    real(RP), intent(out)   :: SFLX_snow(IA,JA)      !! snow flux at surface (used in land and urban model)
    real(DP), intent(in)    :: dt

#ifndef DISABLE_SDM
    ! Work variables
    logical :: lsdmup       ! flag for updating water hydrometeor by SDM

    real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA)
    real(RP) :: crs_dtmp3(KA,IA,JA), crs_dtmp4(KA,IA,JA)
    integer  :: i, j         ! index

    !-------------------------------------------------------------------------------------------------------------------------------

#ifdef _FIPP_
    ! Section specification for fipp profiler
    call fipp_start()
#endif
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_all",0,0)
#endif
    LOG_PROGRESS(*) 'atmosphere / physics / microphysics / SDM'

    ! -----
    if( .not. sdm_calvar(1) .and. &
        .not. sdm_calvar(2) .and. &
        .not. sdm_calvar(3) .and. &
        .not. sdm_calvar(4) .and. &
        .not. sdm_calvar(5)       ) return

    !== run SDM at future ==!
     call sdm_calc(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                   MOMX,MOMY,MOMZ,DENS,RHOT,QTRC,                 &
                   sdm_calvar,sdm_rdnc,sdm_aslset,             &
                   sdm_inisdnc,sdm_sdnmlvol,                   &
                   sdm_mvexchg,   &
                   prr_crs,zph_crs,                      &
                   lsdmup,ni_s2c,nj_s2c,nk_s2c,                   &
                   sdnum_s2c,sdnumasl_s2c,                        &
                   sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,    &
                   sdu_s2c,sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,sdice_s2c,    &
                   sdfmnum_s2c,        sdn_fm,sdliqice_fm,sdx_fm,sdy_fm,sdz_fm,    &
                   sdri_fm,sdrj_fm,sdrk_fm,sdr_fm,sdasl_fm,sdice_fm, &
                   sdrkl_s2c,sdrku_s2c,                           &
                   rng_s2c,rand_s2c,sortid_s2c,sortkey_s2c,       &
                   sortfreq_s2c,sorttag_s2c,                      &
                   bufsiz1,                                       &
                   bufsiz2_r8,bufsiz2_i8,bufsiz2_i2,bufsiz2_i4,   &
                   sdm_itmp1,sdm_itmp2,           &
                   sd_itmp1,sd_itmp2,sd_itmp3,sd_dtmp1,sd_dtmp2,sd_dtmp3,sd_dtmp4,     &
                   crs_dtmp1,crs_dtmp2,crs_dtmp3,crs_dtmp4,       &
                   rbuf_r8,sbuf_r8,rbuf_i8,sbuf_i8,rbuf_i2,sbuf_i2,rbuf_i4,sbuf_i4, &
                   dt )

     !== convert updated contravariant velocity of ==!
     !== super-droplets to {u,v,w} at future       ==!

!     if( sdm_mvexchg>0 ) then

!        call cnt2phy(idsthopt,idtrnopt,idmpopt,idmfcopt,          &
!                     ni,nj,nk,j31,j32,jcb8w,mf,                   &
!                     uf_crs,vf_crs,wcf_crs,wf_crs,                &
!                     rtmp4,rtmp5,rtmp6)
!     end if

! not supported yet. S.Shima
!!$     ! Aerosol formation process of super-droplets
!!$
!!$     if( sdm_aslfmdt>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(sdm_aslfmdt+0.00010_RP))) == 0 ) then
!!$
!!$        call sdm_aslform(DENS,RHOT,QTRC,                             &
!!$                         sdm_calvar,sdm_aslset,                      &
!!$                         sdm_aslfmsdnc,sdm_sdnmlvol,                 &
!!$                         sdm_zupper,sdm_zlower,dtcl,                 &
!!$                         jcb,pbr_crs,ptbr_crs,ppf_crs,               &
!!$                         ptpf_crs,qvf_crs,zph_crs,rhod_crs,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,     &
!!$                         sdy_s2c,sdz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdu_s2c,           &
!!$                         sdv_s2c,sdvz_s2c,sdr_s2c,sdasl_s2c,         &
!!$                         sdfmnum_s2c,sdn_fm,sdx_fm,sdy_fm,sdz_fm,    &
!!$                         sdri_fm,sdrj_fm,sdrk_fm,sdvz_fm,sdr_fm,sdasl_fm,            &
!!$                         ni_s2c,nj_s2c,nk_s2c,                       &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,        &
!!$                         sorttag_s2c,rng_s2c,                        &
!!$!                         sorttag_s2c,                                &
!!$                         sdm_itmp1,sd_itmp1,sd_itmp2,sd_itmp3)
!!$
!!$     end if

! not supported yet. S.Shima
!!$     ! Adjust number of super-droplets
!!$
!!$     if( sdm_nadjdt>0.0_RP .and. &
!!$         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
!!$             int(1.E+3_RP*(sdm_nadjdt+0.00010_RP))) == 0 ) then
!!$
!!$       !== averaged number concentration in a grid ==!
!!$
!!$       sd_nc = sdininum_s2c/real(ni_s2c*nj_s2c*knum_sdm,kind=RP)
!!$
!!$       call sdm_adjsdnum(sdm_nadjvar,ni_s2c,nj_s2c,nk_s2c,          &
!!$                         sdnum_s2c,sdnumasl_s2c,sd_nc,              &
!!$                         sdn_s2c,sdx_s2c,sdy_s2c,sdr_s2c,           &
!!$                         sdasl_s2c,sdvz_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,               &
!!$                         sortid_s2c,sortkey_s2c,sortfreq_s2c,       &
!!$                         sorttag_s2c,rng_s2c,rand_s2c,                &
!!$!                         sorttag_s2c,rand_s2c,                      &
!!$                         sdm_itmp1,sdm_itmp2,sd_itmp1)
!!$
!!$      end if

    do i = IS, IE
    do j = JS, JE
       SFLX_rain(i,j) = prr_crs(i,j,1) * rw ! surface rain rate [kg/m2/s]
       SFLX_snow(i,j) = prr_crs(i,j,3) * rw ! surface snow rate [kg/m2/s]
       !! check the unit needed for SFLX
    enddo
    enddo

    !$omp workshare
    EVAPORATE(:,:,:) = 0.0_RP
    !$omp end workshare

#ifdef _FIPP_
    ! Section specification for fipp profiler
    call fipp_stop()
#endif
#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_all",0,0)
#endif

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_adjustment
  !-----------------------------------------------------------------------------
  !> Calculate Effective Radius
  subroutine ATMOS_PHY_MP_sdm_Effective_Radius( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Re     )
#ifndef DISABLE_SDM
    use m_sdm_coordtrans, only: &
       sdm_x2ri, sdm_y2rj
    use m_sdm_common
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(out) :: Re   (KA,IA,JA,N_HYD) ! effective radius

#ifndef DISABLE_SDM
    real(RP) :: sum3(KA,IA,JA,I_HC:I_HG), sum2(KA,IA,JA,I_HC:I_HG)
    real(RP) :: drate             ! temporary
    integer  :: k, i, j
    integer  :: tlist_l            ! total list number for cloud
    integer  :: lcnt               ! counter
    integer  :: kl                 ! index
    integer  :: ku                 ! index
    integer  :: m, n               ! index
    integer  :: ilist_l(1:sdnum_s2c)
    real(RP), parameter :: m2cm = 100.0_RP
    real(RP) :: r2,r3
    integer  :: ic
    !---------------------------------------------------------------------------

    call sdm_x2ri(IS,IE,sdnum_s2c,sdx_s2c,sdri_s2c,sdrk_s2c)
    call sdm_y2rj(JS,JE,sdnum_s2c,sdy_s2c,sdrj_s2c,sdrk_s2c)

    ! Effective radius is defined by r3m/r2m=1.5/lambda.
    ! dcoef is cancelled by the division of sum3/sum2
!    do j  = JS, JE
!    do i  = IS, IE
!    do k  = KS, KE
!       dcoef(k,i,j) = dxiv_sdm(i) * dyiv_sdm(j) / (zph_crs(k,i,j)-zph_crs(k-1,i,j))
!    enddo
!    enddo
!    enddo

    tlist_l = 0

    lcnt = 0

    do n=1,sdnum_s2c

       if( sdrk_s2c(n)<VALID2INVALID ) cycle

       lcnt = lcnt + 1
       ilist_l(lcnt) = n

    end do

    tlist_l = lcnt

    if( tlist_l>0 ) then

       !-- reset
       !$omp workshare
       sum3(:,:,:,:) = 0.0_RP
       sum2(:,:,:,:) = 0.0_RP
       !$omp end workshare

       do m=1,tlist_l
          n = ilist_l(m)

          i = floor(sdri_s2c(n))+1
          j = floor(sdrj_s2c(n))+1
          k = floor(sdrk_s2c(n))+1

          if ( sdliqice_s2c(n) == STAT_LIQ ) then
             if ( sdr_s2c(n) < sdm_rqc2qr ) then
                ic = I_HC
             else
                ic = I_HR
             end if
             r2 = sdr_s2c(n)**2
             r3 = sdr_s2c(n) * r2
          else if ( sdliqice_s2c(n) == STAT_ICE ) then
             r3 = sdice_s2c%re(n)**2 * sdice_s2c%rp(n)
             r2 = r3**(2.0_RP/3.0_RP)
             if ( sdice_s2c%rho(n) * F_THRD * ONE_PI * r3 * 0.3_RP < sdice_s2c%mrime(n) ) then
                ic = I_HG
             else if ( sdice_s2c%nmono(n) > 10 ) then
                ic = I_HS
             else
                ic = I_HI
             end if
          end if

          r3 = r3 * real(sdn_s2c(n),kind=RP)
          r2 = r2 * real(sdn_s2c(n),kind=RP)

          sum3(k,i,j,ic) = sum3(k,i,j,ic) + r3
          sum2(k,i,j,ic) = sum2(k,i,j,ic) + r2
       end do

       !=== correction at the verical boundary ===!

       !$omp parallel do &
       !$omp private(kl,drate,ku)
       do j=JS, JE
       do i=IS, IE

          !! at lower boundary

          kl    = floor(sdrkl_s2c(i,j))+1
          drate = real(kl,kind=RP) - sdrkl_s2c(i,j)
          if( drate<0.50_RP ) then
             sum3(kl,i,j,:) = 0.0_RP           !! <50% in share
             sum2(kl,i,j,:) = 0.0_RP           !! <50% in share
          else
             sum3(kl,i,j,:) = sum3(kl,i,j,:)/drate
             sum2(kl,i,j,:) = sum2(kl,i,j,:)/drate
          end if

          !! at upper boundary

          ku    = floor(sdrku_s2c(i,j))+1
          drate = sdrku_s2c(i,j) - real(ku-1,kind=RP)

          if( drate<0.50_RP ) then
             sum3(ku,i,j,:) = 0.0_RP           !! <50% in share
             sum2(ku,i,j,:) = 0.0_RP           !! <50% in share
          else
             sum3(ku,i,j,:) = sum3(ku,i,j,:)/drate
             sum2(ku,i,j,:) = sum2(ku,i,j,:)/drate
          end if

       end do
       end do

       !$omp parallel do collapse(3)
       do ic = I_HC, I_HG
       do i = IS, IE
       do j = JS, JE
       do k = KS, KE
          if( sum2(k,i,j,ic) == 0.0_RP ) then
           Re(k,i,j,ic) = 0.0_RP
          else
           Re(k,i,j,ic) = sum3(k,i,j,ic)/sum2(k,i,j,ic)*m2cm
          endif
       end do
       end do
       end do
       end do

    else
       !$omp workshare
       Re(:,:,:,:) = 0.0_RP
       !$omp end workshare
    end if

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_Effective_Radius
  !-----------------------------------------------------------------------------
  !> Calculate Cloud Fraction
  subroutine ATMOS_PHY_MP_sdm_Cloud_Fraction( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       mask_criterion, DENS, &
       cldfrac         )
#ifndef DISABLE_SDM
    use scale_precision
    use m_sdm_common, only: &
         sdm_cold, &
         QHYD_sdm, &
         I_QHC_sdm, &
         I_QHR_sdm, &
         I_QHI_sdm, &
         I_QHS_sdm, &
         I_QHG_sdm, &
         zph_crs,                   &
         sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,        &
         sdr_s2c,sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
         sd_itmp1,sd_itmp2,sd_itmp3
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr,sdm_sd2qiqsqg
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE
    real(RP), intent(in)  :: mask_criterion
    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(out) :: cldfrac(KA,IA,JA)

#ifndef DISABLE_SDM
    real(RP) :: qhydro
    integer  :: k, i, j, iq
    real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
    real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
    real(RP) :: rhoi_scale(KA,IA,JA) ! cloud ice density
    real(RP) :: rhos_scale(KA,IA,JA) ! snow density
    real(RP) :: rhog_scale(KA,IA,JA) ! graupel density
    real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA), crs_dtmp3(KA,IA,JA)
    !---------------------------------------------------------------------------

    !### Diagnose QC and QR from super-droplets ###!
    !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)
    call sdm_sd2qcqr(  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                       DENS,QHYD_sdm(:,:,:,I_QHC_sdm),QHYD_sdm(:,:,:,I_QHR_sdm),          &
                       zph_crs,                   &
                       sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,        &
                       sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
                       rhoc_scale,rhor_scale,                      &
                       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2            )

    !### Diagnose QI, QS, QG from super-droplets ###!
    if( sdm_cold ) then
       call sdm_sd2qiqsqg( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      DENS,QHYD_sdm(:,:,:,I_QHI_sdm),QHYD_sdm(:,:,:,I_QHS_sdm),QHYD_sdm(:,:,:,I_QHG_sdm),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdliqice_s2c, sdx_s2c,sdy_s2c,                     &
                      sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoi_scale,rhos_scale,rhog_scale,                                 &
                      sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3)
    end if

    !$omp parallel do &
    !$omp private(qhydro)
    do k  = KS, KE
    do i  = IS, IE
    do j  = JS, JE
       qhydro = 0.0_RP
       do iq = I_QHC_sdm, I_QHR_sdm
          qhydro = qhydro + QHYD_sdm(k,i,j,iq)
       end do
       if ( sdm_cold ) then
          do iq = I_QHI_sdm, I_QHG_sdm
             qhydro = qhydro + QHYD_sdm(k,i,j,iq)
          end do
       end if
       cldfrac(k,i,j) = 0.5_RP + sign(0.5_RP, qhydro-mask_criterion)
    enddo
    enddo
    enddo

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_Cloud_Fraction
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_qtrc2qhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, &
       Qe     )
#ifndef DISABLE_SDM
    use m_sdm_common, only: &
         sdm_cold, &
         QHYD_sdm, &
         I_QHC_sdm, &
         I_QHR_sdm, &
         I_QHI_sdm, &
         I_QHS_sdm, &
         I_QHG_sdm, &
         zph_crs,                   &
         sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,        &
         sdr_s2c,sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
         sd_itmp1,sd_itmp2,sd_itmp3
    use m_sdm_sd2fluid, only: &
         sdm_sd2qcqr,sdm_sd2qiqsqg
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(out) :: Qe   (KA,IA,JA,N_HYD) ! mixing ratio of each cateory [kg/kg]

#ifndef DISABLE_SDM
    real(RP) :: rhoc_scale(KA,IA,JA) ! cloud water density
    real(RP) :: rhor_scale(KA,IA,JA) ! rain water density
    real(RP) :: rhoi_scale(KA,IA,JA) ! cloud ice density
    real(RP) :: rhos_scale(KA,IA,JA) ! snow density
    real(RP) :: rhog_scale(KA,IA,JA) ! graupel density
    real(RP) :: crs_dtmp1(KA,IA,JA), crs_dtmp2(KA,IA,JA), crs_dtmp3(KA,IA,JA)

    !---------------------------------------------------------------------------

    if ( ( .not. allocated(QHYD_sdm) ) .or. ( .not. initialized ) ) then
       ! for scale_init
       !$omp workshare
       Qe(:,:,:,:) = 0.0_RP
       !$omp end workshare
       return
    end if

    !### Diagnose QC and QR from super-droplets ###!
    !! note that when SDM is used QC := rhoc/(rhod+rhov), QR := rhor/(rhod+rhov)
    call sdm_sd2qcqr(  KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                       DENS,QHYD_sdm(:,:,:,I_QHC_sdm),QHYD_sdm(:,:,:,I_QHR_sdm),          &
                       zph_crs,                   &
                       sdnum_s2c,sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,        &
                       sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,            &
                       rhoc_scale,rhor_scale,                      &
                       sd_itmp1,sd_itmp2,crs_dtmp1,crs_dtmp2            )

    !### Diagnose QI, QS, QG from super-droplets ###!
    if( sdm_cold ) then
       call sdm_sd2qiqsqg( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      DENS,QHYD_sdm(:,:,:,I_QHI_sdm),QHYD_sdm(:,:,:,I_QHS_sdm),QHYD_sdm(:,:,:,I_QHG_sdm),        &
                      zph_crs,                                               &
                      sdnum_s2c,sdn_s2c,sdliqice_s2c, sdx_s2c,sdy_s2c,                     &
                      sdice_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,&
                      rhoi_scale,rhos_scale,rhog_scale,                                 &
                      sd_itmp1,sd_itmp2,sd_itmp3,crs_dtmp1,crs_dtmp2,crs_dtmp3)
    end if

    !$omp workshare
    Qe(:,:,:,I_HC) = QHYD_sdm(:,:,:,I_QHC_sdm)
    Qe(:,:,:,I_HR) = QHYD_sdm(:,:,:,I_QHR_sdm)
    !$omp end workshare
    if( sdm_cold ) then
       !$omp workshare
       Qe(:,:,:,I_HI) = QHYD_sdm(:,:,:,I_QHI_sdm)
       Qe(:,:,:,I_HS) = QHYD_sdm(:,:,:,I_QHS_sdm)
       Qe(:,:,:,I_HG) = QHYD_sdm(:,:,:,I_QHG_sdm)
       !$omp end workshare
    else
       !$omp workshare
       Qe(:,:,:,I_HI) = 0.0_RP
       Qe(:,:,:,I_HS) = 0.0_RP
       Qe(:,:,:,I_HG) = 0.0_RP
       !$omp end workshare
    end if
    !$omp workshare
    Qe(:,:,:,I_HH) = 0.0_RP
    !$omp end workshare

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_qtrc2qhyd
  !-----------------------------------------------------------------------------
  !> Calculate mixing ratio of each category
  subroutine ATMOS_PHY_MP_sdm_qtrc2nhyd( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       Ne     )
#ifndef DISABLE_SDM
    use scale_precision
    use m_sdm_coordtrans, only: &
       sdm_x2ri, sdm_y2rj
    use m_sdm_common
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(out) :: Ne   (KA,IA,JA,N_HYD) ! number concentratio [1/kg]

#ifndef DISABLE_SDM
    real(RP) :: drate
    integer  :: ilist_l(1:sdnum_s2c)
    integer  :: tlist_l, lcnt
    integer  :: kl
    integer  :: ku
    integer  :: n, m, i, j, k, ic
    !---------------------------------------------------------------------------

    call sdm_x2ri(IS,IE,sdnum_s2c,sdx_s2c,sdri_s2c,sdrk_s2c)
    call sdm_y2rj(JS,JE,sdnum_s2c,sdy_s2c,sdrj_s2c,sdrk_s2c)

    tlist_l = 0

    lcnt = 0

    do n=1,sdnum_s2c

       if( sdrk_s2c(n)<VALID2INVALID ) cycle

       lcnt = lcnt + 1
       ilist_l(lcnt) = n

    end do

    tlist_l = lcnt

    !$omp workshare
    Ne(:,:,:,:) = 0
    !$omp end workshare

    if( tlist_l>0 ) then

       do m=1,tlist_l
          n = ilist_l(m)

          i = floor(sdri_s2c(n))+1
          j = floor(sdrj_s2c(n))+1
          k = floor(sdrk_s2c(n))+1

          if ( sdliqice_s2c(n) == STAT_LIQ ) then
             if ( sdr_s2c(n) < 1.0E-6_RP ) cycle ! do not count aerosol particles
             if ( sdr_s2c(n) < sdm_rqc2qr ) then
                ic = I_HC
             else
                ic = I_HR
             end if
          else if ( sdliqice_s2c(n) == STAT_ICE ) then
             if ( (sdice_s2c%re(n) < 1.0E-6_RP) .or. (sdice_s2c%rp(n) < 1.0E-6_RP) ) cycle ! do not count aerosol particles
             if ( sdice_s2c%rho(n) * F_THRD * ONE_PI * sdice_s2c%re(n)**2 * sdice_s2c%rp(n) * 0.3_RP < sdice_s2c%mrime(n) ) then
                ic = I_HG
             else if ( sdice_s2c%nmono(n) > 10 ) then
                ic = I_HS
             else
                ic = I_HI
             end if
          end if

          Ne(k,i,j,ic) = Ne(k,i,j,ic) + sdn_s2c(n)
       end do

       !=== correction at the verical boundary ===!

       !$omp parallel do &
       !$omp private(kl,drate,ku)
       do j = JS, JE
       do i = IS, IE

          !! at lower boundary

          kl    = floor(sdrkl_s2c(i,j))+1
          drate = real(kl,kind=RP) - sdrkl_s2c(i,j)
          if( drate<0.50_RP ) then
             Ne(kl,i,j,:) = 0.0_RP           !! <50% in share
          else
             Ne(kl,i,j,:) = Ne(kl,i,j,:) / drate
          end if

          !! at upper boundary

          ku    = floor(sdrku_s2c(i,j))+1
          drate = sdrku_s2c(i,j) - real(ku-1,kind=RP)
          if( drate<0.50_RP ) then
             Ne(ku,i,j,:) = 0.0_RP           !! <50% in share
          else
             Ne(ku,i,j,:) = Ne(ku,i,j,:) / drate
          end if

       end do
       end do

    end if

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_qtrc2nhyd
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_read( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, RHOT, QTRC )
#ifndef DISABLE_SDM
    use rng_uniform_mt, only: &
         rng_load_state
    use m_sdm_io, only: &
         sdm_sdrestart_read
    use m_sdm_common, only: &
         rng_s2c,RANDOM_IN_BASENAME, &
         sdm_dtcmph,     &
         sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
         sdm_inisdnc,sdm_zlower,             &
         sdm_zupper,sdm_calvar,sd_rest_flg_in
    use m_sdm_ptlformation, only: &
         sdm_iniset
#endif
    implicit none

    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(inout) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(inout) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]

#ifndef DISABLE_SDM
    !### Get random generator seed ###!
    !! Random number generator has already been initialized in scale-les/src/preprocess/mod_mkinit.f90
    !! Be careful. If unit (=fid_random_i) is specified, filename is ignored and the object is initialized by the unit.
    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME))
!    call rng_load_state( rng_s2c, trim(RANDOM_IN_BASENAME), fid_random_i )

    if( sd_rest_flg_in ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_restart_read",*) 'SD read from restart'
       call sdm_sdrestart_read()
    else
       LOG_INFO("ATMOS_PHY_MP_sdm_restart_read",*) 'SD initialized by sdm_iniset'
       call sdm_iniset( &
                      KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                      DENS, RHOT, QTRC,                   &
                      sdm_dtcmph,     &
                      sdm_rdnc,sdm_sdnmlvol,sdm_aslset,   &
                      sdm_inisdnc,sdm_zlower,             &
                      sdm_zupper,sdm_calvar)
    end if

    initialized = .true.

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_restart_read
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_restart_write(otime)
#ifndef DISABLE_SDM
    use scale_time
    use rng_uniform_mt, only: &
         rng_save_state
    use m_sdm_io, only: &
         sdm_sdrestart_write
    use m_sdm_common, only: &
         rng_s2c,RANDOM_OUT_BASENAME,fid_random_o
#endif
    implicit none

    real(DP), intent(in) :: otime

#ifndef DISABLE_SDM
    character(len=17) :: fmt2="(A, '.', A, I*.*)"
    character(len=17) :: fmt3="(3A)"
    character(len=H_LONG) :: ftmp
    character(len=H_LONG) :: basename_random
    character(len=19) :: basename_time

    !--- output restart file of Super Droplet
    call sdm_sdrestart_write(otime)

    !--- output restart file of Random number
    call TIME_gettimelabel(basename_time)
    write(fmt2(14:14),'(I1)') 6
    write(fmt2(16:16),'(I1)') 6
    if( RANDOM_OUT_BASENAME == '' ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_restart_write",*) 'Not found random number output file name. Default used..'
       write(ftmp,fmt3) 'random_number_output', '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    else
       write(ftmp,fmt3) trim(RANDOM_OUT_BASENAME), '_', trim(basename_time)
       write(basename_random,fmt2) trim(ftmp),'pe',mype
    endif

    fid_random_o = IO_get_available_fid()
    LOG_INFO("ATMOS_PHY_MP_sdm_restart_write",*) 'Output random number for SDM ***', trim(basename_random)
!    call rng_save_state( rng_s2c, trim(basename_random))
    call rng_save_state( rng_s2c, trim(basename_random))
!    call rng_save_state( rng_s2c, trim(basename_random), fid_random_o )

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_restart_write
  !-------------------------------------------------------------
  ! generate random number set for SDM
  subroutine ATMOS_PHY_MP_sdm_random_setup
#ifndef DISABLE_SDM
    use rng_uniform_mt,only: &
         rng_init, &
         rng_save_state
    use m_sdm_common, only: &
         rng_s2c
#endif
    implicit none

#ifndef DISABLE_SDM
    character(len=H_LONG) :: basename = ''
    character(len=H_LONG) :: RANDOM_INIT_BASENAME = '' ! name of randon number
    integer :: RANDOM_NUMBER_SEED = 0
    integer :: fid_output, ierr, iseed
    character(len=17) :: fmt="(A, '.', A, I*.*)"

    NAMELIST / PARAM_SDMRANDOM / &
      RANDOM_INIT_BASENAME, &
      RANDOM_NUMBER_SEED

    fid_output = IO_get_available_fid()

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_SDMRANDOM,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("ATMOS_PHY_MP_sdm_random_setup",*) 'Not found namelist. Check!'
       call PRC_abort
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ATMOS_PHY_MP_sdm_random_setup",*) 'Not appropriate names in namelist SDMRANDOM. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_SDMRANDOM)

    if( RANDOM_INIT_BASENAME == '' ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_random_setup",*) 'Not found random number output file name. Default used..'
       write(fmt(14:14),'(I1)') 6
       write(fmt(16:16),'(I1)') 6
       write(basename,fmt) 'random_number_output','pe',mype
    else
       write(fmt(14:14),'(I1)') 6
       write(fmt(16:16),'(I1)') 6
       write(basename,fmt) trim(RANDOM_INIT_BASENAME),'pe',mype
    endif

    iseed = mype+RANDOM_NUMBER_SEED
    call rng_init( rng_s2c, iseed )
    call rng_save_state( rng_s2c, basename )

#endif
    return
  end subroutine ATMOS_PHY_MP_sdm_random_setup
  !-----------------------------------------------------------------------------
  subroutine ATMOS_PHY_MP_sdm_additional_output( &
       KA, KS, KE, IA, IS, IE, JA, JS, JE, &
       DENS, RHOT, QTRC )
#ifndef DISABLE_SDM
    use scale_file_history, only: &
       FILE_HISTORY_query, &
       FILE_HISTORY_in,    &
       FILE_HISTORY_put
    use m_sdm_fluidconv, only: &
       sdm_rhot_qtrc2p_t
    use m_sdm_sd2fluid, only: &
       sdm_sd2rhosd,sdm_sd2rhodropmom
    use m_sdm_io, only: &
       sdm_outasci,sdm_outnetcdf,sdm_outnetcdf_hist
    use m_sdm_idutil, only: &
         sdm_copy_selected_sd
    use m_sdm_coordtrans, only: &
       sdm_rk2z
    use m_sdm_motion, only: &
       sdm_getvz_liq, sdm_getvz_ice
    use m_sdm_common
#endif
    implicit none
    integer, intent(in) :: KA, KS, KE
    integer, intent(in) :: IA, IS, IE
    integer, intent(in) :: JA, JS, JE

    real(RP), intent(in) :: DENS(KA,IA,JA)        !! Density [kg/m3]
    real(RP), intent(in) :: RHOT(KA,IA,JA)        !! DENS * POTT [K*kg/m3]
    real(RP), intent(in) :: QTRC(KA,IA,JA,QA_MP_sdm)    !! Ratio of mass of tracer to total mass[kg/kg]

#ifndef DISABLE_SDM
    real(RP), pointer :: sdx_tmp(:),sdy_tmp(:),sdrk_tmp(:),sdz_tmp(:),sdri_tmp(:),sdrj_tmp(:),sdr_tmp(:),sdvz_tmp(:)
    real(RP), pointer :: sdasl_tmp(:,:)
    integer(i2), pointer :: sdliqice_tmp(:)
    integer(DP), pointer :: sdn_tmp(:)
    type(sdicedef), pointer :: sdice_tmp
    integer :: sdnum_tmp, sdnumasl_tmp
    logical :: do_puthist, do_puthist_0, do_puthist_1, do_puthist_2, do_puthist_3
    integer :: order_n

    real(RP) :: pres_scale(KA,IA,JA) ! Pressure
    real(RP) :: t_scale(KA,IA,JA)    ! Temperature
    real(RP) :: rhosd(KA,IA,JA) ! sd number density
    real(RP) :: rhodropmom0(KA,IA,JA) ! density of the 0th droplet moment
    real(RP) :: rhodropmom1(KA,IA,JA) ! density of the 1st droplet moment
    real(RP) :: rhodropmom2(KA,IA,JA) ! density of the 2nd droplet moment
    real(RP) :: rhodropmom3(KA,IA,JA) ! density of the 3rd droplet moment

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_start("sdm_out",0,0)
#endif

    !!!!!!!!!!!!!!!!!!!!!
    !!!! History output of Accumulated Precipitation Amount

    call FILE_HISTORY_in( prr_crs(:,:,2),       'RAIN_ACC_sd',    'surface rain accumulation',     'm')
    if( sdm_cold ) then
       call FILE_HISTORY_in( prr_crs(:,:,4),       'SNOW_ACC_sd',    'surface snow accumulation',    'm')
    end if
    call FILE_HISTORY_in( prr_crs(:,:,6),       'PREC_ACC_sd',    'surface precipitation accumulation',    'm')

    !!!!!!!!!!!!!!!!!!!!!
    !!!! Output of SD Data

    if( (mod(sdm_dmpvar,10)==1) .and. sdm_dmpitva>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitva+0.00010_RP))) == 0 ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_sd_write",*) 'Output Super Droplet Data in ASCII'
       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z( IS, IE, JS, JE, &
                      sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
       !!! terminal velocity vz
       call sdm_rhot_qtrc2p_t( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                               RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz_liq( KA, IA,IS,IE, JA,JS,JE, &
                           pres_scale,DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
       if( sdm_cold )then
          call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c,sdvz_s2c,  &
                           sd_itmp1,'no_interpolation' )
       end if
       !! Output
       call sdm_outasci(                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdice_s2c, &
                        sdm_dmpnskip)
    end if

    if( ((mod(sdm_dmpvar,100))/10>=1) .and. sdm_dmpitvb>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitvb+0.00010_RP))) == 0 ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_sd_write",*) 'Output Super Droplet Data in NetCDF'
       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z( IS, IE, JS, JE, &
                      sdnum_s2c,sdx_s2c,sdy_s2c,sdrk_s2c,sdz_s2c,sdri_s2c,sdrj_s2c)
       !!! terminal velocity vz
       ! sdvz_s2c(:)=0.0_RP ! temporary for debug
       call sdm_rhot_qtrc2p_t( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                               RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz_liq( KA, IA,IS,IE, JA,JS,JE, &
                           pres_scale,DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdr_s2c,sdvz_s2c,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
       if( sdm_cold )then
          call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sdnum_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdice_s2c,sdvz_s2c,  &
                           sd_itmp1,'no_interpolation' )
       end if
       !! Output
       if( (mod(sdm_dmpvar,100))/10==1) then
          call sdm_outnetcdf(                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdice_s2c, &
                        filetag='all')
       else if( (mod(sdm_dmpvar,100))/10==2) then
          call sdm_outnetcdf_hist(TIME_NOWSEC,                               &
                        sdnum_s2c,sdnumasl_s2c,                    &
                        sdn_s2c,sdliqice_s2c,sdx_s2c,sdy_s2c,sdz_s2c,sdr_s2c,sdasl_s2c,sdvz_s2c, &
                        sdice_s2c, &
                        filetag='all')
       end if
    end if

    if( ((mod(sdm_dmpvar,1000))/100>=1) .and. sdm_dmpitvl>0.0_RP .and. &
         mod(10*int(1.E+2_RP*(TIME_NOWSEC+0.0010_RP)), &
             int(1.E+3_RP*(sdm_dmpitvl+0.00010_RP))) == 0 ) then
       LOG_INFO("ATMOS_PHY_MP_sdm_sd_write",*) ' *** Output Selected Super Droplet Data in NetCDF'

       sdx_tmp  => sd_dtmp1
       sdy_tmp  => sd_dtmp2
       sdrk_tmp => sd_dtmp3
       sdz_tmp  => sd_dtmp4 ! diagnostic variable
       sdri_tmp => sd_dtmp5 ! diagnostic variable (for now)
       sdrj_tmp => sd_dtmp6 ! diagnostic variable (for now)
       sdr_tmp  => sd_dtmp7
       sdvz_tmp => sd_dtmp8 ! diagnostic variable
       sdliqice_tmp => sd_i2tmp1
       sdice_tmp=> sd_icetmp1
       sdn_tmp  => sd_i8tmp1
       sdasl_tmp=> sd_asltmp1

       call sdm_rhot_qtrc2p_t( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                               RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_copy_selected_sd(KA, IA,IS,IE, JA,JS,JE, &
            &                    sdnum_s2c,sdnumasl_s2c,sdn_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c, &
            &                    sdliqice_s2c,sdasl_s2c,sdr_s2c,sdice_s2c,                                  &
            &                    sdnum_tmp,sdnumasl_tmp,sdn_tmp,sdx_tmp,sdy_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp, &
            &                    sdliqice_tmp,sdasl_tmp,sdr_tmp,sdice_tmp,                                  &
            &                    t_scale,sd_itmp1,sdtype='activated') ! options: 'all', 'large', 'activated'

       !! Evaluate diagnostic variables
       !!! z
       call sdm_rk2z( IS, IE, JS, JE, &
                      sdnum_tmp,sdx_tmp,sdy_tmp,sdrk_tmp,sdz_tmp,sdri_tmp,sdrj_tmp)
       !!! terminal velocity vz
       ! sdvz_tmp(:)=0.0_RP ! temporary for debug
       call sdm_rhot_qtrc2p_t( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                               RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_getvz_liq( KA, IA,IS,IE, JA,JS,JE, &
                           pres_scale,DENS,t_scale,            &
                           sdnum_tmp,sdliqice_tmp,sdx_tmp,sdy_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdr_tmp,sdvz_tmp,  &
                           sd_itmp1,sd_itmp2,sd_itmp3,'no_interpolation' )
       if( sdm_cold )then
          call sdm_getvz_ice( KA, IA,IS,IE, JA,JS,JE, &
                           DENS,t_scale,            &
                           sdnum_tmp,sdliqice_tmp,sdx_tmp,sdy_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdice_tmp,sdvz_tmp,  &
                           sd_itmp1,'no_interpolation' )
       end if

       !! Output
       if( (mod(sdm_dmpvar,1000))/100==1) then
          call sdm_outnetcdf(                               &
                        sdnum_tmp,sdnumasl_tmp,                    &
                        sdn_tmp,sdliqice_tmp,sdx_tmp,sdy_tmp,sdz_tmp,sdr_tmp,sdasl_tmp,sdvz_tmp, &
                        sdice_tmp, &
                        filetag='selected')
       else if( (mod(sdm_dmpvar,1000))/100==2) then
          call sdm_outnetcdf_hist(TIME_NOWSEC,                               &
                        sdnum_tmp,sdnumasl_tmp,                    &
                        sdn_tmp,sdliqice_tmp,sdx_tmp,sdy_tmp,sdz_tmp,sdr_tmp,sdasl_tmp,sdvz_tmp, &
                        sdice_tmp, &
                        filetag='selected')
       end if

       nullify(sdx_tmp)
       nullify(sdy_tmp)
       nullify(sdrk_tmp)
       nullify(sdz_tmp)
       nullify(sdri_tmp)
       nullify(sdrj_tmp)
       nullify(sdr_tmp)
       nullify(sdvz_tmp)
       nullify(sdliqice_tmp)
       nullify(sdice_tmp)
       nullify(sdn_tmp)
       nullify(sdasl_tmp)

    end if

    !!!!!!!!!!!!!!!!!!!!!
    !!!! History output of SNC (SD number density).
    !!!! Data is updated only at the time of output

    ! Check whether it is time to input the item
    call FILE_HISTORY_query( histitemid(1), & ! [IN] ! SNC
                             do_puthist     ) ! [OUT]

    if( do_puthist ) then
       call sdm_sd2rhosd(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
            zph_crs,rhosd,sdnum_s2c,sdn_s2c,  &
            sdx_s2c,sdy_s2c,sdr_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c,sdrkl_s2c,sdrku_s2c,      &
            sd_itmp1)
       call FILE_HISTORY_put( histitemid(1), rhosd(:,:,:) )
    endif

    !!!!!!!!!!!!!!!!!!!!!
    !!!! History output of the density of droplet moments.
    !!!! Data is updated only at the time of output

    do_puthist = .false.
    ! Check whether it is time to input the item
    call FILE_HISTORY_query( histitemid(2),  & ! [IN] ! NC
                             do_puthist_0   ) ! [OUT]
    if( do_puthist_0 ) do_puthist = .true.

    ! Check whether the item has been already registered
    do_puthist_1 = .false.
    call FILE_HISTORY_query( histitemid(3), & ! [IN] ! DMOM1
                             do_puthist_1   ) ! [OUT]
    if( do_puthist_1 ) do_puthist = .true.

    ! Check whether it is time to input the item
    do_puthist_2 = .false.
    call FILE_HISTORY_query( histitemid(4), & ! [IN] ! DMOM2
                             do_puthist_2   ) ! [OUT]
    if( do_puthist_2 ) do_puthist = .true.

    ! Check whether it is time to input the item
    do_puthist_3 = .false.
    call FILE_HISTORY_query( histitemid(5), & ! [IN] ! DMOM3
                             do_puthist_3   ) ! [OUT]
    if( do_puthist_3 ) do_puthist = .true.

    if( do_puthist ) then
       sdx_tmp  => sd_dtmp1
       sdy_tmp  => sd_dtmp2
       sdrk_tmp => sd_dtmp3
       sdz_tmp  => sd_dtmp4 ! diagnostic variable
       sdri_tmp => sd_dtmp5 ! diagnostic variable (for now)
       sdrj_tmp => sd_dtmp6 ! diagnostic variable (for now)
       sdr_tmp  => sd_dtmp7
       sdvz_tmp => sd_dtmp8 ! diagnostic variable
       sdliqice_tmp => sd_i2tmp1
       sdice_tmp=> sd_icetmp1
       sdn_tmp  => sd_i8tmp1
       sdasl_tmp=> sd_asltmp1

       call sdm_rhot_qtrc2p_t( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
                               RHOT,QTRC,DENS,pres_scale,t_scale)
       call sdm_copy_selected_sd(KA, IA,IS,IE, JA,JS,JE, &
            &                    sdnum_s2c,sdnumasl_s2c,sdn_s2c,sdx_s2c,sdy_s2c,sdri_s2c,sdrj_s2c,sdrk_s2c, &
            &                    sdliqice_s2c,sdasl_s2c,sdr_s2c,sdice_s2c,                                  &
            &                    sdnum_tmp,sdnumasl_tmp,sdn_tmp,sdx_tmp,sdy_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp, &
            &                    sdliqice_tmp,sdasl_tmp,sdr_tmp,sdice_tmp,                                  &
            &                    t_scale,sd_itmp1,sdtype='activated') ! options: 'all', 'large', 'activated'

       if( do_puthist_0 )then
          order_n = 0
          call sdm_sd2rhodropmom(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               & order_n,zph_crs,rhodropmom0,sdnum_tmp,sdn_tmp,sdliqice_tmp, &
               & sdx_tmp,sdy_tmp,sdr_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdrkl_s2c,sdrku_s2c,    &
               & sd_itmp1)
          call FILE_HISTORY_put( histitemid(2), rhodropmom0(:,:,:) )
       end if

       if( do_puthist_1 )then
          order_n = 1
          call sdm_sd2rhodropmom(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               & order_n,zph_crs,rhodropmom1,sdnum_tmp,sdn_tmp,sdliqice_tmp, &
               & sdx_tmp,sdy_tmp,sdr_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdrkl_s2c,sdrku_s2c,    &
               & sd_itmp1)
          call FILE_HISTORY_put( histitemid(3), rhodropmom1(:,:,:) )
       end if

       if( do_puthist_2 )then
          order_n = 2
          call sdm_sd2rhodropmom(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               & order_n,zph_crs,rhodropmom2,sdnum_tmp,sdn_tmp,sdliqice_tmp, &
               & sdx_tmp,sdy_tmp,sdr_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdrkl_s2c,sdrku_s2c,    &
               & sd_itmp1)
          call FILE_HISTORY_put( histitemid(4), rhodropmom2(:,:,:) )
       end if

       if( do_puthist_3 )then
          order_n = 3
          call sdm_sd2rhodropmom(KA, KS, KE, IA, IS, IE, JA, JS, JE, &
               & order_n,zph_crs,rhodropmom3,sdnum_tmp,sdn_tmp,sdliqice_tmp, &
               & sdx_tmp,sdy_tmp,sdr_tmp,sdri_tmp,sdrj_tmp,sdrk_tmp,sdrkl_s2c,sdrku_s2c,    &
               & sd_itmp1)
          call FILE_HISTORY_put( histitemid(5), rhodropmom3(:,:,:) )
       end if

       nullify(sdx_tmp)
       nullify(sdy_tmp)
       nullify(sdrk_tmp)
       nullify(sdz_tmp)
       nullify(sdri_tmp)
       nullify(sdrj_tmp)
       nullify(sdr_tmp)
       nullify(sdvz_tmp)
       nullify(sdliqice_tmp)
       nullify(sdice_tmp)
       nullify(sdn_tmp)
       nullify(sdasl_tmp)

    endif

#ifdef _FAPP_
    ! Section specification for fapp profiler
    call fapp_stop("sdm_out",0,0)
#endif

#endif
  end subroutine ATMOS_PHY_MP_sdm_additional_output
end module scale_atmos_phy_mp_sdm
!-------------------------------------------------------------------------------
