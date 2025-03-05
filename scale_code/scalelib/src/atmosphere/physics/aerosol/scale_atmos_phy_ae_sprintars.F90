!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prof
  use scale_prc, only: &
       PRC_myrank, &
       PRC_abort
  use scale_statistics, only: &
       STATISTICS_detail
  use scale_atmos_grid_cartesC_index, only: &
       KA, KS, KE, &
       IA, IS, IE, &
       JA, JS, JE
  use scale_file_cartesC, only: &
       FILE_CARTESC_write_var
  use scale_atmos_aerosol, only: &
       N_AE
  use scale_atmos_phy_ae_SPRINTARS_admin, only: &
       ATMOS_PHY_AE_SPRINTARS_AEROSL_setup, &
       ATMOS_PHY_AE_SPRINTARS_AEROSL_finalize, &
       ATMOS_PHY_AE_SPRINTARS_AEROSL
  use scale_atmos_phy_ae_SPRINTARS_common, only: &
       ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA, &
       NSU, NCA, NDU, NSA, NAT, &
       QA_AE_SU, QS_AE_SU, QE_AE_SU, &
       QA_AE_CA, QS_AE_CA, QE_AE_CA, &
       QA_AE_DU, QS_AE_DU, QE_AE_DU, &
       QA_AE_SA, QS_AE_SA, QE_AE_SA, &
       QA_AE_AT, QS_AE_AT, QE_AE_AT
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_tracer_setup
  public :: ATMOS_PHY_AE_SPRINTARS_setup
  public :: ATMOS_PHY_AE_SPRINTARS_finalize
  public :: ATMOS_PHY_AE_SPRINTARS_mkinit
  public :: ATMOS_PHY_AE_SPRINTARS_mkemission
  public :: ATMOS_PHY_AE_SPRINTARS_negative_fixer
  public :: ATMOS_PHY_AE_SPRINTARS_tendency
  
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_AE_SPRINTARS_NAME (:) !> name
  character(len=H_MID)  , public, allocatable :: ATMOS_PHY_AE_SPRINTARS_DESC (:) !> description
  character(len=H_SHORT), public, allocatable :: ATMOS_PHY_AE_SPRINTARS_UNIT (:) !> unit
  real(RP),               public, allocatable :: ATMOS_PHY_AE_SPRINTARS_CV   (:) !> specific heat capacity at constant volume
  real(RP),               public, allocatable :: ATMOS_PHY_AE_SPRINTARS_CP   (:) !> specific heat capacity at constant pressure
  real(RP),               public, allocatable :: ATMOS_PHY_AE_SPRINTARS_R    (:) !> Air constant
  real(RP),               public, allocatable :: ATMOS_PHY_AE_SPRINTARS_ENGI0(:) !> internal energy offset at 0K
  logical,                public, allocatable :: ATMOS_PHY_AE_SPRINTARS_ADVC (:) !> to be advected in the dynamical core
  logical,                public, allocatable :: ATMOS_PHY_AE_SPRINTARS_MASS (:) !> 1 for tracers with mass, otherwise 0
  
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

  
  
  !-----------------------------------------------------------------------------
  !> Tracer setup
  subroutine ATMOS_PHY_AE_SPRINTARS_tracer_setup( QA_AE )
    implicit none

    integer, intent(out) :: QA_AE
    integer   :: iq
    character :: cnumber*2
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_SPRINTARS_tracer_setup",*) 'Setup'

    ! # of tracers
    QA_AE = 0
    if ( NSU /= 0 ) then
       QS_AE_SU = QA_AE + 1
       QA_AE_SU = NSU
       QE_AE_SU = QS_AE_SU + QA_AE_SU - 1
       QA_AE = QA_AE + NSU
    endif
    if ( NCA /= 0 ) then
       QS_AE_CA = QA_AE + 1
       QA_AE_CA = NCA
       QE_AE_CA = QS_AE_CA + QA_AE_CA - 1
       QA_AE = QA_AE + NCA
    endif
    if ( NDU /= 0 ) then
       QS_AE_DU = QA_AE + 1
       QA_AE_DU = NDU
       QE_AE_DU = QS_AE_DU + QA_AE_DU - 1
       QA_AE = QA_AE + NDU
    endif
    if ( NSA /= 0 ) then
       QS_AE_SA = QA_AE + 1
       QA_AE_SA = NSA
       QE_AE_SA = QS_AE_SA + QA_AE_SA - 1
       QA_AE = QA_AE + NSA
    endif
    if ( NAT /= 0 ) then
       QS_AE_AT = QA_AE + 1
       QA_AE_AT = NAT
       QE_AE_AT = QS_AE_AT + QA_AE_AT - 1
       QA_AE = QA_AE + NAT
    endif

    ! assign tracer names
    allocate( ATMOS_PHY_AE_SPRINTARS_NAME (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_DESC (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_UNIT (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_CV   (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_CP   (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_R    (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_ENGI0(QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_ADVC (QA_AE) )
    allocate( ATMOS_PHY_AE_SPRINTARS_MASS (QA_AE) )

    if ( QA_AE_SU /= 0 ) then
       do iq = QS_AE_SU, QE_AE_SU, 1
          write(cnumber,'(i2.2)') iq - QS_AE_SU + 1
          ATMOS_PHY_AE_SPRINTARS_NAME (iq) = 'sulfate' // cnumber
          if     ( iq - QS_AE_SU + 1 == 1 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of sulfate' // cnumber // ' (SO4) aerosol mass to total mass'
          elseif ( iq - QS_AE_SU + 1 == 2 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of sulfate' // cnumber // ' (SO2) mass to total mass'
          elseif ( iq - QS_AE_SU + 1 == 3 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of sulfate' // cnumber // ' (DMS) mass to total mass'
          else
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of sulfate' // cnumber // ' aerosol mass to total mass'
          endif
          ATMOS_PHY_AE_SPRINTARS_UNIT (iq) = 'kg/kg'
          ATMOS_PHY_AE_SPRINTARS_CV   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_CP   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_R    (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ENGI0(iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ADVC (iq) = .true.
          ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .true.
          !ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .false.
       enddo
    endif
    
    if ( QA_AE_CA /= 0 ) then
       do iq = QS_AE_CA, QE_AE_CA, 1
          write(cnumber,'(i2.2)') iq - QS_AE_CA + 1
          ATMOS_PHY_AE_SPRINTARS_NAME (iq) = 'carbon' // cnumber
          if     ( iq - QS_AE_CA + 1 == 1 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (OC;BB+FF;internal mixture) aerosol mass to total mass'
          elseif ( iq - QS_AE_CA + 1 == 2 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (BC;BB+FF;internal mixture) aerosol mass to total mass'
          elseif ( iq - QS_AE_CA + 1 == 3 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (secondary organic carbon (SOC) from NVOC) aerosol mass to total mass'
          elseif ( iq - QS_AE_CA + 1 == 4 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (BC;FF;external mixture) aerosol mass to total mass'
          elseif ( iq - QS_AE_CA + 1 == 5 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (terpene) mass to total mass'
          elseif ( iq - QS_AE_CA + 1 == 6 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' (isoprene) mass to total mass'
          else
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of carbon' // cnumber // ' aerosol mass to total mass'
          endif
          ATMOS_PHY_AE_SPRINTARS_UNIT (iq) = 'kg/kg'
          ATMOS_PHY_AE_SPRINTARS_CV   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_CP   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_R    (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ENGI0(iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ADVC (iq) = .true.
          ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .true.
       enddo
    endif
    
    if ( QA_AE_DU /= 0 ) then
       do iq = QS_AE_DU, QE_AE_DU, 1
          write(cnumber,'(i2.2)') iq - QS_AE_DU + 1
          ATMOS_PHY_AE_SPRINTARS_NAME (iq) = 'dust' // cnumber
          if     ( iq - QS_AE_DU + 1 == 1 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (0.12-0.22um) aerosol mass to total mass'
          elseif ( iq - QS_AE_DU + 1 == 2 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (0.22-0.46um) aerosol mass to total mass'
          elseif ( iq - QS_AE_DU + 1 == 3 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (0.46-1.00um) aerosol mass to total mass'
          elseif ( iq - QS_AE_DU + 1 == 4 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (1.00-2.15um) aerosol mass to total mass'
          elseif ( iq - QS_AE_DU + 1 == 5 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (2.15-4.64um) aerosol mass to total mass'
          elseif ( iq - QS_AE_DU + 1 == 6 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' (4.64-10.00um) aerosol mass to total mass'
          else
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of dust' // cnumber // ' aerosol mass to total mass'
          endif
          ATMOS_PHY_AE_SPRINTARS_UNIT (iq) = 'kg/kg'
          ATMOS_PHY_AE_SPRINTARS_CV   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_CP   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_R    (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ENGI0(iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ADVC (iq) = .true.
          ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .true.
       enddo
    endif
    
    if ( QA_AE_SA /= 0 ) then
       do iq = QS_AE_SA, QE_AE_SA, 1
          write(cnumber,'(i2.2)') iq - QS_AE_SA + 1
          ATMOS_PHY_AE_SPRINTARS_NAME (iq) = 'salt' // cnumber
          if     ( iq - QS_AE_SA + 1 == 1 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of salt' // cnumber // ' (0.100-0.316um) aerosol mass to total mass'
          elseif ( iq - QS_AE_SA + 1 == 2 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of salt' // cnumber // ' (0.316-1.000um) aerosol mass to total mass'
          elseif ( iq - QS_AE_SA + 1 == 3 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of salt' // cnumber // ' (1.000-3.160um) aerosol mass to total mass'
          elseif ( iq - QS_AE_SA + 1 == 4 ) then
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of salt' // cnumber // ' (3.160-10.000um) aerosol mass to total mass'
          else
             ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of salt' // cnumber // ' aerosol mass to total mass'
          endif
          ATMOS_PHY_AE_SPRINTARS_UNIT (iq) = 'kg/kg'
          ATMOS_PHY_AE_SPRINTARS_CV   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_CP   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_R    (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ENGI0(iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ADVC (iq) = .true.
          ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .true.
       enddo
    endif
    
    if ( QA_AE_AT /= 0 ) then
       do iq = QS_AE_AT, QE_AE_AT, 1
          write(cnumber,'(i2.2)') iq - QS_AE_AT + 1
          ATMOS_PHY_AE_SPRINTARS_NAME (iq) = 'tracer' // cnumber
          ATMOS_PHY_AE_SPRINTARS_DESC (iq) = 'Ratio of tracer' // cnumber // ' aerosol mass to total mass'
          ATMOS_PHY_AE_SPRINTARS_UNIT (iq) = 'kg/kg'
          ATMOS_PHY_AE_SPRINTARS_CV   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_CP   (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_R    (iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ENGI0(iq) = 0.0_RP
          ATMOS_PHY_AE_SPRINTARS_ADVC (iq) = .true.
          ATMOS_PHY_AE_SPRINTARS_MASS (iq) = .true.
       enddo
    endif

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_tracer_setup
  !-----------------------------------------------------------------------------

  

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_SPRINTARS_setup
    implicit none
    !----------------------------------------------------
  
    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_SPRINTARS_setup",*) 'Setup'

    call ATMOS_PHY_AE_SPRINTARS_AEROSL_setup
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_setup
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> Setup
  subroutine ATMOS_PHY_AE_SPRINTARS_finalize
    implicit none
    !----------------------------------------------------
  
    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_SPRINTARS_finalize",*) 'Finalize'

    call ATMOS_PHY_AE_SPRINTARS_AEROSL_finalize
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_finalize
  !-----------------------------------------------------------------------------
  
  
  !-----------------------------------------------------------------------------
  !> Mkinit
  subroutine ATMOS_PHY_AE_SPRINTARS_mkinit( QA_AE,    & ! (in)
                                            ccn_init, &
                                            QTRC,     & ! (inout)
                                            CCN       ) ! (inout)
    implicit none

    integer,  intent(in)    :: QA_AE
    real(RP), intent(in)    :: ccn_init
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE)
    real(RP), intent(inout) :: CCN (KA,IA,JA)
    integer :: k, i, j, iq
    !----------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("ATMOS_PHY_AE_SPRINTARS_mkinit",*) 'Mkinit'

    !QTRC
    !---------------------------------
    ! do iq = QS_AE_AT, QE_AE_AT, 1
    !    do j = JS, JE
    !    do i = IS, IE
    !       do k = KS, KE
    !          if ( i >= int(real(IE-IS,kind=RP)/2.0_RP-real(IE-IS,kind=RP)/6.0_RP) .and. &
    !               i <= int(real(IE-IS,kind=RP)/2.0_RP+real(IE-IS,kind=RP)/6.0_RP) .and. &
    !               j >= int(real(JE-JS,kind=RP)/2.0_RP-real(JE-JS,kind=RP)/6.0_RP) .and. &
    !               j <= int(real(JE-JS,kind=RP)/2.0_RP+real(JE-JS,kind=RP)/6.0_RP) ) then
    !             QTRC(k,i,j,iq) = 2.0E-10_RP
    !          else
    !             QTRC(k,i,j,iq) = 0.0E-10_RP
    !          endif
    !       enddo
    !    enddo
    !    enddo
    ! enddo

    do iq = QS_AE_AT, QE_AE_AT, 1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRC(k,i,j,iq) = 0.0_RP
          enddo
       enddo
       enddo
    enddo
    
    !CCN
    !---------------------------------
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          CCN(k,i,j) = ccn_init
       enddo
    enddo
    enddo


    !Negative fixer
    !---------------------------------
    do iq = 1, QA_AE, 1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( QTRC(k,i,j,iq) < 0.0_RP ) QTRC(k,i,j,iq) = 0.0_RP
          enddo
       enddo
       enddo
    enddo
    do j = JS, JE
    do i = IS, IE
       do k = KS, KE
          if ( CCN(k,i,j) < 0.0_RP ) CCN(k,i,j) = 0.0_RP
       enddo
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_mkinit
  !-----------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------
  !> Mkinit
  subroutine ATMOS_PHY_AE_SPRINTARS_mkemission
    use scale_file_cartesC, only: &
       FILE_CARTESC_create, & 
       FILE_CARTESC_def_var, &
       FILE_CARTESC_enddef, &
       FILE_CARTESC_write_var, &
       FILE_CARTESC_close
    use scale_atmos_grid_cartesC_real, only: &
       GDZ      => ATMOS_GRID_CARTESC_REAL_CZ,  & !< GDZ (KA,  IA,JA) : altitude [m]
       GDZM     => ATMOS_GRID_CARTESC_REAL_FZ,  & !< GDZM(0:KA,IA,JA) : altitude (half level) [m]
       LON_REAL => ATMOS_GRID_CARTESC_REAL_LON, & !< longitude [rad,0-2pi]
       LAT_REAL => ATMOS_GRID_CARTESC_REAL_LAT    !< latitude  [rad,-pi,pi]
    use scale_time, only: &
       NOWDATE => TIME_NOWDATE
    use scale_prc, only: &
       PRC_abort
    use scale_atmos_phy_ae_SPRINTARS_common, only: &
       ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA
    implicit none

    ! File name of emission data base
    character(len=H_LONG) :: fname_land = 'sprintars_land_database.nc' !< file name (land)
    character(len=H_LONG) :: fname_chem = 'sprintars_chem_database.nc' !< file name (chem)
    character(len=H_LONG) :: fname_carb = 'sprintars_carb_database.nc' !< file name (carb)
    character(len=H_LONG) :: fname_sulf = 'sprintars_sulf_database.nc' !< file name (sulf)
    character(len=H_LONG) :: basename_emis = 'sprintarts_emission' !< file name of emission file

    character(len=H_LONG) :: type_emission  = 'EMISSION DATA for SCALE-RM with SPRINTARS'
    character(len=H_LONG) :: dtype_emission = 'DEFAULT'
    character(len=H_LONG) :: basename_out_mod

    character(len=H_SHORT) :: varname(176)
    character(len=H_SHORT) :: unit, dimtype(176)
    character(len=H_LONG)  :: longname

    integer :: fid      !< file id for output of emission file of each node
    integer :: vid(176), totalid

    integer, parameter  :: outnstep = 12
    !real(DP), parameter :: time_interval = 2592000.0_RP  ! 30 days (1month)
    real(DP), parameter :: time_interval = 1.0_RP  ! 30 days (1month)
    integer :: it, i, j, k, ierr, iq
    real(RP) :: work2d(IA,JA)

    ! Arrays for input emission data from file
    real(RP), allocatable :: GRLAI_DATA (:,:,:)
    real(RP), allocatable :: GRSNW_DATA (:,:,:,:)
    real(RP), allocatable :: OHRAD_DATA (:,:,:,:)
    real(RP), allocatable :: NO3G_DATA  (:,:,:,:) 
    real(RP), allocatable :: O3ORI_DATA (:,:,:,:)
    real(RP), allocatable :: DAYTIM_DATA(:,:,:)
    real(RP), allocatable :: BBBC_DATA  (:,:,:,:) ! 22 -> 1999-2020
    real(RP), allocatable :: BBOC_DATA  (:,:,:,:) ! 22 -> 1999-2020
    real(RP), allocatable :: ANTBC_DATA (:,:,:,:) ! 22 -> 1999-2020
    real(RP), allocatable :: ANTOC_DATA (:,:,:,:) ! 22 -> 1999-2020
    real(RP), allocatable :: ISOP_DATA  (:,:,:)
    real(RP), allocatable :: TERP_DATA  (:,:,:)
    real(RP), allocatable :: CHLOA_DATA (:,:,:)
    real(RP), allocatable :: DIFATT_DATA(:,:,:)
    real(RP), allocatable :: SULDIO_DATA(:,:,:,:) ! 14 -> 2000-2013
    real(RP), allocatable :: BBSO2_DATA (:,:,:,:) ! 14 -> 2000-2013
    real(RP), allocatable :: VOLSO2_DATA(:,:,:) !
    real(RP), allocatable :: VOLALT_DATA(:,:,:) !
    real(RP), allocatable :: ANTSO4_DATA(:,:,:) !
    real(RP), allocatable :: HOORI_DATA (:,:,:,:)

    namelist / PARAM_ATMOS_PHY_AE_SPRINTARS_INIT_NMDATA / &
         fname_land, &
         fname_chem, &
         fname_carb, &
         fname_sulf, &
         basename_emis

    !PARAM_ATMOS_PHY_AE_SPRINTARS_NMDATA
    !----------------------
    !read config file (land) file
    !-----------------
    rewind( IO_FID_CONF )
    read( IO_FID_CONF, nml=PARAM_ATMOS_PHY_AE_SPRINTARS_INIT_NMDATA, iostat=ierr )
    if( ierr < 0 ) then !--- missing
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_mkemission",*)  'Not found namelist PARAM_ATMOS_PHY_AE_SPRINTARS_INIT_NMDATA. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_mkemission",*) 'Not appropriate names in namelist PARAM_ATMOS_PHY_AE_SPRINTARS_INIT_NMDATA, Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_ATMOS_PHY_AE_SPRINTARS_INIT_NMDATA)

    !----------------------
    !read parameter (land) file
    !-----------------
    !Setup LAI (NETCDF)
    allocate( GRLAI_DATA(IA,JA,12) )
    GRLAI_DATA(:,:,:) = 0.0_RP
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( GRLAI_DATA(:,:,:),                  & ! [OUT]
           trim(adjustl(fname_land)), 'GRLAI', & ! [IN]
           LON_REAL, LAT_REAL                  ) ! [IN]

    !Setup Snow Depth (NETCDF)
    allocate( GRSNW_DATA(IA,JA,12,31) )
    GRSNW_DATA(:,:,:,:) = 0.0_RP
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( GRSNW_DATA(:,:,:,:),                & ! [OUT]
           trim(adjustl(fname_land)), 'GRSNW', & ! [IN]
           LON_REAL, LAT_REAL                  ) ! [IN]
    !-----------------

    !----------------------
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

    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  OHRAD_DATA(:,:,:,:),                & ! [OUT]
            trim(adjustl(fname_chem)), 'OHRAD', & ! [IN]
            LON_REAL, LAT_REAL,                 & ! [IN]
            GDZ, GDZM                           ) ! [IN]

    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  NO3G_DATA(:,:,:,:),                & ! [OUT]
            trim(adjustl(fname_chem)), 'NO3G',  & ! [IN]
            LON_REAL, LAT_REAL,                 & ! [IN]
            GDZ, GDZM                           ) ! [IN]

    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( O3ORI_DATA(:,:,:,:),                & ! [OUT]
           trim(adjustl(fname_chem)), 'OZONE', & ! [IN]
           LON_REAL, LAT_REAL,                 & ! [IN]
           GDZ, GDZM                           ) ! [IN]

    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( DAYTIM_DATA(:,:,:),                  & ! [OUT]
           trim(adjustl(fname_chem)), 'DAYTIM', & ! [IN]
           LON_REAL, LAT_REAL                   ) ! [IN]

    do it = 1, 12
       do j = JA, JA
       do i = IA, IA
          do k = KA, KA
             OHRAD_DATA(k,i,j,it) = max( OHRAD_DATA(k,i,j,it), 0.0_RP )
              NO3G_DATA(k,i,j,it) = max(  NO3G_DATA(k,i,j,it), 0.0_RP )
             O3ORI_DATA(k,i,j,it) = max( O3ORI_DATA(k,i,j,it), 0.0_RP )
          enddo
          DAYTIM_DATA(i,j,it) = max( DAYTIM_DATA(i,j,it), 0.0_RP )
       enddo
       enddo
    enddo
    !-----------------

    !----------------------
    !read emission file
    !-----------------
    !Setup BBBC, BBOC, ANTBC, ANTOC, ISOP, TERP, CHLOA, DIFATT
    allocate(   BBBC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(   BBOC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(  ANTBC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(  ANTOC_DATA(IA,JA,22,12) ) ! 22 -> 1999-2020
    allocate(   ISOP_DATA(IA,JA,   12) )
    allocate(   TERP_DATA(IA,JA,   12) )
    allocate(  CHLOA_DATA(IA,JA,   12) )
    allocate( DIFATT_DATA(IA,JA,   12) )
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (   BBBC_DATA(:,:,:,:), fname_carb, 'BBBC',   LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (   BBOC_DATA(:,:,:,:), fname_carb, 'BBOC',   LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  ANTBC_DATA(:,:,:,:), fname_carb, 'ANTBC',  LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  ANTOC_DATA(:,:,:,:), fname_carb, 'ANTOC',  LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (   ISOP_DATA(:,:,  :), fname_carb, 'ISOP',   LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (   TERP_DATA(:,:,  :), fname_carb, 'TERP',   LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  CHLOA_DATA(:,:,  :), fname_carb, 'CHLOA',  LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( DIFATT_DATA(:,:,  :), fname_carb, 'DIFATT', LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]

    !-----------------
    !Setup SULDIO, BBSO2, VOLSO2, VOLALT, ANTSO4, HOORI
    allocate( SULDIO_DATA   (IA,JA,22,12) ) ! 14 -> 2000-2013
    allocate(  BBSO2_DATA   (IA,JA,22,12) ) ! 14 -> 2000-2013
    allocate( VOLSO2_DATA   (IA,JA,    1) ) !
    allocate( VOLALT_DATA   (IA,JA,    1) ) !
    allocate( ANTSO4_DATA   (IA,JA,    1) ) !
    allocate(  HOORI_DATA(KA,IA,JA,   12) ) !
    SULDIO_DATA  (:,:,:,:) = 0.0_RP
     BBSO2_DATA  (:,:,:,:) = 0.0_RP
    VOLSO2_DATA  (:,:,  :) = 0.0_RP
    VOLALT_DATA  (:,:,  :) = 0.0_RP
    ANTSO4_DATA  (:,:,  :) = 0.0_RP
     HOORI_DATA(:,:,:,  :) = 0.0_RP
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( SULDIO_DATA  (:,:,:,:), fname_sulf, 'ANTSO2', LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  BBSO2_DATA  (:,:,:,:), fname_sulf, 'BBSO2',  LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( VOLSO2_DATA  (:,:,  :), fname_sulf, 'VOLSO2', LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( VOLALT_DATA  (:,:,  :), fname_sulf, 'VOLALT', LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         ( ANTSO4_DATA  (:,:,  :), fname_sulf, 'ANTSO4', LON_REAL, LAT_REAL ) ! [OUT,IN,...,IN]
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA &
         (  HOORI_DATA(:,:,:  ,:), fname_chem,    'H2O2',   LON_REAL, LAT_REAL, GDZ, GDZM ) ! [OUT,IN,...,IN]

    !----------------------
    ! Open or create Emission file for scale-rm
    !-----------------
    basename_out_mod = trim( basename_emis )
    call FILE_CARTESC_create( basename_out_mod, type_emission, dtype_emission, fid )

    !----------------------
    ! Define variables for emission inventory
    !----------------------
    totalid = 0
    ! for GRLAI (XYT)
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'GRLAI'
    write(longname,'(A)') 'Leaf Area Index'
    write(unit,    '(A)') 'm2/s2'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), & ! [IN]
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]
    
    ! for GRSNW 31 set of XYT
    do iq = 1, 31
       totalid = totalid + 1
       write(varname(totalid), '(A,I0)') 'GRSNW_', iq
       write(longname,'(A,I0,A)') 'Snow Depth on day', iq, ' in each month'
       write(unit,    '(A)') 'm'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), & ! [IN]
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo
    
    ! for OHRAD (ZXYT)
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'OHRAD'
    write(longname,'(A)') 'OH Radical'
    write(unit,    '(A)') 'molecules/cm3'
    write(dimtype(totalid), '(A)') 'ZXYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]
    
    ! for NO3G (ZXYT)
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'NO3G'
    write(longname,'(A)') 'NO3 gas'
    write(unit,    '(A)') 'ppt'
    write(dimtype(totalid), '(A)') 'ZXYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for OZONE (ZXYT)
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'O3ORI'
    write(longname,'(A)') 'O3 gas'
    write(unit,    '(A)') 'ppm'
    write(dimtype(totalid), '(A)') 'ZXYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for DAYTIME data (XYT)
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'DAYTIM'
    write(longname,'(A)') 'Day Time'
    write(unit,    '(A)') '1'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for BBBC 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'BBBC_', iq
       write(longname,'(A,i0)') 'BC emission by Biomass Burning on ', 1998+iq
       write(unit,    '(A)') 'kg(BC)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo

    ! for BBOC 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'BBOC_', iq
       write(longname,'(A,i0)') 'OC emission by Biomass Burning on ', 1998+iq
       write(unit,    '(A)') 'kg(OC)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo

    ! for ANTBC 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'ANTBC_', iq
       write(longname,'(A,i0)') 'BC emission by fossil fuel on ', 1998+iq
       write(unit,    '(A)') 'kg(BC)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo

    ! for ANTOC 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'ANTOC_', iq
       write(longname,'(A,i0)') 'OC emission by fossil fuel on ', 1998+iq
       write(unit,    '(A)') 'kg(OC)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo

    ! for ISOP of XYT
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'ISOP'
    write(longname,'(A)') 'Isoplane'
    write(unit,    '(A)') 'kg(C)/m2/s'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for TERP of XYT
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'TERP'
    write(longname,'(A)') 'Terpene'
    write(unit,    '(A)') 'kg(C)/m2/s'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for TERP of XYT
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'CHLOA'
    write(longname,'(A)') 'Chlorophyll-A (CHL-A) concentration'
    write(unit,    '(A)') 'mg/m3'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for DIFATT of XYT
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'DIFATT'
    write(longname,'(A)') 'Diffuse attenuation coef. at 490nm'
    write(unit,    '(A)') '1/m'
    write(dimtype(totalid), '(A)') 'XYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]

    ! for ANTSO2 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'ANTSO2_', iq
       write(longname,'(A,i0)') 'SO2 emission by fossil fuel on ', 1998+iq
       write(unit,    '(A)') 'kg(SO2)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo
   
    ! for BBSO2 22 set of XYT
    do iq = 1, 22
       totalid = totalid + 1
       write(varname(totalid), '(A,i0)') 'BBSO2_', iq
       write(longname,'(A,i0)') 'SO2 emission by Biomass Burning on ', 1998+iq
       write(unit,    '(A)') 'kg(SO2)/m2/s'
       write(dimtype(totalid), '(A)') 'XYT'
       call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                                  trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                                  vid(totalid),                                            & ! [OUT]
                                  "",                                                      & ! [IN:optional]
                                  time_interval, outnstep                                  ) ! [IN:optional]
    enddo

    ! for VOLSO2 of XY
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'VOLSO2'
    write(longname,'(A)') 'SO2 emission from volcano'
    write(unit,    '(A)') 'kg(SO2)/m2/s'
    write(dimtype(totalid), '(A)') 'XY'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid)                                             ) ! [OUT]

    ! for VOLALT of XY
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'VOLALT'
    write(longname,'(A)') 'Altitude of volcano'
    write(unit,    '(A)') 'm'
    write(dimtype(totalid), '(A)') 'XY'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid)                                             ) ! [OUT]

    ! for ANTSO4 of XY
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'ANTSO4'
    write(longname,'(A)') 'SO4 emission from fossil fuel'
    write(unit,    '(A)') 'kg(SO4)/m2/s'
    write(dimtype(totalid), '(A)') 'XY'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid)                                             ) ! [OUT]

    ! for ANTSO4 of XY
    totalid = totalid + 1
    write(varname(totalid), '(A)') 'H2O2'
    write(longname,'(A)') 'H2O2 gas'
    write(unit,    '(A)') 'ppm'
    write(dimtype(totalid), '(A)') 'ZXYT'
    call FILE_CARTESC_def_var( fid, trim(varname(totalid)), trim(longname), trim(unit), &
                               trim(dimtype(totalid)), dtype_emission,                  & ! [IN]
                               vid(totalid),                                            & ! [OUT]
                               "",                                                      & ! [IN:optional]
                               time_interval, outnstep                                  ) ! [IN:optional]


    do iq = 1, totalid
       LOG_INFO("ATMOS_PHY_AE_SPRINTARS_mkemission",'(I5,1X,A,1X,A,1X,I5)') &
                iq, trim(varname(iq)), trim(dimtype(iq)), vid(iq)
    enddo

    !----------------------
    ! End definition of the variables 
    !----------------------
    call FILE_CARTESC_enddef( fid )

    !----------------------
    ! Write to emission file
    !----------------------
    ! GRLAI
    totalid = 1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 GRLAI_DATA(:,:,:),       & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]

    ! 31 set of GRSNW  
    do iq = 1, 31
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                    & ! [IN]
                                    vid(totalid),           & ! [IN]
                                    GRSNW_DATA(:,:,:,iq),   & ! [IN]
                                    trim(varname(totalid)), & ! [IN]
                                    trim(dimtype(totalid)), & ! [IN]
                                    time_interval           ) ! [IN]
    enddo
               
    ! OHRAD                  
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 OHRAD_DATA(:,:,:,:),     & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! NO3G               
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 NO3G_DATA(:,:,:,:),      & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! O3ORI
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 O3ORI_DATA(:,:,:,:),     & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! DAYTIM
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 DAYTIM_DATA(:,:,:),      & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! BBBC
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    BBBC_DATA(:,:,iq,:),     & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! BBOC
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    BBOC_DATA(:,:,iq,:),     & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! ANTBC
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    ANTBC_DATA(:,:,iq,:),    & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! ANTOC
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    ANTOC_DATA(:,:,iq,:),    & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! ISOP
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 ISOP_DATA(:,:,:),        & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    ! TERP
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 TERP_DATA(:,:,:),        & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! CHLOA
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 CHLOA_DATA(:,:,:),       & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! DIFATT
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 DIFATT_DATA(:,:,:),      & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    ! ANTSO2
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    SULDIO_DATA(:,:,iq,:),   & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! BBSO2
    do iq = 1, 22
       totalid = totalid+1
       call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                    vid(totalid),            & ! [IN]
                                    BBSO2_DATA(:,:,iq,:),    & ! [IN]
                                    trim(varname(totalid)),  & ! [IN]
                                    trim(dimtype(totalid)),  & ! [IN]
                                    time_interval            ) ! [IN]
    enddo
       
    ! VOLSO2
    work2d(:,:) = VOLSO2_DATA(:,:,1)
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 work2d(:,:),             & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid))   ) ! [IN]
    
    ! VOLALT
    work2d(:,:) = VOLALT_DATA(:,:,1)
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 work2d(:,:),             & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid))   ) ! [IN]
    
    ! ANTSO4
    work2d(:,:) = ANTSO4_DATA(:,:,1)
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 work2d(:,:),             & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid))   ) ! [IN]
    
    ! HOORI
    totalid = totalid+1
    call FILE_CARTESC_write_var( fid,                     & ! [IN]
                                 vid(totalid),            & ! [IN]
                                 HOORI_DATA(:,:,:,:),     & ! [IN]
                                 trim(varname(totalid)),  & ! [IN]
                                 trim(dimtype(totalid)),  & ! [IN]
                                 time_interval            ) ! [IN]
    
    !----------------------
    ! End definition of the variables 
    !----------------------
    call FILE_CARTESC_close( fid )

    deallocate( GRLAI_DATA )
    deallocate( GRSNW_DATA )
    deallocate( OHRAD_DATA )
    deallocate( NO3G_DATA )
    deallocate( O3ORI_DATA )
    deallocate( DAYTIM_DATA )
    deallocate( BBBC_DATA )
    deallocate( BBOC_DATA )
    deallocate( ANTBC_DATA )
    deallocate( ANTOC_DATA )
    deallocate( ISOP_DATA )
    deallocate( TERP_DATA )
    deallocate( CHLOA_DATA )
    deallocate( DIFATT_DATA )
    deallocate( SULDIO_DATA )
    deallocate( BBSO2_DATA )
    deallocate( VOLSO2_DATA )
    deallocate( VOLALT_DATA )
    deallocate( ANTSO4_DATA )
    deallocate( HOORI_DATA )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_mkemission
  !-----------------------------------------------------------------------------

  
  !-----------------------------------------------------------------------------
  !> Adjustment / negative fixer
  subroutine ATMOS_PHY_AE_SPRINTARS_negative_fixer( QA_AE, & ! [IN]
                                                    QTRC )   ! [INOUT]
    implicit none
    
    integer,  intent(in)    :: QA_AE
    real(RP), intent(inout) :: QTRC(KA,IA,JA,QA_AE)
    integer :: k, i, j, iq
    !----------------------------------------------------
    
    !Negative fixer
    do iq = 1, QA_AE, 1
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             if ( QTRC(k,i,j,iq) < 0.0_RP ) QTRC(k,i,j,iq) = 0.0_RP
          enddo
       enddo
       enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_negative_fixer
  !-----------------------------------------------------------------------------

  
  
  !-----------------------------------------------------------------------------
  !> Aerosol Microphysics 
  subroutine ATMOS_PHY_AE_SPRINTARS_tendency( QA_AE,          & ! [IN]
                                              PRES,           & ! [IN] 3D
                                              TEMP,           & ! [IN]
                                              DENS,           & ! [IN]
                                              W,              & ! [IN]
                                              TKE_BL,         & ! [IN]
                                              QDRY,           & ! [IN]
                                              QV,             & ! [IN]
                                              QC,             & ! [IN]
                                              QI,             & ! [IN]
                                              THETA,          & ! [IN]
                                              RAIN_MP_Re,     & ! [IN]
                                              SNOW_MP_Re,     & ! [IN]
                                              VTERM_MP_RAIN,  & ! [IN]
                                              VTERM_MP_SNOW,  & ! [IN]
                                              RHOQ_t_MP_RAIN, & ! [IN]
                                              RHOQ_t_MP_SNOW, & ! [IN]
                                              RHOQ_t_CP_PREC, & ! [IN]
                                              flux_qr_cp,     & ! [IN]
                                              flux_qs_cp,     & ! [IN]
                                              cldfrac_cp,     & ! [IN]
                                              flux_qr_mp,     & ! [IN]
                                              flux_qs_mp,     & ! [IN]
                                              cldfrac_mp,     & ! [IN]
                                              ATMOS_PHY_MP_TYPE, & ! [IN] 
                                              NUMDEN_MP_CLDW, & ! [IN]
                                              NUMDEN_MP_RAIN, & ! [IN]
                                              NUMDEN_MP_CLDI, & ! [IN]
                                              NUMDEN_MP_SNOW, & ! [IN]
                                              dtau_s,         & ! [IN]
                                              SFC_PRES,       & ! [IN] 2D
                                              T2,             & ! [IN]
                                              U10,            & ! [IN]
                                              V10,            & ! [IN]
                                              U10L,           & ! [IN]
                                              V10L,           & ! [IN]
                                              USTAR,          & ! [IN]
                                              WSTAR,          & ! [IN]
                                              QSTAR,          & ! [IN]
                                              SFC_TEMP,       & ! [IN]
                                              LAND_WATER,     & ! [IN]
                                              LAND_PROPERTY,  & ! [IN]
                                              OCEAN_ICE_FRAC, & ! [IN]
                                              SFCFLX_SW_dn_net, & ! [IN]
                                              RCOSZ,            & ! [IN]
                                              dt_AE,          & ! [IN]
                                              QTRC,           & ! [IN]
                                              RHOQ_t_AE,      & ! [INOUT]
                                              CCN,            & ! [OUT]
                                              CN,              & ! [OUT]
                                              ATMOS_PHY_AE_Re, & ! [OUT]
                                              ATMOS_PHY_AE_Qe  ) ! [OUT]
 
    implicit none

    ! input
    integer,  intent(in)    :: QA_AE

    !3D
    real(RP), intent(in)    :: PRES          (KA,  IA,JA)
    real(RP), intent(in)    :: TEMP          (KA,  IA,JA)
    real(RP), intent(in)    :: DENS          (KA,  IA,JA)
    real(RP), intent(in)    :: W             (KA,  IA,JA)
    real(RP), intent(in)    :: TKE_BL        (KA,  IA,JA)
    real(RP), intent(in)    :: QDRY          (KA,  IA,JA)
    real(RP), intent(in)    :: QV            (KA,  IA,JA)
    real(RP), intent(in)    :: QC            (KA,  IA,JA)
    real(RP), intent(in)    :: QI            (KA,  IA,JA)
    real(RP), intent(in)    :: THETA         (KA,  IA,JA)
    real(RP), intent(in)    :: RAIN_MP_Re    (KA,  IA,JA)
    real(RP), intent(in)    :: SNOW_MP_Re    (KA,  IA,JA)
    real(RP), intent(in)    :: VTERM_MP_RAIN (KA,  IA,JA)
    real(RP), intent(in)    :: VTERM_MP_SNOW (KA,  IA,JA)
    real(RP), intent(in)    :: RHOQ_t_MP_RAIN(KA,  IA,JA)
    real(RP), intent(in)    :: RHOQ_t_MP_SNOW(KA,  IA,JA)
    real(RP), intent(in)    :: RHOQ_t_CP_PREC(KA,  IA,JA)
    real(RP), intent(in)    :: flux_qr_cp    (0:KA,IA,JA)
    real(RP), intent(in)    :: flux_qs_cp    (0:KA,IA,JA)
    real(RP), intent(in)    :: cldfrac_cp    (KA,  IA,JA)
    real(RP), intent(in)    :: flux_qr_mp    (0:KA,IA,JA)
    real(RP), intent(in)    :: flux_qs_mp    (0:KA,IA,JA)
    real(RP), intent(in)    :: cldfrac_mp    (KA,  IA,JA)
    character(len=H_SHORT), intent(in) :: ATMOS_PHY_MP_TYPE
    real(RP), intent(in)    :: NUMDEN_MP_CLDW(KA,  IA,JA)
    real(RP), intent(in)    :: NUMDEN_MP_RAIN(KA,  IA,JA)
    real(RP), intent(in)    :: NUMDEN_MP_CLDI(KA,  IA,JA)
    real(RP), intent(in)    :: NUMDEN_MP_SNOW(KA,  IA,JA)
    real(RP), intent(in)    :: dtau_s        (KA,  IA,JA)
    
    !2D
    real(RP), intent(in)    :: SFC_PRES(IA,JA)
    real(RP), intent(in)    :: T2      (IA,JA)
    real(RP), intent(in)    :: U10     (IA,JA)
    real(RP), intent(in)    :: V10     (IA,JA)
    real(RP), intent(in)    :: U10L    (IA,JA)
    real(RP), intent(in)    :: V10L    (IA,JA)
    real(RP), intent(in)    :: USTAR   (IA,JA)
    real(RP), intent(in)    :: WSTAR   (IA,JA)
    real(RP), intent(in)    :: QSTAR   (IA,JA)
    real(RP), intent(in)    :: SFC_TEMP(IA,JA)
    
    real(RP), intent(in)    :: LAND_WATER    (IA,JA) ! (LIA,LJA)
    real(RP), intent(in)    :: LAND_PROPERTY (IA,JA) ! (LIA,LJA)
    
    real(RP), intent(in)    :: OCEAN_ICE_FRAC(IA,JA) ! (OIA,OJA)
    
    real(RP), intent(in)    :: SFCFLX_SW_dn_net(IA,JA)
    real(RP), intent(in)    :: RCOSZ           (IA,JA)
    
    real(DP), intent(in)    :: dt_AE

    real(RP), intent(in)    :: QTRC     (KA,IA,JA,QA_AE)
    real(RP), intent(inout) :: RHOQ_t_AE(KA,IA,JA,QA_AE)
    real(RP), intent(out)   :: CCN      (KA,IA,JA)
    real(RP), intent(out)   :: CN       (KA,IA,JA)
    real(RP), intent(out)   :: ATMOS_PHY_AE_Re(KA,IA,JA,N_AE)
    real(RP), intent(out)   :: ATMOS_PHY_AE_Qe(KA,IA,JA,N_AE)
    !---------------
    
    integer  :: k, i, j, iq
    real(RP) :: QTRC0(KA,IA,JA,QA_AE)
    logical,  save :: OFIRST = .true.
    
    !----------------------------------------------------

    !copy QTRC -> QTRC0
    !================================================
    do iq = 1, QA_AE
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             QTRC0(k,i,j,iq) = QTRC(k,i,j,iq) ! save
          enddo
       enddo
       enddo
    enddo
    !================================================
    !end copy QTRC -> QTRC0


    !tendency
    !================================================
    call ATMOS_PHY_AE_SPRINTARS_AEROSL( QTRC0     (:,   :,:,1:QA_AE), & ! [INOUT]
                                        CCN       (:,   :,:),         & ! [OUT]
                                        CN        (:,   :,:),         & ! [OUT]
                                        ATMOS_PHY_AE_Re(:,:,:,:),     & ! [OUT]
                                        ATMOS_PHY_AE_Qe(:,:,:,:),     & ! [OUT]
                                        PRES      (:,   :,:),         & ! [IN] 3D
                                        TEMP      (:,   :,:),         & ! [IN]
                                        DENS      (:,   :,:),         & ! [IN]
                                        W         (:,   :,:),         & ! [IN]
                                        TKE_BL    (:,   :,:),         & ! [IN]
                                        QDRY      (:,   :,:),         & ! [IN]
                                        QV        (:,   :,:),         & ! [IN]
                                        QC        (:,   :,:),         & ! [IN]
                                        QI        (:,   :,:),         & ! [IN]
                                        THETA     (:,   :,:),         & ! [IN]
                                        RAIN_MP_Re(:,   :,:),         & ! [IN]
                                        SNOW_MP_Re(:,   :,:),         & ! [IN]
                                        VTERM_MP_RAIN(:,:,:),         & ! [IN]
                                        VTERM_MP_SNOW(:,:,:),         & ! [IN]
                                        RHOQ_t_MP_RAIN(:,:,:),        & ! [IN]
                                        RHOQ_t_MP_SNOW(:,:,:),        & ! [IN]
                                        RHOQ_t_CP_PREC(:,:,:),        & ! [IN]
                                        flux_qr_cp(0:KA,:,:),         & ! [IN]
                                        flux_qs_cp(0:KA,:,:),         & ! [IN]
                                        cldfrac_cp(:,   :,:),         & ! [IN]
                                        flux_qr_mp(0:KA,:,:),         & ! [IN]
                                        flux_qs_mp(0:KA,:,:),         & ! [IN]
                                        cldfrac_mp(:,   :,:),         & ! [IN]
                                        ATMOS_PHY_MP_TYPE,            & ! [IN]
                                        NUMDEN_MP_CLDW(:,:,:),        & ! [IN]
                                        NUMDEN_MP_RAIN(:,:,:),        & ! [IN]
                                        NUMDEN_MP_CLDI(:,:,:),        & ! [IN]
                                        NUMDEN_MP_SNOW(:,:,:),        & ! [IN]
                                        dtau_s    (:,   :,:),         & ! [IN]
                                        SFC_PRES       (:,:),         & ! [IN] 2D
                                        T2             (:,:),         & ! [IN]
                                        U10            (:,:),         & ! [IN]
                                        V10            (:,:),         & ! [IN]
                                        U10L           (:,:),         & ! [IN]
                                        V10L           (:,:),         & ! [IN]
                                        USTAR          (:,:),         & ! [IN]
                                        WSTAR          (:,:),         & ! [IN]
                                        QSTAR          (:,:),         & ! [IN]
                                        SFC_TEMP       (:,:),         & ! [IN]
                                        LAND_WATER     (:,:),         & ! [IN]
                                        LAND_PROPERTY  (:,:),         & ! [IN]
                                        OCEAN_ICE_FRAC (:,:),         & ! [IN]
                                        SFCFLX_SW_dn_net(:,:),        & ! [IN]
                                        RCOSZ           (:,:),        & ! [IN]
                                        ATMOS_PHY_AE_SPRINTARS_NAME(1:QA_AE), & ! [IN]
                                        QA_AE,                        & ! [IN] 
                                        dt_AE                         ) ! [IN]

#if defined(DEBUG) || defined(QUICKDEBUG)
    !check tracer
    call STATISTICS_detail( KA, KS, KE, IA, IS, IE, JA, JS, JE, QA_AE, &
                            ATMOS_PHY_AE_SPRINTARS_NAME(1:QA_AE), QTRC0(:,:,:,1:QA_AE) )
#endif

    do iq = 1, QA_AE
       do j = JS, JE
       do i = IS, IE
          do k = KS, KE
             RHOQ_t_AE(k,i,j,iq) = ( QTRC0(k,i,j,iq) - QTRC(k,i,j,iq) ) * DENS(k,i,j) / real(dt_AE,kind=RP)
          enddo
       enddo
       enddo
    enddo
    !================================================
    !end tendency

    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_tendency
  !-----------------------------------------------------------------------------


  
end module scale_atmos_phy_ae_SPRINTARS
