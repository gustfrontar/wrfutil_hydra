!-------------------------------------------------------------------------------
!> module atmosphere / physics / sprintars - common
!!
!! @par Description
!!          sprintars aerosol microphysics scheme
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_ae_SPRINTARS_common
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
       PI   => CONST_PI,   & !< 
       GRAV => CONST_GRAV, & !< standard acceleration of gravity [m/s2] = 9.80665_RP
       CONST_D2R, &          !< degree to radian
       CONST_UNDEF
  
  ! use scale_comm_cartesC, only: &
  !      COMM_bcast
  
  use scale_interp, only: &
       INTERP_factor2d, &
       INTERP_interp2d, &
       INTERP_domain_compatibility, &
       INTERP_factor3d,             &
       INTERP_interp3d

  use scale_file, only: &
       FILE_open,      &
       FILE_get_shape, &
       FILE_read,      &
       FILE_close

  use scale_atmos_grid_cartesC, only: &
       CX => ATMOS_GRID_CARTESC_CX, &
       CY => ATMOS_GRID_CARTESC_CY
  use scale_mapprojection, only: &
       MAPPROJECTION_lonlat2xy
  use scale_topography, only: &
       GDZS => TOPOGRAPHY_Zsfc
  
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA !netcdf
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY !netcdf
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST_z
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z
  public :: ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC

  interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_1D
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_3D
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_4D
  end interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA

  interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT
  end interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA

  interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D
     module procedure ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D
  end interface ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,  public, parameter :: NSU = 3 !< (3) sulfate   tracer No.
  integer,  public, parameter :: NCA = 6 !< (6) carbon    tracer No.
  integer,  public, parameter :: NDU = 6 !< (6) dust      tracer No.
  integer,  public, parameter :: NSA = 4 !< (4) salt      tracer No.
  integer,  public, parameter :: NAT = 1 !<     arbitrary tracer No.

  integer,  public :: QA_AE_SU = 0, QS_AE_SU = -1, QE_AE_SU = -2 !< QA_AE_SU = NSU, QTRC(:,:,:,QS_AE_SU:QE_AE_SU)
  integer,  public :: QA_AE_CA = 0, QS_AE_CA = -1, QE_AE_CA = -2 !< QA_AE_CA = NCA, QTRC(:,:,:,QS_AE_CA:QE_AE_CA)
  integer,  public :: QA_AE_DU = 0, QS_AE_DU = -1, QE_AE_DU = -2 !< QA_AE_DU = NDU, QTRC(:,:,:,QS_AE_DU:QE_AE_DU)
  integer,  public :: QA_AE_SA = 0, QS_AE_SA = -1, QE_AE_SA = -2 !< QA_AE_SA = NSA, QTRC(:,:,:,QS_AE_SA:QE_AE_SA)
  integer,  public :: QA_AE_AT = 0, QS_AE_AT = -1, QE_AE_AT = -2 !< QA_AE_AT = NAT, QTRC(:,:,:,QS_AE_AT:QE_AE_AT)
  
  integer,  public, parameter :: IWA = 3                         !< number of wavelengths
  integer,  public, parameter :: NRH = 8                         !< number of RH bin
  real(RP), public, parameter :: &
       BRH(NRH) &                                                !< RH bin
       = (/ 0.00E+0_RP, 0.50E+0_RP, 0.70E+0_RP, 0.80E+0_RP, &
            0.90E+0_RP, 0.95E+0_RP, 0.98E+0_RP, 0.99E+0_RP /)

  !real(RP), public, parameter :: VMISS = -9.9999E+30_RP     !< missing value
  real(RP), public, parameter :: NU       =  1.78E-5_RP     !< air viscosity    [kg/m/s] = [Pa*s]
  real(RP), public, parameter :: RADR     =  5.00E-4_RP     !< rain droplet radius         [m]
  real(RP), public, parameter :: VTR      =  0.50E+0_RP     !< rain terminal velocity (CP) [m/s]
  real(RP), public, parameter :: DENSR    =  1.00E+3_RP     !< rain droplet density        [kg/m3]
  real(RP), public, parameter :: DENSW    =  1.00E+3_RP     !< density of water            [kg/m3]
  real(RP), public, parameter :: DENSI    =  9.168E+2_RP    !< density of ice              [kg/m3]
  real(RP), public, parameter :: AVOG     =  6.02214E+23_RP !< Avogadro's number
  real(RP), public, parameter :: THRESRe  =  1.00E-20_RP    !< threshold of rain droplet radius [m] 
  real(RP), public, parameter :: THRESP   =  1.00E-20_RP    !< precipitation threshold [kg/m2/s] (e.g., 0.01 [mm/h] = 2.78E-6 [kg/m2/s] )
  real(RP), public, parameter :: THRESVTR =  1.00E-20_RP    !< terminal velocity threshold for rain(MP) & snow(MP) [m/s]

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
  !> SPRINTARS SET_DAILY (set daily data)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY &
       ( D_DATA,                   &  ! [OUT]
         MD_DATA, M_DATA, MJ_DATA, &  ! [IN]
         nyear, nmonth, nday       )  ! [IN]
    
    implicit none

    real(RP), intent(out) ::  D_DATA(IA,JA)
    real(RP), intent(in)  :: MD_DATA(IA,JA)
    real(RP), intent(in)  ::  M_DATA(IA,JA,12)
    real(RP), intent(in)  :: MJ_DATA(IA,JA)
    integer,  intent(in)  :: nyear, nmonth, nday

    integer :: ndyear
    integer :: smonth, emonth
    integer :: aday, bday
    integer, parameter :: mday(12,2) &
         !              1   2   3   4   5   6   7   8   9  10  11  12
         = reshape( (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31, &
                       31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31  /), &
                    (/ 12, 2 /) )
    integer :: i, j
    !----------------------------------------------------

    ndyear = 1
    if ( mod(nyear,4) == 0 ) then
       ndyear = 2
       if ( mod(nyear,100) == 0 .and. mod(nyear,400) /= 0 ) then
          ndyear = 1
       endif
    endif
    
    if ( nday == 15 ) then
       do j = JS, JE
       do i = IS, IE
          D_DATA(i,j) = M_DATA(i,j,nmonth)
       enddo
       enddo
    
    elseif ( nday < 15 ) then
       smonth = nmonth - 1
       if ( smonth == 0 ) then
          smonth = 12
          aday = mday(smonth,ndyear) - 15 + nday
          bday = mday(smonth,ndyear) - 15 +   15
          do j = JS, JE
          do i = IS, IE
             D_DATA(i,j) = ( M_DATA(i,j,nmonth) - MD_DATA(i,j) ) / bday * aday + MD_DATA(i,j)
          enddo
          enddo
       else
          aday = mday(smonth,ndyear) - 15 + nday
          bday = mday(smonth,ndyear) - 15 +   15
          do j = JS, JE
          do i = IS, IE
             D_DATA(i,j) = ( M_DATA(i,j,nmonth) - M_DATA(i,j,smonth) ) / bday * aday + M_DATA(i,j,smonth)
          enddo
          enddo
       endif
       
    elseif ( nday > 15 ) then
       emonth = nmonth + 1
       if ( emonth == 13 ) then
          emonth = 1
          aday = nday - 15
          bday = mday(nmonth,ndyear) - nday + 15
          do j = JS, JE
          do i = IS, IE
             D_DATA(i,j) = ( MJ_DATA(i,j) - M_DATA(i,j,nmonth) ) / bday * aday + M_DATA(i,j,nmonth)
          enddo
          enddo
       else
          aday = nday - 15
          bday = mday(nmonth,ndyear) - nday + 15
          do j = JS, JE
          do i = IS, IE
             D_DATA(i,j) = ( M_DATA(i,j,emonth) - M_DATA(i,j,nmonth) ) / bday * aday + M_DATA(i,j,nmonth)
          enddo
          enddo
       endif
       
    endif

    return

  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_SET_DAILY
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_NETCDF_DATA (read netCDF data; 1D )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_1D &
       ( DATA, fname, varname ) ! [OUT, IN, IN]
    implicit none
    real(RP),     allocatable, intent(out) :: DATA(:)
    character(*),              intent(in)  :: fname
    character(*),              intent(in)  :: varname
    integer :: fid
    integer :: dims(1)
    real(RP), allocatable :: tmp1d(:)
    logical :: exist
    integer :: m
    !----------------------------------------------------
    inquire( file = trim(adjustl(fname)), exist=exist )
    if ( .not. exist ) then !--- missing
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_1D",*) 'File not found. check! file name :', fname
       call PRC_abort
    endif

    call FILE_open( fname,                & ! (in)
                    fid,                  & ! (out)
                    aggregate = .false.,  & ! (in)
                    postfix   = ""        ) ! (in)
    call FILE_get_shape( fid, varname, dims(:) )
    allocate( DATA (dims(1)) )
    allocate( tmp1d(dims(1)) )
    call FILE_read( fid, varname, tmp1d(:) )
    do m = 1, dims(1)
       DATA(m) = real( tmp1d(m), kind=RP )
    enddo
    deallocate( tmp1d )
    call FILE_close( fid )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_1D
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_NETCDF_DATA (read netCDF data; 3D )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_3D &
       ( DATA, fname, varname ) ! [OUT, IN, IN]
    implicit none
    real(RP),     allocatable, intent(out) :: DATA(:,:,:)
    character(*),              intent(in)  :: fname
    character(*),              intent(in)  :: varname
    integer :: fid
    integer :: dims(3)
    real(RP), allocatable :: tmp3d(:,:,:)
    logical :: exist
    integer :: m, n, t
    !----------------------------------------------------
    inquire( file = trim(adjustl(fname)), exist=exist )
    if ( .not. exist ) then !--- missing
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_3D",*) 'File not found. check! file name :', fname
       call PRC_abort
    endif

    call FILE_open( fname,                & ! (in)
                    fid,                  & ! (out)
                    aggregate = .false.,  & ! (in)
                    postfix   = ""        ) ! (in)
    call FILE_get_shape( fid, varname, dims(:) )
    allocate( DATA (dims(1),dims(2),dims(3)) )
    allocate( tmp3d(dims(1),dims(2),dims(3)) )
    call FILE_read( fid, varname, tmp3d(:,:,:) )
    do t = 1, dims(3)
       do n = 1, dims(2)
          do m = 1, dims(1)
             DATA(m,n,t) = real( tmp3d(m,n,t), kind=RP )
          enddo
       enddo
    enddo
    deallocate( tmp3d )
    call FILE_close( fid )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_3D
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_NETCDF_DATA (read netCDF data; 4D )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_4D &
       ( DATA, fname, varname ) ! [OUT, IN, IN]
    implicit none
    real(RP),     allocatable, intent(out) :: DATA(:,:,:,:)
    character(*),              intent(in)  :: fname
    character(*),              intent(in)  :: varname
    integer :: fid
    integer :: dims(4)
    real(RP), allocatable :: tmp4d(:,:,:,:)
    logical :: exist
    integer :: l, m, n, t
    !----------------------------------------------------
    inquire( file = trim(adjustl(fname)), exist=exist )
    if ( .not. exist ) then !--- missing
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_4D",*) 'File not found. check! file name :', fname
       call PRC_abort
    endif

    call FILE_open( fname,                & ! (in)
                    fid,                  & ! (out)
                    aggregate = .false.,  & ! (in)
                    postfix   = ""        ) ! (in)
    call FILE_get_shape( fid, varname, dims(:) )
    allocate( DATA (dims(1),dims(2),dims(3),dims(4)) )
    allocate( tmp4d(dims(1),dims(2),dims(3),dims(4)) )
    call FILE_read( fid, varname, tmp4d(:,:,:,:) )
    do t = 1, dims(4)
       do n = 1, dims(3)
          do m = 1, dims(2)
             do l = 1, dims(1)
                DATA(l,m,n,t) = real( tmp4d(l,m,n,t), kind=RP )
             enddo
          enddo
       enddo
    enddo
    deallocate( tmp4d )
    call FILE_close( fid )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA_4D
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_DATA ( read data (3D): DATA(IA,JA,month(year)) )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT &
       ( DATA,              & ! [OUT]
         fname, varname,    & ! [IN]
         LON_REAL, LAT_REAL ) ! [IN]
    implicit none
    
    real(RP),     intent(out) :: DATA    (:,:,:) !< data(IA,JA,month)
    character(*), intent(in)  :: fname           !< data file name
    character(*), intent(in)  :: varname         !< variable name
    real(RP),     intent(in)  :: LON_REAL(IA,JA) !< longitude [rad]
    real(RP),     intent(in)  :: LAT_REAL(IA,JA) !< latitude  [rad]

    integer,  parameter   :: nintrp = 4
    real(RP), allocatable :: lon(:,:)
    real(RP), allocatable :: lat(:,:)
    real(RP), allocatable :: data1d_lon(:)
    real(RP), allocatable :: data1d_lat(:)
    real(RP), allocatable :: data3d(:,:,:)
    integer :: nlon, nlat, nmonth
    integer :: ilon, ilat, imonth

    integer  :: idx_i(IA,JA,nintrp)
    integer  :: idx_j(IA,JA,nintrp)
    real(RP) :: hfact(IA,JA,nintrp)
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//varname
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data3d, fname, varname )
    if ( size(data3d,3) /= size(DATA,3) ) then
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'array size is wrong. check! file name:', fname//': '//varname
       call PRC_abort
    endif
    nlon   = size(data3d,1)
    nlat   = size(data3d,2)
    nmonth = size(data3d,3)
!    if ( varname == 'GRLAI' ) then
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
!    else
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon'
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat'
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon' )
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat' )
!    endif
    if ( varname == 'ANTSO4' ) then
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_ANSO4T'
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_ANSO4T'
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_ANSO4T' )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_ANSO4T' )
    else
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
    endif
    allocate( lon(nlon,nlat), lat(nlon,nlat) )
    do ilat = 1, nlat
       do ilon = 1, nlon
          lon(ilon,ilat) = data1d_lon(ilon) * CONST_D2R
          lat(ilon,ilat) = data1d_lat(ilat) * CONST_D2R
       enddo
    enddo
    
    call INTERP_factor2d( nintrp,        & ! [IN]
                          nlon, nlat,    & ! [IN]
                          IA, JA,        & ! [IN]
                          lon(:,:),      & ! [IN]
                          lat(:,:),      & ! [IN]
                          LON_REAL(:,:), & ! [IN]
                          LAT_REAL(:,:), & ! [IN]
                          idx_i(:,:,:),  & ! [OUT]
                          idx_j(:,:,:),  & ! [OUT]
                          hfact(:,:,:)   ) ! [OUT]

    do imonth = 1, nmonth
       call INTERP_interp2d( nintrp,             & ! [IN]
                             nlon, nlat,         & ! [IN]
                             IA, JA,             & ! [IN]
                             idx_i(:,:,:),       & ! [IN]
                             idx_j(:,:,:),       & ! [IN]
                             hfact(:,:,:),       & ! [IN]
                             data3d(:,:,imonth), & ! [IN]
                             DATA  (:,:,imonth)  ) ! [OUT]
    enddo

    deallocate( data1d_lon, data1d_lat, data3d, lon, lat )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_3D_XYT
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_DATA ( read data (4D) : DATA(IA,JA,*,*) )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT &
       ( DATA,              & ! [OUT]
         fname, varname,    & ! [IN]
         LON_REAL, LAT_REAL ) ! [IN]
    implicit none
    
    real(RP),     intent(out) :: DATA    (:,:,:,:) !< data(IA,JA,*,*)
    character(*), intent(in)  :: fname             !< data file name
    character(*), intent(in)  :: varname           !< variable name
    
    real(RP),     intent(in)  :: LON_REAL(IA,JA)   !< longitude [rad]
    real(RP),     intent(in)  :: LAT_REAL(IA,JA)   !< latitude  [rad]

    integer,  parameter   :: nintrp = 4
    real(RP), allocatable :: lon(:,:)
    real(RP), allocatable :: lat(:,:)
    real(RP), allocatable :: data1d_lon(:)
    real(RP), allocatable :: data1d_lat(:)
    real(RP), allocatable :: data4d(:,:,:,:)
    integer :: nlon, nlat, nmonth, nday
    integer :: ilon, ilat, imonth, iday

    integer  :: idx_i(IA,JA,nintrp)
    integer  :: idx_j(IA,JA,nintrp)
    real(RP) :: hfact(IA,JA,nintrp)
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//varname
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d, fname, varname )
    if ( size(data4d,3) /= size(DATA,3) .or. size(data4d,4) /= size(DATA,4) ) then
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'array size is wrong. check! file name :', fname//': '//varname
       call PRC_abort
    endif
    
    nlon   = size(data4d,1)
    nlat   = size(data4d,2)
    nmonth = size(data4d,3)
    nday   = size(data4d,4)
!    if ( varname == 'GRSNW' ) then
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
!    else
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lon'
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lat'
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon' )
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat' )
!    endif
     if( varname == 'ANTSO2' ) then
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_SULDIO'
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_SULDIO'
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_SULDIO' )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_SULDIO' )
     else
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
    endif

    allocate( lon(nlon,nlat), lat(nlon,nlat) )
    do ilat = 1, nlat
       do ilon = 1, nlon
          lon(ilon,ilat) = data1d_lon(ilon) * CONST_D2R
          lat(ilon,ilat) = data1d_lat(ilat) * CONST_D2R
       enddo
    enddo
 
    call INTERP_factor2d( nintrp,        & ! [IN]
                          nlon, nlat,    & ! [IN]
                          IA, JA,        & ! [IN]
                          lon(:,:),      & ! [IN]
                          lat(:,:),      & ! [IN]
                          LON_REAL(:,:), & ! [IN]
                          LAT_REAL(:,:), & ! [IN]
                          idx_i(:,:,:),  & ! [OUT]
                          idx_j(:,:,:),  & ! [OUT]
                          hfact(:,:,:)   ) ! [OUT]

    do imonth = 1, nmonth
       do iday = 1, nday
          call INTERP_interp2d( nintrp,                  & ! [IN]
                                nlon, nlat,              & ! [IN]
                                IA, JA,                  & ! [IN]
                                idx_i(:,:,:),            & ! [IN]
                                idx_j(:,:,:),            & ! [IN]
                                hfact(:,:,:),            & ! [IN]
                                data4d(:,:,imonth,iday), & ! [IN]
                                DATA  (:,:,imonth,iday)  ) ! [OUT]
       enddo
    enddo
    
    deallocate( data1d_lon, data1d_lat, data4d, lon, lat )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_XYTT
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_DATA (read data (4D) : DATA(KA,IA,JA,*) )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT &
       ( DATA,               & ! [OUT]
         fname, varname,     & ! [IN]
         LON_REAL, LAT_REAL, & ! [IN]
         GDZ, GDZM           ) ! [IN]
    implicit none

    real(RP),     intent(out) :: DATA(:,:,:,:)        !< data(KA,IA,JA,year)
    character(*), intent(in)  :: fname                !< data file name
    character(*), intent(in)  :: varname              !< variable name
    real(RP),     intent(in)  :: LON_REAL     (IA,JA) !< longitude [rad]
    real(RP),     intent(in)  :: LAT_REAL     (IA,JA) !< latitude  [rad]
    real(RP),     intent(in)  :: GDZ     (KA,  IA,JA) !< CZ
    real(RP),     intent(in)  :: GDZM    (0:KA,IA,JA) !< FZ
    
    integer,  parameter       :: itp_nh_a = 4
    integer,  parameter       :: nintrp = 4
    real(RP), allocatable     :: lon(:,:)
    real(RP), allocatable     :: lat(:,:)
    real(RP), allocatable     :: z(:,:,:,:)
    real(RP), allocatable     :: data1d_lon(:)
    real(RP), allocatable     :: data1d_lat(:)
    real(RP), allocatable     :: data4d(:,:,:,:)
    real(RP), allocatable     :: data4d2(:,:,:,:)
    real(RP), allocatable     :: data4d_zg(:,:,:,:)
    
    integer :: nlon, nlat, nz, nmonth
    integer :: ilon, ilat, iz, imonth

    integer  :: igrd         (IA,JA,itp_nh_a)
    integer  :: jgrd         (IA,JA,itp_nh_a)
    real(RP) :: hfact        (IA,JA,itp_nh_a)
    integer  :: kgrd    (KA,2,IA,JA,itp_nh_a)
    real(RP) :: vfact   (KA,  IA,JA,itp_nh_a)

    real(RP), allocatable :: X_org(:,:)
    real(RP), allocatable :: Y_org(:,:)
    logical :: zonal, pole
    real(RP) :: wsum(KA,IA,JA)
    real(RP) :: work(KA,IA,JA)
    integer  :: kref

    integer  :: k, i, j
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//varname
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d2, fname, varname )
    if ( size(data4d2,4) /= size(DATA,4) ) then
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'array size is wrong. check! file name:', fname//': '//varname
       call PRC_abort
    endif
    nz     = size(data4d2,1)
    nlon   = size(data4d2,2)
    nlat   = size(data4d2,3)
    nmonth = size(data4d2,4)
    !LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon'
    !LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat'
    !call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon' )
    !call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat' )
    if( varname == 'OZONE' ) then
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_O3ORI'
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_O3ORI'
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'zg'
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_O3ORI' )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_O3ORI' )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d_zg,  fname, 'zg'  )
    else
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'READ DATA, FILE NAME:', fname//': '//'zg'
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d_zg,  fname, 'zg'  )
    endif
    
    allocate( X_org(nlon,nlat), Y_org(nlon,nlat) )
    allocate( data4d(nz+2,nlon,nlat,nmonth) )
    allocate( lon(nlon,nlat), lat(nlon,nlat), z(nz+2,nlon,nlat,nmonth) )
    data4d(:,:,:,:) = CONST_UNDEF
    do ilat = 1, nlat
       do ilon = 1, nlon
          lon(ilon,ilat) = data1d_lon(ilon) * CONST_D2R
          lat(ilon,ilat) = data1d_lat(ilat) * CONST_D2R
       enddo
    enddo
    do iz = 1, nz
       z     (iz+2,:,:,:) = data4d_zg(iz,:,:,:)
       data4d(iz+2,:,:,:) = data4d2  (iz,:,:,:)
    enddo
    
    do imonth = 1, nmonth
       iz = nz+2
       call INTERP_domain_compatibility( lon(:,:),           & ! [IN]
                                         lat(:,:),           & ! [IN]
                                         z  (iz,:,:,imonth), & ! [IN]
                                         LON_REAL(:,:),      & ! [IN]
                                         LAT_REAL(:,:),      & ! [IN]
                                         GDZ  (KE,:,:),      & ! [IN]
                                         GDZM (KE,:,:)       ) ! [IN]
    enddo

    if ( nlon == 1 .or. nlat == 1 ) then
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'LINER interpolation requires nx, ny > 1'
       LOG_ERROR_CONT(*)               'Use "DIST-WEIGHT" as INTRP_TYPE of PARAM_MKINIT_REAL_ATMOS'
       call PRC_abort
    endif

    !$omp parallel do collapse(2)
    do ilat = 1, nlat
    do ilon = 1, nlon
       lat(ilon,ilat) = sign( min( abs(lat(ilon,ilat)), PI * 0.499999_RP ), lat(ilon,ilat) )
    enddo
    enddo

    call MAPPROJECTION_lonlat2xy( nlon, 1, nlon, & ! [IN]
                                  nlat, 1, nlat, & ! [IN]
                                  lon(:,:),      & ! [IN]
                                  lat(:,:),      & ! [IN]
                                  X_org  (:,:),  & ! [OUT]
                                  Y_org  (:,:)   ) ! [OUT]
    
    zonal = ( maxval(lon) - minval(lon) ) > 2.0_RP * PI * 0.9_RP
    pole = ( maxval(lat) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(lat) < - PI * 0.5_RP * 0.9_RP )

    do imonth = 1, nmonth
       
       call INTERP_factor3d( nz+2, 1, nz+2,          & ! [IN]
                             nlon, nlat,             & ! [IN]
                             KA, KS, KE,             & ! [IN]
                             IA, JA,                 & ! [IN]
                             X_org(:,:), Y_org(:,:), & ! [IN]
                             z (:,:,:,imonth),       & ! [IN]
                             CX(:), CY(:),           & ! [IN]
                             GDZ    (:,:,:),         & ! [IN]
                             igrd   (    :,:,:),     & ! [OUT]
                             jgrd   (    :,:,:),     & ! [OUT]
                             hfact  (    :,:,:),     & ! [OUT]
                             kgrd   (:,:,:,:,:),     & ! [OUT]
                             vfact  (:,  :,:,:),     & ! [OUT]
                             flag_extrap = .false.,  & ! [IN]
                             zonal = zonal,          & ! [IN]
                             pole  = pole            ) ! [IN]

       call INTERP_interp3d( itp_nh_a,                    &
                             nz+2, 1, nz+2,               &
                             nlon, nlat,                  &
                             KA, KS, KE,                  &
                             IA, JA,                      &
                             igrd(:,:,:), jgrd(:,:,:),    & ! [IN]
                             hfact(:,:,:),                & ! [IN]
                             kgrd(:,:,:,:,:),             & ! [IN]
                             vfact(:,:,:,:),              & ! [IN]
                             z(:,:,:,imonth), GDZ(:,:,:), & ! [IN]
                             data4d(:,:,:,imonth),        & ! [IN]
                             DATA  (:,:,:,imonth),        & ! [OUT]
                             threshold_undef = 1.0_RP,    & ! [IN]
                             wsum = wsum(:,:,:),          & ! [OUT]
                             val2 = work(:,:,:)           ) ! [OUT]

       do j = JS, JE
       do i = IS, IE

          kref = -1
          do k = KS, KE
             if ( DATA(k,i,j,imonth) /= CONST_UNDEF ) then
                kref = k
                exit
             endif
             if ( k == KE ) then
                LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT",*) 'kref == -1. check!'
                call PRC_abort
             endif
          enddo

          if ( kref >= KS+1 .and. kref <= KE ) then
             do k = kref-1, KS, -1
                ! DATA(k,i,j,imonth) = DATA(k+1,i,j,imonth) * log( ( GDZ(k,i,j) - GDZS(i,j) ) / 1.0E-4_RP ) / log( ( GDZ(k+1,i,j) - GDZS(i,j) ) / 1.0E-4_RP ) * ( 1.0_RP - wsum(k,i,j) ) &
                !                    + work(k,i,j) * wsum(k,i,j)
                DATA(k,i,j,imonth) = DATA(k+1,i,j,imonth)
                ! DATA(k,i,j,imonth) = 0.0_RP
             enddo
          endif
          if ( kref >= KS .and. kref <= KE-1 ) then

             do k = kref+1, KE
                if ( DATA(k,i,j,imonth) == CONST_UNDEF ) DATA(k,i,j,imonth) = DATA(k-1,i,j,imonth)
                ! if ( DATA(k,i,j,imonth) == CONST_UNDEF ) DATA(k,i,j,imonth) = 0.0_RP
             enddo
          endif
          
       enddo
       enddo

    enddo

    deallocate( lon, lat, z, data1d_lon, data1d_lat, data4d, data4d2 )
    deallocate( X_org, Y_org )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_4D_ZXYT
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_DATA_DAILY ( read data (3D): DATA(IA,JA,day) )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D &
       ( DATA,              & ! [OUT]
         fname, varname,    & ! [IN]
         LON_REAL, LAT_REAL ) ! [IN]
    implicit none

    real(RP),     allocatable, intent(out) :: DATA(:,:,:)     !< data(IA,JA,day)
    character(*),              intent(in)  :: fname           !< data file name
    character(*),              intent(in)  :: varname         !< variable name
    real(RP),                  intent(in)  :: LON_REAL(IA,JA) !< longitude [rad]
    real(RP),                  intent(in)  :: LAT_REAL(IA,JA) !< latitude  [rad]
    
    integer,  parameter   :: nintrp = 4
    real(RP), allocatable :: lon(:,:)
    real(RP), allocatable :: lat(:,:)
    real(RP), allocatable :: data1d_lon(:)
    real(RP), allocatable :: data1d_lat(:)
    real(RP), allocatable :: data3d(:,:,:)
    integer  :: nlon, nlat, nday
    integer  :: ilon, ilat, iday
    integer  :: idx_i(IA,JA,nintrp)
    integer  :: idx_j(IA,JA,nintrp)
    real(RP) :: hfact(IA,JA,nintrp)

    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D",*) 'READ DATA, FILE NAME:', fname//': '//varname
    
    !read data3d, data1d_lon, & data1d_lat
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data3d, fname, varname )
    nlon = size(data3d,1)
    nlat = size(data3d,2)
    nday = size(data3d,3)
    allocate( lon(nlon,nlat), lat(nlon,nlat), DATA(nlon,nlat,nday) )
!    if ( varname == 'GRLAI' .or. varname == 'GRSNW' ) then
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
!    else
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D",*) 'READ DATA, FILE NAME:', fname//': '//'lon'
!       LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D",*) 'READ DATA, FILE NAME:', fname//': '//'lat'
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon' )
!       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat' )
!    endif
    
    !set lon & lat
    do ilat = 1, nlat
       do ilon = 1, nlon
          lon(ilon,ilat) = data1d_lon(ilon) * CONST_D2R
          lat(ilon,ilat) = data1d_lat(ilat) * CONST_D2R
       enddo
    enddo

    !set DATA
    call INTERP_factor2d( nintrp,        & ! [IN]
                          nlon, nlat,    & ! [IN]
                          IA, JA,        & ! [IN]
                          lon(:,:),      & ! [IN]
                          lat(:,:),      & ! [IN]
                          LON_REAL(:,:), & ! [IN]
                          LAT_REAL(:,:), & ! [IN]
                          idx_i(:,:,:),  & ! [OUT]
                          idx_j(:,:,:),  & ! [OUT]
                          hfact(:,:,:)   ) ! [OUT]

    do iday = 1, nday
       call INTERP_interp2d( nintrp,           & ! [IN]
                             nlon, nlat,       & ! [IN]
                             IA, JA,           & ! [IN]
                             idx_i(:,:,:),     & ! [IN]
                             idx_j(:,:,:),     & ! [IN]
                             hfact(:,:,:),     & ! [IN]
                             data3d(:,:,iday), & ! [IN]
                             DATA  (:,:,iday)  ) ! [OUT]
    enddo

    deallocate( data1d_lon, data1d_lat, data3d, lon, lat )

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_3D
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS READ_DATA_DAILY (read data (4D) : DATA(KA,IA,JA,day) )
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D &
       ( DATA,               & ! [OUT]
         fname, varname,     & ! [IN]
         LON_REAL, LAT_REAL, & ! [IN]
         GDZ, GDZM           ) ! [IN]
    implicit none

    real(RP),     allocatable, intent(out) :: DATA(:,:,:,:)        !< data(KA,IA,JA,day)
    character(*),              intent(in)  :: fname                !< data file name
    character(*),              intent(in)  :: varname              !< variable name
    real(RP),                  intent(in)  :: LON_REAL     (IA,JA) !< longitude [rad]
    real(RP),                  intent(in)  :: LAT_REAL     (IA,JA) !< latitude  [rad]
    real(RP),                  intent(in)  :: GDZ     (KA,  IA,JA) !< CZ
    real(RP),                  intent(in)  :: GDZM    (0:KA,IA,JA) !< FZ
    
    integer,  parameter   :: itp_nh_a = 4
    integer,  parameter   :: nintrp = 4
    real(RP), allocatable :: lon(:,:)
    real(RP), allocatable :: lat(:,:)
    real(RP), allocatable :: z(:,:,:,:)
    real(RP), allocatable :: data1d_lon(:)
    real(RP), allocatable :: data1d_lat(:)
    real(RP), allocatable :: data1d_zg (:)
    real(RP), allocatable :: data4d(:,:,:,:)
    real(RP), allocatable :: data4d2(:,:,:,:)
    real(RP), allocatable :: data4d_zg(:,:,:,:)
    real(RP), allocatable :: data4d2_zg(:,:,:,:)
    
    integer  :: nlon, nlat, nz, nday
    integer  :: ilon, ilat, iz, iday
    integer  :: igrd         (IA,JA,itp_nh_a)
    integer  :: jgrd         (IA,JA,itp_nh_a)
    real(RP) :: hfact        (IA,JA,itp_nh_a)
    integer  :: kgrd    (KA,2,IA,JA,itp_nh_a)
    real(RP) :: vfact   (KA,  IA,JA,itp_nh_a)
    
    real(RP), allocatable :: X_org(:,:)
    real(RP), allocatable :: Y_org(:,:)
    logical  :: zonal, pole
    real(RP) :: wsum(KA,IA,JA)
    real(RP) :: work(KA,IA,JA)
    integer  :: kref
    integer  :: k, i, j
    !----------------------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//varname

    !read data4d2, data1d_lon, data1d_lat, data1d_zg
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d2, fname, varname )
    nz   = size(data4d2,1)
    nlon = size(data4d2,2)
    nlat = size(data4d2,3)
    nday = size(data4d2,4)
    !LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//'lon'
    !LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//'lat'
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//'lon_'//varname
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//'lat_'//varname
    LOG_INFO("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'READ DATA, FILE NAME:', fname//': '//'zg'
    !call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon' )
    !call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat' )
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lon, fname, 'lon_'//varname )
    call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_lat, fname, 'lat_'//varname )
    if ( varname == 'AIRBC' .or. varname == 'AIROC' .or. varname == 'AIRSO2' ) then
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data1d_zg,  fname, 'rcpz' )
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d2_zg, fname, 'GPH'   )
       allocate( data4d_zg(nz,nlon,nlat,nday) )
       do iday = 1, nday
          do ilat = 1, nlat
             do ilon = 1, nlon
                do iz = 1, nz
                   data4d_zg(iz,ilon,ilat,iday) = data1d_zg(iz) + data4d2_zg(1,ilon,ilat,iday)
                enddo
             enddo
          enddo
       enddo
       deallocate( data1d_zg, data4d2_zg )
    else
       call ATMOS_PHY_AE_SPRINTARS_COMMON_READ_NETCDF_DATA( data4d_zg,  fname, 'GPH'  )
    endif

    allocate( X_org(nlon,nlat), Y_org(nlon,nlat) )
    allocate( data4d(nz+2,nlon,nlat,nday) )
    allocate( lon(nlon,nlat), lat(nlon,nlat), z(nz+2,nlon,nlat,nday) )
    data4d(:,:,:,:) = CONST_UNDEF
    do ilat = 1, nlat
       do ilon = 1, nlon
          lon(ilon,ilat) = data1d_lon(ilon) * CONST_D2R
          lat(ilon,ilat) = data1d_lat(ilat) * CONST_D2R
       enddo
    enddo

    do iz = 1, nz
       z     (iz+2,:,:,:) = data4d_zg(iz,:,:,:)
       data4d(iz+2,:,:,:) = data4d2  (iz,:,:,:)
    enddo

    do iday = 1, nday
       iz = nz+2
       call INTERP_domain_compatibility( lon(:,:),         & ! [IN]
                                         lat(:,:),         & ! [IN]
                                         z  (iz,:,:,iday), & ! [IN]
                                         LON_REAL(:,:),    & ! [IN]
                                         LAT_REAL(:,:),    & ! [IN]
                                         GDZ  (KE,:,:),    & ! [IN]
                                         GDZM (KE,:,:)     ) ! [IN]
    enddo

    if ( nlon == 1 .or. nlat == 1 ) then
       LOG_ERROR("SCALE_ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D",*) 'LINER interpolation requires nx, ny > 1'
       LOG_ERROR_CONT(*)               'Use "DIST-WEIGHT" as INTRP_TYPE of PARAM_MKINIT_REAL_ATMOS'
       call PRC_abort
    endif

    !$omp parallel do collapse(2)
    do ilat = 1, nlat
    do ilon = 1, nlon
       lat(ilon,ilat) = sign( min( abs(lat(ilon,ilat)), PI * 0.499999_RP ), lat(ilon,ilat) )
    enddo
    enddo

    call MAPPROJECTION_lonlat2xy( nlon, 1, nlon, & ! [IN]
                                  nlat, 1, nlat, & ! [IN]
                                  lon(:,:),      & ! [IN]
                                  lat(:,:),      & ! [IN]
                                  X_org(:,:),    & ! [OUT]
                                  Y_org(:,:)     ) ! [OUT]
    
    zonal = ( maxval(lon) - minval(lon) ) > 2.0_RP * PI * 0.9_RP
    pole = ( maxval(lat) > PI * 0.5_RP * 0.9_RP ) .or. ( minval(lat) < - PI * 0.5_RP * 0.9_RP )

    do iday = 1, nday
       
       call INTERP_factor3d( nz+2, 1, nz+2,          & ! [IN]
                             nlon, nlat,             & ! [IN]
                             KA, KS, KE,             & ! [IN]
                             IA, JA,                 & ! [IN]
                             X_org(:,:), Y_org(:,:), & ! [IN]
                             z (:,:,:,iday),         & ! [IN]
                             CX(:), CY(:),           & ! [IN]
                             GDZ    (:,:,:),         & ! [IN]
                             igrd   (    :,:,:),     & ! [OUT]
                             jgrd   (    :,:,:),     & ! [OUT]
                             hfact  (    :,:,:),     & ! [OUT]
                             kgrd   (:,:,:,:,:),     & ! [OUT]
                             vfact  (:,  :,:,:),     & ! [OUT]
                             flag_extrap = .false.,  & ! [IN]
                             zonal = zonal,          & ! [IN]
                             pole  = pole            ) ! [IN]

       call INTERP_interp3d( itp_nh_a,                  &
                             nz+2, 1, nz+2,             &
                             nlon, nlat,                &
                             KA, KS, KE,                &
                             IA, JA,                    &
                             igrd(:,:,:), jgrd(:,:,:),  & ! [IN]
                             hfact(:,:,:),              & ! [IN]
                             kgrd(:,:,:,:,:),           & ! [IN]
                             vfact(:,:,:,:),            & ! [IN]
                             z(:,:,:,iday), GDZ(:,:,:), & ! [IN]
                             data4d(:,:,:,iday),        & ! [IN]
                             DATA  (:,:,:,iday),        & ! [OUT]
                             threshold_undef = 1.0_RP,  & ! [IN]
                             wsum = wsum(:,:,:),        & ! [OUT]
                             val2 = work(:,:,:)         ) ! [OUT]

       do j = 1, JA
       do i = 1, IA
          do k = KS, KA
             if ( DATA(k,i,j,iday) /= CONST_UNDEF ) then
                kref = k
                exit
             endif
          enddo
          do k = kref-1, KS, -1
             DATA(k,i,j,iday) = DATA(k+1,i,j,iday) * log( ( GDZ(k,i,j) - GDZS(i,j) ) / 1.0E-4_RP ) / log( ( GDZ(k+1,i,j) - GDZS(i,j) ) / 1.0E-4_RP ) * ( 1.0_RP - wsum(k,i,j) ) &
                                + work(k,i,j) * wsum(k,i,j)
             ! DATA(k,i,j,iday) = 0.0_RP
          enddo
          do k = kref+1, KE
             if ( DATA(k,i,j,iday) == CONST_UNDEF ) DATA(k,i,j,iday) = DATA(k-1,i,j,iday)
             ! DATA(k,i,j,iday) = 0.0_RP
          enddo
       enddo
       enddo

    enddo

    deallocate( lon, lat, z, data1d_lon, data1d_lat, data4d, data4d2 )
    deallocate( X_org, Y_org )
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_READ_DATA_DAILY_4D
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS INSTAB (instability)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB &
       ( QTRC,             & ! [MODIFIED]
         ONEWQ, DELP, GDPM ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC (KA,  IA,JA) !< aerosol mixing ratio [kg/kg]
    
    ! input
    !-------------------------
    logical,  intent(in)    :: ONEWQ(KA-1,IA,JA) !< switch for instability
    real(RP), intent(in)    :: DELP (KA,  IA,JA) !< GDPM(k-1,i,j)-GDPM(k,i,j)
    real(RP), intent(in)    :: GDPM (0:KA,IA,JA) !< pressure at half levels

    ! internal work
    !-------------------------
    real(RP) :: NEWQTRC
    
    ! others
    !-------------------------
    integer  :: k, i, j

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE-1
          if ( ONEWQ(k,i,j) ) then
             NEWQTRC = ( QTRC(k,i,j) * DELP(k,i,j) + QTRC(k+1,i,j) * DELP(k+1,i,j) ) &
                     / ( GDPM(k-1,i,j) - GDPM(k+1,i,j) )
             QTRC(k  ,i,j) = NEWQTRC
             QTRC(k+1,i,j) = NEWQTRC
          endif
       enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_INSTAB
  !-----------------------------------------------------------------------------

  
  
  !-----------------------------------------------------------------------------
  !> SPRINTARS EVPAER (wet scavenging; re-emission from rain)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER &
       ( QTRCTM, QWTDEP,                & ! [MODIFIED]
         WETF,                          & ! [OUT]
         AERORN, AEROMS, CLTORN,        & ! [IN]
         DENS, FRAIN, DELP, DELZ, dt_AE ) ! [IN]
    
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCTM(KA,  IA,JA) !< outside cloud water [kg/kg]
    real(RP), intent(inout) :: QWTDEP(KA,  IA,JA) !< mixing ratio by wet deposition [kg/kg]
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: WETF       (IA,JA) !< wet deposition flux

    ! input
    !-------------------------
    real(RP), intent(in)    :: AERORN(KA,  IA,JA) !< number density of aerosol entering into raindrops [1/s/m3] <=> collision frequency between aerosol and rain
    real(RP), intent(in)    :: AEROMS             !< aerosol particle mass                             [kg]
    real(RP), intent(in)    :: CLTORN(KA,  IA,JA) !< mixing ratio moving cloud to rain particles       [kg/kg]
    real(RP), intent(in)    :: DENS  (KA,  IA,JA)
    real(RP), intent(in)    :: FRAIN (0:KA,IA,JA)
    real(RP), intent(in)    :: DELP  (KA,  IA,JA)
    real(RP), intent(in)    :: DELZ  (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: FDR(IA,JA)      !< aerosol flux in raindrops                        [kg/m2/s]
    real(RP) :: DF              !< FRAIN(k,i,j) - FRAIN(k-1,i,j)
    real(RP) :: QWTEVP          !< aerosol mass mixing ratio moving raindrop to air [kg/kg]
    integer  :: k, i, j

    !----------------------------------------------------

    FDR (:,:) = 0.0_RP
    WETF(:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE
       do k = KE, KS, -1
          
          FDR(i,j) = FDR(i,j) + DELZ(k,i,j) * (   AERORN(k,i,j) * AEROMS &
                                                + CLTORN(k,i,j) / real(dt_AE,kind=RP) * DENS(k,i,j) )

          if ( FRAIN(k,i,j) > 0.0_RP .and. FRAIN(k,i,j) > FRAIN(k-1,i,j) ) then
             DF     = FRAIN(k,i,j) - FRAIN(k-1,i,j)
             QWTEVP = DF / FRAIN(k,i,j) * FDR(i,j) * real(dt_AE,kind=RP) * GRAV / DELP(k,i,j)
             QTRCTM(k,i,j) = QTRCTM(k,i,j) + QWTEVP
             QWTDEP(k,i,j) = QWTDEP(k,i,j) - QWTEVP
             FDR     (i,j) = FDR(i,j) * FRAIN(k-1,i,j) / FRAIN(k,i,j)
          endif

          WETF(i,j) = WETF(i,j) + QWTDEP(k,i,j) * DELP(k,i,j)

       enddo
    enddo
    enddo

    do j = JS, JE
    do i = IS, IE
       WETF(i,j) = max( WETF(i,j) / GRAV / real(dt_AE,kind=RP), 0.0_RP )
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_EVPAER
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS setup for dry deposition
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET &
       ( QMTX,                   & ! [OUT]
         DENS, DELP, CDVE, dt_AE ) ! [IN]

    implicit none

    
    ! output
    !-------------------------
    real(RP), intent(out) :: QMTX(KA,IA,JA,-1:1)
    
    ! input
    !-------------------------
    real(RP), intent(in)  :: DENS(KA,IA,JA)
    real(RP), intent(in)  :: DELP(KA,IA,JA)
    real(RP), intent(in)  :: CDVE   (IA,JA)
    real(DP), intent(in)  :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: DFZ(KA,IA,JA)
    integer  :: k, i, j
    
    !----------------------------------------------------

    DFZ(:,:,:) = 0.0_RP !vertical diffusion is not calculated by SPRINTARS
    
    do j = JS, JE
    do i = IS, IE

       QMTX(KE,i,j,-1) = - DFZ(KE,i,j) * real(dt_AE,kind=RP)
       QMTX(KE,i,j, 0) = DELP(KE,i,j) / GRAV + DFZ(KE,i,j) * real(dt_AE,kind=RP)
       QMTX(KE,i,j, 1) = 0.0_RP
   
       do k = KS+1, KE-1
          QMTX(k,i,j,-1) = - DFZ(k,i,j) * real(dt_AE,kind=RP)
          QMTX(k,i,j, 0) = DELP(k,i,j) / GRAV + ( DFZ(k,i,j) + DFZ(k+1,i,j) ) * real(dt_AE,kind=RP)
          QMTX(k,i,j, 1) = - DFZ(k,i,j) * real(dt_AE,kind=RP)
       enddo
       
       QMTX(KS,i,j,-1) = 0.0_RP
       QMTX(KS,i,j, 0) = DELP(KS,i,j) / GRAV + DFZ(KS+1,i,j) * real(dt_AE,kind=RP) + DENS(KS,i,j) * CDVE(i,j) * real(dt_AE,kind=RP)
       QMTX(KS,i,j, 1) = - DFZ(KS+1,i,j) * real(dt_AE,kind=RP)

       do k = KE-1, KS, -1
          QMTX(k,i,j, 0) = QMTX(k,i,j, 0) - QMTX(k+1,i,j,-1) / QMTX(k+1,i,j, 0) * QMTX(k,i,j, 1)
       enddo
       
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_DRYSET
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS dry deposition for aerosols
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP &
       ( QTRCA,                        & ! [MODIFIED]
         DRYF,                         & ! [OUT]
         QMTX, DENS, DELP, CDVE, dt_AE ) ! [IN]

    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,IA,JA)
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: DRYF    (IA,JA)
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: QMTX (KA,IA,JA,-1:1)
    real(RP), intent(in)    :: DENS (KA,IA,JA)
    real(RP), intent(in)    :: DELP (KA,IA,JA)
    real(RP), intent(in)    :: CDVE    (IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: DFZ   (KA,  IA,JA)
    real(RP) :: FLUX  (0:KA,IA,JA)
    real(RP) :: GTQ   (KA,  IA,JA)
    real(RP) :: DRYORI     (IA,JA)
    real(RP) :: DRYREM     (IA,JA)
    integer  :: k, i, j
    
    !----------------------------------------------------

    DFZ  (:,:,:) = 0.0_RP !vertical diffusion is not calculated by SPRINTARS
    DRYORI (:,:) = 0.0_RP
    DRYREM (:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE

       FLUX(KE,i,j) = 0.0_RP
       do k = KS, KE-1
          FLUX(k,i,j) = - DFZ(k,i,j) * ( QTRCA(k+1,i,j) - QTRCA(k,i,j) )
       enddo
       FLUX(KS-1,i,j) = - DENS(KS,i,j) * CDVE(i,j) * QTRCA(KS,i,j)

       do k = KS, KE
          DRYORI(i,j) = DRYORI(i,j) + QTRCA(k,i,j) * DELP(k,i,j)
          GTQ (k,i,j) = FLUX(k-1,i,j) - FLUX(k,i,j)
       enddo

       do k = KE-1, KS, -1
          GTQ(k,i,j) = GTQ(k,i,j) - GTQ(k+1,i,j) / QMTX(k+1,i,j, 0) * QMTX(k,i,j, 1)
       enddo

       GTQ(KS,i,j) = GTQ(KS,i,j) / QMTX(KS,i,j, 0)
       do k = KS+1, KE
          GTQ(k,i,j) = ( GTQ(k,i,j) - GTQ(k-1,i,j)*QMTX(k,i,j,-1) ) / QMTX(k,i,j,0)
       enddo
       
       do k = KS, KE
          QTRCA(k,i,j) = QTRCA(k,i,j) + GTQ(k,i,j) * real(dt_AE,kind=RP)
          DRYREM (i,j) = DRYREM (i,j) + QTRCA(k,i,j) * DELP(k,i,j)
       enddo
       DRYF(i,j) = max( ( DRYORI(i,j) - DRYREM(i,j) ) / GRAV / real(dt_AE,kind=RP), 0.0_RP )
              
    enddo
    enddo
    
    return

  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_DRYDEP
  !----------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS gravitational settling for aerosol (original)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST &
       ( QTRCA,                        & ! [MODIFIED]
         QLOAD, GRAVF,                 & ! [OUT]
         GDPM, DELP, DENS, VTER, dt_AE ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,  IA,JA)

    ! output
    !-------------------------
    real(RP), intent(out)   :: QLOAD(KA,  IA,JA)
    real(RP), intent(out)   :: GRAVF     (IA,JA)
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: GDPM (0:KA,IA,JA)
    real(RP), intent(in)    :: DELP (KA,  IA,JA)
    real(RP), intent(in)    :: DENS (KA,  IA,JA)
    real(RP), intent(in)    :: VTER (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: SEM    (KA)
    real(RP) :: SEMSUM
    real(RP) :: QTRCORI(KA)
    real(RP) :: PPLUS  (KA)
    real(RP) :: PMINUS (KA)
    integer  :: IFLAG  (KA)
    integer  :: JPL, JMI
    real(RP) :: GDPMP, GDPMP1
    real(RP) :: GDPMM, GDPMM1
    real(RP) :: FAC1, FAC2
    
    ! others
    !-------------------------
    integer  :: k, i, j
    integer  :: kk

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE

       !-------------------------
       do k = KS, KE
          QTRCORI(k) = QTRCA(k,i,j)
          PPLUS  (k) = GDPM(k,  i,j) + DENS(k,i,j) * GRAV * VTER(k,i,j) * real(dt_AE,kind=RP)
          PMINUS (k) = GDPM(k-1,i,j) + DENS(k,i,j) * GRAV * VTER(k,i,j) * real(dt_AE,kind=RP)
       enddo

       do k = KS, KE
          if ( PPLUS(k) > GDPM(KS-1,i,j) ) then
             IFLAG(k) = 1   ! perfect settling under ground
          elseif ( PMINUS(k) > GDPM(KS-1,i,j) ) then
             IFLAG(k) = 2   ! partial settling under ground
          else
             IFLAG(k) = 0
          endif
       enddo

       do k = KE, KS, -1
          JPL = k
          JMI = k
          GDPMP  = GDPM(k-1,i,j)
          GDPMP1 = GDPM(k,  i,j)
          GDPMM  = GDPM(k-1,i,j)
          GDPMM1 = GDPM(k,  i,j)
          do kk = k, KS+1, -1
             FAC1 = GDPM(kk-1,i,j)
             FAC2 = GDPM(kk-2,i,j)
             if ( PPLUS (k) >= FAC1 ) then
                JPL = kk-1
                GDPMP  = FAC2
                GDPMP1 = FAC1
             endif
             if ( PMINUS(k) >= FAC1 ) then
                JMI = kk-1
                GDPMM  = FAC2
                GDPMM1 = FAC1
             endif
          enddo

          FAC1 = GDPMP - GDPMP1
          if ( IFLAG(k) == 2 ) then
             QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + ( GDPMP - PPLUS(k) ) * QTRCORI(k) / FAC1
          elseif ( IFLAG(k) == 0 ) then
             if ( JPL - JMI > 0 ) then
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + ( GDPMP     - PPLUS(k) ) * QTRCORI(k) / FAC1
                QTRCA(JMI,i,j) = QTRCA(JMI,i,j) + ( PMINUS(k) - GDPMM1   ) * QTRCORI(k) / ( GDPMM - GDPMM1 )
             else
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + DELP(k,i,j) * QTRCORI(k) / FAC1
             endif
          endif

          do kk = KS, k-2
             if (      ( IFLAG(k) == 0 .and. kk >= JMI+1 .and. kk <= JPL+1 ) &
                  .or. ( IFLAG(k) == 2 .and. kk <= JPL+1                   ) ) then
                QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k)
             endif
          enddo

          QTRCA(k,i,j) = - QTRCORI(k) + QTRCA(k,i,j)

          ! do kk = KS, k
          !    FAC1 = QTRCORI(k) * DELP(k,i,j)
          !    if ( IFLAG(k) >= 1 ) then
          !       FAC1 = FAC1 / ( PMINUS(k) - GDPM(k,i,j) )
          !    elseif ( kk >= JMI ) then
          !       FAC1 = FAC1 / ( GDPMM     - GDPM(k,i,j) )
          !    else
          !       FAC1 = 0.0_RP
          !    endif
          !    QLOAD(kk,i,j) = QLOAD(kk,i,j) + FAC1
          ! enddo

          if ( IFLAG(k) == 1 ) then
             SEM(k) = QTRCORI(k) * ( PMINUS(k) - PPLUS(k)       )
          elseif ( IFLAG(k) == 2 ) then
             SEM(k) = QTRCORI(k) * ( PMINUS(k) - GDPM(KS-1,i,j) )
          else
             SEM(k) = 0.0_RP
          endif

       enddo
       
       SEMSUM = 0.0_RP
       do k = KS, KE
          SEMSUM = SEMSUM + SEM(k)
       enddo
       GRAVF(i,j) = max( SEMSUM/real(dt_AE,kind=RP)/GRAV, 0.0_RP )
       !-------------------------
       
    enddo
    enddo
    QLOAD(:,:,:) = QTRCA(:,:,:)
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS gravitational settling for aerosol (original)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST_Z &
       ( QTRCA,                        & ! [MODIFIED]
         QLOAD, GRAVF,                 & ! [OUT]
         GDZM, DELZ, DENS, VTER, dt_AE ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,IA,JA) ! [kg/kg]

    ! output
    !-------------------------
    real(RP), intent(out)   :: QLOAD(KA,IA,JA)
    real(RP), intent(out)   :: GRAVF   (IA,JA)

    ! input
    !-------------------------
    real(RP), intent(in)    :: GDZM  (0:KA,IA,JA)
    real(RP), intent(in)    :: DELZ  (KA,  IA,JA)
    real(RP), intent(in)    :: DENS  (KA,  IA,JA)
    real(RP), intent(in)    :: VTER  (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: QTRCORI(KA)
    real(RP) :: SEMSUM
    real(RP) :: ZPLUS (KA)      !< height of the upper surface 
    real(RP) :: ZMINUS(KA)      !< height of the lower surface 
    integer  :: IFLAG (KA)
    integer  :: JPL             !< z index of upper surface
    integer  :: JMI             !< z index of lower surface
    real(RP) :: GDZMP           !< height of the lower surface of the volume element at JPL
    real(RP) :: GDZMP1          !< height of the upper surface of the volume element at JPL
    real(RP) :: GDZMM           !< height of the lower surface of the volume element at JMI
    real(RP) :: GDZMM1          !< height of the upper surface of the volume element at JMI
    real(RP) :: DZPM            !< ZPLUS(k) - ZMINUS(k)
    
    ! others
    !-------------------------
    integer  :: k, i, j
    integer  :: kk

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE
          
       SEMSUM = 0.0_RP
       
       !-------------------------
       do k = KS, KE
          QTRCORI(k) = QTRCA(k,  i,j)
          ZPLUS  (k) = GDZM (k,  i,j) - VTER(k,i,j) * real(dt_AE,kind=RP)
          ZMINUS (k) = GDZM (k-1,i,j) - VTER(k,i,j) * real(dt_AE,kind=RP)
       enddo

       do k = KS, KE
          if ( ZPLUS(k) < GDZM(KS-1,i,j) ) then
             IFLAG(k) = 1   ! perfect settling under ground
          elseif ( ZMINUS(k) < GDZM(KS-1,i,j) ) then
             IFLAG(k) = 2   ! partial settling under ground
          else
             IFLAG(k) = 0
          endif
       enddo

       do k = KE, KS, -1

          do kk = k, KS, -1
             if ( ZPLUS(k) <= GDZM(kk,i,j) .and. ZPLUS(k) > GDZM(kk-1,i,j) ) then
                JPL = kk
                GDZMP  = GDZM(kk-1,i,j)
                GDZMP1 = GDZM(kk,  i,j)
                exit
             endif
          enddo
          do kk = k, KS, -1
             if ( ZMINUS(k) <= GDZM(kk,i,j) .and. ZMINUS(k) > GDZM(kk-1,i,j) ) then
                JMI = kk
                GDZMM  = GDZM(kk-1,i,j)
                GDZMM1 = GDZM(kk,  i,j)
                exit
             endif
          enddo

          DZPM = ZPLUS(k)-ZMINUS(k)
          if ( IFLAG(k) == 0 ) then
             if ( JPL - JMI > 0 ) then
                do kk = JMI, JPL, 1
                   if ( kk == JMI ) then
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (GDZM(kk,i,j)-ZMINUS(k)) / DZPM
                   elseif ( kk == JPL ) then
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (ZPLUS(k)-GDZM(kk-1,i,j)) / DZPM
                   else
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DZPM
                   endif
                enddo
             else
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(JPL,i,j) * DELZ(k,i,j) / DELZ(JPL,i,j)
             endif
          elseif ( IFLAG(k) == 2 ) then ! partial settling under ground
             do kk = KS, JPL, 1
                if ( kk == JPL ) then
                   QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (ZPLUS(k)-GDZM(kk-1,i,j)) / DZPM
                else
                   QTRCA(kk,i,j) = QTRCA(kk,i,j) + QTRCORI(k) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DZPM
                endif
             enddo
          endif

          QTRCA(k,i,j) = QTRCA(k,i,j) - QTRCORI(k)
          
          if ( IFLAG(k) == 1 ) then     ! perfect settling under ground
             SEMSUM = SEMSUM + QTRCORI(k) * DENS(k,i,j) * DELZ(k,i,j)
          elseif ( IFLAG(k) == 2 ) then ! partial settling under ground
             SEMSUM = SEMSUM + QTRCORI(k) * DENS(k,i,j) * DELZ(k,i,j) * (GDZM(KS-1,i,j)-ZMINUS(k)) / DZPM
          endif

       enddo
       
       GRAVF(i,j) = max( SEMSUM/real(dt_AE,kind=RP), 0.0_RP )
       !-------------------------
       
    enddo
    enddo
    QLOAD(:,:,:) = QTRCA(:,:,:)
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_GRAVST_Z
  !-----------------------------------------------------------------------------

  
      
  !-----------------------------------------------------------------------------
  !> SPRINTARS gravitational settling for aerosol
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST &
       ( QTRCA,                        & ! [MODIFIED]
         QLOAD, GRAVF,                 & ! [OUT]
         GDPM, DELP, DENS, VTER, dt_AE ) ! [IN]
    
    implicit none
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,  IA,JA)

    ! output
    !-------------------------
    real(RP), intent(out)   :: QLOAD(KA,  IA,JA)
    real(RP), intent(out)   :: GRAVF     (IA,JA)
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: GDPM (0:KA,IA,JA)
    real(RP), intent(in)    :: DELP (KA,  IA,JA)
    real(RP), intent(in)    :: DENS (KA,  IA,JA)
    real(RP), intent(in)    :: VTER (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: RDZ
    real(RP) :: GRVEF
    real(RP) :: QTRCB(KA,IA,JA)
    integer  :: k, i, j
    
    !----------------------------------------------------

    GRAVF(:,:) = 0.0_RP
    QTRCB(:,:,:) = QTRCA(:,:,:)

    do j = JS, JE
    do i = IS, IE
       do k = KE, KS, -1
          
          RDZ          = DELP(k,i,j) / GRAV
          QTRCA(k,i,j) = QTRCA(k,i,j) + GRAVF(i,j) / RDZ * real(dt_AE,kind=RP)
          GRVEF        = DENS(k,i,j) * VTER(k,i,j) / RDZ
          GRVEF        = tanh( GRVEF * real(dt_AE,kind=RP) ) / real(dt_AE,kind=RP)
          GRAVF(i,j)   = QTRCB(k,i,j) * GRVEF * RDZ
          QTRCA(k,i,j) = QTRCA(k,i,j) - GRAVF(i,j) / RDZ * real(dt_AE,kind=RP)
          QLOAD(k,i,j) = QTRCA(k,i,j)

       enddo
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_GRAVST
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS rainfall (microphysics scheme) for aerosol
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL &
       ( QTRCTM,                       & ! [MODIFIED]
         WETF,                         & ! [OUT]
         QWTDEP,                       & ! [IN]
         GDPM, DELP, DENS, VTER, dt_AE ) ! [IN]
    
    implicit none
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCTM(KA,  IA,JA) !< outside cloud water [kg/kg]
    
    ! output
    !-------------------------
    real(RP), intent(out)   :: WETF       (IA,JA) !< wet deposition flux
        
    ! input
    !-------------------------
    real(RP), intent(in)    :: QWTDEP(KA,  IA,JA) !< mixing ratio by wet deposition [kg/kg]
    real(RP), intent(in)    :: GDPM  (0:KA,IA,JA)
    real(RP), intent(in)    :: DELP  (KA,  IA,JA)
    real(RP), intent(in)    :: DENS  (KA,  IA,JA)
    real(RP), intent(in)    :: VTER  (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: RDZ
    real(RP) :: WEF
    integer  :: k, i, j
    
    !----------------------------------------------------

    WETF(:,:) = 0.0_RP

    do j = JS, JE
    do i = IS, IE
       do k = KE, KS, -1
          
          RDZ           = DELP(k,i,j) / GRAV
          QTRCTM(k,i,j) = QTRCTM(k,i,j) + WETF(i,j) / RDZ * real(dt_AE,kind=RP)
          WEF           = DENS(k,i,j) * VTER(k,i,j) / RDZ
          WEF           = tanh( WEF * real(dt_AE,kind=RP) ) / real(dt_AE,kind=RP)
          WETF(i,j)     = QWTDEP(k,i,j) * WEF * RDZ
          QTRCTM(k,i,j) = QTRCTM(k,i,j) + QWTDEP(k,i,j)- WETF(i,j) / RDZ * real(dt_AE,kind=RP)

       enddo
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_RAINFALL
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS rainfall (microphysics scheme) for aerosol (original)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL &
       ( QTRCA,                        & ! [MODIFIED]
         WETF,                         & ! [OUT]
         QWTDEP,                       & ! [IN]
         GDPM, DELP, DENS, VTER, dt_AE ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,IA,JA)

    ! output
    !-------------------------
    real(RP), intent(out)   :: WETF    (IA,JA)
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: QWTDEP(KA,  IA,JA) !< mixing ratio by wet deposition [kg/kg]
    real(RP), intent(in)    :: GDPM  (0:KA,IA,JA)
    real(RP), intent(in)    :: DELP  (KA,  IA,JA)
    real(RP), intent(in)    :: DENS  (KA,  IA,JA)
    real(RP), intent(in)    :: VTER  (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: SEM    (KA)
    real(RP) :: SEMSUM
    real(RP) :: PPLUS  (KA)
    real(RP) :: PMINUS (KA)
    integer  :: IFLAG  (KA)
    integer  :: JPL, JMI
    real(RP) :: GDPMP, GDPMP1
    real(RP) :: GDPMM, GDPMM1
    real(RP) :: FAC1, FAC2
    
    ! others
    !-------------------------
    integer  :: k, i, j
    integer  :: kk

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE

       !-------------------------
       do k = KS, KE
          PPLUS  (k) = GDPM(k,  i,j) + DENS(k,i,j) * GRAV * VTER(k,i,j) * real(dt_AE,kind=RP)
          PMINUS (k) = GDPM(k-1,i,j) + DENS(k,i,j) * GRAV * VTER(k,i,j) * real(dt_AE,kind=RP)
       enddo

       do k = KS, KE
          if ( PPLUS(k) > GDPM(KS-1,i,j) ) then
             IFLAG(k) = 1   ! perfect settling under ground
          elseif ( PMINUS(k) > GDPM(KS-1,i,j) ) then
             IFLAG(k) = 2   ! partial settling under ground
          else
             IFLAG(k) = 0
          endif
       enddo

       do k = KE, KS, -1
          JPL = k
          JMI = k
          GDPMP  = GDPM(k-1,i,j)
          GDPMP1 = GDPM(k,  i,j)
          GDPMM  = GDPM(k-1,i,j)
          GDPMM1 = GDPM(k,  i,j)
          do kk = k, KS+1, -1
             if ( PPLUS(k) >= GDPM(kk-1,i,j) .and. PPLUS(k) < GDPM(kk-2,i,j) ) then
                JPL = kk-1
                GDPMP  = GDPM(kk-2,i,j)
                GDPMP1 = GDPM(kk-1,i,j)
                exit
             endif
          enddo
          do kk = k, KS+1, -1
             if ( PMINUS(k) >= GDPM(kk-1,i,j) .and. PMINUS(k) < GDPM(kk-2,i,j) ) then
                JMI = kk-1
                GDPMM  = GDPM(kk-2,i,j)
                GDPMM1 = GDPM(kk-1,i,j)
                exit
             endif
          enddo

          FAC1 = GDPMP - GDPMP1
          if ( IFLAG(k) == 2 ) then
             QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + ( GDPMP - PPLUS(k) ) * QWTDEP(k,i,j) / FAC1
          elseif ( IFLAG(k) == 0 ) then
             if ( JPL - JMI > 0 ) then
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + ( GDPMP     - PPLUS(k) ) * QWTDEP(k,i,j) / FAC1
                QTRCA(JMI,i,j) = QTRCA(JMI,i,j) + ( PMINUS(k) - GDPMM1   ) * QWTDEP(k,i,j) / ( GDPMM - GDPMM1 )
             else
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + DELP(k,i,j) * QWTDEP(k,i,j) / FAC1
             endif
          endif

          do kk = KS, k-2
             if (      ( IFLAG(k) == 0 .and. kk >= JMI+1 .and. kk <= JPL+1 ) &
                  .or. ( IFLAG(k) == 2 .and. kk <= JPL+1                   ) ) then
                QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j)
             endif
          enddo

          ! QTRCA(k,i,j) = - QTRCORI(k) + QTRCA(k,i,j)

          ! do kk = KS, k
          !    FAC1 = QTRCORI(k) * DELP(k,i,j)
          !    if ( IFLAG(k) >= 1 ) then
          !       FAC1 = FAC1 / ( PMINUS(k) - GDPM(k,i,j) )
          !    elseif ( kk >= JMI ) then
          !       FAC1 = FAC1 / ( GDPMM     - GDPM(k,i,j) )
          !    else
          !       FAC1 = 0.0_RP
          !    endif
          !    QLOAD(kk,i,j) = QLOAD(kk,i,j) + FAC1
          ! enddo

          if ( IFLAG(k) == 1 ) then     ! perfect settling under ground
             SEM(k) = QWTDEP(k,i,j) * ( PMINUS(k) - PPLUS(k)       )
          elseif ( IFLAG(k) == 2 ) then ! partial settling under ground
             SEM(k) = QWTDEP(k,i,j) * ( PMINUS(k) - GDPM(KS-1,i,j) )
          else
             SEM(k) = 0.0_RP
          endif

       enddo
       
       SEMSUM = 0.0_RP
       do k = KS, KE
          SEMSUM = SEMSUM + SEM(k)
       enddo
       WETF(i,j) = max( SEMSUM/real(dt_AE,kind=RP)/GRAV, 0.0_RP )
       !-------------------------
       
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS rainfall (microphysics scheme) for aerosol (original)
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z &
       ( QTRCA,                        & ! [MODIFIED]
         WETF,                         & ! [OUT]
         QWTDEP,                       & ! [IN]
         GDZM, DELZ, DENS, VTER, dt_AE ) ! [IN]
    implicit none

    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRCA(KA,IA,JA) ! [kg/kg]

    ! output
    !-------------------------
    real(RP), intent(out)   :: WETF    (IA,JA) ! [kg/m2/s]
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: QWTDEP(KA,  IA,JA) !< mixing ratio by wet deposition [kg/kg]
    real(RP), intent(in)    :: GDZM  (0:KA,IA,JA)
    real(RP), intent(in)    :: DELZ  (KA,  IA,JA)
    real(RP), intent(in)    :: DENS  (KA,  IA,JA)
    real(RP), intent(in)    :: VTER  (KA,  IA,JA)
    real(DP), intent(in)    :: dt_AE

    ! internal work
    !-------------------------
    real(RP) :: SEMSUM
    real(RP) :: ZPLUS (KA)      !< height of the upper surface 
    real(RP) :: ZMINUS(KA)      !< height of the lower surface 
    integer  :: IFLAG (KA)
    integer  :: JPL             !< z index of upper surface
    integer  :: JMI             !< z index of lower surface
    real(RP) :: GDZMP           !< height of the lower surface of the volume element at JPL
    real(RP) :: GDZMP1          !< height of the upper surface of the volume element at JPL
    real(RP) :: GDZMM           !< height of the lower surface of the volume element at JMI
    real(RP) :: GDZMM1          !< height of the upper surface of the volume element at JMI
    real(RP) :: DZPM            !< ZPLUS(k) - ZMINUS(k)
    
    ! others
    !-------------------------
    integer  :: k, i, j
    integer  :: kk

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE
          
       SEMSUM = 0.0_RP
       
       !-------------------------
       do k = KS, KE
          ZPLUS (k) = GDZM(k,  i,j) - VTER(k,i,j) * real(dt_AE,kind=RP)
          ZMINUS(k) = GDZM(k-1,i,j) - VTER(k,i,j) * real(dt_AE,kind=RP)
       enddo

       do k = KS, KE
          if ( ZPLUS(k) < GDZM(KS-1,i,j) ) then
             IFLAG(k) = 1   ! perfect settling under ground
          elseif ( ZMINUS(k) < GDZM(KS-1,i,j) ) then
             IFLAG(k) = 2   ! partial settling under ground
          else
             IFLAG(k) = 0
          endif
       enddo

       do k = KE, KS, -1

          if ( IFLAG(k) == 0 .or. IFLAG(k) == 2 ) then
             do kk = k, KS, -1
                if ( ZPLUS(k) <= GDZM(kk,i,j) .and. ZPLUS(k) > GDZM(kk-1,i,j) ) then
                   JPL = kk
                   GDZMP  = GDZM(kk-1,i,j)
                   GDZMP1 = GDZM(kk,  i,j)
                   exit
                endif
             enddo
          endif

          if ( IFLAG(k) == 0 ) then
             do kk = k, KS, -1
                if ( ZMINUS(k) <= GDZM(kk,i,j) .and. ZMINUS(k) > GDZM(kk-1,i,j) ) then
                   JMI = kk
                   GDZMM  = GDZM(kk-1,i,j)
                   GDZMM1 = GDZM(kk,  i,j)
                   exit
                endif
             enddo
          endif
          
          DZPM = ZPLUS(k)-ZMINUS(k)
          if ( IFLAG(k) == 0 ) then
             if ( JPL - JMI > 0 ) then
                do kk = JMI, JPL, 1
                   if ( kk == JMI ) then
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (GDZM(kk,i,j)-ZMINUS(k)) / DZPM
                   elseif ( kk == JPL ) then
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (ZPLUS(k)-GDZM(kk-1,i,j)) / DZPM
                   else
                      QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DZPM
                   endif
                enddo
             else
                QTRCA(JPL,i,j) = QTRCA(JPL,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(JPL,i,j) * DELZ(k,i,j) / DELZ(JPL,i,j)
             endif
          elseif ( IFLAG(k) == 2 ) then ! partial settling under ground
             do kk = KS, JPL, 1
                if ( kk == JPL ) then
                   QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DELZ(kk,i,j) * (ZPLUS(k)-GDZM(kk-1,i,j)) / DZPM
                else
                   QTRCA(kk,i,j) = QTRCA(kk,i,j) + QWTDEP(k,i,j) * DENS(k,i,j) / DENS(kk,i,j) * DELZ(k,i,j) / DZPM
                endif
             enddo
          endif
          
          if ( IFLAG(k) == 1 ) then     ! perfect settling under ground
             SEMSUM = SEMSUM + QWTDEP(k,i,j) * DENS(k,i,j) * DELZ(k,i,j)
          elseif ( IFLAG(k) == 2 ) then ! partial settling under ground
             SEMSUM = SEMSUM + QWTDEP(k,i,j) * DENS(k,i,j) * DELZ(k,i,j) * (GDZM(KS-1,i,j)-ZMINUS(k)) / DZPM
          endif

       enddo
       
       WETF(i,j) = max( SEMSUM/real(dt_AE,kind=RP), 0.0_RP )
       !-------------------------
       
    enddo
    enddo
    
    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_ORIG_RAINFALL_Z
  !-----------------------------------------------------------------------------


  
  !-----------------------------------------------------------------------------
  !> SPRINTARS in precipitation 
  subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC &
       ( QTRC,                                                         & ! [INOUT]
         QTRCINR, QTRCINS,                                             & ! [OUT]
         QTRCINP,                                                      & ! [IN]
         FRAIN_MP_RAIN1, FRAIN_MP_SNOW1, FRAIN_MP_RAIN, FRAIN_MP_SNOW, & ! [IN]
         DENS1, DENS, THRES                                            ) ! [IN]
    implicit none
    
    ! modified
    !-------------------------
    real(RP), intent(inout) :: QTRC(KA,IA,JA) !< mixing ratio [kg/kg]

    ! output
    !-------------------------
    real(RP), intent(out)   :: QTRCINR(KA,IA,JA) !< inside rain(MP) [kg/kg]
    real(RP), intent(out)   :: QTRCINS(KA,IA,JA) !< inside snow(MP) [kg/kg]
    
    ! input
    !-------------------------
    real(RP), intent(in)    :: QTRCINP       (  KA,IA,JA) !< inside precipitation at 1 step before [kg/kg]
    real(RP), intent(in)    :: FRAIN_MP_RAIN1(0:KA,IA,JA) !< rain(MP) flux at 1 step before        [kg/m2/s]
    real(RP), intent(in)    :: FRAIN_MP_SNOW1(0:KA,IA,JA) !< snow(MP) flux at 1 step before        [kg/m2/s]
    real(RP), intent(in)    :: FRAIN_MP_RAIN (0:KA,IA,JA) !< rain(MP) flux                         [kg/m2/s]
    real(RP), intent(in)    :: FRAIN_MP_SNOW (0:KA,IA,JA) !< snow(MP) flux                         [kg/m2/s]
    real(RP), intent(in)    :: DENS1         (  KA,IA,JA) !< density at 1 step before              [kg/m3]
    real(RP), intent(in)    :: DENS          (  KA,IA,JA) !< density                               [kg/m3]
    real(RP), intent(in)    :: THRES                      !< precipitation threshold               [kg/m2/s]
    
    ! internal work
    !-------------------------
    real(RP) :: FRAIN_MP_PREC1 !< precipitation flux at 1 step before [kg/m2/s]
    real(RP) :: FRAIN_MP_PREC  !< precipitation flux                  [kg/m2/s]
    
    ! others
    !-------------------------
    integer  :: k, i, j

    !----------------------------------------------------

    do j = JS, JE
    do i = IS, IE
       do k = KS, KE

          FRAIN_MP_PREC1 = FRAIN_MP_RAIN1(k,i,j) + FRAIN_MP_SNOW1(k,i,j)
          FRAIN_MP_PREC  = FRAIN_MP_RAIN (k,i,j) + FRAIN_MP_SNOW (k,i,j)
          
          if ( FRAIN_MP_PREC > THRES .and. FRAIN_MP_PREC1 > THRES ) then

             if ( FRAIN_MP_PREC1 >= FRAIN_MP_PREC ) then
                QTRCINR(k,i,j) = FRAIN_MP_RAIN(k,i,j) / FRAIN_MP_PREC1 * QTRCINP(k,i,j) * DENS1(k,i,j) / DENS(k,i,j)
                QTRCINS(k,i,j) = FRAIN_MP_SNOW(k,i,j) / FRAIN_MP_PREC1 * QTRCINP(k,i,j) * DENS1(k,i,j) / DENS(k,i,j)
                if ( QTRC(k,i,j) - ( QTRCINR(k,i,j) + QTRCINS(k,i,j) ) < 0.0_RP ) then
                   QTRCINR(k,i,j) = FRAIN_MP_RAIN(k,i,j) / FRAIN_MP_PREC * QTRC(k,i,j)
                   QTRCINS(k,i,j) = FRAIN_MP_SNOW(k,i,j) / FRAIN_MP_PREC * QTRC(k,i,j)
                   QTRC   (k,i,j) = 0.0_RP
                else
                   QTRC   (k,i,j) = QTRC(k,i,j) - ( QTRCINR(k,i,j) + QTRCINS(k,i,j) )
                endif
                
             else
                QTRCINR(k,i,j) = FRAIN_MP_RAIN(k,i,j) / FRAIN_MP_PREC * QTRCINP(k,i,j) * DENS1(k,i,j) / DENS(k,i,j)
                QTRCINS(k,i,j) = FRAIN_MP_SNOW(k,i,j) / FRAIN_MP_PREC * QTRCINP(k,i,j) * DENS1(k,i,j) / DENS(k,i,j)
                if ( QTRC(k,i,j) - ( QTRCINR(k,i,j) + QTRCINS(k,i,j) ) < 0.0_RP ) then
                   QTRCINR(k,i,j) = FRAIN_MP_RAIN(k,i,j) / FRAIN_MP_PREC * QTRC(k,i,j)
                   QTRCINS(k,i,j) = FRAIN_MP_SNOW(k,i,j) / FRAIN_MP_PREC * QTRC(k,i,j)
                   QTRC   (k,i,j) = 0.0_RP
                else
                   QTRC   (k,i,j) = QTRC(k,i,j) - ( QTRCINR(k,i,j) + QTRCINS(k,i,j) )
                endif
             endif
                
          else
             QTRCINR(k,i,j) = 0.0_RP
             QTRCINS(k,i,j) = 0.0_RP
          endif

       enddo
    enddo
    enddo

    return
  end subroutine ATMOS_PHY_AE_SPRINTARS_COMMON_INPREC
  !-----------------------------------------------------------------------------



end module scale_atmos_phy_ae_SPRINTARS_common
