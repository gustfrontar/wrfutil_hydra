!-------------------------------------------------------------------------------
!> module INITIAL
!!
!! @par Description
!!          common subroutines for preparing initial data
!!
!! @author Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module mod_mkinit_common
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_atmos_grid_cartesC_index
  use scale_cpl_sfc_index

  use scale_prc, only: &
     PRC_abort
  use scale_const, only: &
     PI     => CONST_PI,     &
     GRAV   => CONST_GRAV,   &
     Pstd   => CONST_Pstd,   &
     Rdry   => CONST_Rdry,   &
     Rvap   => CONST_Rvap,   &
     CPdry  => CONST_CPdry,  &
     P00    => CONST_PRE00
  use scale_comm_cartesC, only: &
     COMM_vars8, &
     COMM_wait
  use scale_atmos_grid_cartesC, only: &
     CZ  => ATMOS_GRID_CARTESC_CZ,  &
     CX  => ATMOS_GRID_CARTESC_CX,  &
     CY  => ATMOS_GRID_CARTESC_CY,  &
     FZ  => ATMOS_GRID_CARTESC_FZ,  &
     FX  => ATMOS_GRID_CARTESC_FX,  &
     FY  => ATMOS_GRID_CARTESC_FY,  &
     CXG => ATMOS_GRID_CARTESC_CXG, &
     FXG => ATMOS_GRID_CARTESC_FXG, &
     FYG => ATMOS_GRID_CARTESC_FYG  
  use scale_atmos_grid_cartesC_real, only: &
     REAL_CZ => ATMOS_GRID_CARTESC_REAL_CZ, &
     REAL_FZ => ATMOS_GRID_CARTESC_REAL_FZ, &
     AREA    => ATMOS_GRID_CARTESC_REAL_AREA
  use scale_atmos_hydrometeor, only: &
     HYDROMETEOR_LHV => ATMOS_HYDROMETEOR_LHV       
  use scale_atmos_hydrostatic, only: &
     HYDROSTATIC_buildrho        => ATMOS_HYDROSTATIC_buildrho

  use scale_atmos_dyn_dgm_operator, only: &
     fem_mesh3D => mesh3D, &
     fem_elem3D => elem3D
  use scale_fem_localmesh_3d, only: &
     FEM_LocalMesh3D => LocalMesh3D
  use scale_fem_meshfield_base, only: &
     FEM_MeshField3D => MeshField3D
     
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: MKINIT_common_BUBBLE_setup
  public :: MKINIT_common_BUBBLE_setup_DG
  public :: MKINIT_common_RECT_setup

  public :: MKINIT_common_flux_setup
  public :: MKINIT_common_land_setup
  public :: MKINIT_common_ocean_setup
  public :: MKINIT_common_urban_setup

  public :: MKINIT_common_read_sounding
  public :: MKINIT_common_read_sounding_mpace

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public, parameter           :: THETAstd = 300.0_RP ! [K]

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
  !> Bubble
  subroutine MKINIT_common_BUBBLE_setup( bubble )
    use scale_const, only: &
       CONST_UNDEF8
    implicit none
    real(RP), intent(out) :: bubble(KA,IA,JA)

    ! Bubble
    logical  :: BBL_eachnode = .false.   ! Arrange bubble at each node? [kg/kg]
    real(RP) :: BBL_CZ       =  2.E3_RP  ! center location [m]: z
    real(RP) :: BBL_CX       =  2.E3_RP  ! center location [m]: x
    real(RP) :: BBL_CY       =  2.E3_RP  ! center location [m]: y
    real(RP) :: BBL_RZ       =  0.0_RP   ! bubble radius   [m]: z
    real(RP) :: BBL_RX       =  0.0_RP   ! bubble radius   [m]: x
    real(RP) :: BBL_RY       =  0.0_RP   ! bubble radius   [m]: y
    character(len=H_SHORT) :: BBL_functype = 'COSBELL' ! COSBELL or GAUSSIAN

    namelist / PARAM_BUBBLE / &
       BBL_eachnode, &
       BBL_CZ,       &
       BBL_CX,       &
       BBL_CY,       &
       BBL_RZ,       &
       BBL_RX,       &
       BBL_RY,       &
       BBL_functype

    real(RP) :: CZ_offset
    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: distx, disty, distz

    real(RP) :: Domain_RX, Domain_RY

    logical  :: error
    integer  :: ierr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("BUBBLE_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BUBBLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("BUBBLE_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("BUBBLE_setup",*) 'Not appropriate names in namelist PARAM_BUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_BUBBLE)

    error = .false.

    if ( abs(BBL_RZ*BBL_RX*BBL_RY) <= 0.0_RP ) then
       LOG_INFO("BUBBLE_setup",*) 'no bubble'
       !$acc kernels
       bubble(:,:,:) = 0.0_RP
       !$acc end kernels
    else

       !$acc kernels
       bubble(:,:,:) = CONST_UNDEF8
       !$acc end kernels

       if ( BBL_eachnode ) then
          CZ_offset = CZ(KS)
          CX_offset = CX(IS)
          CY_offset = CY(JS)
          Domain_RX = FX(IE) - FX(IS-1)
          Domain_RY = FY(JE) - FY(JS-1)
       else
          CZ_offset = 0.0_RP
          CX_offset = 0.0_RP
          CY_offset = 0.0_RP
          Domain_RX = FXG(IAG-IHALO) - FXG(IHALO)
          Domain_RY = FYG(JAG-JHALO) - FYG(JHALO)
       endif

       ! make bubble coefficient
       !$acc kernels
       !$acc loop independent collapse(3) reduction(.or.:error)
       do j = 1, JA
       do i = 1, IA
       do k = KS, KE

          distz = ( (CZ(k)-CZ_offset-BBL_CZ)/BBL_RZ )**2

          distx = min( ( (CX(i)-CX_offset-BBL_CX          )/BBL_RX )**2, &
                       ( (CX(i)-CX_offset-BBL_CX-Domain_RX)/BBL_RX )**2, &
                       ( (CX(i)-CX_offset-BBL_CX+Domain_RX)/BBL_RX )**2  )

          disty = min( ( (CY(j)-CY_offset-BBL_CY          )/BBL_RY )**2, &
                       ( (CY(j)-CY_offset-BBL_CY-Domain_RY)/BBL_RY )**2, &
                       ( (CY(j)-CY_offset-BBL_CY+Domain_RY)/BBL_RY )**2  )

          select case(BBL_functype)
          case('COSBELL')
             bubble(k,i,j) = cos( 0.5_RP*PI*sqrt( min(distz+distx+disty,1.0_RP) ) )**2
          case('GAUSSIAN')
             bubble(k,i,j) = exp( -(distz+distx+disty) )
          case default
            LOG_ERROR("BUBBLE_setup",*) 'Not appropriate BBL_functype. Check!', trim(BBL_functype)
#ifdef _OPENACC
            error = .true.
#else
            call PRC_abort
#endif
          end select
       enddo
       enddo
       enddo
       !$acc end kernels
    endif

    if ( error ) call PRC_abort

    return
  end subroutine MKINIT_common_BUBBLE_setup

  !-----------------------------------------------------------------------------
  !> Bubble
  subroutine MKINIT_common_BUBBLE_setup_DG( DG_bubble )
    use scale_const, only: &
       CONST_UNDEF8
    implicit none
    type(FEM_MeshField3D), intent(inout) :: DG_bubble

    ! Bubble
    logical  :: BBL_eachnode = .false.   ! Arrange bubble at each node? [kg/kg]
    real(RP) :: BBL_CZ       =  2.E3_RP  ! center location [m]: z
    real(RP) :: BBL_CX       =  2.E3_RP  ! center location [m]: x
    real(RP) :: BBL_CY       =  2.E3_RP  ! center location [m]: y
    real(RP) :: BBL_RZ       =  0.0_RP   ! bubble radius   [m]: z
    real(RP) :: BBL_RX       =  0.0_RP   ! bubble radius   [m]: x
    real(RP) :: BBL_RY       =  0.0_RP   ! bubble radius   [m]: y
    character(len=H_SHORT) :: BBL_functype = 'COSBELL' ! COSBELL or GAUSSIAN

    namelist / PARAM_BUBBLE / &
       BBL_eachnode, &
       BBL_CZ,       &
       BBL_CX,       &
       BBL_CY,       &
       BBL_RZ,       &
       BBL_RX,       &
       BBL_RY,       &
       BBL_functype

    real(RP) :: CZ_offset
    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP), allocatable :: distx(:), disty(:), distz(:)

    real(RP) :: Domain_RX, Domain_RY

    logical  :: error
    integer  :: ierr
    integer  :: kelem

    integer :: ldomid
    type(FEM_LocalMesh3D), pointer :: lmesh
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("BUBBLE_setup_DG",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_BUBBLE,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("BUBBLE_setup_DG",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("BUBBLE_setup_DG",*) 'Not appropriate names in namelist PARAM_BUBBLE. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_BUBBLE)

    error = .false.

    if ( abs(BBL_RZ*BBL_RX*BBL_RY) <= 0.0_RP ) then
       LOG_INFO("BUBBLE_setup_DG",*) 'no bubble'
       do ldomid=1, fem_mesh3D%LOCAL_MESH_NUM
         DG_bubble%local(ldomid)%val(:,:) = 0.0_RP
       end do       
    else
      do ldomid=1, fem_mesh3D%LOCAL_MESH_NUM
         lmesh => fem_mesh3D%lcmesh_list(ldomid)
         allocate( distx(fem_elem3D%Np), disty(fem_elem3D%Np), distz(fem_elem3D%Np) )

         DG_bubble%local(ldomid)%val(:,:)  = CONST_UNDEF8

         if ( BBL_eachnode ) then
            CZ_offset = CZ(KS)
            CX_offset = CX(IS)
            CY_offset = CY(JS)
            Domain_RX = FX(IE) - FX(IS-1)
            Domain_RY = FY(JE) - FY(JS-1)
         else
            CZ_offset = 0.0_RP
            CX_offset = 0.0_RP
            CY_offset = 0.0_RP
            Domain_RX = FXG(IAG-IHALO) - FXG(IHALO)
            Domain_RY = FYG(JAG-JHALO) - FYG(JHALO)
         endif

         ! make bubble coefficient
         !$omp parallel do private( distx, disty, distz )
         do kelem = lmesh%NeS, lmesh%NeE

            distz(:) = ( (lmesh%zlev(:,kelem)-CZ_offset-BBL_CZ)/BBL_RZ )**2

            distx(:) = min( ( (lmesh%pos_en(:,kelem,1)-CX_offset-BBL_CX          )/BBL_RX )**2, &
                        ( (lmesh%pos_en(:,kelem,1)-CX_offset-BBL_CX-Domain_RX)/BBL_RX )**2, &
                        ( (lmesh%pos_en(:,kelem,1)-CX_offset-BBL_CX+Domain_RX)/BBL_RX )**2  )

            disty(:) = min( ( (lmesh%pos_en(:,kelem,2)-CY_offset-BBL_CY          )/BBL_RY )**2, &
                           ( (lmesh%pos_en(:,kelem,2)-CY_offset-BBL_CY-Domain_RY)/BBL_RY )**2, &
                           ( (lmesh%pos_en(:,kelem,2)-CY_offset-BBL_CY+Domain_RY)/BBL_RY )**2  )

            select case(BBL_functype)
            case('COSBELL')
               DG_bubble%local(ldomid)%val(:,kelem) = cos( 0.5_RP*PI*sqrt( min(distz+distx+disty,1.0_RP) ) )**2
            case('GAUSSIAN')
               DG_bubble%local(ldomid)%val(:,kelem) = exp( -(distz+distx+disty) )
            case default
               LOG_ERROR("BUBBLE_setup_DG",*) 'Not appropriate BBL_functype. Check!', trim(BBL_functype)
#ifdef _OPENACC
               error = .true.
#else
               call PRC_abort
#endif
           end select
         enddo
         deallocate( distx, disty, distz )         
      enddo
     endif

    if ( error ) call PRC_abort

    return
  end subroutine MKINIT_common_BUBBLE_setup_DG

  !-----------------------------------------------------------------------------
  !> Bubble
  subroutine MKINIT_common_RECT_setup( rect )
    use scale_const, only: &
       CONST_UNDEF8
    implicit none
    real(RP), intent(out) :: rect(KA,IA,JA)

    ! Bubble
    logical  :: RCT_eachnode = .false.  ! Arrange rectangle at each node? [kg/kg]
    real(RP) :: RCT_CZ       =  2.E3_RP ! center location [m]: z
    real(RP) :: RCT_CX       =  2.E3_RP ! center location [m]: x
    real(RP) :: RCT_CY       =  2.E3_RP ! center location [m]: y
    real(RP) :: RCT_RZ       =  2.E3_RP ! rectangle z width   [m]: z
    real(RP) :: RCT_RX       =  2.E3_RP ! rectangle x width   [m]: x
    real(RP) :: RCT_RY       =  2.E3_RP ! rectangle y width   [m]: y

    namelist / PARAM_RECT / &
       RCT_eachnode, &
       RCT_CZ,       &
       RCT_CX,       &
       RCT_CY,       &
       RCT_RZ,       &
       RCT_RX,       &
       RCT_RY

    real(RP) :: CZ_offset
    real(RP) :: CX_offset
    real(RP) :: CY_offset
    real(RP) :: dist

    integer  :: ierr
    integer  :: k, i, j
    !---------------------------------------------------------------------------

    LOG_NEWLINE
    LOG_INFO("RECT_setup",*) 'Setup'

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_RECT,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_ERROR("RECT_setup",*) 'Not found namelist. Check!'
       call PRC_abort
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("RECT_setup",*) 'Not appropriate names in namelist PARAM_RECT. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_RECT)

    !$acc kernels
    rect(:,:,:) = CONST_UNDEF8
    !$acc end kernels

    if ( RCT_eachnode ) then
       CZ_offset = CZ(KS)
       CX_offset = CX(IS)
       CY_offset = CY(JS)
    else
       CZ_offset = 0.0_RP
       CX_offset = 0.0_RP
       CY_offset = 0.0_RP
    endif

    !$acc kernels
    do j = 1, JA
    do i = 1, IA
    do k = KS, KE

       ! make tracer rectangle
       dist = 2.0_RP * max( &
            abs(CZ(k) - CZ_offset - RCT_CZ)/RCT_RZ,   &
            abs(CX(i) - CX_offset - RCT_CX)/RCT_RX,   &
            abs(CY(j) - CY_offset - RCT_CY)/RCT_RY    &
            & )
       if ( dist <= 1.0_RP ) then
         rect(k,i,j) = 1.0_RP
       else
         rect(k,i,j) = 0.0_RP
       end if
    enddo
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_common_RECT_setup

  !-----------------------------------------------------------------------------
  !> flux setup
  subroutine MKINIT_common_flux_setup
    use mod_atmos_phy_mp_vars, only: &
       SFLX_rain    => ATMOS_PHY_MP_SFLX_rain, &
       SFLX_snow    => ATMOS_PHY_MP_SFLX_snow
    use mod_atmos_phy_rd_vars, only: &
       SFLX_LW_up   => ATMOS_PHY_RD_SFLX_LW_up, &
       SFLX_LW_dn   => ATMOS_PHY_RD_SFLX_LW_dn, &
       SFLX_SW_up   => ATMOS_PHY_RD_SFLX_SW_up, &
       SFLX_SW_dn   => ATMOS_PHY_RD_SFLX_SW_dn
    implicit none

    ! Flux from Atmosphere
    real(RP) :: FLX_rain   = 0.0_RP ! surface rain flux              [kg/m2/s]
    real(RP) :: FLX_snow   = 0.0_RP ! surface snow flux              [kg/m2/s]
    real(RP) :: FLX_IR_dn  = 0.0_RP ! surface downwad radiation flux [J/m2/s]
    real(RP) :: FLX_NIR_dn = 0.0_RP ! surface downwad radiation flux [J/m2/s]
    real(RP) :: FLX_VIS_dn = 0.0_RP ! surface downwad radiation flux [J/m2/s]

    namelist / PARAM_MKINIT_FLUX / &
       FLX_rain,   &
       FLX_snow,   &
       FLX_IR_dn,  &
       FLX_NIR_dn, &
       FLX_VIS_dn

    integer :: i, j
    integer :: ierr
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_FLUX,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("flux_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("flux_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_FLUX. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_FLUX)

    !$acc kernels
    do j = JSB, JEB
    do i = ISB, IEB
       SFLX_rain (i,j) = FLX_rain
       SFLX_snow (i,j) = FLX_snow

       SFLX_LW_up(i,j) = 0.0_RP
       SFLX_LW_dn(i,j) = FLX_IR_dn
       SFLX_SW_up(i,j) = 0.0_RP
       SFLX_SW_dn(i,j) = FLX_NIR_dn + FLX_VIS_dn
    enddo
    enddo
    !$acc end kernels

    return
  end subroutine MKINIT_common_flux_setup


  !-----------------------------------------------------------------------------
  !> Land setup
  subroutine MKINIT_common_land_setup
    use mod_land_vars, only: &
       LAND_TEMP,       &
       LAND_WATER,      &
       LAND_ICE,        &
       LAND_SFC_TEMP,   &
       LAND_SFC_albedo, &
       SNOW_flag,     &
       SNOW_SFC_TEMP, &
       SNOW_SWE,      &
       SNOW_Depth,    &
       SNOW_Dzero,    &
       SNOW_nosnowsec
    implicit none

    real(RP) :: LND_TEMP                ! land soil temperature      [K]
    real(RP) :: LND_WATER     = 0.15_RP ! land soil moisture         [m3/m3]
    real(RP) :: LND_ICE       = 0.00_RP ! land soil ice              [m3/m3]
    real(RP) :: SFC_TEMP                ! land skin temperature      [K]
    real(RP) :: SFC_albedo_LW = 0.01_RP ! land surface albedo for LW (0-1)
    real(RP) :: SFC_albedo_SW = 0.20_RP ! land surface albedo for SW (0-1)

    namelist / PARAM_MKINIT_LAND / &
       LND_TEMP,      &
       LND_WATER,     &
       LND_ICE,       &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW

    integer  :: ierr
    !---------------------------------------------------------------------------

    LND_TEMP = THETAstd
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_LAND,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("land_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("land_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_LAND. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_LAND)

    !$acc kernels
    LAND_TEMP      (:,:,:)         = LND_TEMP
    LAND_WATER     (:,:,:)         = LND_WATER
    LAND_ICE       (:,:,:)         = LND_ICE
    !$acc end kernels

    !$acc kernels
    LAND_SFC_TEMP  (:,:)           = SFC_TEMP
    !$acc end kernels
    !$acc kernels
    LAND_SFC_albedo(:,:,:,I_R_IR)  = SFC_albedo_LW
    LAND_SFC_albedo(:,:,:,I_R_NIR) = SFC_albedo_SW
    LAND_SFC_albedo(:,:,:,I_R_VIS) = SFC_albedo_SW
    !$acc end kernels

    if ( SNOW_flag ) then
      !!!!! Tentative for snow model !!!!!
       !$acc kernels
       SNOW_SFC_TEMP (:,:) = 273.15_RP
       SNOW_SWE      (:,:) = 0.0_RP
       SNOW_Depth    (:,:) = 0.0_RP
       SNOW_Dzero    (:,:) = 0.0_RP
       SNOW_nosnowsec(:,:) = 0.0_RP
       !$acc end kernels
    end if

    return
  end subroutine MKINIT_common_land_setup

  !-----------------------------------------------------------------------------
  !> Ocean setup
  subroutine MKINIT_common_ocean_setup
    use mod_ocean_vars, only: &
       ICE_flag,         &
       OCEAN_TEMP,       &
       OCEAN_SALT,       &
       OCEAN_UVEL,       &
       OCEAN_VVEL,       &
       OCEAN_OCN_Z0M,    &
       OCEAN_ICE_TEMP,   &
       OCEAN_ICE_MASS,   &
       OCEAN_SFC_TEMP,   &
       OCEAN_SFC_albedo, &
       OCEAN_SFC_Z0M,    &
       OCEAN_SFC_Z0H,    &
       OCEAN_SFC_Z0E
    implicit none

    real(RP) :: OCN_TEMP                 ! ocean temperature                         [K]
    real(RP) :: OCN_SALT      = 0.0_RP   ! ocean salinity                            [psu]
    real(RP) :: OCN_UVEL      = 0.0_RP   ! ocean u-velocity                          [m/s]
    real(RP) :: OCN_VVEL      = 0.0_RP   ! ocean v-velocity                          [m/s]
    real(RP) :: ICE_TEMP                 ! ocean temperature                         [K]
    real(RP) :: ICE_MASS      = 0.0_RP   ! ocean temperature                         [K]
    real(RP) :: SFC_TEMP                 ! ocean skin temperature                    [K]
    real(RP) :: SFC_albedo_LW = 0.04_RP  ! ocean surface albedo for LW               (0-1)
    real(RP) :: SFC_albedo_SW = 0.05_RP  ! ocean surface albedo for SW               (0-1)
    real(RP) :: SFC_Z0M       = 1.E-4_RP ! ocean surface roughness length (momentum) [m]
    real(RP) :: SFC_Z0H       = 1.E-4_RP ! ocean surface roughness length (heat)     [m]
    real(RP) :: SFC_Z0E       = 1.E-4_RP ! ocean surface roughness length (vapor)    [m]

    namelist / PARAM_MKINIT_OCEAN / &
       OCN_TEMP,      &
       OCN_SALT,      &
       OCN_UVEL,      &
       OCN_VVEL,      &
       ICE_TEMP,      &
       ICE_MASS,      &
       SFC_TEMP,      &
       SFC_albedo_LW, &
       SFC_albedo_SW, &
       SFC_Z0M,       &
       SFC_Z0H,       &
       SFC_Z0E

    integer :: ierr
    !---------------------------------------------------------------------------

    OCN_TEMP = THETAstd
    ICE_TEMP = 271.35_RP ! freezing point of the ocean
    SFC_TEMP = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_OCEAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("ocean_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("ocean_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_OCEAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_OCEAN)

    !$acc kernels
    OCEAN_TEMP      (:,:,:) = OCN_TEMP
    OCEAN_SALT      (:,:,:) = OCN_SALT
    OCEAN_UVEL      (:,:,:) = OCN_UVEL
    OCEAN_VVEL      (:,:,:) = OCN_VVEL
    !$acc end kernels
    !$acc kernels
    OCEAN_OCN_Z0M   (:,:)   = SFC_Z0M
    OCEAN_SFC_TEMP  (:,:)           = SFC_TEMP
    !$acc end kernels
    !$acc kernels
    OCEAN_SFC_albedo(:,:,:,I_R_IR)  = SFC_albedo_LW
    OCEAN_SFC_albedo(:,:,:,I_R_NIR) = SFC_albedo_SW
    OCEAN_SFC_albedo(:,:,:,I_R_VIS) = SFC_albedo_SW
    !$acc end kernels
    !$acc kernels
    OCEAN_SFC_Z0M   (:,:)           = SFC_Z0M
    OCEAN_SFC_Z0H   (:,:)           = SFC_Z0H
    OCEAN_SFC_Z0E   (:,:)           = SFC_Z0E
    !$acc end kernels

    if ( ICE_flag ) then
       !$acc kernels
       OCEAN_ICE_TEMP  (:,:)   = ICE_TEMP
       OCEAN_ICE_MASS  (:,:)   = ICE_MASS
       !$acc end kernels
    end if

    return
  end subroutine MKINIT_common_ocean_setup

  !-----------------------------------------------------------------------------
  !> Urban setup
  subroutine MKINIT_common_urban_setup
    use mod_urban_vars, only: &
       URBAN_TR,         &
       URBAN_TB,         &
       URBAN_TG,         &
       URBAN_TC,         &
       URBAN_QC,         &
       URBAN_UC,         &
       URBAN_TRL,        &
       URBAN_TBL,        &
       URBAN_TGL,        &
       URBAN_RAINR,      &
       URBAN_RAINB,      &
       URBAN_RAING,      &
       URBAN_SFC_TEMP,   &
       URBAN_SFC_albedo
    implicit none

    real(RP) :: URB_ROOF_TEMP                 ! Surface temperature of roof           [K]
    real(RP) :: URB_BLDG_TEMP                 ! Surface temperature of building       [K]
    real(RP) :: URB_GRND_TEMP                 ! Surface temperature of ground         [K]
    real(RP) :: URB_CNPY_TEMP                 ! Diagnostic canopy air temperature     [K]
    real(RP) :: URB_CNPY_HMDT       = 0.0_RP  ! Diagnostic canopy humidity            [kg/kg]
    real(RP) :: URB_CNPY_WIND       = 0.0_RP  ! Diagnostic canopy wind                [m/s]
    real(RP) :: URB_ROOF_LAYER_TEMP           ! temperature in layer of roof          [K]
    real(RP) :: URB_BLDG_LAYER_TEMP           ! temperature in layer of building      [K]
    real(RP) :: URB_GRND_LAYER_TEMP           ! temperature in layer of ground        [K]
    real(RP) :: URB_ROOF_RAIN       = 0.0_RP  ! temperature in layer of roof          [kg/m2]
    real(RP) :: URB_BLDG_RAIN       = 0.0_RP  ! temperature in layer of building      [kg/m2]
    real(RP) :: URB_GRND_RAIN       = 0.0_RP  ! temperature in layer of ground        [kg/m2]
    real(RP) :: URB_SFC_TEMP                  ! Grid average of surface temperature   [K]
    real(RP) :: URB_ALB_LW          = 0.10_RP ! Grid average of surface albedo for LW (0-1)
    real(RP) :: URB_ALB_SW          = 0.20_RP ! Grid average of surface albedo for SW (0-1)

    namelist / PARAM_MKINIT_URBAN / &
       URB_ROOF_TEMP,       &
       URB_BLDG_TEMP,       &
       URB_GRND_TEMP,       &
       URB_CNPY_TEMP,       &
       URB_CNPY_HMDT,       &
       URB_CNPY_WIND,       &
       URB_ROOF_LAYER_TEMP, &
       URB_BLDG_LAYER_TEMP, &
       URB_GRND_LAYER_TEMP, &
       URB_ROOF_RAIN,       &
       URB_BLDG_RAIN,       &
       URB_GRND_RAIN,       &
       URB_SFC_TEMP,        &
       URB_ALB_LW,          &
       URB_ALB_SW

    integer :: ierr
    !---------------------------------------------------------------------------

    URB_ROOF_TEMP       = THETAstd
    URB_BLDG_TEMP       = THETAstd
    URB_GRND_TEMP       = THETAstd
    URB_CNPY_TEMP       = THETAstd
    URB_ROOF_LAYER_TEMP = THETAstd
    URB_BLDG_LAYER_TEMP = THETAstd
    URB_GRND_LAYER_TEMP = THETAstd
    URB_SFC_TEMP        = THETAstd

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_URBAN,iostat=ierr)
    if( ierr < 0 ) then !--- missing
       LOG_INFO("urban_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("urban_setup",*) 'Not appropriate names in namelist PARAM_MKINIT_URBAN. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_URBAN)

    !$acc kernels
    URBAN_TRL       (:,:,:)         = URB_ROOF_LAYER_TEMP
    URBAN_TBL       (:,:,:)         = URB_BLDG_LAYER_TEMP
    URBAN_TGL       (:,:,:)         = URB_GRND_LAYER_TEMP
    !$acc end kernels

    !$acc kernels
    URBAN_TR        (:,:)           = URB_ROOF_TEMP
    URBAN_TB        (:,:)           = URB_BLDG_TEMP
    URBAN_TG        (:,:)           = URB_GRND_TEMP
    URBAN_TC        (:,:)           = URB_CNPY_TEMP
    URBAN_QC        (:,:)           = URB_CNPY_HMDT
    URBAN_UC        (:,:)           = URB_CNPY_WIND
    URBAN_RAINR     (:,:)           = URB_ROOF_RAIN
    URBAN_RAINB     (:,:)           = URB_BLDG_RAIN
    URBAN_RAING     (:,:)           = URB_GRND_RAIN
    URBAN_SFC_TEMP  (:,:)           = URB_SFC_TEMP
    !$acc end kernels
    !$acc kernels
    URBAN_SFC_albedo(:,:,:,I_R_IR)  = URB_ALB_LW
    URBAN_SFC_albedo(:,:,:,I_R_NIR) = URB_ALB_SW
    URBAN_SFC_albedo(:,:,:,I_R_VIS) = URB_ALB_SW
    !$acc end kernels

    return
  end subroutine MKINIT_common_urban_setup


  !-----------------------------------------------------------------------------
  !> Read sounding data from file
  subroutine MKINIT_common_read_sounding( &
       DENS, VELX, VELY, POTT, QV )
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP), intent(out) :: DENS(KA)
    real(RP), intent(out) :: VELX(KA)
    real(RP), intent(out) :: VELY(KA)
    real(RP), intent(out) :: POTT(KA)
    real(RP), intent(out) :: QV  (KA)

    real(RP) :: TEMP(KA)
    real(RP) :: PRES(KA)
    real(RP) :: QC  (KA)

    character(len=H_LONG) :: ENV_IN_SOUNDING_file = ''

    integer, parameter :: EXP_klim = 100
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA            ! surface potential temperature [K]
    real(RP) :: SFC_PRES             ! surface pressure [hPa]
    real(RP) :: SFC_QV               ! surface watervapor [g/kg]

    real(RP) :: EXP_z   (EXP_klim+1) ! height      [m]
    real(RP) :: EXP_pott(EXP_klim+1) ! potential temperature [K]
    real(RP) :: EXP_qv  (EXP_klim+1) ! water vapor [g/kg]
    real(RP) :: EXP_u   (EXP_klim+1) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim+1) ! velocity v  [m/s]

    real(RP) :: fact1, fact2
    integer :: k, kref

    integer :: fid
    character(len=H_LONG) :: fname

    logical :: converged

#ifdef _OPENACC
    real(RP) :: work1(KA)
    real(RP) :: work2(KA)
    real(RP) :: work3(KA)
#endif

    integer :: ierr

    namelist / PARAM_MKINIT_SOUNDING / &
       ENV_IN_SOUNDING_file
    
    real(RP) :: pres_sfc, pott_sfc, temp_sfc
    real(RP) :: qv_sfc, qc_sfc
    !---------------------------------------------------------------------------

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SOUNDING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("read_sounding",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("read_sounding",*) 'Not appropriate names in namelist PARAM_MKINIT_SOUNDING. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SOUNDING)

    !--- prepare sounding profile
    fid = IO_get_available_fid()
    call IO_get_fname(fname, ENV_IN_SOUNDING_file)
    LOG_INFO("read_sounding",*) 'Input sounding file:', trim(fname)
    open( fid,                  &
          file   = fname,       &
          form   = 'formatted', &
          status = 'old',       &
          iostat = ierr         )

       if ( ierr /= 0 ) then
          LOG_ERROR("read_sounding",*) '[mod_mkinit/read_sounding] Input file not found!'
       endif

       !--- read sounding file till end
       read(fid,*) SFC_PRES, SFC_THETA, SFC_QV

       LOG_INFO("read_sounding",*) '+ Surface pressure [hPa]',     SFC_PRES
       LOG_INFO("read_sounding",*) '+ Surface pot. temp  [K]',     SFC_THETA
       LOG_INFO("read_sounding",*) '+ Surface water vapor [g/kg]', SFC_QV

       do k = 2, EXP_klim
          read(fid,*,iostat=ierr) EXP_z(k), EXP_pott(k), EXP_qv(k), EXP_u(k), EXP_v(k)
          if ( ierr /= 0 ) exit
       enddo

       EXP_kmax = k - 1
    close(fid)

    ! Boundary
    EXP_z   (1)          = 0.0_RP
    EXP_pott(1)          = SFC_THETA
    EXP_qv  (1)          = SFC_QV
    EXP_u   (1)          = EXP_u   (2)
    EXP_v   (1)          = EXP_v   (2)
    EXP_z   (EXP_kmax+1) = 100.E3_RP
    EXP_pott(EXP_kmax+1) = EXP_pott(EXP_kmax)
    EXP_qv  (EXP_kmax+1) = EXP_qv  (EXP_kmax)
    EXP_u   (EXP_kmax+1) = EXP_u   (EXP_kmax)
    EXP_v   (EXP_kmax+1) = EXP_v   (EXP_kmax)

    do k = 1, EXP_kmax+1
       EXP_qv(k) = EXP_qv(k) * 1.E-3_RP ! [g/kg]->[kg/kg]
    enddo

    ! calc in dry condition
    pres_sfc = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
    pott_sfc = SFC_THETA
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       qv_sfc   = SFC_QV * 1.E-3_RP ! [g/kg]->[kg/kg]
    end if
    qc_sfc = 0.0_RP

    !--- linear interpolate to model grid
    do k = KS, KE
       do kref = 2, EXP_kmax+1
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref  ) ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             pott(k) = EXP_pott(kref-1) * fact1 &
                     + EXP_pott(kref  ) * fact2
             qv  (k) = EXP_qv  (kref-1) * fact1 &
                     + EXP_qv  (kref  ) * fact2
             velx(k) = EXP_u   (kref-1) * fact1 &
                     + EXP_u   (kref  ) * fact2
             vely(k) = EXP_v   (kref-1) * fact1 &
                     + EXP_v   (kref  ) * fact2
          endif
       enddo
    enddo
    if ( ATMOS_HYDROMETEOR_dry ) qv(:) = 0.0_RP

    qc(:) = 0.0_RP

    ! make density & pressure profile in moist condition
    call HYDROSTATIC_buildrho( KA, KS, KE, &
                               pott(:), qv(:), qc(:),               & ! [IN]
                               pres_sfc, pott_sfc, qv_sfc, qc_sfc,  & ! [IN]
                               CZ(:), FZ(:),                        & ! [IN]
#ifdef _OPENACC
                               work1(:), work2(:), work3(:),        & ! [WORK]
#endif
                               DENS(:), temp(:), pres(:), temp_sfc, & ! [OUT]
                               converged                            ) ! [OUT]

    return
  end subroutine MKINIT_common_read_sounding

  !> Read sounding data from file specialized for MPACE project
  subroutine MKINIT_common_read_sounding_mpace( &
       DENS, VELX, VELY, POTT, QV, QCI, QNCI )
    use scale_const, only: &
       Rvap  => CONST_Rvap, &
       PSAT0 => CONST_PSAT0, &
       TEM00 => CONST_TEM00, &
       CL    => CONST_CL,    &
       CPvap => CONST_CPvap, &
       LHV00 => CONST_LHV00
    use scale_atmos_thermodyn, only: &
       ATMOS_THERMODYN_temp_pres2pott
    use scale_atmos_hydrometeor, only: &
       ATMOS_HYDROMETEOR_dry
    implicit none

    real(RP), intent(out) :: DENS(KA)
    real(RP), intent(out) :: VELX(KA)
    real(RP), intent(out) :: VELY(KA)
    real(RP), intent(out) :: POTT(KA)
    real(RP), intent(out) :: QV  (KA)
    real(RP), intent(out) :: QCI (KA)
    real(RP), intent(out) :: QNCI(KA)

    real(RP) :: TEMP(KA)
    real(RP) :: TEMP_OLD(KA)
    real(RP) :: PRES(KA)
    real(RP) :: QC  (KA)
    real(RP) :: POTL(KA)
    real(RP) :: LHV(KA)

    character(len=H_LONG) :: ENV_IN_SOUNDING_file = ''
    logical :: USE_HYDROSTATIC

    integer, parameter :: EXP_klim = 501 ! there was a bug: EXP_klim = 100 (30/3/2020)
    integer            :: EXP_kmax

    real(RP) :: SFC_THETA            ! surface potential temperature [K]
    real(RP) :: SFC_THETAL           ! surface liquid potential temperature [K]
    real(RP) :: SFC_TEMP             ! temperature [K]
    real(RP) :: SFC_PRES             ! surface pressure [hPa]
    real(RP) :: SFC_RH               ! surface relative humidity
    real(RP) :: SFC_QV               ! surface watervapor [g/kg]
    real(RP) :: SFC_QCI              ! surface liquid and ice [g/kg]
    real(RP) :: SFC_QNCI             ! surface liquid and ice conc. [/cm3]

    real(RP) :: EXP_z   (EXP_klim+1) ! height      [m]
    real(RP) :: EXP_pres(EXP_klim+1) ! pressure    [Pa]
    real(RP) :: EXP_temp(EXP_klim+1) ! temperature [K]
    real(RP) :: EXP_pott(EXP_klim+1) ! potential temperature [K]
    real(RP) :: EXP_potl(EXP_klim+1) ! potential temperature [K]
    real(RP) :: EXP_rh  (EXP_klim+1) ! relative humidity
    real(RP) :: EXP_qv  (EXP_klim+1) ! water vapor [g/kg]
    real(RP) :: EXP_qci (EXP_klim+1) ! liquid and ice [g/kg]
    real(RP) :: EXP_qnci(EXP_klim+1) ! liquid and ice conc. [/cm3]
    real(RP) :: EXP_u   (EXP_klim+1) ! velocity u  [m/s]
    real(RP) :: EXP_v   (EXP_klim+1) ! velocity v  [m/s]

    real(RP) :: fact1, fact2, R_gas, CP_gas, qd_assumption
    real(RP) :: CPovR_liq, LovR_liq, RTEM00
    real(RP) :: qdry, Rtot, CPtot

    logical :: converged

    integer :: k, kref
    integer :: fid
    integer :: ierr

    namelist / PARAM_MKINIT_SOUNDING / &
       ENV_IN_SOUNDING_file, &
       USE_HYDROSTATIC

    real(RP) :: pres_sfc, pott_sfc, temp_sfc
    real(RP) :: qv_sfc, qc_sfc
    !---------------------------------------------------------------------------

    ! define necessary constants for evaluating RH
    CPovR_liq = ( CPvap - CL ) / Rvap
    RTEM00 = 1.0_RP / TEM00
    LovR_liq  = LHV00 / Rvap

    !--- read namelist
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_MKINIT_SOUNDING,iostat=ierr)

    if( ierr < 0 ) then !--- missing
       LOG_INFO("read_sounding_mpace",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
       LOG_ERROR("read_sounding_mpace",*) 'Not appropriate names in namelist PARAM_MKINIT_SOUNDING. Check!'
       call PRC_abort
    endif
    LOG_NML(PARAM_MKINIT_SOUNDING)

    !--- prepare sounding profile
    LOG_INFO("read_sounding_mpace",*) 'Input sounding file:', trim(ENV_IN_SOUNDING_file)
    fid = IO_get_available_fid()
    open( fid,                                 &
          file   = trim(ENV_IN_SOUNDING_file), &
          form   = 'formatted',                &
          status = 'old',                      &
          iostat = ierr                        )

    if ( ierr /= 0 ) then
       LOG_ERROR("read_sounding_mpace",*) '[mod_mkinit/read_sounding] Input file not found!'
    endif


    !--- read sounding file till end
    read(fid,*) SFC_PRES, SFC_TEMP, SFC_THETA, SFC_THETAL, SFC_QV, SFC_QCI, SFC_QNCI

    LOG_INFO("read_sounding_mpace",*) '+ Surface pressure [hPa]',                         SFC_PRES
    LOG_INFO("read_sounding_mpace",*) '+ Surface temperature [K]',                        SFC_TEMP
    LOG_INFO("read_sounding_mpace",*) '+ Surface pot. temp  [K]',                         SFC_THETA
    LOG_INFO("read_sounding_mpace",*) '+ Surface water vapor [g/kg]',                     SFC_QV
    LOG_INFO("read_sounding_mpace",*) '+ Surface liquid and ice [g/kg]',                  SFC_QCI
    LOG_INFO("read_sounding_mpace",*) '+ Surface liquid and ice conc. [/cm3]',            SFC_QNCI

    do k = 2, EXP_klim
       read(fid,*,iostat=ierr) EXP_z(k), EXP_pres(k), EXP_temp(k), EXP_pott(k), EXP_potl(k), EXP_qv(k), EXP_u(k), EXP_v(k), EXP_qci(k), EXP_qnci(k)
       if ( ierr /= 0 ) exit
    enddo

    EXP_kmax = k - 1
    close(fid)
       
    ! Boundary
    EXP_z   (1)          = 0.0_RP
    EXP_pres(1)          = SFC_PRES
    EXP_temp(1)          = SFC_TEMP
    EXP_pott(1)          = SFC_THETA
    EXP_potl(1)          = SFC_THETAL
    EXP_qv  (1)          = SFC_QV
    EXP_qci (1)          = SFC_QCI
    EXP_qnci(1)          = SFC_QNCI
    EXP_u   (1)          = EXP_u   (2)
    EXP_v   (1)          = EXP_v   (2)
    EXP_z   (EXP_kmax+1) = 100.E3_RP
    EXP_pres(EXP_kmax+1) = EXP_pres(EXP_kmax)
    EXP_temp(EXP_kmax+1) = EXP_temp(EXP_kmax)
    EXP_pott(EXP_kmax+1) = EXP_pott(EXP_kmax)
    EXP_potl(EXP_kmax+1) = EXP_potl(EXP_kmax)
    EXP_qv  (EXP_kmax+1) = EXP_qv  (EXP_kmax)
    EXP_qci (EXP_kmax+1) = EXP_qci (EXP_kmax)
    EXP_qnci(EXP_kmax+1) = EXP_qnci(EXP_kmax)
    EXP_u   (EXP_kmax+1) = EXP_u   (EXP_kmax)
    EXP_v   (EXP_kmax+1) = EXP_v   (EXP_kmax)
    
    do k = 1, EXP_kmax+1
       EXP_qv(k)   = EXP_qv(k)   * 1.E-3_RP ! [g/kg]->[kg/kg]
       EXP_qci(k)  = EXP_qci(k)  * 1.E-3_RP ! [g/kg]->[kg/kg]
       EXP_pres(k) = EXP_pres(k) * 1.E2_RP  ! [hPa]->[Pa]
    enddo

    ! calc in moist condition
    pres_sfc = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
    pott_sfc = SFC_THETA
    if ( .not. ATMOS_HYDROMETEOR_dry ) then
       qv_sfc = SFC_QV  * 1.E-3_RP ! [g/kg]->[kg/kg]
       qc_sfc = SFC_QCI * 1.E-3_RP ! [g/kg]->[kg/kg]
    end if

    !--- linear interpolate to model grid
    do k = KS, KE
       do kref = 2, EXP_kmax+1
          if (       CZ(k) >  EXP_z(kref-1) &
               .AND. CZ(k) <= EXP_z(kref  ) ) then

             fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
             fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

             PRES(k) = EXP_pres(kref-1) * fact1 &
                     + EXP_pres(kref  ) * fact2
             TEMP(k) = EXP_temp(kref-1) * fact1 &
                     + EXP_temp(kref  ) * fact2
             TEMP_OLD(k) = EXP_temp(kref-1) * fact1 &
                         + EXP_temp(kref  ) * fact2
             POTT(k) = EXP_pott(kref-1) * fact1 &
                     + EXP_pott(kref  ) * fact2
             POTL(k) = EXP_potl(kref-1) * fact1 &
                     + EXP_potl(kref  ) * fact2
             VELX(k) = EXP_u   (kref-1) * fact1 &
                     + EXP_u   (kref  ) * fact2
             VELY(k) = EXP_v   (kref-1) * fact1 &
                     + EXP_v   (kref  ) * fact2
             !QV(k)   = 0.0_RP
             QCI(k)  = 0.0_RP
             QNCI(k)  = 0.0_RP
             QV  (k) = EXP_qv  (kref-1) * fact1 &
                     + EXP_qv  (kref  ) * fact2
             
          endif
       enddo
    enddo

    LOG_INFO("read_sounding_mpace",*) 'Checking initial condition before'
    LOG_INFO("read_sounding_mpace",*) 'Z          D          T          A          P          O          L          V          C'
    do k = KS, KE
       LOG_INFO("read_sounding_mpace",'(9ES15.6)') CZ(k), DENS(k), TEMP(k), TEMP_OLD(k), PRES(k), POTT(k),  POTL(k), QV(k), QCI(k)
    enddo

    if ( ATMOS_HYDROMETEOR_dry ) QV(:) = 0.0_RP

    if ( USE_HYDROSTATIC ) then
       ! make density & pressure profile in dry condition
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  POTL(:), QV(:), QCI(:),              & ! [IN]
                                  pres_sfc, pott_sfc, qv_sfc, qc_sfc,  & ! [IN]
                                  CZ(:), FZ(:),                        & ! [IN]
                                  DENS(:), TEMP(:), PRES(:), temp_sfc, & ! [OUT]
                                  converged                            )

       !--- linear interpolate to model grid
       do k = KS, KE
          do kref = 2, EXP_kmax+1
             if (       CZ(k) >  EXP_z(kref-1) &
                  .AND. CZ(k) <= EXP_z(kref  ) ) then

                fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
                fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

                QCI (k) = EXP_qci (kref-1) * fact1 &
                        + EXP_qci (kref  ) * fact2
                QNCI(k) = EXP_qnci(kref-1) * fact1 &
                        + EXP_qnci(kref  ) * fact2

             endif
          enddo
       enddo

       call HYDROMETEOR_LHV( KA, KS, KE, &
                             TEMP, LHV )

       do k = KS, KE
          TEMP(k) = TEMP(k) + LHV(k) / CPdry * QCI(k)
          qdry = 1.0_RP - QV(k) - QCI(k)
          Rtot = Rdry * qdry + Rvap * QV(k)
          CPtot = CPdry * qdry + CPvap * QV(k) + CL * QCI(k)
          POTT(k) = ( TEMP(k) + LHV(k) / CPdry * QCI(k) ) * ( P00 / PRES(k) )**(Rtot/CPtot)
       enddo

       ! make density & pressure profile in moist condition
       call HYDROSTATIC_buildrho( KA, KS, KE, &
                                  POTT(:), QV(:), QCI(:),              & ! [IN]
                                  pres_sfc, pott_sfc, qv_sfc, qc_sfc,  & ! [IN]
                                  CZ(:), FZ(:),                        & ! [IN]
                                  DENS(:), TEMP(:), PRES(:), temp_sfc, & ! [OUT]
                                  converged                            ) ! [OUT]

       DENS(1:KS-1) = DENS(KS)
       DENS(KE+1:KA) = DENS(KE)
       
    else

       ! calc in moist condition
       pres_sfc = SFC_PRES * 1.E2_RP ! [hPa]->[Pa]
       pott_sfc = SFC_THETA
       if ( .not. ATMOS_HYDROMETEOR_dry ) then
          qv_sfc = SFC_QV  * 1.E-3_RP ! [g/kg]->[kg/kg]
          qc_sfc = SFC_QCI * 1.E-3_RP ! [g/kg]->[kg/kg]
       end if

       !--- linear interpolate to model grid
       do k = KS, KE
          do kref = 2, EXP_kmax+1
             if (       CZ(k) >  EXP_z(kref-1) &
                  .AND. CZ(k) <= EXP_z(kref  ) ) then

                fact1 = ( EXP_z(kref) - CZ(k)   ) / ( EXP_z(kref)-EXP_z(kref-1) )
                fact2 = ( CZ(k) - EXP_z(kref-1) ) / ( EXP_z(kref)-EXP_z(kref-1) )

                QV  (k) = EXP_qv  (kref-1) * fact1 &
                        + EXP_qv  (kref  ) * fact2
                QCI (k) = EXP_qci (kref-1) * fact1 &
                        + EXP_qci (kref  ) * fact2
                QNCI(k) = EXP_qnci(kref-1) * fact1 &
                        + EXP_qnci(kref  ) * fact2
             endif
          enddo
       enddo

       do k = KS, KE
          ! use gas law to deduce density
          ! ASSUMPTION: rho_d = rho
          qd_assumption = 1.0_RP - QV(k) - QCI(k)
          R_gas = qd_assumption * Rdry  + QV(k) * Rvap
          CP_gas = qd_assumption * CPdry + QV(k) * CPvap + QCI(k) * CL
          !R_gas = (1.0_RP - QV(k)) * Rdry  + QV(k) * Rvap
          call ATMOS_THERMODYN_temp_pres2pott(TEMP(k), PRES(k), CP_gas, R_gas, POTT(k))
          DENS(k) = PRES(k) / R_gas / TEMP(k)
          POTT(k) = TEMP(k) * (P00 / PRES(k))**(R_gas/CP_gas)
       enddo
    endif

    LOG_INFO("read_sounding_mpace",*) 'Checking initial condition after'
    LOG_INFO("read_sounding_mpace",*) 'Z          D          T          A          P          O          L          V          C'
    do k = KS, KE
       LOG_INFO("read_sounding_mpace",'(9ES15.6)') CZ(k), DENS(k), TEMP(k), TEMP_OLD(k), PRES(k), POTT(k),  POTL(k), QV(k), QCI(k)
    enddo

    return
  end subroutine MKINIT_common_read_sounding_mpace

end module mod_mkinit_common