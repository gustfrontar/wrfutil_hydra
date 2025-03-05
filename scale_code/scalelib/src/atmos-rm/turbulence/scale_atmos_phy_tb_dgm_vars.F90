!-------------------------------------------------------------------------------
!> module Atmosphere / Physics / Turbulence
!!
!! @par Description
!!          Managing varaibles with DGM for Atmospheric turbulence process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_phy_tb_dgm_vars
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_abort
  use scale_prof
  use scale_tracer
  use scale_const, only: &
    Rdry => CONST_Rdry,   &
    CPdry => CONST_CPdry, &
    CVdry => CONST_CVdry, &
    PRES00 => CONST_PRE00
  use scale_tracer, only: &
    QA, TRACER_DESC
  use scale_fem_element_base, only: &
    ElementBase3D
  use scale_fem_localmesh_2d, only: &
    LocalMesh2D
  use scale_fem_localmesh_3d, only: &
    LocalMesh3D
  use scale_fem_localmeshfield_base, only: &
    LocalMeshFieldBaseList
  use scale_fem_mesh_base3d, only: &
    MeshBase3D,                              &
    DIMTYPE_XYZ  => MeshBase3D_DIMTYPEID_XYZ
  use scale_fem_mesh_cubedom3d, only: &
    MeshCubeDom3D
  use scale_fem_meshfield_base, only: &
    MeshField3D, MeshField2D
  use scale_comm_fem_meshfield_base, only: &
    MeshFieldContainer
  use scale_comm_fem_meshfield_cubedom3d, only: &
    MeshFieldCommCubeDom3D
  use scale_fem_variableinfo, only: &
    VariableInfo
  use scale_file_fem_restart, only: &
    FILE_restart_meshfield_component

  use scale_atmos_phy_tb_common_dgm, only: &
    ATMOS_PHY_TB_AUX_NUM, ATMOS_PHY_TB_AUX_SCALAR_NUM, ATMOS_PHY_TB_AUX_HVEC_NUM, &
    ATMOS_PHY_TB_AUX_HTENSOR_NUM,                                                    &
    ATMOS_PHY_TB_DIAG_NUM,                                                           &
    ATMOS_PHY_TB_TENDS_NUM1, &
    TB_MOMX_t_VID => ATMOS_PHY_TB_MOMX_t_ID,TB_MOMY_t_VID => ATMOS_PHY_TB_MOMY_t_ID, &
    TB_MOMZ_t_VID => ATMOS_PHY_TB_MOMZ_t_ID, TB_RHOT_t_VID => ATMOS_PHY_TB_RHOT_t_ID    

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DFQ3_ID    = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DFQ1_ID    = 2 
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_DFQ2_ID    = 3
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_SCALAR_NUM = 1 
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_HVEC_NUM   = 1
  integer, public, parameter :: ATMOS_PHY_TB_AUXTRC_NUM        = 3
  type(VariableInfo), public :: ATMOS_PHY_TB_AUXTRC_VINFO(ATMOS_PHY_TB_AUXTRC_NUM)
  DATA ATMOS_PHY_TB_AUXTRC_VINFO / &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DFQ3_ID, 'DFQ3', 'Kh * gradient of QTRC (z)',    &
                  'm2/s',  3, 'XYZ',  ''                                           ),  &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DFQ1_ID, 'DFQ1', 'Kh * gradient of QTRC (x)',    &
                  'm2/s',  3, 'XYZ',  ''                                           ),  &
    VariableInfo( ATMOS_PHY_TB_AUXTRC_DFQ2_ID, 'DFQ2', 'Kh * gradient of QTRC (y)',    &
                  'm2/s',  3, 'XYZ',  ''                                           )   /

  !--

  type, public :: AtmosPhyTbDGM_vars
    type(MeshField3D), allocatable :: tends(:)
    integer :: TENDS_NUM_TOT

    type(MeshField3D) :: auxvars(ATMOS_PHY_TB_AUX_NUM)
    type(MeshFieldContainer) :: aux_vars_comm_list(ATMOS_PHY_TB_AUX_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_vars_comm

    type(MeshField3D) :: diagvars(ATMOS_PHY_TB_DIAG_NUM)

    type(MeshField3D) :: auxtrcvars(ATMOS_PHY_TB_AUXTRC_NUM)
    type(MeshFieldContainer) :: aux_trcvars_comm_list(ATMOS_PHY_TB_AUXTRC_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_trcvars_comm

    type(MeshField3D) :: GRAD_DENS(3)
  contains
    procedure :: Init => ATMOS_PHY_TB_DGM_vars_setup
    procedure :: Final => ATMOS_PHY_TB_DGM_vars_finalize
    procedure :: Exchange_auxvars => ATMOS_PHY_TB_DGM_vars_exchange_auxvars
    procedure :: Exchange_auxtrcvars => ATMOS_PHY_TB_DGM_vars_exchange_auxtrcvars
  end type AtmosPhyTbDGM_vars

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
!OCL SERIAL
  subroutine ATMOS_PHY_TB_DGM_vars_setup( this, mesh3D )
    use scale_tracer, only: &
      TRACER_NAME, TRACER_DESC, TRACER_UNIT    
    use scale_atmos_phy_tb_common_dgm, only: &
      atm_phy_tb_common_dgm_get_varinfo
    implicit none
    class(AtmosPhyTbDGM_vars), intent(inout), target :: this
    type(MeshCubeDom3D), intent(in), target :: mesh3D

    integer :: iv, iq
    integer :: idim
    integer :: ldomid
    type(LocalMesh2D), pointer :: lmesh2D
    type(LocalMesh3D), pointer :: lmesh3D

    type(VariableInfo) :: auxvar_info(ATMOS_PHY_TB_AUX_NUM)
    type(VariableInfo) :: diagvar_info(ATMOS_PHY_TB_DIAG_NUM)
    type(VariableInfo) :: tbtend_info(ATMOS_PHY_TB_TENDS_NUM1)    

    type(VariableInfo) :: qtrc_vinfo_tmp    
    !--------------------------------------------------

    this%TENDS_NUM_TOT = ATMOS_PHY_TB_TENDS_NUM1 + QA

    !- Initialize auxiliary and diagnostic variables

    !- Get variable information
    call atm_phy_tb_common_dgm_get_varinfo( auxvar_info, diagvar_info, tbtend_info ) ! (out)

    !- Initialize an object to manage tendencies with turbulent model

    allocate( this%tends(this%TENDS_NUM_TOT) )

    ! Prepare tendency

    do iv=1, ATMOS_PHY_TB_TENDS_NUM1
      call this%tends(iv)%Init( tbtend_info(iv)%NAME, tbtend_info(iv)%UNIT, mesh3D )
    enddo

    qtrc_vinfo_tmp%ndims    = 3
    qtrc_vinfo_tmp%dim_type = 'XYZ'
    qtrc_vinfo_tmp%STDNAME  = ''

    do iq = 1, QA
      iv = ATMOS_PHY_TB_TENDS_NUM1 + iq 
      qtrc_vinfo_tmp%keyID = iv
      qtrc_vinfo_tmp%NAME  = 'TB_'//trim(TRACER_NAME(iq))//'_t'
      qtrc_vinfo_tmp%DESC  = 'tendency of '//trim(TRACER_DESC(iq))//' in TB process'
      qtrc_vinfo_tmp%UNIT  = trim(TRACER_UNIT(iq))//'/s'

      call this%tends(iv)%Init( qtrc_vinfo_tmp%NAME, qtrc_vinfo_tmp%UNIT, mesh3D )
    end do

    !- Initialize an object to manage auxiliary variables with turbulent model
    do iv=1, ATMOS_PHY_TB_AUX_NUM
      call this%auxvars(iv)%Init( auxvar_info(iv)%NAME, auxvar_info(iv)%UNIT, mesh3D )
      this%aux_vars_comm_list(iv)%field3d => this%auxvars(iv)
    enddo
    call this%aux_vars_comm%Init( ATMOS_PHY_TB_AUX_SCALAR_NUM, ATMOS_PHY_TB_AUX_HVEC_NUM, ATMOS_PHY_TB_AUX_HTENSOR_NUM, mesh3D )

    !- Initialize an object to manage diagnostic variables with turbulent model
    do iv=1, ATMOS_PHY_TB_DIAG_NUM
      call this%diagvars(iv)%Init( diagvar_info(iv)%NAME, diagvar_info(iv)%UNIT, mesh3D )
    enddo

    !- Initialize an object to manage tracer variables with turbulent model

    if ( QA > 0 ) then
      do iv=1, ATMOS_PHY_TB_AUXTRC_NUM
        call this%auxtrcvars(iv)%Init( ATMOS_PHY_TB_AUXTRC_VINFO(iv)%NAME, ATMOS_PHY_TB_AUXTRC_VINFO(iv)%UNIT, mesh3D )
        this%aux_trcvars_comm_list(iv)%field3d => this%auxtrcvars(iv)
      enddo
      call this%aux_trcvars_comm%Init( ATMOS_PHY_TB_AUXTRC_SCALAR_NUM, ATMOS_PHY_TB_AUXTRC_HVEC_NUM, 0, mesh3D )

      do idim=1, 3
        call this%GRAD_DENS(idim)%Init( "Grad_DENS", "kg/m4", mesh3D )  
      enddo
    endif

    return
  end subroutine ATMOS_PHY_TB_DGM_vars_setup

!OCL SERIAL
  subroutine ATMOS_PHY_TB_DGM_vars_finalize( this )
    implicit none
    class(AtmosPhyTbDGM_vars), intent(inout) :: this

    integer :: iv
    integer :: idim
    !--------------------------------------------------

    do iv=1, size(this%tends)
      call this%tends(iv)%Final()
    enddo
    deallocate( this%tends )

    call this%aux_vars_comm%Final()
    do iv=1, ATMOS_PHY_TB_AUX_NUM
      call this%auxvars(iv)%Final()
    enddo

    do iv=1, ATMOS_PHY_TB_DIAG_NUM
      call this%diagvars(iv)%Final()
    enddo

    if ( QA > 0 )then
      do idim=1, 3
        call this%GRAD_DENS(idim)%Final()
      enddo
      call this%aux_trcvars_comm%Final()
      do iv=1, ATMOS_PHY_TB_AUXTRC_NUM
        call this%auxtrcvars(iv)%Final()
      enddo
    end if

    return
  end subroutine ATMOS_PHY_TB_DGM_vars_finalize


!OCL SERAIL
  subroutine ATMOS_PHY_TB_DGM_vars_exchange_auxvars( this )
    implicit none
    class(AtmosPhyTbDGM_vars), intent(inout) :: this

    integer :: iv
    !--------------------------------------------------

    call this%aux_vars_comm%Put( this%aux_vars_comm_list, 1 )
    call this%aux_vars_comm%Exchange()
    call this%aux_vars_comm%Get( this%aux_vars_comm_list, 1 )

    return
  end subroutine ATMOS_PHY_TB_DGM_vars_exchange_auxvars

!OCL SERAIL
  subroutine ATMOS_PHY_TB_DGM_vars_exchange_auxtrcvars( this )
    implicit none
    class(AtmosPhyTbDGM_vars), intent(inout) :: this

    integer :: iv
    !--------------------------------------------------

    call this%aux_trcvars_comm%Put( this%aux_trcvars_comm_list, 1 )
    call this%aux_trcvars_comm%Exchange()
    call this%aux_trcvars_comm%Get( this%aux_trcvars_comm_list, 1 )

    return
  end subroutine ATMOS_PHY_TB_DGM_vars_exchange_auxtrcvars

end module scale_atmos_phy_tb_dgm_vars
