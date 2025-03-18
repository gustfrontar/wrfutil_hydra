!-------------------------------------------------------------------------------
!> module Atmosphere / Dynamics / variables
!!
!! @par Description
!!          Managing varaibles with DGM for Atmospheric dynamical process
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
!-------------------------------------------------------------------------------
#include "scalelib.h"
module scale_atmos_dyn_dgm_vars
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

  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: ATMOS_DYN_DGM_vars_setup
  public :: ATMOS_DYN_DGM_vars_finalize

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public, parameter :: PRGVAR_DDENS_ID   = 1
  integer, public, parameter :: PRGVAR_THERM_ID   = 2 ! Variable associated with energy equation (DRHOT or ETOT)
  integer, public, parameter :: PRGVAR_DRHOT_ID   = 2
  integer, public, parameter :: PRGVAR_MOMZ_ID    = 3
  integer, public, parameter :: PRGVAR_MOMX_ID    = 4
  integer, public, parameter :: PRGVAR_MOMY_ID    = 5
  integer, public, parameter :: PRGVAR_SCALAR_NUM = 3
  integer, public, parameter :: PRGVAR_HVEC_NUM   = 1
  integer, public, parameter :: PRGVAR_NUM        = 5  

  integer, public, parameter :: PHYTEND_DENS_ID     = 1
  integer, public, parameter :: PHYTEND_MOMX_ID     = 2
  integer, public, parameter :: PHYTEND_MOMY_ID     = 3
  integer, public, parameter :: PHYTEND_MOMZ_ID     = 4
  integer, public, parameter :: PHYTEND_RHOT_ID     = 5
  integer, public, parameter :: PHYTEND_RHOH_ID     = 6
  integer, public, parameter :: PHYTEND_NUM         = 6
  integer, public :: PHYTEND_NUM_TOT

  !-
  integer, public, parameter :: AUXVAR_PRESHYDRO_ID     = 1
  integer, public, parameter :: AUXVAR_DENSHYDRO_ID     = 2
  integer, public, parameter :: AUXVAR_PRES_ID          = 3
  integer, public, parameter :: AUXVAR_PT_ID            = 4
  integer, public, parameter :: AUXVAR_Rtot_ID          = 5
  integer, public, parameter :: AUXVAR_CVtot_ID         = 6
  integer, public, parameter :: AUXVAR_CPtot_ID         = 7
  integer, public, parameter :: AUXVAR_Qdry_ID          = 8
  integer, public, parameter :: AUXVAR_PRESHYDRO_REF_ID = 9
  integer, public, parameter :: AUXVAR_NUM              = 9

  integer, public, parameter :: AUX_DYNVAR_NUM          = 1

  !-
  integer, public, parameter :: TRCVARS3D_TRCADV_ID   = 1
  integer, public, parameter :: TRCVARS3D_DENS_ID     = 2
  integer, public, parameter :: TRCVARS3D_DENS0_ID    = 3
  integer, public, parameter :: TRCVARS3D_NUM         = 3

  integer, public, parameter :: AUXTRCVARS3D_FCTCOEF_ID  = 1
  integer, public, parameter :: AUXTRCVARS3D_NUM         = 1

  integer, public, parameter :: MASSFLX_Z_ID    = 1  
  integer, public, parameter :: MASSFLX_X_ID    = 2
  integer, public, parameter :: MASSFLX_Y_ID    = 3
  integer, public, parameter :: MASS_FLUX_NUM   = 3

  !--

  type, public :: AtmosDynDGM_vars
    type(MeshField3D) :: prg_vars(PRGVAR_NUM)
    type(MeshFieldContainer) :: prg_vars_comm_list(PRGVAR_NUM)
    type(MeshFieldCommCubeDom3D) :: prg_vars_comm
    integer :: prgvars_hist_id(PRGVAR_NUM)

    type(MeshField3D) :: aux_vars(AUXVAR_NUM)
    type(MeshFieldContainer) :: aux_vars_comm_list(AUXVAR_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_vars_comm
    integer :: auxvars_hist_id(AUXVAR_NUM)

    !-
    type(MeshField3D), allocatable :: qtrc_vars(:)
    type(MeshFieldContainer), allocatable :: qtrc_vars_comm_list(:)
    type(MeshFieldCommCubeDom3D) :: qtrc_vars_comm
    integer, allocatable :: qtrcvars_hist_id(:)

    !-
    type(MeshField3D) :: trc_vars(TRCVARS3D_NUM)
    type(MeshFieldContainer) :: trc_vars_comm_list(1)
    type(MeshFieldCommCubeDom3D) :: trc_vars_comm

    type(MeshField3D) :: aux_trcvars(AUXTRCVARS3D_NUM)
    type(MeshFieldContainer) :: aux_trcvars_comm_list(AUXTRCVARS3D_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_trcvars_comm

    type(MeshField3D) :: aux_trcflux(MASS_FLUX_NUM)
    type(MeshFieldContainer) :: aux_trcflux_comm_list(MASS_FLUX_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_trcflux_comm

    type(MeshField3D) :: alphaDensM, alphaDensP ! coeffcient with stabilizion terms in numerical flux (element boundary data)

    !-
    type(MeshField3D) :: DPRES
    type(MeshFieldContainer) :: aux_dynvars_comm_list(AUX_DYNVAR_NUM)
    type(MeshFieldCommCubeDom3D) :: aux_dynvars_comm

    !-
    type(MeshField3D), allocatable :: phy_tend(:)
    integer, allocatable :: phy_tend_hist_id(:)

    type(MeshField2D) :: Coriolis
    class(MeshCubeDom3D), pointer :: mesh3D

    type(FILE_restart_meshfield_component) :: restart_file
    !-
  contains
    procedure :: Init => ATMOS_DYN_DGM_vars_setup
    procedure :: Final => ATMOS_DYN_DGM_vars_finalize
    procedure :: SetVars_from_FV => ATMOS_DYN_DGM_vars_setvars_from_FV
    procedure :: SetVars_from_FV_lc => ATMOS_DYN_DGM_vars_setvars_from_FV_lc
    procedure :: GetVarsPtr => ATMOS_DYN_DGM_vars_getvars_ptr
    procedure :: GetPhyTendPtr => ATMOS_DYN_DGM_vars_getphytend_ptr
    procedure :: GetTRCVarsPtr => ATMOS_DYN_DGM_trcvars_getvars_ptr
    procedure :: GetTRCAuxVarsPtr => ATMOS_DYN_DGM_trcauxvars_getvars_ptr
    procedure :: GetAuxVarsPtr => ATMOS_DYN_DGM_vars_getauxvars_ptr
    procedure :: Calcuate_pressure => ATMOS_DYN_DGM_vars_calc_DPRES
    procedure :: Calculate_specific_heat => ATMOS_DYN_calc_specific_heat
    procedure :: Calculate_AuxVars => ATMOS_DYN_DGM_vars_calc_axuvars
    procedure :: Fix_negative_value_qtrc => ATMOS_DYN_DGM_vars_negative_fixer_qtrc
    procedure :: Get_vel_fv => ATMOS_DYN_DGM_vars_get_vel_fv
    procedure :: Get_prgvars_fv1 => ATMOS_DYN_DGM_vars_get_prgvar_fv1
    procedure :: Get_prgvars_fv2 => ATMOS_DYN_DGM_vars_get_prgvar_fv2
    generic :: Get_prgvars_fv => Get_prgvars_fv1, Get_prgvars_fv2
    procedure :: Flush_phytend => ATMOS_DYN_DGM_vars_flush_phytend
    procedure :: Set_FVphytend => ATMOS_DYN_DGM_vars_set_phytend_fv
    procedure :: Exchange_prgvars => ATMOS_DYN_DGM_vars_exchange_prgvars
    procedure :: Exchange_massflux => ATMOS_DYN_DGM_vars_exchange_massflux
    procedure :: Exchange_qtrcvars => ATMOS_DYN_DGM_vars_exchange_qtrcvars
    procedure :: Exchange_trcvar => ATMOS_DYN_DGM_vars_exchange_trcvar
    procedure :: Exchange_qtrcauxvars => ATMOS_DYN_DGM_vars_exchange_qtrcauxvars
    procedure :: Exchange_auxvars => ATMOS_DYN_DGM_vars_exchange_auxvars
    procedure :: Setup_hist => ATMOS_DYN_DGM_vars_setup_history
    procedure :: Put_hist => ATMOS_DYN_DGM_vars_put_history
    procedure :: Put_hist_phytend => ATMOS_DYN_DGM_vars_put_history_phytend
    procedure :: Setup_restart => ATMOS_DYN_DGM_vars_setup_restart
    procedure :: Open_restart => ATMOS_DYN_DGM_vars_open_restart
    procedure :: Create_restart => ATMOS_DYN_DGM_vars_create_restart
    procedure :: Read_restart => ATMOS_DYN_DGM_vars_read_restart
    procedure :: Defvar_restart => ATMOS_DYN_DGM_vars_defvar_restart
    procedure :: Enddef_restart => ATMOS_DYN_DGM_vars_enddef_restart
    procedure :: Write_restart => ATMOS_DYN_DGM_vars_write_restart
    procedure :: Close_restart => ATMOS_DYN_DGM_vars_close_restart
  end type AtmosDynDGM_vars

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------

  type(VariableInfo) :: PRGVAR_VARINFO(PRGVAR_NUM)
  DATA PRGVAR_VARINFO / &
    VariableInfo( PRGVAR_DDENS_ID, 'DG_DDENS', 'deviation of density',         &
                    'kg/m3',  3, 'FEM_XYZ',  'air_density'              ),     &
    VariableInfo( PRGVAR_THERM_ID, 'DG_THERM', 'THERM',                        &
                    'kg/m3.K', 3, 'FEM_XYZ',  ''                            ), &
    VariableInfo( PRGVAR_MOMZ_ID , 'DG_MOMZ', 'momentum z',                    &
                    'kg/m2/s', 3, 'FEM_XYZ', 'northward_mass_flux_of_air'   ), &
    VariableInfo( PRGVAR_MOMX_ID , 'DG_MOMX', 'momentum x',                    &
                    'kg/m2/s', 3, 'FEM_XYZ', 'upward_mass_flux_of_air'      ), &
    VariableInfo( PRGVAR_MOMY_ID , 'DG_MOMY', 'momentum y',                    &
                    'kg/m2/s', 3, 'FEM_XYZ', 'eastward_mass_flux_of_air'    )  /    

  type(VariableInfo) :: AUXVAR_VARINFO(AUXVAR_NUM)  
  DATA AUXVAR_VARINFO / &
    VariableInfo( AUXVAR_PRESHYDRO_ID, 'DG_PRES_hyd', 'hydrostatic part of pressure',          &
                    'Pa', 3, 'FEM_XYZ', ''                                                  ), &
    VariableInfo( AUXVAR_DENSHYDRO_ID, 'DG_DENS_hyd', 'hydrostatic part of density',           &
                  'kg/m3', 3, 'FEM_XYZ', ''                                                 ), &
    VariableInfo( AUXVAR_PRES_ID     ,     'DG_PRES', 'pressure',                              &
                    'Pa', 3, 'FEM_XYZ', 'air_pressure'                                      ), &
    VariableInfo( AUXVAR_PT_ID       ,       'DG_PT', 'potential temperature',                 &
                      'K', 3, 'FEM_XYZ', 'potential_temperature'                            ), &
    VariableInfo( AUXVAR_Rtot_ID     ,     'DG_RTOT', 'Total gas constant',                    &
                      'J/kg/K', 3, 'FEM_XYZ', ''                                            ), &
    VariableInfo( AUXVAR_CVtot_ID    ,    'DG_CVTOT', 'Total heat capacity',                   &
                      'J/kg/K', 3, 'FEM_XYZ', ''                                            ), &
    VariableInfo( AUXVAR_CPtot_ID    ,    'DG_CPTOT', 'Total heat capacity',                   &
                      'J/kg/K', 3, 'FEM_XYZ', ''                                            ), &
    VariableInfo( AUXVAR_QDRY_ID     ,     'DG_QDRY', 'dry air',                               &
                      'kg/kg', 3, 'FEM_XYZ', ''                                             ), &
    VariableInfo( AUXVAR_PRESHYDRO_REF_ID, 'DG_PRES_hyd_REF', 'hydrostatic reference pressure',&
                    'Pa', 3, 'FEM_XYZ', ''                                                    )/

  type(VariableInfo), public :: ATMOS_DYN_TRCVARS3D_VINFO(TRCVARS3D_NUM)
  DATA ATMOS_DYN_TRCVARS3D_VINFO / &
    VariableInfo( TRCVARS3D_TRCADV_ID, 'TRCADV', '',     '1',  3, 'XYZ',  '' ), &
    VariableInfo( TRCVARS3D_DENS_ID  ,   'DENS', '', 'kg/m3',  3, 'XYZ',  '' ), & 
    VariableInfo( TRCVARS3D_DENS0_ID  ,  'DENS0', '', 'kg/m3',  3, 'XYZ',  '' ) / 

  type(VariableInfo), public :: ATMOS_DYN_MASS_FLUX_VINFO(MASS_FLUX_NUM)
  DATA ATMOS_DYN_MASS_FLUX_VINFO / &
    VariableInfo( MASSFLX_Z_ID, 'MASSFLX_Z', 'flux in z-direction',  &
                  'kg/s/m2',  3, 'FEM_XYZ',  ''                     ),   &
    VariableInfo( MASSFLX_X_ID, 'MASSFLX_X', 'flux in x-direction',  &
                  'kg/s/m2',  3, 'FEM_XYZ',  ''                     ),   & 
    VariableInfo( MASSFLX_Y_ID, 'MASSFLX_Y', 'flux in y-direction',  &
                  'kg/s/m2',  3, 'FEM_XYZ',  ''                     )    /     

  type(VariableInfo), public :: ATMOS_DYN_AUXTRCVARS3D_VINFO(AUXTRCVARS3D_NUM)
  DATA ATMOS_DYN_AUXTRCVARS3D_VINFO / &
    VariableInfo( AUXTRCVARS3D_FCTCOEF_ID, 'TRCADV_FCTCOEF', '',  &
                  '1',  3, 'FEM_XYZ',  ''                          )  /   
                  
    type(VariableInfo) :: ATMOS_DYN_PHYTEND_VARINFO(PHYTEND_NUM)
    DATA ATMOS_DYN_PHYTEND_VARINFO / &
      VariableInfo( PHYTEND_DENS_ID, 'DG_DENS_tp', 'DENS_tp',                        &
                    'kg/m3/s',  3, 'FEM_XYZ',  'tendency of physical process for DENS' ),   &
      VariableInfo( PHYTEND_MOMX_ID, 'DG_MOMX_tp', 'MOMX_tp',                        &
                    'kg/m2/s',  3, 'FEM_XYZ',  'tendency of physical process for MOMX' ),   &
      VariableInfo( PHYTEND_MOMY_ID, 'DG_MOMY_tp', 'MOMY_tp',                        &
                    'kg/m2/s',  3, 'FEM_XYZ',  'tendency of physical process for MOMY' ),   &
      VariableInfo( PHYTEND_MOMZ_ID, 'DG_MOMZ_tp', 'MOMZ_tp',                        &
                    'kg/m2/s',  3, 'FEM_XYZ',  'tendency of physical process for MOMZ' ),   &
      VariableInfo( PHYTEND_RHOT_ID, 'DG_RHOT_tp', 'RHOT_tp',                        &
                    'kg/m3.K/s',  3, 'FEM_XYZ',  'tendency of physical process for RHOT' ), &
      VariableInfo( PHYTEND_RHOH_ID,  'DG_RHOH_p',  'RHOH_p',                        &
                    'kg/m3.J/s',  3, 'FEM_XYZ',  'heating of physical process for THERM' )   /
                  
contains
!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_setup( this, mesh3D )
    use scale_atmos_grid_cartesC, only: &
      DOMAIN_CENTER_Y => ATMOS_GRID_CARTESC_DOMAIN_CENTER_Y    
    use scale_coriolis, only: &
      CORIOLIS_calc_f
    use scale_tracer, only: &
      TRACER_NAME, TRACER_UNIT 
    use scale_atmos_hydrometeor, only: &
      ATMOS_HYDROMETEOR_dry       
    use scale_fem_localmeshfield_base, only: LOCAL_MESHFIELD_TYPE_NODES_FACEVAL
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    type(MeshCubeDom3D), intent(in), target :: mesh3D

    integer :: iv, iq
    integer :: ldomid
    type(LocalMesh2D), pointer :: lmesh2D
    type(LocalMesh3D), pointer :: lmesh3D

    type(VariableInfo) :: qtrc_vinfo_tmp  
    !--------------------------------------------------

    ! Prepare prognostic variables
    do iv=1, PRGVAR_NUM
      call this%prg_vars(iv)%Init( PRGVAR_VARINFO(iv)%NAME, PRGVAR_VARINFO(iv)%UNIT, mesh3D )
      this%prg_vars_comm_list(iv)%field3d => this%prg_vars(iv)
    enddo
    call this%prg_vars_comm%Init( PRGVAR_SCALAR_NUM, PRGVAR_HVEC_NUM, 0, mesh3D )

    ! Prepare tracer variables
    if ( QA > 0 ) then
      allocate( this%qtrc_vars(QA), this%qtrc_vars_comm_list(QA) )
      allocate( this%qtrcvars_hist_id(QA) )
      do iv=1, QA
        call this%qtrc_vars(iv)%Init( 'DG_'//trim(TRACER_NAME(iv)), TRACER_UNIT(iv), mesh3D )
        this%qtrc_vars_comm_list(iv)%field3d => this%qtrc_vars(iv)
      enddo
      call this%qtrc_vars_comm%Init( QA, 0, 0, mesh3D )

      do iv=1, TRCVARS3D_NUM
        call this%trc_vars(iv)%Init( ATMOS_DYN_TRCVARS3D_VINFO(iv)%NAME, ATMOS_DYN_TRCVARS3D_VINFO(iv)%UNIT, mesh3D )
      enddo
      this%trc_vars_comm_list(1)%field3d => this%trc_vars(TRCVARS3D_TRCADV_ID)
      call this%trc_vars_comm%Init( 1, 0, 0, mesh3D )

      do iv=1, MASS_FLUX_NUM
        call this%aux_trcflux(iv)%Init( ATMOS_DYN_MASS_FLUX_VINFO(iv)%NAME, ATMOS_DYN_MASS_FLUX_VINFO(iv)%UNIT, mesh3D )
        this%aux_trcflux_comm_list(iv)%field3d => this%aux_trcflux(iv)
      enddo
      call this%aux_trcflux_comm%Init( 1, 1, 0, mesh3D )

      do iv=1, AUXTRCVARS3D_NUM
        call this%aux_trcvars(iv)%Init( ATMOS_DYN_AUXTRCVARS3D_VINFO(iv)%NAME, ATMOS_DYN_AUXTRCVARS3D_VINFO(iv)%UNIT, mesh3D )
        this%aux_trcvars_comm_list(iv)%field3d => this%aux_trcvars(iv)
      enddo
      call this%aux_trcvars_comm%Init( AUXTRCVARS3D_NUM, 0, 0, mesh3D )
    end if

    ! Prepare auxilirary variables

    this%aux_dynvars_comm_list(1)%field3d => this%DPRES    
    call this%aux_dynvars_comm%Init( 1, 0, 0, mesh3D )

    do iv=1, AUXVAR_NUM
      call this%aux_vars(iv)%Init( AUXVAR_VARINFO(iv)%NAME, AUXVAR_VARINFO(iv)%UNIT, mesh3D )
      this%aux_vars_comm_list(iv)%field3d => this%aux_vars(iv)
    enddo
    call this%aux_vars_comm%Init( AUXVAR_NUM, 0, 0, mesh3D )

    ! Initialize the tendency of physical processes

    PHYTEND_NUM_TOT = PHYTEND_NUM + QA
    allocate( this%phy_tend(PHYTEND_NUM_TOT) )
    allocate( this%phy_tend_hist_id(PHYTEND_NUM_TOT) )

    do iv=1, PHYTEND_NUM
      call this%phy_tend(iv)%Init( ATMOS_DYN_PHYTEND_VARINFO(iv)%NAME, ATMOS_DYN_PHYTEND_VARINFO(iv)%UNIT, mesh3D, fill_value=0.0_RP )
    enddo
    if (QA > 0) then
      call this%alphaDensM%Init( "alphaDensM", "kg/m3.m/s", mesh3D, LOCAL_MESHFIELD_TYPE_NODES_FACEVAL )
      call this%alphaDensP%Init( "alphaDensP", "kg/m3.m/s", mesh3D, LOCAL_MESHFIELD_TYPE_NODES_FACEVAL )

      qtrc_vinfo_tmp%ndims    = 3
      qtrc_vinfo_tmp%dim_type = 'FEM_XYZ'
      qtrc_vinfo_tmp%STDNAME  = ''

      do iq=1, QA
        iv = PHYTEND_NUM + iq
        qtrc_vinfo_tmp%keyID = iv
        qtrc_vinfo_tmp%NAME  = 'DG_'//trim(TRACER_NAME(iq))//'_tp'
        qtrc_vinfo_tmp%DESC  = 'tendency of physical process for '//trim(TRACER_DESC(iq))
        qtrc_vinfo_tmp%UNIT  = trim(TRACER_UNIT(iq))//'/s'

        call this%phy_tend(iv)%Init( qtrc_vinfo_tmp%NAME, qtrc_vinfo_tmp%UNIT, mesh3D, fill_value=0.0_RP )
      enddo
    end if

    !
    call this%DPRES%Init( "DPRES", "Pa", mesh3D )
    call this%Coriolis%Init( "Coriolis", "s-1", mesh3D%mesh2D )

    do ldomid=1, mesh3D%LOCAL_MESH_NUM
      lmesh3D => mesh3D%lcmesh_list(ldomid)
#if defined(__GFORTRAN__) && __GNUC__ < 5
      LOG_ERROR('ATMOS_DYN_DGM_vars_setup',*) 'GCC version should be >= 5! STOP.'
      call PRC_abort
#else
      lmesh2D => lmesh3D%lcmesh2D
#endif

      call CORIOLIS_calc_f( 1, lmesh2D%refElem2D%Np * lmesh2D%Ne,                            & ! [in]
        lmesh2D%pos_en(:,lmesh2D%NeS:lmesh2D%NeE,2), lmesh2D%lat(:,lmesh2D%NeS:lmesh2D%NeE), & ! [in]
        this%Coriolis%local(ldomid)%val(:,lmesh2D%NeS:lmesh2D%NeE) ) ! [out]
      
      this%aux_vars(AUXVAR_CPtot_ID)%local(ldomid)%val(:,:) = CPdry
      this%aux_vars(AUXVAR_CVtot_ID)%local(ldomid)%val(:,:) = CVdry
      this%aux_vars(AUXVAR_Rtot_ID )%local(ldomid)%val(:,:) = Rdry
      this%aux_vars(AUXVAR_PRESHYDRO_REF_ID)%local(ldomid)%val(:,:) = 0.0_RP
    enddo

    this%mesh3D => mesh3D
    return
  end subroutine ATMOS_DYN_DGM_vars_setup

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_setup_history( this )
    use scale_file_fem_history, only: &
      FILE_HISTORY_meshfield_setup, &
      FILE_HISTORY_meshfield_put
    use scale_file_history, only: FILE_HISTORY_reg
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv, iq
    type(VariableInfo) :: qtrc_vinfo_tmp    
    !--------------------------------------------------

    call FILE_HISTORY_meshfield_setup( mesh3D_=this%mesh3D, skip_hist_setup=.true. )

    do iv=1, PRGVAR_NUM
      call FILE_HISTORY_reg( PRGVAR_VARINFO(iv)%NAME, PRGVAR_VARINFO(iv)%DESC, PRGVAR_VARINFO(iv)%unit, this%prgvars_hist_id(iv), &
        dim_type=PRGVAR_VARINFO(iv)%dim_type )    
    enddo

    do iv=1, QA
      call FILE_HISTORY_reg( this%qtrc_vars(iv)%varname, TRACER_DESC(iv), this%qtrc_vars(iv)%unit, this%qtrcvars_hist_id(iv), &
        dim_type='FEM_XYZ' )    
    enddo

    do iv=1, AUXVAR_NUM
      call FILE_HISTORY_reg( AUXVAR_VARINFO(iv)%NAME, AUXVAR_VARINFO(iv)%DESC, AUXVAR_VARINFO(iv)%unit, this%auxvars_hist_id(iv), &
        dim_type=AUXVAR_VARINFO(iv)%dim_type )    
    enddo

    do iv=1, PHYTEND_NUM
      call FILE_HISTORY_reg( ATMOS_DYN_PHYTEND_VARINFO(iv)%NAME, ATMOS_DYN_PHYTEND_VARINFO(iv)%DESC, ATMOS_DYN_PHYTEND_VARINFO(iv)%unit, this%phy_tend_hist_id(iv), &
        dim_type=ATMOS_DYN_PHYTEND_VARINFO(iv)%dim_type )
    enddo
    qtrc_vinfo_tmp%ndims    = 3
    qtrc_vinfo_tmp%dim_type = 'FEM_XYZ'
    qtrc_vinfo_tmp%STDNAME  = ''
    do iq=1, QA
      iv = PHYTEND_NUM + iq
      qtrc_vinfo_tmp%keyID = iv
      qtrc_vinfo_tmp%NAME  = 'DG_'//trim(TRACER_NAME(iq))//'_tp'
      qtrc_vinfo_tmp%DESC  = 'tendency of physical process for '//trim(TRACER_DESC(iq))
      qtrc_vinfo_tmp%UNIT  = trim(TRACER_UNIT(iq))//'/s'
      call FILE_HISTORY_reg( qtrc_vinfo_tmp%NAME, qtrc_vinfo_tmp%DESC, qtrc_vinfo_tmp%unit, this%phy_tend_hist_id(iv), &
        dim_type=qtrc_vinfo_tmp%dim_type )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_setup_history

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_put_history( this )
    use scale_file_fem_history, only: &
      FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv
    !--------------------------------------------------

    do iv=1, PRGVAR_NUM
      call FILE_HISTORY_meshfield_put( this%prgvars_hist_id(iv), this%prg_vars(iv) )  
    enddo

    do iv=1, QA
      call FILE_HISTORY_meshfield_put( this%qtrcvars_hist_id(iv), this%qtrc_vars(iv) )  
    enddo

    do iv=1, AUXVAR_NUM
      call FILE_HISTORY_meshfield_put( this%auxvars_hist_id(iv), this%aux_vars(iv) )  
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_put_history

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_put_history_phytend( this )
    use scale_file_fem_history, only: &
      FILE_HISTORY_meshfield_put
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv
    !--------------------------------------------------

    do iv=1, PHYTEND_NUM_TOT
      call FILE_HISTORY_meshfield_put( this%phy_tend_hist_id(iv), this%phy_tend(iv) )  
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_put_history_phytend

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_setup_restart( this, & 
    in_basename, in_postfix_timelabel,                &
    out_basename, out_postfix_timelabel,              &
    out_dtype, out_title                              )

    use scale_file_fem_restart, only: &
      FILE_restart_meshfield_setup
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    character(*), intent(in) :: in_basename
    logical, intent(in) :: in_postfix_timelabel
    character(*), intent(in) :: out_basename
    logical, intent(in) :: out_postfix_timelabel
    character(*), intent(in) :: out_title
    character(*), intent(in) :: out_dtype
    !--------------------------------------------------

    call FILE_restart_meshfield_setup

    call this%restart_file%Init( 'ATM_DYN_DGM', &
      in_basename, in_postfix_timelabel,                         &
      out_basename, out_postfix_timelabel, out_dtype, out_title, &
      PRGVAR_NUM + 2 + QA, mesh3D=this%mesh3D)
    
    return
  end subroutine ATMOS_DYN_DGM_vars_setup_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_open_restart( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    !--------------------------------------------------
    call this%restart_file%Open()
    return
  end subroutine ATMOS_DYN_DGM_vars_open_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_read_restart( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv
    !--------------------------------------------------

    do iv=1, PRGVAR_NUM
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%prg_vars(iv)%varname, this%prg_vars(iv) )
    enddo
    do iv=1, QA
      call this%restart_file%Read_Var( DIMTYPE_XYZ, this%qtrc_vars(iv)%varname, this%qtrc_vars(iv) )
    enddo
    do iv=1, AUXVAR_DENSHYDRO_ID
      call this%restart_file%Read_var( DIMTYPE_XYZ, this%aux_vars(iv)%varname, this%aux_vars(iv) )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_read_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_create_restart( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    !--------------------------------------------------
    call this%restart_file%Create()
    return
  end subroutine ATMOS_DYN_DGM_vars_create_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_defvar_restart( this )
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv, rf_vid
    !--------------------------------------------------

    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Def_var( this%prg_vars(iv), &
        PRGVAR_VARINFO(iv)%DESC, rf_vid, DIMTYPE_XYZ )
    enddo
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Def_var( this%aux_vars(iv), &
        AUXVAR_VARINFO(iv)%DESC, rf_vid, DIMTYPE_XYZ )
    enddo
    do iv=1, QA
      rf_vid = PRGVAR_NUM + 2 + iv
      call this%restart_file%Def_var( this%qtrc_vars(iv), &
        TRACER_DESC(iv), rf_vid, DIMTYPE_XYZ )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_defvar_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_enddef_restart( this )
    class(AtmosDynDGM_vars), intent(inout), target :: this
    !--------------------------------------------------
    call this%restart_file%End_def()
    return
  end subroutine ATMOS_DYN_DGM_vars_enddef_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_write_restart( this )
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv, rf_vid
    !--------------------------------------------------

    do iv=1, PRGVAR_NUM
      rf_vid = iv
      call this%restart_file%Write_var(rf_vid, this%prg_vars(iv) )
    enddo
    do iv=1, AUXVAR_DENSHYDRO_ID
      rf_vid = PRGVAR_NUM + iv
      call this%restart_file%Write_var(rf_vid, this%aux_vars(iv) )
    enddo
    do iv=1, QA
      rf_vid = PRGVAR_NUM + 2 + iv
      call this%restart_file%Write_var(rf_vid, this%qtrc_vars(iv) )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_write_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_close_restart( this )
    class(AtmosDynDGM_vars), intent(inout), target :: this
    !--------------------------------------------------
    call this%restart_file%Close()
    return
  end subroutine ATMOS_DYN_DGM_vars_close_restart

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_setvars_from_FV( this, &
    DENS_fv, VELX_fv, VELY_fv, VELZ_fv, POTT_fv, QTRC_fv, &
    interp_ord )
    use scale_atmos_grid_cartesC_index, only: &
      IA, JA, KA, IS, IE, IHALO, JS, JE, JHALO, KS, KE, KHALO
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    real(RP), intent(in) :: DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: VELX_fv(KA,IA,JA)
    real(RP), intent(in) :: VELY_fv(KA,IA,JA)
    real(RP), intent(in) :: VELZ_fv(KA,IA,JA)
    real(RP), intent(in) :: POTT_fv(KA,IA,JA)
    real(RP), intent(in) :: QTRC_fv(KA,IA,JA,QA)
    integer, intent(in) :: interp_ord

    integer :: domid, iq
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D
    type(LocalMeshFieldBaseList) :: lc_QTRC_dg(QA)
    !--------------------------------------------------

    domid = 1
    lcmesh3D => this%mesh3D%lcmesh_list(domid)
    elem3D => lcmesh3D%refElem3D

    do iq=1, QA
      lc_QTRC_dg(iq)%ptr => this%qtrc_vars(iq)%local(domid)
    end do
    call this%SetVars_from_FV_lc( &
      DENS_fv, VELX_fv, VELY_fv, VELZ_fv, POTT_fv, QTRC_fv,                                                     &
      this%prg_vars(PRGVAR_DDENS_ID)%local(domid)%val, this%prg_vars(PRGVAR_MOMX_ID)%local(domid)%val,          &
      this%prg_vars(PRGVAR_MOMY_ID)%local(domid)%val, this%prg_vars(PRGVAR_MOMZ_ID)%local(domid)%val,           &
      this%prg_vars(PRGVAR_THERM_ID)%local(domid)%val, lc_QTRC_dg,                                              & 
      this%aux_vars(AUXVAR_DENSHYDRO_ID)%local(domid)%val, this%aux_vars(AUXVAR_PRESHYDRO_ID)%local(domid)%val, &
      interp_ord, lcmesh3D, elem3D )

    return
  end subroutine ATMOS_DYN_DGM_vars_setvars_from_FV

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_setvars_from_FV_lc( this, &
    DENS_fv, VELX_fv, VELY_fv, VELZ_fv, POTT_fv, QTRC_fv,   &
    DDENS_dg, MOMX_dg, MOMY_dg, MOMZ_dg, DRHOT_dg, QTRC_dg, &
    DENS_hyd_dg, PRES_hyd_dg, &
    interp_ord, lcmesh3D, elem3D )
    use scale_atmos_grid_cartesC_index, only: &
      IA, JA, KA, IS, IE, IHALO, JS, JE, JHALO, KS, KE, KHALO
    use scale_statistics, only: &
      STATISTICS_horizontal_mean
    use scale_interp_vert, only: &
      INTERP_VERT_xi2z
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    use scale_fem_meshfield_util, only: &
      MeshFieldUtil_interp_FVtoDG_3D
    use scale_atmos_dyn_dgm_hydrostatic, only: &
      hydrostatic_calc_basicstate_constPT, &
      hydrostaic_build_rho_XYZ
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(in) :: DENS_fv(KA,IA,JA)
    real(RP), intent(in) :: VELX_fv(KA,IA,JA)
    real(RP), intent(in) :: VELY_fv(KA,IA,JA)
    real(RP), intent(in) :: VELZ_fv(KA,IA,JA)
    real(RP), intent(in) :: POTT_fv(KA,IA,JA)
    real(RP), intent(in) :: QTRC_fv(KA,IA,JA,QA)
    real(RP), intent(out) :: DDENS_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: MOMX_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: MOMY_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: MOMZ_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: DRHOT_dg(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(inout) :: QTRC_dg(QA)
    real(RP), intent(out) :: DENS_hyd_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(out) :: PRES_hyd_dg(elem3D%Np,lcmesh3D%NeA)
    integer, intent(in) :: interp_ord

    integer :: domid, kelem, kelem2D, ke_z, pz
    integer :: i, j, k, ke_x, ke_y
    integer :: iq

    real(RP) :: DG_VELX(elem3D%Np,lcmesh3D%NeA), DG_VELY(elem3D%Np,lcmesh3D%NeA), DG_VELZ(elem3D%Np,lcmesh3D%NeA)
    real(RP) :: DG_DENS(elem3D%Np,lcmesh3D%NeA), DG_POTT(elem3D%Np,lcmesh3D%NeA)
    real(RP) :: DDENS_hyd(elem3D%Np,lcmesh3D%NeA)
    real(RP) :: POTT_z(elem3D%Np,lcmesh3D%NeZ,lcmesh3D%NeX,lcmesh3D%NeY)
    real(RP) :: Rtot_z(elem3D%Np,lcmesh3D%NeZ,lcmesh3D%NeX,lcmesh3D%NeY)
    real(RP) :: CPtot_ov_CVtot_z(elem3D%Np,lcmesh3D%NeZ,lcmesh3D%NeX,lcmesh3D%NeY)
    real(RP) :: bnd_SFC_pres(lcmesh3D%lcmesh2D%refElem2D%Np,lcmesh3D%lcmesh2D%NeA)
    real(RP) :: RHOT_hyd(elem3D%Np)

    real(RP) :: r_h1(elem3D%Np), r_h2(elem3D%Np)

    real(RP) :: q_tmp(elem3D%Np,QA)
    real(RP) :: Qdry_el(elem3D%Np), Rtot_el(elem3D%Np), CPtot_el(elem3D%Np), CVtot_el(elem3D%Np)
    !--------------------------------------------------

    call MeshFieldUtil_interp_FVtoDG_3D( DG_DENS, &
      DENS_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord )
    call MeshFieldUtil_interp_FVtoDG_3D( DG_VELX, &
      VELX_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord )
    call MeshFieldUtil_interp_FVtoDG_3D( DG_VELY, &
      VELY_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord )
    call MeshFieldUtil_interp_FVtoDG_3D( DG_VELZ, &
      VELZ_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord )
    call MeshFieldUtil_interp_FVtoDG_3D( DG_POTT, &
      POTT_fv, lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
      interp_ord )
    do iq=1, QA
      call MeshFieldUtil_interp_FVtoDG_3D( QTRC_dg(iq)%ptr%val, &
        QTRC_fv(:,:,:,iq), lcmesh3D, elem3D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, KS, KE, KA, KHALO, &
        interp_ord )
    end do    

    !$omp parallel do private( ke_x, ke_y, ke_z, kelem, iq, &
    !$omp q_tmp, Qdry_el, Rtot_el, CVtot_el, CPtot_el) collapse(3)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      do iq = 1, QA
        q_tmp(:,iq) = QTRC_dg(iq)%ptr%val(:,kelem)
      end do
      call ATMOS_THERMODYN_specific_heat( &
        elem3D%Np, 1, elem3D%Np, QA,                                         & ! (in)
        q_tmp(:,:), TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
        Qdry_el(:), Rtot_el(:), CVtot_el(:), CPtot_el(:)                     ) ! (out)

      POTT_z(:,ke_z,ke_x,ke_y) = DG_POTT(:,kelem)
      CPtot_ov_CVtot_z(:,ke_z,ke_x,ke_y) = CPtot_el(:) / CVtot_el(:)
      Rtot_z(:,ke_z,ke_x,ke_y) = Rtot_el(:)
    end do
    end do
    end do

    ! Construct hydrostatic state
    
    ! allocate( FZ(elem3D%Nnode_v,lcmesh3D%NeZ) )
    ! allocate( work(elem3D%Nnode_v*))
    ! do ke_z=1, lcmesh3D%NeZ
    !   FZ(:,ke_z) = lcmesh3D%pos_en(elem3D%Colmask(:,1),1+(ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY,3)
    ! end do
    ! call INTERP_VERT_xi2z( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
    !   FZ(:,:), lcmesh3D%zlev(:,:), DG_POTT(:,lcmesh3D%NeS:lcmesh3D%NeS), & ! [IN]
    !   work(:,:,:)                         ) ! [OUT]
    ! call STATISTICS_horizontal_mean( KA, KS, KE, IA, IS, IE, JA, JS, JE, &
    !   work(:,:,:), area(:,:),  & ! [IN]
    !   ATMOS_REFSTATE1D_pott(:) ) ! [OUT]

    call hydrostatic_calc_basicstate_constPT( DENS_hyd_dg, PRES_hyd_dg,                      &
      300.0_RP, 1E5_RP, lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2), lcmesh3D%zlev(:,:),  &
      lcmesh3D, elem3D )

    !$omp parallel private(kelem)
    !$omp do collapse(2)
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX
      bnd_SFC_pres(:,kelem) = PRES00 * ( DG_DENS(elem3D%Hslice(:,1),kelem) * Rtot_z(elem3D%Hslice(:,1),1,ke_x,ke_y) * POTT_z(elem3D%Hslice(:,1),1,ke_x,ke_y) / PRES00 )**CPtot_ov_CVtot_z(elem3D%Hslice(:,1),1,ke_x,ke_y)
    end do
    end do
    !$omp end parallel

    call hydrostaic_build_rho_XYZ( DDENS_hyd, &
      DENS_hyd_dg, PRES_hyd_dg, POTT_z, Rtot_z, CPtot_ov_CVtot_z, &
      lcmesh3D%pos_en(:,:,1), lcmesh3D%pos_en(:,:,2), lcmesh3D%pos_en(:,:,3), &
      lcmesh3D, elem3D, bnd_SFC_PRES=bnd_SFC_pres)
    
    !$omp parallel do private(kelem, kelem2D, RHOT_hyd, r_h1, r_h2)
    do kelem=lcmesh3D%NeS, lcmesh3D%NeE
      kelem2D = lcmesh3D%EMap3Dto2D(kelem)
      r_h1(:) = sqrt( lcmesh3D%G_ij(elem3D%IndexH2Dto3D,kelem2D,1,1) )
      r_h2(:) = sqrt( lcmesh3D%G_ij(elem3D%IndexH2Dto3D,kelem2D,2,2) )

      DG_DENS(:,kelem) = DENS_hyd_dg(:,kelem) + DDENS_hyd(:,kelem)
      DDENS_dg(:,kelem) = DG_DENS(:,kelem) - DENS_hyd_dg(:,kelem)
      MOMX_dg(:,kelem) = DG_DENS(:,kelem) * r_h1(:) * DG_VELX(:,kelem)
      MOMY_dg(:,kelem) = DG_DENS(:,kelem) * r_h2(:) * DG_VELY(:,kelem)
      MOMZ_dg(:,kelem) = DG_DENS(:,kelem) * DG_VELZ(:,kelem)

      RHOT_hyd(:) = PRES00 / Rdry * ( PRES_hyd_dg(:,kelem) / PRES00 )**(CvDry/CpDry)
      DRHOT_dg(:,kelem) = DG_DENS(:,kelem) * DG_POTT(:,kelem) - RHOT_hyd(:)
    enddo
    
    return
  end subroutine ATMOS_DYN_DGM_vars_setvars_from_FV_lc  

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_finalize( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this

    integer :: iv
    !--------------------------------------------------

    call this%prg_vars_comm%Final()
    do iv=1, PRGVAR_NUM
      call this%prg_vars(iv)%Final()
    enddo

    call this%aux_vars_comm%Final()
    do iv=1, AUXVAR_NUM
      call this%aux_vars(iv)%Final()
    enddo

    do iv=1, PHYTEND_NUM_TOT
      call this%phy_tend(iv)%Final()
    enddo
    deallocate( this%phy_tend_hist_id )

    if (QA > 0) then
      call this%qtrc_vars_comm%Final()
      do iv=1, QA
        call this%qtrc_vars(iv)%Final()
      enddo
      deallocate( this%qtrc_vars, this%qtrc_vars_comm_list )
      deallocate( this%qtrcvars_hist_id )

      call this%trc_vars_comm%Final()
      do iv=1, TRCVARS3D_NUM
        call this%trc_vars(iv)%Final()
      enddo

      call this%aux_trcflux_comm%Final()
      do iv=1, MASS_FLUX_NUM
        call this%aux_trcflux(iv)%Final()
      enddo

      call this%aux_trcvars_comm%Final()
      do iv=1, AUXTRCVARS3D_NUM
        call this%aux_trcvars(iv)%Final()
      enddo
    end if

    call this%aux_dynvars_comm%Final()
    call this%DPRES%Final()

    call this%Coriolis%Final()

    return
  end subroutine ATMOS_DYN_DGM_vars_finalize

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_flush_phytend( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    integer :: iv
    integer :: ldomid
    !----------------------------------------

    !$omp parallel do private(ldomid)
    do iv=1, PHYTEND_NUM_TOT
    do ldomid=1, this%mesh3D%LOCAL_MESH_NUM
      this%phy_tend(iv)%local(ldomid)%val(:,:) = 0.0_RP
    enddo
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_flush_phytend

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_getvars_ptr( this, &
    DDENS, MOMX, MOMY, MOMZ, THERM, &
    DENS_hyd, PRES_hyd, PRES_hyd_ref )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    type(MeshField3D), pointer, intent(out) :: DDENS, MOMX, MOMY, MOMZ, THERM
    type(MeshField3D), pointer, intent(out) :: DENS_hyd, PRES_hyd
    type(MeshField3D), pointer, intent(out), optional :: PRES_hyd_ref
    !--------------------------------------------------

    DDENS => this%prg_vars(PRGVAR_DDENS_ID)
    MOMX => this%prg_vars(PRGVAR_MOMX_ID)
    MOMY => this%prg_vars(PRGVAR_MOMY_ID)
    MOMZ => this%prg_vars(PRGVAR_MOMZ_ID)
    THERM => this%prg_vars(PRGVAR_THERM_ID)

    DENS_hyd => this%aux_vars(AUXVAR_DENSHYDRO_ID)
    PRES_hyd => this%aux_vars(AUXVAR_PRESHYDRO_ID)
    if ( present(PRES_hyd_ref) ) PRES_hyd_ref => this%aux_vars(AUXVAR_PRESHYDRO_REF_ID)

    return
  end subroutine ATMOS_DYN_DGM_vars_getvars_ptr

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_getphytend_ptr( this, &
    DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, THERM_tp, RHOH_tp )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    type(MeshField3D), pointer, intent(out) :: DENS_tp, MOMX_tp, MOMY_tp, MOMZ_tp, THERM_tp, RHOH_tp
    !--------------------------------------------------

    DENS_tp => this%phy_tend(PHYTEND_DENS_ID)
    MOMX_tp => this%phy_tend(PHYTEND_MOMX_ID)
    MOMY_tp => this%phy_tend(PHYTEND_MOMY_ID)
    MOMZ_tp => this%phy_tend(PHYTEND_MOMZ_ID)
    THERM_tp => this%phy_tend(PHYTEND_RHOT_ID)
    RHOH_tp => this%phy_tend(PHYTEND_RHOH_ID)

    return
  end subroutine ATMOS_DYN_DGM_vars_getphytend_ptr

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_trcauxvars_getvars_ptr( this, &
    QTRC_tmp, FCT_coef, DDENS, DDENS0, DENS_hyd,   &
    mflx_x, mflx_y, mflx_z, alphDENS_M, alphDENS_P )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    type(MeshField3D), pointer, intent(out) :: QTRC_tmp
    type(MeshField3D), pointer, intent(out) :: FCT_coef
    type(MeshField3D), pointer, intent(out) :: DDENS, DDENS0, DENS_hyd
    type(MeshField3D), pointer, intent(out) :: mflx_x, mflx_y, mflx_z
    type(MeshField3D), pointer, intent(out) :: alphDENS_M, alphDENS_P
    !--------------------------------------------------

    QTRC_tmp => this%trc_vars(TRCVARS3D_TRCADV_ID)
    FCT_coef => this%aux_trcvars(AUXTRCVARS3D_FCTCOEF_ID)
    DDENS => this%trc_vars(TRCVARS3D_DENS_ID)
    DDENS0 => this%trc_vars(TRCVARS3D_DENS0_ID)
    DENS_hyd => this%aux_vars(AUXVAR_DENSHYDRO_ID)

    mflx_x => this%aux_trcflux(MASSFLX_X_ID)
    mflx_y => this%aux_trcflux(MASSFLX_Y_ID)
    mflx_z => this%aux_trcflux(MASSFLX_Z_ID)

    alphDENS_M => this%alphaDensM
    alphDENS_P => this%alphaDensP

    return
  end subroutine ATMOS_DYN_DGM_trcauxvars_getvars_ptr  

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_trcvars_getvars_ptr( this, iq, &
    QTRC, RHOQ_tp )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    integer, intent(in) :: iq
    type(MeshField3D), pointer, intent(out) :: QTRC
    type(MeshField3D), pointer, intent(out) :: RHOQ_tp 
    !--------------------------------------------------

    QTRC => this%qtrc_vars(iq)
    RHOQ_tp => this%phy_tend(PHYTEND_NUM + iq)
    return
  end subroutine ATMOS_DYN_DGM_trcvars_getvars_ptr  

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_getauxvars_ptr( this, &
    Rtot, CVtot, CPtot, PT, PRES )
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    type(MeshField3D), pointer, intent(out) :: Rtot, CVtot, CPtot
    type(MeshField3D), pointer, intent(out), optional :: PT, PRES
    !--------------------------------------------------

    Rtot => this%aux_vars(AUXVAR_Rtot_ID)
    CVtot => this%aux_vars(AUXVAR_CVtot_ID)
    CPtot => this%aux_vars(AUXVAR_CPtot_ID)
    if ( present(PT) ) PT => this%aux_vars(AUXVAR_PT_ID)
    if ( present(PRES) ) PRES => this%aux_vars(AUXVAR_PRES_ID)

    return
  end subroutine ATMOS_DYN_DGM_vars_getauxvars_ptr

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_calc_DPRES( this )
    class(AtmosDynDGM_vars), intent(inout) :: this

    integer :: ldomid
    class(LocalMesh3D), pointer :: lcmesh3D
    !---------------------------

    do ldomid=1, this%mesh3D%LOCAL_MESH_NUM
      lcmesh3D => this%mesh3D%lcmesh_list(ldomid)
      call calc_DRHOT2PRES( this%aux_vars(AUXVAR_PRES_ID)%local(ldomid)%val, this%DPRES%local(ldomid)%val,                                   &
        this%prg_vars(PRGVAR_THERM_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_PRESHYDRO_ID)%local(ldomid)%val,                                              &
        this%aux_vars(AUXVAR_RTOT_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_CVTOT_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_CPTOT_ID)%local(ldomid)%val, &
        lcmesh3D, lcmesh3D%refElem3D )
    end do

    return
  end subroutine ATMOS_DYN_DGM_vars_calc_DPRES
  
!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_calc_axuvars( this )
    class(AtmosDynDGM_vars), intent(inout) :: this

    integer :: ldomid
    class(LocalMesh3D), pointer :: lcmesh3D
    !---------------------------

    call this%Calculate_specific_heat()

    do ldomid=1, this%mesh3D%LOCAL_MESH_NUM
      lcmesh3D => this%mesh3D%lcmesh_list(ldomid)
      call calc_auxvars( this%aux_vars(AUXVAR_PRES_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_PT_ID)%local(ldomid)%val,  &
        this%prg_vars(PRGVAR_DDENS_ID)%local(ldomid)%val, this%prg_vars(PRGVAR_THERM_ID)%local(ldomid)%val,               &
        this%aux_vars(AUXVAR_DENSHYDRO_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_PRESHYDRO_ID)%local(ldomid)%val,       &
        this%aux_vars(AUXVAR_RTOT_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_CVTOT_ID)%local(ldomid)%val, this%aux_vars(AUXVAR_CPTOT_ID)%local(ldomid)%val, &
        lcmesh3D, lcmesh3D%refElem3D )
    end do

    return
  end subroutine ATMOS_DYN_DGM_vars_calc_axuvars

!OCL SERIAL  
  subroutine ATMOS_DYN_calc_specific_heat( this )
    use scale_const, only: &
      Rdry => CONST_Rdry,      &
      CPdry => CONST_CPdry,    &
      CVdry => CONST_CVdry,    &
      PRES00 => CONST_PRE00
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP    
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this

    class(LocalMesh3D), pointer :: lcmesh3D
    integer :: lcdomid
    integer :: varid
    integer :: kelem
    integer :: iq

    class(ElementBase3D), pointer :: elem3D

    real(RP), allocatable :: q_tmp(:,:)
    !-------------------------------------------------------

    ! Calculate specific heat
    do lcdomid=1, this%aux_vars(1)%mesh%LOCAL_MESH_NUM
      lcmesh3D => this%aux_vars(1)%mesh%lcmesh_list(lcdomid)
      elem3D => lcmesh3D%refElem3D
      allocate( q_tmp(elem3D%Np,QA) )

      !$omp parallel do private(kelem, iq, q_tmp)
      do kelem = lcmesh3D%NeS, lcmesh3D%NeE
        do iq = 1, QA
          q_tmp(:,iq) = this%qtrc_vars(iq)%local(lcdomid)%val(:,kelem)
        end do
        call ATMOS_THERMODYN_specific_heat( &
          elem3D%Np, 1, elem3D%Np, QA,                                         & ! (in)
          q_tmp(:,:), TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
          this%aux_vars(AUXVAR_QDRY_ID )%local(lcdomid)%val(:,kelem),            & ! (out)
          this%aux_vars(AUXVAR_Rtot_ID )%local(lcdomid)%val(:,kelem),            & ! (out)
          this%aux_vars(AUXVAR_CVtot_ID)%local(lcdomid)%val(:,kelem),            & ! (out)
          this%aux_vars(AUXVAR_CPtot_ID)%local(lcdomid)%val(:,kelem)             ) ! (out)
      end do
      deallocate(q_tmp)
    end do
  end subroutine ATMOS_DYN_calc_specific_heat

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_negative_fixer_qtrc( this )
    use scale_atmos_hydrometeor, only: &
      QLA, QIA
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    
    integer :: lcdomid
    integer :: iq

    class(LocalMesh3D), pointer :: lcmesh3D
    type(LocalMeshFieldBaseList) :: lc_qtrc(QA)
    !----------------------------------------------------

    do lcdomid=1, this%mesh3D%LOCAL_MESH_NUM      
      lcmesh3D => this%mesh3D%lcmesh_list(lcdomid)

      do iq=1, QA
        lc_qtrc(iq)%ptr => this%qtrc_vars(iq)%local(lcdomid)
      enddo
      call negative_fixer_qtrc( lc_qtrc, &
        this%prg_vars(PRGVAR_DDENS_ID)%local(lcdomid)%val, this%aux_vars(AUXVAR_PRES_ID)%local(lcdomid)%val,      &
        this%aux_vars(AUXVAR_CVtot_ID)%local(lcdomid)%val, this%aux_vars(AUXVAR_CPtot_ID)%local(lcdomid)%val,     &
        this%aux_vars(AUXVAR_Rtot_ID )%local(lcdomid)%val, this%aux_vars(AUXVAR_DENSHYDRO_ID)%local(lcdomid)%val, &
        this%aux_vars(AUXVAR_PRESHYDRO_ID)%local(lcdomid)%val, &
        lcmesh3D, lcmesh3D%refElem3D, QA, QLA, QIA, this%prg_vars(PRGVAR_DRHOT_ID)%local(lcdomid)%val )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_negative_fixer_qtrc


!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_get_vel_fv( this, U_fv, V_fv, W_fv )
    use scale_atmos_grid_cartesC_index, only: &
      IA, JA, KA
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    real(RP), intent(out) :: U_fv(KA,IA,JA)
    real(RP), intent(out) :: V_fv(KA,IA,JA)
    real(RP), intent(out) :: W_fv(KA,IA,JA)

    integer :: lcdomid
    class(LocalMesh3D), pointer :: lcmesh3D
    !----------------------------------------------------

    do lcdomid=1, this%mesh3D%LOCAL_MESH_NUM      
      lcmesh3D => this%mesh3D%lcmesh_list(lcdomid)
      call cal_vel_fv( U_fv, V_fv, W_fv, &
        this%prg_vars(PRGVAR_DDENS_ID)%local(lcdomid)%val, this%prg_vars(PRGVAR_MOMX_ID)%local(lcdomid)%val, &
        this%prg_vars(PRGVAR_MOMY_ID)%local(lcdomid)%val, this%prg_vars(PRGVAR_MOMZ_ID)%local(lcdomid)%val,  &
        this%aux_vars(AUXVAR_DENSHYDRO_ID)%local(lcdomid)%val, lcmesh3D, lcmesh3D%refElem3D                  )
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_get_vel_fv

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_get_prgvar_fv1( this, DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv, &
    dyn_bnd, FILL_BND, qtrc_flag )
    use scale_atmos_grid_cartesC_index, only: &
      IA, JA, KA
    use scale_atmos_dyn_dgm_bnd, only: &
      AtmDynBnd      
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    real(RP), intent(inout) :: DENS_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_fv(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_fv(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_fv(KA,IA,JA,QA)
    class(AtmDynBnd), intent(in) :: dyn_bnd
    logical, intent(in), optional :: FILL_BND
    logical, intent(in), optional :: qtrc_flag(QA)

    integer :: lcdomid
    class(LocalMesh3D), pointer :: lcmesh3D

    type(LocalMeshFieldBaseList) :: lc_QTRC_fv(QA)
    integer :: iq
    !----------------------------------------------------

    do lcdomid=1, this%mesh3D%LOCAL_MESH_NUM      
      lcmesh3D => this%mesh3D%lcmesh_list(lcdomid)
      do iq=1, QA
        lc_QTRC_fv(iq)%ptr => this%qtrc_vars(iq)%local(lcdomid)
      end do
      call ATMOS_DYN_DGM_vars_get_prgvar_fv2( this, DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv,             & ! (out)
        this%prg_vars(PRGVAR_DDENS_ID)%local(lcdomid)%val, this%prg_vars(PRGVAR_MOMX_ID)%local(lcdomid)%val,          & ! (in)
        this%prg_vars(PRGVAR_MOMY_ID)%local(lcdomid)%val, this%prg_vars(PRGVAR_MOMZ_ID)%local(lcdomid)%val,           & ! (in)
        this%prg_vars(PRGVAR_THERM_ID)%local(lcdomid)%val, lc_QTRC_fv,                                                & ! (in)
        this%aux_vars(AUXVAR_DENSHYDRO_ID)%local(lcdomid)%val, this%aux_vars(AUXVAR_PRESHYDRO_ID)%local(lcdomid)%val, & ! (in)
        lcmesh3D, lcmesh3D%refElem3D, dyn_bnd, FILL_BND, qtrc_flag                                                    ) ! (in)
    enddo

    return
  end subroutine ATMOS_DYN_DGM_vars_get_prgvar_fv1

!OCL SERIAL
  subroutine ATMOS_DYN_DGM_vars_get_prgvar_fv2( this, DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv, &
    lc_DENS_dg, lc_MOMX_dg, lc_MOMY_dg, lc_MOMZ_dg, lc_RHOT_dg, lc_QTRC_dg, lc_DENS_hyd_dg, lc_PRES_hyd_dg, lcmesh3D, elem3D, &
    dyn_bnd, FILL_BND, qtrc_flag )
    use scale_atmos_grid_cartesC_index, only: &
      IA, JA, KA
    use scale_atmos_dyn_dgm_bnd, only: &
      AtmDynBnd
    implicit none
    class(AtmosDynDGM_vars), intent(inout), target :: this
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: DENS_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_fv(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_fv(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_fv(KA,IA,JA,QA)
    real(RP), intent(in) :: lc_DENS_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: lc_MOMX_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: lc_MOMY_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: lc_MOMZ_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: lc_RHOT_dg(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: lc_QTRC_dg(QA)
    real(RP), intent(in) :: lc_DENS_hyd_dg(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: lc_PRES_hyd_dg(elem3D%Np,lcmesh3D%NeA)
    class(AtmDynBnd), intent(in) :: dyn_bnd
    logical, intent(in), optional :: FILL_BND
    logical, intent(in), optional :: qtrc_flag(QA)
    !----------------------------------------------------

    call cal_progvars_fv( DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv,     & ! (inout)
      lc_DENS_dg, lc_MOMX_dg, lc_MOMY_dg, lc_MOMZ_dg, lc_RHOT_dg, lc_QTRC_dg,       & ! (in)
      lc_DENS_hyd_dg, lc_PRES_hyd_dg,                                               & ! (in)
      lcmesh3D, elem3D, dyn_bnd%BND_W, dyn_bnd%BND_E, dyn_bnd%BND_S, dyn_bnd%BND_N, & ! (in)
      FILL_BND, qtrc_flag )

    return
  end subroutine ATMOS_DYN_DGM_vars_get_prgvar_fv2

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_prgvars( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%prg_vars_comm%Put( this%prg_vars_comm_list, 1 )
    call this%prg_vars_comm%Exchange()
    call this%prg_vars_comm%Get( this%prg_vars_comm_list, 1 )

    call this%aux_dynvars_comm%Put( this%aux_dynvars_comm_list, 1 )
    call this%aux_dynvars_comm%Exchange()
    call this%aux_dynvars_comm%Get( this%aux_dynvars_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_prgvars

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_qtrcvars( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%qtrc_vars_comm%Put( this%qtrc_vars_comm_list, 1 )
    call this%qtrc_vars_comm%Exchange()
    call this%qtrc_vars_comm%Get( this%qtrc_vars_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_qtrcvars

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_trcvar( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%trc_vars_comm%Put( this%trc_vars_comm_list, 1 )
    call this%trc_vars_comm%Exchange()
    call this%trc_vars_comm%Get( this%trc_vars_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_trcvar

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_qtrcauxvars( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%aux_trcvars_comm%Put( this%aux_trcvars_comm_list, 1 )
    call this%aux_trcvars_comm%Exchange()
    call this%aux_trcvars_comm%Get( this%aux_trcvars_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_qtrcauxvars

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_massflux( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%aux_trcflux_comm%Put( this%aux_trcflux_comm_list, 1 )
    call this%aux_trcflux_comm%Exchange()
    call this%aux_trcflux_comm%Get( this%aux_trcflux_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_massflux

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_exchange_auxvars( this )
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    !--------------------------------------------------

    call this%aux_vars_comm%Put( this%aux_vars_comm_list, 1 )
    call this%aux_vars_comm%Exchange()
    call this%aux_vars_comm%Get( this%aux_vars_comm_list, 1 )

    return
  end subroutine ATMOS_DYN_DGM_vars_exchange_auxvars

!OCL SERAIL
  subroutine ATMOS_DYN_DGM_vars_set_phytend_fv( this, &
    DENS_tp_fv, RHOU_tp_fv, RHOV_tp_fv, MOMX_tp_fv, MOMY_tp_fv, MOMZ_tp_fv, RHOT_tp_fv, &
    MAPF_fv )
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO, I_XY
    use scale_atmos_dyn_dgm_operator, only: &
      set_fv_phytend => ATMOS_DYN_DGM_operator_set_fv_phytend
    implicit none
    class(AtmosDynDGM_vars), intent(inout) :: this
    real(RP), intent(in) :: DENS_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: RHOU_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: RHOV_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMX_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMY_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: MOMZ_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: RHOT_tp_fv(KA,IA,JA)
    real(RP), intent(in) :: MAPF_fv(IA,JA,2,4)  !< map factor

    integer :: lcdomid
    class(LocalMesh3D), pointer :: lcmesh3D
    class(ElementBase3D), pointer :: elem3D

    real(RP) :: RHOU_tp_fv_tmp(KA,IA,JA)
    real(RP) :: RHOV_tp_fv_tmp(KA,IA,JA)
    real(RP) :: RHOW_tp_fv_tmp(KA,IA,JA)

    integer :: i, j, k
    integer :: kelem, ke2D
    !--------------------------------------------------

    lcdomID=1
    lcmesh3D => this%mesh3D%lcmesh_list(lcdomID)
    elem3D => lcmesh3D%refElem3D

    !$omp parallel do private(i,j,k) collapse(2)
    do j=JS, JE
    do i=IS, IE
    do k=KS, KE
      RHOU_tp_fv_tmp(k,i,j) = ( 0.5_RP * ( MOMX_tp_fv(k,i-1,j) + MOMX_tp_fv(k,i,j) ) + RHOU_tp_fv(k,i,j) ) * MAPF_fv(i,j,1,I_XY)
      RHOV_tp_fv_tmp(k,i,j) = ( 0.5_RP * ( MOMY_tp_fv(k,i,j-1) + MOMY_tp_fv(k,i,j) ) + RHOV_tp_fv(k,i,j) ) * MAPF_fv(i,j,2,I_XY)
      RHOW_tp_fv_tmp(k,i,j) = 0.5_RP * ( MOMZ_tp_fv(k-1,i,j) + MOMZ_tp_fv(k,i,j) )
    enddo
    enddo
    enddo

    call set_fv_phytend( this%phy_tend(PHYTEND_DENS_ID)%local(lcdomid)%val, DENS_tp_fv, lcmesh3D, lcmesh3D%refElem3D )
    call set_fv_phytend( this%phy_tend(PHYTEND_MOMX_ID)%local(lcdomid)%val, RHOU_tp_fv_tmp, lcmesh3D, lcmesh3D%refElem3D )
    call set_fv_phytend( this%phy_tend(PHYTEND_MOMY_ID)%local(lcdomid)%val, RHOV_tp_fv_tmp, lcmesh3D, lcmesh3D%refElem3D )
    call set_fv_phytend( this%phy_tend(PHYTEND_MOMZ_ID)%local(lcdomid)%val, RHOW_tp_fv_tmp, lcmesh3D, lcmesh3D%refElem3D )
    call set_fv_phytend( this%phy_tend(PHYTEND_RHOT_ID)%local(lcdomid)%val, RHOT_tp_fv, lcmesh3D, lcmesh3D%refElem3D )

    return
  end subroutine ATMOS_DYN_DGM_vars_set_phytend_fv

!-- private -------------------

!OCL SERIAL
  subroutine calc_DRHOT2PRES( PRES, DPRES, &
    DRHOT, PRES_hyd, Rtot, CVtot, CPtot,                            &
    lcmesh, elem3D                                                  )

    use scale_const, only: &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(out) :: PRES(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(out) :: DPRES(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA) 

    integer :: ke
    real(RP) :: RHOT(elem3D%Np)
    real(RP) :: rP0
    !---------------------------------------------------------------

    rP0 = 1.0_RP / PRES00
    !$omp parallel do private( RHOT )
    do ke=lcmesh%NeS, lcmesh%NeE
      RHOT(:) = PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CvDry/CpDry) + DRHOT(:,ke)

      PRES(:,ke) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT(:) )**( CPtot(:,ke) / CVtot(:,ke) )
      DPRES(:,ke) = PRES(:,ke) - PRES_hyd(:,ke)
    end do

    return
  end subroutine calc_DRHOT2PRES

!OCL SERIAL
  subroutine calc_auxvars( PRES, PT, & 
    DDENS, DRHOT, DENS_hyd, PRES_hyd, &   
    Rtot, CVtot, CPtot,               &
    lcmesh, elem3D                    )

    use scale_const, only: &
      Rdry => CONST_Rdry,   &
      CPdry => CONST_CPdry, &
      CVdry => CONST_CVdry, &
      PRES00 => CONST_PRE00
        
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(out) :: PRES(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(out) :: PT(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh%NeA) 
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: Rtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CVtot(elem3D%Np,lcmesh%NeA)
    real(RP), intent(in) :: CPtot(elem3D%Np,lcmesh%NeA)
    
    integer :: ke
    real(RP) :: RHOT(elem3D%Np), DENS(elem3D%Np)
    real(RP) :: rP0
    !---------------------------------------------------------------

    rP0 = 1.0_RP / PRES00
    !$omp parallel do private( RHOT, DENS )
    do ke=lcmesh%NeS, lcmesh%NeE
      RHOT(:) = PRES00 / Rdry * ( PRES_hyd(:,ke) / PRES00 )**(CvDry/CpDry) + DRHOT(:,ke)
      DENS(:) = DDENS(:,ke) + DENS_hyd(:,ke)

      PT(:,ke) = RHOT(:) / DENS(:)
      PRES(:,ke) = PRES00 * ( Rtot(:,ke) * rP0 * RHOT(:) )**( CPtot(:,ke) / CVtot(:,ke) )
    end do

    return
  end subroutine calc_auxvars

!OCL SERIAL
  subroutine negative_fixer_qtrc( &
    QTRC, DDENS, PRES,                       &
    CVtot, CPtot, Rtot,                      &
    DENS_hyd, PRES_hyd,                      &
    lmesh, elem, QA, QLA, QIA,               &
    DRHOT                                    )

    use scale_const, only: &
      CVdry => CONST_CVdry,  &
      CPdry => CONST_CPdry,  &
      Rdry => CONST_Rdry,    &
      Grav => CONST_GRAV
    use scale_tracer, only: &
      TRACER_MASS, TRACER_R, TRACER_CV, TRACER_CP
    use scale_atmos_thermodyn, only: &
      ATMOS_THERMODYN_specific_heat
    use scale_prc
    implicit none

    class(LocalMesh3D), intent(in) :: lmesh
    class(ElementBase3D), intent(in) :: elem
    integer, intent(in) :: QA
    type(LocalMeshFieldBaseList), intent(inout) :: QTRC(QA)
    real(RP), intent(inout) :: DDENS(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: PRES(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: CVtot(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: CPtot(elem%Np,lmesh%NeA)
    real(RP), intent(inout) :: Rtot(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: DENS_hyd(elem%Np,lmesh%NeA)
    real(RP), intent(in) :: PRES_hyd(elem%Np,lmesh%NeA)
    integer, intent(in) :: QLA, QIA
    real(RP), intent(inout), optional :: DRHOT(elem%Np,lmesh%NeA)

    integer :: kelem
    integer :: iq

    real(RP) ::  int_w(elem%Np)

    real(RP) :: DENS(elem%Np)
    real(RP) :: DDENS0(elem%Np)

    real(RP) :: TRCMASS0(elem%Np), TRCMASS1(elem%Np,QA)
    real(RP) :: MASS0_elem, MASS1_elem
    real(RP) :: IntEn0_elem, IntEn_elem

    real(RP) :: QTRC_tmp(elem%Np,QA), Qdry(elem%Np)
    real(RP) :: CVtot_old(elem%Np), CPtot_old(elem%Np), Rtot_old(elem%Np)
    real(RP) :: InternalEn(elem%Np), InternalEn0(elem%Np), TEMP(elem%Np)
    real(RP) :: RHOT_hyd(elem%Np)
#ifdef SINGLE
    real(RP), parameter :: TRC_EPS = 1E-32_RP
#else
    real(RP), parameter :: TRC_EPS = 1E-128_RP
#endif
    !------------------------------------------------

    !$omp parallel do private( &
    !$omp kelem, iq, DENS, DDENS0, InternalEn, InternalEn0, TEMP, QTRC_tmp, &
    !$omp TRCMASS0, TRCMASS1, MASS0_elem, MASS1_elem, IntEn0_elem, IntEn_elem, &
    !$omp Qdry, CVtot_old, CPtot_old, Rtot_old,                  &
    !$omp int_w, RHOT_hyd )
    do kelem = lmesh%NeS, lmesh%NeE

      do iq = 1, QA
        QTRC_tmp(:,iq) = QTRC(iq)%ptr%val(:,kelem)
      end do
      call ATMOS_THERMODYN_specific_heat( & 
        elem%Np, 1, elem%Np, QA,                                           & ! (in)
        QTRC_tmp, TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
        Qdry, Rtot_old, CVtot_old, CPtot_old                               ) ! (out)

      DENS(:) = DENS_hyd(:,kelem) + DDENS(:,kelem)
      DDENS0(:) = DDENS(:,kelem)
      
      InternalEn0(:) = CVtot(:,kelem) * PRES(:,kelem) / Rtot(:,kelem)
      TEMP(:) = InternalEn0(:) / ( DENS(:) * CVtot(:,kelem) )
      InternalEn(:) = InternalEn0(:)

      int_w(:) = lmesh%Gsqrt(:,kelem) * lmesh%J(:,kelem) * elem%IntWeight_lgl(:)
      do iq = 1, QA !1 + QLA + QIA  
        TRCMASS0(:) = DENS(:) * QTRC_tmp(:,iq)
        TRCMASS1(:,iq) = max( TRC_EPS, TRCMASS0(:) )

        MASS0_elem = sum( int_w(:) * TRCMASS0(:)    )
        MASS1_elem = sum( int_w(:) * TRCMASS1(:,iq) )
        TRCMASS1(:,iq) = max(MASS0_elem, 0.0E0_RP) / MASS1_elem * TRCMASS1(:,iq)
        
        if ( TRACER_MASS(iq) > 0.0_RP ) then
          DDENS(:,kelem) = DDENS(:,kelem) &
            + TRACER_MASS(iq) * ( TRCMASS1(:,iq) - TRCMASS0(:) )
          InternalEn(:) = InternalEn(:) &
            + TRACER_MASS(iq) * ( TRCMASS1(:,iq) - TRCMASS0(:) ) * TRACER_CV(iq) * TEMP(:)
        endif
      end do

      !--

      DENS(:) = DENS_hyd(:,kelem) + DDENS(:,kelem)
      do iq = 1, QA
        QTRC_tmp(:,iq) = TRCMASS1(:,iq) / DENS(:)
        QTRC(iq)%ptr%val(:,kelem) = QTRC_tmp(:,iq)
      end do

      IntEn0_elem = sum( int_w(:) * InternalEn0(:) )
      IntEn_elem  = sum( int_w(:) * InternalEn (:) )
      InternalEn(:) = IntEn0_elem / IntEn_elem * InternalEn(:)

      call ATMOS_THERMODYN_specific_heat( &
        elem%Np, 1, elem%Np, QA,                                           & ! (in)
        QTRC_tmp, TRACER_MASS(:), TRACER_R(:), TRACER_CV(:), TRACER_CP(:), & ! (in)
        Qdry, Rtot(:,kelem), CVtot(:,kelem), CPtot(:,kelem)                ) ! (out)

      InternalEn(:) = InternalEn(:) - ( DDENS(:,kelem) - DDENS0(:) ) * Grav * lmesh%zlev(:,kelem)
      PRES(:,kelem) = InternalEn(:) * Rtot(:,kelem) / CVtot(:,kelem)

      if ( present(DRHOT) ) then
        DRHOT(:,kelem) = PRES00 / Rtot(:,kelem) * ( PRES(:,kelem) / PRES00 )**( CVtot(:,kelem) / CPtot(:,kelem) ) &
                       - PRES00 / Rdry * ( PRES_hyd(:,kelem) / PRES00 )**( CVdry / CPdry )
      end if
    end do

    return
  end subroutine negative_fixer_qtrc

!OCL SERIAL
  subroutine cal_vel_fv( U_fv, V_fv, W_fv, &
    DDENS, MOMX, MOMY, MOMZ, DENS_hyd,     &
    lcmesh3D, elem3D )
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO
    use scale_comm_cartesC, only: &
      COMM_vars8, &
      COMM_wait      
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(out) :: U_fv(KA,IA,JA)
    real(RP), intent(out) :: V_fv(KA,IA,JA)
    real(RP), intent(out) :: W_fv(KA,IA,JA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh3D%NeA)

    integer :: ke_x, ke_y, ke_z
    integer :: i, j, k
    integer :: kelem, kelem2D
    real(RP) :: IntWeight(elem3D%Np)
    real(RP) :: vol
    real(RP) :: r_dens
    !--------------------------------------------------------------------

    !$omp parallel private( ke_x, ke_y, ke_z, i, j, k, kelem, kelem2D, IntWeight, vol, r_dens ) 

    !$omp do collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      kelem2D = lcmesh3D%EMap3Dto2D(kelem)
      i = IHALO + ke_x
      j = JHALO + ke_y
      k = KHALO + ke_z

      IntWeight(:) = elem3D%IntWeight_lgl(:) * lcmesh3D%J(:,kelem) * lcmesh3D%Gsqrt(:,kelem)
      vol = sum(IntWeight(:))
      IntWeight(:) = IntWeight(:) / vol

      r_dens = 1.0_RP / sum( IntWeight(:) * ( DENS_hyd(:,kelem) + DDENS(:,kelem) ) )
      U_fv(k,i,j) = r_dens * sum( IntWeight(:) * MOMX(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,1,1)) )
      V_fv(k,i,j) = r_dens * sum( IntWeight(:) * MOMY(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,2,2)))
      W_fv(k,i,j) = r_dens * sum( IntWeight(:) * MOMZ(:,kelem) )
    enddo
    enddo
    enddo
    !$omp end do
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS, IE
      do k=1, KS-1
        U_fv(k,i,j) = U_fv(KS,i,j)
        V_fv(k,i,j) = V_fv(KS,i,j)
        W_fv(k,i,j) = W_fv(KS,i,j)
      enddo
      do k=KE+1, KA
        U_fv(k,i,j) = U_fv(KE,i,j)
        V_fv(k,i,j) = V_fv(KE,i,j)
        W_fv(k,i,j) = W_fv(KE,i,j)
      enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel

    call COMM_vars8( W_fv(:,:,:), 1 )
    call COMM_vars8( U_fv(:,:,:), 2 )
    call COMM_vars8( V_fv(:,:,:), 3 )
    call COMM_wait ( W_fv(:,:,:), 1, .false. )
    call COMM_wait ( U_fv(:,:,:), 2, .false. )
    call COMM_wait ( V_fv(:,:,:), 3, .false. )    

    return
  end subroutine cal_vel_fv

!OCL SERIAL
  subroutine cal_progvars_fv( DENS_fv, MOMX_fv, MOMY_fv, MOMZ_fv, RHOT_fv, QTRC_fv, &
    DDENS, MOMX, MOMY, MOMZ, DRHOT, QTRC, DENS_hyd, PRES_hyd, &
    lcmesh3D, elem3D, BND_W, BND_E, BND_S, BND_N, &
    FILL_BND, qtrc_flag )
    use scale_atmos_grid_cartesC_index, only: &
      KS, KE, KA, IS, IE, IA, JS, JE, JA, &
      IHALO, JHALO, KHALO
    use scale_comm_cartesC, only: &
       COMM_vars8, &
       COMM_wait    
    implicit none
    class(LocalMesh3D), intent(in) :: lcmesh3D
    class(ElementBase3D), intent(in) :: elem3D
    real(RP), intent(inout) :: DENS_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMX_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMY_fv(KA,IA,JA)
    real(RP), intent(inout) :: MOMZ_fv(KA,IA,JA)
    real(RP), intent(inout) :: RHOT_fv(KA,IA,JA)
    real(RP), intent(inout) :: QTRC_fv(KA,IA,JA,QA)
    real(RP), intent(in) :: DDENS(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMX(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMY(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: MOMZ(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: DRHOT(elem3D%Np,lcmesh3D%NeA)
    type(LocalMeshFieldBaseList), intent(in) :: QTRC(QA)
    real(RP), intent(in) :: DENS_hyd(elem3D%Np,lcmesh3D%NeA)
    real(RP), intent(in) :: PRES_hyd(elem3D%Np,lcmesh3D%NeA)
    logical,  intent(in)    :: BND_W
    logical,  intent(in)    :: BND_E
    logical,  intent(in)    :: BND_S
    logical,  intent(in)    :: BND_N
    logical, intent(in), optional :: FILL_BND
    logical, intent(in), optional :: qtrc_flag(QA)

    integer :: ke_x, ke_y, ke_z
    integer :: i, j, k
    integer :: iq
    integer :: kelem, kelem2D
    real(RP) :: IntWeight(elem3D%Np)
    real(RP) :: vol(lcmesh3D%Ne)
    real(RP) :: RHOT_hyd(elem3D%Np)
    real(RP) :: r_dens
    real(RP) :: MOMX_c(KA,IA,JA)
    real(RP) :: MOMY_c(KA,IA,JA)
    real(RP) :: MOMZ_c(KA,IA,JA)

    integer :: IS_, IE_, JS_, JE_

    logical :: fill_bnd_
    logical :: qtrc_flag_(QA)
    !--------------------------------------------------------------------

    fill_bnd_ = .false.
    qtrc_flag_(:) = .true.
    if ( present(FILL_BND) ) fill_bnd_ = FILL_BND
    if ( present(qtrc_flag) ) qtrc_flag_(:) = qtrc_flag(:)

    !$omp parallel private( &
    !$omp ke_x, ke_y, ke_z, i, j, k, kelem, kelem2D, iq, IntWeight, RHOT_hyd ) 
    
    !$omp do collapse(2)
    do ke_z=1, lcmesh3D%NeZ
    do ke_y=1, lcmesh3D%NeY
    do ke_x=1, lcmesh3D%NeX
      kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
      kelem2D = lcmesh3D%EMap3Dto2D(kelem)
      i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

      IntWeight(:) = elem3D%IntWeight_lgl(:) * lcmesh3D%J(:,kelem) * lcmesh3D%Gsqrt(:,kelem)
      vol(kelem) = sum(IntWeight(:))
      IntWeight(:) = IntWeight(:) / vol(kelem)

      DENS_fv(k,i,j) = sum( IntWeight(:) * ( DENS_hyd(:,kelem) + DDENS(:,kelem) ) )
      MOMX_c(k,i,j) = sum( IntWeight(:) * MOMX(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,1,1)) )
      MOMY_c(k,i,j) = sum( IntWeight(:) * MOMY(:,kelem) * sqrt(lcmesh3D%GIJ(elem3D%IndexH2Dto3D(:),kelem2D,2,2)) )
      MOMZ_c(k,i,j) = sum( IntWeight(:) * MOMZ(:,kelem) )

      RHOT_hyd(:) = PRES00 / Rdry * ( PRES_hyd(:,kelem) / PRES00 )**(CvDry/CpDry)
      RHOT_fv(k,i,j) = sum( IntWeight(:) * ( RHOT_hyd(:) + DRHOT(:,kelem) ) )
    enddo
    enddo
    enddo

    do iq=1, QA
      if ( qtrc_flag_(iq) ) then
        !$omp do collapse(2)
        do ke_z=1, lcmesh3D%NeZ
        do ke_y=1, lcmesh3D%NeY
        do ke_x=1, lcmesh3D%NeX
          kelem = ke_x + (ke_y-1)*lcmesh3D%NeX + (ke_z-1)*lcmesh3D%NeX*lcmesh3D%NeY
          i = IHALO + ke_x; j = JHALO + ke_y; k = KHALO + ke_z

          IntWeight(:) = elem3D%IntWeight_lgl(:) * lcmesh3D%J(:,kelem) * lcmesh3D%Gsqrt(:,kelem) &
                      / vol(kelem)
          QTRC_fv(k,i,j,iq) = sum( IntWeight(:) * ( DENS_hyd(:,kelem) + DDENS(:,kelem) ) * QTRC(iq)%ptr%val(:,kelem) ) &
                            / DENS_fv(k,i,j)
        enddo
        enddo
        enddo
      end if
    enddo

    !$omp end parallel

    call COMM_vars8( DENS_fv(:,:,:), 1 )
    call COMM_vars8( MOMX_c(:,:,:), 2 )
    call COMM_vars8( MOMY_c(:,:,:), 3 )
    call COMM_vars8( RHOT_fv(:,:,:), 4 )

    call COMM_wait ( DENS_fv(:,:,:), 1, fill_bnd_  )
    call COMM_wait ( MOMX_c(:,:,:), 2, fill_bnd_  )
    call COMM_wait ( MOMY_c(:,:,:), 3, fill_bnd_  )
    call COMM_wait ( RHOT_fv(:,:,:), 4, fill_bnd_ )

    do iq=1, QA
      if ( qtrc_flag_(iq) ) then
        call COMM_vars8( QTRC_fv(:,:,:,iq), iq )
      end if
    enddo

    if ( BND_W ) then
      IS_ = IS
    else
      IS_ = IS-1
    end if
    if ( BND_E ) then
      IE_ = IE-1
    else
      IE_ = IE
    end if
    if ( BND_S ) then
      JS_ = JS
    else
      JS_ = JS-1
    end if
    if ( BND_N ) then
      JE_ = JE-1
    else
      JE_ = JE
    end if

    !$omp parallel private(i,j,k)
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS_, IE_
    do k=KS, KE
      MOMX_fv(k,i,j) = 0.5_RP * ( MOMX_c(k,i,j) + MOMX_c(k,i+1,j) )
    enddo
    enddo
    enddo

    !$omp do collapse(2)
    do j=JS_, JE_
    do i=IS, IE
    do k=KS, KE
      MOMY_fv(k,i,j) = 0.5_RP * ( MOMY_c(k,i,j) + MOMY_c(k,i,j+1) )
    enddo
    enddo
    enddo
    !$omp do collapse(2)
    do j=JS, JE
    do i=IS, IE
      do k=KS, KE-1
        MOMZ_fv(k,i,j) = 0.5_RP * ( MOMZ_c(k,i,j) + MOMZ_c(k+1,i,j) )
      enddo
      MOMZ_fv(1:KS-1,i,j) = 0.0_RP; MOMZ_fv(KE:KA,i,j) = 0.0_RP
    enddo
    enddo 
    !$omp end parallel

    do iq = 1, QA
      if ( qtrc_flag_(iq) ) then
        call COMM_wait ( QTRC_fv(:,:,:,iq), iq, fill_bnd_ )
      end if
    enddo
    
    call COMM_vars8( MOMZ_fv(:,:,:), 1 )
    call COMM_vars8( MOMX_fv(:,:,:), 2 )
    call COMM_vars8( MOMY_fv(:,:,:), 3 )

    call COMM_wait ( MOMZ_fv(:,:,:), 1, fill_bnd_ )
    call COMM_wait ( MOMX_fv(:,:,:), 2, fill_bnd_ )
    call COMM_wait ( MOMY_fv(:,:,:), 3, fill_bnd_ )

    return
  end subroutine cal_progvars_fv

end module scale_atmos_dyn_dgm_vars
