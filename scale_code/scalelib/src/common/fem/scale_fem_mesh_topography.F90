!-------------------------------------------------------------------------------
!> module common / FEM mesh / Topography
!!
!! @par Description
!!         Manage topography data with FEM mesh
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scalelib.h"
module scale_fem_mesh_topography

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_fem_localmesh_3d, only: LocalMesh3D
  use scale_fem_localmesh_2d, only: LocalMesh2D
  use scale_fem_mesh_base2d, only: MeshBase2D
  use scale_fem_mesh_rectdom2d, only: MeshRectDom2D
  use scale_fem_mesh_base3d, only: MeshBase3D
  use scale_fem_element_base, only: ElementBase3D
  use scale_fem_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_comm_fem_meshfield_base, only: &
    MeshFieldCommBase, MeshFieldContainer

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, public :: MeshTopography
    type(MeshField2D) :: topo 
  contains
    procedure :: Init => MeshTopography_Init
    procedure :: Final => MeshTopography_Final
    procedure :: SetVCoordinate => MeshTopography_set_vcoordinate
    procedure :: Read => MeshTopography_read
    procedure :: Write => MeshTopography_write
    procedure :: Set_from_FVtopo => MeshTopography_set_from_FVtopo
  end type MeshTopography

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

contains
!OCL SERIAL
  subroutine MeshTopography_Init( this, varname, mesh )
    implicit none

    class(MeshTopography), intent(inout) :: this
    character(len=*), intent(in) :: varname
    class(MeshBase2D), intent(in), target :: mesh

    integer :: n
    !-----------------------------------------------------------------------------

    call this%topo%Init( varname, "m", mesh )

    do n=1, mesh%LOCAL_MESH_NUM
!$omp parallel workshare
      this%topo%local(n)%val(:,:) = 0.0_RP
!$omp end parallel workshare
    end do

    return
  end subroutine MeshTopography_Init

!OCL SERIAL
  subroutine MeshTopography_Final( this )
    implicit none

    class(MeshTopography), intent(inout) :: this
    !-----------------------------------------------------------------------------
          
    call this%topo%Final()

    return
  end subroutine MeshTopography_Final
  
!OCL SERIAL
  subroutine MeshTopography_set_vcoordinate( this, mesh3D, &
      vcoord_id, zTop, comm3D, comm2D                      )
    
    use scale_sparsemat, only: SparseMat
    use scale_fem_meshutil_vcoord, only: MeshUtil_VCoord_GetMetric
    implicit none

    class(MeshTopography), target, intent(inout) :: this
    class(MeshBase3D), intent(inout), target :: mesh3D
    integer, intent(in) :: vcoord_id
    real(RP), intent(in) :: zTop
    class(MeshFieldCommBase), intent(inout) :: comm3D
    class(MeshFieldCommBase), intent(inout) :: comm2D

    type(SparseMat) :: Dx2D, Dy2D, Lift2D
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase3D), pointer :: elem3D

    integer :: n
    integer :: ke, ke2D

    type(MeshFieldContainer) :: comm2d_varlist(1)
    type(MeshFieldContainer) :: comm3d_varlist(4)

    type(MeshField3D), target :: zlev, GsqrtV, G13, G23
    type(MeshField3D), target :: tmp_G13, tmp_G23

    logical :: flag_covariantvec
    !-------------------------------------------------------------

    lcmesh => mesh3D%lcmesh_list(1)
    lcmesh2D => lcmesh%lcmesh2D
    call Dx2D  %Init( lcmesh2D%refElem2D%Dx1  )
    call Dy2D  %Init( lcmesh2D%refElem2D%Dx2  )
    call Lift2D%Init( lcmesh2D%refElem2D%Lift )

    ! Exchange topography data to fill halo

    comm2d_varlist(1)%field2d => this%topo
    call comm2D%Put(comm2d_varlist, 1)
    call comm2D%Exchange()
    call comm2D%Get(comm2d_varlist, 1)

    ! Calculate metric factors associated with general vertical coordinates

    call zlev%Init( "zlev", "", mesh3D )
    call GsqrtV%Init( "GsqrtV", "", mesh3D )
    call G13%Init( "G13", "", mesh3D )
    call G23%Init( "G23", "", mesh3D )

    call tmp_G13%Init( "tmp_G13", "", mesh3D )
    call tmp_G23%Init( "tmp_G23", "", mesh3D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      lcmesh2D => lcmesh%lcmesh2D
      elem3D => lcmesh%refElem3D

      call MeshUtil_VCoord_GetMetric( &
        G13%local(n)%val, G23%local(n)%val, zlev%local(n)%val,  & ! (out)
        GsqrtV%local(n)%val,                                    & ! (out)
        this%topo%local(n)%val(:,:), zTop, vcoord_id,           & ! (in)
        lcmesh, lcmesh%refElem3D, lcmesh2D, lcmesh2D%refElem2D, & ! (in)
        Dx2D, Dy2D, Lift2D                                      ) ! (in)

      !$omp parallel do private(ke2D)
      do ke=lcmesh%NeS, lcmesh%NeE
        ke2D = lcmesh%EMap3Dto2D(ke)
        tmp_G13%local(n)%val(:,ke) = &
          lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,1,1) * G13%local(n)%val(:,ke) &
        + lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,1,2) * G23%local(n)%val(:,ke)
        tmp_G23%local(n)%val(:,ke) = &
          lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,2,1) * G13%local(n)%val(:,ke) &
        + lcmesh%GIJ(elem3D%IndexH2Dto3D(:),ke2D,2,2) * G23%local(n)%val(:,ke)
      end do        
    end do

    ! Exchange metric data to fill halo

    comm3d_varlist(1)%field3d => zlev
    comm3d_varlist(2)%field3d => GsqrtV
    comm3d_varlist(3)%field3d => tmp_G13
    comm3d_varlist(4)%field3d => tmp_G23

    flag_covariantvec = .false.
    ! select type(comm3D)
    ! type is (MeshFieldCommCubedSphereDom3D)
    !   call comm3D%SetCovariantVec( 1, G13, G23 )
    !   flag_covariantvec = .true.
    ! end select

    call comm3D%Put(comm3d_varlist, 1)
    call comm3D%Exchange()
    call comm3D%Get(comm3d_varlist, 1)

    ! Update metric data managed by local mesh

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      elem3D => lcmesh%refElem3D

      if ( flag_covariantvec ) then
        call mesh3D%Set_geometric_with_vcoord( n, &
          GsqrtV%local(n)%val, zlev%local(n)%val, G13%local(n)%val, G23%local(n)%val )
      else
        call mesh3D%Set_geometric_with_vcoord( n, &
          GsqrtV%local(n)%val, zlev%local(n)%val, tmp_G13%local(n)%val, tmp_G23%local(n)%val )
      end if
    end do

    !---
    call zlev%Final()
    call GsqrtV%Final()
    call G13%Final()
    call G23%Final()
    call tmp_G13%Final()
    call tmp_G23%Final()

    call Dx2D%Final()
    call Dy2D%Final()
    call Lift2D%Final()

    return
  end subroutine MeshTopography_set_vcoordinate

!---
  subroutine MeshTopography_write( this )
    use scale_prc, only: &
      PRC_myrank, PRC_abort
    use scale_time, only: &
      NOWDATE   => TIME_NOWDATE,   &
      NOWSUBSEC => TIME_NOWSUBSEC, &
      NOWDAYSEC => TIME_NOWDAYSEC
    use scale_file_fem_base, only: FILE_base_meshfield
    use scale_fem_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY

    implicit none
    class(MeshTopography), intent(inout) :: this

    character(len=H_LONG) :: TOPOGRAPHY_FEM_IN_BASENAME = ''                     !< basename of the input  file
    character(len=H_LONG) :: TOPOGRAPHY_FEM_IN_VARNAME  = 'topo'                 !< variable name of topo in the input  file
    character(len=H_LONG) :: TOPOGRAPHY_FEM_OUT_BASENAME = ''                          !< basename of the output file
    character(len=H_MID)  :: TOPOGRAPHY_FEM_OUT_TITLE    = 'SCALE-RM TOPOGRAPHY (FEM)' !< title    of the output file
    character(len=H_SHORT) :: TOPOGRAPHY_FEM_OUT_DTYPE    = 'DEFAULT'                   !< REAL4 or REAL8
    namelist / PARAM_TOPOGRAPHY_FEM / &
      TOPOGRAPHY_FEM_IN_BASENAME,  &
      TOPOGRAPHY_FEM_IN_VARNAME,   &
      TOPOGRAPHY_FEM_OUT_BASENAME, &
      TOPOGRAPHY_FEM_OUT_TITLE,    &
      TOPOGRAPHY_FEM_OUT_DTYPE   

    integer :: ierr

    class(MeshBase2D), pointer :: mesh2D

    type(FILE_base_meshfield) :: file
    logical :: file_existed
    integer, parameter :: vid_topo = 1
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("MeshTopography (FEM)",*) 'Write'
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPOGRAPHY_FEM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("ATMOS_DYN_DGM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("ATMOS_DYN_DGM_setup",*) 'Not appropriate names in namelist PARAM_TOPOHGRAPHY_FEM. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_TOPOGRAPHY_FEM)
    !--------------

    if ( TOPOGRAPHY_FEM_OUT_BASENAME /= '' .and. TOPOGRAPHY_FEM_OUT_BASENAME /= TOPOGRAPHY_FEM_IN_BASENAME ) then
      select type(mesh2D => this%topo%mesh)
      type is (MeshRectDom2D)
        call file%Init( 1, mesh2D=mesh2D )
      end select
      call file%Create( TOPOGRAPHY_FEM_OUT_BASENAME, TOPOGRAPHY_FEM_OUT_TITLE, TOPOGRAPHY_FEM_OUT_DTYPE, &
                        file_existed,                       &
                        myrank=PRC_myrank                   )

      if ( .not. file_existed ) then
        call file%Put_GlobalAttribute_time( NOWDATE, NOWSUBSEC )
      end if
      
      call file%Def_Var( this%topo, "TOPOGRAPHY", vid_topo, MFTYPE2D_XY, 'XY' )
      call file%End_def()
      call file%Write_var2D( vid_topo, this%topo, NOWDAYSEC, NOWDAYSEC )
      call file%Close()
      call file%Final()
    end if
    return
  end subroutine MeshTopography_write

  subroutine MeshTopography_read( this )
    use scale_prc, only: &
      PRC_myrank, PRC_abort
    use scale_topography, only: &
      TOPOGRAPHY_exist
    use scale_file_fem_base, only: FILE_base_meshfield
    use scale_fem_mesh_base2d, only: &
      MFTYPE2D_XY => MeshBase2D_DIMTYPEID_XY
    implicit none
    class(MeshTopography), intent(inout) :: this

    character(len=H_LONG) :: TOPOGRAPHY_FEM_IN_BASENAME = ''                     !< basename of the input  file
    character(len=H_LONG) :: TOPOGRAPHY_FEM_IN_VARNAME  = 'topo'                 !< variable name of topo in the input  file
    character(len=H_LONG) :: TOPOGRAPHY_FEM_OUT_BASENAME = ''                          !< basename of the output file
    character(len=H_MID)  :: TOPOGRAPHY_FEM_OUT_TITLE    = 'SCALE-RM TOPOGRAPHY (FEM)' !< title    of the output file
    character(len=H_SHORT) :: TOPOGRAPHY_FEM_OUT_DTYPE    = 'DEFAULT'                   !< REAL4 or REAL8
    namelist / PARAM_TOPOGRAPHY_FEM / &
      TOPOGRAPHY_FEM_IN_BASENAME,  &
      TOPOGRAPHY_FEM_IN_VARNAME,   &
      TOPOGRAPHY_FEM_OUT_BASENAME, &
      TOPOGRAPHY_FEM_OUT_TITLE,    &
      TOPOGRAPHY_FEM_OUT_DTYPE   

    integer :: ierr

    class(MeshBase2D), pointer :: mesh2D

    type(FILE_base_meshfield) :: file
    !---------------------------------------
    
    LOG_NEWLINE
    LOG_INFO("MeshTopography (FEM)",*) 'Write'
    
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=PARAM_TOPOGRAPHY_FEM,iostat=ierr)
    if( ierr < 0 ) then !--- missing
        LOG_INFO("ATMOS_DYN_DGM_setup",*) 'Not found namelist. Default used.'
    elseif( ierr > 0 ) then !--- fatal error
        LOG_ERROR("ATMOS_DYN_DGM_setup",*) 'Not appropriate names in namelist PARAM_TOPOHGRAPHY_FEM. Check!'
        call PRC_abort
    endif
    LOG_NML(PARAM_TOPOGRAPHY_FEM)
    !--------------

    if ( TOPOGRAPHY_FEM_IN_BASENAME /= '' ) then
      select type(mesh2D => this%topo%mesh)
      type is (MeshRectDom2D)
        call file%Init( 1, mesh2D=mesh2D )
      end select

      call file%Open( TOPOGRAPHY_FEM_IN_BASENAME, myrank=PRC_myrank)
      call file%Read_Var( MFTYPE2D_XY, TOPOGRAPHY_FEM_IN_VARNAME, this%topo )

      call file%Close()
      call file%Final()
    else
       if ( TOPOGRAPHY_exist ) then
          LOG_ERROR("MeshTopography_read",*) "Topography data for FEM is necessary"
          call PRC_abort
       end if
    end if

    return
  end subroutine MeshTopography_read

!> Interpolate FV topography data to DG nodes
!!
  subroutine MeshTopography_set_from_FVtopo( this, &
    TOPO_fv, FVtoDG_interp_order,                &
    IA, IS, IE, IHALO, JA, JS, JE, JHALO )
    use scale_fem_meshfield_util, only: &
      MeshFieldUtil_interp_FVtoDG_2D
    implicit none
    class(MeshTopography), intent(inout), target :: this
    integer, intent(in) :: IA, IS, IE, IHALO, JA, JS, JE, JHALO
    real(RP), intent(in) :: TOPO_fv(IA,JA)
    integer, intent(in) :: FVtoDG_interp_order

    class(LocalMesh2D), pointer :: lcmesh2D
    !---------------------------------------------

    lcmesh2D => this%topo%mesh%lcmesh_list(1)
    call MeshFieldUtil_interp_FVtoDG_2D( this%topo%local(1)%val,    & ! (out)
      TOPO_fv, lcmesh2D, lcmesh2D%refElem2D, IS, IE, IA, IHALO, JS, JE, JA, JHALO, & ! (in)
      FVtoDG_interp_order )                                                          ! (in)

    
    return
  end subroutine MeshTopography_set_from_FVtopo

end module scale_fem_mesh_topography
