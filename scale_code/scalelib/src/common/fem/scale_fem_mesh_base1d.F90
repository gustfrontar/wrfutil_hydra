!> module common / FEM mesh / base 1D
!!
!! @par Description
!!           A module for managing mesh with FEM 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scalelib.h"
module scale_fem_mesh_base1d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision

  use scale_fem_localmesh_1d, only: &
    LocalMesh1D, LocalMesh1D_Init, LocalMesh1D_Final
  
  use scale_fem_mesh_base, only: MeshBase, &
    MeshBase_Init, MeshBase_Final, &
    MeshBase_setGeometricInfo
  
  use scale_fem_element_base, only: ElementBase1d
  use scale_fem_element_line, only: LineElement

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, abstract, public, extends(MeshBase) :: Meshbase1d
    type(LocalMesh1D), allocatable :: lcmesh_list(:)
    class(ElementBase1D), pointer :: refElem1D
    
    integer, public :: NeG
    integer, public :: Nprc
    real(RP), public :: xmin_gl, xmax_gl
    real(RP), public, allocatable :: FX(:)
  contains
    procedure(Meshbase1d_generate), deferred :: Generate 
    procedure :: GetLocalMesh => MeshBase1D_get_localmesh
  end type Meshbase1d

  interface 
    subroutine Meshbase1d_generate(this)
      import Meshbase1d
      class(Meshbase1d), intent(inout), target :: this
    end subroutine Meshbase1d_generate
  end interface

  public :: Meshbase1d_Init, Meshbase1d_Final
  public :: Meshbase1d_setGeometricInfo, Meshbase1D_assignDomID, MeshBase1D_setupLocalDom

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer, public :: MeshBase1D_DIMTYPE_NUM   = 2
  integer, public :: MeshBase1D_DIMTYPEID_X   = 1
  integer, public :: MeshBase1D_DIMTYPEID_XT  = 2
  
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !

contains
  subroutine Meshbase1d_Init( this,         &
    NeG,                                    &
    dom_xmin, dom_xmax,                     &
    refElem, NLocalMeshPerPrc,              &
    nprocs, myrank, FX                      )
    
    implicit none

    class(Meshbase1d), intent(inout) :: this
    integer, intent(in) :: NeG
    real(RP), intent(in) :: dom_xmin
    real(RP), intent(in) :: dom_xmax   
    class(ElementBase1d), intent(in), target :: refElem
    integer, intent(in) :: NLocalMeshPerPrc
    integer, intent(in), optional :: nprocs
    integer, intent(in), optional :: myrank
    real(RP), intent(in), optional :: FX(NeG+1)

    integer :: n
    integer :: k
    real(RP) :: dx
    !-----------------------------------------------------------------------------
    
    this%NeG = NeG

    this%xmin_gl       = dom_xmin
    this%xmax_gl       = dom_xmax

    !- Fx
    allocate( this%FX(NeG+1) )
    if ( present(FX) ) then
      this%FX(:) = FX(:)
    else
      this%FX(1    ) = dom_xmin
      this%FX(NeG+1) = dom_xmax
      dx = (dom_xmax - dom_xmin) / dble(NeG)
      do k=2, NeG
        this%FX(k) = this%FX(k-1) + dx
      end do
    end if

    this%refElem1D => refElem
    call MeshBase_Init( this,          &
      MeshBase1D_DIMTYPE_NUM, refElem, &
      NLocalMeshPerPrc, 2,             &
      nprocs                           )

    this%Nprc = this%PRC_NUM
         
    allocate( this%lcmesh_list(this%LOCAL_MESH_NUM) )
    do n=1, this%LOCAL_MESH_NUM
      call LocalMesh1D_Init( this%lcmesh_list(n), n, refElem, myrank )
    end do

    call this%SetDimInfo( MeshBase1D_DIMTYPEID_X, "FEM_x", "m", "X-coordinate" )
    call this%SetDimInfo( MeshBase1D_DIMTYPEID_XT, "FEM_xt", "m", "X-coordinate" )

    return
  end subroutine Meshbase1d_Init

  subroutine Meshbase1d_Final( this )
    implicit none    
    class(Meshbase1d), intent(inout) :: this

    integer :: n
    !-----------------------------------------------------------------------------
  
    if ( allocated ( this%lcmesh_list ) ) then 
      do n=1, this%LOCAL_MESH_NUM
        call LocalMesh1D_Final( this%lcmesh_list(n), this%isGenerated )
      end do
  
      deallocate( this%lcmesh_list )
    end if
    
    call MeshBase_Final(this)

    return
  end subroutine Meshbase1d_Final
  
  subroutine MeshBase1D_get_localmesh( this, id, ptr_lcmesh )
    use scale_fem_localmesh_base, only: LocalMeshBase
    implicit none

    class(MeshBase1D), target, intent(in) :: this
    integer, intent(in) :: id
    class(LocalMeshBase), pointer, intent(out) :: ptr_lcmesh
    !-------------------------------------------------------------

    ptr_lcmesh => this%lcmesh_list(id)
    return
  end subroutine MeshBase1D_get_localmesh

  subroutine Meshbase1d_setGeometricInfo( lcmesh )
    implicit none
    
    type(LocalMesh1D), intent(inout) :: lcmesh

    class(Elementbase1D), pointer :: refElem
    integer :: kelem
    integer :: f
    integer :: i

    integer :: node_ids(lcmesh%refElem%Nv)
    real(RP) :: vx(lcmesh%refElem%Nv)
    real(RP) :: xr(lcmesh%refElem%Np)
    real(DP) :: Escale(1,1,lcmesh%refElem%Np)
    integer :: fmask(lcmesh%refElem%NfpTot)
    integer :: fid(lcmesh%refElem1D%Nfp,lcmesh%refElem1D%Nfaces)

  !-----------------------------------------------------------------------------

    call MeshBase_setGeometricInfo(lcmesh, 1)
    refElem => lcmesh%refElem1D
    
    fmask(:) = reshape(refElem%Fmask, shape(fmask))
    !$omp parallel private(f,i,kelem, node_ids,vx,xr)
    !$omp do
    do f=1, refElem%Nfaces
    do i=1, refElem%Nfp
       fid(i,f) = i + (f-1)*refElem%Nfp
    end do
    end do
    !$omp end do
    !$omp do
    do kelem=1, lcmesh%Ne
      node_ids(:) = lcmesh%EToV(kelem,:)
      vx(:) = lcmesh%pos_ev(node_ids(:),1)
      lcmesh%pos_en(:,kelem,1) = vx(1) + 0.5_RP*(refElem%x1(:) + 1.0_RP)*(vx(2) - vx(1))
      
      xr(:) = 0.5_RP*(vx(2) - vx(1)) !matmul(refElem%Dr,mesh%x(:,kelem))

      lcmesh%J(:,kelem) = xr
      lcmesh%Escale(:,kelem,1,1) =   1.0_RP/lcmesh%J(:,kelem)

      !* Face

      !
      !mesh%fx(:,n) = mesh%x(fmask(:),n)
      !mesh%fy(:,n) = mesh%y(fmask(:),n)

      ! Calculate normal vectors
      lcmesh%normal_fn(fid(:,1),kelem,1) = - 1.0_RP
      lcmesh%normal_fn(fid(:,2),kelem,1) = + 1.0_RP
      lcmesh%sJ(:,kelem) = 1.0_RP
      lcmesh%Fscale(:,kelem) = lcmesh%sJ(:,kelem)/lcmesh%J(fmask(:),kelem)
      lcmesh%Gsqrt(:,kelem) = 1.0_RP
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine Meshbase1d_setGeometricInfo

  subroutine Meshbase1D_assignDomID( this, &
    tileID_table, panelID_table,        &
    pi_table )
  
    use scale_fem_meshutil_1d, only: &       
      MeshUtil1D_buildGlobalMap
    implicit none
    
    class(Meshbase1d), intent(inout) :: this    
    integer, intent(out) :: tileID_table(this%LOCAL_MESH_NUM, this%PRC_NUM)
    integer, intent(out) :: panelID_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    integer, intent(out) :: pi_table(this%LOCAL_MESH_NUM*this%PRC_NUM)
    
    integer :: ldomID
    integer :: p
    integer :: tileID
    !-----------------------------------------------------------------------------
    
    call MeshUtil1D_buildGlobalMap( &
      panelID_table, pi_table,                                                      & ! (out)
      this%tileID_globalMap, this%tileFaceID_globalMap, this%tilePanelID_globalMap, & ! (out)
      this%LOCAL_MESH_NUM_global )                                                    ! (in)

    !$omp parallel do private(ldomID, tileID)
    do p=1, this%PRC_NUM
    do ldomID=1, this%LOCAL_MESH_NUM
      tileID = ldomID + (p-1)*this%LOCAL_MESH_NUM

      tileID_table(ldomID,p)              = tileID
      this%tileID_global2localMap(tileID) = ldomID
      this%PRCRank_globalMap(tileID)      = p - 1
    end do
    end do

    return
  end subroutine Meshbase1D_assignDomID

  subroutine MeshBase1D_setupLocalDom( lcmesh,  &
    tileID, panelID,                            &
    i, Nprc,                                    &
    dom_xmin, dom_xmax,                         &
    Ne, FX )

    use scale_fem_meshutil_1d, only: &
      MeshUtil1D_genConnectivity,   &
      MeshUtil1D_genLineDomain,     &
      MeshUtil1D_BuildInteriorMap,  &
      MeshUtil1D_genPatchBoundaryMap

    use scale_fem_localmesh_base, only: BCTYPE_INTERIOR
    implicit none
      
    type(LocalMesh1D), intent(inout) :: lcmesh
    integer, intent(in) :: tileID
    integer, intent(in) :: panelID
    integer, intent(in) :: i
    integer, intent(in) :: Nprc
    real(RP) :: dom_xmin, dom_xmax
    integer, intent(in) ::Ne
    real(RP), intent(in) :: FX(Ne*Nprc+1)
    
    class(ElementBase1D), pointer :: elem
    real(RP) :: delx
    real(RP) :: FX_lc(Ne+1)    
    !-----------------------------------------------------------------------------

    elem => lcmesh%refElem1D
    lcmesh%tileID = tileID
    lcmesh%panelID = panelID
    
    !--

    lcmesh%Ne  = Ne
    lcmesh%Nv  = Ne + 1
    lcmesh%NeS = 1
    lcmesh%NeE = lcmesh%Ne
    lcmesh%NeA = lcmesh%Ne + 2

    !delx = (dom_xmax - dom_xmin)/dble(Nprc)
    FX_lc(:) = Fx((i-1)*Ne+1:i*Ne)
    lcmesh%xmin = FX_lc(1)
    lcmesh%xmax = FX_lc(Ne+1)

    allocate(lcmesh%pos_ev(lcmesh%Nv,1))
    allocate( lcmesh%EToV(lcmesh%Ne,2) )
    allocate( lcmesh%EToE(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%EToF(lcmesh%Ne,elem%Nfaces) )
    allocate( lcmesh%BCType(lcmesh%refElem%Nfaces,lcmesh%Ne) )
    allocate( lcmesh%VMapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%VMapP(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapM(elem%NfpTot, lcmesh%Ne) )
    allocate( lcmesh%MapP(elem%NfpTot, lcmesh%Ne) )

    lcmesh%BCType(:,:) = BCTYPE_INTERIOR

    !----

    call MeshUtil1D_genLineDomain( lcmesh%pos_ev, lcmesh%EToV,   & ! (out)
        lcmesh%Ne, lcmesh%xmin, lcmesh%xmax, FX=FX_lc   )            ! (in)


    !---
    call MeshBase1D_setGeometricInfo( lcmesh )
    
    !---

    call MeshUtil1D_genConnectivity( lcmesh%EToE, lcmesh%EToF, & ! (out)
        & lcmesh%EToV, lcmesh%Ne, elem%Nfaces )                  ! (in)

    !---
    call MeshUtil1D_BuildInteriorMap( lcmesh%VmapM, lcmesh%VMapP, lcmesh%MapM, lcmesh%MapP,  &
      & lcmesh%pos_en, lcmesh%pos_ev, lcmesh%EToE, lcmesh%EtoF, lcmesh%EtoV,                   &
      & elem%Fmask, lcmesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, lcmesh%Nv )

    call MeshUtil1D_genPatchBoundaryMap( lcmesh%VMapB, lcmesh%MapB, lcmesh%VMapP, &
      & lcmesh%pos_en, lcmesh%xmin, lcmesh%xmax,                                  &
      & elem%Fmask, lcmesh%Ne, elem%Np, elem%Nfp, elem%Nfaces, lcmesh%Nv)
    
    return
  end subroutine MeshBase1D_setupLocalDom

end module scale_fem_mesh_base1d