!> module common / FEM mesh / Map projection
!!
!! @par Description
!!      Mangage map projection
!! @author Yuta Kawai, Team SCALE
#include "scalelib.h"
module scale_fem_mesh_mapprojection

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_sparsemat, only: SparseMat
  use scale_fem_localmesh_3d, only: LocalMesh3D
  use scale_fem_localmesh_2d, only: LocalMesh2D
  use scale_fem_mesh_base2d, only: MeshBase2D
  use scale_fem_mesh_rectdom2d, only: MeshRectDom2D
  use scale_fem_mesh_base3d, only: MeshBase3D
  use scale_fem_element_base, only: &
    ElementBase2D, ElementBase3D
  use scale_fem_meshfield_base, only: &
    MeshField2D, MeshField3D
  use scale_comm_fem_meshfield_base, only: &
    MeshFieldCommBase, MeshFieldContainer

  use scale_mapprojection, only: &
    MAPPROJECTION_basepoint_lon, &
    MAPPROJECTION_basepoint_lat, &
    MAPPROJECTION_xy2lonlat,     &
    MAPPROJECTION_mapfactor    
    
  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type MapProjection_local_data
    real(RP), allocatable :: Gam(:,:,:,:,:) !< Christoffel symbols of the second kind
  end type
  type, public :: MeshMapProjection
    type(MapProjection_local_data), allocatable :: local(:)
  contains
    procedure :: Init => MeshMapProjection_Init
    procedure :: Final => MeshMapProjection_Final
    procedure :: SetHCoordinate => MeshMapProjection_set_hcoordinate
  end type MeshMapProjection

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
  subroutine MeshMapProjection_Init( this, mesh )
    implicit none

    class(MeshMapProjection), intent(inout) :: this
    class(MeshBase2D), intent(in), target :: mesh

    class(LocalMesh2D), pointer :: lcmesh2D
    integer :: n
    !-----------------------------------------------------------------------------

    allocate( this%local(mesh%LOCAL_MESH_NUM) )
    do n=1, mesh%LOCAL_MESH_NUM
      lcmesh2D => mesh%lcmesh_list(n)
      allocate( this%local(n)%Gam(lcmesh2D%refElem2D%Np,2,2,lcmesh2D%Ne,2) )
    end do
    return
  end subroutine MeshMapProjection_Init

!OCL SERIAL
  subroutine MeshMapProjection_Final( this )
    implicit none
    class(MeshMapProjection), intent(inout) :: this

    integer :: n    
    !-----------------------------------------------------------------------------

    do n=1, size(this%local)
      deallocate( this%local(n)%Gam )
    end do
    deallocate( this%local )

    return
  end subroutine MeshMapProjection_Final
  
!OCL SERIAL
  subroutine MeshMapProjection_set_hcoordinate( this, mesh3D, &
    comm2D                                    )    
    use scale_fem_meshutil_vcoord, only: MeshUtil_VCoord_GetMetric
    implicit none

    class(MeshMapProjection), target, intent(inout) :: this
    class(MeshBase3D), intent(inout), target :: mesh3D
    class(MeshFieldCommBase), intent(inout) :: comm2D

    type(SparseMat) :: Dx2D, Dy2D, Lift2D
    class(LocalMesh3D), pointer :: lcmesh
    class(LocalMesh2D), pointer :: lcmesh2D
    class(ElementBase2D), pointer :: elem2D
    class(MeshBase2D), pointer :: mesh2D

    integer :: n
    integer :: ke, ke2D, p

    type(MeshFieldContainer) :: comm2d_varlist(2)

    type(MeshField2D), target :: h1, h2
    real(RP), allocatable :: m1(:,:), m2(:,:)
    !-------------------------------------------------------------

    lcmesh => mesh3D%lcmesh_list(1)
    lcmesh2D => lcmesh%lcmesh2D
    call Dx2D  %Init( lcmesh2D%refElem2D%Dx1  )
    call Dy2D  %Init( lcmesh2D%refElem2D%Dx2  )
    call Lift2D%Init( lcmesh2D%refElem2D%Lift )

    call mesh3D%GetMesh2D( mesh2D )

    ! Calculate metric factors associated with map projection

    call h1%Init( "h1", "", mesh2D )
    call h2%Init( "h2", "", mesh2D )

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh => mesh3D%lcmesh_list(n)
      lcmesh2D => lcmesh%lcmesh2D
      elem2D => lcmesh2D%refElem2D

      !$omp parallel do private(ke2D, p) collapse(2)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
      do p=1, elem2D%Np
        call MAPPROJECTION_xy2lonlat( lcmesh2D%pos_en(p,ke2D,1), lcmesh2D%pos_en(p,ke2D,2), & ! (in)
          lcmesh2D%lon(p,ke2D), lcmesh2D%lat(p,ke2D) ) ! (out)
      end do
      end do

      allocate( m1(1:elem2D%Np,1:lcmesh2D%Ne), m2(1:elem2D%Np,1:lcmesh2D%Ne) )
      call MAPPROJECTION_mapfactor( elem2D%Np, 1, elem2D%Np, lcmesh2D%Ne, 1, lcmesh2D%Ne, lcmesh2D%lat(:,:), & ! (in)
        m1(:,:), m2(:,:) ) ! (out)
      
      !$omp parallel do private(ke2D)
      do ke2D=lcmesh2D%NeS, lcmesh2D%NeE
        h1%local(n)%val(:,ke2D) = 1.0_RP / m1(:,ke2D)
        h2%local(n)%val(:,ke2D) = 1.0_RP / m2(:,ke2D)
      end do
      deallocate( m1, m2)
    end do

    ! Exchange scale factor to fill halo

    comm2d_varlist(1)%field2d => h1
    comm2d_varlist(2)%field2d => h2
    call comm2D%Put( comm2d_varlist, 1 )
    call comm2D%Exchange()
    call comm2D%Get( comm2d_varlist, 1 )

    ! Set metric data managed by local mesh

    do n=1, mesh3D%LOCAL_MESH_NUM
      lcmesh2D => lcmesh%lcmesh2D
      call cal_Christoffel_symbols( this%local(n)%Gam, &
        h1%local(n)%val, h2%local(n)%val, lcmesh2D, lcmesh2D%refElem2D, &
        Dx2D, Dy2D, Lift2D )

      call mesh3D%Set_geometric_with_hcoord( n, h1%local(n)%val, h2%local(n)%val )
    end do

    !---
    call h1%Final(); call h2%Final()

    call Dx2D%Final()
    call Dy2D%Final()
    call Lift2D%Final()

    return
  end subroutine MeshMapProjection_set_hcoordinate

!-- private

!OCL SERIAL
  subroutine cal_Christoffel_symbols( Gam, &
    h1, h2, lcmesh2D, elem2D,  &
    Dx2D, Dy2D, Lift2D )
    use scale_sparsemat, only: sparsemat_matmul
    implicit none

    class(LocalMesh2D), intent(in) :: lcmesh2D
    class(ElementBase2D), intent(in) :: elem2D
    real(RP), intent(out) :: Gam(elem2D%Np,2,2,lcmesh2D%Ne,2)
    real(RP), intent(in) :: h1(elem2D%Np,lcmesh2D%NeA)
    real(RP), intent(in) :: h2(elem2D%Np,lcmesh2D%NeA)
    type(SparseMat), intent(in) :: Dx2D
    type(SparseMat), intent(in) :: Dy2D
    type(SparseMat), intent(in) :: Lift2D

    integer :: ke, ke2D
    real(RP) :: del_flux_h1(elem2D%NfpTot,lcmesh2D%Ne,2)
    real(RP) :: del_flux_h2(elem2D%NfpTot,lcmesh2D%Ne,2)
    real(RP) :: Fx2D(elem2D%Np), Fy2D(elem2D%Np), LiftDelFlux2D(elem2D%Np)
    real(RP) :: Grad_h1(elem2D%Np,lcmesh2D%Ne,2)
    real(RP) :: Grad_h2(elem2D%Np,lcmesh2D%Ne,2)
    !------------------------------------------------


    call cal_del_flux( del_flux_h1, del_flux_h2, &
        h1, h2, lcmesh2D%normal_fn(:,:,1), lcmesh2D%normal_fn(:,:,2), &
        lcmesh2D%VMapM, lcmesh2D%VMapP, lcmesh2D, elem2D              )
    
    !$omp parallel private(ke2D, ke,        &
    !$omp Fx2D, Fy2D, LiftDelFlux2D )

    !$omp do
    do ke2D=1, lcmesh2D%Ne
      call sparsemat_matmul( Dx2D, h1(:,ke2D), Fx2D )
      call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux_h1(:,ke2D,1), LiftDelFlux2D)
      Grad_h1(:,ke2D,1) = lcmesh2D%Escale(:,ke2D,1,1) * Fx2D(:) + LiftDelFlux2D(:)

      call sparsemat_matmul( Dy2D, h1(:,ke2D), Fy2D )
      call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux_h1(:,ke2D,2), LiftDelFlux2D)
      Grad_h1(:,ke2D,2) = lcmesh2D%Escale(:,ke2D,2,2) * Fy2D(:) + LiftDelFlux2D(:)

      call sparsemat_matmul( Dx2D, h2(:,ke2D), Fx2D )
      call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux_h2(:,ke2D,1), LiftDelFlux2D)
      Grad_h2(:,ke2D,1) = lcmesh2D%Escale(:,ke2D,1,1) * Fx2D(:) + LiftDelFlux2D(:)

      call sparsemat_matmul( Dy2D, h2(:,ke2D), Fy2D )
      call sparsemat_matmul( Lift2D, lcmesh2D%Fscale(:,ke2D) * del_flux_h2(:,ke2D,2), LiftDelFlux2D)
      Grad_h2(:,ke2D,2) = lcmesh2D%Escale(:,ke2D,2,2) * Fy2D(:) + LiftDelFlux2D(:)
    end do
    !$omp end do

    !$omp do
    do ke2D=1, lcmesh2D%Ne
      Gam(:,1,1,ke2D,1) = 1.0_RP / h1(:,ke2D) * Grad_h1(:,ke2D,1)
      Gam(:,2,1,ke2D,1) = 1.0_RP / h1(:,ke2D) * Grad_h1(:,ke2D,2)
      Gam(:,1,2,ke2D,1) = Gam(:,2,1,ke2D,1) 
      Gam(:,2,2,ke2D,1) = - 1.0_RP / h1(:,ke2D)**2 * ( h2(:,ke2D) * Grad_h2(:,ke2D,1) )

      Gam(:,1,1,ke2D,2) = - 1.0_RP / h2(:,ke2D)**2 * ( h1(:,ke2D) * Grad_h1(:,ke2D,2) )
      Gam(:,2,1,ke2D,2) = 1.0_RP / h2(:,ke2D) * Grad_h2(:,ke2D,1)
      Gam(:,1,2,ke2D,2) = Gam(:,2,1,ke2D,2) 
      Gam(:,2,2,ke2D,2) = 1.0_RP / h2(:,ke2D) * Grad_h2(:,ke2D,2)
    end do
    !$omp end do
    !$omp end parallel

    return
  end subroutine cal_Christoffel_symbols

!OCL SERIAL
  subroutine cal_del_flux( del_flux_h1, del_flux_h2, &
    h1, h2, nx, ny, vmapM, vmapP, lmesh, elem )
    implicit none
    type(LocalMesh2D), intent(in) :: lmesh
    type(ElementBase2D), intent(in) :: elem
    real(RP), intent(out) :: del_flux_h1(elem%NfpTot*lmesh%Ne,2)
    real(RP), intent(out) :: del_flux_h2(elem%NfpTot*lmesh%Ne,2)
    real(RP), intent(in) :: h1(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: h2(elem%Np*lmesh%NeA)
    real(RP), intent(in) :: nx(elem%NfpTot*lmesh%Ne)
    real(RP), intent(in) :: ny(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapM(elem%NfpTot*lmesh%Ne)
    integer, intent(in) :: vmapP(elem%NfpTot*lmesh%Ne)

    integer :: i
    integer :: iP, iM
    real(RP) :: del
    !-------------------------------------

    !$omp parallel do private( i, iM, iP, del )
    do i=1, elem%NfpTot * lmesh%Ne
      iM = vmapM(i); iP = vmapP(i)

      del = 0.5_RP * ( h1(iP) - h1(iM) )
      del_flux_h1(i,1) = del * nx(i)
      del_flux_h1(i,2) = del * ny(i)

      del = 0.5_RP * ( h2(iP) - h2(iM) )
      del_flux_h2(i,1) = del * nx(i)
      del_flux_h2(i,2) = del * ny(i)
    end do

    return
  end subroutine cal_del_flux

end module scale_fem_mesh_mapprojection