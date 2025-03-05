!> module common / FEM local mesh / base
!!
!! @par Description
!!           A module for managing local mesh 3D with FEM 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scalelib.h"
module scale_fem_localmesh_3d

  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use scale_precision
  use scale_io
  use scale_prc, only: &
    PRC_abort
  use scale_fem_localmesh_base, only: &
    LocalMeshBase, LocalMeshBase_Init, LocalMeshBase_Final
  use scale_fem_localmesh_2d, only: LocalMesh2D
  use scale_fem_element_base, only: elementbase, elementbase3D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh3D

    type(ElementBase3D), pointer :: refElem3D
    real(RP) :: xmin, xmax
    real(RP) :: ymin, ymax
    real(RP) :: zmin, zmax
    integer :: NeX
    integer :: NeY
    integer :: NeZ
    integer :: Ne2D
    integer :: Ne2DA

    real(RP), allocatable :: Sz(:,:)
    real(RP), allocatable :: zS(:,:)
    real(RP), allocatable :: GI3(:,:,:)    !< The contravariant component of metric tensor with vertical general coordinate 
    real(RP), allocatable :: GsqrtH(:,:)   !< The Jacobian of horizontal transformation in the computational coordinate
    real(RP), allocatable :: zlev(:,:)
    real(RP), allocatable :: gam(:,:)      !< Factor for approximation with spherical shell domain (= r/a)

    class(LocalMesh2D), pointer :: lcmesh2D
    real(RP), allocatable :: lon2D(:,:)     
    real(RP), allocatable :: lat2D(:,:)

    integer, allocatable :: EMap3Dto2D(:)  
  contains
    procedure :: SetLocalMesh2D => LocalMesh3D_setLocalMesh2D  
  end type LocalMesh3D

  public :: LocalMesh3D_Init
  public :: LocalMesh3D_Final

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
  subroutine LocalMesh3D_Init( this, &
    lcdomID, refElem, myrank )
    implicit none

    class(LocalMesh3D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase3D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

#if defined(__GFORTRAN__) && __GNUC__ < 5
    LOG_ERROR('LocalMesh3D_Init',*) 'GCC version should be >= 5! STOP.'
    call PRC_abort
#else
    this%refElem3D  => refElem
#endif
    nullify( this%lcmesh2D )

    call LocalMeshBase_Init( this, lcdomID, refElem, 3, myrank )

    return
  end subroutine LocalMesh3D_Init

!OCL SERIAL
  subroutine LocalMesh3D_Final( this, is_generated )
    implicit none
    type(LocalMesh3D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    call LocalMeshBase_Final( this, is_generated )
    if (is_generated) then
      deallocate( this%zS, this%Sz )
      deallocate( this%GI3, this%GsqrtH )
      deallocate( this%zlev )
      deallocate( this%gam )
      deallocate( this%lon2D, this%lat2D )
      deallocate( this%EMap3Dto2D )
    end if

    return
  end subroutine LocalMesh3D_Final
  
!OCL SERIAL
  subroutine LocalMesh3D_setLocalMesh2D( this, lcmesh2D )
    implicit none
    class(LocalMesh3D), intent(inout) :: this
    class(LocalMesh2D), intent(in), target :: lcmesh2D
    !-------------------------------------------------

    this%lcmesh2D => lcmesh2D

    return
  end subroutine LocalMesh3D_setLocalMesh2D

end module scale_fem_localmesh_3d