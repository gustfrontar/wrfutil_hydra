!> module common / FEM local mesh / base
!!
!! @par Description
!!           A module for managing local mesh 2D with FEM 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scalelib.h"
module scale_fem_localmesh_2d

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
  use scale_fem_element_base, only: ElementBase, ElementBase2D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh2D

    type(ElementBase2D), pointer :: refElem2D
    real(RP) :: xmin, xmax
    real(RP) :: ymin, ymax
    integer :: NeX, NeY

    real(RP), allocatable :: lon(:,:)     
    real(RP), allocatable :: lat(:,:)     
  end type LocalMesh2D

  public :: LocalMesh2D_Init, LocalMesh2D_Final

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
  subroutine LocalMesh2D_Init( this, &
    lcdomID, refElem, myrank )
    
    implicit none

    class(LocalMesh2D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase2D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

#if defined(__GFORTRAN__) && __GNUC__ < 5
    LOG_ERROR('LocalMesh2D_Init',*) 'GCC version should be >= 5! STOP.'
    call PRC_abort
#else
    this%refElem2D  => refElem
#endif
    call LocalMeshBase_Init(this, lcdomID, refElem, 2, myrank)

    return
  end subroutine LocalMesh2D_Init

  subroutine LocalMesh2D_Final( this, is_generated )
    implicit none
    type(LocalMesh2D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    if (is_generated) then
      deallocate( this%lon, this%lat )
    end if
    call LocalMeshBase_Final( this, is_generated )

    return
  end subroutine LocalMesh2D_Final
  
end module scale_fem_localmesh_2d