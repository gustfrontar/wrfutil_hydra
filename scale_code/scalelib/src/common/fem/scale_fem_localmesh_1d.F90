!> module common / FEM local mesh / base
!!
!! @par Description
!!           A module for managing local mesh 1D with FEM 
!!
!! @author Yuta Kawai, Team SCALE
!!
!<
#include "scalelib.h"
module scale_fem_localmesh_1d

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
  use scale_fem_element_base, only: ElementBase, ElementBase1D

  !-----------------------------------------------------------------------------
  implicit none
  private

  !-----------------------------------------------------------------------------
  !
  !++ Public type & procedure
  ! 
  type, extends(LocalMeshBase), public :: LocalMesh1D

    type(ElementBase1D), pointer :: refElem1D
    real(RP) :: xmin, xmax 
  end type LocalMesh1D

  public :: LocalMesh1D_Init, LocalMesh1D_Final

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
  subroutine LocalMesh1D_Init( this, &
    lcdomID, refElem, myrank )
    
    implicit none

    type(LocalMesh1D), intent(inout) :: this
    integer, intent(in) :: lcdomID
    class(ElementBase1D), intent(in), target :: refElem
    integer, intent(in), optional :: myrank
    !-------------------------------------------------

#if defined(__GFORTRAN__) && __GNUC__ < 5
    LOG_ERROR('LocalMesh1D_Init',*) 'GCC version should be >= 5! STOP.'
    call PRC_abort
#else
    this%refElem1D => refElem
#endif
    call LocalMeshBase_Init(this, lcdomID, refElem, 1, myrank)

    return
  end subroutine LocalMesh1D_Init

  subroutine LocalMesh1D_Final( this, is_generated )
    implicit none

    type(LocalMesh1D), intent(inout) :: this
    logical, intent(in) :: is_generated
    !-------------------------------------------------

    call LocalMeshBase_Final( this, is_generated )

    return
  end subroutine LocalMesh1D_Final
  
end module scale_fem_localmesh_1d