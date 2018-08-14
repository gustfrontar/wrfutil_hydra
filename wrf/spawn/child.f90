Program child
!=====================================================================
!
!       test_comm_spawn
!
! A program to test MPI_Comm_spawn and _multiple
!
!=====================================================================
  Implicit None
  Integer :: ierr
  
!---------------------------------------------------------------------


  Call MPI_Init(ierr)

    WRITE(*,*)"Hello world" 

 
  Call MPI_Finalize(ierr)
  
End Program child



