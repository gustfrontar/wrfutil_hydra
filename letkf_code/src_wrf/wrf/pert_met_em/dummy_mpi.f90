Program  mpi_serial_wrapper
!=====================================================================
!
!       mpi_serial_wrapper
!
!  A tool that can be use in combination with the spawn.exe to distribute
!  serial programs into different nodes. 
!
!
!=====================================================================
  Implicit None
  include 'mpif.h'
  Integer :: ierr , myrank , nranks
  Character(500) :: command
!---------------------------------------------------------------------
  
  !Get the serial command to be executed
  CALL GETARG ( 1, command )    !Command to be executed.
  Call MPI_Init(ierr)

  Call MPI_COMM_SIZE(MPI_COMM_WORLD,nranks,ierr)
  Call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

  !Command should be a serial task which is run from rank=0

  If ( myrank == 0 ) then
     Call system( command )
  Endif
  Call MPI_Finalize(ierr)

End Program mpi_serial_wrapper

