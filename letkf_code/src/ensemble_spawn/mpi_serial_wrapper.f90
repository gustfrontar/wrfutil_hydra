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
  Character(50)  :: command
  Character(50)  :: workdir1,workdir2
  Character(100) :: workdir
!---------------------------------------------------------------------
  
  !Get the serial command to be executed
  CALL GETARG ( 1, command )    !Command to be executed.
  CALL GETARG ( 2, workdir1 )   !Working directory to run the command (part1)
  CALL GETARG ( 3, workdir2 )   !Working directory to run the command (part2)
  !Join the work directory all together.
  workdir = TRIM(workdir1)//TRIM(workdir2) 

  Call MPI_Init(ierr)

  Call MPI_COMM_SIZE(MPI_COMM_WORLD,nranks,ierr)
  Call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

  !Command should be a serial task which is run from rank=0
  If ( myrank == 0 ) then
     WRITE(*,*)'Launching the program ',command,' at dir ',workdir
     Call Chdir( workdir )
     Call system( 'cd '//workdir )
     Call system( command )
  Endif
  Call MPI_Finalize(ierr)
  STOP 0 

End Program mpi_serial_wrapper

