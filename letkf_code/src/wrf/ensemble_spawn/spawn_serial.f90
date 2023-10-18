Program spawn_serial
!=====================================================================
!
! A program to spawn multiple programs using MPI_Comm_spawn and _multiple
!
!=====================================================================

implicit none
include 'mpif.h'

Integer,allocatable ::  errcode(:)
Character(len=100)  ::  arg
!Character(len=10) , dimension(3) :: args_to_wrapper !The size of this string is limited to 50 
                                                    !setting it to higher values results in 
                                                    !segmentation fault (with ifort)
character(len=20)   :: args_to_wrapper(1) = [character(len=20) :: ""]
character(len=100)  ::  command
Character(len=100)  ::  ibasepath
Character(len=100)  ::  tmpdir
Character(len=5)    ::  ensmember
Integer             ::  nprocs , nchilds , inimem , endmem
Integer             ::  ierr , i
Integer             ::  nranks , myrank
Integer      :: info
Integer      :: intercomm

!---------------------------------------------------------------------

CALL GETARG ( 1, arg )
command=TRIM(arg)    !The serial command to be executed.
CALL GETARG ( 2, arg )
READ(arg,*)nprocs         !The number of procs to be allocated (serial command will use only one of those)
CALL GETARG ( 3, arg )
ibasepath=arg             !The basepath where the serial comand and the mpi_serial_wrapper will be run.
CALL GETARG ( 4, arg )
READ(arg,*)inimem         !The initial ensemble member for this job range.
CALL GETARG ( 5, arg )
READ(arg,*)endmem         !The final ensemble member for this job range.

Call MPI_Init(ierr)
Call MPI_COMM_SIZE(MPI_COMM_WORLD,nranks,ierr)
Call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

allocate(errcode(nprocs))

nchilds=endmem-inimem+1

DO i=1,nchilds
    if ( myrank == 0 ) then
       WRITE(*,*) "JOB", i
       call MPI_Info_create(info,ierr)
       WRITE(ensmember,'(I2.2)')i+inimem-1
       tmpdir=TRIM(ibasepath)//'/'//TRIM(ensmember)//'/'

       call MPI_Info_set(info,"wdir",tmpdir,ierr)
       !args_to_wrapper(2)=TRIM(tmpdir) 
       args_to_wrapper(1) = TRIM(command) 
       !args_to_wrapper(2) = TRIM(tmpdir)
       !args_to_wrapper=[ TRIM(command) , TRIM(tmpdir) , '' ]
    endif   
    !WRITE(*,*)intercomm
   
    Call MPI_Comm_spawn('./mpi_serial_wrapper.exe',command,nprocs,info,0,MPI_COMM_WORLD, &
    &      intercomm, errcode, ierr)
 
    !Call MPI_Comm_spawn('./mpi_serial_wrapper.exe',MPI_ARGV_NULL,nprocs,info,0,MPI_COMM_WORLD, &
    !&      intercomm, errcode, ierr)

!    Call MPI_Comm_spawn('./mpi_serial_wrapper.exe',command,nprocs,info,0,MPI_COMM_WORLD, &
!    &      intercomm, errcode, ierr)
    WRITE(*,*) "IERR", ierr
    if (ierr /= 0) then
       WRITE(*,*) "MPI_Comm_spawn failure",i
       STOP 1
    endif

ENDDO

call MPI_COMM_FREE(intercomm,ierr)
Call MPI_Finalize(ierr)

STOP 0

End Program spawn_serial
