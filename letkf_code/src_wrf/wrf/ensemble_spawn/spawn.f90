Program spawn
!=====================================================================
!
! A program to spawn multiple programs using MPI_Comm_spawn and _multiple
!
!=====================================================================

implicit none
include 'mpif.h'

Integer,allocatable ::  errcode(:)
Character(len=100)  ::  arg
Character(len=100)  ::  command
Character(len=100)  ::  ibasepath
Character(len=100)  ::  tmpdir
Character(len=10)   ::  runtime_flags(1)
Character(len=5)    ::  ensmember
Integer             ::  nprocs , nchilds , inimem , endmem
Integer             ::  ierr , i
Integer      :: info
Integer      :: intercomm

!---------------------------------------------------------------------

CALL GETARG ( 1, arg )
command=arg               !MPI program to be executed.
CALL GETARG ( 2, arg )
READ(arg,*)nprocs         !Number of procs assigned to each program.
CALL GETARG ( 3, arg )
ibasepath=arg             !Basepath where the program is.
CALL GETARG ( 4, arg )
READ(arg,*)inimem         !Initial ensemble member for this job range
CALL GETARG ( 5, arg )
READ(arg,*)endmem         !Final ensemble member for this job range
CALL GETARG ( 6, arg )    
runtime_flags(1) = arg    !Runtime flags to be passed to the spawned executable.


Call MPI_Init(ierr)

allocate(errcode(nprocs))
errcode=0

nchilds=endmem-inimem+1

DO i=1,nchilds

    WRITE(*,*) "JOB", i
    call MPI_Info_create(info,ierr)
    WRITE(ensmember,'(I2.2)')i+inimem-1
    tmpdir=TRIM(ibasepath)//'/'//ensmember
    WRITE(*,*)tmpdir
    call MPI_Info_set(info,"wdir",tmpdir,ierr)
    if ( TRIM( runtime_flags(1)) == 'NULL' )then
       Call MPI_Comm_spawn(command,MPI_ARGV_NULL,nprocs,info,0,MPI_COMM_WORLD, &
    &      intercomm, errcode, ierr)
    else 
       Call MPI_Comm_spawn(command,runtime_flags,nprocs,info,0,MPI_COMM_WORLD, &
    &      intercomm, errcode, ierr)
    endif

    WRITE(*,*) "IERR", ierr
    if (ierr /= 0) then
       WRITE(*,*) "MPI_Comm_spawn failure",i
       STOP 1
    endif

ENDDO

call MPI_COMM_FREE(intercomm,ierr)
Call MPI_Finalize(ierr)

End Program spawn
