Program spawn
!=====================================================================
!
!       test_comm_spawn
!
! A program to test MPI_Comm_spawn and _multiple
!
!=====================================================================
  implicit none
  include 'mpif.h'

  Integer,allocatable ::  errcode(:)
  Character(len=100)  ::  arg
  Character(len=100)  ::  icommand
  Character(len=100)  ::  ibasepath
  Character(len=100)  ::  tmpdir
  Character(len=5)    ::  ensmember
  Integer             ::  nprocs , nchilds , inimem , endmem 
  Integer             ::  ierr , i 
  Integer      :: info  
  Integer      :: intercomm
  
!---------------------------------------------------------------------

  CALL GETARG ( 1, arg )
  WRITE(*,*)arg
  icommand=arg
  CALL GETARG ( 2, arg )
  WRITE(*,*)arg
  READ(arg,*)nprocs
  CALL GETARG ( 3, arg )
  WRITE(*,*)arg
  ibasepath=arg
  CALL GETARG ( 4, arg )
  WRITE(*,*)arg
  READ(arg,*)inimem
  CALL GETARG ( 5, arg )
  WRITE(*,*)arg
  READ(arg,*)endmem

Call MPI_Init(ierr)

  allocate(errcode(nprocs))
  errcode=0

  nchilds=endmem-inimem+1

  DO i=1,nchilds
    call MPI_Info_create(info,ierr)
    WRITE(ensmember,'(I5.5)')i+inimem-1
    tmpdir=TRIM(ibasepath)//'/WRF'//ensmember
    call MPI_Info_set(info,"wdir",tmpdir,ierr)
 
  WRITE(*,*) errcode
  Call MPI_Comm_spawn(icommand,MPI_ARGV_NULL,nprocs,info,0,MPI_COMM_WORLD, &
 &      intercomm, errcode, ierr)
    if (ierr /= 0) then
       WRITE(*,*) "MPI_Comm_spawn failure",i
       STOP
    endif
 
  ENDDO  
  call MPI_COMM_FREE(intercomm,ierr) 
  Call MPI_Finalize(ierr)
  
End Program spawn



