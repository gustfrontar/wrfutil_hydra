module covariance_matrix_tools
!=======================================================================
!
! [PURPOSE:] Common procedure for verification against grided data
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_verification
  USE map_utils
  USE ifport
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  !NAMELIST INPUT
  INTEGER :: nbv=1  !Number of ensemble members.
  integer, parameter :: max_npoints = 20
  integer :: npoints
  real(r_size) :: plon(max_npoints)  , plat(max_npoints) , plev(max_npoints) 
  character(50) :: pvarname(max_npoints)
  real(r_size) :: dep(max_npoints) , error(max_npoints)
  LOGICAL :: bootstrap=.true.
  INTEGER :: bootstrap_samples=20

  INTEGER :: skipx=1,skipy=1,skipz=1
  INTEGER :: nignore=0
  INTEGER , PARAMETER :: max_nignore = 20
  character(50) :: ignorevarname(max_nignore)
  
  character(50) :: input_endian='little_endian'
  character(50) :: output_endian='big_endian'

  CHARACTER(LEN=100) :: NAMELIST_FILE='./covariance_matrix.namelist'

  !GRID INFO
  TYPE(proj_info)  :: projection
  TYPE(ctl_info)   :: ctl

!  integer :: nv !Total number of variables (vars times levs)


CONTAINS

!-----------------------------------------------------------------------
! Read namelist
!-----------------------------------------------------------------------

SUBROUTINE read_namelist
IMPLICIT NONE
LOGICAL  :: file_exist
INTEGER  :: ierr

plon=undef
plat=undef
plev=undef
dep=undef
error=undef

!In the current version PLEV represents not only the level but also the
!variable. 

NAMELIST / GENERAL / nbv , npoints , plon , plat , pvarname , plev , &
                     dep , error , bootstrap , bootstrap_samples   , &
                     nignore , ignorevarname , skipx , skipy , skipz,&
                     input_endian , output_endian

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)

END SUBROUTINE read_namelist

!-----------------------------------------------------------------------
! File I/O
!-----------------------------------------------------------------------

!-- Read a grid file ---------------------------------------------------
SUBROUTINE read_grd(filename,nx,ny,nz,var,undefmask,undefbin)
  IMPLICIT NONE
  INTEGER,INTENT(IN)      :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(OUT) :: var(nx,ny,nz) 
  LOGICAL     ,INTENT(OUT) :: undefmask(nx,ny,nz)
  REAL(r_sngl) :: buf4(nx,ny)
  REAL(r_sngl),INTENT(IN) :: undefbin
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength,ios ,ii,jj
  LOGICAL :: file_exist
  

  undefmask=.true.

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength, &
              CONVERT=input_endian )
  ELSE
   WRITE(*,*)"[Warning]: Could not finde file! ",filename
   undefmask=.false.
   var=0.0d0
   RETURN
  ENDIF

  DO n=1,nz
    READ(iunit,rec=n,IOSTAT=ios)buf4
     IF( ios == 0 )THEN
      WHERE( buf4 == undefbin )
       undefmask(:,:,n)=.false.
      ENDWHERE
       var(:,:,n) = REAL(buf4,r_size)
     ELSE !We reached the end of the file.
       IF( n==1)THEN
         WRITE(*,*)"[Warning]: Error reading, this might be a new file "
         WRITE(*,*)filename
         var=0.0d0     
         undefmask=.false.
       ELSE
         WRITE(*,*)"[Error]: Unexpectively reached end of file"
         WRITE(*,*)filename
         STOP
       ENDIF
     ENDIF
  END DO

  CLOSE(iunit)

  RETURN
END SUBROUTINE read_grd

!-- Write a grid file -------------------------------------------------
SUBROUTINE write_grd(filename,nx,ny,nz,var,undefmask,undefbin)
  IMPLICIT NONE
  INTEGER, INTENT(IN)     :: nx,ny,nz
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_sngl),INTENT(IN) :: var(nx,ny,nz)
  LOGICAL     ,INTENT(IN) :: undefmask(nx,ny,nz)
  REAL(r_sngl),INTENT(IN) :: undefbin
  REAL(r_sngl) :: buf4(nx,ny) 
  INTEGER :: iunit,iolen
  INTEGER :: k,n,reclength

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


  iunit=55
  INQUIRE(IOLENGTH=iolen) iolen
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength, &
             CONVERT=output_endian )

  DO n=1,nz
   buf4 = REAL(var(:,:,n),r_sngl)
   WHERE( .NOT. undefmask(:,:,n) )
    buf4=undefbin
   ENDWHERE     
   
   WRITE(iunit,rec=n) buf4
  ENDDO

  CLOSE(iunit)

  RETURN
END SUBROUTINE write_grd

SUBROUTINE generate_sample(sampledindex,n)
!Resample randomly picking with substitution.
IMPLICIT NONE
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(OUT) :: sampledindex(n)
REAL(r_size) :: randomloc(n)
INTEGER :: intind(n) , ii

 CALL com_rand(n,randomloc)
 sampledindex=INT(1+randomloc*(n-1))

END SUBROUTINE generate_sample


END MODULE covariance_matrix_tools
