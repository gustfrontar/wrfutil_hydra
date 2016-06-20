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
!  USE common_smooth2d
!  USE ifport
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

  LOGICAL  :: computemoments=.true.
  INTEGER  :: max_moments=4

  LOGICAL      :: smoothcov=.false.
  REAL(r_size) :: smoothcovlength=1.0d5  !Lanczos filter length scale in the same unit as dx.
  REAL(r_size) :: smoothdx=1.0d3         !Data resolution 

  INTEGER :: skipx=1,skipy=1,skipz=1
  INTEGER :: nignore=0
  INTEGER , PARAMETER :: max_nignore = 20
  character(50) :: ignorevarname(max_nignore)
  character(50) :: inputendian='little_endian'
  character(50) :: outputendian='big_endian'

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

!In the current version PLEV represents not only the level but also the
!variable. 

NAMELIST / GENERAL / nbv , npoints , plon , plat , pvarname , plev,   &
                     dep , error , bootstrap , bootstrap_samples  ,   &
                     nignore , ignorevarname , skipx , skipy , skipz, &
                     inputendian , outputendian , smoothcov ,         &
                     smoothdx , smoothcovlength , computemoments ,    &
                     max_moments

plon=undef
plat=undef
plev=undef
dep=undef
error=undef


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
  INTEGER :: k,n,reclength,ios
  LOGICAL :: file_exist
  

  undefmask=.true.

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nx * ny


 iunit=33
 INQUIRE(FILE=filename,EXIST=file_exist)
  IF(file_exist) THEN
   OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=inputendian)
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
  OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=outputendian)

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

!Compute the firtst N moments of the PDF centered around the mean.
SUBROUTINE compute_moments(ensemble,nx,ny,nz,nbv,nmoments,undefmask,undefbin)
IMPLICIT NONE
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
INTEGER, INTENT(IN)       :: nmoments     !Number of moments to be computed.
LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(IN)  :: undefbin
REAL(r_sngl)              :: moments(nx,ny,nz,nmoments)
REAL(r_sngl)              :: tmp(nbv)
INTEGER                   :: ii , jj , kk , im
CHARACTER(20)             :: moment_out='momentXXX.grd'

moments=0.0e0

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,im,tmp)
  DO ii=1,nx
   DO im=1,nmoments
    DO jj=1,ny
     DO kk=1,nz
       tmp=( ( ensemble(ii,jj,kk,:)-moments(ii,jj,kk,1) ) )**im
       CALL com_mean_sngl(nbv,tmp,moments(ii,jj,kk,im))
     ENDDO
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO
  
 DO im=1,nmoments
   WRITE(moment_out(7:9),'(I3.3)')im
   CALL write_grd(moment_out,nx,ny,nz,moments(:,:,:,im),undefmask,undefbin)
 ENDDO

END SUBROUTINE compute_moments


!Perform a 2D smoothing of the ensemble using a Lanczos filter.
SUBROUTINE smooth_2d(mydata,nx,ny,nz,dx,xfs,undefmask)
IMPLICIT NONE
INTEGER      , INTENT(IN)       :: nx,ny,nz               !Grid dimensions
REAL(r_sngl) , INTENT(INOUT)    :: mydata(nx,ny,nz)       !Model data
REAL(r_size) , INTENT(IN)       :: dx      !Grid resolution.
REAL(r_size) , INTENT(IN)       :: xfs     !Filter scale in x and y.
LOGICAL      , INTENT(IN)       :: undefmask(nx,ny,nz)
INTEGER                         :: kk
INTEGER                         :: integermask(nx,ny,nz)
INTEGER                         :: filter_size_x
REAL(r_size)                    :: tmpdata(nx,ny,nz)

integermask=0
WHERE( undefmask )
  integermask=1
END WHERE

filter_size_x = NINT( xfs / dx )     

tmpdata=REAL(mydata,r_size)

!$OMP PARALLEL DO PRIVATE(kk)
DO kk=1,nz  
    CALL lanczos_2d(tmpdata(:,:,kk),tmpdata(:,:,kk),integermask(:,:,kk),filter_size_x,nx,ny)
ENDDO
!$OMP END PARALLEL DO

mydata=REAL(tmpdata,r_sngl)

END SUBROUTINE smooth_2d


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



SUBROUTINE lanczos_2d(inputvar2d,outputvar2d,mask,lambda,nx,ny)
!This is a single-routine lanczos filter. This routine has been prepared in order to be used
!with openmp (common_smooth_2d module shares variables among different routines and can not be
!safely used with openmp)
  IMPLICIT NONE
  INTEGER,INTENT(IN)      :: nx,ny 
  REAL(r_size) , INTENT(IN) :: inputvar2d(nx,ny) 
  REAL(r_size) , INTENT(OUT):: outputvar2d(nx,ny)
  INTEGER      , INTENT(IN) :: mask(nx,ny)
  INTEGER      , INTENT(IN) :: lambda
  REAL(r_size)                         :: fn, rspval        ! cutoff freq/wavelength, spval
  REAL(r_size), DIMENSION(0:2*lambda)  :: dec, de           ! weight in r8, starting index 0:nband
  INTEGER                              :: ji, jj, jmx, jkx, jk, jt, jvar !  dummy loop index
  INTEGER                              :: ik1x, ik2x, ikkx
  INTEGER                              :: ifrst=0
  INTEGER                              :: inxmin, inxmaxi
  INTEGER                              :: inymin, inymaxi
  REAL(r_size), DIMENSION(nx,ny)       :: dl_tmpx, dl_tmpy
  REAL(r_size)                         :: dl_yy, dl_den, dl_pi, dl_ey, dl_coef

  !   PRINT *,'        ncut     : number of grid step to be filtered'
  !  remark: for a spatial filter, fn=dx/lambda where dx is spatial step, lamda is cutting wavelength

  fn    = 1./lambda

  !Filter init -------------------------------------  
 
    dl_pi   = ACOS(-1.d0)
    dl_coef = 2.0d0*dl_pi*fn

    de(0) = 2.d0*fn
    DO  ji=1,2*lambda
       de(ji) = SIN(dl_coef*ji)/(dl_pi*ji)
    END DO
    !
    dec(0) = 2.d0*fn
    DO ji=1,2*lambda
       dl_ey   = dl_pi*ji/(2*lambda)
       dec(ji) = de(ji)*SIN(dl_ey)/dl_ey
    END DO

  !End of filter init --------------------------------
  
  !FILTER-----------------------------------------------

   IF ( lambda /= 0 )THEN

    inxmin   =  2*lambda
    inxmaxi  =  nx-2*lambda+1
    inymin   =  2*lambda
    inymaxi  =  ny-2*lambda+1


    DO jj=1,ny
       DO  jmx=1,nx
          ik1x = -2*lambda
          ik2x =  2*lambda
          !
          IF (jmx <= inxmin ) ik1x = 1-jmx
          IF (jmx >= inxmaxi) ik2x = nx-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (mask(jkx+jmx,jj)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*inputvar2d(jkx+jmx,jj)
             END IF
          END DO
          !
          dl_tmpx(jmx,jj)=dl_yy/dl_den
       END DO
    END DO

    DO ji=1,nx
       DO  jmx=1,ny
          ik1x = -2*lambda
          ik2x =  2*lambda
          !
          IF (jmx <= inymin ) ik1x = 1-jmx
          IF (jmx >= inymaxi) ik2x = ny-jmx
          !
          dl_yy  = 0.d0
          dl_den = 0.d0
          !
          DO jkx=ik1x,ik2x
             ikkx=ABS(jkx)
             IF (mask(ji,jkx+jmx)  ==  1) THEN
                dl_den = dl_den + dec(ikkx)
                dl_yy  = dl_yy  + dec(ikkx)*dl_tmpx(ji,jkx+jmx)
             END IF
          END DO
          outputvar2d(ji,jmx)=0.
          IF (dl_den /=  0.) outputvar2d(ji,jmx) = dl_yy/dl_den
       END DO
    END DO
    !


   ENDIF 

  !END OF FILTER-----------------------------------------


   IF ( lambda == 0 ) outputvar2d=inputvar2d

   WHERE( mask == 0 )outputvar2d = inputvar2d


  END SUBROUTINE lanczos_2d

 



END MODULE covariance_matrix_tools





 

  

  
  




















