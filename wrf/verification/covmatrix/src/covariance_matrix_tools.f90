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
  integer, parameter :: max_npoints = 20 , max_vars = 20
  integer :: npoints
  real(r_size) :: plon(max_npoints)  , plat(max_npoints) , plev(max_npoints) 
  character(50) :: pvarname(max_npoints)
  real(r_size) :: dep(max_npoints) , error(max_npoints)
  LOGICAL :: bootstrap=.true.
  INTEGER :: bootstrap_samples=20
  LOGICAL :: computehistogram=.true.
  INTEGER :: max_histogram=50            !Number of bins

  LOGICAL  :: computemoments=.true.
  INTEGER  :: max_moments=4

  LOGICAL      :: smoothcov=.false.
  REAL(r_size) :: smoothcovlength=1.0d5  !Lanczos filter length scale in the same unit as dx.
  REAL(r_size) :: smoothdx=1.0d3         !Data resolution 

  LOGICAL :: computecovar=.true.
  INTEGER :: skipx=1,skipy=1,skipz=1
  INTEGER :: nignore=0
  INTEGER , PARAMETER :: max_nignore = 20
  character(50) :: ignorevarname(max_nignore)

  LOGICAL :: computeindex=.true.       !Compute or not covariance and correlation index.
  INTEGER :: delta = 1
  INTEGER :: skip  = 1                 !Skip grid points in covariance strenght computation.

  CHARACTER(clen) :: BASE_VARS(max_vars)  !This is a list of reference or observed variables.
  CHARACTER(clen) :: COV_VARS(max_vars)   !This is the variables whose correlation / covariance with the observed variables will be measured.
  INTEGER :: NBASE_VARS=0 , NCOV_VARS=0   !Number of base vars and cov vars.
  REAL(r_sngl)    :: tr_rain = 5.0e-4 , tr_norain = 1.0e-5 !Condensate thresholds to separate rainy grid points.

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
                     smoothdx , smoothcovlength , computemoments ,   &
                     max_moments , computeindex , delta , computecovar, &
                     computehistogram , max_histogram , base_vars ,  &
                     cov_vars , nbase_vars , ncov_vars , tr_rain , tr_norain

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
  INTEGER :: k,n,reclength,ios ,ii,jj
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

SUBROUTINE write_profile(filename,profile,delta,nz,ncat)
IMPLICIT NONE
INTEGER,      INTENT(IN) :: delta , nz , ncat
CHARACTER(*), INTENT(IN) :: filename
REAL(r_sngl), INTENT(IN) :: profile(delta,nz,ncat)
INTEGER                  :: ii,jj,kk,iunit=66,reclen
 
 reclen=delta * 20

 OPEN(iunit,FILE=filename,FORM='formatted',RECL=reclen)
 WRITE(iunit,*)delta,nz,ncat
 DO jj = 1, ncat
  DO kk=1,nz
    WRITE(iunit,*)(profile(ii,kk,jj),ii=1,delta)
  ENDDO
 ENDDO     
 
 CLOSE(iunit)

END SUBROUTINE write_profile

!Compute the firtst N moments of the PDF centered around the mean.
SUBROUTINE compute_moments(ensemble,nx,ny,nz,nbv,nmoments,moments,undefmask,undefbin)
IMPLICIT NONE
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
INTEGER, INTENT(IN)       :: nmoments     !Number of moments to be computed.
LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(IN)  :: undefbin
REAL(r_sngl), INTENT(OUT) :: moments(nx,ny,nz,nmoments)
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
       IF ( im == 1 )THEN
        moments(ii,jj,kk,im) = SUM( tmp )/REAL( nbv , r_sngl ) 
       ELSE
        moments(ii,jj,kk,im) = SUM( tmp )/REAL( nbv -1 , r_sngl )
       ENDIF
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


SUBROUTINE compute_histogram(ensemble,nx,ny,nz,nbv,nbins,undefmask,undef,varmin,varmax,histogram)
IMPLICIT NONE
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv           !Ensemble dimensions.
REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
INTEGER, INTENT(IN)       :: nbins                  !Number of moments to be computed.
LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz)    !Valid grids.
LOGICAL                   :: tmpundefmask(nx,ny,nz)
REAL(r_sngl), INTENT(IN)  :: undef
INTEGER*2, INTENT(OUT)    :: histogram(nx,ny,nz,nbins)
INTEGER                   :: tmpindex , reclength , iunit
REAL(r_sngl),INTENT(OUT)  :: varmin(nx,ny,nz),varmax(nx,ny,nz)
REAL(r_sngl)              :: tmp(nbv) , dvar
INTEGER                   :: ii , jj , kk , im , irec

histogram=0

tmpundefmask=.true. !Do not mask output fields.

!$OMP PARALLEL DO PRIVATE(ii,jj,kk)
  DO ii=1,nx
   DO jj=1,ny
    DO kk=1,nz
       varmin(ii,jj,kk)=MINVAL( ensemble(ii,jj,kk,:) )
       varmax(ii,jj,kk)=MAXVAL( ensemble(ii,jj,kk,:) )
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,im,tmp,dvar,tmpindex)
  DO ii=1,nx
   DO jj=1,ny
    DO kk=1,nz
      if( .not. undefmask(ii,jj,kk) ) then

       dvar=( varmax(ii,jj,kk)-varmin(ii,jj,kk) )/real(nbins,r_sngl)

       IF( dvar > 0 )THEN
         DO im=1,nbv
           tmpindex=floor( ( ensemble(ii,jj,kk,im) - varmin(ii,jj,kk) )/dvar ) + 1
           IF( tmpindex > nbins )tmpindex=nbins
           IF( tmpindex < 1     )tmpindex=1
           !WRITE(*,*)histogram(ii,jj,kk,tmpindex),tmpindex
           histogram(ii,jj,kk,tmpindex)=histogram(ii,jj,kk,tmpindex) + 1
         ENDDO
       ENDIF
      else
         histogram(ii,jj,kk,:)=undef
      endif
    ENDDO
   ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_histogram


SUBROUTINE compute_kld(ensemble,nx,ny,nz,nbv,nbins,kld,undefmask,undef)
IMPLICIT NONE
!REAL(r_sngl),PARAMETER :: pi=3.1415926535e0
INTEGER, INTENT(IN)       :: nbins
INTEGER, INTENT(IN)       :: nx,ny,nz,nbv !Ensemble dimensions.
REAL(r_sngl), INTENT(IN)  :: ensemble(nx,ny,nz,nbv) !Ensemble data.
REAL(r_sngl), INTENT(IN)  :: undef        !Missing data code.
LOGICAL, INTENT(IN)       :: undefmask(nx,ny,nz) !Valid grids.
REAL(r_sngl), INTENT(OUT) :: kld(nx,ny,nz)
REAL(r_sngl)              :: tmp(nbv)
INTEGER                   :: ii , jj , kk , ib
INTEGER*2                 :: histogram(nx,ny,nz,nbins)
REAL(r_sngl)              :: varmin(nx,ny,nz) , varmax(nx,ny,nz)
REAL(r_sngl)              :: moments(nx,ny,nz,2)

REAL(r_sngl)              :: p(nbins) , q(nbins) , dvar , gconst , gconst2 , x

kld=0.0e0

 !Compute the histogram.
 CALL compute_histogram(ensemble,nx,ny,nz,nbv,nbins,undefmask,undef,varmin,varmax,histogram)

 !Compute the mean and standard deviation.
 CALL compute_moments(ensemble,nx,ny,nz,nbv,2,moments,undefmask,undef)

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,ib,x,dvar,p,q,gconst,gconst2)
  DO ii=1,nx
    DO jj=1,ny
     DO kk=1,nz
      if( .not. undefmask(ii,jj,kk) )then
       dvar=( varmax(ii,jj,kk)-varmin(ii,jj,kk) )/real(nbins,r_sngl)
       gconst=sqrt(1./(pi*2*moments(ii,jj,kk,2)))
       gconst2=(1./(2*moments(ii,jj,kk,2)))
       kld(ii,jj,kk)=0.0e0
       DO ib=1,nbins
         x=varmin(ii,jj,kk)+(ib-0.5)*dvar
         p(ib)=real(histogram(ii,jj,kk,ib),r_sngl)/real(nbv,r_sngl)            !Ensemble pdf
         q(ib)=dvar*gconst*exp(-((x-moments(ii,jj,kk,1))**2)*gconst2)          !Gaussian pdf.

         if( p(ib) > 0 )then
           kld(ii,jj,kk)=kld(ii,jj,kk)+p(ib)*( log(p(ib)) - log(q(ib) ) )
         endif
       ENDDO
      else
       kld(ii,jj,kk)=undef
      endif
     ENDDO
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE compute_kld

!Compute covariance strhength and covariance profiles among two variables.

SUBROUTINE covariance_strenght(var1_ens,var2_ens,var_qhyd,covst,corrst,corrstdist,cov_profile,corr_profile, & 
                               num_profile,nx,ny,nz,nbv,delta,skip_cov,undefmask,undefbin)

IMPLICIT NONE
INTEGER, INTENT(IN) :: nx,ny,nz,nbv,delta,skip_cov
REAL(r_sngl), INTENT(IN) :: var1_ens(nx,ny,nz,nbv)  , undefbin  !Base variable
REAL(r_sngl), INTENT(IN) :: var2_ens(nx,ny,nz,nbv)              !Covariance variable
REAL(r_sngl), INTENT(IN) :: var_qhyd(nx,ny,nz)                  !Hydrometeor content.
LOGICAL ,     INTENT(IN) :: undefmask(nx,ny,nz) 
REAL(r_sngl), INTENT(OUT):: covst(nx,ny,nz),corrst(nx,ny,nz),corrstdist(nx,ny,nz) !Covariance and correlation strhenght index.
REAL(r_sngl)             :: var1_m(nx,ny,nz,2) , var2_m(nx,ny,nz,2)
INTEGER :: imin , imax , jmin , jmax , ii , jj , kk , iii , jjj , contador
REAL(r_sngl)             :: pensemble(nbv) , tmpcov , tmpcorr , tmpdist
REAL(r_sngl)             :: tmp_cov_profile(delta) , tmp_corr_profile(delta) !, dist_profile(delta)
INTEGER                  :: tmp_num_profile(delta) , current_index
REAL(r_sngl),INTENT(OUT) :: cov_profile(delta,nz,2) , corr_profile(delta,nz,2)
INTEGER     ,INTENT(OUT) :: num_profile(delta,nz,2)

covst=0
corrst=0
corrstdist=0

cov_profile=0.0e0     !Note that if OPENMP is used this variable is initialized to 0
corr_profile=0.0e0    !Note that if OPENMP is used this variable is initialized to 0
num_profile=0.0e0     !Note that if OPENMP is used this variable is initialized to 0

!tr_norain=1.0e-5
!tr_rain  =1.0e-3

CALL compute_moments(var1_ens,nx,ny,nz,nbv,2,var1_m,undefmask,undefbin)
CALL compute_moments(var2_ens,nx,ny,nz,nbv,2,var2_m,undefmask,undefbin)

var1_m(:,:,:,2)=sqrt(var1_m(:,:,:,2)) !Get standard deviation.
var2_m(:,:,:,2)=sqrt(var2_m(:,:,:,2)) !Get standard deviation.

!$OMP PARALLEL DO PRIVATE(imin,imax,jmin,jmax,ii,jj,kk,iii,jjj,pensemble,tmpcov,tmpcorr,contador,tmpdist, &
!$OMP tmp_cov_profile,tmp_corr_profile,tmp_num_profile,current_index)                                     &
!$OMP reduction(+:cov_profile,corr_profile,num_profile)

DO ii=1,nx
 !WRITE(*,*)ii
 DO jj=1,ny
  DO kk=1,nz
   IF( undefmask(ii,jj,kk) )THEN

    imin=ii-delta
    imax=ii+delta
    if( imin < 1 )imin=1
    if( imax > nx)imax=nx
    jmin=jj-delta
    jmax=jj+delta
    if( jmin < 1 )jmin=1
    if( jmax > ny)jmax=ny

    pensemble=var1_ens(ii,jj,kk,:)

    tmp_cov_profile=0.0e0
    tmp_corr_profile=0.0e0
    !dist_profile=0.0e0
    tmp_num_profile=0


    DO iii=imin,imax,skip_cov
     DO jjj=jmin,jmax,skip_cov

      tmpdist= sqrt(real(iii-ii,r_sngl)**2 + real(jjj-jj,r_sngl)**2) !/real(delta,r_sngl)

      IF( undefmask(iii,jjj,kk) .and. tmpdist <= real(delta,r_size) )THEN
   
       IF( var1_m(ii,jj,kk,2) > 0 .and. var2_m(iii,jjj,kk,2) > 0 )THEN

        CALL com_covar_sngl(nbv,var2_ens(iii,jjj,kk,:),pensemble,tmpcov)
  
        tmpcorr=tmpcov/( var1_m(ii,jj,kk,2) * var2_m(iii,jjj,kk,2) ) 

        current_index=FLOOR( tmpdist ) + 1

        IF( current_index >= delta) current_index=delta

        tmp_cov_profile(current_index)=tmp_cov_profile(current_index)+ABS(tmpcov) 
  
        tmp_corr_profile(current_index)=tmp_corr_profile(current_index)+ABS(tmpcorr)

        !dist_profile(current_index)=dist_profile(current_index)+tmpdist*ABS(tmpcorr)

        tmp_num_profile(current_index)=tmp_num_profile(current_index) + 1
    
       ENDIF
  
      ENDIF
       
     ENDDO
    ENDDO

    IF( var_qhyd(ii,jj,kk) >= tr_rain )THEN
      !Accumulate the profiles for rainy grid points
      cov_profile(:,kk,1) =cov_profile(:,kk,1) + tmp_cov_profile 
      corr_profile(:,kk,1)=corr_profile(:,kk,1) + tmp_corr_profile
      num_profile(:,kk,1)=num_profile(:,kk,1) + tmp_num_profile
    ELSEIF( var_qhyd(ii,jj,kk) <= tr_norain )THEN
      !Accumulate the profiles for no rain grid points
      cov_profile(:,kk,2) =cov_profile(:,kk,2) + tmp_cov_profile
      corr_profile(:,kk,2)=corr_profile(:,kk,2) + tmp_corr_profile
      num_profile(:,kk,2)=num_profile(:,kk,2) + tmp_num_profile
    ENDIF
  
    contador=0 
    DO iii=1,delta
        IF(  tmp_num_profile(iii) >= 1 )THEN
          tmp_cov_profile(iii) = tmp_cov_profile(iii) / tmp_num_profile(iii)
          tmp_corr_profile(iii) = tmp_corr_profile(iii) / tmp_num_profile(iii)

          covst(ii,jj,kk) = covst(ii,jj,kk) + tmp_cov_profile(iii) 
          corrst(ii,jj,kk) = corrst(ii,jj,kk) + tmp_corr_profile(iii)
          corrstdist(ii,jj,kk) = corrstdist(ii,jj,kk) + tmp_corr_profile(iii) * iii
 
          contador=contador + 1
        ENDIF
    ENDDO

    IF( contador >= 0 )THEN
      corrstdist(ii,jj,kk)=corrstdist(ii,jj,kk)/corrst(ii,jj,kk)
      covst(ii,jj,kk)=covst(ii,jj,kk)/REAL(contador,r_sngl)
      corrst(ii,jj,kk)=corrst(ii,jj,kk)/REAL(contador,r_sngl)
    ELSE
      covst(ii,jj,kk)=undefbin
      corrst(ii,jj,kk)=undefbin
      corrstdist(ii,jj,kk)=undefbin
    ENDIF
   
   ELSE
 
    covst(ii,jj,kk)=undefbin
    corrst(ii,jj,kk)=undefbin
    corrstdist(ii,jj,kk)=undefbin

   ENDIF

  ENDDO
 ENDDO
ENDDO

!$OMP END PARALLEL DO

  DO ii=1,delta
   DO kk=1,nz
    DO jj=1,2
     IF( num_profile(ii,kk,jj) >= 1 )THEN
       cov_profile(ii,kk,jj)=cov_profile(ii,kk,jj) / real( num_profile(ii,kk,jj) , r_sngl )
       corr_profile(ii,kk,jj)=corr_profile(ii,kk,jj) / real( num_profile(ii,kk,jj) , r_sngl )
     ELSE
       cov_profile(ii,kk,jj)=undefbin
       corr_profile(ii,kk,jj)=undefbin
     ENDIF
    ENDDO
   ENDDO
  ENDDO

END SUBROUTINE covariance_strenght

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





 

  

  
  




















