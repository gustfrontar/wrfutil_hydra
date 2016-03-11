!===============================================================================
! covariance matrix: main program
!-------------------------------------------------------------------------------
!===============================================================================

program covariance_matrix
!-------------------------------------------------------------------------------

  use common
  use covariance_matrix_tools
  implicit none

  !Original grid files and variables
  character(100) :: fcstfile = 'fcst00000.grd'
  character(100) :: covfile  = 'covmat_xXXXyXXXzXXX.grd'
  character(100) :: oimfile  = 'impact_xXXXyXXXzXXX.grd'
  character(100) :: bssfile  = 'bssprd_xXXXyXXXzXXX.grd' !Bootstrap spread
  character(100) :: bsmfile  = 'bsmean_xXXXyXXXzXXX.grd' !Bootstrap mean
  character(100) :: ctl_file = 'input.ctl'

  real(r_sngl) ,ALLOCATABLE :: ensemble(:,:,:,:) !nlon,nlat,nv,nens
  real(r_sngl) ,ALLOCATABLE :: pensemble(:)        !nens
  real(r_sngl) :: stdens

  !Covariance between and observation at a certain location and the 
  !state vector. 
  real(r_sngl) ,ALLOCATABLE :: covar(:,:,:)        !nlon,nlat,nv
  !Obs impact how much impact an observation located at a certain location
  !will have on the state vector.
  real(r_sngl) ,ALLOCATABLE :: obsimpact(:,:,:)    !nlon,nlat,nv
  real(r_sngl) ,ALLOCATABLE :: bssprd(:,:,:) !Bootstrap covariance spread
  real(r_sngl) ,ALLOCATABLE :: bsmean(:,:,:) !Bottstrap covariance mean
  real(r_sngl) ,ALLOCATABLE :: tmpsample(:,:),tmp(:,:)  !Temporary buffer for data resampling.
  real(r_sngl) :: cov , covmean , covsprd


  LOGICAL , ALLOCATABLE :: undefmask(:,:,:) , totalundefmask(:,:,:)

  integer :: i, ii , jj , kk  , is
 
  real(r_size) ri(1) , rj(1)

  integer :: iunit, iolen, irec , ip , gridi , gridj , gridk

!-------------------------------------------------------------------------------

!===============================================================================

  call read_namelist !ensemble size and general options.

  call parse_ctl(ctl_file) !Get number of variables, X, Y and Z

  ALLOCATE( ensemble(nlon,nlat,nv,nbv) )
  ALLOCATE( pensemble(nbv) )  
  ALLOCATE( covar(nlon,nlat,nv) )
  ALLOCATE( obsimpact(nlon,nlat,nv) )

  ALLOCATE( undefmask(nlon,nlat,nv) , totalundefmask(nlon,nlat,nv) )

  if( bootstrap )then
    ALLOCATE( bssprd(nlon,nlat,nv) , bsmean(nlon,nlat,nv) )
    ALLOCATE( tmpsample(nbv,2) , tmp(nbv,2) )
  endif

  !First read all the ensemble and store it in memory (if this is not possible
  !then the code should be modified).

   totalundefmask=.true.

DO i=1,nbv

     WRITE( fcstfile(5:9), '(I5.5)' )i
     call read_grd(fcstfile,nlon,nlat,nv,ensemble(:,:,:,i),undefmask)

     WHERE( .NOT. undefmask )
        totalundefmask=.false.
     END WHERE

ENDDO ![End do over ensemble members]

   !Loop over points.

DO ip=1,npoints
     write(*,*)"Processing point ",ip," of ",npoints
     write(*,*)"Lat ",plat(ip)," Lon ",plon(ip)," Var ",plev(ip)

     !Get the closest grid point to the selected location.
     CALL com_pos2ij(1,nlon,nlat,lon,lat,1,plon(ip),plat(ip),ri(1),rj(1))
     gridi=NINT(ri(1))
     gridj=NINT(rj(1))
     gridk=NINT(plev(ip))

     if( gridi > nlon .or. gridi < 1 .or. gridj > nlat .or. gridj < 1 .or. gridk > nv )then
       write(*,*)"[Warning]: Requested coordinates are outside current grid."
       cycle
     endif

     !Get ensemble values for this grid location.
     pensemble=ensemble(gridi,gridj,gridk,:) 
 
     !Compute spread at the selected location/variable
     CALL com_stdev_sngl(nbv,pensemble,stdens)  
 
 
     !Compute the covariance between this location/variable and the rest
     !of the state space.

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,cov)
     DO ii = 1,nlon
      DO jj = 1,nlat
       DO kk = 1,nv
        IF( totalundefmask(ii,jj,kk) ) then
        CALL com_covar_sngl(nbv,ensemble(ii,jj,kk,:),pensemble,cov)
        covar(ii,jj,kk)=REAL(cov,r_sngl)
        obsimpact(ii,jj,kk)=covar(ii,jj,kk) * dep(ip) /( stdens + error(ip) )
        ENDIF
       ENDDO
      ENDDO
     ENDDO 
!$OMP END PARALLEL DO


   if( bootstrap )then
     write(*,*)"Using bootstrap, with ",bootstrap_samples," samples"
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,is,cov,tmp,tmpsample,covmean,covsprd)
     DO ii = 1,nlon
      tmp(:,1)=pensemble
      DO jj = 1,nlat
       DO kk = 1,nv
        IF( totalundefmask(ii,jj,kk) ) then
           covmean=0.0d0
           covsprd=0.0d0
           DO is = 1,bootstrap_samples
             tmp(:,2)=ensemble(ii,jj,kk,:)
             !tmpsample=tmp
             CALL generate_sample(tmp,tmpsample,nbv,2)
             CALL com_covar_sngl(nbv,tmpsample(:,1),tmpsample(:,2),cov)
             covmean=covmean+cov
             covsprd=covsprd+cov**2
           ENDDO
           bsmean(ii,jj,kk)=covmean/REAL(bootstrap_samples,r_sngl)
           bssprd(ii,jj,kk)=sqrt( covsprd/REAL(bootstrap_samples,r_sngl)-( covmean/REAL(bootstrap_samples,r_sngl) )**2 )
        ENDIF
       ENDDO
      ENDDO
     ENDDO
!$OMP END PARALLEL DO
   endif


   !Write the results to a file.
   WRITE(covfile(9 :11),'(I3.3)')gridi
   WRITE(covfile(13:15),'(I3.3)')gridj
   WRITE(covfile(17:19),'(I3.3)')gridk

   CALL write_grd(covfile,nlon,nlat,nv,covar,totalundefmask)

  if( bootstrap ) then

   WRITE(bsmfile(9 :11),'(I3.3)')gridi
   WRITE(bsmfile(13:15),'(I3.3)')gridj
   WRITE(bsmfile(17:19),'(I3.3)')gridk

   CALL write_grd(bsmfile,nlon,nlat,nv,bsmean,totalundefmask)

   WRITE(bssfile(9 :11),'(I3.3)')gridi
   WRITE(bssfile(13:15),'(I3.3)')gridj
   WRITE(bssfile(17:19),'(I3.3)')gridk

   CALL write_grd(bssfile,nlon,nlat,nv,bssprd,totalundefmask)

  endif

   WRITE(oimfile(9 :11),'(I3.3)')gridi
   WRITE(oimfile(13:15),'(I3.3)')gridj
   WRITE(oimfile(17:19),'(I3.3)')gridk
 
   CALL write_grd(oimfile,nlon,nlat,nv,obsimpact,totalundefmask)

ENDDO ![End do over points]



end program covariance_matrix
!===============================================================================
