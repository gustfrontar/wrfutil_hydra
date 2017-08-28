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
  character(100) :: hstfile  = 'histgr_xXXXyXXXzXXX.grd' !Histogram.
  !character(100) :: covindexfile =      'covariance_index.grd'
  !character(100) :: corrindexfile =     'correlation_index.grd'
  character(100) :: tmpfilename 
  character(100) :: ctl_file = 'input.ctl'

  real(r_sngl) ,ALLOCATABLE :: ensemble(:,:,:,:) !nlon,nlat,nfields,nens
  real(r_sngl) ,ALLOCATABLE :: pensemble(:)      !nens
  integer      ,ALLOCATABLE :: sampleindex(:)    !bootstrap sample index array
  real(r_sngl) :: stdens

  !Covariance between and observation at a certain location and the 
  !state vector. 
  real(r_sngl) ,ALLOCATABLE :: var1(:,:,:) , var2(:,:,:) , mean_cond(:,:,:)
  real(r_sngl) ,ALLOCATABLE :: covar(:,:,:)        !nlon,nlat,nfields
  real(r_sngl) ,ALLOCATABLE :: localcov(:,:,:)     !nlon,nlat,nfields

  real(r_sngl) ,ALLOCATABLE :: cov_profile(:,:,:) , corr_profile(:,:,:) 
  integer      ,ALLOCATABLE ::  num_profile(:,:,:)


  !Obs impact how much impact an observation located at a certain location
  !will have on the state vector.
  real(r_sngl) ,ALLOCATABLE :: obsimpact(:,:,:)         !nlon,nlat,nfields
  real(r_sngl) ,ALLOCATABLE :: bssprd(:,:,:)            !Bootstrap covariance spread
  real(r_sngl) ,ALLOCATABLE :: bsmean(:,:,:)            !Bottstrap covariance mean
  real(r_sngl) ,ALLOCATABLE :: tmpsample(:,:),tmp(:,:)  !Temporary buffer for data resampling.
  real(r_sngl) ,ALLOCATABLE :: moments(:,:,:,:)         !PDF moments
  real(r_sngl) ,ALLOCATABLE :: covariance_index(:,:,:)      !Covariance strength index
  real(r_sngl) ,ALLOCATABLE :: correlation_index(:,:,:)     !Correlation strength index
  real(r_sngl) ,ALLOCATABLE :: correlationdist_index(:,:,:) !Correlation and distance index
  integer*2    ,ALLOCATABLE :: histogram(:,:,:,:)           !Histogram data.
  real(r_sngl) ,ALLOCATABLE :: varmin(:,:,:) , varmax(:,:,:) 
  real(r_sngl) :: cov , covmean , covsprd


  LOGICAL , ALLOCATABLE :: undefmask(:,:,:) , totalundefmask(:,:,:) , local_undefmask(:,:,:)

  integer :: i, j , ii , jj , kk  , is , ilev
 
  real(r_size) ri(1) , rj(1)

  integer :: iunit, iolen, irec , ip , gridi , gridj , gridk

  integer :: cond_i , var1_i , var2_i 
  integer :: cond_si , var1_si , var2_si 
  integer :: cond_ei , var1_ei , var2_ei
  integer :: local_nzc , local_nz1 , local_nz2

  integer :: reclength

!-------------------------------------------------------------------------------

!===============================================================================

  call read_namelist !ensemble size and general options.

  call parse_ctl(ctl_file , ctl ) !Get ctl info (varnames, levs, grid size, etc)

  ALLOCATE( ensemble(ctl%nlon,ctl%nlat,ctl%nfields,nbv) )
  ALLOCATE( pensemble(nbv) )  
  ALLOCATE( covar(ctl%nlon,ctl%nlat,ctl%nfields) )
  ALLOCATE( obsimpact(ctl%nlon,ctl%nlat,ctl%nfields) )

  ALLOCATE( undefmask(ctl%nlon,ctl%nlat,ctl%nfields) , totalundefmask(ctl%nlon,ctl%nlat,ctl%nfields) )

  covar=ctl%undefbin
  ensemble=0.0e0
  pensemble=0.0e0
  obsimpact=ctl%undefbin
  

!First read all the ensemble and store it in memory (if this is not possible
!then the code should be modified).

   totalundefmask=.true.

DO i=1,nbv

     WRITE( fcstfile(5:9), '(I5.5)' )i
     WRITE(*,*)"Reading file ",fcstfile
     call read_grd(fcstfile,ctl%nlon,ctl%nlat,ctl%nfields,ensemble(:,:,:,i),undefmask,ctl%undefbin)

     DO ii=1,ctl%nlon
      DO jj=1,ctl%nlat
       DO kk=1,ctl%nfields
        IF( .not. undefmask(ii,jj,kk) )totalundefmask(ii,jj,kk)=.false.
       ENDDO
      ENDDO
     ENDDO

     IF( smoothcov )THEN
         CALL smooth_2d(ensemble(:,:,:,i),ctl%nlon,ctl%nlat,ctl%nfields,smoothdx,smoothcovlength,undefmask)
         IF( i==1)THEN 
          CALL write_grd('testsmooth.grd',ctl%nlon,ctl%nlat,ctl%nfields,ensemble(:,:,:,i),totalundefmask,ctl%undefbin)
         ENDIF
     ENDIF

ENDDO ![End do over ensemble members]

!Apply variable and grid filter.
ilev=1
DO i=1,ctl%nfields
 DO j=1,nignore
    IF( ignorevarname(j) == ctl%varnameall(i) )THEN
      !This field won't be prcoessed
      write(*,*)"[Warning]: Field ",TRIM(ctl%varnameall(i))," at lev ",ctl%levall(i)," will be ignored"
      totalundefmask(:,:,i)=.false. 
    ENDIF
 ENDDO
      
 IF( i > 1 )THEN
   IF(  ctl%varnameall(i) /= ctl%varnameall(i-1) )ilev=1
 ENDIF
 IF( MOD(ilev-1,skipz) /= 0)THEN
   totalundefmask(:,:,i)=.false.
   write(*,*)"[Warning]: Field ",TRIM(ctl%varnameall(i))," at lev ",ctl%levall(i)," will be ignored"
 ENDIF
 ilev=ilev+1
ENDDO
DO i=1,ctl%nlon
  IF( MOD(i-1,skipx) /= 0)THEN
    totalundefmask(i,:,:)=.false.
    write(*,*)"[Warning]: Longitude ",i," will be ignored"
  ENDIF
ENDDO
DO i=1,ctl%nlat
  IF( MOD(i-1,skipy) /= 0)THEN
    totalundefmask(:,i,:)=.false.
    write(*,*)"[Warning]: Latitude ",i," will be ignored"
  ENDIF
ENDDO

!PDF Moments computation

IF( computemoments )THEN
   !Compute moments of the PDF at every grid point
   ALLOCATE( moments(ctl%nlon,ctl%nlat,ctl%nfields,max_moments) )
   CALL compute_moments(ensemble,ctl%nlon,ctl%nlat,ctl%nfields,nbv,max_moments,moments,totalundefmask,ctl%undefbin)
   DEALLOCATE( moments )
ENDIF

IF( computehistogram )THEN
   !Compute histograms
   WRITE(*,*)ctl%nlon,ctl%nlat,ctl%nfields,max_histogram
   ALLOCATE( histogram(ctl%nlon,ctl%nlat,ctl%nfields,max_histogram) )
   ALLOCATE( varmin(ctl%nlon,ctl%nlat,ctl%nfields) , varmax(ctl%nlon,ctl%nlat,ctl%nfields) )
   CALL compute_histogram(ensemble,ctl%nlon,ctl%nlat,ctl%nfields,nbv,max_histogram,totalundefmask,        &
                          ctl%undefbin,varmin,varmax,histogram)

  INQUIRE(IOLENGTH=reclength) reclength
  reclength = ctl%nlon * ctl%nlat !* reclength

  iunit=55
  tmpfilename='histogram.grd'
  !OPEN(iunit,FILE=histogram_out,FORM='unformatted',ACCESS='direct',RECL=reclength,CONVERT=outputendian)
  OPEN(iunit,FILE=tmpfilename,FORM='unformatted',ACCESS='sequential',CONVERT=outputendian)

  !irec=1
  DO ii=1,ctl%nfields
   DO jj=1,max_histogram
    WRITE(iunit)histogram(:,:,ii,jj)
   ENDDO
  ENDDO

  CLOSE(iunit)
 
  tmpfilename='minvar.grd'
  CALL write_grd(tmpfilename,ctl%nlon,ctl%nlat,ctl%nfields,varmin,totalundefmask,ctl%undefbin)
  tmpfilename='maxvar.grd'
  CALL write_grd(tmpfilename,ctl%nlon,ctl%nlat,ctl%nfields,varmax,totalundefmask,ctl%undefbin)

  DEALLOCATE( histogram , varmin , varmax )

ENDIF



!Covariance computation
!Loop over points.

IF( computecovar )THEN

DO ip=1,npoints
     write(*,*)"Processing point ",ip," of ",npoints
     write(*,*)"Lat ",plat(ip)," Lon ",plon(ip)," Var ",plev(ip)

     !Get the closest grid point to the selected location.
     CALL com_pos2ij(1,ctl%nlon,ctl%nlat,ctl%lon,ctl%lat,1,plon(ip),plat(ip),ri(1),rj(1))
     gridi=NINT(ri(1))
     gridj=NINT(rj(1))
     !Get the file index corresponding to the selected variable and level.
     CALL get_var_index(ctl,pvarname(ip),plev(ip),gridk)

     if( gridi > ctl%nlon .or. gridi < 1 .or. gridj > ctl%nlat .or. gridj < 1 .or. gridk > ctl%nfields )then
       write(*,*)"[Warning]: Requested coordinates are outside current grid."
       cycle
     endif
     if( gridk == 0 )then
       write(*,*)"[Warning]: Requested variable ",pvarname(ip)," was not found at level ",plev(ip)
       cycle
     endif

     !Get ensemble values for this grid location.
     pensemble=ensemble(gridi,gridj,gridk,:) 
 
     !Compute spread at the selected location/variable
     CALL com_stdev_sngl(nbv,pensemble,stdens)  

 
     !Compute the covariance between this location/variable and the rest
     !of the state space.

!$OMP PARALLEL DO PRIVATE(ii,jj,kk,cov)
     DO ii = 1,ctl%nlon
      DO jj = 1,ctl%nlat
       DO kk = 1,ctl%nfields
        IF( totalundefmask(ii,jj,kk) ) then
        CALL com_covar_sngl(nbv,ensemble(ii,jj,kk,:),pensemble,covar(ii,jj,kk))
        obsimpact(ii,jj,kk)=covar(ii,jj,kk) * dep(ip) /( stdens + error(ip) )
        ENDIF
       ENDDO
      ENDDO
     ENDDO 
!$OMP END PARALLEL DO



   if( bootstrap )then
     ALLOCATE( bssprd(ctl%nlon,ctl%nlat,ctl%nfields) , bsmean(ctl%nlon,ctl%nlat,ctl%nfields) )
     write(*,*)"Using bootstrap, with ",bootstrap_samples," samples"
     bsmean=ctl%undefbin
     bssprd=ctl%undefbin
!$OMP PARALLEL DO PRIVATE(ii,jj,kk,is,cov,covmean,covsprd,localcov,sampleindex) 
    DO is = 1,bootstrap_samples 
     ALLOCATE( localcov(ctl%nlon,ctl%nlat,ctl%nfields) , sampleindex(nbv) )
     CALL generate_sample(sampleindex,nbv)
     localcov=0.0d0
     DO ii = 1,ctl%nlon
      DO jj = 1,ctl%nlat
       DO kk = 1,ctl%nfields
        IF( totalundefmask(ii,jj,kk) ) then
          CALL com_covar_sngl_sample(nbv,pensemble,ensemble(ii,jj,kk,:),sampleindex,localcov(ii,jj,kk))
        ENDIF
       ENDDO
      ENDDO
     ENDDO
!$OMP CRITICAL     
      bsmean=bsmean+localcov
      bssprd=bssprd+localcov**2
!$OMP END CRITICAL
     DEALLOCATE(localcov,sampleindex)
   ENDDO
!$OMP END PARALLEL DO
   bsmean=bsmean/REAL(bootstrap_samples,r_sngl)
   bssprd=sqrt( bssprd/REAL(bootstrap_samples,r_sngl)-( bsmean )**2 )
   endif


   !Write the results to a file.
   WRITE(covfile(9 :11),'(I3.3)')gridi
   WRITE(covfile(13:15),'(I3.3)')gridj
   WRITE(covfile(17:19),'(I3.3)')gridk

   CALL write_grd(covfile,ctl%nlon,ctl%nlat,ctl%nfields,covar,totalundefmask,ctl%undefbin)

  if( bootstrap ) then

   WRITE(bsmfile(9 :11),'(I3.3)')gridi
   WRITE(bsmfile(13:15),'(I3.3)')gridj
   WRITE(bsmfile(17:19),'(I3.3)')gridk

   CALL write_grd(bsmfile,ctl%nlon,ctl%nlat,ctl%nfields,bsmean,totalundefmask,ctl%undefbin)

   WRITE(bssfile(9 :11),'(I3.3)')gridi
   WRITE(bssfile(13:15),'(I3.3)')gridj
   WRITE(bssfile(17:19),'(I3.3)')gridk

   CALL write_grd(bssfile,ctl%nlon,ctl%nlat,ctl%nfields,bssprd,totalundefmask,ctl%undefbin)

   deallocate( bsmean , bssprd ) 
  endif

   WRITE(oimfile(9 :11),'(I3.3)')gridi
   WRITE(oimfile(13:15),'(I3.3)')gridj
   WRITE(oimfile(17:19),'(I3.3)')gridk
 
   CALL write_grd(oimfile,ctl%nlon,ctl%nlat,ctl%nfields,obsimpact,totalundefmask,ctl%undefbin)

ENDDO ![End do over points]

ENDIF ![Endif over compuation of covariance]

IF( computeindex )THEN

 !Get the mean condensate.
 CALL get_var_startend( ctl , TRIM('qhydro') ,cond_si,cond_ei,cond_i)
 local_nzc=ctl%varlev( cond_i )
 ALLOCATE( mean_cond( ctl%nlon , ctl%nlat , local_nzc ) )
 IF( cond_i /= 0 )THEN
   CALL compute_moments(ensemble(:,:,cond_si:cond_ei,:),ctl%nlon,ctl%nlat,local_nzc,nbv,1,mean_cond, &
                           totalundefmask(:,:,cond_si:cond_ei),ctl%undefbin)
 ELSE
   mean_cond=0.0e0
   WRITE(*,*)"Warning: Could not find condesate variable. ALl points will be considered as no rain points"
 ENDIF

 DO ii=1,NBASE_VARS !Loop over base vars
  DO jj=1,NCOV_VARS !Loop over cov vars

    CALL get_var_startend( ctl , TRIM(BASE_VARS(ii)) , var1_si,var1_ei,var1_i)
    CALL get_var_startend( ctl , TRIM(COV_VARS(jj))  , var2_si,var2_ei,var2_i)

    local_nz1=ctl%varlev( var1_i )
    local_nz2=ctl%varlev( var2_i )

    WRITE(*,*)"Processing covariance strength between ",TRIM(BASE_VARS(ii))," and ",TRIM(COV_VARS(jj))

    IF( ( local_nz1 /= local_nz2 ) .OR. ( local_nz1 /= local_nzc ) )THEN
      WRITE(*,*)"Error: BASE variable and COV variable must have the same number of vertical levels"
      WRITE(*,*)TRIM(BASE_VARS(ii)),TRIM(COV_VARS(jj)),local_nz1,local_nz2,var1_i,var2_i
      STOP
    ENDIF
    
    ALLOCATE( covariance_index(ctl%nlon,ctl%nlat,local_nz1) , correlation_index(ctl%nlon,ctl%nlat,local_nz1) &
          , correlationdist_index(ctl%nlon,ctl%nlat,local_nz1) )

    ALLOCATE( local_undefmask( ctl%nlon , ctl%nlat , local_nz1)  )
 
    ALLOCATE( cov_profile(delta,local_nz1,2) , corr_profile(delta,local_nz1,2) , num_profile(delta,local_nz1,2) )

    local_undefmask=totalundefmask(:,:,var1_si:var1_ei)
    where( .not. totalundefmask(:,:,var2_si:var2_ei) )
       local_undefmask=.false.
    endwhere

    CALL covariance_strenght(ensemble(:,:,var1_si:var1_ei,:),ensemble(:,:,var2_si:var2_ei,:),mean_cond,covariance_index,  &
                             correlation_index,correlationdist_index,cov_profile,corr_profile,num_profile,ctl%nlon,       &
                             ctl%nlat,local_nz1,nbv,delta,skip,local_undefmask,ctl%undefbin)

    tmpfilename='covindex_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.grd'

    CALL write_grd(tmpfilename,ctl%nlon,ctl%nlat,local_nz1,covariance_index,totalundefmask,ctl%undefbin)

    tmpfilename='corrindex_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.grd'

    CALL write_grd(tmpfilename,ctl%nlon,ctl%nlat,local_nz1,correlation_index,totalundefmask,ctl%undefbin)

    tmpfilename='corrdistindex_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.grd'

    CALL write_grd(tmpfilename,ctl%nlon,ctl%nlat,local_nz1,correlationdist_index,totalundefmask,ctl%undefbin)

    tmpfilename='cov_profile_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.txt'

    CALL write_profile(tmpfilename,cov_profile,delta,local_nz1,2)

    tmpfilename='corr_profile_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.txt'

    CALL write_profile(tmpfilename,corr_profile,delta,local_nz1,2)

    tmpfilename='num_profile_' // TRIM(BASE_VARS(ii)) // '_' // TRIM(COV_VARS(jj)) // '.txt'

    CALL write_profile(tmpfilename,REAL(num_profile,r_sngl),delta,local_nz1,2)

    DEALLOCATE( cov_profile , corr_profile , num_profile )

    DEALLOCATE( local_undefmask )

    DEALLOCATE( covariance_index , correlation_index , correlationdist_index)

  ENDDO !End loop over cov vars
 ENDDO  !End loop over base vars
ENDIF ![Endif over computatio of the covariance and correlation index]

end program covariance_matrix
!===============================================================================
