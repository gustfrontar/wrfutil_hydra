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

  call parse_ctl(ctl_file , ctl ) !Get ctl info (varnames, levs, grid size, etc)

  ALLOCATE( ensemble(ctl%nlon,ctl%nlat,ctl%nfields,nbv) )
  ALLOCATE( pensemble(nbv) )  
  ALLOCATE( covar(ctl%nlon,ctl%nlat,ctl%nfields) )
  ALLOCATE( obsimpact(ctl%nlon,ctl%nlat,ctl%nfields) )

  ALLOCATE( undefmask(ctl%nlon,ctl%nlat,ctl%nfields) , totalundefmask(ctl%nlon,ctl%nlat,ctl%nfields) )

  if( bootstrap )then
    ALLOCATE( bssprd(ctl%nlon,ctl%nlat,ctl%nfields) , bsmean(ctl%nlon,ctl%nlat,ctl%nfields) )
  endif

  !First read all the ensemble and store it in memory (if this is not possible
  !then the code should be modified).

   totalundefmask=.true.

DO i=1,nbv

     WRITE( fcstfile(5:9), '(I5.5)' )i
     call read_grd(fcstfile,ctl%nlon,ctl%nlat,ctl%nfields,ensemble(:,:,:,i),undefmask,ctl%undefbin)

     WHERE( .NOT. undefmask )
        totalundefmask=.false.
     END WHERE

ENDDO ![End do over ensemble members]


   !Loop over points.

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
     DO ii = 1,ctl%nlon
      ALLOCATE( tmpsample(nbv,2) , tmp(nbv,2) )
      tmp(:,1)=pensemble
      DO jj = 1,ctl%nlat
       DO kk = 1,ctl%nfields
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
      DEALLOCATE(tmpsample,tmp) 
     ENDDO
!$OMP END PARALLEL DO
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

  endif

   WRITE(oimfile(9 :11),'(I3.3)')gridi
   WRITE(oimfile(13:15),'(I3.3)')gridj
   WRITE(oimfile(17:19),'(I3.3)')gridk
 
   CALL write_grd(oimfile,ctl%nlon,ctl%nlat,ctl%nfields,obsimpact,totalundefmask,ctl%undefbin)

ENDDO ![End do over points]



end program covariance_matrix
!===============================================================================
