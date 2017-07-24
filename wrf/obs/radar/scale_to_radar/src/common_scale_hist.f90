MODULE common_scale
!=======================================================================
!
! [PURPOSE:] Common Information for SCALE
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   10/03/2012 Guo-Yuan Lien     modified for GFS model
!   07/24/2014 Guo-Yuan Lien     modified for SCALE model
!   13/02/2017 Juan Ruiz         adapted  for SCALE breeding 
!   11/07/2017 Juan Ruiz         adapted  for SCALE_TO_RADAR
!
!=======================================================================
!!$USE OMP_LIB
  use common
  use netcdf , only : NF90_WRITE , NF90_NOWRITE
  use common_ncio
  use common_namelist
  use scale_const
  use scale_mapproj

  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nv3d=12    ! 
  INTEGER,PARAMETER :: nv2d=2     ! 
  INTEGER,PARAMETER :: iv3d_q=1
  INTEGER,PARAMETER :: iv3d_qc=2
  INTEGER,PARAMETER :: iv3d_qr=3
  INTEGER,PARAMETER :: iv3d_qi=4
  INTEGER,PARAMETER :: iv3d_qs=5
  INTEGER,PARAMETER :: iv3d_qg=6
  INTEGER,PARAMETER :: iv3d_u=7
  INTEGER,PARAMETER :: iv3d_v=8
  INTEGER,PARAMETER :: iv3d_w=9
  INTEGER,PARAMETER :: iv3d_t=10
  INTEGER,PARAMETER :: iv3d_p=11
  INTEGER,PARAMETER :: iv3d_hgt=12

  INTEGER,PARAMETER :: iv2d_lon=1
  INTEGER,PARAMETER :: iv2d_lat=2

!  INTEGER,SAVE :: nlon  ! # grids in I-direction [subdomain]
!  INTEGER,SAVE :: nlat  ! # grids in J-direction [subdomain]
  INTEGER,SAVE :: nlev  ! # grids in K-direction with halo
  INTEGER,SAVE :: nlonh ! # grids in I-direction with halo [global domain]
  INTEGER,SAVE :: nlath ! # grids in J-direction with halo [global domain]
  INTEGER,SAVE :: ntime ! # of times in file.

  INTEGER,SAVE :: nlong2
  INTEGER,SAVE :: nlatg2
  INTEGER,SAVE :: nlev2 !Number of levs + 2*KHALO
!  INTEGER,SAVE :: nij0
!  INTEGER,SAVE :: nlevall
!  INTEGER,SAVE :: nlevalld
  INTEGER,SAVE :: ngpv
  INTEGER,SAVE :: ngpvd

  REAL(r_size) , allocatable  :: topo(:,:)  !#Model topography.

!  INTEGER,SAVE :: PRC_NUM_X , PRC_NUM_Y !Number of subdomains that we have in X and Y.

  REAL(r_size) :: dx , dy , dz ! Horizontal and vertical resolution

  INTEGER,PARAMETER :: vname_max = 10
  CHARACTER(vname_max),SAVE :: v3d_name(nv3d)
  CHARACTER(vname_max),SAVE :: v2d_name(nv2d)

  LOGICAL :: from_file_3d(nv3d) , from_file_2d(nv2d) 

  REAL(r_size) , allocatable :: CXG(:) , CYG(:) !Global domain x and y position
  REAL(r_size) , allocatable :: FXG(:) , FYG(:) !Global grid face possition.
  REAL(r_size) , allocatable :: CZ(:),CX(:),CY(:) !Location of vertical levels + halo.
  REAL(r_size) , allocatable :: levels(:)

  REAL(r_size) :: GRID_DOMAIN_CENTER_X , GRID_DOMAIN_CENTER_Y

  INTEGER :: KHALO , IHALO , JHALO
  INTEGER :: NSUB , TMPINT

CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE set_common_scale(file_prefix , topo_file_prefix )

  IMPLICIT NONE
  character(len=*), intent(in)    :: file_prefix , topo_file_prefix
  character(len=100)              :: filename
  integer                         :: ncid , ierr
  integer                         :: ii , jj , kk , isub
  integer                         :: nlons , nlats  
  character(len=12)               :: filesuffix='.pe000000.nc'

  WRITE(6,'(A)') 'Hello from set_common_scale'

  !Initialize model constants used in the variable transform operations.
  call const_init 

  isub=0
  DO 
     write (filesuffix(4:9),'(I6.6)') isub
     filename = trim( file_prefix ) // filesuffix
     open(unit=55,file=filename , status='old', iostat = ierr )
 
     if( ierr == 0 )then 
       write(*,*)"We found the file ",filename
       isub=isub+1
       close(55)
     else
       write(*,*)"That was the last file"
       exit
     endif
  ENDDO
  
  nsub=isub !Set the number of subdomains that we have.
 
  write (filesuffix(4:9),'(I6.6)') 0
  filename = trim( file_prefix ) // filesuffix

  !GET THE DIMENSIONS FROM THE NETCDF FILE AND NAMELIST
  write(*,*)filename
  call ncio_open(TRIM( filename ),NF90_NOWRITE,ncid)


  call ncio_get_dim(ncid, 'x' , nlons )
  call ncio_get_dim(ncid, 'y' , nlats )
  call ncio_get_dim(ncid, 'z' , nlev )
  call ncio_get_dim(ncid, 'time', ntime )

  call ncio_get_dim(ncid,'CXG', nlonh )
  call ncio_get_dim(ncid,'CYG', nlath )

  !PRC_NUM_X = NINT( REAL(nlonh,r_size) / REAL(nlons,r_size) )
  !PRC_NUM_Y = NINT( REAL(nlath,r_size) / REAL(nlats,r_size) )

  !  nsub= PRC_NUM_X * PRC_NUM_Y

    write(*,*)"--------------------------------------------------------------------"
    write(*,*)" NSUB =",nsub
    write(*,*)"--------------------------------------------------------------------"

    !nij0 = nlon * nlat
    !nlevall  = nlev * nv3d  + nv2d
    !ngpv  = nij0 * nlevall

  if( allocated(cxg) )deallocate(cxg)
  if( allocated(cyg) )deallocate(cyg)

  allocate(cxg(nlonh),cyg(nlath))
  call ncio_read(ncid, trim('CXG'),nlonh, 1, cxg )
  call ncio_read(ncid, trim('CYG'),nlath, 1, cyg )

  !Get the horizontal resolution computing the distance
  !between the two first grid points.
  dx=cxg(2)-cxg(1)
  dy=cyg(2)-cyg(1)

  !write(*,*)nlon,nlat,nlev

  if(allocated(levels) )deallocate(levels)

  allocate( levels(nlev) )
  call ncio_read_1d_r8(ncid, trim('z'),nlev,1,levels)

  !Compute average vertical resolution.
  dz = 0.0d0
  DO kk = 1 , nlev - 1
   dz = dz + levels(kk+1) - levels(kk)
  ENDDO
  dz = dz / REAL( nlev - 1 , r_size ) 

  call ncio_get_dim(ncid,'FXG', nlong2 )
  allocate( FXG(nlong2) )
  call ncio_read(ncid, trim('FXG'),nlong2, 1, FXG )

  call ncio_get_dim(ncid,'FYG', nlatg2 )
  if( allocated(FYG) )deallocate(FYG)
  allocate( FYG(nlatg2) )
  call ncio_read(ncid, trim('FYG'),nlatg2, 1, FYG )

  GRID_DOMAIN_CENTER_X = 0.5d0 * ( FXG(0) + FXG(nlong2) )
  GRID_DOMAIN_CENTER_Y = 0.5d0 * ( FYG(0) + FYG(nlatg2) )

  !Initialize projection routines.
  CALL MPRJ_setup( GRID_DOMAIN_CENTER_X, GRID_DOMAIN_CENTER_Y )

  call ncio_get_dim(ncid,'CZ',nlev2)  !Number of levels + 2*KHALO
  if( allocated(CZ) )deallocate(CZ)
  allocate(CZ(nlev2))
  call ncio_read(ncid,trim('CZ'),nlev2,1,CZ )
  KHALO=( nlev2-nlev ) / 2.0d0 

  call ncio_get_dim(ncid,'CX', tmpint )  !Number of longitudes + 2*IHALO
  IHALO=( tmpint-nlons ) / 2.0d0

  call ncio_get_dim(ncid,'CY', tmpint )  !Number of longitudes + 2*JHALO
  JHALO=( tmpint-nlats ) / 2.0d0

  call ncio_close(ncid) 

    !
    ! Variable names (same as in the NetCDF file)
    !

    !V3D
    v3d_name(iv3d_q)    = 'QV'
    v3d_name(iv3d_qc)   = 'QC'
    v3d_name(iv3d_qr)   = 'QR'
    v3d_name(iv3d_qi)   = 'QI'
    v3d_name(iv3d_qs)   = 'QS'
    v3d_name(iv3d_qg)   = 'QG'
    v3d_name(iv3d_u)    = 'U'
    v3d_name(iv3d_v)    = 'V'
    v3d_name(iv3d_w)    = 'W'
    v3d_name(iv3d_t)    = 'T'
    v3d_name(iv3d_hgt)  = 'height'
    v3d_name(iv3d_p)    = 'PRES'
    !V2D
    v2d_name(iv2d_lon)  = 'lon'
    v2d_name(iv2d_lat)  = 'lat'
    !
    ! state variables (in 'history' files, for LETKF)
    from_file_3d(iv3d_q)    =  .true.
    from_file_3d(iv3d_qc)   =  .true.
    from_file_3d(iv3d_qr)   =  .true.
    from_file_3d(iv3d_qi)   =  .true.
    from_file_3d(iv3d_qs)   =  .true.
    from_file_3d(iv3d_qg)   =  .true.
    from_file_3d(iv3d_u)    =  .true.
    from_file_3d(iv3d_v)    =  .true.
    from_file_3d(iv3d_w)    =  .true.
    from_file_3d(iv3d_t)    =  .true.
    from_file_3d(iv3d_hgt)  =  .true.
    from_file_3d(iv3d_p)    =  .true.

    from_file_2d(iv2d_lon)  = .true.
    from_file_2d(iv2d_lat)  = .true.


    !Read the topography.
 
    allocate( topo( nlonh , nlath ) )
    write(*,*)"Reading the topography"
    CALL read_topo(topo_file_prefix,topo)



  RETURN
END SUBROUTINE set_common_scale

!!-----------------------------------------------------------------------
!! File I/O
!!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE read_history_subdomain(filename,v3dg,v2dg,rtime)
  IMPLICIT NONE
  INTEGER     ,INTENT(IN) :: rtime !Read time.
  INTEGER                 :: nlons , nlats
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(INOUT) :: v3dg(nlev,nlonh,nlath,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2dg(nlonh,nlath,nv2d)

  REAL(r_size),ALLOCATABLE :: v3dstmp(:,:,:)

  REAL(r_size),ALLOCATABLE :: v3ds(:,:,:,:)
  REAL(r_size),ALLOCATABLE :: v2ds(:,:,:)
  REAL(r_size),ALLOCATABLE :: cxs(:),cys(:)
  integer :: iv3d,iv2d,ncid,istart,iend,jstart,jend
  integer :: ii,jj,i,j,k
  real(r_size) :: dummy = 0.0d0

  !Read data from a subdomain file.

  call ncio_open(trim(filename),NF90_NOWRITE,ncid)

     write(*,*)"Reading ",trim( filename)
     call ncio_get_dim(ncid, 'x' , nlons )
     call ncio_get_dim(ncid, 'y' , nlats )

     allocate( cxs(nlons) , cys(nlats) )

     call ncio_read_1d_r8(ncid, trim('x'),nlons, 1, cxs )
     call ncio_read_1d_r8(ncid, trim('y'),nlats, 1, cys )
     !call ncio_close(ncid)

     allocate( v3ds(nlev,nlons,nlats,nv3d) , v2ds(nlons,nlats,nv2d) )
     allocate( v3dstmp(nlons,nlats,nlev)  )


     do iv3d = 1, nv3d
      if( from_file_3d( iv3d ) ) then
        write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
        call ncio_read_3d_r8(ncid, trim(v3d_name(iv3d)), nlons, nlats, nlev, rtime, v3dstmp )
        forall (i=1:nlons, j=1:nlats, k=1:nlev) v3ds(k,i,j,iv3d) = v3dstmp(i,j,k) ! use FORALL to change order of dimensions
      endif
     end do

     do iv2d = 1, nv2d
      if( from_file_2d( iv2d ) )then
       write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
       call ncio_read_2d_r8(ncid, trim(v2d_name(iv2d)), nlons , nlats , rtime, v2ds(:,:,iv2d) )
      endif
     end do

     istart=0
     iend=0
     jstart=0
     jend=0

     DO ii = 1 , nlonh
        if( abs( cxg( ii ) - cxs(1) )    < tiny(dummy) )istart = ii
        if( abs( cxg( ii ) - cxs(nlons) ) < tiny(dummy) )iend   = ii
     ENDDO
     DO jj = 1 , nlath
        if( abs( cyg( jj ) - cys(1)  )   < tiny(dummy) )jstart = jj
        if( abs( cyg( jj ) - cys(nlats) ) < tiny(dummy) )jend   = jj
     ENDDO

     if( istart == 0 .or. jstart == 0 .or. iend == 0 .or. jend == 0 )then
       WRITE(6,*)"Error: Could not find subdomain location in subroutine read_restart"
     endif

     v3dg(:,istart:iend,jstart:jend,:)=v3ds(:,:,:,:)
     v2dg(  istart:iend,jstart:jend,:)=v2ds(:,:,:)

     deallocate( v3ds , v2ds , cxs , cys)



  call ncio_close(ncid)

  

  RETURN
END SUBROUTINE read_history_subdomain



SUBROUTINE read_history(filenameprefix,v3dg,v2dg,rtime)

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filenameprefix
  INTEGER     ,INTENT(IN) :: rtime          !Time to be read (history files may contain several times)
  REAL(r_size),INTENT(OUT) :: v3dg(nlev,nlonh,nlath,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dg(nlonh,nlath,nv2d)
  REAL(r_size) :: dummy
 
  character(len=12) :: filesuffix = '.pe000000.nc'

  integer :: iv3d,iv2d,isub,istart,iend,jstart,jend,ii,jj,ncid

  dummy=0.0d0

  DO isub = 1 , nsub

     write (filesuffix(4:9),'(I6.6)') isub-1

     write(*,*)"Reading ",trim( filenameprefix ) // filesuffix

     call read_history_subdomain( trim( filenameprefix ) // filesuffix , v3dg , v2dg , rtime )

  ENDDO

END SUBROUTINE read_history


SUBROUTINE read_topo(filenameprefix,topo)

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filenameprefix
  REAL(r_size),INTENT(OUT) :: topo(nlonh,nlath)
  REAL(r_size) , allocatable :: cxs(:) , cys(:)
  REAL(r_size) , allocatable :: topos(:,:)
  INTEGER                    :: nlons , nlats !Subdomain dimension
  REAL(r_size) :: dummy
  character(len=100) :: filename
 
  character(len=12) :: filesuffix = '.pe000000.nc'

  integer :: iv3d,iv2d,isub,istart,iend,jstart,jend,ii,jj,ncid

  dummy=0.0d0

  DO isub = 1 , nsub

     write (filesuffix(4:9),'(I6.6)') isub-1

     filename=trim( filenameprefix ) // filesuffix

     write(*,*)"Reading ",filename

     call ncio_open(filename,NF90_NOWRITE,ncid)

     call ncio_get_dim(ncid, 'x' , nlons )
     call ncio_get_dim(ncid, 'y' , nlats )
 
     allocate( cxs(nlons) , cys(nlats) )

     call ncio_read_1d_r8(ncid, trim('x'),nlons, 1, cxs )
     call ncio_read_1d_r8(ncid, trim('y'),nlats, 1, cys )
     !call ncio_close(ncid)
     allocate( topos(nlons,nlats) )
      
     !call ncio_open(trim(filename),NF90_NOWRITE,ncid)
     call ncio_read_2d_r8(ncid,'TOPO', nlons, nlats, 1, topos )
     call ncio_close(ncid)

     !Get the subdomain location with respect to the global domain.

     istart=0
     iend=0
     jstart=0
     jend=0

     DO ii = 1 , nlonh
        if( abs( cxg( ii ) - cxs(1) )    < tiny(dummy) )istart = ii
        if( abs( cxg( ii ) - cxs(nlons) ) < tiny(dummy) )iend   = ii
     ENDDO
     DO jj = 1 , nlath
        if( abs( cyg( jj ) - cys(1)  )   < tiny(dummy) )jstart = jj
        if( abs( cyg( jj ) - cys(nlats) ) < tiny(dummy) )jend   = jj
     ENDDO

     if( istart == 0 .or. jstart == 0 .or. iend == 0 .or. jend == 0 )then 
       WRITE(6,*)"Error: Could not find subdomain location in subroutine read_restart"
     endif

     topo(istart:iend,jstart:jend)=topos
     
     deallocate( topos , cxs , cys)

  ENDDO

END SUBROUTINE read_topo


END MODULE common_scale
