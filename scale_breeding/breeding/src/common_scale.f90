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
!
!=======================================================================
!!$USE OMP_LIB
  use common
  use netcdf , only : NF90_WRITE , NF90_NOWRITE
  use common_ncio
  use common_namelist
  use scale_const


  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: nv3d=17    ! 
  INTEGER,PARAMETER :: nv2d=2     ! 
  INTEGER,PARAMETER :: iv3d_rho=1  !-- State in restart files
  INTEGER,PARAMETER :: iv3d_rhou=2 !
  INTEGER,PARAMETER :: iv3d_rhov=3 !
  INTEGER,PARAMETER :: iv3d_rhow=4 !
  INTEGER,PARAMETER :: iv3d_rhot=5 !
  INTEGER,PARAMETER :: iv3d_q=6
  INTEGER,PARAMETER :: iv3d_qc=7
  INTEGER,PARAMETER :: iv3d_qr=8
  INTEGER,PARAMETER :: iv3d_qi=9
  INTEGER,PARAMETER :: iv3d_qs=10
  INTEGER,PARAMETER :: iv3d_qg=11

! For norm computation
  INTEGER,PARAMETER :: iv3d_u=12
  INTEGER,PARAMETER :: iv3d_v=13
  INTEGER,PARAMETER :: iv3d_w=14
  INTEGER,PARAMETER :: iv3d_t=15
  INTEGER,PARAMETER :: iv3d_hgt=16
  INTEGER,PARAMETER :: iv3d_p=17

  INTEGER,PARAMETER :: iv2d_lon=1
  INTEGER,PARAMETER :: iv2d_lat=2

!  INTEGER,SAVE :: nlon  ! # grids in I-direction [subdomain]
!  INTEGER,SAVE :: nlat  ! # grids in J-direction [subdomain]
  INTEGER,SAVE :: nlev  ! # grids in K-direction
  INTEGER,SAVE :: nlong ! # grids in I-direction [global domain]
  INTEGER,SAVE :: nlatg ! # grids in J-direction [global domain]
!  INTEGER,SAVE :: nij0
!  INTEGER,SAVE :: nlevall
!  INTEGER,SAVE :: nlevalld
  INTEGER,SAVE :: ngpv
  INTEGER,SAVE :: ngpvd

!  INTEGER,SAVE :: PRC_NUM_X , PRC_NUM_Y !Number of subdomains that we have in X and Y.

  REAL(r_size) :: dx , dz ! Horizontal and vertical resolution

  INTEGER,PARAMETER :: vname_max = 10
  CHARACTER(vname_max),SAVE :: v3d_name(nv3d)
  CHARACTER(vname_max),SAVE :: v2d_name(nv2d)

  LOGICAL :: from_file_3d(nv3d) , from_file_2d(nv2d) 

  REAL(r_size) , allocatable :: CXG(:) , CYG(:) !Global domain x and y position
  REAL(r_size) , allocatable :: levels(:)

  INTEGER :: NSUB


CONTAINS

!-----------------------------------------------------------------------
SUBROUTINE set_common_scale(file_prefix)

  IMPLICIT NONE
  character(len=*), intent(in)    :: file_prefix
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
  call ncio_open(TRIM( filename ),NF90_NOWRITE,ncid)

  call ncio_get_dim(ncid, 'x' , nlons )
  call ncio_get_dim(ncid, 'y' , nlats )
  call ncio_get_dim(ncid, 'z' , nlev )

  call ncio_get_dim(ncid,'CXG', nlong )
  call ncio_get_dim(ncid,'CYG', nlatg )

  !PRC_NUM_X = NINT( REAL(nlong,r_size) / REAL(nlons,r_size) )
  !PRC_NUM_Y = NINT( REAL(nlatg,r_size) / REAL(nlats,r_size) )

  !  nsub= PRC_NUM_X * PRC_NUM_Y

    write(*,*)"--------------------------------------------------------------------"
    write(*,*)" NSUB =",nsub
    write(*,*)"--------------------------------------------------------------------"

    !nij0 = nlon * nlat
    !nlevall  = nlev * nv3d  + nv2d
    !ngpv  = nij0 * nlevall

  if( allocated(cxg) )deallocate(cxg)
  if( allocated(cyg) )deallocate(cyg)

  allocate(cxg(nlong),cyg(nlatg))
  call ncio_read(ncid, trim('CXG'),nlong, 1, cxg )
  call ncio_read(ncid, trim('CYG'),nlatg, 1, cyg )

  !Get the horizontal resolution computing the distance
  !between the two first grid points.
  dx=cxg(2)-cxg(1)

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


 
  call ncio_close(ncid) 

    !
    ! Variable names (same as in the NetCDF file)
    !
    ! state variables (in 'restart' files, for LETKF)
    v3d_name(iv3d_rho)  = 'DENS'
    v3d_name(iv3d_rhou) = 'MOMX'
    v3d_name(iv3d_rhov) = 'MOMY'
    v3d_name(iv3d_rhow) = 'MOMZ'
    v3d_name(iv3d_rhot) = 'RHOT'
    v3d_name(iv3d_q)    = 'QV'
    v3d_name(iv3d_qc)   = 'QC'
    v3d_name(iv3d_qr)   = 'QR'
    v3d_name(iv3d_qi)   = 'QI'
    v3d_name(iv3d_qs)   = 'QS'
    v3d_name(iv3d_qg)   = 'QG'
    !
    ! diagnostic variables (for norm computation)
    v3d_name(iv3d_u)    = 'U'
    v3d_name(iv3d_v)    = 'V'
    v3d_name(iv3d_w)    = 'W'
    v3d_name(iv3d_t)    = 'T'
    v3d_name(iv3d_hgt)  = 'height'
    v3d_name(iv3d_p)    = 'P'
    v2d_name(iv2d_lon)  = 'lon'
    v2d_name(iv2d_lat)  = 'lat'
    !
    ! state variables (in 'restart' files, for LETKF)
    from_file_3d(iv3d_rho)  =  .true.
    from_file_3d(iv3d_rhou) =  .true.
    from_file_3d(iv3d_rhov) =  .true.
    from_file_3d(iv3d_rhow) =  .true.
    from_file_3d(iv3d_rhot) =  .true.
    from_file_3d(iv3d_q)    =  .true.
    from_file_3d(iv3d_qc)   =  .true.
    from_file_3d(iv3d_qr)   =  .true.
    from_file_3d(iv3d_qi)   =  .true.
    from_file_3d(iv3d_qs)   =  .true.
    from_file_3d(iv3d_qg)   =  .true.
    !
    ! diagnostic variables (for norm computation)
    from_file_3d(iv3d_u)    = .false.
    from_file_3d(iv3d_v)    = .false.
    from_file_3d(iv3d_w)    = .false.
    from_file_3d(iv3d_hgt)  = .false.
    from_file_3d(iv3d_p)    = .false.

    from_file_2d(iv2d_lon)  = .true.
    from_file_2d(iv2d_lat)  = .true.


  RETURN
END SUBROUTINE set_common_scale



!!-----------------------------------------------------------------------
!! File I/O
!!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

SUBROUTINE read_restart_subdomain(filename,v3ds,v2ds,nlons,nlats)
  IMPLICIT NONE
  INTEGER     ,INTENT(IN) :: nlons , nlats
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(OUT) :: v3ds(nlev,nlons,nlats,nv3d)
  REAL(r_size),INTENT(OUT) :: v2ds(nlons,nlats,nv2d)
  integer :: iv3d,iv2d,ncid


  !Read data from a subdomain file.

  call ncio_open(trim(filename),NF90_NOWRITE,ncid)

  do iv3d = 1, nv3d
    if( from_file_3d( iv3d ) ) then
      write(6,'(1x,A,A15)') '*** Read 3D var: ', trim(v3d_name(iv3d))
      call ncio_read_3d_r8(ncid, trim(v3d_name(iv3d)), nlev, nlons, nlats, 1, v3ds(:,:,:,iv3d) )
    endif
  end do

  do iv2d = 1, nv2d
    if( from_file_2d( iv2d ) )then
      write(6,'(1x,A,A15)') '*** Read 2D var: ', trim(v2d_name(iv2d))
      call ncio_read_2d_r8(ncid, trim(v2d_name(iv2d)), nlons, nlats, 1, v2ds(:,:,iv2d))
    endif
  end do


  call ncio_close(ncid)

  

  RETURN
END SUBROUTINE read_restart_subdomain

SUBROUTINE write_restart_subdomain(filename,v3ds,v2ds,nlons,nlats)
  implicit none
  INTEGER , INTENT(IN)    :: nlons , nlats 
  CHARACTER(*),INTENT(IN) :: filename
  REAL(r_size),INTENT(IN) :: v3ds(nlev,nlons,nlats,nv3d)
  REAL(r_size),INTENT(IN) :: v2ds(nlons,nlats,nv2d)
  integer :: iv3d,iv2d,ncid

  !Write output to a subdomain file.

  call ncio_open(filename, NF90_WRITE, ncid)

  do iv3d = 1, nv3d
   if( from_file_3d( iv3d ) )then

    write(6,'(1x,A,A15)') '*** Write 3D var: ', trim(v3d_name(iv3d))
    call ncio_write(ncid, trim(v3d_name(iv3d)), nlev , nlons, nlats, 1, v3ds(:,:,:,iv3d))
   endif
  end do

  do iv2d = 1, nv2d
   if( from_file_2d( iv2d ) )then
    write(6,'(1x,A,A15)') '*** Write 2D var: ', trim(v2d_name(iv2d))
    call ncio_write(ncid, trim(v2d_name(iv2d)), nlons, nlats, 1, v2ds(:,:,iv2d))
   endif
  end do

  call ncio_close(ncid)

  RETURN
END SUBROUTINE write_restart_subdomain


SUBROUTINE read_restart(filenameprefix,v3dg,v2dg)

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filenameprefix
  REAL(r_size),INTENT(OUT) :: v3dg(nlev,nlong,nlatg,nv3d)
  REAL(r_size),INTENT(OUT) :: v2dg(nlong,nlatg,nv2d)
  REAL(r_size) , allocatable :: cxs(:) , cys(:)
  REAL(r_size) , allocatable :: v3ds(:,:,:,:)
  REAL(r_size) , allocatable :: v2ds(:,:,:)
  INTEGER                    :: nlons , nlats !Subdomain dimension
  REAL(r_size) :: dummy
 
  character(len=12) :: filesuffix = '.pe000000.nc'

  integer :: iv3d,iv2d,isub,istart,iend,jstart,jend,ii,jj,ncid

  dummy=0.0d0

  DO isub = 1 , nsub

     write (filesuffix(4:9),'(I6.6)') isub-1

     write(*,*)"Reading ",trim( filenameprefix ) // filesuffix

     call ncio_open( trim( filenameprefix ) // filesuffix ,NF90_NOWRITE,ncid)
     call ncio_get_dim(ncid, 'x' , nlons )
     call ncio_get_dim(ncid, 'y' , nlats )
 
     allocate( cxs(nlons) , cys(nlats) )

     call ncio_read_1d_r8(ncid, trim('x'),nlons, 1, cxs )
     call ncio_read_1d_r8(ncid, trim('y'),nlats, 1, cys )
     call ncio_close(ncid)
     allocate( v3ds(nlev,nlons,nlats,nv3d) , v2ds(nlons,nlats,nv2d) )
     call read_restart_subdomain( trim( filenameprefix ) // filesuffix , v3ds , v2ds , nlons , nlats )

     !Get the subdomain location with respect to the global domain.

     istart=0
     iend=0
     jstart=0
     jend=0

     DO ii = 1 , nlong
        if( abs( cxg( ii ) - cxs(1) )    < tiny(dummy) )istart = ii
        if( abs( cxg( ii ) - cxs(nlons) ) < tiny(dummy) )iend   = ii
     ENDDO
     DO jj = 1 , nlatg
        if( abs( cyg( jj ) - cys(1)  )   < tiny(dummy) )jstart = jj
        if( abs( cyg( jj ) - cys(nlats) ) < tiny(dummy) )jend   = jj
     ENDDO

     if( istart == 0 .or. jstart == 0 .or. iend == 0 .or. jend == 0 )then 
       WRITE(6,*)"Error: Could not find subdomain location in subroutine read_restart"
     endif

     v3dg(:,istart:iend,jstart:jend,:)=v3ds(:,:,:,:)
     v2dg(  istart:iend,jstart:jend,:)=v2ds(:,:,:)

     deallocate( v3ds , v2ds , cxs , cys)

  ENDDO

END SUBROUTINE read_restart


SUBROUTINE write_restart(filenameprefix,v3dg,v2dg)

  IMPLICIT NONE
  CHARACTER(*),INTENT(IN) :: filenameprefix
  REAL(r_size),INTENT(IN) :: v3dg(nlev,nlong,nlatg,nv3d)
  REAL(r_size),INTENT(IN) :: v2dg(nlong,nlatg,nv2d)
  REAL(r_size) , allocatable :: cxs(:) , cys(:)
  REAL(r_size) , allocatable :: v3ds(:,:,:,:)
  REAL(r_size) , allocatable :: v2ds(:,:,:)
  REAL(r_size) :: dummy
  INTEGER      :: nlons , nlats

  character(len=12) :: filesuffix = '.pe000000.nc'

  integer :: iv3d,iv2d,isub,istart,iend,jstart,jend,ii,jj,ncid

  dummy=0.0d0

  DO isub = 1 , nsub 

     write (filesuffix(4:9),'(I6.6)') isub - 1
     write(*,*)trim( filenameprefix ) // filesuffix

     call ncio_open( trim( filenameprefix ) // filesuffix ,NF90_NOWRITE,ncid)
     call ncio_get_dim(ncid, 'x' , nlons )
     call ncio_get_dim(ncid, 'y' , nlats )

     allocate( cxs(nlons) , cys(nlats) )

     call ncio_read(ncid, trim('x'),nlons, 1, cxs )
     call ncio_read(ncid, trim('y'),nlats, 1, cys )
     call ncio_close(ncid)

     !Get the subdomain location with respect to the global domain.

     istart=0
     iend=0
     jstart=0
     jend=0

     DO ii = 1 , nlong
        if( abs( cxg( ii ) - cxs(1) ) < tiny(dummy) )    istart = ii
        if( abs( cxg( ii ) - cxs(nlons) ) < tiny(dummy) ) iend   = ii
     ENDDO
     DO jj = 1 , nlatg
        if( abs( cyg( jj ) - cys(1) ) < tiny(dummy) )    jstart = jj
        if( abs( cyg( jj ) - cys(nlats) ) < tiny(dummy) ) jend   = jj
     ENDDO

     if( istart == 0 .or. iend == 0 .or. jstart == 0 .or. jend == 0 )then
       WRITE(6,*)"Error: Could not find subdomain location in subroutine read_restart"
     endif

     allocate( v3ds(nlev,nlons,nlats,nv3d) , v2ds(nlons,nlats,nv2d) )
     v3ds=0.0d0
     v2ds=0.0d0


     v3ds(:,:,:,:)=v3dg(:,istart:iend,jstart:jend,:)
     v2ds(:,:,:)  =v2dg(  istart:iend,jstart:jend,:)
     
     call write_restart_subdomain( trim( filenameprefix ) // filesuffix , v3ds , v2ds , nlons , nlats )

     deallocate ( v3ds , v2ds , cxs , cys )

  ENDDO

END SUBROUTINE write_restart


!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------


SUBROUTINE state_trans(v3dg)


  IMPLICIT NONE

  REAL(r_size),INTENT(INOUT) :: v3dg(nlev,nlong,nlatg,nv3d)
  REAL(r_size) :: rho,pres,temp
  real(r_size) :: qdry,CVtot,Rtot,CPovCV
  integer :: i,j,k,iv3d


!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CPovCV,rho,pres,temp) COLLAPSE(2)
  do j = 1, nlatg
    do i = 1, nlong
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, iv3d_qg ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
       enddo
         CVtot = v3dg(k,i,j,iv3d_q)*CONST_CVvap + v3dg(k,i,j,iv3d_qr)*CONST_CL + v3dg(k,i,j,iv3d_qc)*CONST_CL +  &
                 v3dg(k,i,j,iv3d_qi)*CONST_CI   + v3dg(k,i,j,iv3d_qg)*CONST_CI + v3dg(k,i,j,iv3d_qs)*CONST_CI +  & 
                 qdry * CONST_CVdry

         Rtot  = CONST_Rdry  * qdry + CONST_Rvap * v3dg(k,i,j,iv3d_q)
         CPovCV = ( CVtot + Rtot ) / CVtot

       rho = v3dg(k,i,j,iv3d_rho)
       pres = CONST_PRE00 * ( v3dg(k,i,j,iv3d_rhot) * Rtot / CONST_PRE00 )**CPovCV
       temp = pres / ( rho * Rtot )

       v3dg(k,i,j,iv3d_u) = v3dg(k,i,j,iv3d_rhou) / rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_v) = v3dg(k,i,j,iv3d_rhov) / rho !!!!!!
       v3dg(k,i,j,iv3d_w) = v3dg(k,i,j,iv3d_rhow) / rho !!!!!!
       v3dg(k,i,j,iv3d_t) = temp
       v3dg(k,i,j,iv3d_p) = pres

      enddo
    enddo
  enddo


  RETURN
END SUBROUTINE state_trans


SUBROUTINE state_trans_inv(v3dg)
  IMPLICIT NONE
  REAL(r_size),INTENT(INOUT) :: v3dg(nlev,nlong,nlatg,nv3d)
  REAL(r_size) :: rho,rhot
  real(r_size) :: qdry,CVtot,Rtot,CVovCP
  integer :: i,j,k,iv3d

  !Check that nor the vapor, nor the condensates are negative.
  do iv3d = 1, nv3d
    if ( iv3d == iv3d_q .or. iv3d == iv3d_qc .or. iv3d == iv3d_qr .or.  &
         iv3d == iv3d_qi .or. iv3d == iv3d_qs .or. iv3d == iv3d_qg ) then
      v3dg(:,:,:,iv3d) = max(v3dg(:,:,:,iv3d), 0.0d0)
    end if
  end do

!!$OMP PARALLEL DO PRIVATE(i,j,k,iv3d,qdry,CVtot,Rtot,CVovCP,rho,rhot) COLLAPSE(2)
  do j = 1, nlatg
    do i = 1, nlong
      do k = 1, nlev
       qdry  = 1.0d0
       CVtot = 0.0d0
       do iv3d = iv3d_q, iv3d_qg ! loop over all moisture variables
         qdry  = qdry - v3dg(k,i,j,iv3d)
       enddo
 
         CVtot = v3dg(k,i,j,iv3d_q)*CONST_CVvap + v3dg(k,i,j,iv3d_qr)*CONST_CL + v3dg(k,i,j,iv3d_qc)*CONST_CL +  &
               v3dg(k,i,j,iv3d_qi)*CONST_CI   + v3dg(k,i,j,iv3d_qg)*CONST_CI + v3dg(k,i,j,iv3d_qs)*CONST_CI   +  &
               qdry * CONST_CVdry

         Rtot  = CONST_Rdry  * qdry + CONST_Rvap * v3dg(k,i,j,iv3d_q)
         CVovCP = CVtot / ( CVtot + Rtot )

       rho = v3dg(k,i,j,iv3d_p) / (Rtot * v3dg(k,i,j,iv3d_t))
       rhot = CONST_PRE00 / Rtot * (v3dg(k,i,j,iv3d_p) / CONST_PRE00)**CVovCP

       v3dg(k,i,j,iv3d_rhot) = rhot
       v3dg(k,i,j,iv3d_rhow) = v3dg(k,i,j,iv3d_w) * rho !!!!!! inaccurate! do not consider staggered grid !!!!!!
       v3dg(k,i,j,iv3d_rhov) = v3dg(k,i,j,iv3d_v) * rho !!!!!!
       v3dg(k,i,j,iv3d_rhou) = v3dg(k,i,j,iv3d_u) * rho !!!!!!
       v3dg(k,i,j,iv3d_rho) = rho
      enddo
    enddo
  enddo
!!$OMP END PARALLEL DO

  RETURN
END SUBROUTINE state_trans_inv

END MODULE common_scale
