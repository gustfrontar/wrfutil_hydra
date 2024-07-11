MODULE common_mpi_wrf
!=======================================================================
!
! [PURPOSE:] MPI procedures
!
! [ATTENTION:]
!   DO NOT COMPILE WITH BOTH INLINE EXPANSION AND OMP OPTIONS TOGETHER
!    (Use ONE if you want, but DON'T USE BOTH AT THE SAME TIME)
!
! [HISTORY:]
!   01/23/2009 Takemasa Miyoshi  created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  USE common_mpi
  USE common_wrf
  USE common_smooth2d
  USE common_namelist
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: mpibufsize=10000 !200
  INTEGER,SAVE :: nij1
  INTEGER,SAVE :: nij1max
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lon1(:),lat1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonu1(:),latu1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonv1(:),latv1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: ri1(:),rj1(:)
  REAL(r_size),ALLOCATABLE,SAVE :: landmask1(:)

CONTAINS
SUBROUTINE set_common_mpi_wrf
  REAL(r_size),ALLOCATABLE :: ri(:,:),rj(:,:)
  INTEGER :: i,j,n,ierr

  WRITE(6,'(A)') 'Hello from set_common_mpi_wrf'

  !BROADCAST P00 AND DNW
  CALL MPI_BCAST(p0,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(dnw,nlev,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  CALL MPI_BCAST(znu,nlev,MPI_REAL,0,MPI_COMM_WORLD,ierr)


  !COMPUTE SIZE OF LOCAL DOMAINS
  i = MOD(nlon*nlat,nprocs)
  nij1max = (nlon*nlat - i)/nprocs + 1
  IF(myrank < i) THEN
    nij1 = nij1max
  ELSE
    nij1 = nij1max - 1
  END IF
  WRITE(6,'(A,I3.3,A,I6)') 'MYRANK ',myrank,' number of grid points: nij1= ',nij1
  ALLOCATE(nij1node(nprocs))
  DO n=1,nprocs
    IF(n-1 < i) THEN
      nij1node(n) = nij1max
    ELSE
      nij1node(n) = nij1max - 1
    END IF
  END DO

  ALLOCATE(phi1(nij1))
  ALLOCATE(landmask1(nij1))
  ALLOCATE(lon1(nij1))
  ALLOCATE(lat1(nij1))
  ALLOCATE(lonu1(nij1))
  ALLOCATE(latu1(nij1))
  ALLOCATE(lonv1(nij1))
  ALLOCATE(latv1(nij1))
  ALLOCATE(ri1(nij1))
  ALLOCATE(rj1(nij1))

  IF(myrank .eq. 0)THEN
  ALLOCATE(ri(nlon,nlat))
  ALLOCATE(rj(nlon,nlat))
  DO j=1,nlat
    DO i=1,nlon
      ri(i,j) = REAL(i,r_sngl)
      rj(i,j) = REAL(j,r_sngl)
    END DO
  END DO
  ENDIF
  CALL scatter_ngrd_mpi(0,SNGL(lon),lon1,1)
  CALL scatter_ngrd_mpi(0,SNGL(lat),lat1,1)
  CALL scatter_ngrd_mpi(0,SNGL(lonu),lonu1,1)
  CALL scatter_ngrd_mpi(0,SNGL(latu),latu1,1)
  CALL scatter_ngrd_mpi(0,SNGL(lonv),lonv1,1)
  CALL scatter_ngrd_mpi(0,SNGL(latv),latv1,1)
  CALL scatter_ngrd_mpi(0,SNGL(ri),ri1,1)
  CALL scatter_ngrd_mpi(0,SNGL(rj),rj1,1)
  CALL scatter_ngrd_mpi(0,SNGL(phi0),phi1,1)
  CALL scatter_ngrd_mpi(0,SNGL(landmask),landmask1,1)

  IF(myrank.eq.0)DEALLOCATE(ri,rj)  

  RETURN
END SUBROUTINE set_common_mpi_wrf
!-----------------------------------------------------------------------
! Scatter gridded data to processes (nrank -> all)
!-----------------------------------------------------------------------
!SUBROUTINE scatter_grd_mpi(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
!
!  !IF(mpibufsize > nij1max) THEN
!  !  WRITE(6,*)"SCATTER GRID MPI FAST",mpibufsize,nij1max
!    CALL scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!  !ELSE
!  !  WRITE(6,*)"SCATTER GRID MPI SAFE",mpibufsize,nij1max
!  !  CALL scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!  !END IF
!
!  RETURN
!END SUBROUTINE scatter_grd_mpi
!
!SUBROUTINE scatter_grd_mpi_safe(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
!  REAL(r_sngl) :: tmp(nij1max,nprocs)
!  REAL(r_sngl) :: bufs(mpibufsize,nprocs)
!  REAL(r_sngl) :: bufr(mpibufsize)
!  INTEGER :: i,j,k,n,ierr,ns,nr
!  INTEGER :: iter,niter
!
!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))
!
!  DO n=1,nv3d
!    DO k=1,nlev
!      IF(myrank == nrank) CALL grd_to_buf(v3dg(:,:,k,n),tmp)
!      DO iter=1,niter
!        IF(myrank == nrank) THEN
!          i = mpibufsize * (iter-1)
!          DO j=1,mpibufsize
!            i=i+1
!            IF(i > nij1max) EXIT
!            bufs(j,:) = tmp(i,:)
!          END DO
!        END IF
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                       & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1) EXIT
!          v3d(i,k,n) = REAL(bufr(j),r_size)
!        END DO
!      END DO
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    IF(myrank == nrank) CALL grd_to_buf(v2dg(:,:,n),tmp)
!    DO iter=1,niter
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          bufs(j,:) = tmp(i,:)
!        END DO
!      END IF
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        v2d(i,n) = REAL(bufr(j),r_size)
!      END DO
!    END DO
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  RETURN
!END SUBROUTINE scatter_grd_mpi_safe

!SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
!  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
!  REAL(r_sngl) :: bufr(nij1max,nlevall)
!  INTEGER :: j,k,n,ierr,ns,nr
!
!  ns = nij1max * nlevall
!  nr = ns
!  IF(myrank == nrank) THEN
!    j=0
!    DO n=1,nv3d
!      DO k=1,nlev
!        j = j+1
!        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
!      END DO
!    END DO
!
!    DO n=1,nv2d
!      j = j+1
!      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
!    END DO
!  END IF
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!
!  j=0
!  DO n=1,nv3d
!    DO k=1,nlev
!      j = j+1
!      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    j = j+1
!    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  RETURN
!END SUBROUTINE scatter_grd_mpi_fast


!-----------------------------------------------------------------------
! Gather gridded data (all -> nrank)
!-----------------------------------------------------------------------
!SUBROUTINE gather_grd_mpi(nrank,v3d,v2d,v3dg,v2dg)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
!  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!
!  !IF(mpibufsize > nij1max) THEN
!
!    CALL gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
!  !ELSE
!  !  WRITE(*,*)"GATHER GRID MPI SAFE",mpibufsize,nij1max
!  !  CALL gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
!  !END IF
!
!  RETURN
!END SUBROUTINE gather_grd_mpi
!
!SUBROUTINE gather_grd_mpi_safe(nrank,v3d,v2d,v3dg,v2dg)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
!  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_sngl) :: tmp(nij1max,nprocs)
!  REAL(r_sngl) :: bufs(mpibufsize)
!  REAL(r_sngl) :: bufr(mpibufsize,nprocs)
!  INTEGER :: i,j,k,n,ierr,ns,nr
!  INTEGER :: iter,niter
!
!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))
!
!  DO n=1,nv3d
!    DO k=1,nlev
!      DO iter=1,niter
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1) EXIT
!          bufs(j) = REAL(v3d(i,k,n),r_sngl)
!        END DO
!        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!        CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                      & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!        IF(myrank == nrank) THEN
!          i = mpibufsize * (iter-1)
!          DO j=1,mpibufsize
!            i=i+1
!            IF(i > nij1max) EXIT
!            tmp(i,:) = bufr(j,:)
!          END DO
!        END IF
!      END DO
!      IF(myrank == nrank) CALL buf_to_grd(tmp,v3dg(:,:,k,n))
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    DO iter=1,niter
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        bufs(j) = REAL(v2d(i,n),r_sngl)
!      END DO
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          tmp(i,:) = bufr(j,:)
!        END DO
!      END IF
!    END DO
!    IF(myrank == nrank) CALL buf_to_grd(tmp,v2dg(:,:,n))
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  RETURN
!END SUBROUTINE gather_grd_mpi_safe

!SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
!  INTEGER,INTENT(IN) :: nrank
!  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
!  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
!  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
!  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
!  REAL(r_sngl) :: bufs(nij1max,nlevall)
!  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
!  INTEGER :: j,k,n,ierr,ns,nr
!
!  ns = nij1max * nlevall
!  nr = ns
!  j=0
!  DO n=1,nv3d
!    DO k=1,nlev
!      j = j+1
!      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
!    END DO
!  END DO
!
!  DO n=1,nv2d
!    j = j+1
!    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!
!  IF(myrank == nrank) THEN
!    j=0
!    DO n=1,nv3d
!      DO k=1,nlev
!        j = j+1
!        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
!      END DO
!    END DO
!
!    DO n=1,nv2d
!      j = j+1
!      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
!    END DO
!  END IF
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!  RETURN
!END SUBROUTINE gather_grd_mpi_fast

!-----------------------------------------------------------------------
!Scatter_ngrd_mpi to scatter n fields with the safe scatter
!-----------------------------------------------------------------------
!This routine can be used to improve memory performance.
SUBROUTINE scatter_ngrd_mpi(nrank,vg,v,ndim)
  INTEGER,INTENT(IN) :: nrank,ndim
  REAL(r_sngl),INTENT(IN) :: vg(nlon,nlat,ndim)
  REAL(r_size),INTENT(OUT) :: v(nij1,ndim)
  REAL(r_sngl) :: tmp(nij1max,nprocs)
  REAL(r_sngl) :: bufs(nij1max,ndim,nprocs)
  REAL(r_sngl) :: bufr(nij1max,ndim)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  v=0.0d0
  bufs=0.0e0
  bufr=0.0e0


!SAFE
!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))
!
!  DO n=1,ndim
!    IF(myrank == nrank) CALL grd_to_buf(vg(:,:,n),tmp)
!    DO iter=1,niter
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          bufs(j,:) = tmp(i,:)
!        END DO
!      END IF
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
!                     & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        v(i,n) = REAL(bufr(j),r_size)
!      END DO
!    END DO
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!FAST
  ns = nij1max * ndim
  nr = ns
  IF(myrank == nrank) THEN
    DO n=1,ndim
        CALL grd_to_buf(vg(:,:,n),bufs(:,n,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  DO n=1,ndim
      v(:,n) = REAL(bufr(1:nij1,n),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_ngrd_mpi

!==========================================
!Gather n fields using mpi safe
!==========================================

SUBROUTINE gather_ngrd_mpi(nrank,v,vg,ndim)
  INTEGER,INTENT(IN) :: nrank,ndim
  REAL(r_size),INTENT(IN) :: v(nij1,ndim)
  REAL(r_sngl),INTENT(OUT) :: vg(nlon,nlat,ndim)
  REAL(r_sngl) :: tmp(nij1max,nprocs)
  REAL(r_sngl) :: bufs(nij1max,ndim)
  REAL(r_sngl) :: bufr(nij1max,ndim,nprocs)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

!SAFE
!  ns = mpibufsize
!  nr = ns
!  niter = CEILING(REAL(nij1max)/REAL(mpibufsize))
!
!  DO n=1,ndim
!    DO iter=1,niter
!      i = mpibufsize * (iter-1)
!      DO j=1,mpibufsize
!        i=i+1
!        IF(i > nij1) EXIT
!        bufs(j) = REAL(v(i,n),r_sngl)
!      END DO
!      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!      CALL MPI_GATHER(bufs,ns,MPI_REAL,&
!                    & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
!      IF(myrank == nrank) THEN
!        i = mpibufsize * (iter-1)
!        DO j=1,mpibufsize
!          i=i+1
!          IF(i > nij1max) EXIT
!          tmp(i,:) = bufr(j,:)
!        END DO
!      END IF
!    END DO
!    IF(myrank == nrank) CALL buf_to_grd(tmp,vg(:,:,n))
!  END DO
!
!  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)


!FAST
  ns = nij1max * ndim
  nr = ns
  DO n=1,ndim
      bufs(1:nij1,n) = REAL(v(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    DO n=1,ndim
        CALL buf_to_grd(bufr(:,n,:),vg(:,:,n))
    END DO

  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_ngrd_mpi

!-----------------------------------------------------------------------
! Read ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ens_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file !prefix
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: fieldg(nlon,nlat,nlev)
  INTEGER :: l,n,ll,im,im2,i,iunit,iv,ilev
  CHARACTER(20) :: filename='xxxx00000'
  INTEGER(4) :: ncid(member)



  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:9),'(A4,I5.5)') file,im
      !OPEN NC FILE
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
      CALL open_wrf_file(filename,'ro',ncid(im))
    END IF
  END DO

  !Read one field at a time and scatter to process.
  DO iv=1,nv3d
   !READ THE DATA 1 LEVEL AT A TIME
   !DO ilev=1,nlev
      DO l=1,ll 
      im = myrank+1 + (l-1)*nprocs
      IF( im <= member)THEN
          CALL read_var_wrf(ncid(im),iv,1,nlev,fieldg,'3d')
      ENDIF

     DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg,v3d(:,ilev,im2,iv),nlev)
      END IF
     END DO
    END DO
   !END DO  !End do over levels
  END DO  !End do over variables


  !Read one field at a time and scatter to process.
  DO iv=1,nv2d
    DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF( im <= member)THEN
       CALL read_var_wrf(ncid(im),iv,1,1,fieldg,'2d')
    ENDIF

     DO n=0,nprocs-1
       im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg,v2d(:,im2,iv),1)
       END IF
      END DO
    END DO
  END DO  !End do over variables

  DO l=1,ll
   im = myrank+1 + (l-1)*nprocs
   IF(im <= member) THEN
      CALL close_wrf_file(ncid(im))
    END IF
  END DO


!
  RETURN
END SUBROUTINE read_ens_mpi

!-----------------------------------------------------------------------
! Read bias data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_bias_mpi(filename,v3d,v2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl) :: fieldg(nlon,nlat)
  INTEGER :: l,n,ll,im,im2,i,iunit,iv,ilev
  CHARACTER(20),INTENT(IN) :: filename
  INTEGER(4) :: reclength,irec
  LOGICAL    :: ex
  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nlon * nlat
  IF(myrank .EQ. 0) THEN

 INQUIRE(FILE=filename,EXIST=ex)
 IF( ex )THEN
     OPEN(unit=33,FILE=filename,form='unformatted',access='direct',recl= nlon*nlat*reclength )

      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
      OPEN(unit=33,FILE=filename,form='unformatted',access='direct',recl=nlon*nlat*reclength)
  END IF
  !Read one field at a time and scatter to process.
  DO iv=1,nv3d
   !READ THE DATA 1 LEVEL AT A TIME
   DO ilev=1,nlev
          IF(myrank .EQ. 0 )THEN
           irec=(iv-1)*nlev + ilev
           READ(33,rec=irec)fieldg 
          ENDIF
          CALL scatter_ngrd_mpi(0,fieldg,v3d(:,ilev,iv),1)
   END DO  !End do over levels
  END DO  !End do over variables
  !Read one field at a time and scatter to process.
  DO iv=1,nv2d
       IF( myrank .EQ. 0 )THEN
           irec=nv3d*nlev + iv 
           READ(33,rec=irec)fieldg
       ENDIF
       CALL scatter_ngrd_mpi(0,fieldg,v2d(:,iv),1)
  END DO  !End do over variables
   IF( myrank .EQ. 0 ) THEN
      CLOSE(33)
   END IF

 ELSE !IF we could not found the file... then
   v3d=0.0d0
   v2d=0.0d0
   IF(myrank .EQ. 0)THEN
   WRITE(6,*)"WARNING: FILE ",filename," NOT FOUND."
   WRITE(6,*)"BIAS IS SET TO 0"
   ENDIF

 END IF


  RETURN
END SUBROUTINE read_bias_mpi

!-----------------------------------------------------------------------
! Gather bias and write it to a file.
!-----------------------------------------------------------------------
SUBROUTINE write_bias_mpi(filename,v3d,v2d)
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl) :: fieldg(nlon,nlat)
  INTEGER :: l,n,ll,im,im2,i,iunit,iv,ilev
  CHARACTER(20),INTENT(IN) :: filename
  INTEGER(4) :: reclength,irec
  INQUIRE(IOLENGTH=reclength) reclength
  reclength = reclength * nlon * nlat
  IF(myrank .EQ. 0) THEN
      OPEN(unit=33,FILE=filename,form='unformatted',access='direct',recl=nlon*nlat*reclength)
  END IF
  !Read one field at a time and scatter to process.
  DO iv=1,nv3d
   !READ THE DATA 1 LEVEL AT A TIME
   DO ilev=1,nlev
          CALL gather_ngrd_mpi(0,v3d(:,ilev,iv),fieldg,1)
          IF(myrank .EQ. 0 )THEN
           irec=(iv-1)*nlev+ilev
           WRITE(33,rec=irec)fieldg
          ENDIF
   END DO  !End do over levels
  END DO  !End do over variables
  !Read one field at a time and scatter to process.
  DO iv=1,nv2d
       CALL gather_ngrd_mpi(0,v2d(:,iv),fieldg,1)
       IF(myrank .EQ. 0)THEN
         irec=nv3d*nlev+iv
         WRITE(33,rec=irec)fieldg
       ENDIF
  END DO  !End do over variables
   IF( myrank .EQ. 0 ) THEN
      CLOSE(33)
   END IF
  RETURN
END SUBROUTINE write_bias_mpi

!-----------------------------------------------------------------------
! Read parameter ensemble data and distribute to processes
!-----------------------------------------------------------------------
SUBROUTINE read_ensp_mpi(file,member,vp2d)

CHARACTER(4),INTENT(IN) :: file !prefix
INTEGER,INTENT(IN) :: member
REAL(r_size),INTENT(OUT) :: vp2d(nij1,member,np2d)
REAL(r_sngl) :: fieldg(nlon,nlat)
INTEGER :: l,n,ll,im,im2,i,iunit,iv,ilev
CHARACTER(20) :: filename='xxxx00000'
INTEGER(4) :: ncid(member)

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:9),'(A4,I5.5)') file,im
      !OPEN NC FILE
      CALL open_wrf_file(filename,'ro',ncid(im))
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
    END IF
  END DO

 DO iv=1,np2d
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      CALL read_var_wrf(ncid(im),iv,1,1,fieldg,'pa')
    END IF
    DO n=0,nprocs-1
       im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg,vp2d(:,im2,iv),1)
      END IF
    END DO
  END DO

 END DO !END DO OVER THE PARAMETERS

  !CLOSE FILES
  DO l=1,ll
   im = myrank+1 + (l-1)*nprocs
   IF(im <= member) THEN
      CALL close_wrf_file(ncid(im))
    END IF
  END DO


RETURN
END SUBROUTINE read_ensp_mpi

!-----------------------------------------------------------------------
! Write parameter ensemble data and distribute to processes
!-----------------------------------------------------------------------

SUBROUTINE write_ensp_mpi(file,member,vp2d)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: vp2d(nij1,member,nv2d)
  REAL(r_sngl) :: fieldg(nlon,nlat)
  INTEGER :: l,n,ll,im,im2,iunit,iv,ilev,ncid(member),ierr,j,k
  CHARACTER(20) :: filename='xxxx00000'

  !OPEN FILES AND READ CONSTANTS.
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      IF(member > 1)THEN
      WRITE(filename(1:9),'(A4,I5.5)') file,im
      ELSE
      WRITE(filename(1:9),'(A4,I5.5)') file,nbv+1
      ENDIF
      !OPEN NC FILE
      CALL open_wrf_file(filename,'rw',ncid(im))
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    END IF
  END DO

  !WRITE 2D VARIABLES
  !1 GRID AT A TIME
  DO iv=1,np2d

   DO l=1,ll
    DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL gather_ngrd_mpi(n,vp2d(:,im2,iv),fieldg,1)
      END IF
    END DO

     im = myrank+1 + (l-1)*nprocs
     IF (im <= member)THEN
        CALL write_var_wrf(ncid(im),iv,1,1,fieldg,'pa')
     ENDIF
    END DO
  ENDDO

  !CLOSE FILES
  DO l=1,ll
   im = myrank+1 + (l-1)*nprocs
   IF(im <= member) THEN
      CALL close_wrf_file(ncid(im))
    END IF
  END DO

END SUBROUTINE write_ensp_mpi
!-----------------------------------------------------------------------
! Write ensemble data after collecting data from processes
!-----------------------------------------------------------------------
SUBROUTINE write_ens_mpi(file,member,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(4) :: sdmd, s1md
  REAL(r_sngl) :: fieldg(nlon,nlat,nlev)
  INTEGER :: l,n,ll,im,im2,iunit,iv,ilev,ncid(member),ierr,j,k
  CHARACTER(20) :: filename='xxxx00000'
  REAL(r_size),ALLOCATABLE ::  mu(:,:),ps(:,:),qv(:,:,:),t(:,:,:)
 
  ALLOCATE(mu(nij1,member),ps(nij1,member),qv(nij1,nlev,member))

  !OPEN FILES AND READ CONSTANTS.
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      IF(member > 1)THEN
        WRITE(filename(1:9),'(A4,I5.5)')file,im
      ELSE
        WRITE(filename(1:9),'(A4,I5.5)')file,nbv+1
      ENDIF
      !OPEN NC FILE
      CALL open_wrf_file(filename,'rw',ncid(im))
      WRITE(6,'(A,I3.3,2A)')'MYRANK ',myrank,' is writing a file ',filename
    END IF
  END DO

! FIRST READ SOME VARIABLES FROM THE GUES FILE AND COMPUTE INCREMENT IN MU
   !READ GUES QV
   !DO ilev=1,nlev
      DO l=1,ll
      im = myrank+1 + (l-1)*nprocs
      IF( im <= member)THEN
        CALL read_var_wrf(ncid(im),iv3d_qv,1,nlev,fieldg,'3d')
      ENDIF
     DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg,qv(:,:,im2),nlev)
      END IF
     END DO
    END DO
   !END DO  !End do over levels
   !READ MU
   DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF( im <= member)THEN
      CALL read_var_wrf(ncid(im),iv2d_mu,1,1,fieldg(:,:,1),'2d')
    ENDIF
    DO n=0,nprocs-1
       im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg(:,:,1),mu(:,im2),1)
       END IF
    END DO
   END DO
   !READ PSFC
   DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF( im <= member)THEN
      CALL read_var_wrf(ncid(im),iv2d_ps,1,1,fieldg(:,:,1),'2d')
    ENDIF
    DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg(:,:,1),ps(:,im2),1)
       END IF
    END DO
   END DO

   !COMPUTE UPDATED MU IN LOCAL VARIABLES

   !WRITE(*,*)mu(10,1),v3d(10,10,1,iv3d_qv),qv(10,10,1),dnw(k)
   DO j=1,nij1
    DO im=1,member
      sdmd=0.0
      s1md=0.0
      DO k=1,nlev-1
        IF(v3d(j,k,im,iv3d_qv) >= 0.0) THEN
          sdmd = sdmd + (v3d(j,k,im,iv3d_qv) - qv(j,k,im)) * REAL(dnw(k),r_size)     
        ELSE
          sdmd = sdmd + (0.0       - qv(j,k,im)) * REAL(dnw(k),r_size)
        END IF
        s1md = s1md + (1.0 + qv(j,k,im)) * REAL(dnw(k),r_size)
      END DO
      mu(j,im) = mu(j,im) - ((v2d(j,im,iv2d_ps) - ps(j,im)) + mu(j,im) * sdmd) / s1md
    END DO
   END DO
   !WRITE(*,*)mu(10,1),v3d(10,10,1,iv3d_qv),qv(10,10,1),dnw(k)

    !WRITE MU
    DO l=1,ll
     DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL gather_ngrd_mpi(n,mu(:,im2),fieldg(:,:,1),1)
      END IF
     END DO

     im = myrank+1 + (l-1)*nprocs
     IF (im <= member)THEN
       CALL write_var_wrf(ncid(im),iv2d_mu,1,1,fieldg(:,:,1),'2d')
     ENDIF
    END DO
    DEALLOCATE(ps,qv,mu)

! SECOND WRITE THE ENSEMBLE ONE FIELD AT A TIME

  !WRITE 3D VARIABLES
  !GATHER AND WRITE ONE GRID AT A TIME
  DO iv=1,nv3d
   !From T to theta'
   IF( iv == iv3d_t)THEN
    ALLOCATE(t(nij1,nlev,member))
    t=0.0d0
    DO k=1,nlev-1
      t(:,k,:) = v3d(:,k,:,iv3d_t) *(p0 /  v3d(:,k,:,iv3d_p)) ** (rd / cp) - t0
    ENDDO
   ENDIF

   !DO ilev=1,nlev
    DO l=1,ll
     DO n=0,nprocs-1
       im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        IF ( iv == iv3d_t)THEN
           CALL gather_ngrd_mpi(n,t(:,:,im2),fieldg,nlev)
        ELSE
           CALL gather_ngrd_mpi(n,v3d(:,:,im2,iv),fieldg,nlev)
        ENDIF
       END IF
     END DO
     im = myrank+1 + (l-1)*nprocs
     IF( im <= member)THEN
      IF(iv == iv3d_w)THEN
        CALL damp_latbnd(spec_bdy_width,nlon,nlat,nlev,fieldg)
      ENDIF
      CALL write_var_wrf(ncid(im),iv,1,nlev,fieldg,'3d')
     ENDIF

    END DO

   !END DO  !End do over levels

   IF( iv == iv3d_t )DEALLOCATE(t)

  END DO   !End do over variables

  !WRITE 2D VARIABLES
  !1 GRID AT A TIME
  DO iv=1,nv2d

   DO l=1,ll
    DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL gather_ngrd_mpi(n,v2d(:,im2,iv),fieldg(:,:,1),1)
      END IF
    END DO
    
     im = myrank+1 + (l-1)*nprocs 
     IF (im <= member)THEN
        CALL write_var_wrf(ncid(im),iv,1,1,fieldg(:,:,1),'2d')
     ENDIF
    END DO
  ENDDO


  !CLOSE FILES
  DO l=1,ll
   im = myrank+1 + (l-1)*nprocs
   IF(im <= member) THEN
      CALL close_wrf_file(ncid(im))
    END IF
  END DO


  RETURN
END SUBROUTINE write_ens_mpi
!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  INTEGER :: i,j,m,ilon,ilat


  buf(nij1max,:)=0.0
  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      buf(i,m) = grd(ilon,ilat)
    END DO
  END DO

  RETURN
END SUBROUTINE grd_to_buf
!-----------------------------------------------------------------------
! buffer -> gridded data
!-----------------------------------------------------------------------
SUBROUTINE buf_to_grd(buf,grd)
  REAL(r_sngl),INTENT(IN) :: buf(nij1max,nprocs)
  REAL(r_sngl),INTENT(OUT) :: grd(nlon,nlat)
  INTEGER :: i,j,m,ilon,ilat

  DO m=1,nprocs
    DO i=1,nij1node(m)
      j = m-1 + nprocs * (i-1)
      ilon = MOD(j,nlon) + 1
      ilat = (j-ilon+1) / nlon + 1
      grd(ilon,ilat) = buf(i,m)
    END DO
  END DO

  RETURN
END SUBROUTINE buf_to_grd
!-----------------------------------------------------------------------
! STORING DATA (ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_ensmspr_mpi(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_size) :: v3dm(nij1,nlev,nv3d),tmpv3dm(nij1,nlev,1,nv3d)
  REAL(r_size) :: v2dm(nij1,nv2d),tmpv2dm(nij1,1,nv2d)
  REAL(r_size) :: v3ds(nij1,nlev,nv3d)
  REAL(r_size) :: v2ds(nij1,nv2d)
  REAL(r_sngl) :: fieldg(nlon,nlat)
  INTEGER :: i,k,m,n,iunit
  CHARACTER(20) :: filename='           '

  !COMPUTE AND WRITE ENSEMBLE MEAN (IN NETCDF FORMAT)
  CALL ensmean_grd(member,nij1,v3d,v2d,v3dm,v2dm)
  tmpv3dm(:,:,1,:)=v3dm
  tmpv2dm(:,1,:)=v2dm
  CALL write_ens_mpi(file,1,tmpv3dm,tmpv2dm)

  !COMPUTE AND WRITE ENSEMBLE SPREAD (IN BINARY FORMAT).
  DO n=1,nv3d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO k=1,nlev
      DO i=1,nij1
        v3ds(i,k,n) = (v3d(i,k,1,n)-v3dm(i,k,n))**2
        DO m=2,member
          v3ds(i,k,n) = v3ds(i,k,n) + (v3d(i,k,m,n)-v3dm(i,k,n))**2
        END DO
        v3ds(i,k,n) = SQRT(v3ds(i,k,n) / REAL(member-1,r_size))
      END DO
    END DO
!!$OMP END PARALLEL DO
  END DO
  DO n=1,nv2d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO i=1,nij1
      v2ds(i,n) = (v2d(i,1,n)-v2dm(i,n))**2
      DO m=2,member
        v2ds(i,n) = v2ds(i,n) + (v2d(i,m,n)-v2dm(i,n))**2
      END DO
      v2ds(i,n) = SQRT(v2ds(i,n) / REAL(member-1,r_size))
    END DO
!!$OMP END PARALLEL DO
  END DO

  IF(myrank .EQ. 0)THEN
    iunit=33
    WRITE(filename(1:9),'(A4,A3)') file,'_sp  '
    WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename
    OPEN(unit=iunit,file=filename,form='unformatted',access='sequential')
  ENDIF

  DO n=1,nv3d
    DO k=1,nlev
      CALL gather_ngrd_mpi(0,v3ds(:,k,n),fieldg,1)
      IF(myrank .EQ. 0)WRITE(iunit)fieldg
    ENDDO
  ENDDO
  DO n=1,nv2d
      CALL gather_ngrd_mpi(0,v2ds(:,n),fieldg,1)
      IF(myrank .EQ. 0)WRITE(iunit)fieldg
  ENDDO

  IF(myrank .EQ. 0)CLOSE(iunit)

  RETURN
END SUBROUTINE write_ensmspr_mpi

!-----------------------------------------------------------------------
! STORING DATA (parameter ensemble mean and spread)
!-----------------------------------------------------------------------
SUBROUTINE write_enspmspr_mpi(file,member,vp2d)
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: vp2d(nij1,member,np2d)
  REAL(r_size) :: vp2dm(nij1,np2d)
  REAL(r_size) :: vp2ds(nij1,np2d)
  REAL(r_sngl) :: vp2dg(nlon,nlat)
  INTEGER :: i,k,m,n
  CHARACTER(14) :: filename


  CALL ensmean_ngrd(member,nij1,vp2d,vp2dm,np2d)
  CALL write_ensp_mpi(file,1,vp2d)

  DO n=1,np2d
!!$OMP PARALLEL DO PRIVATE(i,k,m)
    DO i=1,nij1
      vp2ds(i,n) = (vp2d(i,1,n)-vp2dm(i,n))**2
      DO m=2,member
        vp2ds(i,n) = vp2ds(i,n) + (vp2d(i,m,n)-vp2dm(i,n))**2
      END DO
      vp2ds(i,n) = SQRT(vp2ds(i,n) / REAL(member-1,r_size))
    END DO
!!$OMP END PARALLEL DO
  END DO

  DO n=1,np2d
  CALL gather_ngrd_mpi(0,vp2ds(:,n),vp2dg,1)
  IF(myrank == 0) THEN
    filename= file // '_sp'
    OPEN(37,FILE=filename,FORM='unformatted',ACCESS='sequential',POSITION='append')
    WRITE(6,*) 'MYRANK ',myrank,' is appending the parameters to file ',filename
    WRITE(37)vp2dg
  END IF
  ENDDO
  CLOSE(37)

  RETURN
END SUBROUTINE write_enspmspr_mpi


!-----------------------------------------------------------------------
! MPI_ALLREDUCE of hdxf
!-----------------------------------------------------------------------
SUBROUTINE allreduce_obs_mpi(n,nbv,hdxf)
  INTEGER,INTENT(IN) :: n
  INTEGER,INTENT(IN) :: nbv
  REAL(r_size),INTENT(INOUT) :: hdxf(n,nbv)
  REAL(r_size) :: bufs(mpibufsize)
  REAL(r_size) :: bufr(mpibufsize)
  REAL(r_size),ALLOCATABLE :: tmp(:,:)
  INTEGER :: i,j,k
  INTEGER :: iter,niter
  INTEGER :: ierr

  niter = CEILING(REAL(n*nbv)/REAL(mpibufsize))
  ALLOCATE(tmp(mpibufsize,niter))
  bufs=0.0d0
  i=1
  iter=1
  DO k=1,nbv
    DO j=1,n
      bufs(i) = hdxf(j,k)
      i=i+1
      IF(i > mpibufsize) THEN
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_ALLREDUCE(bufs,bufr,mpibufsize,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,ierr)
        tmp(:,iter) = bufr
        i=1
        iter=iter+1
      END IF
    END DO
  END DO
  IF(iter == niter) THEN
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    CALL MPI_ALLREDUCE(bufs,bufr,mpibufsize,MPI_DOUBLE_PRECISION,MPI_SUM,&
      & MPI_COMM_WORLD,ierr)
    tmp(:,iter) = bufr

  END IF

  i=1
  iter=1
  DO k=1,nbv
    DO j=1,n
      hdxf(j,k) = tmp(i,iter)
      i=i+1
      IF(i > mpibufsize) THEN
        i=1
        iter=iter+1
      END IF
    END DO
  END DO
  DEALLOCATE(tmp)

  RETURN
END SUBROUTINE allreduce_obs_mpi



!-----------------------------------------------------------------------
! SMOOTH PARAMETER UPDATE
!-----------------------------------------------------------------------
SUBROUTINE smooth_par_update(p2da,p2dg,member,smooth_par)
  REAL(r_size),INTENT(INOUT) :: p2da(nij1,member,np2d)
  REAL(r_size),INTENT(IN)    :: p2dg(nij1,member,np2d)
  INTEGER     ,INTENT(IN)    :: smooth_par(np2d)       !Wheter parameter will be smoothed or not.
  REAL(r_sngl) :: p2dga(nlon,nlat)
  REAL(r_sngl) :: p2dgg(nlon,nlat)
  REAL(r_size) :: field(nlon,nlat)
  INTEGER      :: mask(nlon,nlat)
  INTEGER :: l,n,ll,im,im2,ip,ns,nr
  INTEGER,INTENT(IN)         :: member     !Ensemble size.
  INTEGER :: iter,niter
  REAL(r_sngl) :: bufs(nij1max,np2d)
  REAL(r_sngl) :: bufr(nij1max,np2d,nprocs)
  INTEGER :: j,k,ierr,np
  CHARACTER(256) :: filter_type
  INTEGER :: filter_option
  INTEGER :: lambda
  lambda=30  !Waves with length less than 30 will be seriously filtered
             !filter length is authomatically computed as 2*lambda.
             !No boundary conditions are assumed. (point beyond boundaries 
             !are treated as missing values)
             !This values is in grid points.
  filter_type='Lanczos'
  mask=1
  mask(nlon,:)=0
  mask(:,nlat)=0

  DO np=1,np2d
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    DO n=0,nprocs-1
      im = n+1 + (l-1)*nprocs
      IF(im <= member) THEN
        CALL gather_ngrd_mpi(n,p2dg(:,im,np),p2dgg,1)
        CALL gather_ngrd_mpi(n,p2da(:,im,np),p2dga,1)
      END IF
    END DO

    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
       field=REAL(p2dga-p2dgg,r_size)   !Parameter update
        IF( smooth_par(np) == 1)THEN
         CALL filter_2d(field,field,mask,lambda,filter_type,nlon,nlat)
        ENDIF
       p2dga=p2dgg+REAL(field,r_sngl) !Parameter smoothed analysis.
    END IF

    DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,p2dga,p2da(:,im2,np),1)
      END IF
    END DO

  END DO
  END DO !End do over the parameters
RETURN
END SUBROUTINE smooth_par_update

!-----------------------------------------------------------------------
! Initialize parameter ensemble.
!----------------------------------------------------------------------
SUBROUTINE init_par_mpi(analp2d,member,update_parameter,param_default_value,param_sprd_ini)
IMPLICIT NONE
INTEGER     ,INTENT(IN)    :: member,update_parameter(np2d)
REAL(r_size),INTENT(IN)    :: param_default_value(np2d),param_sprd_ini(np2d)
REAL(r_size),INTENT(INOUT) :: analp2d(nij1,member,np2d)
INTEGER                    :: ip,ij,ierr,ll,l,im,jj,n
REAL(r_size)               :: MAX_LOCAL_SPREAD,TMP_SPRD,MAX_GLOBAL_SPREAD,TMP_MEAN
REAL(r_size)               :: MIN_LOCAL_SPREAD,MIN_GLOBAL_SPREAD
REAL(r_size)               :: tmpfield(nlon,nlat)
INTEGER                    :: mask(nlon,nlat)
INTEGER                    :: lambda
CHARACTER(256) :: filter_type
  lambda=30  !Waves with length less than 10 will be seriously filtered
             !filter length is authomatically computed as 2*lambda.
             !No boundary conditions are assumed. (point beyond boundaries 
             !are treated as missing values)
             !This values is in grid points.
  filter_type='Lanczos'
  mask=1
  mask(nlon,:)=0
  mask(:,nlat)=0

!DO over the parameters
DO ip=1,np2d
IF(update_parameter(ip)==1)THEN

 !Check the maximum parameter sprad.
  MAX_LOCAL_SPREAD=0.0d0
  MIN_LOCAL_SPREAD=UNDEF
  DO ij=1,nij1
    IF( analp2d(ij,1,ip) .NE. REAL(REAL(UNDEF,r_sngl),r_size) )THEN
    CALL com_stdev(member,analp2d(ij,:,ip),TMP_SPRD)
    IF(TMP_SPRD > MAX_LOCAL_SPREAD)MAX_LOCAL_SPREAD=TMP_SPRD
    IF(TMP_SPRD < MIN_LOCAL_SPREAD)MIN_LOCAL_SPREAD=TMP_SPRD
    ENDIF
  ENDDO
  CALL MPI_ALLREDUCE(MAX_LOCAL_SPREAD,MAX_GLOBAL_SPREAD,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(MIN_LOCAL_SPREAD,MIN_GLOBAL_SPREAD,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierr)

  WRITE(6,*)"MAX / MIN PARAMETER SPREAD IS:",MAX_GLOBAL_SPREAD," ",MIN_GLOBAL_SPREAD

  IF(MAX_GLOBAL_SPREAD .LE. 1.0d-10 )THEN

      WRITE(6,*)"PARAMETERS HAS NOT BEEN INITIALIZED."

      !Initialize parameter perturbations using random filtered fields.
      ll = CEILING(REAL(member)/REAL(nprocs))
      DO l=1,ll
       im = myrank+1 + (l-1)*nprocs
        IF(im <= member) THEN
         tmpfield=0.0d0
          DO jj=1,nlat-1
           CALL com_randn2(nlon,tmpfield(:,jj),myrank)
          ENDDO
           CALL filter_2d(tmpfield,tmpfield,mask,lambda,filter_type,nlon,nlat)

        END IF

         DO n=0,nprocs-1
          im = n+1 + (l-1)*nprocs
          IF(im <= member) THEN
           CALL scatter_ngrd_mpi(n,REAL(tmpfield,r_sngl),analp2d(:,im,ip),1)
          END IF
        END DO
      END DO

  !Adjust the parameter ensemble spread. 
  !WRITE(*,*)myrank,analp2d(10,:,ip)
     DO ij=1,nij1
      CALL com_stdev(member,analp2d(ij,:,ip),TMP_SPRD)
      CALL com_mean(member,analp2d(ij,:,ip),TMP_MEAN)
       IF(TMP_SPRD .ne. 0)analp2d(ij,:,ip)=(analp2d(ij,:,ip)-TMP_MEAN)*param_sprd_ini(ip)/TMP_SPRD 
        analp2d(ij,:,ip)=analp2d(ij,:,ip)+param_default_value(ip)
        !write(*,*)TMP_SPRD,param_sprd_ini(ip)
     END DO
       !WRITE(*,*)myrank,analp2d(10,:,ip)

     !ELSE
     !  !THIS PARAMETER IS NOT BEING ESTIMATED THEN SET TO DEFALUT VALUE.
     !   analp2d(:,:,ip)=param_defalut_value(ip)

  ENDIF
ENDIF
ENDDO

END SUBROUTINE init_par_mpi

!----------------------------------------------------------------------------------
! Gaussian filter
!----------------------------------------------------------------------------------

SUBROUTINE gaussian_filter_mpi(fieldl,fieldls,sigma,enssize)
IMPLICIT NONE
REAL(r_size) :: sigma
INTEGER      :: enssize
REAL(r_size) :: fieldl(nij1,enssize),fieldls(nij1,enssize)
REAL(r_size) :: field(nlon,nlat,enssize),weigth(nlon,nlat)
REAL(r_size) :: dist,mean
REAL(r_sngl) :: bufl(nij1),bufg(nlon,nlat),bufg2(nlon,nlat)
INTEGER      :: i,ierr,ilat,ilon,iens
!Sigma should be given in meters.

field=0.0d0
fieldls=0.0d0
weigth=0.0d0

IF(myrank == 0)bufg=REAL(lon,r_sngl)
CALL MPI_BCAST(bufg,nlon*nlat,MPI_REAL,0,MPI_COMM_WORLD,ierr)
lon=REAL(bufg,r_size)

!CALL MPI_BARRIER()

IF(myrank == 0)bufg=REAL(lat,r_sngl)
CALL MPI_BCAST(bufg,nlon*nlat,MPI_REAL,0,MPI_COMM_WORLD,ierr)
lat=REAL(bufg,r_size)

!CALL MPI_BARRIER()

DO i=1,nij1 !Loop over local grid points.
  !WRITE(*,*)i
  DO ilat=1,nlat-1
    DO ilon=1,nlon-1
     IF( ABS( lat1(i) - lat(ilon,ilat)) < 30)THEN
       CALL com_distll_1(lon1(i),lat1(i),lon(ilon,ilat),lat(ilon,ilat),dist)
       !IF(lon(ilon,ilat) < 0)WRITE(*,*)lon1(i),lat1(i),lon(ilon,ilat),lat(ilon,ilat),dist
       !Compute Gaussian weigth.
       IF(dist .LE. 3*sigma)THEN
         weigth(ilon,ilat)  =  weigth(ilon,ilat)+exp(-dist**2/(2*sigma**2))
         field(ilon,ilat,:)   =  field(ilon,ilat,:)+exp(-dist**2/(2*sigma**2))*fieldl(i,:)
       ENDIF
     ENDIF
    ENDDO
  ENDDO
ENDDO

 DO iens=1,enssize
   bufg=REAL(field(:,:,iens),r_sngl)
   CALL MPI_REDUCE(bufg,bufg2,nlon*nlat,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
   IF(myrank == 0)field(:,:,iens)=REAL(bufg2,r_size)
 ENDDO
 bufg=REAL(weigth,r_sngl)
 CALL MPI_REDUCE(bufg,bufg2,nlon*nlat,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
 IF (myrank== 0)weigth=REAL(bufg2,r_size)
 

 !Compute the weigthed average. 
 DO iens=1,enssize
    IF(myrank==0)THEN
     DO ilon=1,nlon-1
       DO ilat=1,nlat-1
       field(ilon,ilat,iens)=field(ilon,ilat,iens)/weigth(ilon,ilat)
       ENDDO 
     ENDDO
     !Compute field at nlon and nlat.
     mean=0.0d0
     DO ilon=1,nlon-1
       mean=mean + field(ilon,nlat-1,iens)
     ENDDO
       mean=mean/REAL(nlon-1,r_size)
       field(:,nlat,iens)=mean
       !Use global ciclic boundary condition.
       field(nlon,:,iens)=0.5*(field(nlon-1,:,iens)+field(1,:,iens)) 

    bufg=REAL(field(:,:,iens),r_sngl)
    ENDIF

    CALL MPI_SCATTER(bufg,nlat*nlon,MPI_REAL,bufl,nij1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    fieldls(:,iens)=REAL(bufl,r_size)
 ENDDO

 !For debug
 IF(myrank == 0)THEN
 OPEN(33,FILE='smoth_norm.grd',STATUS='unknown',FORM='unformatted',ACCESS='sequential')
 DO iens=1,enssize
    WRITE(33)field(:,:,iens)
 ENDDO
 CLOSE(33)
 OPEN(44,FILE='weight.grd',STATUS='unknown',FORM='unformatted',ACCESS='sequential')
 WRITE(44)weigth
 CLOSE(44)
 ENDIF
 !For debug
  
 !CALL MPI_BARRIER(MPI_COMM_WORLD)

END SUBROUTINE gaussian_filter_mpi

!-----------------------------------------------------------------------
! Write obsope output
!-----------------------------------------------------------------------

SUBROUTINE obsope_io(member,nobs,option,oelm,olon,olat,olev,odat,oerr,otyp,olot  &
                        ,oi,oj,odep,ohdxf )
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: member
  INTEGER , INTENT(IN) :: nobs
  REAL(r_size) , INTENT(INOUT) :: oelm(nobs) , olon(nobs) , olat(nobs) , olev(nobs)
  REAL(r_size) , INTENT(INOUT) :: odat(nobs) , oerr(nobs) , otyp(nobs) , olot(nobs)
  REAL(r_size) , INTENT(INOUT) :: oi(nobs)   , oj(nobs)   , odep(nobs) , ohdxf(nobs,nbv)
  INTEGER                   :: ii , dummy
  INTEGER      , INTENT(IN) :: option ! 1 - read , 2 - write

  REAL(r_sngl) :: fieldg(nlon,nlat)
  INTEGER :: l,n,ll,im,im2,iunit,iv,ilev,ncid(member),ierr,j,k 
  CHARACTER(20) :: filename='obsop00000'

  IF( option == 1 )THEN
    oelm=0.0d0
    olon=0.0d0
    olat=0.0d0
    olev=0.0d0
    odat=0.0d0
    oerr=0.0d0
    otyp=0.0d0
    olot=0.0d0
    oi=0.0d0
    oj=0.0d0
    odep=0.0d0
    ohdxf=0.0d0
  ENDIF
  
  iunit=99
  !OPEN FILES AND READ CONSTANTS.
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(6:10),'(I5.5)') im
      OPEN(iunit,FILE=filename,FORM='unformatted',ACCESS='sequential')
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is writing a file ',filename

         IF( option == 1 )THEN
          READ(iunit)dummy
          DO ii=1,nobs
          READ(iunit)oelm(ii),olon(ii),olat(ii),olev(ii),odat(ii),oerr(ii),otyp(ii),olot(ii),oi(ii) &
           ,oj(ii),odep(ii),ohdxf(ii,im)
          ENDDO
         ELSEIF( option == 2 )THEN
          WRITE(iunit)nobs
          DO ii=1,nobs
          WRITE(iunit)oelm(ii),olon(ii),olat(ii),olev(ii),odat(ii),oerr(ii),otyp(ii),olot(ii),oi(ii) &
           ,oj(ii),odep(ii),ohdxf(ii,im)
          ENDDO
         ENDIF
  
     ENDIF
  ENDDO

  !CLOSE FILES
  CLOSE(iunit)

END SUBROUTINE obsope_io


END MODULE common_mpi_wrf
