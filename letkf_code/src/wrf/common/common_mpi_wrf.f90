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

  CALL scatter_sgrd_mpi(0,SNGL(lon),lon1)
  CALL scatter_sgrd_mpi(0,SNGL(lat),lat1)
  CALL scatter_sgrd_mpi(0,SNGL(lonu),lonu1)
  CALL scatter_sgrd_mpi(0,SNGL(latu),latu1)
  CALL scatter_sgrd_mpi(0,SNGL(lonv),lonv1)
  CALL scatter_sgrd_mpi(0,SNGL(latv),latv1)
  CALL scatter_sgrd_mpi(0,SNGL(ri),ri1)
  CALL scatter_sgrd_mpi(0,SNGL(rj),rj1)
  CALL scatter_sgrd_mpi(0,SNGL(phi0),phi1)
  CALL scatter_sgrd_mpi(0,SNGL(landmask),landmask1)

  IF(myrank.eq.0)DEALLOCATE(ri,rj)  

  RETURN
END SUBROUTINE set_common_mpi_wrf

SUBROUTINE scatter_grd_mpi_fast(nrank,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      v3d(:,k,n) = REAL(bufr(1:nij1,j),r_size)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    v2d(:,n) = REAL(bufr(1:nij1,j),r_size)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_fast


SUBROUTINE gather_grd_mpi_fast(nrank,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: j,k,n,ierr,ns,nr

  ns = nij1max * nlevall
  nr = ns
  j=0
  DO n=1,nv3d
    DO k=1,nlev
      j = j+1
      bufs(1:nij1,j) = REAL(v3d(:,k,n),r_sngl)
    END DO
  END DO

  DO n=1,nv2d
    j = j+1
    bufs(1:nij1,j) = REAL(v2d(:,n),r_sngl)
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_GATHER(bufs,ns,MPI_REAL,&
                & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)

  IF(myrank == nrank) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO

    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE gather_grd_mpi_fast

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

!FAST
  ns = nij1max * ndim
  nr = ns
  IF(myrank == nrank) THEN
    DO n=1,ndim
        CALL grd_to_buf(vg(:,:,n),bufs(:,n,:))
    END DO
  END IF

  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  DO n=1,ndim
      v(:,n) = REAL(bufr(1:nij1,n),r_size)
  END DO

  !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_ngrd_mpi


!-----------------------------------------------------------------------
!Scatter_sgrd_mpi to scatter 1 field with the safe scatter
!-----------------------------------------------------------------------
!This routine can be used to improve memory performance.
SUBROUTINE scatter_sgrd_mpi(nrank,vg,v)
  INTEGER,INTENT(IN) :: nrank
  REAL(r_sngl),INTENT(IN) :: vg(nlon,nlat)
  REAL(r_size),INTENT(OUT) :: v(nij1)
  REAL(r_sngl) :: tmp(nij1max,nprocs)
  REAL(r_sngl) :: bufs(nij1max,nprocs)
  REAL(r_sngl) :: bufr(nij1max)
  INTEGER :: i,j,k,n,ierr,ns,nr
  INTEGER :: iter,niter

  v=0.0d0
  bufs=0.0e0
  bufr=0.0e0

  ns = nij1max 
  nr = ns
  IF(myrank == nrank) THEN
    CALL grd_to_buf(vg(:,:),bufs(:,:))
  END IF

  CALL MPI_SCATTER(bufs,ns,MPI_REAL,&
                 & bufr,nr,MPI_REAL,nrank,MPI_COMM_WORLD,ierr)
  v(1:nij1) = REAL(bufr(1:nij1),r_size)

  RETURN
END SUBROUTINE scatter_sgrd_mpi


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

  ns = nij1max * ndim
  nr = ns
  DO n=1,ndim
      bufs(1:nij1,n) = REAL(v(:,n),r_sngl)
  END DO

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
          CALL read_var_wrf(ncid(im),iv,1,nlev,1,fieldg,'3d')
      ENDIF

     DO n=0,nprocs-1
      im2 = n+1 + (l-1)*nprocs
      IF(im2 <= member) THEN
        CALL scatter_ngrd_mpi(n,fieldg,v3d(:,1:nlev,im2,iv),nlev)
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
       CALL read_var_wrf(ncid(im),iv,1,1,1,fieldg,'2d')
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

SUBROUTINE read_ens_mpi_fast(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file !prefix
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: fieldg3(nlon,nlat,nlev,nv3d)
  REAL(r_sngl) :: fieldg2(nlon,nlat,nv2d)
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
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF( im <= member)THEN

      DO iv=1,nv3d
         CALL read_var_wrf(ncid(im),iv,1,nlev,1,fieldg3(:,:,:,iv),'3d')
      ENDDO

      DO iv=1,nv2d
         CALL read_var_wrf(ncid(im),iv,1,1,1,fieldg2(:,:,iv),'2d')
      ENDDO
    ENDIF
  END DO  !End do over variables

  DO l=1,ll
     DO n=0,nprocs-1
       im2 = n+1 + (l-1)*nprocs
       IF(im2 <= member) THEN
        CALL scatter_grd_mpi_fast(n,fieldg3,fieldg2,v3d(:,:,im2,:),v2d(:,im2,:))
       END IF
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
END SUBROUTINE read_ens_mpi_fast


SUBROUTINE read_ens_mpi_alltoall(file,member,v3d,v2d)
  CHARACTER(4),INTENT(IN) :: file !prefix
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(OUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) , ALLOCATABLE :: field3dg(:,:,:,:)
  REAL(r_sngl) , ALLOCATABLE :: field2dg(:,:,:)
  !REAL(r_sngl) :: fieldg3(nlon,nlat,nlev,nv3d)
  !REAL(r_sngl) :: fieldg2(nlon,nlat,nv2d)
  INTEGER :: l,n,ll,im,im2,i,ierr,iv,ilev,mstart,mend,skip
  !REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  !REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  !REAL(r_size),INTENT(INOUT) :: v3d(nij1,nlev,member,nv3d)
  !REAL(r_size),INTENT(INOUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) , ALLOCATABLE :: bufs(:,:,:)
  REAL(r_sngl) , ALLOCATABLE :: bufr(:,:,:)
  !REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  !REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: k,n,j,m,mcount,p
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  CHARACTER(20) :: filename='xxxx00000'
  INTEGER(4) :: ncid

  ll = CEILING(REAL(member)/REAL(nprocs))
  skip = FLOOR( REAL(nprocs) / REAL(member) )
  if ( skip == 0 )skip = 1
  DO l=1,ll
    IF ( MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip + (l-1)*nprocs <= member )THEN
      !This is a reading rank.
      !Read the data from the netcdf file.
      im = (myrank+1)/skip + (l-1)*nprocs
      ALLOCATE( field3dg(nlon,nlat,nlev,nv3d) , field2dg(nlon,nlat,nv2d) )
      WRITE(filename(1:9),'(A4,I5.5)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
      CALL open_wrf_file(filename,'ro',ncid)
      DO iv=1,nv3d
          CALL read_var_wrf(ncid,iv,1,nlev,1,field3dg(:,:,:,iv),'3d')
      END DO
      DO iv=1,nv2d
          CALL read_var_wrf(ncid,iv,1,1,1,field2dg(:,:,iv),'2d')
      END DO
      CALL close_wrf_file(ncid) 
      !Store data into the send bufr
      ALLOCATE( bufs(nij1max,nlevall,nprocs) )
      j=0
      DO n=1,nv3d
         DO k=1,nlev
            j = j+1
            CALL grd_to_buf(field3dg(:,:,k,n),bufs(:,j,:))
         END DO
      END DO
      DO n=1,nv2d
         j = j+1
         CALL grd_to_buf(field2dg(:,:,n),bufs(:,j,:))
      END DO
      DEALLOCATE( field3dg , field2dg )
    ELSE 
      ALLOCATE( bufs(1,1,1) )
    END IF

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    !Set the send and recieve bufr size for each rank
    IF (mstart <= mend) THEN
       mcount = mend - mstart + 1
       ALLOCATE( bufr(nij1max,nlevall,mcount ) )
       IF(mcount > nprocs .OR. mcount <= 0) STOP
       nr  = 0
       nrt = 0
       ns  = 0
       nst = 0
       DO p=1,mcount
         nr(p*skip) = nij1max*nlevall
       ENDDO
       IF (MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount )THEN
          ns(:) = nij1max*nlevall
       END IF
       DO p=2,nprocs
          nrt(p) = nrt(p-1) + nr(p-1)
          nst(p) = nst(p-1) + ns(p-1)
       END DO
       CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                          bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)

       DO m = mstart,mend
          j=0
          DO n=1,nv3d
            DO k=1,nlev
              j = j+1
              v3d(:,k,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)
           END DO
         END DO
         DO n=1,nv2d
           j = j+1
           v2d(:,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)
         END DO
       END DO

    END IF   !END IF mstart < mend
    DEALLOCATE( bufr , bufs )

  END DO  !END DO l
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)



  RETURN
END SUBROUTINE read_ens_mpi_alltoall

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
      CALL read_var_wrf(ncid(im),iv,1,1,1,fieldg,'pa')
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
      WRITE(filename(1:9),'(A4,A5)') file,'emean'
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

SUBROUTINE write_ens_mpi_alltoall(file,member,v3d,v2d)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) , ALLOCATABLE :: field3dg(:,:,:,:)
  REAL(r_sngl) , ALLOCATABLE :: field2dg(:,:,:)
  INTEGER :: ll,im,iv,ncid,mstart,mend
  CHARACTER(20) :: filename='xxxx00000'
  REAL(r_sngl) , ALLOCATABLE :: bufs(:,:,:)
  REAL(r_sngl) , ALLOCATABLE :: bufr(:,:,:)
  !REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  !REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: l,k,n,j,m,mcount,skip,ierr,p
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)


  skip = FLOOR( REAL(nprocs) / REAL(member) )
  if ( skip == 0 )skip = 1
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll

    !Gather the data
    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    mcount = mend - mstart + 1
    IF(mcount > nprocs .OR. mcount <= 0) STOP
    !All member will send data (allocate bufers)
    ALLOCATE( bufs(nij1max,nlevall,mcount ) )
    DO m = mstart,mend
       j=0
       DO n=1,nv3d
         DO k=1,nlev
           j = j+1
           bufs(1:nij1,j,m-mstart+1) = REAL(v3d(:,k,m,n),r_sngl)
         END DO
       END DO
       DO n=1,nv2d
         j = j+1
         bufs(1:nij1,j,m-mstart+1) = REAL(v2d(:,m,n),r_sngl)
       END DO
    END DO
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    nr  = 0
    nrt = 0
    ns  = 0
    nst = 0
    DO p=1,mcount
      ns(p*skip) = nij1max*nlevall
    ENDDO
    IF (MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount)THEN
      ALLOCATE( bufr(nij1max,nlevall,nprocs ) ) 
      nr(:) = nij1max*nlevall
    ELSE
      ALLOCATE( bufr( 1,1,1 ) )
    END IF
    DO p=2,nprocs
      nrt(p) = nrt(p-1) + nr(p-1)
      nst(p) = nst(p-1) + ns(p-1)
    END DO
    WRITE(*,*)'skip',skip
    WRITE(*,*)'bufs',bufs(1,1,:)
    WRITE(*,*)'bufr',bufr(1,1,:)
    WRITE(*,*)'nr',nr
    WRITE(*,*)'ns',ns

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

    CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                       bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
    WRITE(*,*)ierr
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)          
    WRITE(*,*)'bufs',bufs(1,1,:)
    WRITE(*,*)'bufr',bufr(1,1,:)

    DEALLOCATE( bufs )
    IF(MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount) THEN
      !From bufr to grid
      ALLOCATE( field3dg(nlon,nlat,nlev,nv3d) , field2dg(nlon,nlat,nv2d) )
      j=0
      DO n=1,nv3d
        DO k=1,nlev
          j = j+1
          CALL buf_to_grd(bufr(:,j,:),field3dg(:,:,k,n))
        END DO
      END DO
      DO n=1,nv2d
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),field2dg(:,:,n))
      END DO

      !Write data to file
      im = (myrank+1)/skip + (l-1)*nprocs
      IF(member > 1)THEN
        WRITE(filename(1:9),'(A4,I5.5)')file,im
      ELSE
        WRITE(filename(1:9),'(A4,A5)')file,'emean'
      ENDIF
      !OPEN NC FILE
      WRITE(6,'(A,I3.3,2A)')'MYRANK ',myrank,' is writing a file ',filename
      CALL open_wrf_file(filename,'rw',ncid)

      DO iv=1,nv3d
        IF(iv == iv3d_w)THEN
          CALL damp_latbnd(spec_bdy_width,nlon,nlat,nlev,field3dg(:,:,:,iv))
        ENDIF
        CALL write_var_wrf(ncid,iv,1,nlev,field3dg(:,:,:,iv),'3d')
      END DO !End do over 3D variables
  
      DO iv=1,nv2d
        CALL write_var_wrf(ncid,iv,1,1,field2dg(:,:,iv),'2d')
      END DO !End do over 2D variables
      CALL close_wrf_file(ncid)
      DEALLOCATE( field3dg , field2dg , bufr )
    ELSE
      DEALLOCATE( bufr )
    END IF
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END DO

  RETURN
END SUBROUTINE write_ens_mpi_alltoall


!-----------------------------------------------------------------------
! gridded data -> buffer
!-----------------------------------------------------------------------
SUBROUTINE grd_to_buf(grd,buf)
  REAL(r_sngl),INTENT(IN) :: grd(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: buf(nij1max,nprocs)
  !REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
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
  REAL(r_size) :: v3dm(nij1,nlev,1,nv3d)
  REAL(r_size) :: v2dm(nij1,1,nv2d)

  !COMPUTE AND WRITE ENSEMBLE MEAN (IN NETCDF FORMAT)
  v3dm(:,:,1,:)=SUM( v3d , DIM=3 ) / REAL( member , r_size )
  v2dm(:,1,:)=SUM( v2d , DIM=2 ) / REAL( member , r_size )
  CALL write_ens_mpi_alltoall(file,1,v3dm,v2dm)

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
SUBROUTINE init_par_mpi(p2d,member,update_parameter_2d,update_parameter_0d,param_default_value,param_sprd_ini)
IMPLICIT NONE
INTEGER     ,INTENT(IN)    :: member,update_parameter_2d(np2d),update_parameter_0d(np2d)
REAL(r_size),INTENT(IN)    :: param_default_value(np2d),param_sprd_ini(np2d)
REAL(r_size),INTENT(INOUT) :: p2d(nij1,member,np2d)
INTEGER                    :: ip,ij,ierr,ll,l,im,jj,n , i , m
REAL(r_size)               :: MAX_LOCAL_SPREAD,TMP_SPRD,MAX_GLOBAL_SPREAD,TMP_MEAN
REAL(r_size)               :: MIN_LOCAL_SPREAD,MIN_GLOBAL_SPREAD
REAL(r_size)               :: tmpfield(nlon,nlat) , tmpp0d(member) , tmpsprdp0d , tmpmeanp0d
INTEGER                    :: mask(nlon,nlat)
INTEGER                    :: lambda
REAL(r_size)               :: tmpwork(member)
INTEGER                    :: tmpiwork , npoints
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
IF( update_parameter_2d(ip) == 1 )THEN

 !Check the maximum parameter sprad.
  MAX_LOCAL_SPREAD=0.0d0
  MIN_LOCAL_SPREAD=UNDEF
  DO ij=1,nij1
    IF( p2d(ij,1,ip) .NE. REAL(REAL(UNDEF,r_sngl),r_size) )THEN
    CALL com_stdev(member,p2d(ij,:,ip),TMP_SPRD)
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
           CALL scatter_ngrd_mpi(n,REAL(tmpfield,r_sngl),p2d(:,im,ip),1)
          END IF
        END DO
      END DO

  !Adjust the parameter ensemble spread. 
     DO ij=1,nij1
      CALL com_stdev(member,p2d(ij,:,ip),TMP_SPRD)
      CALL com_mean(member,p2d(ij,:,ip),TMP_MEAN)
       IF(TMP_SPRD .ne. 0)p2d(ij,:,ip)=(p2d(ij,:,ip)-TMP_MEAN)*param_sprd_ini(ip)/TMP_SPRD
        p2d(ij,:,ip)=p2d(ij,:,ip)+param_default_value(ip)
     END DO

  ENDIF

ELSEIF( update_parameter_0d(ip) == 1 )THEN

!  !First spatially average 2d parameters.
  tmpp0d=0.0d0
  npoints=0
  DO ij=1,nij1
    IF( landmask1(ij) == 0 .and. p2d(ij,1,ip) /= REAL(REAL(undef,r_sngl),r_size) )THEN
      tmpp0d(:)=tmpp0d(:)+p2d(ij,:,ip)
      npoints=npoints+1
    ENDIF
  ENDDO
!
!  !All the processes will have the same spatially averaged parameters.
!  
  CALL MPI_ALLREDUCE(tmpp0d,tmpwork,member,MPI_DOUBLE_PRECISION,MPI_SUM,&
          & MPI_COMM_WORLD,ierr)
  CALL MPI_ALLREDUCE(npoints,tmpiwork,1,MPI_INTEGER,MPI_SUM,&
          & MPI_COMM_WORLD,ierr)
  !Get the parameter mean.
  tmpp0d=tmpwork/REAL(tmpiwork,r_size)
  CALL com_stdev(member,tmpp0d,tmpsprdp0d)

   WRITE(6,*)tmpsprdp0d

   IF ( tmpsprdp0d <= 1.0d-10 )THEN
      WRITE(6,*)"INITIALIZING 0D PARAMETERS"
      IF( myrank == 0 )then
       CALL com_randn(member,tmpp0d)
       CALL com_stdev(member,tmpp0d,tmpsprdp0d)
       CALL com_mean(member,tmpp0d,tmpmeanp0d)
       DO i=1,member
         tmpp0d(i)=(tmpp0d(i)-tmpmeanp0d)*param_sprd_ini(ip)/tmpsprdp0d
         tmpp0d(i)=tmpp0d(i)+param_default_value(ip)
       ENDDO
      ENDIF

      CALL MPI_BCAST( tmpp0d, member, MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         DO m=1,member
            p2d(:,m,ip)=tmpp0d(m)
        ENDDO

  ENDIF ![End if over tmpsprdp0d == 0 initializatio of 0d parameters ]
   WRITE(6,*)tmpsprdp0d

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

  INTEGER :: l,n,ll,im,im2,iunit,iv,ilev,ncid(member),ierr,j,k 
  CHARACTER(20) :: filename='obsopeMMMMM'

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
      WRITE(filename(7:11),'(I5.5)') im
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


SUBROUTINE scatter_grd_mpi_alltoall(member,mstart,mend,v3dg,v2dg,v3d,v2d)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_size),INTENT(INOUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(INOUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: k,n,j,m,mcount,ierr,p,skip
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  skip = FLOOR( REAL(nprocs) / REAL(nbv) )
  if ( skip == 0 )skip = 1
  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP 

  IF ( MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount )THEN
  !IF(myrank < mcount) THEN  !Ver cual seria la condicion aca.
  !WRITE(*,*)'Estoy en ',myrank,'completo el bufr'
    j=0  
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL grd_to_buf(v3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      CALL grd_to_buf(v2dg(:,:,n),bufs(:,j,:))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  nr  = 0
  nrt = 0
  ns  = 0
  nst = 0
  DO p=1,mcount
    nr(p*skip) = nij1max*nlevall
  ENDDO
  IF (MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount )THEN
    !WRITE(*,*)'Estoy en ',myrank,' y seteo todo para enviar datos'
    !IF(myrank+1 == p) THEN   !TODO ver cual seria la condicion aca.
    ns(:) = nij1max*nlevall
  END IF
  DO p=2,nprocs
    nrt(p) = nrt(p-1) + nr(p-1)
    nst(p) = nst(p-1) + ns(p-1)
  END DO
  !WRITE(6,*)'nr',nr
  !WRITE(6,*)'ns',ns
      
  !CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs,nr,nrt,ns,nst)
  CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                     bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)

  DO m = mstart,mend
    !WRITE(6,*)'bufr',bufr(1,1,(m-mstart+1)*skip)
    !WRITE(6,*)'bufrcompleto',bufr(1,1,:)
    j=0  
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        v3d(:,k,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)  
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      v2d(:,m,n) = REAL(bufr(1:nij1,j,m-mstart+1),r_size)
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-----------------------------------------------------------------------

SUBROUTINE gather_grd_mpi_alltoall(member,mstart,mend,v3d,v2d,v3dg,v2dg)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_size),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_size),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs) 
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs) 
  INTEGER :: k,n,j,m,mcount,skip,ierr,p
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  skip = FLOOR( REAL(nprocs) / REAL(member) )
  if ( skip == 0 )skip = 1
  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP

  DO m = mstart,mend
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = REAL(v3d(:,k,m,n),r_sngl)
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      bufs(1:nij1,j,m-mstart+1) = REAL(v2d(:,m,n),r_sngl)
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  nr  = 0
  nrt = 0
  ns  = 0
  nst = 0
  DO p=1,mcount
    ns(p*skip) = nij1max*nlevall
  ENDDO
  IF (MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount)THEN
    nr(:) = nij1max*nlevall
  END IF
  DO p=2,nprocs
    nrt(p) = nrt(p-1) + nr(p-1)
    nst(p) = nst(p-1) + ns(p-1)
  END DO
  !WRITE(*,*)'nr',nr
  !WRITE(*,*)'ns',ns

  CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                     bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)

  !WRITE(*,*)'bufr',bufr(1,1,:)
  IF(MOD( myrank + 1 , skip ) == 0 .AND. (myrank+1)/skip <= mcount) THEN
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),v3dg(:,:,k,n))
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      CALL buf_to_grd(bufr(:,j,:),v2dg(:,:,n))
    END DO
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
!SUBROUTINE set_alltoallv_counts(mcount,ngpblock,np,n_ens,nt_ens,n_mem,nt_mem)
!  INTEGER,INTENT(IN) :: mcount,ngpblock
!  INTEGER,INTENT(IN) :: np
!  INTEGER,INTENT(OUT) :: n_ens(np),nt_ens(np),n_mem(np),nt_mem(np)
!  INTEGER :: p

!  n_ens = 0
!  nt_ens = 0
!  n_mem = 0
!  nt_mem = 0
!  DO p=1,mcount
!    n_ens(p) = ngpblock
!    IF(myrank+1 == p) THEN
!      n_mem(:) = ngpblock
!    END IF
!  END DO
!  DO p=2,np
!    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
!    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
!  END DO
!
!  RETURN
!END SUBROUTINE set_alltoallv_counts

END MODULE common_mpi_wrf
