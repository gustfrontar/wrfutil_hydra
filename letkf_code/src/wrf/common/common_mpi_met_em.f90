MODULE common_mpi_met_em
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
  USE common_met_em
  USE common_namelist_met_em
  IMPLICIT NONE
  PUBLIC

  INTEGER,PARAMETER :: mpibufsize=10000 !200
  INTEGER,SAVE :: nij1
  INTEGER,SAVE :: nij1max
  INTEGER,ALLOCATABLE,SAVE :: nij1node(:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: ri1(:),rj1(:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: wmatrix(:,:)

CONTAINS
SUBROUTINE set_common_mpi_met_em
  INTEGER :: i,n,ierr

  WRITE(6,'(A)') 'Hello from set_common_mpi_met_em'

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

  ALLOCATE(ri1(nij1))
  ALLOCATE(rj1(nij1))

  CALL scatter_sgrd_mpi(0,SNGL(ri),ri1)
  CALL scatter_sgrd_mpi(0,SNGL(rj),rj1)

  IF(myrank.eq.0)DEALLOCATE(ri,rj)  

  ALLOCATE( wmatrix(nbv_ori,nbv_tar) )
  if( myrank .eq. 0 ) then
     wmatrix = 0.0e0
     !TODO call the transformation matrix generation.
     !CALL get_weigths( nbv_ori , nbv_tar , method , wmatrix )
  endif
  CALL MPI_BCAST(wmatrix,nbv_ori*nbv_tar,MPI_REAL,0,MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE set_common_mpi_met_em

!-----------------------------------------------------------------------
!Scatter_sgrd_mpi to scatter 1 field with the safe scatter
!-----------------------------------------------------------------------
!This routine can be used to improve memory performance.
SUBROUTINE scatter_sgrd_mpi(nrank,vg,v)
  INTEGER,INTENT(IN)       :: nrank
  REAL(r_sngl),INTENT(IN)  :: vg(nlon,nlat)
  REAL(r_sngl),INTENT(OUT) :: v(nij1)
  REAL(r_sngl) :: bufs(nij1max,nprocs)
  REAL(r_sngl) :: bufr(nij1max)
  INTEGER :: ierr,ns,nr

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
  v(1:nij1) = bufr(1:nij1)

  RETURN
END SUBROUTINE scatter_sgrd_mpi

SUBROUTINE read_ens_mpi_alltoall(file,member,v3d,v2d,vs3d)
  CHARACTER(4),INTENT(IN) :: file !prefix
  INTEGER,INTENT(IN) :: member
  REAL(r_sngl),INTENT(OUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_sngl),INTENT(OUT) :: vs3d(nij1,nlev_soil,member,ns3d)
  REAL(r_sngl),INTENT(OUT) :: v2d(nij1,member,nv2d)

  REAL(r_sngl)  :: fields3(nlon,nlat,nlev_soil,ns3d)
  REAL(r_sngl)  :: fieldg3(nlon,nlat,nlev,nv3d)
  REAL(r_sngl)  :: fieldg2(nlon,nlat,nv2d)
  INTEGER       :: l,ll,im,mstart,mend
  CHARACTER(20) :: filename='xxxx00000'

  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:9),'(A4,I5.5)') file,im
      WRITE(6,'(A,I3.3,2A)') 'MYRANK ',myrank,' is reading a file ',filename
      !Read the ncfile
      CALL read_grd(filename,fieldg3,fieldg2,fields3)
    END IF
  END DO  !End do over variables

  !! Agregado por MPM
  DO l=1,ll 
    mstart = 1 + (l-1)*nprocs
    mend = MIN(l*nprocs, member)
    if (mstart <= mend) then
      CALL scatter_grd_mpi_alltoall(member,mstart,mend,fieldg3,fieldg2,fields3,v3d,v2d,vs3d)
    end if
  END DO

  RETURN
END SUBROUTINE read_ens_mpi_alltoall


SUBROUTINE write_ens_mpi_alltoall(file,member,v3d,v2d,vs3d)
  IMPLICIT NONE
  CHARACTER(4),INTENT(IN) :: file
  INTEGER,INTENT(IN) :: member
  REAL(r_sngl),INTENT(IN) :: vs3d(nij1,nlev_soil,member,ns3d)
  REAL(r_sngl),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl)  :: fields3(nlon,nlat,nlev_soil,ns3d)
  REAL(r_sngl)  :: fieldg3(nlon,nlat,nlev,nv3d)
  REAL(r_sngl)  :: fieldg2(nlon,nlat,nv2d)
  INTEGER       :: l,ll,im,mstart,mend
  CHARACTER(20) :: filename='xxxx00000'

! Junta de todos los nodos de procesos a los nodos que escriben todas las 
! variables y todos los niveles
  ll = CEILING(REAL(member)/REAL(nprocs))
  DO l=1,ll
      mstart = 1 + (l-1)*nprocs
      mend = MIN(l*nprocs, member)
      if (mstart <= mend) then
        CALL gather_grd_mpi_alltoall(member,mstart,mend,v3d,v2d,vs3d,fieldg3,fieldg2,fields3) !cambiamos tmpv3d
      end if
  END DO

  !Write the data to the ensemble files.
  DO l=1,ll
    im = myrank+1 + (l-1)*nprocs
    IF(im <= member) THEN
      WRITE(filename(1:9),'(A4,I5.5)')file,im
      !Write a netcdf file
      WRITE(6,'(A,I3.3,2A)')'MYRANK ',myrank,' is writing a file ',filename
      CALL write_grd(filename,fieldg3,fieldg2,fields3)
    END IF
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

SUBROUTINE scatter_grd_mpi_alltoall(member,mstart,mend,v3dg,v2dg,vs3dg,v3d,v2d,vs3d)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_sngl),INTENT(IN) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(IN) :: vs3dg(nlon,nlat,nlev_soil,ns3d)
  REAL(r_sngl),INTENT(IN) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl),INTENT(INOUT) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_sngl),INTENT(INOUT) :: vs3d(nij1,nlev_soil,member,ns3d)
  REAL(r_sngl),INTENT(INOUT) :: v2d(nij1,member,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP 

  IF(myrank < mcount) THEN 
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
    DO n=1,ns3d
      DO k=1,nlev_soil
         j = j+1
         CALL grd_to_buf(vs3dg(:,:,k,n),bufs(:,j,:))
      END DO
    END DO
  END IF


  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF(mcount == nprocs) THEN 
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL, &
                      bufr, nij1max*nlevall, MPI_REAL, MPI_COMM_WORLD, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs,nr,nrt,ns,nst)
    CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                       bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)

  END IF

  DO m = mstart,mend
    j=0  
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        v3d(:,k,m,n) = bufr(1:nij1,j,m-mstart+1)
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      v2d(:,m,n) = bufr(1:nij1,j,m-mstart+1)
    END DO
    DO n=1,ns3d
      DO k=1,nlev_soil
        j = j+1
        vs3d(:,k,m,n) = bufr(1:nij1,j,m-mstart+1)
      END DO
    END DO
  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

  RETURN
END SUBROUTINE scatter_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Gather gridded data using MPI_ALLTOALL(V) (all -> mstart~mend)
!-----------------------------------------------------------------------

SUBROUTINE gather_grd_mpi_alltoall(member,mstart,mend,v3d,v2d,vs3d,v3dg,v2dg,vs3dg)
  INTEGER,INTENT(IN) :: mstart,mend,member
  REAL(r_sngl),INTENT(IN) :: v3d(nij1,nlev,member,nv3d)
  REAL(r_sngl),INTENT(IN) :: vs3d(nij1,nlev_soil,member,ns3d)
  REAL(r_sngl),INTENT(IN) :: v2d(nij1,member,nv2d)
  REAL(r_sngl),INTENT(OUT) :: v3dg(nlon,nlat,nlev,nv3d)
  REAL(r_sngl),INTENT(OUT) :: vs3dg(nlon,nlat,nlev_soil,ns3d)
  REAL(r_sngl),INTENT(OUT) :: v2dg(nlon,nlat,nv2d)
  REAL(r_sngl) :: bufs(nij1max,nlevall,nprocs)
  REAL(r_sngl) :: bufr(nij1max,nlevall,nprocs)
  INTEGER :: k,n,j,m,mcount,ierr
  INTEGER :: ns(nprocs),nst(nprocs),nr(nprocs),nrt(nprocs)

  mcount = mend - mstart + 1
  IF(mcount > nprocs .OR. mcount <= 0) STOP

  DO m = mstart,mend
    j=0
    DO n=1,nv3d
      DO k=1,nlev
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = v3d(:,k,m,n)
      END DO
    END DO
    DO n=1,nv2d
      j = j+1
      bufs(1:nij1,j,m-mstart+1) = v2d(:,m,n)
    END DO
    DO n=1,ns3d
      DO k=1,nlev_soil
        j = j+1
        bufs(1:nij1,j,m-mstart+1) = vs3d(:,k,m,n)
      END DO
    END DO

  END DO

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF(mcount == nprocs) THEN
    CALL MPI_ALLTOALL(bufs, nij1max*nlevall, MPI_REAL, &
                      bufr, nij1max*nlevall, MPI_REAL, MPI_COMM_WORLD, ierr)
  ELSE
    CALL set_alltoallv_counts(mcount,nij1max*nlevall,nprocs,ns,nst,nr,nrt)
    CALL MPI_ALLTOALLV(bufs, ns, nst, MPI_REAL, &
                       bufr, nr, nrt, MPI_REAL, MPI_COMM_WORLD, ierr)
  END IF


  IF(myrank < mcount) THEN
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
    DO n=1,ns3d
      DO k=1,nlev_soil
        j = j+1
        CALL buf_to_grd(bufr(:,j,:),vs3dg(:,:,k,n))
      END DO
    END DO    
  END IF

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  RETURN
END SUBROUTINE gather_grd_mpi_alltoall

!-----------------------------------------------------------------------
! Set the send/recieve counts of MPI_ALLTOALLV
!-----------------------------------------------------------------------
SUBROUTINE set_alltoallv_counts(mcount,ngpblock,np,n_ens,nt_ens,n_mem,nt_mem)
  INTEGER,INTENT(IN) :: mcount,ngpblock
  INTEGER,INTENT(IN) :: np
  INTEGER,INTENT(OUT) :: n_ens(np),nt_ens(np),n_mem(np),nt_mem(np)
  INTEGER :: p

  n_ens = 0
  nt_ens = 0
  n_mem = 0
  nt_mem = 0
  DO p=1,mcount
    n_ens(p) = ngpblock
    IF(myrank+1 == p) THEN
      n_mem(:) = ngpblock
    END IF
  END DO
  DO p=2,np
    nt_ens(p) = nt_ens(p-1) + n_ens(p-1)
    nt_mem(p) = nt_mem(p-1) + n_mem(p-1)
  END DO

  RETURN
END SUBROUTINE set_alltoallv_counts

END MODULE common_mpi_met_em
