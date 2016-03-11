PROGRAM  merge_wrfinput
!Copy the atmosphere data of one wrf input into another wrfinput.
!The grid (horizontal and vertical) of both wrfinput must be the same.
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4) :: nlon,nlat,nlev
  INTEGER(4) :: ncid1,ncid2,varid,dimids(3),shape(3)
  INTEGER(4) :: reclength,irec
  INTEGER(4) :: i,j,k
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  REAL(4),ALLOCATABLE :: dnw(:) 
  REAL(4),ALLOCATABLE :: buf2d(:,:),buf3d(:,:,:)
  CHARACTER(9) :: ncin1  = 'input1.nc'
  CHARACTER(9) :: ncin2  = 'input2.nc'
!
! initialization
!
  CALL check_io(NF_OPEN(ncin1,NF_NOWRITE,ncid1))
  CALL check_io(NF_OPEN(ncin2,NF_WRITE,ncid2))

! Get dimensions from file 1 (we will assume that file 2 has the 
! same dimensions!!!)
  CALL check_io(NF_INQ_VARID(ncid1,'U',varid))
  CALL check_io(NF_INQ_VARDIMID(ncid1,varid,dimids))
  DO i=1,3
    CALL check_io(NF_INQ_DIMLEN(ncid1,dimids(i),shape(i)))
  END DO
  nlon = shape(1)
  nlat = shape(2) + 1
  nlev = shape(3) + 1
  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev 
!
! ALLOCATE valiables
!
  ALLOCATE(buf2d(nlon,nlat))  
  ALLOCATE(buf3d(nlon,nlat,nlev))

!
! COPY netcdf file (2-D valiables)
!
  ALLOCATE(start(3),count(3))
  start = (/ 1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (2D-VAR) ***'
  !!! PSFC
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'PSFC',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'PSFC',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!! MU
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'MU',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'MU',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!! MUB
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'MUB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'MUB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!! T2
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'T2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'T2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!!Q2
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'Q2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'Q2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!!U10
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'U10',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'U10',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!!V10
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'V10',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'V10',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  !!!TH2
  count = (/ nlon-1,nlat-1,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'TH2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'TH2',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon-1,1:nlat-1)))


  WRITE(6,'(a)') '***  END  : READ NETCDF (2D-VAR) ***'
  DEALLOCATE(start,count)

!
! COPY netcdf file (3-D valiables)
!
  ALLOCATE(start(4),count(4))
  start = (/ 1,1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (3D-VAR) ***'

  !!! U
  count = (/ nlon,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'U',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon,1:nlat-1,1:nlev-1)))

  !!! V
  count = (/ nlon-1,nlat,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'V',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat,1:nlev-1)))

  !!! W
  count = (/ nlon-1,nlat-1,nlev,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'W',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  CALL check_io(NF_INQ_VARID(ncid2,'W',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  !!! PB
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'PB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'PB',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! P
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'P',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'P',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! T
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'T',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'T',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! PHB
  count = (/ nlon-1,nlat-1,nlev,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'PHB',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  CALL check_io(NF_INQ_VARID(ncid2,'PHB',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  !!! PH
  count = (/ nlon-1,nlat-1,nlev,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'PH',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  CALL check_io(NF_INQ_VARID(ncid2,'PH',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev)))

  !!! QVAPOR
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'QVAPOR',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'QVAPOR',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! QCLOUD
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'QCLOUD',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'QCLOUD',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! QRAIN
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'QRAIN',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'QRAIN',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! QICE
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'QICE',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'QICE',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  !!! QSNOW
  count = (/ nlon-1,nlat-1,nlev-1,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'QSNOW',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))

  CALL check_io(NF_INQ_VARID(ncid2,'QSNOW',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon-1,1:nlat-1,1:nlev-1)))



  WRITE(6,'(a)') '***  END  : READ NETCDF (3D-VAR) ***'
  DEALLOCATE(start,count)
  CALL check_io(NF_CLOSE(ncid1))
  CALL check_io(NF_CLOSE(ncid2))



  STOP
END PROGRAM merge_wrfinput

!-----------------------------------------------------------------------
SUBROUTINE check_io(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'  
  INTEGER(4),INTENT(IN) :: status
  
  IF(status /= nf_noerr) THEN
    WRITE(6,*) trim(nf_strerror(status))
    STOP 10
  ENDIF
  
  RETURN
END SUBROUTINE check_io
