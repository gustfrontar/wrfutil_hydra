PROGRAM  merge_metem
!Complete met_em files with soil data. (Letkf is not currently using this, so we can take it from gdas)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4) :: nlon,nlat,nlev
  INTEGER(4) :: ncid1,ncid2,varid,dimids(3),shape(3)
  INTEGER(4) :: reclength,irec
  INTEGER(4) :: i,j,k
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  REAL(4),ALLOCATABLE :: dnw(:) 
  REAL(4),ALLOCATABLE :: buf2d(:,:),buf3d(:,:,:)
  CHARACTER(9) :: ncin1  = 'input'
  CHARACTER(9) :: ncin2  = 'output'
!
! initialization
!
  CALL check_io(NF_OPEN(ncin1,NF_NOWRITE,ncid1))
  CALL check_io(NF_OPEN(ncin2,NF_WRITE,ncid2))

! Get dimensions from file 1 (we will assume that file 2 has the 
! same dimensions!!!)
  CALL check_io(NF_INQ_VARID(ncid1,'ST',varid))
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


  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'GHT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'GHT',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SNOWH',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SNOWH',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SNOW',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SNOW',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SKINTEMP',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SKINTEMP',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SOILHGT',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SOILHGT',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'LANDSEA',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'LANDSEA',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SEAICE',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SEAICE',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'ST100200',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'ST100200',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'ST040100',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'ST040100',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'ST010040',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'ST010040',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'ST000010',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'ST000010',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SM100200',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SM100200',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SM040100',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SM040100',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SM010040',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SM010040',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  count = (/ nlon,nlat,1 /)
  buf2d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SM000010',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf2d(1:nlon,1:nlat)))

  CALL check_io(NF_INQ_VARID(ncid2,'SM000010',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf2d(1:nlon,1:nlat)))

  WRITE(6,'(a)') '***  END  : READ NETCDF (2D-VAR) ***'
  DEALLOCATE(start,count)

!
! COPY netcdf file (3-D valiables)
!
  ALLOCATE(start(4),count(4))
  start = (/ 1,1,1,1 /)
  WRITE(6,'(a)') '*** START : READ NETCDF (3D-VAR) ***'

  count = (/ nlon,nlat,nlev,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'ST',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon,1:nlat,1:nlev)))

  CALL check_io(NF_INQ_VARID(ncid2,'ST',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon,1:nlat,1:nlev)))


  count = (/ nlon,nlat,nlev,1 /)
  buf3d = 0.0
  CALL check_io(NF_INQ_VARID(ncid1,'SM',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid1,varid,start,count,buf3d(1:nlon,1:nlat,1:nlev)))

  CALL check_io(NF_INQ_VARID(ncid2,'SM',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid2,varid,start,count,buf3d(1:nlon,1:nlat,1:nlev)))


  WRITE(6,'(a)') '***  END  : READ NETCDF (3D-VAR) ***'
  DEALLOCATE(start,count)
  CALL check_io(NF_CLOSE(ncid1))
  CALL check_io(NF_CLOSE(ncid2))



STOP
END PROGRAM merge_metem

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
