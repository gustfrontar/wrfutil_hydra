PROGRAM test
!=======================================================================
!
! [HISTORY:]
!   12/02/2017 Juan Ruiz created
!=======================================================================
USE netcdf


IMPLICIT NONE
INTEGER :: nlon , nlat , nlev
INTEGER(4) :: ncid , varid , status , dimid , dimids(3) 
INTEGER :: i
INTEGER(4) ::  shape(3) 
INTEGER,PARAMETER :: r_size=kind(0.0d0)
INTEGER,PARAMETER :: r_dble=kind(0.0d0)
INTEGER,PARAMETER :: r_sngl=kind(0.0e0)
REAL(r_size) , allocatable :: var(:,:,:) 
Character(len=200) :: filename='pp.pe000000.nc' , varname='RHOT' , tmpchar


  nlon=62
  nlat=62
  nlev=60

  status=nf90_open(TRIM(filename),NF90_NOWRITE, ncid)


  allocate(var(nlev,nlon,nlat))

  var=0.0d0

  status=nf90_inq_varid(ncid,TRIM(varname),varid)
  status=nf90_get_var(ncid,varid,var,(/1,1,1,1/),(/nlev,nlon,nlat,1/)) 
  write(*,*)"Read status ",status


  write(*,*)"var in ",var(10,10,10)

  status=nf90_close(ncid)

  status=nf90_open(TRIM(filename),NF90_WRITE, ncid)
  write(*,*)"status open",status,NF90_NOERR

  

  status=nf90_put_var(ncid, varid, var,         &
                               start = (/ 1, 1, 1, 1 /), &
                               count = (/ nlev, nlon, nlat, 1 /))

  write(*,*)"status write",status,NF90_NOERR
 

  status=nf90_close(ncid) 
 
  write(*,*)"status close",status,NF90_NOERR

  write(*,*)"Read the data again"

  status=nf90_open(TRIM(filename),NF90_NOWRITE, ncid)

  var=0.0d0

  status=nf90_inq_varid(ncid,TRIM(varname),varid)
  status=nf90_get_var(ncid,varid,var,(/1,1,1,1/),(/nlev,nlon,nlat,1/))
  write(*,*)"Read status ",status


  write(*,*)"var in ",var(10,10,10)




END PROGRAM test


!SUBROUTINE check_io(status)
!  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
!  INTEGER(4),INTENT(IN) :: status
!
!  IF(status /= nf_noerr) THEN
!    WRITE(6,*) TRIM(nf_strerror(status))
!    STOP 10
!  ENDIF
!
!  RETURN
!END SUBROUTINE check_io

