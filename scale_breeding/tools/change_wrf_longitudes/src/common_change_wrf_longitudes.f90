MODULE common_change_wrf_longitudes
!=======================================================================
!
! [PURPOSE:] Common Information for WRF
!
! [HISTORY:]
!   10/15/2004 Takemasa Miyoshi  created
!   01/23/2009 Takemasa Miyoshi  modified
!   07/14/2010 Takemasa Miyoshi  modified for WRF
!   TODO: This routine is for NCEP data, other input sources may produce
!   different variables inside met_em files. In that case the routine 
!   has to be able to process the data including the new variables and
!   automatically detecting the abssence of some other variables.
!=======================================================================
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
  INTEGER,PARAMETER :: r_size=kind(0.0d0)
  INTEGER,PARAMETER :: r_dble=kind(0.0d0)
  INTEGER,PARAMETER :: r_sngl=kind(0.0e0)

  !MODULE VARIABLES (NOT READ FROM NAMELIST)
  INTEGER, PARAMETER :: nv3d=15
  INTEGER, PARAMETER :: nv2d=3
  INTEGER, PARAMETER :: np2d=3
  INTEGER, PARAMETER :: nid_obs=10

  REAL(r_sngl) :: tmpatt
  INTEGER,SAVE :: dx,dy ! grid spacing [m]
  REAL(r_size),SAVE :: pdx,pdy
  INTEGER,SAVE :: nlon
  INTEGER,SAVE :: nlat
  INTEGER,SAVE :: nlev

  REAL(4),SAVE      :: p0 != 1.0e+5
  REAL(4),PARAMETER :: t0 = 300.0
  REAL(r_size)      :: p_top
  INTEGER(4),PARAMETER :: imt_grd = 50
  INTEGER(4),PARAMETER :: spec_bdy_width = 5
  REAL(r_size),ALLOCATABLE,SAVE :: lon(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lat(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonu(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latu(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: lonv(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: latv(:,:)
!  REAL(r_size),ALLOCATABLE,SAVE :: fcori(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: phi0(:,:)
  REAL(r_size),ALLOCATABLE,SAVE :: landmask(:,:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: dnw(:)
  REAL(r_sngl),ALLOCATABLE,SAVE :: znu(:)

  CHARACTER(19) :: hdate !Date in input files.

CONTAINS

!-----------------------------------------------------------------------
! Read in longitudes
!-----------------------------------------------------------------------
SUBROUTINE set_common_wrf(inputfile)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  REAL(r_sngl),ALLOCATABLE :: buf4(:,:)
  INTEGER :: i
  INTEGER(4) :: ncid,varid,dimids(3),shape(3)
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  CHARACTER(*)           :: inputfile
  INTEGER :: ierr

  WRITE(6,'(A)') 'Hello from set_common_wrf'

  WRITE(6,'(A)') 'READING: ',inputfile
  CALL check_io(NF_OPEN(inputfile,NF_NOWRITE,ncid))
  CALL check_io(NF_INQ_VARID(ncid,'U',varid))
  CALL check_io(NF_INQ_VARDIMID(ncid,varid,dimids))
  DO i=1,3
    CALL check_io(NF_INQ_DIMLEN(ncid,dimids(i),shape(i)))
  END DO
  nlon = shape(1)
  nlat = shape(2) + 1
  nlev = shape(3) + 1
  WRITE(6,'(A)') '*** grid information ***'
  WRITE(6,'(3(2X,A,I5))') 'nlon =',nlon,'nlat =',nlat,'nlev =',nlev

  ALLOCATE(buf4(nlon,nlat))
  ALLOCATE(lon(nlon,nlat))
  ALLOCATE(lonu(nlon,nlat))
  ALLOCATE(lonv(nlon,nlat))

  !!! Time
  ALLOCATE(start(3),count(3))

  start = (/ 1,1,1 /)
  !!! XLONG
  count = (/ nlon-1,nlat-1,1 /)
  lon=0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  lon(1:nlon-1,1:nlat-1)=REAL(buf4(1:nlon-1,1:nlat-1),r_size)
  !!! XLONG_U
  count = (/ nlon,nlat-1,1 /)
  lonu = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_U',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon,1:nlat-1)))
  lonu(1:nlon,1:nlat-1)=REAL(buf4(1:nlon,1:nlat-1),r_size)
  !!! XLONG_V
  count = (/ nlon-1,nlat,1 /)
  lonv = 0.0d0
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_V',varid))
  CALL check_io(NF_GET_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat)))
  lonv(1:nlon-1,1:nlat)=REAL(buf4(1:nlon-1,1:nlat),r_size)

  DEALLOCATE(start,count)

  WRITE(6,'(a)') 'DOMAIN CORNERS'
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,nlat-1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,nlat-1)


  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'

  CALL check_io(NF_CLOSE(ncid))

  RETURN
END SUBROUTINE set_common_wrf

!-----------------------------------------------------------------------
! Modify longitudes
!-----------------------------------------------------------------------
SUBROUTINE modify_longitudes(inputfile)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  REAL(r_sngl),ALLOCATABLE :: buf4(:,:)
  INTEGER :: i
  INTEGER(4) :: ncid,varid
  INTEGER(4),ALLOCATABLE :: start(:),count(:)
  CHARACTER(*)           :: inputfile
  INTEGER :: ierr , ii , jj

  DO ii=1,nlon
    DO jj=1,nlat
       if( lon(ii,jj) < 0.0d0 )then
         lon(ii,jj)=360.0d0 + lon(ii,jj)
       endif
       if( lonu(ii,jj) < 0.0d0 )then
         lonu(ii,jj)=360.0d0 + lonu(ii,jj)
       endif
       if( lonv(ii,jj) < 0.0d0 )then
         lonv(ii,jj)=360.0d0 + lonv(ii,jj)
       endif
    ENDDO
  ENDDO


  WRITE(6,'(A)') 'READING: ',inputfile
  CALL check_io(NF_OPEN(inputfile,NF_WRITE,ncid))

  ALLOCATE(buf4(nlon,nlat))


  ALLOCATE(start(3),count(3))

  start = (/ 1,1,1 /)
  !!! XLONG
  count = (/ nlon-1,nlat-1,1 /)
  buf4=REAL(lon(1:nlon-1,1:nlat-1),r_sngl)
  CALL check_io(NF_INQ_VARID(ncid,'XLONG',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat-1)))
  !!! XLONG_U
  count = (/ nlon,nlat-1,1 /)
  buf4=REAL(lonu(1:nlon,1:nlat-1),r_sngl)
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_U',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,buf4(1:nlon,1:nlat-1)))
  !!! XLONG_V
  count = (/ nlon-1,nlat,1 /)
  buf4 = REAL(lonv(1:nlon-1,1:nlat),r_sngl)
  CALL check_io(NF_INQ_VARID(ncid,'XLONG_V',varid))
  CALL check_io(NF_PUT_VARA_REAL(ncid,varid,start,count,buf4(1:nlon-1,1:nlat)))

  DEALLOCATE(start,count)

  WRITE(6,'(a)') 'DOMAIN CORNERS'
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(1,nlat-1)
  WRITE(6,'(A4,F7.2,A4,F7.2)')"LON",lon(nlon-1,nlat-1)


  WRITE(6,'(a)') '***  END  : READ NETCDF (CNST) ***'

  CALL check_io(NF_CLOSE(ncid))

  RETURN
END SUBROUTINE modify_longitudes


SUBROUTINE check_io(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER(4),INTENT(IN) :: status

  IF(status /= nf_noerr) THEN
    WRITE(6,*) TRIM(nf_strerror(status))
    STOP 10
  ENDIF

  RETURN
END SUBROUTINE check_io

SUBROUTINE open_wrf_file(filename,mode,ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(OUT)   :: ncid
  CHARACTER(*) , INTENT(IN) :: filename,mode
  INCLUDE 'netcdf.inc'

  IF(mode .eq. 'ro' )THEN
  CALL check_io(NF_OPEN(TRIM(filename),NF_NOWRITE,ncid))
  ELSEIF(mode .eq. 'rw' )THEN
  CALL check_io(NF_OPEN(TRIM(filename),NF_WRITE,ncid))
  END IF

  RETURN
END SUBROUTINE open_wrf_file

SUBROUTINE close_wrf_file(ncid)
  IMPLICIT NONE
  INTEGER(4), INTENT(IN)   :: ncid
  INCLUDE 'netcdf.inc'

  CALL  check_io(NF_CLOSE(ncid))

  RETURN
END SUBROUTINE close_wrf_file

!-----------------------------------------------------------------------
! Get time nc
!-----------------------------------------------------------------------
SUBROUTINE get_date(inputfile,date)
IMPLICIT NONE
CHARACTER(19),INTENT(OUT) :: date
CHARACTER(*) ,INTENT(IN)  :: inputfile
INTEGER :: ncid , start(2) , count(2) , varid
INCLUDE 'netcdf.inc'

  CALL check_io(NF_OPEN(inputfile,NF_NOWRITE,ncid))

  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  start = (/  1,1 /)
  count = (/ 19,1 /)
  CALL check_io(NF_GET_VARA_TEXT(ncid,varid,start,count,date))

  CALL check_io(NF_CLOSE(ncid))

END SUBROUTINE get_date

!-----------------------------------------------------------------------
! Put time nc
!-----------------------------------------------------------------------

SUBROUTINE put_date(inputfile,date)
IMPLICIT NONE
CHARACTER(19), INTENT(IN) :: DATE
CHARACTER(*) , INTENT(IN) :: inputfile
INTEGER :: ncid , rstart(2) , rend(2) ,  varid
INCLUDE 'netcdf.inc'


  !OPEN NC FILE
  CALL open_wrf_file(inputfile,'rw',ncid)
  
  rstart = (/ 1,1 /)
  rend   = (/ 19,1/)

  CALL check_io(NF_INQ_VARID(ncid,'Times',varid))
  CALL check_io(NF_PUT_VARA_TEXT(ncid,varid,rstart,rend,DATE))


  !Modify input date 
  CALL check_io(NF_REDEF(ncid))
  CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'START_DATE',19,date))
  CALL check_io(NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'SIMULATION_START_DATE',19,date))

  !CLOSE FILES
  CALL close_wrf_file(ncid)

RETURN

END SUBROUTINE put_date



END MODULE common_change_wrf_longitudes
