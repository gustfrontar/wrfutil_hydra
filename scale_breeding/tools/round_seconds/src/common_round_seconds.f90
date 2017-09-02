MODULE common_round_seconds
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

CONTAINS


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



END MODULE common_round_seconds
