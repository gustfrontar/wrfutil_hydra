MODULE COMMON_RADAR_TOOLS
!=======================================================================
!
! [PURPOSE:] Common radar procedures and variables.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!
!=======================================================================
!$USE OMP_LIB
  USE common
  IMPLICIT NONE
  PUBLIC
!-----------------------------------------------------------------------
! General parameters
!-----------------------------------------------------------------------
! Define the radar type

  REAL(r_size), PARAMETER :: minz = 0.01d0 !Minimum radar power.


  !Binary format input/output
  TYPE RADAR
  INTEGER :: na,nr,ne,nt
  INTEGER :: radar_type  !Code to identify the radar and the procedures

  REAL(r_size), ALLOCATABLE :: azimuth(:),rrange(:),elevation(:),time(:)
  CHARACTER(23)             :: start_time , end_time
  REAL(r_size)              :: res_azimuth , res_range !Resolutions in azimuth and range.
  REAL(r_size), ALLOCATABLE :: local_elevation(:,:)   !Corrected elevation (function of range and elevation)
  REAL(r_size), ALLOCATABLE :: beam_wid_h(:) , beam_wid_v(:) !Beam width in horizontal and vertical direction (degrees)
  REAL(r_size), ALLOCATABLE :: radarv3d(:,:,:)        !Observations
  REAL(r_size), ALLOCATABLE :: radarv3d_model(:,:,:)  !Model derived quantities in radar obs space
  REAL(r_size), ALLOCATABLE :: oerror(:,:,:)
  INTEGER     , ALLOCATABLE :: qcflag(:,:,:)
  REAL(r_size), ALLOCATABLE :: attenuation(:,:)
  REAL(r_size)              :: syear,smonth,sday,shour,sminute,ssecond
  REAL(r_size)              :: eyear,emonth,eday,ehour,eminute,esecond
  REAL(r_size)  :: missing
  REAL(r_size), ALLOCATABLE :: lat(:,:) ! Grid latitude
  REAL(r_size), ALLOCATABLE :: lon(:,:) ! Grid longitude
  REAL(r_size), ALLOCATABLE :: distance_to_radar(:,:)
  REAL(r_size), ALLOCATABLE :: z(:,:)   ! Heigth over sea level.
  REAL(r_size), ALLOCATABLE :: i(:,:) , j(:,:) , k(:,:)
  REAL(r_size)  :: lon0 , lat0   !Radara location
  REAL(r_size)  :: z0            !Radar height
  REAL(r_size)  :: frequency
  REAL(r_size), ALLOCATABLE :: sounding_t(:)
  REAL(r_size), ALLOCATABLE :: sounding_q(:)
  REAL(r_size), ALLOCATABLE :: sounding_z(:)
  INTEGER                   :: nsounding   !number of sounding levels.
  INTEGER  ::   nv3d     = 3
  INTEGER  ::   iv3d_ref = 1
  INTEGER  ::   iv3d_v   = 2
  INTEGER  ::   iv3d_vda = 3

  CHARACTER(20) :: element(3)=(/'dBZ','V','Vda'/)

 
  END TYPE 


CONTAINS

!-----------------------------------------------------------------------
! Set the parameters only reads radar information
!-----------------------------------------------------------------------
SUBROUTINE radar_set_common( myradar , myfile  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT)       :: myradar
  CHARACTER(LEN=*) , INTENT(IN)    :: input_file
  INTEGER      :: iunit
  INTEGER      :: ncid , varid
  INTEGER      :: varid,dimids,dimsize
  INCLUDE 'netcdf.inc'

  CALL check_io(NF90_OPEN(myfile,NF_NOWRITE,ncid))

  CALL get_cfradial_dates(ncid,myradar)

  !Radar loction
  CALL check_io(NF90_INQ_VARID(ncid,'latitude',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,myradar%lat0))
  CALL check_io(NF90_INQ_VARID(ncid,'longitude',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,myradar%lon0))
  CALL check_io(NF90_INQ_VARID(ncid,'altitude',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,myradar%z0))

  !Get beam dimensions
  CALL check_io(NF90_INQ_VARID(ncid,'radar_beam_width_h',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,myradar%beam_wid_h))
  CALL check_io(NF90_INQ_VARID(ncid,'radar_beam_width_v',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,myradar%beam_wid_v))

  !Get radar frequency
  CALL check_io(NF90_INQ_VARID(ncid,'frequency',varid))
  CALL check_io(NF90_GET_VAR_REAL(ncid,varid,frequency))

  !Get volume dimensions
  CALL check_io(NF90_INQ_VARID(ncid,'azimuth',varid))
  CALL check_io(NF90_INQ_VARDIMID(ncid,varid,dimid))
  CALL check_io(NF90_INQ_DIMLEN(ncid,dimid,dimsize))
  myradar%na = dimsize
  CALL check_io(NF90_INQ_VARID(ncid,'range',varid))
  CALL check_io(NF90_INQ_VARDIMID(ncid,varid,dimid))
  CALL check_io(NF90_INQ_DIMLEN(ncid,dimid,dimsize))
  myradar%nr = dimsize
  CALL check_io(NF90_INQ_VARID(ncid,'sweep_number',varid))
  CALL check_io(NF90_INQ_VARDIMID(ncid,varid,dimid))
  CALL check_io(NF90_INQ_DIMLEN(ncid,dimid,dimsize))
  myradar%ne = dimsize
  CALL check_io(NF90_INQ_VARID(ncid,'time',varid))
  CALL check_io(NF90_INQ_VARDIMID(ncid,varid,dimid))
  CALL check_io(NF90_INQ_DIMLEN(ncid,dimid,dimsize))
  myradar%nt = dimsize
 
  !Get grid coordinates
  ALLOCATE( azimuth(na),rrange(nr),elevation(nt),time(nt) )

  ALLOCATE(start(2),count(2),dnw(nlev),znu(nlev))

  CALL check_io(NF90_INQ_VARID(ncid,'azimuth',varid))
  CALL check_io(NF90_GET_VARA_REAL(ncid,varid,1,na,azimuth))
  CALL check_io(NF90_INQ_VARID(ncid,'range',varid))
  CALL check_io(NF90_GET_VARA_REAL(ncid,varid,1,nr,rrange))
  CALL check_io(NF90_INQ_VARID(ncid,'elevation',varid))
  CALL check_io(NF90_GET_VARA_REAL(ncid,varid,1,nt,elevation))
  CALL check_io(NF90_INQ_VARID(ncid,'time',varid))
  CALL check_io(NF90_GET_VARA_REAL(ncid,varid,1,nt,time))

  WRITE(6,'(a)') '============================================='
  WRITE(6,'(a)') '               RADAR PARAMETERS'
  WRITE(6,'(a)') ' -------------------------------------------'
  WRITE(6,'(a, 3f12.3)') '  location     :', myradar%lon0, myradar%lat0 , myradar%z0
  WRITE(6,'(a, 3i)') '      grid size    :', myradar%na ,  myradar%nr , myradar%ne

 
  RETURN
END SUBROUTINE radar_set_common


!-----------------------------------------------------------------------
! Compute the radar grid in LAT , LON , Z
!-----------------------------------------------------------------------
SUBROUTINE radar_georeference( myradar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT) :: myradar
REAL(r_size)  :: ke 
REAL(r_size)  :: tmp1 , tmp2 
INTEGER       :: it , ir 

ke=(4d0/3d0)

ALLOCATE( input_radar%z( myradar%nt , myradar%nr  ) )
ALLOCATE( input_radar%lat( myradar%nt, myradar%nr ) )
ALLOCATE( input_radar%lon( myradar%nt, myradar%nr ) )
ALLOCATE( input_radar%local_elevation( myradar%nt , myradar%nr ) )
ALLOCATE( input_radar%distance_to_radar( myradar%nt , myradar%nr ) )

!Two possibilityes:
!1) if we have a vertical profile of T and Q we can use them to compute elevation
!for each radar beam.
!2) if we do not have information about the vertical profile of T or Q then
!proced with the standard computation.

  IF ( ALLOCATED( myradar%sounding_t ) .AND. ALLOCATED( myradar%sounding_q ) .AND. ALLOCATED( myradar%sounding_z) ) THEN
  !Perform the accurate height estimation based on the input sounding.
  

  !TODO Code this option.


  ELSE

  !Perform standard height beam heigth computation.
   DO ir=1,myradar%nr
     DO it=1,myradar%nt

      tmp1= myradar%rrange(ir)**2 + (ke*Re)**2 + 2*myradar%rrange(ir)*ke*Re*sin( myradar%elevation(it)*deg2rad )
      myradar%z(it,ir)=myradar%z0 + SQRT( tmp1 ) - ke*Re
      !Compute the local elevation angle tacking into account the efective earth radius.
      tmp1= myradar%rrange(ir) * cos( myradar%elevation(it)*deg2rad )
      tmp2= myradar%rrange(ir) * sin( myradar%elevation(it)*deg2rad ) + ke*Re
      myradar%local_elevation(it,ir) = myradar%elevation(it) + atan(tmp1/tmp2)
      
     ENDDO
   ENDDO

  ENDIF

  !Compute the latitude and longitude corresponding to each radar point.
   DO ir=1,myradar%nr
    DO it=1,myradar%nt
       myradar%distance_to_radar(it,ir)=ke*Re*asin( myradar%rrange(ir) * cos( myradar%elevation(it)*deg2rad) / ( myradar%z(it,ir) + ke*Re) )
       
       CALL com_ll_arc_distance(myradar%lon0,myradar%lat0,                                          &
             myradar%distance_to_radar(it,ir),myradar%azimuth(it),myradar%lon(it,ir), &
             myradar%lat(it,ir))
    ENDDO
   ENDDO


END SUBROUTINE radar_georeference



SUBROUTINE get_cfradial_dates(ncid,myradar)
IMPLICIT NONE
INTEGER, INTENT(IN)         :: ncid
TYPE(RADAR) , intent(INOUT) :: myradar

INCLUDE 'netcdf.inc'

  !OPEN NC FILE
  CALL open_wrf_file(inputfile,'ro',ncid)

  !READ START AND END DATES
  CALL check_io(NF90_GET_ATT(ncid, NF90_GLOBAL,'start_time',myradar%start_time))
  CALL check_io(NF90_GET_ATT(ncid, NF90_GLOBAL,'end_time',myradar%end_time))

  !DECOMPOSE DATES
  myradar%syear   = read(myradar%start_time(1:4))
  myradar%smonth  = read(myradar%start_time(6:7))
  myradar%sday    = read(myradar%start_time(9:10))
  myradar%shour   = read(myradar%start_time(12:13))
  myradar%sminute = read(myradar%start_time(15:16))
  myradar%ssecond = read(myradar%start_time(18:23))
  
  myradar%eyear   = read(myradar%start_time(1:4))
  myradar%emonth  = read(myradar%start_time(6:7))
  myradar%eday    = read(myradar%start_time(9:10))
  myradar%ehour   = read(myradar%start_time(12:13))
  myradar%eminute = read(myradar%start_time(15:16))
  myradar%esecond = read(myradar%start_time(18:23))
  
  !CLOSE FILES
  CALL close_wrf_file(ncid)

RETURN

END SUBROUTINE get_cfradial_dates


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


END MODULE COMMON_RADAR_TOOLS
