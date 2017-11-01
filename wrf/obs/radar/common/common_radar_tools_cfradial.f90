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
  USE netcdf
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
  CHARACTER(23)             :: start_time 
  REAL(r_size)              :: res_azimuth , res_range !Resolutions in azimuth and range.
  REAL(r_size), ALLOCATABLE :: local_elevation(:,:)   !Corrected elevation (function of range and elevation)
  REAL(r_size)              :: beam_wid_h , beam_wid_v !Beam width in horizontal and vertical direction (degrees)

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
  REAL(r_size)  :: frequency , lambda
  REAL(r_size), ALLOCATABLE :: sounding_t(:)
  REAL(r_size), ALLOCATABLE :: sounding_q(:)
  REAL(r_size), ALLOCATABLE :: sounding_z(:)
  INTEGER                   :: nsounding   !number of sounding levels.
  INTEGER  ::   nv3d     = 3
  INTEGER  ::   iv3d_ref = 1
  INTEGER  ::   iv3d_rv  = 2
  INTEGER  ::   iv3d_rvda = 3

  INTEGER  ::   nv3d_model = 2
  INTEGER  ::   iv3d_ref_model = 1
  INTEGER  ::   iv3d_rv_model  = 2

  CHARACTER(20) :: element(3)=(/'dBZ','V','Vda'/)

  CHARACTER(20) :: element_model(2)=(/'dBZ_model','V_model'/)
 
  END TYPE 


CONTAINS

!-----------------------------------------------------------------------
! Set the parameters only reads radar information
!-----------------------------------------------------------------------
SUBROUTINE radar_set_common( myradar , myfile  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT)       :: myradar
  CHARACTER(LEN=*) , INTENT(IN)    :: myfile 
  INTEGER      :: ncid,varid,dimid,dimsize
  CHARACTER(LEN=100) :: dummychar
!  INCLUDE 'netcdf.inc'

  CALL check_io(NF90_OPEN(myfile,NF90_NOWRITE,ncid))

  CALL get_cfradial_dates(ncid,myradar)

  !Radar loction
  CALL check_io(NF90_INQ_VARID(ncid,'latitude',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%lat0))
  CALL check_io(NF90_INQ_VARID(ncid,'longitude',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%lon0))
  CALL check_io(NF90_INQ_VARID(ncid,'altitude',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%z0))

  !Get beam dimensions
  CALL check_io(NF90_INQ_VARID(ncid,'radar_beam_width_h',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%beam_wid_h))
  CALL check_io(NF90_INQ_VARID(ncid,'radar_beam_width_v',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%beam_wid_v))

  !Get radar frequency
  CALL check_io(NF90_INQ_VARID(ncid,'frequency',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%frequency))

  myradar%lambda=100*clight*2*pi/myradar%frequency  !Assuming frequency is in seconds and lambda is in cm.

  !Get volume dimensions
  CALL check_io(NF90_INQ_DIMID(ncid,'range',dimid))
  CALL check_io(NF90_INQUIRE_DIMENSION(ncid,dimid,dummychar,myradar%nr))
  CALL check_io(NF90_INQ_DIMID(ncid,'sweep',dimid))
  CALL check_io(NF90_INQUIRE_DIMENSION(ncid,dimid,dummychar,myradar%ne))
  CALL check_io(NF90_INQ_DIMID(ncid,'time',dimid))
  CALL check_io(NF90_INQUIRE_DIMENSION(ncid,dimid,dummychar,myradar%nt))

  !The number of azimuths is the same as the number of times. 
  myradar%na=myradar%nt

  !Get grid coordinates
  ALLOCATE( myradar%azimuth(myradar%na),myradar%rrange(myradar%nr),myradar%elevation(myradar%nt),myradar%time(myradar%nt) )

  CALL check_io(NF90_INQ_VARID(ncid,'azimuth',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%azimuth))
  CALL check_io(NF90_INQ_VARID(ncid,'range',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%rrange))
  CALL check_io(NF90_INQ_VARID(ncid,'elevation',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%elevation))
  CALL check_io(NF90_INQ_VARID(ncid,'time',varid))
  CALL check_io(NF90_GET_VAR(ncid,varid,myradar%time))

  WRITE(6,'(a)') '============================================='
  WRITE(6,'(a)') '               RADAR PARAMETERS'
  WRITE(6,'(a)') ' -------------------------------------------'
  WRITE(6,'(a, 3f12.3)') '  location     :', myradar%lon0, myradar%lat0 , myradar%z0
  WRITE(6,'(a, 3i)') '      grid size    :', myradar%na ,  myradar%nr , myradar%ne
  WRITE(6,*)' frequency: ', myradar%frequency
  WRITE(6,*)' lambda:    ', myradar%lambda

  CALL check_io(NF90_CLOSE(ncid))
 
  RETURN
END SUBROUTINE radar_set_common

!-----------------------------------------------------------------------
! Write model data.
!-----------------------------------------------------------------------
SUBROUTINE radar_write_model( myradar , myfile  )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(IN)          :: myradar
  CHARACTER(LEN=*) , INTENT(IN)    :: myfile
  INTEGER      :: ncid , varid , ncstatus , iv
  INTEGER      :: dimid,dimsize
  INTEGER      :: timedimid , rangedimid 
!  INCLUDE 'netcdf.inc'

  CALL check_io(NF90_OPEN(myfile,NF90_WRITE,ncid))

  !Get volume dimensions
  CALL check_io(NF90_INQ_DIMID(ncid,'range',rangedimid))
  CALL check_io(NF90_INQ_DIMID(ncid,'time',timedimid))

  CALL check_io(NF90_REDEF(ncid))

 DO iv=1,myradar%nv3d_model

  NCSTATUS = NF90_INQ_VARID(ncid,myradar%element_model(iv),varid)

  IF( NCSTATUS /= NF90_NOERR) THEN !Variable is not present, then I have to create it.

    !Create the variable
    CALL check_io(NF90_DEF_VAR(ncid,myradar%element_model(iv),NF90_FLOAT, (/rangedimid,timedimid/),varid))

    !Add attributes
    IF( iv == myradar%iv3d_ref_model )THEN
   
      CALL check_io(NF90_PUT_ATT(ncid,varid,'long_name','reflectivity_from_horizontal_polarization'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'standard_name','equivalent_reflectivity_factor'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'units','dBZ'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'sampling_ratio',1.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'_FillValue',undefs))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'scale_factor',1.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'add_offset',0.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'grid_mapping','grid_mapping'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'coordinates','time range'))

    ELSEIF( iv == myradar%iv3d_rv_model )THEN

      CALL check_io(NF90_PUT_ATT(ncid,varid,'long_name','radial_velocity_from_horizontal_polarization'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'standard_name','radial_velocity_of_scatterers_away_from_instrument'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'units','m/s'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'sampling_ratio',1.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'_FillValue',undefs))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'scale_factor',1.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'add_offset',0.0e0))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'grid_mapping','grid_mapping'))
      CALL check_io(NF90_PUT_ATT(ncid,varid,'coordinates','time range'))

    ENDIF

  ENDIF
 
 ENDDO

 !End define mode
 CALL check_io(NF90_ENDDEF(ncid))
 
 !Write the data
 DO iv=1,myradar%nv3d_model

   IF ( iv == myradar%iv3d_ref_model .or.  iv == myradar%iv3d_rv_model )THEN

     write(*,*)"Writing ",myradar%element_model(iv)

     CALL check_io(NF90_INQ_VARID(ncid,myradar%element_model(iv),varid))

     CALL check_io(NF90_PUT_VAR(ncid,varid,myradar%radarv3d_model(:,:,iv),(/1,1/),(/myradar%nt , myradar%nr/)))

   ENDIF

 ENDDO


  CALL check_io(NF90_CLOSE(ncid))

END SUBROUTINE  radar_write_model


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

ALLOCATE( myradar%z( myradar%nt , myradar%nr  ) )
ALLOCATE( myradar%lat( myradar%nt, myradar%nr ) )
ALLOCATE( myradar%lon( myradar%nt, myradar%nr ) )
ALLOCATE( myradar%local_elevation( myradar%nt , myradar%nr ) )
ALLOCATE( myradar%distance_to_radar( myradar%nt , myradar%nr ) )

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
INTEGER                     :: varid
TYPE(RADAR) , intent(INOUT) :: myradar

!INCLUDE 'netcdf.inc'

  !READ START TIME
  CALL check_io(NF90_INQ_VARID(ncid,'time',varid))
  CALL check_io(NF90_GET_ATT(ncid,varid,'units',myradar%start_time)) 
  myradar%start_time=myradar%start_time(15:34)
  
  !CALL check_io(NF90_GET_ATT(ncid, NF90_GLOBAL,'start_time',myradar%start_time))
  !CALL check_io(NF90_GET_ATT(ncid, NF90_GLOBAL,'end_time',myradar%end_time))

          !DECOMPOSE DATES
  write(*,*)myradar%start_time

  read(myradar%start_time(1:4),*) myradar%syear
  read(myradar%start_time(6:7),*) myradar%smonth
  read(myradar%start_time(9:10),*) myradar%sday
  read(myradar%start_time(12:13),*) myradar%shour
  read(myradar%start_time(15:16),*) myradar%sminute
  read(myradar%start_time(18:19),*) myradar%ssecond

RETURN

END SUBROUTINE get_cfradial_dates


SUBROUTINE check_io(status)
  IMPLICIT NONE
!  INCLUDE 'netcdf.inc'
  INTEGER(4),INTENT(IN) :: status

  IF(status /= NF90_NOERR) THEN
    WRITE(6,*) TRIM(NF90_STRERROR(status))
    STOP 10
  ENDIF

  RETURN
END SUBROUTINE check_io


END MODULE COMMON_RADAR_TOOLS
