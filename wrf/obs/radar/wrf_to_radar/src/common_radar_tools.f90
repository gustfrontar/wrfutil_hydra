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

  TYPE RADAR
  INTEGER :: na,nr,ne,nv3d
  INTEGER :: radar_type  !Code to identify the radar and the procedures
                         !that will be used (some of them radar dependent)
                         ! 1 - PAWR
                         ! 2 - LIDAR
                         ! 3 - Add more radars here
  REAL(r_size), ALLOCATABLE :: azimuth(:),rrange(:),elevation(:)
  REAL(r_size), ALLOCATABLE :: local_elevation(:,:)   !Corrected elevation (function of range and elevation)
  REAL(r_size)              :: range_resolution
  REAL(r_size), ALLOCATABLE :: beam_wid_h(:) , beam_wid_v(:) !Beam width in horizontal and vertical direction (degrees)
  REAL(r_size), ALLOCATABLE :: radarv3d(:,:,:,:)        !Observations
  REAL(r_size), ALLOCATABLE :: radarv3d_model(:,:,:,:)  !Model derived quantities in radar obs space
  REAL(r_size), ALLOCATABLE :: oerror(:,:,:,:)
  REAL(r_size), ALLOCATABLE :: qcflag(:,:,:)
  REAL(r_size), ALLOCATABLE :: attenuation(:,:,:)
  REAL(r_size)              :: year,month,day,hour,minute,second
  REAL(r_size)  :: missing
  REAL(r_size), ALLOCATABLE :: lat(:,:,:) ! Grid latitude
  REAL(r_size), ALLOCATABLE :: lon(:,:,:) ! Grid longitude
  REAL(r_size), ALLOCATABLE :: z(:,:,:)   ! Heigth over sea level.
  REAL(r_size), ALLOCATABLE :: p(:,:,:)   ! Pressure (for vertical localization)
  REAL(r_size), ALLOCATABLE :: i(:,:,:) , j(:,:,:) , k(:,:,:)
  REAL(r_size)  :: lon0 , lat0   !Radara location
  REAL(r_size)  :: z0            !Radar height
  REAL(r_size)  :: lambda
  REAL(r_size), ALLOCATABLE :: sounding_t(:)
  REAL(r_size), ALLOCATABLE :: sounding_q(:)
  REAL(r_size), ALLOCATABLE :: sounding_z(:)
  INTEGER                   :: nsounding   !number of sounding levels.
  INTEGER  :: iv3d_wind , iv3d_ref , iv3d_prh , iv3d_windwidth , iv3d_snr 
  REAL(r_size)        :: total_attenuation_factor
 
  END TYPE 

  !This value will be used in case of very low reflectivityes. This will be the minimum 
  !reflectivity that will be used either in the model or in the radar.
  REAL(r_size)      :: minz=0.01d0 
  !These two flags indicates if the "offline"  QCFLAG and the attenuation of 
  !the radar data will be used in the quality control of the data.
  LOGICAL           :: USE_QCFLAG      =.FALSE.
  LOGICAL           :: USE_ATTENUATION =.FALSE.

  !IF USE_ATTENUATION == TRUE, then gates with estimated attenuation 
  !greather than the threshold will be rejected. (this does not affect
  !the computation of attenuation in the forward operator)
  REAL(r_size)      :: ATTENUATION_THRESHOLD=0.01 !0.01 means -10 dBz

  !Constants for the computation of observation error.
  REAL(r_size)      :: CONSTANT_REFLECTIVITY_ERROR=3 !Reflectivity error in dBz
  REAL(r_size)      :: CONSTANT_VR_ERROR=1           !Radial wind error in m/s 
  REAL(r_size)      :: CONSTANT_PSEUDORH_ERROR=0.2   !Pseudo relative humidity error. 

  !Threshold for pseudo hr conversion. If reflectivity is larger than this threshold
  !then a pseudo relative humidity of 100% will be assumed.
  REAL(r_size)      :: PRH_Z_THRESHOLD=0.0d0

CONTAINS
!-----------------------------------------------------------------------
! Set the parameters only reads radar information
!-----------------------------------------------------------------------
SUBROUTINE radar_set_common( input_radar , input_file )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT)       :: input_radar
  CHARACTER(LEN=*) , INTENT(IN)    :: input_file
  INTEGER      :: iunit
  INTEGER      :: ia , ir , ie
  REAL(r_sngl) ,ALLOCATABLE :: buf4(:)

  !INITIALIZE VARIABLES ID  (-9 means variable not avaliable in this radar)
  input_radar%iv3d_ref=-9
  input_radar%iv3d_wind=-9
  input_radar%iv3d_prh=-9
  input_radar%iv3d_windwidth=-9
  input_radar%iv3d_snr=-9


  IF ( input_radar % radar_type .EQ. 1 )THEN
  iunit=99

  !Define data order for the PAWR data 
  input_radar%iv3d_ref=1    !Reflectivity
  input_radar%iv3d_wind=2   !Radial velocity

  !Get radar dimension from the input file.
  !Skip data only read array dimensions and grid information.
  OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential')
  
  !READ TIME HEADER
  ALLOCATE( buf4( 6 ) )
  READ(iunit)buf4
  input_radar%year  =REAL(buf4(1),r_size)
  input_radar%month =REAL(buf4(2),r_size)
  input_radar%day   =REAL(buf4(3),r_size)
  input_radar%hour  =REAL(buf4(4),r_size)
  input_radar%minute=REAL(buf4(5),r_size)
  input_radar%second=REAL(buf4(6),r_size)
  DEALLOCATE( buf4)

  !READ RADAR LOCATION AND CHARACTERISTICS HEADER
  ALLOCATE( buf4 ( 8 ) )
  READ(iunit) buf4
  input_radar%lon0            =REAL(buf4(1),r_size)
  input_radar%lat0            =REAL(buf4(2),r_size)
  input_radar%z0              =REAL(buf4(3),r_size)
  ALLOCATE( input_radar%beam_wid_h( input_radar%na ) )
  ALLOCATE( input_radar%beam_wid_v( input_radar%ne ) )
  input_radar%beam_wid_h      =REAL(buf4(4),r_size) !In PAWR this is a constant.
  input_radar%beam_wid_v      =REAL(buf4(5),r_size) !In PAWR this is a constant.
  input_radar%range_resolution=REAL(buf4(6),r_size)
  input_radar%lambda          =REAL(buf4(7),r_size)
  input_radar%missing         =REAL(buf4(8),r_size)
  DEALLOCATE( buf4 )

  !READ DATA DIMESION HEADER
  READ(iunit)input_radar%na , input_radar%nr , input_radar%ne , input_radar%nv3d

  !READ DIMENSIONS DATA
  ALLOCATE( input_radar % azimuth  ( input_radar % na ) )
  ALLOCATE( input_radar % rrange   ( input_radar % nr ) )
  ALLOCATE( input_radar % elevation( input_radar % ne ) )

  ALLOCATE( buf4( input_radar % na ) ) 
  READ(iunit)buf4
  input_radar%azimuth=REAL(buf4,r_size)
  DEALLOCATE(buf4)
  ALLOCATE(buf4 ( input_radar % nr ) )
  READ(iunit)buf4
  input_radar%rrange=REAL(buf4,r_size)
  DEALLOCATE(buf4)
  ALLOCATE(buf4 ( input_radar % ne ) )
  READ(iunit)buf4
  input_radar%elevation=REAL(buf4,r_size)
  DEALLOCATE(buf4)

  !READ TOTAL ATTENUATION FACTOR
  ALLOCATE(buf4(1))
  READ(iunit)buf4
  input_radar%total_attenuation_factor=REAL(buf4(1),r_size)
  DEALLOCATE(buf4)


  WRITE(6,'(a)') '============================================='
  WRITE(6,'(a)') '               RADAR PARAMETERS'
  WRITE(6,'(a)') ' -------------------------------------------'
  WRITE(6,'(a, 3f12.3)') '  location     :', input_radar%lon0, input_radar%lat0 ,input_radar%z0
  WRITE(6,'(a, 3i)') '  grid size :', input_radar%na , input_radar%nr , input_radar%ne

  CLOSE(iunit) 
 
  ELSEIF( input_radar % radar_type .EQ. 2)THEN
  !LIDAR DATA
    iunit=99

    !Define data order for the LIDAR data.
    input_radar%iv3d_snr=1    !signal to noise ratio.
    input_radar%iv3d_wind=2   !Radial velocity
    input_radar%iv3d_windwidth=3 !Velocity spectral width.

    !Get radar dimension from the input file.
    !Skip data only read array dimensions and grid information.
    OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential')
  
    !READ TIME HEADER
    ALLOCATE( buf4( 6 ) )
    READ(iunit)buf4
    input_radar%year  =REAL(buf4(1),r_size)
    input_radar%month =REAL(buf4(2),r_size)
    input_radar%day   =REAL(buf4(3),r_size)
    input_radar%hour  =REAL(buf4(4),r_size)
    input_radar%minute=REAL(buf4(5),r_size)
    input_radar%second=REAL(buf4(6),r_size)
    DEALLOCATE( buf4)

    !READ RADAR LOCATION AND CHARACTERISTICS HEADER
    ALLOCATE( buf4 ( 6 ) )
    READ(iunit) buf4
    input_radar%lon0            =REAL(buf4(1),r_size)
    input_radar%lat0            =REAL(buf4(2),r_size)
    input_radar%z0              =REAL(buf4(3),r_size)
    input_radar%range_resolution=REAL(buf4(4),r_size)
    input_radar%lambda          =REAL(buf4(5),r_size)
    input_radar%missing         =REAL(buf4(6),r_size)
    DEALLOCATE( buf4 )
  
    !READ DATA DIMESION HEADER
    READ(iunit)input_radar%na , input_radar%nr , input_radar%ne , input_radar%nv3d
 
    !READ DIMENSIONS DATA
    ALLOCATE( input_radar % azimuth  ( input_radar % na ) )
    ALLOCATE( input_radar % rrange   ( input_radar % nr ) )
    ALLOCATE( input_radar % elevation( input_radar % ne ) )
    ALLOCATE( input_radar%beam_wid_h ( input_radar%na ) )
    ALLOCATE( input_radar%beam_wid_v ( input_radar%ne ) )

    ALLOCATE( buf4( input_radar % na ) )
    READ(iunit)buf4
    input_radar%azimuth=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % nr ) )
    READ(iunit)buf4
    input_radar%rrange=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % ne ) )
    READ(iunit)buf4
    input_radar%elevation=REAL(buf4,r_size)
    DEALLOCATE(buf4)

    ALLOCATE(buf4 ( input_radar % na ) )
    READ(iunit)buf4
    input_radar%beam_wid_h=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % ne ) )
    READ(iunit)buf4
    input_radar%beam_wid_v=REAL(buf4,r_size)
    DEALLOCATE(buf4)

    WRITE(6,'(a)') '============================================='
    WRITE(6,'(a)') '               LIDAR PARAMETERS'
    WRITE(6,'(a)') ' -------------------------------------------'
    WRITE(6,'(a, 3f12.3)') '  location     :', input_radar%lon0, input_radar%lat0 ,input_radar%z0
    WRITE(6,'(a, 3i)') '  grid size :', input_radar%na , input_radar%nr , input_radar%ne

    CLOSE(iunit)

  ELSE
   WRITE(6,*)'ERROR: Not recognized radar type in set_common_radar'
   STOP

  ENDIF !END OF RADAR DATA INITIALIZATION

  RETURN
END SUBROUTINE radar_set_common

!-----------------------------------------------------------------------
! Write a radar file
!-----------------------------------------------------------------------
SUBROUTINE radar_write_file( input_radar , input_data , nvar , output_file )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT) :: input_radar
  INTEGER,     INTENT(IN)    :: nvar
  CHARACTER*20, INTENT(IN)   :: output_file
  REAL(r_size), INTENT(IN)   :: input_data(input_radar%na,input_radar%nr,input_radar%ne,nvar)
  INTEGER      :: iunit , ivar
  INTEGER      :: ia , ir , ie

  IF ( input_radar % radar_type .EQ. 1 )THEN
    iunit=99
    !Write data in PAWR data format.
    OPEN(iunit,FILE=output_file,FORM='unformatted',ACCESS='sequential')
  
    !WRITE TIME HEADER
    WRITE(iunit)REAL(input_radar % year   , r_sngl), &
                REAL(input_radar % month  , r_sngl), &
                REAL(input_radar % day    , r_sngl), &
                REAL(input_radar % hour   , r_sngl), &
                REAL(input_radar % minute , r_sngl), &
                REAL(input_radar % second , r_sngl)
  
    !WRITE RADAR CHARACTERISTICS HEADER
    WRITE(iunit)  REAL( input_radar%lon0        , r_sngl), &
                  REAL( input_radar%lat0        , r_sngl), &
                  REAL( input_radar%z0          , r_sngl), &
                  REAL( input_radar%beam_wid_h(1), r_sngl), &
                  REAL( input_radar%beam_wid_v(1), r_sngl), &
                  REAL( input_radar%range_resolution, r_sngl), &
                  REAL( input_radar%lambda      , r_sngl), &
                  REAL( input_radar%missing , r_sngl)

    !READ DATA DIMESION HEADER
    WRITE(iunit)input_radar%na , input_radar%nr , input_radar%ne , nvar
  
  
    !READ DIMENSIONS DATA
    WRITE(iunit)REAL( input_radar%azimuth  , r_sngl )
    WRITE(iunit)REAL( input_radar%rrange   , r_sngl )
    WRITE(iunit)REAL( input_radar%elevation, r_sngl )

    !READ TOTAL ATTENUATION FACTOR
    WRITE(iunit)REAL( input_radar%total_attenuation_factor , r_sngl)

  
    DO ivar=1,nvar
      DO ie=1,input_radar % ne
        WRITE(iunit) REAL( input_data(:,:,ie,ivar) , r_sngl )
      ENDDO
    ENDDO

    !QC FLAG
    DO ie=1,input_radar % ne
      WRITE(iunit) REAL( input_radar%qcflag(:,:,ie) , r_sngl )
    ENDDO
    !ATTENUATION
    DO ie=1,input_radar % ne
       WRITE(iunit) REAL( input_radar%attenuation(:,:,ie), r_sngl )
    ENDDO


  
    CLOSE(iunit)

  ELSEIF ( input_radar % radar_type .EQ. 2 )THEN
  !THIS IS LIDAR DATA FORMAT

    iunit=99
    !Write data in LIDAR data format.
    OPEN(iunit,FILE=output_file,FORM='unformatted',ACCESS='sequential')

    !WRITE TIME HEADER
    WRITE(iunit)REAL(input_radar % year   , r_sngl), &
                REAL(input_radar % month  , r_sngl), &
                REAL(input_radar % day    , r_sngl), &
                REAL(input_radar % hour   , r_sngl), &
                REAL(input_radar % minute , r_sngl), &
                REAL(input_radar % second , r_sngl)

    !WRITE RADAR CHARACTERISTICS HEADER
    WRITE(iunit)  REAL( input_radar%lon0        , r_sngl), &
                  REAL( input_radar%lat0        , r_sngl), &
                  REAL( input_radar%z0          , r_sngl), &
                  REAL( input_radar%range_resolution, r_sngl), &
                  REAL( input_radar%lambda      , r_sngl), &
                  REAL( input_radar%missing     , r_sngl) 

    !READ DATA DIMESION HEADER
    WRITE(iunit)input_radar%na , input_radar%nr , input_radar%ne , nvar

    !READ TOTAL ATTENUATION FACTOR
    WRITE(iunit)REAL( input_radar%total_attenuation_factor , r_sngl)

    !READ DIMENSIONS DATA
    WRITE(iunit)REAL( input_radar%azimuth  , r_sngl )
    WRITE(iunit)REAL( input_radar%rrange   , r_sngl )
    WRITE(iunit)REAL( input_radar%elevation, r_sngl )
    WRITE(iunit)REAL( input_radar%beam_wid_h, r_sngl)
    WRITE(iunit)REAL( input_radar%beam_wid_v, r_sngl)

    DO ivar=1,input_radar % nv3d
      DO ie=1,input_radar % ne
        WRITE(iunit) REAL( input_data(:,:,ie,ivar) , r_sngl )
      ENDDO
    ENDDO



    CLOSE(iunit)


  ELSE
    WRITE(6,*)'ERROR: I can not write this radar because I dont recognize the radar type'

  ENDIF !END OF FILE WRITING

  RETURN

END SUBROUTINE radar_write_file


!-----------------------------------------------------------------------
! Read radar data
! This function reads one volume of radar data.
!-----------------------------------------------------------------------
SUBROUTINE radar_read_data( input_radar , input_file )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT)       :: input_radar
CHARACTER(LEN=*) ,INTENT(IN)     :: input_file
INTEGER                    :: iunit , ivar , inttmp
INTEGER                    :: ia , ir , ie 
REAL(r_sngl),ALLOCATABLE   :: buf4(:),buf2d4(:,:)

!Reads radar information and data at the same time

iunit=98
!First select procedure according to the radar type.

 IF ( input_radar % radar_type .EQ. 1 )THEN
  OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential')

  !Define data order for the PAWR data 
  input_radar%iv3d_ref=1    !Reflectivity
  input_radar%iv3d_wind=2   !Radial velocity

  !READ TIME HEADER
  ALLOCATE( buf4( 6 ) )
  READ(iunit)buf4
  input_radar%year  =REAL(buf4(1),r_size)
  input_radar%month =REAL(buf4(2),r_size)
  input_radar%day   =REAL(buf4(3),r_size)
  input_radar%hour  =REAL(buf4(4),r_size)
  input_radar%minute=REAL(buf4(5),r_size)
  input_radar%second=REAL(buf4(6),r_size)
  DEALLOCATE( buf4)

  !READ RADAR LOCATION AND CHARACTERISTICS HEADER
  ALLOCATE( buf4 ( 8 ) )
  READ(iunit) buf4
  input_radar%lon0            =REAL(buf4(1),r_size)
  input_radar%lat0            =REAL(buf4(2),r_size)
  input_radar%z0              =REAL(buf4(3),r_size)

  ALLOCATE( input_radar%beam_wid_h( input_radar%na ) )
  ALLOCATE( input_radar%beam_wid_v( input_radar%ne ) )

  input_radar%beam_wid_h      =REAL(buf4(4),r_size)  !In PAWR this is a constant.
  input_radar%beam_wid_v      =REAL(buf4(5),r_size)  !In PAWR this is a constant.
  input_radar%range_resolution=REAL(buf4(6),r_size)
  input_radar%lambda          =REAL(buf4(7),r_size)
  input_radar%missing         =REAL(buf4(8),r_size)
  DEALLOCATE( buf4 )

  !READ DATA DIMESION HEADER
  READ(iunit)input_radar%na , input_radar%nr , input_radar%ne , input_radar%nv3d

  !READ DIMENSIONS DATA
  ALLOCATE( input_radar % azimuth  ( input_radar % na ) )
  ALLOCATE( input_radar % rrange   ( input_radar % nr ) )
  ALLOCATE( input_radar % elevation( input_radar % ne ) )

  ALLOCATE( buf4( input_radar % na ) )
  READ(iunit)buf4
  input_radar%azimuth=REAL(buf4,r_size)
  DEALLOCATE(buf4)
  ALLOCATE(buf4 ( input_radar % nr ) )
  READ(iunit)buf4
  input_radar%rrange=REAL(buf4,r_size)
  DEALLOCATE(buf4)
  ALLOCATE(buf4 ( input_radar % ne ) )
  READ(iunit)buf4
  input_radar%elevation=REAL(buf4,r_size)
  DEALLOCATE(buf4)

  !READ TOTAL ATTENUATION FACTOR
  ALLOCATE(buf4(1))
  READ(iunit)buf4
  input_radar%total_attenuation_factor=REAL(buf4(1),r_size)
  DEALLOCATE(buf4)

  ALLOCATE(input_radar % radarv3d(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d))
  ALLOCATE(input_radar % qcflag(input_radar%na,input_radar%nr,input_radar%ne))
  ALLOCATE(input_radar % attenuation(input_radar%na,input_radar%nr,input_radar%ne)) 

  ALLOCATE( buf2d4 ( input_radar%na , input_radar%nr ) )
  DO ivar=1,input_radar % nv3d
    DO ie=1,input_radar % ne
      READ(iunit) buf2d4
      input_radar%radarv3d(:,:,ie,ivar)=REAL(buf2d4,r_size) 
    ENDDO
  ENDDO

    DO ie=1,input_radar % ne
      READ(iunit) buf2d4
      input_radar%qcflag(:,:,ie)=REAL(buf2d4,r_size)
    ENDDO

    DO ie=1,input_radar % ne
      READ(iunit) buf2d4
      input_radar%attenuation(:,:,ie)=REAL(buf2d4,r_size)
    ENDDO
  
  DEALLOCATE( buf2d4 )


  CLOSE(iunit)

 ELSEIF( input_radar % radar_type .EQ. 2)THEN
 !LIDAR DATA

  !LIDAR DATA
    iunit=99

    !Define data order for the PAWR data 
    input_radar%iv3d_snr=1    !signal to noise ratio.
    input_radar%iv3d_wind=2   !Radial velocity
    input_radar%iv3d_windwidth=3 !Velocity spectral width.

    !Get radar dimension from the input file.
    !Skip data only read array dimensions and grid information.
    OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential')

    !READ TIME HEADER
    ALLOCATE( buf4( 6 ) )
    READ(iunit)buf4
    input_radar%year  =REAL(buf4(1),r_size)
    input_radar%month =REAL(buf4(2),r_size)
    input_radar%day   =REAL(buf4(3),r_size)
    input_radar%hour  =REAL(buf4(4),r_size)
    input_radar%minute=REAL(buf4(5),r_size)
    input_radar%second=REAL(buf4(6),r_size)
    DEALLOCATE( buf4)

    !READ RADAR LOCATION AND CHARACTERISTICS HEADER
    ALLOCATE( buf4 ( 6 ) )
    READ(iunit) buf4
    input_radar%lon0            =REAL(buf4(1),r_size)
    input_radar%lat0            =REAL(buf4(2),r_size)
    input_radar%z0              =REAL(buf4(3),r_size)
    input_radar%range_resolution=REAL(buf4(4),r_size)
    input_radar%lambda          =REAL(buf4(5),r_size)
    input_radar%missing         =REAL(buf4(6),r_size)
    DEALLOCATE( buf4 )

    !READ DATA DIMESION HEADER
    READ(iunit)input_radar%na , input_radar%nr , input_radar%ne , input_radar%nv3d

    !READ DIMENSIONS DATA
    ALLOCATE( input_radar % azimuth  ( input_radar % na ) )
    ALLOCATE( input_radar % rrange   ( input_radar % nr ) )
    ALLOCATE( input_radar % elevation( input_radar % ne ) )
    ALLOCATE( input_radar%beam_wid_h ( input_radar%na ) )
    ALLOCATE( input_radar%beam_wid_v ( input_radar%ne ) )

    ALLOCATE( buf4( input_radar % na ) )
    READ(iunit)buf4
    input_radar%azimuth=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % nr ) )
    READ(iunit)buf4
    input_radar%rrange=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % ne ) )
    READ(iunit)buf4
    input_radar%elevation=REAL(buf4,r_size)
    DEALLOCATE(buf4)

    ALLOCATE(buf4 ( input_radar % na ) )
    READ(iunit)buf4
    input_radar%beam_wid_h=REAL(buf4,r_size)
    DEALLOCATE(buf4)
    ALLOCATE(buf4 ( input_radar % ne ) )
    READ(iunit)buf4
    input_radar%beam_wid_v=REAL(buf4,r_size)
    DEALLOCATE(buf4)

    ALLOCATE(input_radar % radarv3d(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d))

    ALLOCATE( buf2d4 ( input_radar%na , input_radar%nr ) )

    WRITE(*,*)input_radar%na,input_radar%nr
    DO ivar=1,input_radar % nv3d
      DO ie=1,input_radar % ne
       READ(iunit) buf2d4
       input_radar%radarv3d(:,:,ie,ivar)=REAL(buf2d4,r_size)
      ENDDO
    ENDDO

    DEALLOCATE( buf2d4 )

    CLOSE(iunit)

   ELSE

    
   WRITE(6,*)'ERROR: Not recognized radar type in radar_read_data'
   STOP

 ENDIF  ! END 

END SUBROUTINE radar_read_data

!-----------------------------------------------------------------------
! Compute radar pseudo rh (for convection initiation / supression )
!-----------------------------------------------------------------------
SUBROUTINE radar_pseudorh ( input_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT)  :: input_radar
INTEGER                     :: ia , ir , ie , iv
REAL(r_size),ALLOCATABLE    :: tmparray(:,:,:,:)
LOGICAL                     :: isalloc

!This procedure will derive psedo hr if reflectivity is present.
ISALLOC=ALLOCATED(input_radar%radarv3d)
IF( input_radar%iv3d_ref .GT. 0 .AND. ISALLOC .AND. input_radar%iv3d_prh .LT. 0 )THEN

  !Allocate the temporal array
  ALLOCATE(tmparray(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d))
  tmparray=input_radar%radarv3d !store data in the temporary array
  !Increas the number of radar variables by one
  input_radar % nv3d = input_radar % nv3d + 1
  !Increment the size of input_radar%radarv3d
  DEALLOCATE( input_radar%radarv3d)
  ALLOCATE( input_radar%radarv3d(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d))
  input_radar%radarv3d(:,:,:,1:input_radar%nv3d-1)=tmparray
  input_radar%iv3d_prh=input_radar%nv3d 

  input_radar%radarv3d(:,:,:,input_radar%iv3d_prh)=input_radar%missing

  !Compute input radar pseudo relative humidity.
  DO ia=1,input_radar % na
   DO ir=1,input_radar % nr
    DO ie=1,input_radar % ne
      IF ( input_radar%radarv3d(ia,ir,ie,input_radar%iv3d_ref) .EQ. input_radar%missing) CYCLE
        !At this point pseudo hr will be 100% where we have reflectivities
        !over 0dbz and 0 elsewhere. Of course the 0 values will be changed later
        !to something closer to the rh forecasted by the model in order to dry it
        !if necessary but not too much.
       IF ( input_radar%radarv3d(ia,ir,ie,input_radar%iv3d_ref) .GT. PRH_Z_THRESHOLD )THEN
         input_radar%radarv3d(ia,ir,ie,input_radar%iv3d_prh)= 1.0d0
       ELSE
         input_radar%radarv3d(ia,ir,ie,input_radar%iv3d_prh)=0.0d0
       ENDIF
     ENDDO
    ENDDO
   ENDDO

ELSEIF( input_radar%iv3d_ref .GT. 0 .AND. .NOT. ISALLOC .AND. input_radar%iv3d_prh .LT. 0)THEN
    !Data is not allocated, we just find the place to put this variable.
    input_radar % nv3d = input_radar % nv3d + 1
    input_radar%iv3d_prh=input_radar % nv3d
ENDIF

END SUBROUTINE radar_pseudorh

!-----------------------------------------------------------------------
! Apply radar QC data
!-----------------------------------------------------------------------
SUBROUTINE radar_apply_qc( input_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT) :: input_radar
REAL ( r_size )            :: mindbz
INTEGER                    :: ia , ir , ie , iv
LOGICAL                    :: isalloc

mindbz=10.0d0*log10(minz)

!In PAWR zero reflectivity appear with a very low value.
!Values lower to the minimum dbz value defined in this module will be set
!to the minimum dbz value. The same value will be used to limit the lower bound
!of model derived reflectivities.

IF( input_radar % iv3d_ref .GT. 0 ) THEN

  DO ia=1,input_radar%na
    DO ir=1,input_radar%nr
      DO ie=1,input_radar%ne
        IF( input_radar%radarv3d(ia,ir,ie,input_radar%iv3d_ref) .LT. mindbz )THEN
          input_radar % radarv3d(ia,ir,ie,input_radar%iv3d_ref)=mindbz
        ENDIF
      ENDDO
    ENDDO
  ENDDO

ENDIF

!QC FLAG is applied to all variables.
ISALLOC=ALLOCATED( input_radar % qcflag )
IF( ISALLOC .AND. USE_QCFLAG ) THEN

  DO ia=1,input_radar%na
    DO ir=1,input_radar%nr
      DO ie=1,input_radar%ne
        IF( input_radar % qcflag(ia,ir,ie) .GT. 900 )THEN
          input_radar % radarv3d(ia,ir,ie,:)=input_radar% missing
        ENDIF
      ENDDO
    ENDDO
  ENDDO

ENDIF

!Attenuation only affects reflectivity.
ISALLOC=ALLOCATED( input_radar % attenuation )
IF( ISALLOC .AND. USE_ATTENUATION .AND. input_radar%iv3d_ref .GT. 0)THEN
  DO ia=1,input_radar%na
    DO ir=1,input_radar%nr
      DO ie=1,input_radar%ne
        IF( input_radar % attenuation(ia,ir,ie) < ATTENUATION_THRESHOLD )THEN
          input_radar % radarv3d(ia,ir,ie,input_radar%iv3d_ref)=input_radar% missing
        ENDIF
      ENDDO
    ENDDO
  ENDDO

ENDIF


END SUBROUTINE radar_apply_qc
!-----------------------------------------------------------------------
! Compute the radar grid in LAT , LON , Z
!-----------------------------------------------------------------------
SUBROUTINE radar_georeference( input_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT) :: input_radar
REAL(r_size)  :: ke 
REAL(r_size)  :: tmp1 , tmp2 
REAL(r_size)  :: distance_to_radar
INTEGER       :: ia , ir , ie

ke=(4d0/3d0)

ALLOCATE( input_radar%z( input_radar%na , input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%lat(input_radar%na, input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%lon(input_radar%na, input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%local_elevation( input_radar%nr , input_radar%ne    ) )

!Two possibilityes:
!1) if we have a vertical profile of T and Q we can use them to compute elevation
!for each radar beam.
!2) if we do not have information about the vertical profile of T or Q then
!proced with the standard computation.

  IF ( ALLOCATED( input_radar%sounding_t ) .AND. ALLOCATED( input_radar%sounding_q ) .AND. ALLOCATED( input_radar%sounding_z) ) THEN
  !Perform the accurate height estimation based on the input sounding.
  

  !TODO Code this option.


  ELSE

  !Perform standard height beam heigth computation.
   DO ir=1,input_radar%nr
     DO ie=1,input_radar%ne

      tmp1= input_radar%rrange(ir)**2 + (ke*Re)**2 + 2*input_radar%rrange(ir)*ke*Re*sin( input_radar%elevation(ie)*deg2rad )
      input_radar%z(:,ir,ie)=input_radar%z0 + SQRT( tmp1 )-ke*Re;
      !Compute the local elevation angle tacking into account the efective earth radius.
      tmp1= input_radar%rrange(ir) * cos( input_radar%elevation(ie)*deg2rad )
      tmp2= input_radar%rrange(ir) * sin( input_radar%elevation(ie)*deg2rad ) + ke*Re
      input_radar%local_elevation(ir,ie) = input_radar%elevation(ie) + atan(tmp1/tmp2)
      
     ENDDO
   ENDDO

  ENDIF

  !Compute the latitude and longitude corresponding to each radar point.
   DO ir=1,input_radar%nr
    DO ie=1,input_radar%ne
       distance_to_radar=ke*Re*asin( input_radar%rrange(ir) * cos(input_radar%elevation(ie)*deg2rad) / (ke*Re) )
       
       DO ia=1,input_radar%na
        CALL com_ll_arc_distance(input_radar%lon0,input_radar%lat0,               &
             distance_to_radar,input_radar%azimuth(ia),input_radar%lon(ia,ir,ie), &
             input_radar%lat(ia,ir,ie))
       ENDDO
    ENDDO
   ENDDO

   !DEBUG
   !DO ir=1,input_radar%nr
   !
   !    WRITE(6,*)input_radar%lat(40,ir,1) , input_radar%lon(40,ir,1)
   !ENDDO
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE radar_georeference

SUBROUTINE radar_deallocate ( input_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT) :: input_radar
LOGICAL                    :: ISALLOC

ISALLOC = ALLOCATED(input_radar % azimuth)
IF(ISALLOC)DEALLOCATE( input_radar % azimuth )
ISALLOC = ALLOCATED(input_radar % rrange)
IF(ISALLOC)DEALLOCATE( input_radar % rrange  )
ISALLOC = ALLOCATED(input_radar % elevation) 
IF(ISALLOC)DEALLOCATE( input_radar % elevation  )
ISALLOC = ALLOCATED(input_radar % local_elevation) 
IF(ISALLOC)DEALLOCATE( input_radar % local_elevation  )
ISALLOC = ALLOCATED(input_radar % radarv3d) 
IF(ISALLOC)DEALLOCATE( input_radar % radarv3d  )
ISALLOC = ALLOCATED(input_radar % radarv3d_model) 
IF(ISALLOC)DEALLOCATE( input_radar % radarv3d_model  )
ISALLOC = ALLOCATED(input_radar % qcflag) 
IF(ISALLOC)DEALLOCATE( input_radar % qcflag  )
ISALLOC = ALLOCATED(input_radar % attenuation) 
IF(ISALLOC)DEALLOCATE( input_radar % attenuation  )
ISALLOC = ALLOCATED(input_radar % lat) 
IF(ISALLOC)DEALLOCATE( input_radar % lat  )
ISALLOC = ALLOCATED(input_radar % lon) 
IF(ISALLOC)DEALLOCATE( input_radar % lon  )
ISALLOC = ALLOCATED(input_radar % z) 
IF(ISALLOC)DEALLOCATE( input_radar % z  )
ISALLOC = ALLOCATED(input_radar % sounding_t) 
IF(ISALLOC)DEALLOCATE( input_radar % sounding_t  )
ISALLOC = ALLOCATED(input_radar % sounding_q) 
IF(ISALLOC)DEALLOCATE( input_radar % sounding_q  )
ISALLOC = ALLOCATED(input_radar % sounding_z) 
IF(ISALLOC)DEALLOCATE( input_radar % sounding_z  )
ISALLOC = ALLOCATED(input_radar % oerror )
IF(ISALLOC)DEALLOCATE( input_radar % oerror )
ISALLOC = ALLOCATED(input_radar % beam_wid_h )
IF(ISALLOC)DEALLOCATE( input_radar % beam_wid_h )
ISALLOC = ALLOCATED(input_radar % beam_wid_v )
IF(ISALLOC)DEALLOCATE( input_radar % beam_wid_v )


END SUBROUTINE  radar_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE ERROR IN RADAR OBSERVATIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE radar_observation_error(input_radar)

IMPLICIT NONE
TYPE(RADAR)  :: input_radar

!Currently a simple constant error model is implemented.
!In the future error can be computed as a function of attenuation factor or distance from
!gate to radar as many authors suggested.

ALLOCATE( input_radar%oerror(input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )

!WRITE(6,*)input_radar%iv3d_ref , input_radar%iv3d_wind , input_radar%iv3d_prh,input_radar%nv3d

IF( input_radar%iv3d_ref .GT. 0 )THEN
input_radar%oerror(:,:,:,input_radar%iv3d_ref)=CONSTANT_REFLECTIVITY_ERROR
ENDIF
IF( input_radar%iv3d_wind .GT. 0 )THEN
input_radar%oerror(:,:,:,input_radar%iv3d_wind)=CONSTANT_VR_ERROR
ENDIF
IF( input_radar%iv3d_prh .GT. 0 )THEN
input_radar%oerror(:,:,:,input_radar%iv3d_prh)=CONSTANT_PSEUDORH_ERROR
ENDIF

END SUBROUTINE radar_observation_error

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DUPLICATE A RADAR ARRAY STRUCTURE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE radar_copy( input_radar , copy_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(IN) :: input_radar
TYPE(RADAR), INTENT(INOUT) :: copy_radar
!All the data in input radar will be transfered to copy_radar

copy_radar%na = input_radar%na
copy_radar%ne = input_radar%ne
copy_radar%nr = input_radar%nr
copy_radar%range_resolution= input_radar%range_resolution
copy_radar%radar_type = input_radar%radar_type
!copy_radar%beam_wid_h=input_radar%beam_wid_h
!copy_radar%beam_wid_v=input_radar%beam_wid_v
copy_radar%year=input_radar%year
copy_radar%month=input_radar%month
copy_radar%day=input_radar%day
copy_radar%hour=input_radar%hour
copy_radar%minute=input_radar%minute
copy_radar%second=input_radar%second
copy_radar%missing=input_radar%missing
copy_radar%lon0=input_radar%lon0
copy_radar%lat0=input_radar%lat0
copy_radar%z0=input_radar%z0
copy_radar%lambda=input_radar%lambda
copy_radar%total_attenuation_factor=input_radar%total_attenuation_factor
copy_radar%iv3d_wind=input_radar%iv3d_wind
copy_radar%iv3d_ref=input_radar%iv3d_ref


IF( ALLOCATED( input_radar % beam_wid_h ) )THEN
  ALLOCATE(copy_radar % beam_wid_h(input_radar%na ) )
  copy_radar%beam_wid_h=input_radar%beam_wid_h
ENDIF
IF( ALLOCATED( input_radar % beam_wid_v ) )THEN
  ALLOCATE(copy_radar % beam_wid_v(input_radar%ne ) )
  copy_radar%beam_wid_v=input_radar%beam_wid_v
ENDIF

IF( ALLOCATED( input_radar % azimuth ) )THEN
  ALLOCATE(copy_radar % azimuth(input_radar%na ) )
  copy_radar%azimuth=input_radar%azimuth
ENDIF
IF( ALLOCATED( input_radar % rrange ) )THEN
  ALLOCATE(copy_radar % rrange( input_radar%nr ) )
  copy_radar%rrange=input_radar%rrange
ENDIF
IF( ALLOCATED( input_radar % elevation ) )THEN
  ALLOCATE(copy_radar % elevation(input_radar%ne ) )
  copy_radar%elevation=input_radar%elevation
ENDIF
IF( ALLOCATED( input_radar % local_elevation ) )THEN
  ALLOCATE(copy_radar % local_elevation( input_radar%nr, input_radar%ne ) )
  copy_radar%local_elevation=input_radar%local_elevation
ENDIF
IF( ALLOCATED( input_radar % radarv3d ) )THEN
  ALLOCATE(copy_radar% radarv3d( input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )
  copy_radar%radarv3d=input_radar%radarv3d
ENDIF
IF( ALLOCATED( input_radar % radarv3d_model ) )THEN
  ALLOCATE(copy_radar% radarv3d_model( input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )
  copy_radar%radarv3d_model=input_radar%radarv3d_model
ENDIF
IF( ALLOCATED( input_radar % oerror ) )THEN
  ALLOCATE(copy_radar%oerror( input_radar%na,input_radar%nr,input_radar%ne,input_radar%nv3d) )
  copy_radar%oerror=input_radar%oerror
ENDIF
IF( ALLOCATED( input_radar % qcflag ) )THEN
  ALLOCATE(copy_radar%qcflag( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%qcflag=input_radar%qcflag
ENDIF
IF( ALLOCATED( input_radar % attenuation ) )THEN
  ALLOCATE(copy_radar%attenuation( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%attenuation=input_radar%attenuation
ENDIF
IF( ALLOCATED( input_radar % lat ) )THEN
  ALLOCATE(copy_radar%lat( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%lat=input_radar%lat
ENDIF
IF( ALLOCATED( input_radar % lon ) )THEN
  ALLOCATE(copy_radar%lon( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%lon=input_radar%lon
ENDIF
IF( ALLOCATED( input_radar % z ) )THEN
  ALLOCATE(copy_radar%z( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%z=input_radar%z
ENDIF
IF( ALLOCATED( input_radar % p ) )THEN
  ALLOCATE(copy_radar%p( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%p=input_radar%p
ENDIF
IF( ALLOCATED( input_radar % i ) )THEN
  ALLOCATE(copy_radar%i( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%i=input_radar%i
ENDIF
IF( ALLOCATED( input_radar % j ) )THEN
  ALLOCATE(copy_radar%j( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%j=input_radar%j
ENDIF
IF( ALLOCATED( input_radar % k ) )THEN
  ALLOCATE(copy_radar%k( input_radar%na,input_radar%nr,input_radar%ne) )
  copy_radar%k=input_radar%k
ENDIF
IF( ALLOCATED( input_radar % sounding_t ) )THEN
  ALLOCATE(copy_radar%sounding_t( input_radar%nsounding) )
  copy_radar%sounding_t = input_radar%sounding_t
ENDIF
IF( ALLOCATED( input_radar % sounding_q ) )THEN
  ALLOCATE(copy_radar%sounding_q( input_radar%nsounding) )
  copy_radar%sounding_q = input_radar%sounding_q
ENDIF
IF( ALLOCATED( input_radar % sounding_z ) )THEN
  ALLOCATE(copy_radar%sounding_z( input_radar%nsounding) )
  copy_radar%sounding_z = input_radar%sounding_z
ENDIF

END SUBROUTINE radar_copy

        END MODULE COMMON_RADAR_TOOLS
