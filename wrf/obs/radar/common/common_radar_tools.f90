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
  INTEGER     , ALLOCATABLE :: qcflag(:,:,:)
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
  INTEGER  :: iv3d_vr , iv3d_ref , iv3d_prh , iv3d_vrwidth , iv3d_snr 
  REAL(r_size)        :: total_attenuation_factor
  REAL(r_size), ALLOCATABLE :: distance_to_radar(:,:)

 
  END TYPE 

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
  input_radar%iv3d_vr=-9
  input_radar%iv3d_prh=-9
  input_radar%iv3d_vrwidth=-9
  input_radar%iv3d_snr=-9


  IF ( input_radar % radar_type .EQ. 1 )THEN
  iunit=99

  !Define data order for the PAWR data 
  input_radar%iv3d_ref=1    !Reflectivity
  input_radar%iv3d_vr=2   !Radial velocity

  !Get radar dimension from the input file.
  !Skip data only read array dimensions and grid information.

  !Input from PAWR data is in little_endian.
  OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential', &
       CONVERT='little_endian')
  
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
    input_radar%iv3d_vr=2   !Radial velocity
    input_radar%iv3d_vrwidth=3 !Velocity spectral width.

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
SUBROUTINE radar_write_file( input_radar , reflectivity , wind , qcflag , attenuation ,  output_file )
  IMPLICIT NONE
  TYPE(RADAR), INTENT(INOUT) :: input_radar
  INTEGER                    :: nvar
  CHARACTER(LEN=*),INTENT(IN):: output_file
  REAL(r_sngl)               :: buf4(input_radar%na,input_radar%nr)
  REAL(r_size), INTENT(IN)   :: reflectivity(input_radar%na,input_radar%nr,input_radar%ne)
  REAL(r_size), INTENT(IN)   :: wind(input_radar%na,input_radar%nr,input_radar%ne)
  REAL(r_size), INTENT(IN)   :: attenuation(input_radar%na,input_radar%nr,input_radar%ne)
  INTEGER, INTENT(IN)        :: qcflag(input_radar%na,input_radar%nr,input_radar%ne)
  INTEGER      :: iunit , ivar
  INTEGER      :: ia , ir , ie

  nvar=2
  input_radar%total_attenuation_factor=0.0d0

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
                  REAL( undefs  )

    !READ DATA DIMESION HEADER
    WRITE(iunit)input_radar%na , input_radar%nr , input_radar%ne , nvar


    !READ DIMENSIONS DATA
    WRITE(iunit)REAL( input_radar%azimuth  , r_sngl )
    WRITE(iunit)REAL( input_radar%rrange   , r_sngl )
    WRITE(iunit)REAL( input_radar%elevation, r_sngl )

    !READ TOTAL ATTENUATION FACTOR
    WRITE(iunit)REAL( input_radar%total_attenuation_factor , r_sngl)

    !REFLECTIVITY
    DO ie=1,input_radar % ne
        buf4=REAL( reflectivity(:,:,ie) , r_sngl )
        WHERE( reflectivity(:,:,ie) == UNDEF )buf4=UNDEFS
        WRITE(iunit)buf4
    ENDDO

    !RADIAL VELOCITY
    DO ie=1,input_radar % ne
        buf4=REAL( wind(:,:,ie) , r_sngl )
        WHERE( wind(:,:,ie) == UNDEF )buf4=UNDEFS
        WRITE(iunit)buf4
    ENDDO

    !QC FLAG
    DO ie=1,input_radar % ne
      WRITE(iunit) REAL( qcflag(:,:,ie) , r_sngl )
    ENDDO
    !ATTENUATION
    DO ie=1,input_radar % ne
       buf4=REAL( attenuation(:,:,ie) , r_sngl )
       WHERE( attenuation(:,:,ie) == UNDEF )buf4=UNDEFS
       WRITE(iunit)buf4
    ENDDO

    CLOSE(iunit)

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
  OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='sequential', &
       CONVERT='little_endian')

  !Define data order for the PAWR data 
  input_radar%iv3d_ref=1    !Reflectivity
  input_radar%iv3d_vr=2   !Radial velocity

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

  ALLOCATE( input_radar%beam_wid_h( 1 ) )
  ALLOCATE( input_radar%beam_wid_v( 1 ) )

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
    input_radar%iv3d_vr=2   !Radial velocity
    input_radar%iv3d_vrwidth=3 !Velocity spectral width.

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
! Read radar data
! This function reads one volume of radar data.
!-----------------------------------------------------------------------
SUBROUTINE radar_read_data_old( input_radar , input_file )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT)       :: input_radar
CHARACTER(LEN=*) ,INTENT(IN)     :: input_file
INTEGER                    :: iunit , ivar , inttmp
INTEGER                    :: ia , ir , ie
REAL(r_sngl)               :: buf4
REAL(r_size)               :: mindbz
INTEGER                    :: record
REAL(r_size),ALLOCATABLE   :: tmp_data(:,:,:,:) , tmp_elev(:)
INTEGER                    :: total_lev

mindbz=10.0d0*log10(minz)

input_radar%nv3d=1    !Only reflectivity will be read.
input_radar%iv3d_ref=1
!Reads radar information and data at the same time

record=1
iunit=98
!First select procedure according to the radar type.

  OPEN(iunit,FILE=input_file,FORM='unformatted',ACCESS='stream')

  IF(ALLOCATED(input_radar%beam_wid_h) )DEALLOCATE( input_radar%beam_wid_h)
  ALLOCATE( input_radar%beam_wid_h( 1 ) )
  IF(ALLOCATED(input_radar%beam_wid_v) )DEALLOCATE( input_radar%beam_wid_v)
  ALLOCATE( input_radar%beam_wid_v( 1 ) )

  !READ RADAR LOCATION AND CHARACTERISTICS HEADER

  READ(iunit) buf4
  input_radar%lon0            =REAL(buf4,r_size)
  READ(iunit) buf4
  input_radar%lat0            =REAL(buf4,r_size)
  READ(iunit) buf4
  input_radar%z0              =REAL(buf4,r_size)
  READ(iunit) buf4
  input_radar%beam_wid_h      =REAL(buf4,r_size)  !In PAWR this is a constant.
  READ(iunit) buf4
  input_radar%beam_wid_v      =REAL(buf4,r_size)  !In PAWR this is a constant.


  !READ DATA DIMESION HEADER
  READ(iunit)input_radar%na
  READ(iunit)input_radar%nr
  READ(iunit)input_radar%ne

  !WRITE(6,*) input_radar%na , input_radar%nr , input_radar%ne ,
  !input_radar%lon0 , input_radar%lat0 , input_radar%z0 , input_radar%beam_wid_h
  !, input_radar%beam_wid_v
  !STOP

  IF(ALLOCATED(input_radar%radarv3d) )DEALLOCATE(input_radar%radarv3d)
  ALLOCATE( input_radar%radarv3d(input_radar%na , input_radar%nr ,input_radar%ne , input_radar%nv3d ) )

    DO ie=1,input_radar % ne
     DO ir=1,input_radar % nr
      DO ia=1,input_radar % na
      READ(iunit) buf4
      input_radar%radarv3d(ia,ir,ie,1)=REAL(buf4,r_size)
      ENDDO
     ENDDO
    ENDDO

  READ(iunit)input_radar%na
  READ(iunit)input_radar%ne

  !READ DIMENSIONS DATA
  IF(ALLOCATED(input_radar%azimuth) )DEALLOCATE( input_radar%azimuth)
  ALLOCATE( input_radar % azimuth  ( input_radar % na ) )
  IF(ALLOCATED(input_radar%rrange)  )DEALLOCATE( input_radar%rrange)
  ALLOCATE( input_radar % rrange   ( input_radar % nr ) )
  IF(ALLOCATED(input_radar%elevation))DEALLOCATE( input_radar%elevation)
  ALLOCATE( input_radar % elevation( input_radar % ne ) )

   DO ia=1,input_radar % na
    DO ie=1,input_radar % ne
      READ(iunit) buf4
      IF( ie == 1 ) input_radar%azimuth(ia)=REAL(buf4,r_size)
    ENDDO
   ENDDO

   DO ia=1,input_radar % na
    DO ie=1,input_radar % ne
      READ(iunit) buf4
      IF( ia == 1 ) input_radar%elevation(ie)=REAL(buf4,r_size)
    ENDDO
   ENDDO


  DO ir=1,input_radar%nr
   input_radar%rrange(ir)=(ir-1)*100
  ENDDO

  CLOSE(iunit)


  !Some levels are repeated in PAWR data, the reflectivity corresponding to
  !these levels
  !will first be averaged and then the extra levels will be eliminated.

  total_lev=1
  DO ie=2,input_radar%ne
    IF( input_radar%elevation(ie) /= input_radar%elevation(ie-1) ) total_lev=total_lev + 1
  ENDDO

  ALLOCATE(tmp_data(input_radar%na,input_radar%nr,total_lev,input_radar%nv3d) )
  ALLOCATE(tmp_elev( total_lev ) )

  total_lev=1
  tmp_data(:,:,1,:)=input_radar%radarv3d(:,:,1,:)
  tmp_elev(1)=input_radar%elevation(1)
  DO ie=2,input_radar%ne
     IF( input_radar%elevation(ie) == input_radar%elevation(ie-1) )THEN
      tmp_data(:,:,total_lev,:)=0.5*( input_radar%radarv3d(:,:,ie,:)+input_radar%radarv3d(:,:,ie-1,:) )
     ELSE
      total_lev=total_lev+1
      tmp_data(:,:,total_lev,:)=input_radar%radarv3d(:,:,ie,:)
      tmp_elev(total_lev)=input_radar%elevation(ie)
     ENDIF
  ENDDO

  input_radar%ne = total_lev

  DEALLOCATE( input_radar%radarv3d , input_radar%elevation )
  ALLOCATE( input_radar%radarv3d(input_radar%na, input_radar%nr , input_radar%ne , input_radar%nv3d ) )
  ALLOCATE( input_radar%elevation( input_radar%ne ) )


  input_radar%radarv3d=tmp_data
  input_radar%elevation=tmp_elev

  DEALLOCATE( tmp_data , tmp_elev )

END SUBROUTINE radar_read_data_old

!-----------------------------------------------------------------------
! Compute the radar grid in LAT , LON , Z
!-----------------------------------------------------------------------
SUBROUTINE radar_georeference( input_radar )
IMPLICIT NONE
TYPE(RADAR), INTENT(INOUT) :: input_radar
REAL(r_size)  :: ke 
REAL(r_size)  :: tmp1 , tmp2 
INTEGER       :: ia , ir , ie

ke=(4d0/3d0)

ALLOCATE( input_radar%z( input_radar%na , input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%lat(input_radar%na, input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%lon(input_radar%na, input_radar%nr , input_radar%ne ) )
ALLOCATE( input_radar%local_elevation( input_radar%nr , input_radar%ne    ) )
ALLOCATE( input_radar%distance_to_radar( input_radar%nr , input_radar%ne  ) )
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
       input_radar%distance_to_radar(ir,ie)=ke*Re*asin( input_radar%rrange(ir) * cos(input_radar%elevation(ie)*deg2rad) / (ke*Re) )
       
       DO ia=1,input_radar%na
        CALL com_ll_arc_distance(input_radar%lon0,input_radar%lat0,                                  &
             input_radar%distance_to_radar(ir,ie),input_radar%azimuth(ia),input_radar%lon(ia,ir,ie), &
             input_radar%lat(ia,ir,ie))
       ENDDO
    ENDDO
   ENDDO


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
copy_radar%iv3d_vr=input_radar%iv3d_vr
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
