PROGRAM WRF_TO_RADAR
!=======================================================================
!
! [PURPOSE:] Main program of WRF_TO_RADAR
! This program interpolates WRF model data to a radar grid.
! The output is reflectivity factor and doppler wind.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!=======================================================================

USE common
USE common_wrf
USE common_wrf_to_radar
USE common_radar_tools

 IMPLICIT NONE

 REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:)
 REAL(r_size),ALLOCATABLE :: gues2d(:,:,:)

 REAL(r_size),ALLOCATABLE :: radar3d(:,:,:,:)
 REAL(r_size)             :: rtimer00 , rtimer
 INTEGER                  :: i

 CHARACTER*20  :: model_file_name='input_model.nc'
 CHARACTER*20  :: inputfile='input_radar.grd'    !Radar data input
 CHARACTER*20  :: outputfile='output_radar.grd'        !output files in radar format

 LOGICAL   ::   FAKE_RADAR=.TRUE. !We will create a non existant radar.

 REAL(r_size) :: tmp
 !--------------------------------------------FAKE RADAR PARAMETERS
 INTEGER , PARAMETER  :: radar_cen_x = 65 , radar_cen_y = 65 !Grid point position of radar.
 REAL(r_size), PARAMETER :: radar_z = 0.0d0
 REAL(r_size), PARAMETER :: radar_az_res = 1.2d0 , radar_r_res = 100.0d0 , radar_el_res = 1.0d0
 REAL(r_size), PARAMETER :: radar_min_az = 0.0d0 , radar_max_az = 360.0d0
 REAL(r_size), PARAMETER :: radar_min_r  = 0.0d0 , radar_max_r  = 60000
 REAL(r_size), PARAMETER :: radar_min_el = 0.0d0 , radar_max_el = 89.0d0


 !------------------------------------------------------------------


 INTEGER       :: ia , ir , ie

 TYPE(RADAR) :: RADAR_1  !We can have more than one radar 

CALL CPU_TIME(rtimer00)


CALL set_common_wrf(model_file_name)

ALLOCATE( gues3d(nlon,nlat,nlev,nv3d),gues2d(nlon,nlat,nv2d) )


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(SET MODEL):',rtimer,rtimer-rtimer00
rtimer00=rtimer


CALL read_grd(model_file_name,gues3d,gues2d)

gues3d(:,:,:,iv3d_ph)=gues3d(:,:,:,iv3d_ph)/gg

!-------------- DEBUG
!  gues3d(:,:,:,iv3d_u)=20*sin(gues3d(:,:,:,iv3d_ph)*pi/20e3 )
!  gues3d(:,:,:,iv3d_v)=0.0d0
!  gues3d(:,:,:,iv3d_w)=0.0d0
!  gues3d(:,:,:,iv3d_qr)=1e-4*(gues3d(:,:,:,iv3d_ph)**2 / 20e3 ** 2)
!  gues3d(:,:,:,iv3d_qs)=0.0d0
!  gues3d(:,:,:,iv3d_qg)=0.0d0
!  WRITE(6,*)gues3d(10,10,:,iv3d_u),gues3d(10,10,:,iv3d_ph)
!---------------------



CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(READ MODEL):',rtimer,rtimer-rtimer00
rtimer00=rtimer


!!!!! RADAR DATA INPUT !!!!!!!!!!
IF ( .NOT. FAKE_RADAR )THEN   
!Dealing with a real radar, lets read the radar location and charatheristics.
 inputfile='input_radar.grd'  
 RADAR_1 % radar_type =1
 CALL radar_read_data( RADAR_1 , inputfile ) 
ELSE
!Dealing with a fake radar, lets create the radar location and characteristics.
 RADAR_1 % lon0 = lon(radar_cen_x,radar_cen_y)
 RADAR_1 % lat0 = lat(radar_cen_x,radar_cen_y)

 RADAR_1 % nv3d = 2
 RADAR_1 % radar_type = 1
 
  tmp= ( radar_max_az-radar_min_az ) / radar_az_res
  RADAR_1 % na = NINT(tmp)
  tmp= ( radar_max_r-radar_min_r ) / radar_r_res
  RADAR_1 % nr = NINT(tmp)
  tmp= ( radar_max_el-radar_min_el ) / radar_el_res
  RADAR_1 % ne = NINT(tmp)

  ALLOCATE( RADAR_1 % azimuth( RADAR_1 % na) )
  ALLOCATE( RADAR_1 % rrange( RADAR_1 % nr) )
  ALLOCATE( RADAR_1 % elevation( RADAR_1 % ne) )

  DO i=1,RADAR_1 % na
     RADAR_1 % azimuth (i) = (i-1)*radar_az_res + radar_min_az
  ENDDO
  DO i=1,RADAR_1 % nr
     RADAR_1 % rrange (i) = (i-1)*radar_r_res + radar_min_r
  ENDDO
  DO i=1,RADAR_1 % ne
     RADAR_1 % elevation (i) = (i-1)*radar_el_res + radar_min_el
  ENDDO



  ALLOCATE(RADAR_1 % radarv3d(RADAR_1%na,RADAR_1%nr,RADAR_1%ne,RADAR_1%nv3d))
  ALLOCATE(RADAR_1 % qcflag(RADAR_1%na,RADAR_1%nr,RADAR_1%ne))
  ALLOCATE(RADAR_1 % attenuation(RADAR_1%na,RADAR_1%nr,RADAR_1%ne))
 
  RADAR_1 % radarv3d    = 0.0d0
  RADAR_1 % qcflag      = 0.0d0
  RADAR_1 % attenuation = 0.0d0 


  RADAR_1 % range_resolution = radar_r_res

  RADAR_1 %year  =2014.0
  RADAR_1 %month =7.0
  RADAR_1 %day   =16.0
  RADAR_1 %hour  =0.0
  RADAR_1 %minute=0.0
  RADAR_1 %second=0.0
  
  ALLOCATE(RADAR_1%beam_wid_h(RADAR_1%na),RADAR_1%beam_wid_v(RADAR_1%ne) )

  RADAR_1 %beam_wid_h     =0.1  !Assume this constant.
  RADAR_1%beam_wid_v      =0.1  !Assume this constant.
  RADAR_1%lambda          =3.3
  RADAR_1%missing         =9.99e9

  RADAR_1%total_attenuation_factor=0.0

  RADAR_1%nv3d = 2

  RADAR_1%iv3d_ref=1    !Reflectivity
  RADAR_1%iv3d_wind=2   !Radial velocity

  WRITE(6,*)' FAKE RADAR SUMMARY      '
  WRITE(6,*)' NA   ',RADAR_1%na
  WRITE(6,*)' NR   ',RADAR_1%nr
  WRITE(6,*)' NE   ',RADAR_1%ne
ENDIF


 WRITE(6,*)'RADAR LOCATION'
 WRITE(6,*)'RADAR LAT= ',RADAR_1 % lat0
 WRITE(6,*)'RADAR LON= ',RADAR_1 % lon0


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(SET RADAR):',rtimer,rtimer-rtimer00
rtimer00=rtimer

!Generate latitude, longitude an height for radar.
CALL radar_georeference( RADAR_1 )


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(GEOREFERENCE RADAR):',rtimer,rtimer-rtimer00
rtimer00=rtimer

!!!!!! END OF RADAR DATA INPUT !!!!!!!

CALL model_to_radar( RADAR_1 , gues3d , gues2d )


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(FORWARD OPERATOR):',rtimer,rtimer-rtimer00
rtimer00=rtimer


CALL radar_write_file( RADAR_1 , RADAR_1 % radarv3d_model , RADAR_1 % nv3d , outputfile )

CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(WRITE DATA):',rtimer,rtimer-rtimer00
rtimer00=rtimer


END PROGRAM WRF_TO_RADAR

