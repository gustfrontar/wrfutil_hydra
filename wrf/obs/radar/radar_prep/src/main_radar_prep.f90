PROGRAM MAIN_RADAR_PREP
!=======================================================================
!
! [PURPOSE:] Main program of  RADAR_PREP
! This program interpolates radar data to model grid using box averaging.
! In this version spatially correlated random errors are generated to 
! simulate observational error. (For OSSE experiments)
! TODO: 
! -Implement radar thinning in this same program.
! -Implement superobbing to the model grid (this will requiere
!  model specific routines)
! -Currently ontly PAWR IO is supported. Implement IO for other radar tyeps
!
! [HISTORY:]
!   07/10/2014 Juan Ruiz  created
!   15/02/2016 Juan Ruiz modified, added namelist input.
!=======================================================================

USE common
USE common_superobbing
USE common_radar_tools
USE common_namelist

 IMPLICIT NONE

 REAL(r_size)             :: rtimer00 , rtimer

 CHARACTER*20  :: file_name_radar
 INTEGER       :: ia , ir , ie , iobs , kk , ii , jj

 REAL(r_sngl) wk(7)

 TYPE(RADAR) :: PAWR  !We can have more than one radar 

CALL CPU_TIME(rtimer00)

CALL read_namelist()

WRITE(*,*)"Welcome to radar prep routine "
WRITE(*,*)"---------------------------------------"
WRITE(*,*)"Here are the most important parameters "
WRITE(*,*)"---------------------------------------"
WRITE(*,*)"DX=                           ",dx
WRITE(*,*)"DZ=                           ",dz
WRITE(*,*)"MAXRANGE=                     ",maxrange
WRITE(*,*)"MAXZ=                         ",maxz
WRITE(*,*)"DBZ_TO_POWER=                 ",dbz_to_power
WRITE(*,*)"ERROR_REF=                    ",error_ref
WRITE(*,*)"ERROR_VR =                    ",error_vr
WRITE(*,*)"id_ref_obs=                   ",id_ref_obs
WRITE(*,*)"id_ref_vr=                    ",id_vr_obs
WRITE(*,*)"DEBUG_OUTPUT=                 ",debug_output
if( osse_exp) THEN !Add osse configuration variables.
WRITE(*,*)"This is an OSSE experiment"
WRITE(*,*)"SIMULATE_ERROR=               ",simulate_error
WRITE(*,*)"ERROR_REF_SCLX=               ",error_ref_sclx
WRITE(*,*)"ERROR_VR_SCLX=                ",error_vr_sclx
WRITE(*,*)"ERROR_REF_SCLZ=               ",error_ref_sclz
WRITE(*,*)"ERROR_VR_SCLZ=                ",error_vr_sclz

endif



!!!!! RADAR DATA INPUT !!!!!!!!!!
 file_name_radar='radar_input.grd'  
 PAWR % radar_type = 1 

CALL radar_read_data( PAWR , file_name_radar , endian )

!do ii=1,PAWR%na
!  do jj=1,PAWR%nr
!    do kk=1,PAWR%ne
!      if( PAWR%attenuation(ii,jj,kk) < attenuation_threshold )then
!       WRITE(*,*)PAWR%attenuation(ii,jj,kk)
!      endif
!    enddo
!  enddo
!enddo

CALL radar_georeference( PAWR )

WRITE(*,*)"---------------------------------------"
WRITE(*,*)"Procesing data for date: "
WRITE(*,*)PAWR%year,PAWR%month,PAWR%day,PAWR%hour,PAWR%minute,PAWR%second
WRITE(*,*)"---------------------------------------"

!Use radar location to define the grid  for superobbing
CALL define_grid( PAWR )


CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(READ RADAR):',rtimer,rtimer-rtimer00
rtimer00=rtimer
!!!!!! END OF RADAR DATA INPUT !!!!!!!

CALL radar_superobbing( PAWR )

!
CALL CPU_TIME(rtimer)
WRITE(6,'(A,2F10.2)') '### TIMER(THINNING AND WRITING):',rtimer,rtimer-rtimer00
rtimer00=rtimer
   


END PROGRAM MAIN_RADAR_PREP

