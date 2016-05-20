MODULE common_namelist
!=======================================================================
!
! [PURPOSE:] Namelist input for LETKF analysis
!
! [HISTORY:]
!   25/09/2009 Juan Ruiz  created
!
!=======================================================================
  USE common
  IMPLICIT NONE
  PUBLIC


 !NAMELIST VARIABLES FOR SUPEROBBING

 !Superobbing grid parameters
  REAL(r_size) :: DX=1000.0d0    !X resolution
  REAL(r_size) :: DZ=1000.0d0    !Y resolution
  REAL(r_size) :: MAXRANGE=60d3  !Max horizontal range
  REAL(r_size) :: MAXZ=20d3      !Max vertical range

  LOGICAL :: DBZ_TO_POWER = .false. !Whether reflectivity is dbz or not.
  LOGICAL :: DEBUG_OUTPUT=.false.
  
  REAL(r_size)  :: error_ref=5.0d0     !Reflectivity error in dBz (for output file)
  REAL(r_size)  :: error_vr=1.0d0      !Doppler error in m/s (for output file)

  INTEGER :: id_ref_obs=4001 !Id for reflectivity in ouput file.
  INTEGER :: id_vr_obs=4002  !Id for radial velocity in output file.

  LOGICAL :: latlon_coord=.true.    !Whether obs location will be provided in 
                                    !lat,lon,height or az,r,elev
  LOGICAL :: use_attenuation=.true. !Consider attenuation in superobbing
  LOGICAL :: use_qcflag=.true.      !Consider or not qc flag.
  LOGICAL :: use_vr_std=.true.           !If we are going to use the wind std threshold within each box.
  REAL(r_size) :: vr_std_threshold=5.0d0 !If wind variability within each superob is greather than this threshold the box is rejected.
  !IF USE_ATTENUATION == TRUE, then gates with estimated attenuation
  !greather than the threshold will be rejected. (this does not affect
  !the computation of attenuation in the forward operator)
  REAL(r_size)      :: attenuation_threshold=0.01 !0.1 is 10dbz, 0.5 is aprox 5 dbz.

  !Parameters for osse experiments
  LOGICAL  :: osse_exp = .false.
  LOGICAL  :: SIMULATE_ERROR=.false.  !Wheter we are going to add noise
  REAL(r_size) :: ERROR_REF_SCLX=5000d0   !SPATIAL SCALE OF ERROR IN METERS
  REAL(r_size) :: ERROR_VR_SCLX=5000d0    !SPATIAL SCALE OF ERRORS IN METERS
  REAL(r_size) :: ERROR_REF_SCLZ=5000d0   !VERTICAL SCALE OF ERRORS IN METERS
  REAL(r_size) :: ERROR_VR_SCLZ=5000d0    !VERTICLA SCALE OF ERRORS IN METERS

  CHARACTER(1)      :: endian='b' ! 'b' -> big_endian , 'l' -> little_endian

  CHARACTER(LEN=20) :: NAMELIST_FILE='./radarprep.namelist'

CONTAINS

SUBROUTINE READ_NAMELIST

IMPLICIT NONE
INTEGER :: IERR
LOGICAL :: file_exist
!Namelist declaration

NAMELIST / GENERAL / dx,dz,maxrange,maxz,dbz_to_power,error_ref,error_vr,                 &
                     id_ref_obs,id_vr_obs,debug_output , latlon_coord , use_attenuation , &
                     attenuation_threshold , use_qcflag
NAMELIST / OSSE    / osse_exp,simulate_error,      &
                     error_ref_sclx , error_vr_sclx , error_ref_sclz , error_vr_sclz , &
                     error_ref , error_vr 

INQUIRE(FILE=NAMELIST_FILE, EXIST=file_exist)

IF( .NOT. file_exist )THEN
  WRITE(*,*)"ERROR: Could not find namelist"
ENDIF

OPEN(54,FILE=NAMELIST_FILE)
READ(54,NML=GENERAL,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at GENERAL section"
WRITE(*,*)"Using default values"
ENDIF
REWIND(54)
READ(54,NML=OSSE,IOSTAT=IERR)
IF(IERR /=0)THEN
WRITE(*,*)"Warning!! Error during namelist reading at OSSE section"
WRITE(*,*)"Using default values"
ENDIF

CLOSE(54)


END SUBROUTINE READ_NAMELIST

END MODULE common_namelist
