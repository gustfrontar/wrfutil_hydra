PROGRAM CHANGE_WRF_LONGITUDES
!=======================================================================
!
! [PURPOSE:] Generate an ensemble of perturbed WRF runs.
!   Change WRF longitude convention from -180 to 180 to 0 360.
! [HISTORY:]
!   21/10/2016 Juan Ruiz created
!=======================================================================
USE common_change_wrf_longitudes

 IMPLICIT NONE
 CHARACTER(15)            :: input_file='input_file.nc'


 CALL set_common_wrf(input_file)

 CALL modify_longitudes(input_file)


END PROGRAM CHANGE_WRF_LONGITUDES

