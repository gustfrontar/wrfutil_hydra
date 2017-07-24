PROGRAM SCALE_TO_RADAR
!=======================================================================
!
! [PURPOSE:] Main program of SCALE_TO_RADAR
! This program interpolates SCALE model data to a radar grid.
! The output is reflectivity factor and doppler wind.
! IO is performed in CFRADIAL format.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!   03/17/2016 Juan Ruiz modified the code to add multiple radars and 
!                        and model files. A namelist file to control the
!                        interpolation process was also added.
!   10/07/2017 Juan Ruiz Adapted to SCALE model.

! TODO LIST:
! -Incorporate a C-band observation operator.
! -Incorporate the computation of attenuation.
! -Incorporate the effect of anomalous propagation.
! -Incorporate beam broadening effect. 
!=======================================================================

USE common
USE common_scale
USE common_scale_to_radar
USE common_radar_tools


implicit none

TYPE(RADAR) :: test_radar

!------------------------------------------------------------------

call get_namelist_vars

WRITE(6,*)'Interpolationg model data to a real radar location and geometry'
CALL interp_to_real_radar


END PROGRAM SCALE_TO_RADAR

