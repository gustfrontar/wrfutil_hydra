PROGRAM obs_gen
!=======================================================================
!
! [PURPOSE:] Create observations from a "nature" run.
!
! [HISTORY:]
!   01/01/2013 Juan Ruiz  created based on the code provided by Takemasa
!   Miyoshi.
!
!=======================================================================
  USE common
  USE common_wrf
  USE common_namelist
  USE common_wrf_to_bufr

  IMPLICIT NONE


  !Get configuration flags.
  CALL get_namelist_vars()

!-----------------------------------------------------------------------
! From model to obs steps. 
!-----------------------------------------------------------------------
!
  !Get model grid info.
  CALL set_common_wrf(guesfile)

  !Get observation number.
  CALL get_obs_num()

  ALLOCATE(v3d(nlon,nlat,nlev,nv3d))
  ALLOCATE(v2d(nlon,nlat,nv2d))

  !Read model data.
  CALL read_grd(guesfile,v3d,v2d)

  !Get observation location.
  CALL get_obs_ijk()

  !Apply operator from model space to observation space.
  CALL model_to_obs()

  !Write observations in LETKF format.
  CALL write_obs(outfile)

STOP
END PROGRAM obs_gen



