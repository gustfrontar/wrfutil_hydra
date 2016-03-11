PROGRAM COMPUTE_PERTURBATION
!=======================================================================
!
! [PURPOSE:] Generate an ensemble of perturbed WRF runs.
!
! [HISTORY:]
!   02/19/2014 Juan Ruiz created
!   13/02/2016 Major bug corrected SLP was not being perturbed.
!=======================================================================
USE common
USE common_wrf
USE common_perturb_ensemble

 IMPLICIT NONE
 REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:,:)
 REAL(r_size),ALLOCATABLE :: gues2d(:,:,:,:)
 REAL(r_size), ALLOCATABLE :: levels(:) 
 CHARACTER(15)            :: inputfile='wrfoutXXX.nc'    !Where perturbation will be dadded
 INTEGER :: ilat, ilon ,ilev ,ivar , kk  , nlev_pert , it
 REAL(r_size) :: rmse

 REAL(r_size) :: tmp
 CHARACTER(10) :: input_par


!-----------------------------------------------------------------------------
! RANDOM PERTURBATION PARAMETERS
!-----------------------------------------------------------------------------

 CALL read_namelist()

!--------------------------------------------------------------------------------

!Get model domain from the first ensemble member.
write(inputfile(7:9),'(I3.3)')1
CALL set_common_wrf(inputfile)


allocate(gues3d(nlon,nlat,nlev,nv3d,ntimes),gues2d(nlon,nlat,nv2d,ntimes)) 
!Read the first time (we need it to filter the perturbations
DO it=1,ntimes
 write(inputfile(7:9),'(I3.3)')it
 !OPEN CONTROL FILE
  CALL read_grd(inputfile,gues3d(:,:,:,:,it),gues2d(:,:,:,it))
ENDDO

 !ADD A RANDOM PERTURBATION TO THE TIME INTERPOLATED BALANCED PERTURBATION.
 CALL PERTURB_MEMBER_RANDOM(gues3d,gues2d)

 !Check that positive variables are stil possitive
 CALL check_pert(gues3d,gues2d)

DO it=1,ntimes
 write(inputfile(7:9),'(I3.3)')it
 CALL write_grd(inputfile,gues3d(:,:,:,:,it),gues2d(:,:,:,it))
ENDDO

END PROGRAM COMPUTE_PERTURBATION

