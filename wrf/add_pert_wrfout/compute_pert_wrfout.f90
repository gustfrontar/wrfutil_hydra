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
 REAL(r_size),ALLOCATABLE :: pert3da(:,:,:,:)
 REAL(r_size),ALLOCATABLE :: pert2da(:,:,:)
 REAL(r_size),ALLOCATABLE :: pert3db(:,:,:,:)
 REAL(r_size),ALLOCATABLE :: pert2db(:,:,:)
 REAL(r_size),ALLOCATABLE :: pert3d_final(:,:,:,:)
 REAL(r_size),ALLOCATABLE :: rpert2d(:,:,:),rpert3d(:,:,:,:)
 REAL(r_size), ALLOCATABLE :: levels(:) , levels_final(:)
 CHARACTER(15)            :: met_ema1='input_filea1.nc',met_ema2='input_filea2.nc'
 CHARACTER(15)            :: met_emb1='input_fileb1.nc',met_emb2='input_fileb2.nc'
 CHARACTER(15)            :: ctrl_met_em='ctrl_met_em.nc'    !Where perturbation will be dadded
 LOGICAL                  :: USE_RANDOM_PERT
 INTEGER :: ilat, ilon ,ilev ,ivar , kk  , nlev_pert
 REAL(r_size) :: rmse

 REAL(r_size) :: tmp
 CHARACTER(10) :: input_par
 INTEGER       :: current_time , max_time
 REAL(r_size)  :: time_factor

!-----------------------------------------------------------------------------
! GET PERTURBATION AMPLITUDE (thes parameters should go to a namelist)
!-----------------------------------------------------------------------------

  CALL GETARG ( 1, input_par )
  READ(input_par,*)amp_factor
  WRITE(*,*)"AMP_FACTOR is ",amp_factor
  CALL GETARG ( 2, input_par )
  READ(input_par,*)random_amp_factor
  WRITE(*,*)"RANDOM AMP_FACTOR is ",random_amp_factor
  CALL GETARG ( 3, input_par )
  READ(input_par,*)current_time
  WRITE(*,*)"CURRENT TIME is ",current_time
  CALL GETARG ( 4, input_par )
  READ(input_par,*)max_time
  WRITE(*,*)"MAX TIME is ",max_time

  IF( current_time .GT. max_time )THEN
    write(*,*)"[Error] Max time = ",max_time," and Current time = ",current_time
    stop
  ENDIF  

!-----------------------------------------------------------------------------
! RANDOM PERTURBATION PARAMETERS (this should go to a namelist )
!-----------------------------------------------------------------------------

!Wheter variables will be perturbed
PERTURB_T = .TRUE.
PERTURB_RH = .TRUE.
PERTURB_WIND = .TRUE.

!Perturbation amplitude.
PERTURB_T_AMP = 1.0d0
PERTURB_RH_AMP = 10.0d0
PERTURB_WIND_AMP = 1.0d0

!Perturbation lenght scale (horizontal)
PERTURB_T_SCLH = 10000d0
PERTURB_RH_SCLH = 10000d0
PERTURB_WIND_SCLH = 10000d0

!Perturbation length scale (vertical)
PERTURB_T_SCLV = 10000d0
PERTURB_RH_SCLV = 10000d0
PERTURB_WIND_SCLV = 10000d0

!-----------------------------------------------------------------------------
! BALANCED PERTURBATION PARAMETERS (this should go to a namelist )
!-----------------------------------------------------------------------------

PERTURB_ATMOSPHERE= .TRUE.  !If atmospheric variables will be perturbed
PERTURB_SST       = .TRUE.  !If sst will be perturbed
PERTURB_SOIL      = .TRUE.  !If soiltemp and soilmoisture will be perturbed

!--------------------------------------------------------------------------------

!Get model domain from the first ensemble member.
CALL set_common_wrf(met_ema1)

ALLOCATE( gues3d(nlon,nlat,nlev,nv3d,2),gues2d(nlon,nlat,nv2d,2) )
ALLOCATE( pert3da(nlon,nlat,nlev,nv3d),pert2da(nlon,nlat,nv2d) )
ALLOCATE( pert3db(nlon,nlat,nlev,nv3d),pert2db(nlon,nlat,nv2d) )

ALLOCATE( levels(nlev) )

!Get perturbation at the initial time in the window.
 CALL read_grd(met_ema1,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 CALL read_grd(met_ema2,gues3d(:,:,:,:,2),gues2d(:,:,:,2))

 levels=gues3d(1,1,:,iv3d_pres,1)  !Get input levels.

 CALL PERTURB_MEMBER_BALANCED(gues3d(:,:,:,:,1),gues2d(:,:,:,1),   &
                              gues3d(:,:,:,:,2),gues2d(:,:,:,2),   &
                              pert3da(:,:,:,:),pert2da(:,:,:)  )

!Get perturbation at the final time in the window.
 CALL read_grd(met_emb1,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 CALL read_grd(met_emb2,gues3d(:,:,:,:,2),gues2d(:,:,:,2))


 levels=gues3d(1,1,:,iv3d_pres,1)  !Get input levels.

 CALL PERTURB_MEMBER_BALANCED(gues3d(:,:,:,:,1),gues2d(:,:,:,1),   &
                              gues3d(:,:,:,:,2),gues2d(:,:,:,2),   &
                              pert3da(:,:,:,:),pert2da(:,:,:)  )

 !Interpolate perturbation linearly in time
 time_factor=( REAL(current_time,r_size) / REAL(max_time,r_size) )
 pert3da=(1.0d0-time_factor)*pert3da + (time_factor)*pert3db
 pert2da=(1.0d0-time_factor)*pert2da + (time_factor)*pert2db


 !ADD A RANDOM PERTURBATION TO THE TIME INTERPOLATED BALANCED PERTURBATION.
 IF( random_amp_factor > 0.0d0 )THEN
  ALLOCATE(rpert3d(nlon,nlat,nlev,nv3d),rpert2d(nlon,nlat,nv2d))
  CALL PERTURB_MEMBER_RANDOM(gues3d(:,:,:,:,1),gues2d(:,:,:,1),rpert3d,rpert2d)
  pert3da=pert3da+random_amp_factor*rpert3d
  pert2da=pert2da+random_amp_factor*rpert2d
  DEALLOCATE(rpert3d,rpert2d)
 ENDIF

 nlev_pert=nlev

 !OPEN CONTROL FILE 
 CALL set_common_wrf(ctrl_met_em)
 ALLOCATE( pert3d_final(nlon,nlat,nlev,nv3d)  ) 
 ALLOCATE( levels_final(nlev) )
 DEALLOCATE( gues3d , gues2d)
 ALLOCATE( gues3d(nlon,nlat,nlev,nv3d,1),gues2d(nlon,nlat,nv2d,1) )

 CALL read_grd(ctrl_met_em,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 
 !Lets keep only the levels that are the same as the final levels.
 !Always keep the first level for all the variables.

 levels_final=gues3d(1,1,:,iv3d_pres,1)  !Get final levels.

 pert3d_final=0.0d0

 !First level corresponds to surface pressure and surface variables.
 pert3d_final(:,:,1,:)=pert3da(:,:,1,:)

 WRITE(*,*)"Levels final",levels_final
 WRITE(*,*)"Levels ",levels
 

 DO ilev=2,nlev
  DO kk=2,nlev_pert
    !WRITE(*,*)levels_final(ilev),levels(kk),levels_final(ilev)-levels(kk),tiny(levels(kk))
    IF( abs(levels(kk) - levels_final(ilev)) <= tiny(levels(kk)) )THEN
    !  WRITE(*,*)"Match!",ilev,kk,levels_final(ilev),levels(kk)
      pert3d_final(:,:,ilev,:)=pert3da(:,:,kk,:)
      exit
    ENDIF 
    IF( kk == nlev_pert )THEN
      WRITE(*,*)"FATAL!!!: REQUESTED LEVEL NOT FOUND IN INPUT DATA!"
      WRITE(*,*)"Current level ", levels_final(ilev)
      WRITE(*,*)"Levels list",levels 
      STOP
    ENDIF
  ENDDO
 ENDDO

 gues3d(:,:,:,:,1)=gues3d(:,:,:,:,1) + pert3d_final  
 gues2d(:,:,:,1)  =gues2d(:,:,:,1)   + pert2da

 !Check that positive variables are stil possitive
 CALL check_pert(gues3d,gues2d)


 CALL write_grd(ctrl_met_em,gues3d(:,:,:,:,1),gues2d(:,:,:,1))

END PROGRAM COMPUTE_PERTURBATION

