PROGRAM TIME_INTERP_MET_EM
!=======================================================================
!
! [PURPOSE:] Linearly interpolate two met_ems in time.
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
 CHARACTER(15)            :: met_em1='input_file1.nc',met_em2='input_file2.nc'
 CHARACTER(15)            :: ctrl_met_em='ctrl_met_em.nc'    !Where perturbation will be dadded
 LOGICAL                  :: USE_RANDOM_PERT
 CHARACTER(19) :: cdate , date1 , date2
 INTEGER :: ilat, ilon ,ilev ,ivar , kk  , nlev_pert
 REAL(r_size) :: rmse

 REAL(r_size) :: tmp
 CHARACTER(10) :: input_par
 INTEGER       :: current_time , max_time
 REAL(r_size)  :: time_factor

 INTEGER :: iy , im , id , ih
 REAL(r_size) :: sec , tai1 , tai2 , ctai

!-----------------------------------------------------------------------------
! GET INPUT PARAMETERS (thes parameters should go to a namelist)
!-----------------------------------------------------------------------------


  CALL GETARG ( 1, input_par )
  READ(input_par,*)cdate
  WRITE(*,*)"CURRENT TIME is ",current_time  !YYYY-MM-DD_HH:MM:SS

!--------------------------------------------------------------------------------

!Get model domain from the first ensemble member.
CALL set_common_wrf(met_ema1)

ALLOCATE( gues3d(nlon,nlat,nlev,nv3d,2),gues2d(nlon,nlat,nv2d,2) )


ALLOCATE( levels(nlev) )

!Get perturbation at the initial time in the window.
 CALL read_grd(met_em1,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 CALL read_grd(met_em2,gues3d(:,:,:,:,2),gues2d(:,:,:,2))

 CALL get_date(met_em1,date1)
 CALL get_date(met_em2,date2)

!Convert UTC to tai  
 READ(date1(1:4),*)iy
 READ(date1(6:7),*)im
 READ(date1(9:10),*)id
 READ(date1(12:13),*)ih
 READ(date1(15:16),*)im
 READ(date1(18:19),*)sec
 CALL com_utc2tai(iy,im,id,ih,imin,sec,tai1)

!Convert UTC to tai  
 READ(date2(1:4),*)iy
 READ(date2(6:7),*)im
 READ(date2(9:10),*)id
 READ(date2(12:13),*)ih
 READ(date2(15:16),*)im
 READ(date2(18:19),*)sec
 CALL com_utc2tai(iy,im,id,ih,imin,sec,tai2)

!Convert UTC to tai  
 READ(cdate(1:4),*)iy
 READ(cdate(6:7),*)im
 READ(cdate(9:10),*)id
 READ(cdate(12:13),*)ih
 READ(cdate(15:16),*)im
 READ(cdate(18:19),*)sec
 CALL com_utc2tai(iy,im,id,ih,imin,sec,ctai)

 !Interpolate perturbation linearly in time
 IF ( ctai < ctai2 .and. ctai > ctai1)
     IF( tai2 - tai1 .NE. 0.0d0 )THEN
      time_factor= (ctai-tai1) / (tai2 - tai1)
     ELSE
      time_factor=1.0d0
     ENDIF
 ELSE
     WRITE(*,*)"[Error]: The requested date is not within the time range, "
     STOP
 END
 gues3d(:,:,:,:,1)=(1.0d0-time_factor)*gues3d(:,:,:,:,1) + (time_factor)*gues3d(:,:,:,:,2)
 gues2d(:,:,:,1)=(1.0d0-time_factor)*gues2d(:,:,:,1) + (time_factor)*gues2d(:,:,:,2)

 CALL write_grd(ctrl_met_em,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 
 CALL put_date(ctrl_met_em, cdate)

END PROGRAM TIME_INTERP_MET_EM

