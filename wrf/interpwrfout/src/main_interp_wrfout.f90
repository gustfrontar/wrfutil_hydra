PROGRAM INTERP_WRFOUT
!=======================================================================
!
! [PURPOSE:] Interpolate two wrfouts linearly in time.
! [USE:] interp_wrfout.exe YYYY MM DD HH MM SS
!
! [HISTORY:]
!   23/02/2016 Juan Ruiz created
!=======================================================================
USE common
USE common_wrf

 IMPLICIT NONE
 REAL(r_size),ALLOCATABLE :: gues3d(:,:,:,:,:)
 REAL(r_size),ALLOCATABLE :: gues2d(:,:,:,:)
 CHARACTER(20)            :: wout1='input_file1.nc',wout2='input_file2.nc'
 CHARACTER(20)            :: woutint='input_fileint.nc'    !Where interpolated fields will be stored.
 LOGICAL                  :: USE_RANDOM_PERT
 INTEGER :: ilat, ilon ,ilev ,ivar , kk  , nlev_pert
 REAL(r_size) :: rmse

 REAL(r_size) :: tmp
 CHARACTER(14) :: input_par

 INTEGER      :: intyear , intmonth , intday , inthour , intminute
 REAL(r_size) :: intsecond

 INTEGER      :: year1 , month1 , day1 , hour1 , minute1
 REAL(r_size) :: second1

 INTEGER      :: year2 , month2 , day2 , hour2 , minute2
 REAL(r_size) :: second2

 REAL(r_size) :: tai1 , tai2 , inttai !Tai times
 
 REAL(r_size) :: w1,w2 !Interpolation weigths

!-----------------------------------------------------------------------------
! GET PERTURBATION AMPLITUDE (these parameters should go to a namelist)
!-----------------------------------------------------------------------------

  CALL GETARG ( 1, input_par )
  READ(input_par(1:4),'(I4)')intyear
  READ(input_par(5:6),'(I2)')intmonth
  READ(input_par(7:8),'(I2)')intday
  READ(input_par(9:10),'(I2)')inthour
  READ(input_par(11:12),'(I2)')intminute
  READ(input_par(13:14),'(F2.0)')intsecond
   
  WRITE(*,*)"The input date is : "
  WRITE(*,'(I5,I3.2,I3.2,I3.2,I3.2,F4.1)')intyear,intmonth,intday,inthour,intminute,intsecond

!--------------------------------------------------------------------------------

!Get model domain from the first ensemble member.
CALL set_common_wrf(wout1)

!Get the dates from the two input files and set the interpolation weights.

CALL read_date_wrf(wout1,year1,month1,day1,hour1,minute1,second1)
CALL read_date_wrf(wout2,year2,month2,day2,hour2,minute2,second2)

!Get tai time
CALL com_utc2tai(year1,month1,day1,hour1,minute1,second1,tai1)
CALL com_utc2tai(year2,month2,day2,hour2,minute2,second2,tai2)

CALL com_utc2tai(intyear,intmonth,intday,inthour,intminute,intsecond,inttai)

if( inttai <= tai2 .AND. inttai >= tai1 )then
  WRITE(*,*)"We will interpolate data using the following times:"
  WRITE(*,*)"Initial time:"
  WRITE(*,*)year1,month1,day1,hour1,minute1,second1
  WRITE(*,*)"Final time:"
  WRITE(*,*)year2,month2,day2,hour2,minute2,second2
else
  WRITE(*,*)"[Error]: The requested time is outside the time window."
  WRITE(*,*)"Initial time:"
  WRITE(*,*)year1,month1,day1,hour1,minute1,second1
  WRITE(*,*)"Final time:"
  WRITE(*,*)year2,month2,day2,hour2,minute2,second2
  STOP
endif

  w1=0.0d0
  w2=0.0d0

  if ( inttai - tai1 /= 0 )then
   w2 =  ( inttai - tai1 ) / ( tai2 - tai1 )
  else
   w2 = 0.0d0
  endif

  w1=1.0d0-w2

  WRITE(*,*)"Weigth for time 1 is : ",w1
  WRITE(*,*)"Weigth for time 2 is : ",w2

ALLOCATE( gues3d(nlon,nlat,nlev,nv3d,2),gues2d(nlon,nlat,nv2d,2) )

 !Get perturbation at the initial time in the window.
 CALL read_grd(wout1,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 CALL read_grd(wout2,gues3d(:,:,:,:,2),gues2d(:,:,:,2))

 !Time linear interpolation.
 gues3d(:,:,:,:,1)=gues3d(:,:,:,:,1)*w1+gues3d(:,:,:,:,2)*w2
 gues2d(:,:,:,1)  =gues2d(:,:,:,1)*w1  +gues2d(:,:,:,2)*w2

 CALL write_grd(woutint,gues3d(:,:,:,:,1),gues2d(:,:,:,1))
 CALL write_date_wrf(woutint,intyear,intmonth,intday,inthour,intminute,intsecond)

END PROGRAM INTERP_WRFOUT

