PROGRAM update_wrf_datetime
!=======================================================================
!
! [PURPOSE:] Update times in wrf files.
!
! [HISTORY:]
!   09/29/2023 Takemasa Miyoshi  created
!
!=======================================================================
  USE common
  USE common_wrf

  IMPLICIT NONE
  INTEGER :: ierr
  CHARACTER(100) ::wrf_file
  CHARACTER(19) :: new_date
  CHARACTER(100) ::arg
!-----------------------------------------------------------------------
! Updating the date
!-----------------------------------------------------------------------
  CALL GETARG ( 1, arg )    !Wrf_file to be updated
  wrf_file = arg
  CALL GETARG ( 2, arg )    !New date in the format YYYY-MM-DD_HH:MN:SS
  new_date=  TRIM(arg)

  WRITE(*,*)'Updating times for file ',wrf_file ,' new date will be ',new_date

  CALL write_date_wrf(wrf_file,new_date)

STOP
END PROGRAM update_wrf_datetime
