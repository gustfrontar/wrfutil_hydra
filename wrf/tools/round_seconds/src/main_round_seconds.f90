PROGRAM ROUND_SECONDS
!=======================================================================
!
! [PURPOSE:] Generate an ensemble of perturbed WRF runs.
!   Round the date within WRF output files to the nearest integer minute.
!   We round to the nearest lower minute.
! [HISTORY:]
!   21/10/2016 Juan Ruiz created
!=======================================================================
USE common_round_seconds

 IMPLICIT NONE
 CHARACTER(15)            :: input_file='input_file.nc'
 CHARACTER(19)            :: my_date
 INTEGER                  :: ncid
 INTEGER                  :: iy,im,id,ih,imn
 REAL(r_size)             :: sec


 CALL get_date(input_file,my_date)

 READ(my_date(1:4),*)iy
 READ(my_date(6:7),*)im
 READ(my_date(9:10),*)id
 READ(my_date(12:13),*)ih
 READ(my_date(15:16),*)imn
 READ(my_date(18:19),*)sec

 WRITE(my_date(18:19),*)"00"

 CALL put_date(input_file,my_date)

END PROGRAM ROUND_SECONDS

